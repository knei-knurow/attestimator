// Package attestimator is for Szym.
package attestimator

import "math"

const (
	// If an acc measurement has norm-squared less than this, then it is considered to be faulty and the measurement is discarded.
	AccToleranceSquared = 1e-12 * 1e-12

	// If an acc-only generated quaternion `Qy` has norm-squared less than this, then the quaternion is considered to be zero and a fallback solution is used.
	QyNormToleranceSquared = 1e-12 * 1e-12

	// If a supposedly near-unit quaternion has norm-squared less than this during normalisation, then we declare a general state of emergency and completely reset the estimator.
	QhatNormToleranceSquared = 1e-12 * 1e-12

	// If the norm-squared of the calculated `xG` and/or `yG` basis vectors (`= norm(m_magCalib)*norm(mag)*sin(angle(acc,mag))`) is less than this prior to basis normalisation, then the calculation is assumed to be faulty due to acc/mag collinearity, and the mag measurement is discarded. Note however that this constant is also used for various checks in the `updateQy()` special cases.
	XGYGNormToleranceSquared = 1e-12 * 1e-12

	// If the norm-squared of the `Qhat` vector with its `x` and `y` components zeroed out is less than this, then normalisation is avoided and a fallback method is used instead to generate a unit pure yaw quaternion in the fused yaw acc-only case.
	WEZENormToleranceSquared = 1e-12 * 1e-12

	// If the absolute value of components of the `zGhat` vector are less than this, then they are considered to be zero (used deep inside the near-impossible special cases only).
	ZGHATAbsTolerance = 1e-12
)

type AccMethod int

// Acc-only resolution method enumeration.
const (
	// Resolve acc-only cases into full 3D orientations using the relative fused yaw method (default).
	MethodFusedYaw AccMethod = iota

	// Resolve acc-only cases into full 3D orientations using the absolute fused yaw method.
	MethodAbsFusedYaw

	// Resolve acc-only cases into full 3D orientations using the relative ZYX yaw method.
	MethodZYXyaw

	// Total number of acc-only resolution methods.
	MethodCount

	// Default acc-only resolution method (MethodFusedYaw)
	MethodDefault = MethodFusedYaw
)

type AttitudeEstimator struct {
	// Acc-only measurement resolution method
	accMethod AccMethod // The method to use if we cannot use a mag measurement to resolve the corresponding acc measurement into a full 3D orientation

	// Configuration variables
	m_Kp      float64 // Kp for standard operation (at lambda = 1)
	m_Ti      float64 // Ti for standard operation (at lambda = 1)
	m_KpQuick float64 // Kp for quick learning (at lambda = 0)
	m_TiQuick float64 // Ti for quick learning (at lambda = 0)
	m_QLTime  float64 // The nominal duration of quick learning when lambda is reset to zero

	// Magnetometer calibration
	m_magCalib [3]float64 // The mag value corresponding to/calibrated at the identity orientation

	// Internal variables
	m_lambda                             float64    // The quick learning parameter (0 => Quick learning, 1 => Normal operation, In-between => Faded PI gains with m_lambda as the interpolation constant)
	m_Qy, m_Qtilde, m_Qhat               [4]float64 // Quaternions: Format is (q0,qvec) = (w,x,y,z), m_Qhat must *always* be a unit quaternion (or within a few eps of it)!
	m_dQ, m_dQold                        [4]float64 // Quaternion derivatives: Format is (dw,dx,dy,dz)
	m_w, m_wold, m_bhat, m_omega, m_base [3]float64 // 3D vectors: Format is (x,y,z)
	m_Ry                                 [9]float64 // Rotation matrix: The numbering of the 3x3 matrix entries is left to right, top to bottom
	m_Ehat                               [3]float64 // Euler angles: Format is (psi,theta,phi) = (yaw,pitch,roll) and follows the ZYX convention
	m_Fhat                               [3]float64 // Fused angles: Format is (psi,theta,phi) = (yaw,pitch,roll) and follows the standard definition of fused angles
	m_FhatHemi                           bool       // Fused angles: Extra fourth hemisphere parameter, true is taken to mean 1 (positive z hemisphere), and false is taken to mean -1 (negative z hemisphere)
	m_eulerValid, m_fusedValid           bool       // Flags specifying whether the currently stored Euler and fused angle representations are up to date with the quaternion stored in m_Qhat
}

// Update updates the current attitude estimate based on new sensor measurements
// and the amount of elapsed time since the last update.
//
// This function requires (in a logically grouped sense) four inputs:
//
// - gyro: New gyroscope reading (rad/s)
//
// - dt: The time since the last call to `update()` (s) (in the first call should just be the nominal loop
// - period, must be in seconds)
//
// - acc: New accelerometer reading (can be in any self-consistent units, preferably ms)
//
// - mag: New magnetometer reading (can be in any self-consistent units, preferably gauss
// - and preferably the same units as `MagCalib`, but this is not required)
//
// The algorithm at work is based on the Passive Complementary Filter described in:
//
// - R. Mahoney, T. Hamel, and J.-M. Pflimlin, "Nonlinear complementary filters on the special orthogonal group",
//
// - IEEE Transactions on Automatic Control, vol. 53, no. 5, pp. 1203-1218.
//
// This implementation has been made ultimately robust, with detailed consideration and mathematical analysis having been performed to avoid
// all possible numerical issues, and to deal with all possible special cases. Nevertheless, the implementation has simultaneously been
// optimised to use the least number of floating point operations as possible. The source code of this and the subordinate functions thereof
// has been heavily documented.
func (e *AttitudeEstimator) Update(dt, gyroX, gyroY, gyroZ, accX, accY, accZ, magX, magY, magZ float64) {
	// Do not modify *anything* unless you know *exactly* what you're doing (even if you think you're not changing anything significant),
	// or it could affect the unconditional and unadulterated robustness of the entire estimator. This applies to updateQy() as well.

	// Update lambda
	if e.m_lambda < 1 {
		e.m_lambda += dt / e.m_QLTime
		if e.m_lambda <= 0 {
			e.m_lambda = 0
		}
		if e.m_lambda >= 1 {
			e.m_lambda = 1
		}
	}

	// Calculate the required filter gains for this update (Note: Ki = Kp / Ti)
	Kp := e.m_lambda*e.m_Kp + (1-e.m_lambda)*e.m_KpQuick
	Ti := e.m_lambda*e.m_Ti + (1-e.m_lambda)*e.m_TiQuick

	// Save the old values of the required variables
	e.m_wold[0] = e.m_w[0]
	e.m_wold[1] = e.m_w[1]
	e.m_wold[2] = e.m_w[2]
	e.m_dQold[0] = e.m_dQ[0]
	e.m_dQold[1] = e.m_dQ[1]
	e.m_dQold[2] = e.m_dQ[2]
	e.m_dQold[3] = e.m_dQ[3]

	// Calculate Qy, the current quaternion orientation measurement, based on the acc and mag readings
	e.updateQy(accX, accY, accZ, magX, magY, magZ) // Writes to m_Qy internally, and never fails...

	// Calculate the rotational error between the current Qhat and the new measured Qy
	e.m_Qtilde[0] = e.m_Qhat[0]*e.m_Qy[0] + e.m_Qhat[1]*e.m_Qy[1] + e.m_Qhat[2]*e.m_Qy[2] + e.m_Qhat[3]*e.m_Qy[3]
	e.m_Qtilde[1] = e.m_Qhat[0]*e.m_Qy[1] - e.m_Qhat[1]*e.m_Qy[0] - e.m_Qhat[2]*e.m_Qy[3] + e.m_Qhat[3]*e.m_Qy[2]
	e.m_Qtilde[2] = e.m_Qhat[0]*e.m_Qy[2] + e.m_Qhat[1]*e.m_Qy[3] - e.m_Qhat[2]*e.m_Qy[0] - e.m_Qhat[3]*e.m_Qy[1]
	e.m_Qtilde[3] = e.m_Qhat[0]*e.m_Qy[3] - e.m_Qhat[1]*e.m_Qy[2] + e.m_Qhat[2]*e.m_Qy[1] - e.m_Qhat[3]*e.m_Qy[0]

	// Calculate the angular velocity feedback term required to act in the direction of reducing Qtilde
	wscale := 2 * Kp * e.m_Qtilde[0]
	e.m_w[0] = wscale * e.m_Qtilde[1]
	e.m_w[1] = wscale * e.m_Qtilde[2]
	e.m_w[2] = wscale * e.m_Qtilde[3]

	// Update the estimated gyro bias (trapezoidal integration of -Ki*w)
	bscale := 0.5 * dt / Ti
	e.m_bhat[0] -= bscale * (e.m_w[0] + e.m_wold[0])
	e.m_bhat[1] -= bscale * (e.m_w[1] + e.m_wold[1])
	e.m_bhat[2] -= bscale * (e.m_w[2] + e.m_wold[2])

	// Calculate the required (combined predictive/corrective) angular velocity to apply to our current attitude estimate
	e.m_omega[0] = gyroX - e.m_bhat[0] + e.m_w[0]
	e.m_omega[1] = gyroY - e.m_bhat[1] + e.m_w[1]
	e.m_omega[2] = gyroZ - e.m_bhat[2] + e.m_w[2]

	// Convert the calculated angular velocity into a quaternion velocity (the missing factor of 0.5 here has been taken into dscale below)
	e.m_dQ[0] = -e.m_Qhat[1]*e.m_omega[0] - e.m_Qhat[2]*e.m_omega[1] - e.m_Qhat[3]*e.m_omega[2]
	e.m_dQ[1] = e.m_Qhat[0]*e.m_omega[0] + e.m_Qhat[2]*e.m_omega[2] - e.m_Qhat[3]*e.m_omega[1]
	e.m_dQ[2] = e.m_Qhat[0]*e.m_omega[1] - e.m_Qhat[1]*e.m_omega[2] + e.m_Qhat[3]*e.m_omega[0]
	e.m_dQ[3] = e.m_Qhat[0]*e.m_omega[2] + e.m_Qhat[1]*e.m_omega[1] - e.m_Qhat[2]*e.m_omega[0]

	// Update the attitude estimate using the calculated quaternion velocity (trapezoidal integration of dQ)
	dscale := 0.25 * dt // The extra factor of 0.5 here comes from the omission thereof in the calculation above
	e.m_Qhat[0] += dscale * (e.m_dQ[0] + e.m_dQold[0])
	e.m_Qhat[1] += dscale * (e.m_dQ[1] + e.m_dQold[1])
	e.m_Qhat[2] += dscale * (e.m_dQ[2] + e.m_dQold[2])
	e.m_Qhat[3] += dscale * (e.m_dQ[3] + e.m_dQold[3])

	// Renormalise the current attitude estimate
	qscale := e.m_Qhat[0]*e.m_Qhat[0] + e.m_Qhat[1]*e.m_Qhat[1] + e.m_Qhat[2]*e.m_Qhat[2] + e.m_Qhat[3]*e.m_Qhat[3]
	if qscale < QhatNormToleranceSquared {
		// The quaternion is so far away from being normalised (almost zero norm when it should be 1)
		// that something must be dreadfully wrong... (avoid potential division by zero below)
		e.Reset(true, true)
		return
	}
	qscale = 1 / math.Sqrt(qscale)
	e.m_Qhat[0] *= qscale
	e.m_Qhat[1] *= qscale
	e.m_Qhat[2] *= qscale
	e.m_Qhat[3] *= qscale

	// Reset the alternative representation validity flags
	e.m_eulerValid = false
	e.m_fusedValid = false
}

// updateQy calculates a 'measured' orientation based on given sensor measurements.
func (e *AttitudeEstimator) updateQy(accX, accY, accZ, magX, magY, magZ float64) {
	// This function generates a quaternion Qy (saved in the class member variable m_Qy) representing
	// the measured orientation of the body (i.e. the body-fixed coordinate system) relative to the
	// global coordinate system. It is referred to as the 'measured' orientation of the body, as the
	// quaternion derives from the given acc and mag measurements. The acc vector is guaranteed to
	// be respected, no matter what the inputs to this function and the current state estimate are.
	// That is to say, the generated Qy is guaranteed to be such that if the robot is exactly at the
	// rotation Qy, and gravity points down the negative global z-axis, then the given value of acc
	// is analytically expected to be measured in terms of body-fixed coordinates. The mag vector is
	// respected as long as it can be used to resolve the acc vector into a full 3D orientation, given
	// the current magnetometer calibration vector.
	//
	// This function uses the value of m_magCalib = (mtx,mty,mtz). Only the projection of this vector
	// onto the global xG-yG plane is considered however, which is equivalent to (mtx,mty,0). As such,
	// only the values of mtx and mty are used by this function, and m_magCalib is effectively taken
	// to be (mtx,mty,0), irrespective of what the value of mtz is. The only difference this makes is
	// in the comments, where for example ||m_magCalib|| is used as shorthand for sqrt(mtx^2 + mty^2).
	// In actual fact, all that really matters about the m_magCalib vector is its heading, because in
	// the calculation of Qy that's all it's effectively reduced to. More precisely stated, only the
	// value of atan2(mty,mtx) is important. The magnetometer readings do not affect any component of
	// the resulting rotation other than the yaw about the measured acceleration vector!
	//
	// The output may be Qy or -Qy. That is, Qy isn't forced to live in any particular half of
	// the quaternion space, and any code that uses the outputs of this function must be able to deal
	// with that. This is usually not a problem, as both quaternion multiplication and quaternion-
	// vector rotation implicitly handle this ambiguity.
	//
	// The rotation matrix m_Ry (if used) is stored as a 1D array of 9 doubles, arranged as follows:
	//        G    | Ry[0] Ry[1] Ry[2] |   | xBx yBx zBx |   | xGx xGy xGz |   B
	//   Ry =  R = | Ry[3] Ry[4] Ry[5] | = | xBy yBy zBy | = | yGx yGy yGz | =  R'
	//        B    | Ry[6] Ry[7] Ry[8] |   | xBz yBz zBz |   | zGx zGy zGz |   G
	//
	// Thus the standard orthonormal basis vectors of the global coordinate system written in the
	// body-fixed coordinate system are (G => Global coords, B => Body-fixed coords):
	//   xGhat = G(1,0,0) = B(xGx, xGy, xGz) = B(Ry[0], Ry[1], Ry[2])
	//   yGhat = G(0,1,0) = B(yGx, yGy, yGz) = B(Ry[3], Ry[4], Ry[5])
	//   zGhat = G(0,0,1) = B(zGx, zGy, zGz) = B(Ry[6], Ry[7], Ry[8])
	// Most of the algorithm cases in this function work by constructing the xGhat, yGhat and zGhat
	// vectors (in body-fixed coordinates). These are then used to populate the rotation matrix Ry,
	// which is then converted into the required quaternion Qy.
	//
	// We define the coordinate system H to be the coordinate system obtained when rotating the body-
	// fixed coordinate system B by the inverse of m_Qhat. In other words, H is our current estimate
	// of how the global coordinate system is oriented, without having considered the sensor
	// measurements of this cycle yet. The coordinate system G is defined as our 'measured' orientation
	// of the global csys, based on the sensor measurements of this cycle. The csys G *always* has zGhat
	// as its z-axis (i.e. respects the acc measurement). The orientations of the xG and yG axes are
	// then calculated either using the mag measurement, or m_Qhat and some assumption about the yaw.
	// Letting Qy represent the quaternion rotation from G to B (our 'measured' orientation of B for
	// this cycle, towards which m_Qhat limits via the filter equations) we have the following; taking
	// B to be our reference csys, G is the csys obtained by the inverse rotation of Qy, and H is the
	// csys obtained by the inverse rotation of m_Qhat.

	var dot, nxG, nyG float64

	// Calculate the norm (squared) of the acc measurement
	naccsq := accX*accX + accY*accY + accZ*accZ // = ||acc||^2

	// If the acc square norm is too close to zero then we actually have no idea where we are, so just set Qy to the identity quaternion and return
	if naccsq < AccToleranceSquared {
		e.m_Qy[0] = 1 // Qy -> w
		e.m_Qy[1] = 0 // Qy -> x
		e.m_Qy[2] = 0 // Qy -> y
		e.m_Qy[3] = 0 // Qy -> z
		return
	}

	// Define zGhat as the unit vector pointing in the direction of acc (the opposite direction to which gravity is pointing)
	// Note: To machine precision the resulting zGhat will have unit norm, i.e. ||zGhat|| == 1
	naccinv := 1 / math.Sqrt(naccsq)
	e.m_Ry[6] = accX * naccinv // zGhat -> x = zGx
	e.m_Ry[7] = accY * naccinv // zGhat -> y = zGy
	e.m_Ry[8] = accZ * naccinv // zGhat -> z = zGz

	// Project mag into the xG-yG plane (in body-fixed coordinates) using the formula mtilde := mag - dot(mag,zGhat)*zGhat
	// Note: To machine precision the resulting mtilde is perpendicular to zGhat, and ||mtilde|| = ||mag||*sin(angle(acc,mag)),
	//       where the angle is taken to be in the range [0,pi].
	dot = magX*e.m_Ry[6] + magY*e.m_Ry[7] + magZ*e.m_Ry[8] // dot(mag,zGhat)
	magX -= dot * e.m_Ry[6]                                // mtilde -> x
	magY -= dot * e.m_Ry[7]                                // mtilde -> y
	magZ -= dot * e.m_Ry[8]                                // mtilde -> z

	// Generate a second orthogonal basis vector for the xG-yG plane, complementary to mtilde, using the cross product (defines the csys X = m_base / Y = mtilde / Z = zGhat)
	// Note: To machine precision ||m_base|| == ||mtilde|| and dot(m_base,mtilde) == 0
	e.m_base[0] = magY*e.m_Ry[8] - magZ*e.m_Ry[7] // m_base -> x
	e.m_base[1] = magZ*e.m_Ry[6] - magX*e.m_Ry[8] // m_base -> y
	e.m_base[2] = magX*e.m_Ry[7] - magY*e.m_Ry[6] // m_base -> z

	// Calculate orthogonal xG and yG such that mtilde is collinear with m_magCalib (a vector in xG-yG-zG coordinates) projected onto the xG-yG plane
	// Note: To machine precision ||xG|| == ||yG|| == ||m_magCalib||*||mtilde|| == ||m_magCalib||*||mag||*sin(angle(acc,mag)), where the
	//       angle is taken to be in the range [0,pi]. Zero xG/yG arise iff m_magCalib is zero, or acc and mag are collinear.
	e.m_Ry[0] = e.m_magCalib[1]*e.m_base[0] + e.m_magCalib[0]*magX // xG -> x
	e.m_Ry[1] = e.m_magCalib[1]*e.m_base[1] + e.m_magCalib[0]*magY // xG -> y
	e.m_Ry[2] = e.m_magCalib[1]*e.m_base[2] + e.m_magCalib[0]*magZ // xG -> z
	e.m_Ry[3] = e.m_magCalib[1]*magX - e.m_magCalib[0]*e.m_base[0] // yG -> x
	e.m_Ry[4] = e.m_magCalib[1]*magY - e.m_magCalib[0]*e.m_base[1] // yG -> y
	e.m_Ry[5] = e.m_magCalib[1]*magZ - e.m_magCalib[0]*e.m_base[2] // yG -> z

	// Calculate the xG and yG vector (square) norms
	// Note: The calculated nxG and nyG should theoretically be identical.
	nxG = e.m_Ry[0]*e.m_Ry[0] + e.m_Ry[1]*e.m_Ry[1] + e.m_Ry[2]*e.m_Ry[2] // = ||xG||^2
	nyG = e.m_Ry[3]*e.m_Ry[3] + e.m_Ry[4]*e.m_Ry[4] + e.m_Ry[5]*e.m_Ry[5] // = ||yG||^2

	// Check whether the basis vector generation was successful (non-zero xG/yG)
	// Note: We check both nxG and nyG, even though to machine precision we should have nxG == nyG.
	if nxG < XGYGNormToleranceSquared || nyG < XGYGNormToleranceSquared {
		// This IF block executes if the algorithm was unable to resolve the acc and mag measurements into a consistent and unambiguous
		// 3D rotation. Assuming the m_magCalib is not the zero vector, this occurs iff acc and mag are collinear. In this case the mag
		// vector gives us no valid information and we are forced to discard it. We must now generate a full 3D rotation from just the
		// acc vector (i.e. from zGhat), using some additional assumption in order to be able to resolve the global yaw. This additional
		// assumption should ideally make the 'measured' orientation match as closely as possible to the current orientation estimate,
		// while respecting the zGhat vector.

		// Use the required acc-only resolution method to generate a full 3D orientation measurement Qy from zGhat (=> m_Ry[6-8]) and m_Qhat
		if e.accMethod == MethodFusedYaw {
			//
			// Fused Yaw Method
			//

			// We attempt to use the current orientation estimate m_Qhat to resolve the zGhat vector into a complete 3D rotation that agrees
			// "as much as possible" in a relative sense with m_Qhat. The code in this IF block assumes that m_Qhat is a unit quaternion.
			// Within reason this code can deal with numeric deviations thereof, but eventually these act to reduce the numerical accuracy of
			// the computed values. The term "fused yaw" refers to the yaw (about the z-axis) as defined by the first component of the fused
			// angles representation.

			// If the resolution of zGhat using m_Qhat is successful, then a quaternion is generated directly and returned. If on the other
			// hand the resolution attempt fails, a fallback solution is used, and this IF block outputs orthogonal xG (stored in m_Ry[0-2])
			// and yG (stored in m_Ry[3-5]), and the respective square norms of each, nxG = ||xG||^2 and nyG = ||yG||^2, such that xG-yG-zGhat
			// defines the 'measured' global coordinate system relative to the body-fixed coordinate system (the inverse rotation of m_Qy).

			// Rotate zGhat = BzG (in B coordinates) by m_Qhat = HqB to get the same vector HzG, but expressed in H coordinates
			// Note: HzG = m_Qhat * BzG * conj(m_Qhat), and assuming m_Qhat is a unit quaternion we have ||HzG|| = ||BzG|| as expected.
			perpx := 2 * (e.m_Qhat[2]*e.m_Ry[8] - e.m_Qhat[3]*e.m_Ry[7])                  // perp -> x
			perpy := 2 * (e.m_Qhat[3]*e.m_Ry[6] - e.m_Qhat[1]*e.m_Ry[8])                  // perp -> y
			perpz := 2 * (e.m_Qhat[1]*e.m_Ry[7] - e.m_Qhat[2]*e.m_Ry[6])                  // perp -> z
			HzGx := e.m_Ry[6] + e.m_Qhat[0]*perpx - e.m_Qhat[3]*perpy + e.m_Qhat[2]*perpz // HzG -> x
			HzGy := e.m_Ry[7] + e.m_Qhat[3]*perpx + e.m_Qhat[0]*perpy - e.m_Qhat[1]*perpz // HzG -> y
			HzGz := e.m_Ry[8] - e.m_Qhat[2]*perpx + e.m_Qhat[1]*perpy + e.m_Qhat[0]*perpz // HzG -> z

			// Calculate the Qytilde, the unnormalised quaternion rotation from G to B, such that G and H have no fused yaw relative to one another
			// Note: We have that the required (unnormalised) GqHtilde = (1+HzGz, HzGy, -HzGx, 0), which can easily be seen to have no fused yaw.
			//       We then have Qytilde = GqBtilde = GqHtilde * HqB = (1+HzGz, HzGy, -HzGx, 0) * m_Qhat.
			HzGztilde := 1 + HzGz
			e.m_Qy[0] = e.m_Qhat[0]*HzGztilde - e.m_Qhat[1]*HzGy + e.m_Qhat[2]*HzGx // Qytilde -> w
			e.m_Qy[1] = e.m_Qhat[1]*HzGztilde + e.m_Qhat[0]*HzGy - e.m_Qhat[3]*HzGx // Qytilde -> x
			e.m_Qy[2] = e.m_Qhat[2]*HzGztilde - e.m_Qhat[3]*HzGy - e.m_Qhat[0]*HzGx // Qytilde -> y
			e.m_Qy[3] = e.m_Qhat[3]*HzGztilde + e.m_Qhat[2]*HzGy + e.m_Qhat[1]*HzGx // Qytilde -> z

			// Calculate the square norm of the generated quaternion Qy
			// Note: This should only ever be zero if m_Qhat is zero, or if HzG = (0,0,-1).
			nQysq := e.m_Qy[0]*e.m_Qy[0] + e.m_Qy[1]*e.m_Qy[1] + e.m_Qy[2]*e.m_Qy[2] + e.m_Qy[3]*e.m_Qy[3] // = |Qy|^2

			// Normalise and return the generated quaternion, if possible, or use a fallback solution to calculate an alternative 'measured' acc-only quaternion Qy
			if nQysq >= QyNormToleranceSquared {
				// Normalise Qy and return
				nQyinv := 1 / math.Sqrt(nQysq) // = 1/|Qy|
				e.m_Qy[0] *= nQyinv            // Qy -> w
				e.m_Qy[1] *= nQyinv            // Qy -> x
				e.m_Qy[2] *= nQyinv            // Qy -> y
				e.m_Qy[3] *= nQyinv            // Qy -> z
				return
			} else {
				// This case only executes if the z-axes of G and H are exactly opposite, hence bringing the fused yaw representation of the relative
				// rotation between G and H to a singularity point. Said in another way, the above generation of m_Qy only fails if m_Qhat is zero
				// (assumed not to ever be the case), acc is zero (in which case the above code wouldn't have executed anyway because that case is
				// checked earlier), and/or HzG = (0,0,-1) (implying that zGhat and zH point in exactly opposite directions). So in essence, this case
				// only executes if the acc vector points *exactly* in the wrong direction according to the current orientation estimate m_Qhat.

				// Calculate orthogonal xG and yG (orthogonal to each other and zGhat) such that the ZYX yaw of H with respect to G is zero (as opposed to the fused yaw now)
				// Note: The critical underlying observation behind this efficient calculation is that the ZYX yaw of a rotation rotates the global
				//       x-axis so that it is collinear with the [global] projection of the body-fixed (i.e. rotated) x-axis into the [global] xy-plane.
				//       As such, the problem becomes finding an xG such that xG projected into the xG-yG plane (== xG) is equal to xh (the global x-axis
				//       expressed in body-fixed coordinates if the body is at an orientation of m_Qhat) projected into the xG-yG plane. Hence we never
				//       actually need to calculate a ZYX yaw, and can instead just set xG to be exactly the projection of xh into the plane perpendicular
				//       to zGhat. If m_Qhat is a unit quaternion, to machine precision ||0.5*xh|| = 0.5 and ||xG|| = ||yG|| = 0.5*sin(angle(xh,zGhat)),
				//       where the angle is taken to be in the range [0,pi]. Zero xG/yG arise iff xh and zGhat/acc are collinear. yG is calculated as the
				//       cross product of zGhat and xG.
				e.m_Ry[0] = 0.5 - (e.m_Qhat[2]*e.m_Qhat[2] + e.m_Qhat[3]*e.m_Qhat[3]) // 0.5*xh -> x
				e.m_Ry[1] = e.m_Qhat[1]*e.m_Qhat[2] - e.m_Qhat[3]*e.m_Qhat[0]         // 0.5*xh -> y
				e.m_Ry[2] = e.m_Qhat[1]*e.m_Qhat[3] + e.m_Qhat[2]*e.m_Qhat[0]         // 0.5*xh -> z
				dot = e.m_Ry[0]*e.m_Ry[6] + e.m_Ry[1]*e.m_Ry[7] + e.m_Ry[2]*e.m_Ry[8] // = dot(0.5*xh,zGhat)
				e.m_Ry[0] -= dot * e.m_Ry[6]                                          // xG -> x
				e.m_Ry[1] -= dot * e.m_Ry[7]                                          // xG -> y
				e.m_Ry[2] -= dot * e.m_Ry[8]                                          // xG -> z
				e.m_Ry[3] = e.m_Ry[2]*e.m_Ry[7] - e.m_Ry[1]*e.m_Ry[8]                 // yG -> x
				e.m_Ry[4] = e.m_Ry[0]*e.m_Ry[8] - e.m_Ry[2]*e.m_Ry[6]                 // yG -> y
				e.m_Ry[5] = e.m_Ry[1]*e.m_Ry[6] - e.m_Ry[0]*e.m_Ry[7]                 // yG -> z

				// Calculate the xG and yG vector (square) norms
				// Note: The calculated nxG and nyG should theoretically be identical.
				nxG = e.m_Ry[0]*e.m_Ry[0] + e.m_Ry[1]*e.m_Ry[1] + e.m_Ry[2]*e.m_Ry[2] // = ||xG||^2
				nyG = e.m_Ry[3]*e.m_Ry[3] + e.m_Ry[4]*e.m_Ry[4] + e.m_Ry[5]*e.m_Ry[5] // = ||yG||^2

				// Check whether the basis vector generation was successful (non-zero xG/yG)
				// Note: We check both nxG and nyG, even though to machine precision we should have nxG == nyG.
				if nxG < XGYGNormToleranceSquared || nyG < XGYGNormToleranceSquared {
					// The fallback ZYX yaw based method only fails if m_Qhat is such that xh is collinear with acc. That is to say, zGhat is collinear
					// with xH. In this case it is impossible however for zGhat to also be antiparallel to zH, and so in all such cases we should have
					// that the initial fused yaw method solution worked. As such, this case should be impossible to reach. We handle it anyway though
					// so that there is no conceivable way that a division by zero can occur.

					// Set Qy to the identity quaternion and return
					e.m_Qy[0] = 1 // Qy -> w
					e.m_Qy[1] = 0 // Qy -> x
					e.m_Qy[2] = 0 // Qy -> y
					e.m_Qy[3] = 0 // Qy -> z
					return
				}
			}
		} else if e.accMethod == MethodAbsFusedYaw {
			//
			// Absolute Fused Yaw Method
			//

			// We attempt to use the current orientation estimate m_Qhat to resolve the zGhat vector into a complete 3D rotation that agrees
			// "as much as possible" in an absolute sense with m_Qhat. The code in this IF block assumes that m_Qhat is a unit quaternion.
			// Within reason this code can deal with numeric deviations thereof, but eventually these act to reduce the numerical accuracy of
			// the computed values. The term "fused yaw" refers to the yaw (about the z-axis) as defined by the first component of the fused
			// angles representation.

			// If the resolution of zGhat using m_Qhat is successful, then a quaternion is generated directly and returned. If on the other
			// hand the resolution attempt fails, a fallback solution is used, and this IF block outputs orthogonal xG (stored in m_Ry[0-2])
			// and yG (stored in m_Ry[3-5]), and the respective square norms of each, nxG = ||xG||^2 and nyG = ||yG||^2, such that xG-yG-zGhat
			// defines the 'measured' global coordinate system relative to the body-fixed coordinate system (the inverse rotation of m_Qy).

			var wEhat, zEhat float64

			// Calculate the quaternion that corresponds to the pure fused yaw component of m_Qhat
			// Note: If we let m_Qhat = (wE,xE,yE,zE), then the fused yaw component quaternion is just (wE,0,0,zE) normalised to unit magnitude.
			//       The resulting normalised quaternion is then (wEhat,0,0,zEhat). It is important to stress that this shortcut for extracting
			//       the yaw is in fact true and proven for fused yaw, but is *not* true for ZYX yaw, as is sometimes erroneously thought. Care
			//       needs to be taken in the vicinity of wE = zE = 0, to avoid divisions by zero.
			zEsq := e.m_Qhat[3] * e.m_Qhat[3]         // zE^2
			nwEzEsq := e.m_Qhat[0]*e.m_Qhat[0] + zEsq // wE^2 + zE^2 = ||(we,0,0,zE)||^2
			if nwEzEsq >= WEZENormToleranceSquared {
				nwEzEinv := 1 / math.Sqrt(nwEzEsq)
				wEhat = e.m_Qhat[0] * nwEzEinv
				zEhat = e.m_Qhat[3] * nwEzEinv
			} else {
				// <-- Too close to wE = zE = 0 to just apply normalisation => We are in the vicinity of fused yaw
				// instability of m_Qhat, and must resort to using the (in that region) stable ZYX yaw instead
				hZYXyawE := 0.5 * math.Atan2(e.m_Qhat[0]*e.m_Qhat[3]+e.m_Qhat[1]*e.m_Qhat[2], 0.5-(e.m_Qhat[2]*e.m_Qhat[2]+zEsq)) // Half of the ZYX yaw of m_Qhat
				wEhat = math.Cos(hZYXyawE)
				zEhat = math.Sin(hZYXyawE)
			}

			// Calculate a quaternion that when normalised corresponds to a Qy that has the same fused yaw as m_Qhat (more precisely, the same fused yaw as (wEhat,0,0,zEhat))
			// Note: Taking zGhat = (zGx,zGy,zGz), the quaternion that respects this zGhat and has zero fused yaw is given by the normalisation
			//       of (1+zGz,zGy,-zGx,0). Thus, the quaternion that respects the given zGhat, and has fused yaw equal to that of the pure yaw
			//       quaternion (wEhat,0,0,zEhat), is given by the normalisation of the quaternion product (wEhat,0,0,zEhat)*(1+zGz,zGy,-zGx,0).
			//       This expands to Qytilde = (wEhat*(1+zGz), zGx*zEhat+zGy*wEhat, zGy*zEhat-zGx*wEhat, zEhat*(1+zGz)).
			zGztilde := 1 + e.m_Ry[8]                     // 1 + zGz
			e.m_Qy[0] = wEhat * zGztilde                  // Qytilde -> w (unnormalised)
			e.m_Qy[1] = wEhat*e.m_Ry[7] + zEhat*e.m_Ry[6] // Qytilde -> x (unnormalised)
			e.m_Qy[2] = zEhat*e.m_Ry[7] - wEhat*e.m_Ry[6] // Qytilde -> y (unnormalised)
			e.m_Qy[3] = zEhat * zGztilde                  // Qytilde -> z (unnormalised)

			// Calculate the square norm of the generated quaternion
			// Note: We should have nQysq = ||Qytilde||^2 = 1 + 2*zGz + zGx^2 + zGy^2 + zGz^2 = 2*(1+zGz) as zGhat is a unit vector.
			nQysq := e.m_Qy[0]*e.m_Qy[0] + e.m_Qy[1]*e.m_Qy[1] + e.m_Qy[2]*e.m_Qy[2] + e.m_Qy[3]*e.m_Qy[3]

			// Normalise and return the generated quaternion, if possible, or use a fallback solution to calculate an alternative 'measured' acc-only quaternion Qy
			if nQysq >= QyNormToleranceSquared {
				// Normalise m_Qy and return
				nQyinv := 1 / math.Sqrt(nQysq)
				e.m_Qy[0] *= nQyinv // Qy -> w
				e.m_Qy[1] *= nQyinv // Qy -> x
				e.m_Qy[2] *= nQyinv // Qy -> y
				e.m_Qy[3] *= nQyinv // Qy -> z
				return
			} else {
				// This case executes if the m_Qy that was generated using the method above was essentially zero, and thus couldn't be normalised
				// to obtain a unit 'measured' orientation as required. The square norm of Qytilde is expected to be 2*(1+zGz), and thus this case
				// invokes when zGz is very close to -1, that is, when zGhat is very close to (0,0,-1). We identify this situation as being exactly
				// when m_Qy is close to the fused yaw singularity. As the fused yaw goes unstable near this singularity, we calculate a fallback
				// solution by matching ZYX yaws (instead of matching fused yaws like before). That is, we calculate a rotation matrix m_Ry that
				// respects zGhat and has the same ZYX yaw as the quaternion (wEhat,0,0,zEhat) (but not in general the same ZYX yaw as m_Qhat!).
				// We define xh to be the projection of xBhat = B(1,0,0) onto the global xG-yG plane. This is the plane perpendicular to the known
				// vector zGhat. The angle from xG to xh about the zG axis is equivalent to the ZYX yaw of the body-fixed frame (i.e. of m_Qy).

				// Calculate xh in both the global and body-fixed coordinate systems
				// Note: We have that xh = G(cos(psiE),sin(psiE),0), where psiE is the (independent of definition) yaw of (wEhat,0,0,zEhat).
				//       We also have xh = xBhat - dot(xBhat,zGhat)*zGhat = B(1-zGx^2,-zGx*zGy,-zGx*zGz) and yh = cross(zGhat,xBhat) = B(0,zGz,-zGy).
				//       The magnitudes of xh and yh should to machine precision be ||xh||^2 = ||yh||^2 = zGy^2 + zGz^2, assuming that zGhat is a
				//       unit vector. The vectors xh, yh and zGhat together form an orthogonal basis. The double angle formulas for sin and cos are
				//       used to evaluate cos(psiE) and sin(psiE), using that we know that wEhat = cos(psiE/2) and zEhat = sin(psiE/2). Note that
				//       yh isn't actually explicitly calculated here, it is used implicitly in the next step.
				cpsiE := wEhat*wEhat - zEhat*zEhat // cos(psiE) = cos(psiE/2)^2 - sin(psiE/2)^2
				spsiE := 2 * wEhat * zEhat         // sin(psiE) = 2*sin(psiE/2)*cos(psiE/2)
				xhx := 1 - e.m_Ry[6]*e.m_Ry[6]     // xh -> x = 1 - zGx^2
				xhy := -e.m_Ry[6] * e.m_Ry[7]      // xh -> y = -zGx*zGy
				xhz := -e.m_Ry[6] * e.m_Ry[8]      // xh -> z = -zGx*zGz

				// Calculate orthogonal xG and yG such that the global and body-fixed representations of xh match (i.e. xh is xG yawed by psiE)
				// Note: The required equations are xG = xh*cos(psiE) - yh*sin(psiE) and yG = cross(zGhat,xG), with neither xG nor yG necessarily
				//       coming out a unit vector. Important however is that ||xh||^2 = ||yh||^2, implying that ||xG||^2 = A. As ||yG||^2 = ||xG||^2
				//       this means that to machine precision we have ||xG||^2 = ||yG||^2 = ||xh||^2 = ||yh||^2 = zGy^2 + zGz^2.
				e.m_Ry[0] = xhx * cpsiE                               // xG -> x
				e.m_Ry[1] = xhy*cpsiE - e.m_Ry[8]*spsiE               // xG -> y
				e.m_Ry[2] = xhz*cpsiE + e.m_Ry[7]*spsiE               // xG -> z
				e.m_Ry[3] = e.m_Ry[7]*e.m_Ry[2] - e.m_Ry[1]*e.m_Ry[8] // yG -> x
				e.m_Ry[4] = e.m_Ry[8]*e.m_Ry[0] - e.m_Ry[2]*e.m_Ry[6] // yG -> y
				e.m_Ry[5] = e.m_Ry[6]*e.m_Ry[1] - e.m_Ry[0]*e.m_Ry[7] // yG -> z

				// Calculate the xG and yG vector (square) norms
				// Note: The calculated nxG and nyG should theoretically be identical.
				nxG = e.m_Ry[0]*e.m_Ry[0] + e.m_Ry[1]*e.m_Ry[1] + e.m_Ry[2]*e.m_Ry[2] // = ||xG||^2
				nyG = e.m_Ry[3]*e.m_Ry[3] + e.m_Ry[4]*e.m_Ry[4] + e.m_Ry[5]*e.m_Ry[5] // = ||yG||^2

				// Check whether the basis vector generation was successful (non-zero xG/yG)
				// Note: We check both nxG and nyG, even though to machine precision we should have nxG == nyG.
				if nxG < XGYGNormToleranceSquared || nyG < XGYGNormToleranceSquared {

					// This case should be physically impossible, as for the matching of fused yaws to have failed we must have zGz very close to -1,
					// and for the matching of ZYX yaws to have failed we must have zGy^2 + zGz^2 very close to 0. It is clearly not possible for
					// both of these conditions to be true at the same time. To cover all bases though, and avoid absolutely all remotely conceivable
					// possibilities of division by zero, we handle this dud case separately.

					// Set Qy to the identity quaternion and return
					e.m_Qy[0] = 1 // Qy -> w
					e.m_Qy[1] = 0 // Qy -> x
					e.m_Qy[2] = 0 // Qy -> y
					e.m_Qy[3] = 0 // Qy -> z
					return
				}
			}
		} else if e.accMethod == MethodZYXyaw {
			//
			// ZYX Yaw Method
			//

			// We attempt to use the current orientation estimate m_Qhat to resolve the zGhat vector into a complete 3D rotation that agrees
			// "as much as possible" with m_Qhat. The code in this IF block assumes that m_Qhat is a unit quaternion. Within reason this code
			// can deal with numeric deviations thereof, but eventually these act to reduce the numerical accuracy of the computed values.
			// The terms "ZYX yaw" and "ZXY yaw" in the comments below refer to the yaw of a rotation (about the z-axis) in the case of ZYX
			// and ZXY Euler angles respectively.
			//
			// The output of this IF block in all cases is orthogonal xG (stored in m_Ry[0-2]) and yG (stored in m_Ry[3-5]), and the respective
			// square norms of each, nxG = ||xG||^2 and nyG = ||yG||^2, such that xG-yG-zGhat defines the 'measured' global coordinate system
			// relative to the body-fixed coordinate system (the inverse rotation of m_Qy).

			// Calculate orthogonal xG and yG (orthogonal to each other and zGhat) such that the ZYX yaw of H with respect to G is zero
			// Note: The critical underlying observation behind this efficient calculation is that the ZYX yaw of a rotation rotates the global
			//       x-axis so that it is collinear with the [global] projection of the body-fixed (i.e. rotated) x-axis into the [global] xy-plane.
			//       As such, the problem becomes finding an xG such that xG projected into the xG-yG plane (== xG) is equal to xh (the global x-axis
			//       expressed in body-fixed coordinates if the body is at an orientation of m_Qhat) projected into the xG-yG plane. Hence we never
			//       actually need to calculate a ZYX yaw, and can instead just set xG to be exactly the projection of xh into the plane perpendicular
			//       to zGhat. If m_Qhat is a unit quaternion, to machine precision ||0.5*xh|| = 0.5 and ||xG|| = ||yG|| = 0.5*sin(angle(xh,zGhat)),
			//       where the angle is taken to be in the range [0,pi]. Zero xG/yG arise iff xh and zGhat/acc are collinear. yG is calculated as the
			//       cross product of zGhat and xG.
			e.m_Ry[0] = 0.5 - (e.m_Qhat[2]*e.m_Qhat[2] + e.m_Qhat[3]*e.m_Qhat[3]) // 0.5*xh -> x
			e.m_Ry[1] = e.m_Qhat[1]*e.m_Qhat[2] - e.m_Qhat[3]*e.m_Qhat[0]         // 0.5*xh -> y
			e.m_Ry[2] = e.m_Qhat[1]*e.m_Qhat[3] + e.m_Qhat[2]*e.m_Qhat[0]         // 0.5*xh -> z
			dot = e.m_Ry[0]*e.m_Ry[6] + e.m_Ry[1]*e.m_Ry[7] + e.m_Ry[2]*e.m_Ry[8] // = dot(0.5*xh,zGhat)
			e.m_Ry[0] -= dot * e.m_Ry[6]                                          // xG -> x
			e.m_Ry[1] -= dot * e.m_Ry[7]                                          // xG -> y
			e.m_Ry[2] -= dot * e.m_Ry[8]                                          // xG -> z
			e.m_Ry[3] = e.m_Ry[2]*e.m_Ry[7] - e.m_Ry[1]*e.m_Ry[8]                 // yG -> x
			e.m_Ry[4] = e.m_Ry[0]*e.m_Ry[8] - e.m_Ry[2]*e.m_Ry[6]                 // yG -> y
			e.m_Ry[5] = e.m_Ry[1]*e.m_Ry[6] - e.m_Ry[0]*e.m_Ry[7]                 // yG -> z

			// Calculate the xG and yG vector (square) norms
			// Note: The calculated nxG and nyG should theoretically be identical.
			nxG = e.m_Ry[0]*e.m_Ry[0] + e.m_Ry[1]*e.m_Ry[1] + e.m_Ry[2]*e.m_Ry[2] // = ||xG||^2
			nyG = e.m_Ry[3]*e.m_Ry[3] + e.m_Ry[4]*e.m_Ry[4] + e.m_Ry[5]*e.m_Ry[5] // = ||yG||^2

			// Check whether the basis vector generation was successful (non-zero xG/yG)
			// Note: We check both nxG and nyG, even though to machine precision we should have nxG == nyG.
			if nxG < XGYGNormToleranceSquared || nyG < XGYGNormToleranceSquared {
				// This IF block executes if the mag vector had to be discarded (i.e. if m_magCalib is zero, or acc and mag are collinear) and
				// m_Qhat is such that xh is collinear with acc. That is, mag was discarded and our current attitude estimate places the global
				// x-axis in the same (or 180 degree opposite) direction as the acc we just measured. Assuming m_Qhat is a unit quaternion, this
				// can only happen if our estimate is *exactly* out by 90 degrees, and perchance in *exactly* the right direction. If m_Qhat
				// deviates from unit norm by more than a negligible eps-ish amount however (which it should never do), then funkier and slightly
				// counterintuitive things can happen that produce further special cases below. The xh (ZYX) and yh (ZXY) methods produce
				// identical results for collinear zGhat and zh. As the angle between these two vectors increases, the xh and yh methods
				// continuously (but only gradually) start to differ more and more in their output (in terms of global yaw only though, the zGhat
				// vector always points directly opposite to the measured acc). As the angle between zGhat (opposite of measured acc) and zh
				// (current estimate of the up direction in body-fixed coordinates) approaches 90 degrees, there are two singularities
				// (corresponding to gimbal lock of the two different Euler angle conventions) that make the output of the xh and yh methods
				// start to differ vastly and in an unstable way. For example, close to zGhat == xh the xh output becomes extremely sensitive to
				// small variations in zGhat (and nxG/nyG tend to zero, which is the condition that is checked for by this IF block), while the
				// yh output remains stable, and vice versa for zGhat == yh. We wish to have a continuous xG/yG output for as much of the domain
				// as possible though, so we nominally always use xh (as it corresponds to matching global yaw in a sense consistent with our
				// nominal ZYX Euler angle convention) unless it produces an effectively zero output, in which case we switch to yh, which in
				// this neighbourhood is guaranteed to be stable (with the unit quaternion technical caveat). The switching is of course a
				// locally discontinuous operation though, but we really have no choice here due to the nasty directional dependence of the
				// singularities. So in summary, the important points to take home from all this are:
				// 1) This yh method/IF block will practically never execute, but it cannot be neglected.
				// 2) The singularities and regions of progressive instability in xG/yG do NOT depend on the absolute orientation of acc and
				//    m_Qhat at all (i.e. it doesn't matter if we're upright, upside-down or sideways), it depends only on the *relative*
				//    orientation/rotation between acc and m_Qhat. What's more, it can only even possibly be a problem if acc and m_Qhat
				//    disagree in what direction is up by almost exactly 90 degrees, and by coincidence *exactly* in the wrong direction.

				// Calculate orthogonal xG and yG so that the ZXY yaw of m_Qy, the rotation defined by the resulting orthogonal basis xG-yG-zGhat, is equal to the ZXY yaw of m_Qhat
				// Note: The critical underlying observation behind this efficient calculation is that the ZXY yaw of a rotation rotates the global
				//       y-axis so that it is collinear with the [global] projection of the body-fixed (i.e. rotated) y-axis into the [global] xy-plane.
				//       As such, the problem becomes finding a yG such that yG projected into the xG-yG plane (== yG) is equal to yh (the global y-axis
				//       expressed in body-fixed coordinates if the body is at an orientation of m_Qhat) projected into the xG-yG plane. Hence we never
				//       actually need to calculate a ZXY yaw, and can instead just set yG to be exactly the projection of yh into the plane perpendicular
				//       to zGhat. If m_Qhat is a unit quaternion, to machine precision ||0.5*yh|| = 0.5 and ||yG|| = ||xG|| = 0.5*sin(angle(yh,zGhat)),
				//       where the angle is taken to be in the range [0,pi]. Zero xG/yG arise iff yh and zGhat/acc are collinear. xG is calculated as the
				//       cross product of yG and zGhat.
				e.m_Ry[3] = e.m_Qhat[1]*e.m_Qhat[2] + e.m_Qhat[3]*e.m_Qhat[0]         // 0.5*yh -> x
				e.m_Ry[4] = 0.5 - (e.m_Qhat[1]*e.m_Qhat[1] + e.m_Qhat[3]*e.m_Qhat[3]) // 0.5*yh -> y
				e.m_Ry[5] = e.m_Qhat[2]*e.m_Qhat[3] - e.m_Qhat[1]*e.m_Qhat[0]         // 0.5*yh -> z
				dot = e.m_Ry[3]*e.m_Ry[6] + e.m_Ry[4]*e.m_Ry[7] + e.m_Ry[5]*e.m_Ry[8] // = dot(0.5*yh,zGhat)
				e.m_Ry[3] -= dot * e.m_Ry[6]                                          // yG -> x
				e.m_Ry[4] -= dot * e.m_Ry[7]                                          // yG -> y
				e.m_Ry[5] -= dot * e.m_Ry[8]                                          // yG -> z
				e.m_Ry[0] = e.m_Ry[4]*e.m_Ry[8] - e.m_Ry[5]*e.m_Ry[7]                 // xG -> x
				e.m_Ry[1] = e.m_Ry[5]*e.m_Ry[6] - e.m_Ry[3]*e.m_Ry[8]                 // xG -> y
				e.m_Ry[2] = e.m_Ry[3]*e.m_Ry[7] - e.m_Ry[4]*e.m_Ry[6]                 // xG -> z

				// Calculate the xG and yG vector (square) norms
				// Note: The calculated nxG and nyG should theoretically be identical.
				nxG = e.m_Ry[0]*e.m_Ry[0] + e.m_Ry[1]*e.m_Ry[1] + e.m_Ry[2]*e.m_Ry[2] // = ||xG||^2
				nyG = e.m_Ry[3]*e.m_Ry[3] + e.m_Ry[4]*e.m_Ry[4] + e.m_Ry[5]*e.m_Ry[5] // = ||yG||^2

				// Check whether the basis vector generation was successful (non-zero xG/yG)
				// Note: We check both nxG and nyG, even though to machine precision we should have nxG == nyG.
				if nxG < XGYGNormToleranceSquared || nyG < XGYGNormToleranceSquared {
					// Ok, so you're asking me why I'm even checking this case, seeing as it's impossible anyway, right? Right...?
					// Well, almost. Somewhat surprisingly it *is* actually possible for the calculated xh and yh to turn out collinear, even though
					// the corresponding expressions used to calculate them are orthogonal - as long as m_Qhat is a unit quaternion! That exactly is
					// the catch. This entire function was written with the mindset that it should never ever break, even if it receives completely
					// rubbish inputs. In that case it should simply output whatever it thinks is most appropriate, and *never* start performing
					// divisions by zero or such. No matter what happens, this function *must* be guaranteed to return a valid non-Inf/NaN quaternion,
					// or the effects of one bad call could ripple through and break the entire attitude estimator, which is not robust. I've managed
					// to prove mathematically that the previously calculated xh and yh are collinear iff m_Qhat is a pure vector quaternion of norm
					// 1/sqrt(2). If, furthermore, the resulting collinear xh and yh also happen to be collinear with the measured acc, AND with the
					// measured mag (or if m_magCalib is zero), then this case is invoked! If this is the case, then we know that m_Qhat is useless,
					// and so we are forced to discard it also. In order to still be able to resolve our single acc measurement into a full 3D
					// orientation, some extra (arbitrary) assumption is required. The two assumptions in use below amount to assumptions of zero yaw
					// in terms of either ZYX yaw or ZXY yaw, with preference to the (nominal) former convention for Euler angles.

					// Check whether zGhat is collinear with (1,0,0), and produce an appropriate output either way, with the assumption of zero ZYX yaw or ZXY yaw
					// Note: If zGhat is collinear with (1,0,0), then the code in the first case doesn't work because (1,0,0) is the ZYX gimbal lock situation.
					//       In this case however we know that zGhat can't also be collinear with (0,1,0), the ZXY gimbal lock situation, so it is safe to
					//       calculate a new xG/yG pair based on a zero ZXY yaw assumption without any further checks.
					if math.Abs(e.m_Ry[7]) >= ZGHATAbsTolerance || math.Abs(e.m_Ry[8]) >= ZGHATAbsTolerance {
						// If zGhat is *not* collinear with (1,0,0)...
						// Assume zero ZYX yaw: xG is the projection of (1,0,0) into the plane perpendicular to zGhat, yG is the cross product of zGhat and xG
						// Note: To machine precision ||xG|| = ||yG|| = sin(angle(zGhat,(1,0,0))), where the angle is taken to be in the range [0,pi].
						e.m_Ry[0] = 1 - e.m_Ry[6]*e.m_Ry[6] // xG -> x
						e.m_Ry[1] = -e.m_Ry[6] * e.m_Ry[7]  // xG -> y
						e.m_Ry[2] = -e.m_Ry[6] * e.m_Ry[8]  // xG -> z
						e.m_Ry[3] = 0                       // yG -> x
						e.m_Ry[4] = e.m_Ry[8]               // yG -> y
						e.m_Ry[5] = -e.m_Ry[7]              // yG -> z
					} else {
						// If zGhat *is* collinear with (1,0,0), and hence not collinear with (0,1,0) (as it is a unit vector)...
						// Assume zero ZXY yaw: yG is the projection of (0,1,0) into the plane perpendicular to zGhat, xG is the cross product of yG and zGhat
						// Note: To machine precision ||xG|| = ||yG|| = sin(angle(zGhat,(0,1,0))), where the angle is taken to be in the range [0,pi].
						//       This case is only invoked if m_Qhat is either exactly (0,+-1/sqrt(2),0,0) or (0,0,0,+-1/sqrt(2)), acc is non-zero and
						//       along the body-fixed x-axis, and mag is collinear with acc or m_magCalib is zero.
						e.m_Ry[0] = e.m_Ry[8]                 // xG -> x
						e.m_Ry[1] = 0.0                       // xG -> y
						e.m_Ry[2] = -e.m_Ry[6]                // xG -> z
						e.m_Ry[3] = -e.m_Ry[7] * e.m_Ry[6]    // yG -> x
						e.m_Ry[4] = 1.0 - e.m_Ry[7]*e.m_Ry[7] // yG -> y
						e.m_Ry[5] = -e.m_Ry[7] * e.m_Ry[8]    // yG -> z
					}

					// Calculate the xG and yG vector (square) norms
					// Note: The calculated nxG and nyG should theoretically be identical.
					nxG = e.m_Ry[0]*e.m_Ry[0] + e.m_Ry[1]*e.m_Ry[1] + e.m_Ry[2]*e.m_Ry[2] // = ||xG||^2
					nyG = e.m_Ry[3]*e.m_Ry[3] + e.m_Ry[4]*e.m_Ry[4] + e.m_Ry[5]*e.m_Ry[5] // = ||yG||^2

					// Check whether the basis vector generation was successful (non-zero xG/yG)
					// Note: We check both nxG and nyG, even though to machine precision we should have nxG == nyG.
					if nxG < XGYGNormToleranceSquared || nyG < XGYGNormToleranceSquared {
						// This case should be physically impossible, i.e. both theoretically and numerically. To cover all bases though,
						// and avoid absolutely all remotely conceivable possibilities of division by zero, we handle this dud case separately.

						// Set Qy to the identity quaternion and return
						e.m_Qy[0] = 1 // Qy -> w
						e.m_Qy[1] = 0 // Qy -> x
						e.m_Qy[2] = 0 // Qy -> y
						e.m_Qy[3] = 0 // Qy -> z
						return
					}
				}
			}
		} else {
			// <-- Unrecognised acc-only resolution method, should never happen!
			// Set Qy to the identity quaternion and return
			e.m_Qy[0] = 1 // Qy -> w
			e.m_Qy[1] = 0 // Qy -> x
			e.m_Qy[2] = 0 // Qy -> y
			e.m_Qy[3] = 0 // Qy -> z
			return
		}
	}

	// Normalise xG and yG to obtain an orthonormal basis (in conjunction with zGhat) that forms the rows of the orthogonal rotation matrix m_Ry
	// Note: The calculated orthonormal basis <xGhat, yGhat, zGhat> is placed in the rows of this matrix (as opposed to the columns) as we want
	//       the inverse (i.e. transpose) rotation for m_Ry. That is, the global to body-fixed frame rotation, not the body-fixed to global rotation.
	nxG = 1 / math.Sqrt(nxG)
	nyG = 1 / math.Sqrt(nyG)
	e.m_Ry[0] *= nxG // xGhat -> x = xGx
	e.m_Ry[1] *= nxG // xGhat -> y = xGy
	e.m_Ry[2] *= nxG // xGhat -> z = xGz
	e.m_Ry[3] *= nyG // yGhat -> x = yGx
	e.m_Ry[4] *= nyG // yGhat -> y = yGy
	e.m_Ry[5] *= nyG // yGhat -> z = yGz

	// Declare variables
	var r, s, t float64

	// Convert the rotation matrix m_Ry into the quaternion m_Qy
	// Note: m_Qy and -m_Qy are both valid and completely equivalent outputs here, so we have
	//       the freedom to arbitrarily choose the sign of *one* of the quaternion parameters.
	t = e.m_Ry[0] + e.m_Ry[4] + e.m_Ry[8]
	if t >= 0.0 {
		// Option 1: Centred at identity rotation... [Condition ensures |w| >= 0.5, WLOG choose the sign w >= 0.5]
		r = math.Sqrt(1 + t)                    // = 2*|w|
		s = 0.5 / r                             // = 1/(4*|w|)
		e.m_Qy[0] = 0.5 * r                     // = |w|           = w*sgn(w) = w
		e.m_Qy[1] = s * (e.m_Ry[7] - e.m_Ry[5]) // = (4xw)/(4*|w|) = x*sgn(w) = x
		e.m_Qy[2] = s * (e.m_Ry[2] - e.m_Ry[6]) // = (4yw)/(4*|w|) = y*sgn(w) = y
		e.m_Qy[3] = s * (e.m_Ry[3] - e.m_Ry[1]) // = (4zw)/(4*|w|) = z*sgn(w) = z
	} else if (e.m_Ry[8] >= e.m_Ry[4]) && (e.m_Ry[8] >= e.m_Ry[0]) {
		// Option 2: Centred at 180 deg z-rotation... [Conditions ensure |z| > 0.5, WLOG choose the sign z > 0.5]
		r = math.Sqrt(1 - (e.m_Ry[0] + e.m_Ry[4] - e.m_Ry[8])) // = 2*|z|
		s = 0.5 / r                                            // = 1/(4*|z|)
		e.m_Qy[0] = s * (e.m_Ry[3] - e.m_Ry[1])                // = (4zw)/(4*|z|) = w*sgn(z) = w
		e.m_Qy[1] = s * (e.m_Ry[2] + e.m_Ry[6])                // = (4xz)/(4*|z|) = x*sgn(z) = x
		e.m_Qy[2] = s * (e.m_Ry[7] + e.m_Ry[5])                // = (4yz)/(4*|z|) = y*sgn(z) = y
		e.m_Qy[3] = 0.5 * r                                    // = |z|           = z*sgn(z) = z
	} else if e.m_Ry[4] >= e.m_Ry[0] {
		// Option 3: Centred at 180 deg y-rotation... [Conditions ensure |y| > 0.5, WLOG choose the sign y > 0.5]
		r = math.Sqrt(1 - (e.m_Ry[0] - e.m_Ry[4] + e.m_Ry[8])) // = 2*|y|
		s = 0.5 / r                                            // = 1/(4*|y|)
		e.m_Qy[0] = s * (e.m_Ry[2] - e.m_Ry[6])                // = (4yw)/(4*|y|) = w*sgn(y) = w
		e.m_Qy[1] = s * (e.m_Ry[3] + e.m_Ry[1])                // = (4xy)/(4*|y|) = x*sgn(y) = x
		e.m_Qy[2] = 0.5 * r                                    // = |y|           = y*sgn(y) = y
		e.m_Qy[3] = s * (e.m_Ry[7] + e.m_Ry[5])                // = (4yz)/(4*|y|) = z*sgn(y) = z
	} else {
		// Option 4: Centred at 180 deg x-rotation... [Conditions ensure |x| > 0.5, WLOG choose the sign x > 0.5]
		r = math.Sqrt(1 + (e.m_Ry[0] - e.m_Ry[4] - e.m_Ry[8])) // = 2*|x|
		s = 0.5 / r                                            // = 1/(4*|x|)
		e.m_Qy[0] = s * (e.m_Ry[7] - e.m_Ry[5])                // = (4xw)/(4*|x|) = w*sgn(x) = w
		e.m_Qy[1] = 0.5 * r                                    // = |x|           = x*sgn(x) = x
		e.m_Qy[2] = s * (e.m_Ry[3] + e.m_Ry[1])                // = (4xy)/(4*|x|) = y*sgn(x) = y
		e.m_Qy[3] = s * (e.m_Ry[2] + e.m_Ry[6])                // = (4xz)/(4*|x|) = z*sgn(x) = z
	}

	// Any deviations from being a unit quaternion (that might be experienced here due to the inaccuracies
	// of floating point arithmetic) are pretty much irrelevant. This is firstly because they will only
	// ever be extremely minor eps deviations, given the mathematical correctness of this algorithm, but also
	// because any scaling in Qy is swallowed up by the Kp in the expression for m_w in the update()
	// function anyway. This is why no quaternion normalisation step has been added here.
}
