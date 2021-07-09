// Package attestimator is for Szym.
package attestimator

import "math"

const (
	// If an acc measurement has norm-squared less than this, then it is considered to be faulty and the measurement is discarded.
	AccToleranceSquared = 1e-12 * 1e-12 //

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

// Acc-only resolution method enumeration.
const (
	// Default acc-only resolution method (`ME_FUSED_YAW`).
	MethodDefault = 0

	// Resolve acc-only cases into full 3D orientations using the relative fused yaw method (default).
	MethodFusedYaw = MethodDefault

	// Resolve acc-only cases into full 3D orientations using the absolute fused yaw method.
	MethodAbsFusedYaw = 1

	// Resolve acc-only cases into full 3D orientations using the relative ZYX yaw method.
	MethodZYXyaw = 2

	// Total number of acc-only resolution methods.
	MethodCount = 3
)

type AttitudeEstimator struct {
	// Acc-only measurement resolution method
	AccMethodEnum int // The method to use if we cannot use a mag measurement to resolve the corresponding acc measurement into a full 3D orientation

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

func (e *AttitudeEstimator) Reset(quickLearn bool, resetGyroBias bool) {
	e.ResetState(resetGyroBias)

	if quickLearn {
		e.SetLambda(0) // reset lambda
	} else {
		e.SetLambda(1)
	}
}

func (e *AttitudeEstimator) ResetState(resetGyroBias bool) {
	var i int

	// Initialise the attitude estimate
	e.SetAttitude(1, 0, 0, 0) // Resets m_Qhat, m_eulerValid and m_fusedValid internally

	// Update the alternative attitude estimate representations
	e.updateEuler() // Resets m_Ehat and m_eulerValid internally
	e.updateFused() // Resets m_Fhat, m_FhatHemi and m_fusedValid internally

	// Initialise the gyro bias estimate
	if resetGyroBias {
		e.SetGyroBias(0, 0, 0) // Resets m_bhat internally
	}

	// Initialise the remaining size 3 internal variables
	for i = 0; i < 3; i++ {
		e.m_w[i] = 0
		e.m_wold[i] = 0
		e.m_omega[i] = 0
		e.m_base[i] = 0
	}

	// Initialise the remaining size 4 internal variables
	for i = 0; i < 4; i++ {
		e.m_Qy[i] = 0
		e.m_Qtilde[i] = 0
		e.m_dQ[i] = 0
		e.m_dQold[i] = 0
	}
	e.m_Qy[0] = 1
	e.m_Qtilde[0] = 1

	// Initialise the remaining size 9 internal variables
	for i = 0; i < 9; i++ {
		e.m_Ry[i] = 0
	}

	e.m_Ry[0] = 1
	e.m_Ry[4] = 1
	e.m_Ry[8] = 1
}

// SetLambda sets estimator's lambda to the desired value, by default this is `1.0`, which turns
// quick learning off. Values outside of `[0,1]` are coerced.
func (e *AttitudeEstimator) SetLambda(value float64) {
	if value >= 1 {
		e.m_lambda = 1
	} else {
		if value <= 0 {
			e.m_lambda = 0
		} else {
			e.m_lambda = value
		}
	}
}

// SetAttitude sets the current attitude estimate.
func (e *AttitudeEstimator) SetAttitude(w, x, y, z float64) {
	var qscale float64 = w*w + x*x + y*y + z*z // calculate the quaternion square norm

	// update the current attitude estimate
	if qscale < QhatNormToleranceSquared {
		// reset the attitude estimate to the identity orientation if the norm is too close to zero
		e.m_Qhat[0] = 1
		e.m_Qhat[1] = 0
		e.m_Qhat[2] = 0
		e.m_Qhat[3] = 0
	} else {
		qscale = 1.0 / math.Sqrt(qscale)
		e.m_Qhat[0] = qscale * w
		e.m_Qhat[1] = qscale * x
		e.m_Qhat[2] = qscale * y
		e.m_Qhat[3] = qscale * z
	}

	// Reset the alternative representation validity flags
	e.m_eulerValid = false
	e.m_fusedValid = false
}

// SetAttitudeEuler sets the current attitude estimate to a particular set of ZYX Euler angles.
func (e *AttitudeEstimator) SetAttitudeEuler(yaw, pitch, roll float64) {
	// halve the yaw, pitch and roll values (for calculation purposes only)
	yaw *= 0.5
	pitch *= 0.5
	roll *= 0.5

	// precalculate the required sin and cos values

	var (
		cpsi = math.Cos(yaw)
		spsi = math.Sin(yaw)
		cth  = math.Cos(pitch)
		sth  = math.Sin(pitch)
		cphi = math.Cos(roll)
		sphi = math.Sin(roll)
	)

	// calculate the required quaternion components

	var (
		w = cpsi*cth*cphi + spsi*sth*sphi
		x = cpsi*cth*sphi - spsi*sth*cphi
		y = cpsi*sth*cphi + spsi*cth*sphi
		z = spsi*cth*cphi - cpsi*sth*sphi
	)

	// Set the current attitude estimate to the calculated quaternion orientation
	e.SetAttitude(w, x, y, z)
}

func (e *AttitudeEstimator) SetAttitudeFused(yaw, pitch, roll float64, hemi bool) {
	// precalculate the sin values
	var (
		sth  = math.Sin(pitch)
		sphi = math.Sin(roll)
	)

	// Calculate the sine sum criterion
	crit := sth*sth + sphi*sphi

	// Calculate the tilt angle alpha
	var alpha float64
	if crit >= 1.0 {
		alpha = math.Pi * 2
	} else {
		if hemi {
			alpha = math.Acos(math.Sqrt(1 - crit))
		} else {
			alpha = math.Acos(-math.Sqrt(1 - crit))
		}
	}

	// Calculate the tilt axis gamma
	gamma := math.Atan2(sth, sphi)

	// Evaluate the required intermediate angles
	var (
		halpha  = 0.5 * alpha
		hpsi    = 0.5 * yaw
		hgampsi = gamma + hpsi
	)

	// Precalculate trigonometric terms involved in the quaternion expression
	var (
		chalpha  = math.Cos(halpha)
		shalpha  = math.Sin(halpha)
		chpsi    = math.Cos(hpsi)
		shpsi    = math.Sin(hpsi)
		chgampsi = math.Cos(hgampsi)
		shgampsi = math.Sin(hgampsi)
	)

	// Calculate the required quaternion orientation and set it as the current attitude estimate
	e.SetAttitude(chalpha*chpsi, shalpha*chgampsi, shalpha*shgampsi, chalpha*shpsi)
}

func (e *AttitudeEstimator) updateEuler() {
	// These calculations rely on the assumption that m_Qhat is a unit quaternion!
	//
	// The output ranges are:
	//   Yaw:    psi  = m_Ehat[0] is in (-pi,pi]
	//   Pitch: theta = m_Ehat[1] is in [-pi/2,pi/2]
	//   Roll:   phi  = m_Ehat[2] is in (-pi,pi]

	// Calculate pitch
	stheta := 2.0 * (e.m_Qhat[0]*e.m_Qhat[2] - e.m_Qhat[3]*e.m_Qhat[1])

	// Coerce stheta to [-1,1]
	if stheta >= 1 {
		stheta = 1
	} else if stheta <= -1 {
		stheta = -1
	}

	e.m_Ehat[1] = math.Asin(stheta)

	// Calculate yaw and roll
	ysq := e.m_Qhat[2] * e.m_Qhat[2]
	e.m_Ehat[0] = math.Atan2(e.m_Qhat[0]*e.m_Qhat[3]+e.m_Qhat[1]*e.m_Qhat[2], 0.5-(ysq+e.m_Qhat[3]*e.m_Qhat[3]))
	e.m_Ehat[2] = math.Atan2(e.m_Qhat[0]*e.m_Qhat[1]+e.m_Qhat[2]*e.m_Qhat[3], 0.5-(ysq+e.m_Qhat[1]*e.m_Qhat[1]))

	// Set the Euler angles valid flag
	e.m_eulerValid = true
}

func (e *AttitudeEstimator) updateFused() {
	// These calculations rely on the assumption that m_Qhat is a unit quaternion!
	//
	// The output ranges are:
	//   Fused yaw:    psi  = m_Fhat[0]  is in (-pi,pi]
	//   Fused pitch: theta = m_Fhat[1]  is in [-pi/2,pi/2]
	//   Fused roll:   phi  = m_Fhat[2]  is in [-pi/2,pi/2]
	//   Hemisphere:    h   = m_FhatHemi is in {-1,1} (stored as {false,true} respectively)

	// Calculate and wrap the fused yaw
	e.m_Fhat[0] = 2 * math.Atan2(e.m_Qhat[3], e.m_Qhat[0]) // Output of atan2 is [-pi,pi], so this expression is in [-2*pi,2*pi]
	if e.m_Fhat[0] > math.Pi {
		// Fhat[0] is now in [-2*pi,pi]
		e.m_Fhat[0] -= math.Pi * 2
	}
	if e.m_Fhat[0] <= -math.Pi {
		// Fhat[0] is now in (-pi,pi]
		e.m_Fhat[0] += math.Pi * 2
	}

	// Calculate the fused pitch and roll
	stheta := 2.0 * (e.m_Qhat[2]*e.m_Qhat[0] - e.m_Qhat[1]*e.m_Qhat[3])
	sphi := 2.0 * (e.m_Qhat[2]*e.m_Qhat[3] + e.m_Qhat[1]*e.m_Qhat[0])

	// Coerce stheta to [-1,1]
	if stheta >= 1 {
		stheta = 1
	} else if stheta <= -1 {
		stheta = -1
	}

	// Coerce sphi to [-1,1]
	if sphi >= 1 {
		sphi = 1
	} else if sphi <= -1 {
		sphi = -1
	}

	e.m_Fhat[1] = math.Asin(stheta)
	e.m_Fhat[2] = math.Asin(sphi)

	// Calculate the hemisphere of the rotation
	e.m_FhatHemi = (0.5-(e.m_Qhat[1]*e.m_Qhat[1]+e.m_Qhat[2]*e.m_Qhat[2]) >= 0)

	// Set the fused angles valid flag
	e.m_fusedValid = true
}
