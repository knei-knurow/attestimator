package attestimator

import "math"

// GetAttitude gets the current attitude estimate.
func (e *Estimator) GetAttitude() (w, x, y, z float64) {
	return e.m_Qhat[0], e.m_Qhat[1], e.m_Qhat[2], e.m_Qhat[3]
}

// SetAttitude sets the current attitude estimate.
func (e *Estimator) SetAttitude(w, x, y, z float64) {
	qscale := w*w + x*x + y*y + z*z // calculate the quaternion square norm

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
func (e *Estimator) SetAttitudeEuler(yaw, pitch, roll float64) {
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

func (e *Estimator) SetAttitudeFused(yaw, pitch, roll float64, hemi bool) {
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

func (e *Estimator) updateEuler() {
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

func (e *Estimator) updateFused() {
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
