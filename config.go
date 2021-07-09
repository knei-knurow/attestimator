package attestimator

func (e *AttitudeEstimator) Reset(quickLearn bool, resetGyroBias bool) {
	e.resetState(resetGyroBias)

	if quickLearn {
		e.SetLambda(0) // reset lambda
	} else {
		e.SetLambda(1)
	}
}

func (e *AttitudeEstimator) ResetAll(quickLearn bool) {
	// Initialise the acc-only resolution method
	e.SetAccMethod(MethodDefault)

	// Initialise the configuration variables
	e.SetPIGains(2.20, 2.65, 10.0, 1.25)
	e.SetQLTime(3)

	// Initialise the magnetometer calibration
	e.SetMagCalib(1, 0, 0)

	// Reset the attitude estimator
	e.Reset(true, true)
}

func (e *AttitudeEstimator) resetState(resetGyroBias bool) {
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

// SetAccMethod sets the acc-only measurement resolution method to use.
func (e *AttitudeEstimator) SetAccMethod(method AccMethod) {
	if method < MethodDefault || method >= MethodCount {
		e.accMethod = MethodDefault
	} else {
		e.accMethod = method
	}
}

func (e *AttitudeEstimator) PIGains() (Kp, Ti, KpQuick, TiQuick float64) {
	return e.m_Kp, e.m_Ti, e.m_KpQuick, e.m_TiQuick
}

// SetPIGains sets the PI gains to use in both normal situations and for quick learning.
//
// Each `Kp/Ti` pair must be a pair of positive values or the
// respective pair is not updated.
func (e *AttitudeEstimator) SetPIGains(Kp, Ti, KpQuick, TiQuick float64) {
	// Set the standard PI gains
	if Kp > 0 && Ti > 0 {
		e.m_Kp = Kp
		e.m_Ti = Ti
	}

	// Set the quick learning PI gains
	if KpQuick > 0 && TiQuick > 0 {
		e.m_KpQuick = KpQuick
		e.m_TiQuick = TiQuick
	}
}

// QLTime Returns the current quick learning time. A rate at which to
// auto-increment lambda is calculated so that quick learning fades over into
// standard operation in exactly `QLTime` seconds.
func (e *AttitudeEstimator) QLTime() float64 { return e.m_QLTime }

// SetQLTime sets the quick learning time to use.
//
// Non-positive values of QLtime are ignored by this function.
func (e *AttitudeEstimator) SetQLTime(QLtime float64) {
	if QLtime > 0 {
		e.m_QLTime = QLtime
	}
}
