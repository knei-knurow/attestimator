package attestimator

// Reset resets the entire class, except for the magnetometer calibration,
// the acc-only resolution method, and the configuration variables
// (i.e. the PI gains and the quick learn time).
func (e *Estimator) Reset(quickLearn bool, resetGyroBias bool) {
	e.ResetState(resetGyroBias)

	if quickLearn {
		e.SetLambda(0) // reset lambda
	} else {
		e.SetLambda(1)
	}
}

// ResetAll resets the entire class, including all variables not reset by the `reset()` function.
func (e *Estimator) ResetAll(quickLearn bool) {
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

// ResetState is equivalent to Reset(), but also leaves the lambda value untouched.
func (e *Estimator) ResetState(resetGyroBias bool) {
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

// AccMethod returns the currently selected acc-only measurement resolution method.
func (e *Estimator) AccMethod() AccMethod {
	return e.accMethod
}

// SetAccMethod sets the acc-only measurement resolution method to use.
func (e *Estimator) SetAccMethod(method AccMethod) {
	if method < MethodDefault || method >= MethodCount {
		e.accMethod = MethodDefault
	} else {
		e.accMethod = method
	}
}

// PIGains returns the current attitude estimator PI gains. Refer to `getLambda()` for
// more information on the PI gains.
func (e *Estimator) PIGains() (Kp, Ti, KpQuick, TiQuick float64) {
	return e.m_Kp, e.m_Ti, e.m_KpQuick, e.m_TiQuick
}

// SetPIGains sets the PI gains to use in both normal situations and for quick learning.
//
// Each `Kp/Ti` pair must be a pair of positive values or the
// respective pair is not updated.
func (e *Estimator) SetPIGains(Kp, Ti, KpQuick, TiQuick float64) {
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
func (e *Estimator) QLTime() float64 {
	return e.m_QLTime
}

// SetQLTime sets the quick learning time to use.
//
// Non-positive values of QLtime are ignored by this function.
func (e *Estimator) SetQLTime(QLtime float64) {
	if QLtime > 0 {
		e.m_QLTime = QLtime
	}
}
