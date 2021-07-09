package attestimator

// SetMagCalib sets the magnetometer calibration vector to use in the update functions.
//
// This should be the value of `(magX,magY,magZ)` that corresponds to a true
// orientation of identity.
func (e *AttitudeEstimator) SetMagCalib(mtx, mty, mtz float64) {
	e.m_magCalib[0] = mtx
	e.m_magCalib[1] = mty
	e.m_magCalib[2] = mtz
}

// MagCalib returns the current magnetometer calibration.
func (e *AttitudeEstimator) MagCalib() [3]float64 {
	return e.m_magCalib
}
