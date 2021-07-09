package attestimator

// GyroBias returns a vector representing the current guroscope bias.
func (e *AttitudeEstimator) GyroBias() (bx, by, bz float64) {
	return e.m_bhat[0], e.m_bhat[1], e.m_bhat[2]
}

// SetGyroBias sets the current estimated gyroscope bias to a particular vector value.
func (e *AttitudeEstimator) SetGyroBias(bx, by, bz float64) {
	e.m_bhat[0] = bx
	e.m_bhat[1] = by
	e.m_bhat[2] = bz
}
