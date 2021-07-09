package attestimator


func (e *AttitudeEstimator) SetGyroBias(bx, by, bz float64) {
	e.m_bhat[0] = bx
	e.m_bhat[1] = by
	e.m_bhat[2] = bz
}
