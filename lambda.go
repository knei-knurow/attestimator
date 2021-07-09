package attestimator

// Lambda returns the current lambda paramater.
func (e *Estimator) Lambda() float64 {
	return e.m_lambda
}

// SetLambda sets the estimator's lambda.
// The value of 1 turns quick learning off.
// Values outside of the [0,1] range are coerced.
func (e *Estimator) SetLambda(value float64) {
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
