package attestimator

func (e *AttitudeEstimator) Lambda() float64 {
	return e.m_lambda
}

// SetLambda sets estimator's lambda to the desired value.
// The value of 1 turns quick learning off. Values outside of `[0,1]` are coerced.
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

// ResetLambda sets lambda to 0.
func (e *AttitudeEstimator) ResetLambda() {
	e.SetLambda(0)
}
