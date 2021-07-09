package attestimator

// EulerYaw returns the ZYX Euler yaw of the current attitude estimate (1st of the three ZYX Euler angles).
func (e *AttitudeEstimator) EulerYaw() float64 {
	if !e.m_eulerValid {
		e.updateEuler()
	}
	return e.m_Ehat[0]
}

// EulerPitch returns the ZYX Euler pitch of the current attitude estimate (2nd of the three ZYX Euler angles).
func (e *AttitudeEstimator) EulerPitch() float64 {
	if !e.m_eulerValid {
		e.updateEuler()
	}
	return e.m_Ehat[1]
}

// EulerRoll returns the ZYX Euler roll of the current attitude estimate (3rd of the three ZYX Euler angles).
func (e *AttitudeEstimator) EulerRoll() float64 {
	if !e.m_eulerValid {
		e.updateEuler()
	}
	return e.m_Ehat[2]
}

// FusedYaw returns the fused yaw of the current attitude estimate (1st of the fused angles).
func (e *AttitudeEstimator) FusedYaw() float64 {
	if !e.m_fusedValid {
		e.updateFused()
	}
	return e.m_Fhat[0]
}

// FusedPitch returns the fused pitch of the current attitude estimate (2nd of the fused angles).
func (e *AttitudeEstimator) FusedPitch() float64 {
	if !e.m_fusedValid {
		e.updateFused()
	}
	return e.m_Fhat[1]
}

// FusedRoll returns the fused roll of the current attitude estimate (3rd of the fused angles).
func (e *AttitudeEstimator) FusedRoll() float64 {
	if !e.m_fusedValid {
		e.updateFused()
	}
	return e.m_Fhat[2]
}
