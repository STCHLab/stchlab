import numpy as np

from .CEF import CEF, SIGMA_MO2


class CeO2(CEF):
    name = "CeO2"
    sigma_model = SIGMA_MO2
    delta_max = 0.25
    coef = np.array([
        -0.34172328711020095,
        -32.393301703343276,
        89.77868949465461,
        50.17315692774026,
        -101.99944927577718,
        309.4130186362184,
        0,
        0,
        0,
    ])


name = CeO2.name
coef = CeO2.coef
delta_max = CeO2.delta_max
Cp = CeO2.Cp
G_R = CeO2.G_R
H_R = CeO2.H_R
