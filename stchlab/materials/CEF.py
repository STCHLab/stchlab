"""Compound Energy Formalism (CEF) helpers for metal oxide models."""

import math

import numpy as np


SIGMA_MO2 = "MO2"
SIGMA_PEROVSKITE_B = "AA'B(1-y)B'yO3"

DEFAULT_CP = 80
COEF_COUNT = 9
R_GAS = 8.3144598
T_REF = 1073.15
P_REF = 21223.5
LOG_1E5_OVER_P_REF = math.log(1e5 / P_REF)

MU_O2_C0 = 23.957412864121114
MU_O2_C1 = 2.7594606723207535
MU_O2_C2 = 0.05818937063497993
MU_O2_CLOG = -4.912356983241369
MU_O2_CBETA_LOG = -0.7950813082655013

SIGMA_ALIASES = {
    "mo2": SIGMA_MO2,
    "perovskite": SIGMA_PEROVSKITE_B,
    "aab1ybyo3": SIGMA_PEROVSKITE_B,
}


def _normalize_sigma_model(sigma_model):
    if sigma_model in {SIGMA_MO2, SIGMA_PEROVSKITE_B}:
        return sigma_model

    key = "".join(ch for ch in str(sigma_model).lower() if ch.isalnum())
    try:
        return SIGMA_ALIASES[key]
    except KeyError:
        raise ValueError(f"Unknown Compound Energy Formalism sigma model: {sigma_model!r}") from None


def _require_perovskite_y(y):
    if y is None:
        raise ValueError("Perovskite CEF materials must define y.")
    return y


def _default_delta_max(y, sigma_model):
    if sigma_model == SIGMA_MO2:
        return 0.25
    y = _require_perovskite_y(y)
    return min((1 - y) / 2, 1 / 3)


def _as_coef9(coef):
    coef = np.asarray(coef, dtype=float)
    if coef.shape != (COEF_COUNT,):
        raise ValueError(f"CEF coef must contain exactly {COEF_COUNT} coefficients.")
    return coef


def _has_p2(coef6, coef7, coef8):
    return not (coef6 == 0.0 and coef7 == 0.0 and coef8 == 0.0)


def _sigma_mo2(delta, log):
    # MO2 fluorite model, used by CeO2.
    return 2 * log(1 - 2 * delta) + log(2 - delta) - 3 * log(delta)


def _sigma_perovskite_b(delta, y, log):
    # Perovskite AA'B(1-y)B'yO3 model, where B cation changes oxidation state.
    return (
        2 * log(1 - y - 2 * delta)
        - log(delta / 3)
        - 2 * log(2 * delta)
        + log(1 - delta / 3)
    )


def _mu_o2_bar(beta, log_beta):
    return (
        MU_O2_C0
        + MU_O2_C1 * beta
        + MU_O2_C2 * beta * beta
        + MU_O2_CLOG * log_beta
        + MU_O2_CBETA_LOG * beta * log_beta
    )


def _dmu_o2_bar_d_beta(beta, log_beta):
    return (
        MU_O2_C1
        + 2 * MU_O2_C2 * beta
        + MU_O2_CLOG / beta
        + MU_O2_CBETA_LOG * (log_beta + 1)
    )


def make_cef_functions(coef, y=None, cp=DEFAULT_CP, sigma_model=SIGMA_PEROVSKITE_B, delta_max=None):
    """Build fast module-level functions for a Compound Energy Formalism material."""
    sigma_model = _normalize_sigma_model(sigma_model)
    coef0, coef1, coef2, coef3, coef4, coef5, coef6, coef7, coef8 = _as_coef9(coef)

    if sigma_model == SIGMA_PEROVSKITE_B:
        y = _require_perovskite_y(y)
    if delta_max is None:
        delta_max = _default_delta_max(y, sigma_model)

    is_mo2 = sigma_model == SIGMA_MO2
    log = math.log
    has_p2 = _has_p2(coef6, coef7, coef8)
    coef4_half = coef4 / 2
    coef5_third = coef5 / 3
    coef7_half = coef7 / 2
    coef8_third = coef8 / 3

    def Cp(T):
        # J/mol K
        return cp

    def G_R(delta, T, pO2):
        beta = T_REF / T
        delta2 = delta * delta
        log_beta = log(beta)

        p0 = coef0 + coef1 * delta + coef2 * delta2
        p1 = coef3 + coef4 * delta + coef5 * delta2
        sigma = _sigma_mo2(delta, log) if is_mo2 else _sigma_perovskite_b(delta, y, log)

        if has_p2:
            p2 = coef6 + coef7 * delta + coef8 * delta2
            dmu_sol = p0 - beta * p1 - log_beta * p2 + sigma
        else:
            dmu_sol = p0 - beta * p1 + sigma

        muO2ref = _mu_o2_bar(beta, log_beta) + LOG_1E5_OVER_P_REF
        u = 0.5 * log(P_REF / pO2)
        equilibriumR = dmu_sol - (-0.5 * muO2ref - u)
        return -R_GAS * T * equilibriumR

    def H_R(delta, T):
        beta = T_REF / T
        delta2 = delta * delta
        delta3 = delta2 * delta
        log_beta = log(beta)

        ip1 = coef3 * delta + coef4_half * delta2 + coef5_third * delta3
        dmuO2bar = _dmu_o2_bar_d_beta(beta, log_beta)

        if has_p2:
            ip2 = coef6 * delta + coef7_half * delta2 + coef8_third * delta3
            return -R_GAS * T_REF * (-ip1 - ip2 / beta + 0.5 * dmuO2bar * delta)
        return -R_GAS * T_REF * (-ip1 + 0.5 * dmuO2bar * delta)

    return delta_max, Cp, G_R, H_R


class CEF:
    """Base class for Compound Energy Formalism (CEF) thermodynamic models."""

    name = None
    coef = None
    y = None
    cp = DEFAULT_CP
    delta_max = None
    sigma_model = SIGMA_PEROVSKITE_B

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cls.sigma_model = _normalize_sigma_model(cls.sigma_model)

        if cls.coef is None:
            return

        cls.coef = _as_coef9(cls.coef)
        cls.delta_max, Cp, G_R, H_R = make_cef_functions(
            cls.coef,
            cls.y,
            cp=cls.cp,
            sigma_model=cls.sigma_model,
            delta_max=cls.delta_max,
        )
        cls.Cp = staticmethod(Cp)
        cls.G_R = staticmethod(G_R)
        cls.H_R = staticmethod(H_R)
