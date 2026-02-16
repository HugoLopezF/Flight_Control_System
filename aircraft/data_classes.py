from dataclasses import dataclass


@dataclass
class Geometry:
    S: float
    c: float
    b: float
    eps_s: float


@dataclass
class MassProperties:
    W: float
    I_xx: float
    I_yy: float
    I_zz: float
    I_xz: float


@dataclass
class ConditionCoefficients:
    CL_s: float
    CD_s: float
    CT_s: float
    Cm: float
    Cm_T: float


@dataclass
class FlightCondition:
    rho: float
    u_s: float
    theta_s: float
    h: float
    M: float
    cg: float
    q: float


@dataclass
class LongitudinalCoefficients:
    CL_0: float
    CL_alpha: float
    CL_q: float
    CL_alpha_dot: float
    CL_u: float
    CL_delta_e: float
    CD_0: float
    CD_alpha: float
    CD_u: float
    CD_delta_e: float
    Cm_0: float
    Cm_alpha: float
    Cm_q: float
    Cm_alpha_dot: float
    Cm_u: float
    Cm_Tu: float
    Cm_T_alpha: float
    Cm_delta_e: float
    CT_u: float


@dataclass
class LateralDirectionalCoefficients:
    CY_beta: float
    CY_p: float
    CY_r: float
    CY_delta_a: float
    CY_delta_r: float
    Cl_beta: float
    Cl_p: float
    Cl_r: float
    Cl_delta_a: float
    Cl_delta_r: float
    Cn_beta: float
    Cn_p: float
    Cn_r: float
    Cn_delta_a: float
    Cn_delta_r: float 


@dataclass
class StabilityCoefficients:
    long: LongitudinalCoefficients
    latdir: LateralDirectionalCoefficients


@dataclass
class LongitudinalStabilityDerivatives:
    Xu: float
    Xw: float
    Xdelta_e: float
    Zu: float
    Zw: float
    Zw_dot: float
    Zq: float
    Zdelta_e: float
    Mu: float
    Mw: float
    Mw_dot: float
    Mq: float
    Mdelta_e: float


@dataclass
class LateralDirectionalStabilityDerivatives:
    Yv: float
    Yp: float
    Yr: float
    Ydelta_r: float
    Lv: float
    Lp: float
    Lr: float
    Ldelta_a: float
    Ldelta_r: float
    Nv: float
    Np: float
    Nr: float
    Ndelta_a: float
    Ndelta_r: float


@dataclass(frozen=True)
class LinearizationParameters:
    geom: Geometry
    mass_prop: MassProperties
    flight_cond: FlightCondition
    cond_coeffs: ConditionCoefficients
    stab_coeffs: StabilityCoefficients


@dataclass
class StabilityDerivatives:
    long: LongitudinalStabilityDerivatives
    latdir: LateralDirectionalStabilityDerivatives