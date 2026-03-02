from dataclasses import dataclass


@dataclass
class Geometry:
    """
    Aircraft geometry.

    Attributes
    ----------
    S : float
        Wing surface area.
    c : float
        Wing chord.
    b : float
        Wingspan.
    eps_s : float
        Thrust pitch angle with respect to aircraft's x axis.
    """

    S: float
    c: float
    b: float
    eps_s: float


@dataclass
class MassProperties:
    """
    Aircraft mass properties.

    Attributes
    ----------
    W : float
        Aircraft's mass.
    I_xx : float
        Aircraft's moment of inertia with respect to x axis.
    I_yy : float
        Aircraft's moment of inertia with respect to y axis.
    I_zz : float
        Aircraft's moment of inertia with respect to z axis.
    I_xz : float
        Aircraft's moment of inertia with respect to x-z plane.
    """

    W: float
    I_xx: float
    I_yy: float
    I_zz: float
    I_xz: float


@dataclass
class ConditionCoefficients:
    """
    Aircraft flight condition coefficients.

    Attributes
    ----------
    CL_s : float
        Lift coefficient in stability axes.
    CD_s : float
        Drag coefficient in stability axes.
    CT_s : float
        Thrust coefficient in stability axes.
    Cm : float
        Aerodynamic moment coefficient in stability axes.
    Cm_T : float
        Thrust moment coefficient in stability axes.
    """

    CL_s: float
    CD_s: float
    CT_s: float
    Cm: float
    Cm_T: float


@dataclass
class FlightCondition:
    """
    Aircraft flight condition data.

    Attributes
    ----------
    rho : float
        Air density.
    u_s : float
        Aerodynamic speed.
    h : float
        Altitude.
    M : float
        Mach number.
    cg : float
        Gravity center position (% chord).
    q : float
        Dynamic pressure.
    """

    rho: float
    u_s: float
    theta_s: float
    h: float
    M: float
    cg: float
    q: float


@dataclass
class LongitudinalCoefficients:
    """
    Aircraft longitudinal stability coefficients data.

    Attributes
    ----------
    CL_0 : float
        Zero angle of attack lift coefficient.
    CL_alpha : float
        Derivative of lift coefficient with respect to angle of attack.
    CL_q : float
        Derivative of lift coefficient with respect to pitch rate.
    CL_alpha_dot : float
        Derivative of lift coefficient with respect to angle of attack rate.
    CL_u : float
        Derivative of lift coefficient with respect to x-axis aerodynamic speed.
    CL_delta_e : float
        Derivative of lift coefficient with respect to elevator deflection.
    CD_0 : float
        Zero angle of attack drag coefficient.
    CD_alpha : float
        Derivative of drag coefficient with respect to angle of attack.
    CD_u : float
        Derivative of drag coefficient with respect to x-axis aerodynamic speed.
    CD_delta_e : float
        Derivative of drag coefficient with respect to elevator deflection.
    Cm_0 : float
        Zero angle of attack pitch moment coefficient.
    Cm_alpha : float
        Derivative of pitch moment coefficient with respect to angle of attack.
    Cm_q : float
        Derivative of pitch moment coefficient with respect to pitch rate.
    Cm_alpha_dot : float
        Derivative of pitch moment coefficient with respect to angle of attack rate.
    Cm_u : float
        Derivative of pitch moment coefficient with respect to x-axis aerodynamic speed.
    Cm_Tu : float
        Derivative of thrust moment coefficient with respect to x-axis aerodynamic speed.
    Cm_T_alpha : float
        Derivative of thrust moment coefficient with respect to angle of attack.
    Cm_delta_e : float
        Derivative of pitch moment coefficient with respect to elevator deflection.
    CT_u : float
        Derivative of thrust coefficient with respect to x-axis aerodynamic speed.
    """

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
    """
    Aircraft lateral-directional stability coefficients data.

    Attributes
    ----------
    CY_beta : float
        Derivative of lateral force coefficient with respect to sideslip angle.
    CY_p : float
        Derivative of lateral force coefficient with respect to roll rate.
    CY_r : float
        Derivative of lateral force coefficient with respect to yaw rate.
    CY_delta_a : float
        Derivative of lateral force coefficient with respect to aileron deflection.
    CY_delta_r : float
        Derivative of lateral force coefficient with respect to rudder deflection.
    Cl_beta : float
        Derivative of roll moment coefficient with respect to sideslip angle.
    Cl_p : float
        Derivative of roll moment coefficient with respect to roll rate.
    Cl_r : float
        Derivative of roll moment coefficient with respect to yaw rate.
    Cl_delta_a : float
        Derivative of roll moment coefficient with respect to aileron deflection.
    Cl_delta_r : float
        Derivative of roll moment coefficient with respect to rudder deflection.
    Cn_beta : float
        Derivative of yaw moment coefficient with respect to sideslip angle.
    Cn_p : float
        Derivative of yaw moment coefficient with respect to roll rate.
    Cn_r : float
        Derivative of yaw moment coefficient with respect to yaw rate.
    Cn_delta_a : float
        Derivative of yaw moment coefficient with respect to aileron deflection.
    Cn_delta_r : float
        Derivative of yaw moment coefficient with respect to rudder deflection.
    """

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
    """
    Aircraft stability coefficients.

    Attributes
    ----------
    long : LongitudinalCoefficients
        Aircraft longitudinal stability coefficients.
    latdir : LateralDirectionalCoefficients
        Aircraft lateral-directional stability coefficients.
    """

    long: LongitudinalCoefficients
    latdir: LateralDirectionalCoefficients


@dataclass
class LongitudinalStabilityDerivatives:
    """
    Aircraft longitudinal stability derivatives.

    Attributes
    ----------
    Xu : float
        Derivative of longitudinal force with respect to x-axis aerodynamic speed.
    Xw : float
        Derivative of longitudinal force with respect to z-axis aerodynamic speed.
    Xdelta_e : float
        Derivative of longitudinal force with respect to elevator deflection.
    Zu : float
        Derivative of vertical force with respect to x-axis aerodynamic speed.
    Zw : float
        Derivative of vertical force with respect to z-axis aerodynamic speed.
    Zw_dot : float
        Derivative of vertical force with respect to z-axis aerodynamic acceleration.
    Zq : float
        Derivative of vertical force with respect to pitch rate.
    Mu : float
        Derivative of pitch moment with respect to x-axis aerodynamic speed.
    Mw : float
        Derivative of pitch moment with respect to z-axis aerodynamic speed.
    Mw_dot : float
        Derivative of pitch moment with respect to z-axis aerodynamic acceleration.
    Mq : float
        Derivative of pitch moment with respect to pitch rate.
    Mdelta_e : float
        Derivative of pitch moment with respect to elevator deflection.
    """

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
    """
    Aircraft lateral-directional stability derivatives.

    Attributes
    ----------
    Yv : float
        Derivative of lateral force with respect to y-axis aerodynamic speed.
    Yp : float
        Derivative of lateral force with respect to roll rate.
    Yr : float
        Derivative of lateral force with respect to yaw rate.
    Ydelta_r : float
        Derivative of lateral force with respect to rudder deflection.
    Lv : float
        Derivative of roll moment with respect to y-axis aerodynamic speed.
    Lp : float
        Derivative of roll moment with respect to roll rate.
    Lr : float
        Derivative of roll moment with respect to yaw rate.
    Ldelta_a : float
        Derivative of roll moment with respect to aileron deflection.
    Ldelta_r : float
        Derivative of roll moment with respect to rudder deflection.
    Nv : float
        Derivative of yaw moment with respect to y-axis aerodynamic speed.
    Np : float
        Derivative of yaw moment with respect to roll rate.
    Nr : float
        Derivative of yaw moment with respect to yaw rate.
    Ndelta_a : float
        Derivative of yaw moment with respect to aileron deflection.
    Ndelta_r : float
        Derivative of yaw moment with respect to rudder deflection.
    """

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
    """
    Aircraft linearized system parameters.

    Attributes
    ----------
    geom : Geometry
        Aircraft geometry.
    mass_prop : MassProperties
        Aircraft mass properties.
    flight_cond : FlightCondition
        Aircraft flight condition data.
    cond_coeffs : ConditionCoefficients
        Aircraft flight condition coefficients.
    stab_coeffs : StabilityCoefficients
        Aircraft stability coefficients.
    """

    geom: Geometry
    mass_prop: MassProperties
    flight_cond: FlightCondition
    cond_coeffs: ConditionCoefficients
    stab_coeffs: StabilityCoefficients


@dataclass
class StabilityDerivatives:
    """
    Aircraft stability derivatives.

    Attributes
    ----------
    long : Geometry
        Aircraft longitudinal stability derivatives.
    latdir : MassProperties
        Aircraft lateral-directional stability derivatives.
    """

    long: LongitudinalStabilityDerivatives
    latdir: LateralDirectionalStabilityDerivatives