from .data_classes import LongitudinalStabilityDerivatives, LateralDirectionalStabilityDerivatives, LinearizationParameters
from utilities.constants import g
from math import sin, cos


class StabilityDerivativesCalculator:
        
    @staticmethod
    def calculate_longitudinal(params: LinearizationParameters) -> LongitudinalStabilityDerivatives:
        """
        Calculate all longitudinal stability derivatives.

        """
        return LongitudinalStabilityDerivatives(
            Xu = StabilityDerivativesCalculator.calculate_Xu(params),
            Xw = StabilityDerivativesCalculator.calculate_Xw(params),
            Xdelta_e = StabilityDerivativesCalculator.calculate_Xdelta_e(params),
            Zu = StabilityDerivativesCalculator.calculate_Zu(params),
            Zw = StabilityDerivativesCalculator.calculate_Zw(params),
            Zw_dot = StabilityDerivativesCalculator.calculate_Zw_dot(params),
            Zq = StabilityDerivativesCalculator.calculate_Zq(params),
            Zdelta_e = StabilityDerivativesCalculator.calculate_Zdelta_e(params),
            Mu = StabilityDerivativesCalculator.calculate_Mu(params),
            Mw = StabilityDerivativesCalculator.calculate_Mw(params),
            Mw_dot = StabilityDerivativesCalculator.calculate_Mw_dot(params),
            Mq = StabilityDerivativesCalculator.calculate_Mq(params),
            Mdelta_e = StabilityDerivativesCalculator.calculate_Mdelta_e(params)
        )
    
    @staticmethod
    def calculate_lateral_directional(params: LinearizationParameters) -> LongitudinalStabilityDerivatives:
        """
        Calculate all lateral directional stability derivatives.

        """
        return LateralDirectionalStabilityDerivatives(
            Yv = StabilityDerivativesCalculator.calculate_Yv(params),
            Yp = StabilityDerivativesCalculator.calculate_Yp(params),
            Yr = StabilityDerivativesCalculator.calculate_Yr(params),
            Ydelta_r = StabilityDerivativesCalculator.calculate_Ydelta_r(params),
            Lv = StabilityDerivativesCalculator.calculate_Lv(params),
            Lp = StabilityDerivativesCalculator.calculate_Lp(params),
            Lr = StabilityDerivativesCalculator.calculate_Lr(params),
            Ldelta_a = StabilityDerivativesCalculator.calculate_Ldelta_a(params),
            Ldelta_r = StabilityDerivativesCalculator.calculate_Ldelta_r(params),
            Nv = StabilityDerivativesCalculator.calculate_Nv(params),
            Np = StabilityDerivativesCalculator.calculate_Np(params),
            Nr = StabilityDerivativesCalculator.calculate_Nr(params),
            Ndelta_a = StabilityDerivativesCalculator.calculate_Ndelta_a(params),
            Ndelta_r = StabilityDerivativesCalculator.calculate_Ndelta_r(params)
        )
    
    # Longitudinal derivatives
    def calculate_Xu(params) -> float:
        """
        Calculate Xu derivative.

        :rtype: float
        """
        W = params.mass_prop.W
        theta_s = params.flight_cond.theta_s
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S
        CT_u = params.stab_coeffs.long.CT_u
        CD_u = params.stab_coeffs.long.CD_u
        eps_s = params.geom.eps_s

        CX_s = W * g * sin(theta_s) / (0.5 * rho * S * u_s ** 2)
        CX_u = CT_u * cos(eps_s) - CD_u

        Xu = rho * S * u_s * (CX_s + 0.5 * CX_u)
        return Xu
    
    def calculate_Xw(params) -> float:
        """
        Calculate Xw derivative.

        :rtype: float
        """
        CL_s = params.cond_coeffs.CL_s
        CD_alpha = params.stab_coeffs.long.CD_alpha
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S

        CX_alpha = CL_s - CD_alpha

        Xw = 0.5 * rho * S * u_s * CX_alpha
        return Xw
    
    def calculate_Xdelta_e(params) -> float:
        """
        Calculate Xdelta_e derivative.

        :rtype: float
        """
        CD_delta_e = params.stab_coeffs.long.CD_delta_e
        CX_delta_e = - CD_delta_e
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S

        Xdelta_e = 0.5 * rho * S * u_s ** 2 * CX_delta_e
        return Xdelta_e
    
    def calculate_Zu(params) -> float:
        """
        Calculate Zu derivative.

        :rtype: float
        """
        W = params.mass_prop.W
        theta_s = params.flight_cond.theta_s
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S
        CT_u = params.stab_coeffs.long.CT_u
        CL_u = params.stab_coeffs.long.CL_u
        eps_s = params.geom.eps_s

        CZ_s = - W * g * cos(theta_s) / (0.5 * rho * S * u_s ** 2)
        CZ_u = - CT_u * sin(eps_s) - CL_u

        Zu = rho * S * u_s * (CZ_s + 0.5 * CZ_u)
        return Zu
    
    def calculate_Zw(params) -> float:
        """
        Calculate Zw derivative.

        :rtype: float
        """
        CL_alpha = params.stab_coeffs.long.CL_alpha
        CD_s =  params.cond_coeffs.CD_s
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S

        CZ_alpha = - CL_alpha - CD_s

        Zw = 0.5 * rho * S * u_s * CZ_alpha
        return Zw
    
    def calculate_Zw_dot(params) -> float:
        """
        Calculate Zw_dot derivative.

        :rtype: float
        """
        CL_alpha_dot = params.stab_coeffs.long.CL_alpha_dot
        rho = params.flight_cond.rho
        S = params.geom.S
        c = params.geom.c

        CZ_alpha_dot = - CL_alpha_dot

        Zw_dot = 0.25 * rho * S * c * CZ_alpha_dot
        return Zw_dot
    
    def calculate_Zq(params) -> float:
        """
        Calculate Zq derivative.

        :rtype: float
        """
        CL_q = params.stab_coeffs.long.CL_q
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S
        c = params.geom.c

        CZ_q = - CL_q

        Zq = 0.25 * rho * S * c * u_s * CZ_q
        return Zq
    
    def calculate_Zdelta_e(params) -> float:
        """
        Calculate Zdelta_e derivative.

        :rtype: float
        """
        CL_delta_e = params.stab_coeffs.long.CL_delta_e
        CZ_delta_e = - CL_delta_e
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S

        Zdelta_e = 0.5 * rho * S * u_s ** 2 * CZ_delta_e
        return Zdelta_e
    
    def calculate_Mu(params) -> float:
        """
        Calculate Mu derivative.

        :rtype: float
        """
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S
        c = params.geom.c
        Cm_u = params.stab_coeffs.long.Cm_u

        Mu = 0.5* rho * S * c * u_s * Cm_u
        return Mu
    
    def calculate_Mw(params) -> float:
        """
        Calculate Mw derivative.

        :rtype: float
        """
        Cm_alpha = params.stab_coeffs.long.Cm_alpha
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S
        c = params.geom.c

        Mw = 0.5 * rho * S * c * u_s * Cm_alpha
        return Mw
    
    def calculate_Mw_dot(params) -> float:
        """
        Calculate Mw_dot derivative.

        :rtype: float
        """
        Cm_alpha_dot = params.stab_coeffs.long.Cm_alpha_dot
        rho = params.flight_cond.rho
        S = params.geom.S
        c = params.geom.c

        Mw_dot = 0.25 * rho * S * c ** 2 * Cm_alpha_dot
        return Mw_dot
    
    def calculate_Mq(params) -> float:
        """
        Calculate Mq derivative.

        :rtype: float
        """
        Cm_q = params.stab_coeffs.long.Cm_q
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S
        c = params.geom.c

        Mq = 0.25 * rho * S * c ** 2 * u_s * Cm_q
        return Mq
    
    def calculate_Mdelta_e(params) -> float:
        """
        Calculate Mdelta_e derivative.

        :rtype: float
        """
        Cm_delta_e = params.stab_coeffs.long.Cm_delta_e
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S
        c = params.geom.c

        Mdelta_e = 0.5 * rho * S * c * u_s ** 2 * Cm_delta_e
        return Mdelta_e
    
    # Lateral-directional derivatives
    def calculate_Yv(params) -> float:
        """
        Calculate Yv derivative.

        :rtype: float
        """
        CY_beta = params.stab_coeffs.latdir.CY_beta
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S

        Yv = 0.5 * rho * S * u_s * CY_beta
        return Yv
    
    def calculate_Yp(params) -> float:
        """
        Calculate Yp derivative.

        :rtype: float
        """
        CY_p = params.stab_coeffs.latdir.CY_p
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S
        b = params.geom.b

        Yp = 0.25 * rho * S * b * u_s * CY_p
        return Yp
    
    def calculate_Yr(params) -> float:
        """
        Calculate Yr derivative.

        :rtype: float
        """
        CY_r = params.stab_coeffs.latdir.CY_r
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S
        b = params.geom.b

        Yr = 0.25 * rho * S * b * u_s * CY_r
        return Yr
    
    def calculate_Ydelta_r(params) -> float:
        """
        Calculate Y_delta_r derivative.

        :rtype: float
        """
        CY_delta_r = params.stab_coeffs.latdir.CY_delta_r
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S

        Ydelta_r = 0.5 * rho * S * u_s ** 2 * CY_delta_r
        return Ydelta_r
    
    def calculate_Lv(params) -> float:
        """
        Calculate Lv derivative.

        :rtype: float
        """
        Cl_beta = params.stab_coeffs.latdir.Cl_beta
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S
        b = params.geom.b

        Lv = 0.5 * rho * S * b * u_s * Cl_beta
        return Lv
    
    def calculate_Lp(params) -> float:
        """
        Calculate Lp derivative.

        :rtype: float
        """
        Cl_p = params.stab_coeffs.latdir.Cl_p
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S
        b = params.geom.b

        Lp = 0.25 * rho * S * b ** 2 * u_s * Cl_p
        return Lp
    
    def calculate_Lr(params) -> float:
        """
        Calculate Lr derivative.

        :rtype: float
        """
        Cl_r = params.stab_coeffs.latdir.Cl_r
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S
        b = params.geom.b

        Lr = 0.25 * rho * S * b ** 2 * u_s * Cl_r
        return Lr
    
    def calculate_Ldelta_a(params) -> float:
        """
        Calculate Ldelta_a derivative.

        :rtype: float
        """
        Cl_delta_a = params.stab_coeffs.latdir.Cl_delta_a
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S
        b = params.geom.b

        Ldelta_a = 0.5 * rho * S * b * u_s ** 2 * Cl_delta_a
        return Ldelta_a
    
    def calculate_Ldelta_r(params) -> float:
        """
        Calculate Ldelta_r derivative.

        :rtype: float
        """
        Cl_delta_r = params.stab_coeffs.latdir.Cl_delta_r
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S
        b = params.geom.b

        Ldelta_r = 0.5 * rho * S * b * u_s ** 2 * Cl_delta_r
        return Ldelta_r
    
    def calculate_Nv(params) -> float:
        """
        Calculate Nv derivative.

        :rtype: float
        """
        Cn_beta = params.stab_coeffs.latdir.Cn_beta
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S
        b = params.geom.b

        Nv = 0.5 * rho * S * b * u_s * Cn_beta
        return Nv
    
    def calculate_Np(params) -> float:
        """
        Calculate Np derivative.

        :rtype: float
        """
        Cn_p = params.stab_coeffs.latdir.Cn_p
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S
        b = params.geom.b

        Np = 0.25 * rho * S * b ** 2 * u_s * Cn_p
        return Np
    
    def calculate_Nr(params) -> float:
        """
        Calculate Nr derivative.

        :rtype: float
        """
        Cn_r = params.stab_coeffs.latdir.Cn_r
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S
        b = params.geom.b

        Nr = 0.25 * rho * S * b ** 2 * u_s * Cn_r
        return Nr
    
    def calculate_Ndelta_a(params) -> float:
        """
        Calculate Ndelta_a derivative.

        :rtype: float
        """
        Cn_delta_a = params.stab_coeffs.latdir.Cn_delta_a
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S
        b = params.geom.b

        Ndelta_a = 0.5 * rho * S * b * u_s ** 2 * Cn_delta_a
        return Ndelta_a
    
    def calculate_Ndelta_r(params) -> float:
        """
        Calculate Ndelta_r derivative.

        :rtype: float
        """
        Cn_delta_r = params.stab_coeffs.latdir.Cn_delta_r
        rho = params.flight_cond.rho
        u_s = params.flight_cond.u_s
        S = params.geom.S
        b = params.geom.b

        Ndelta_r = 0.5 * rho * S * b * u_s ** 2 * Cn_delta_r
        return Ndelta_r