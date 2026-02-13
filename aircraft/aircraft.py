# -*- coding: utf-8 -*-
"""
Aircraft
=========

Aircraft parameters and calculations.

"""

from pathlib import Path
import os
import json
from math import sin, cos
from utilities.constants import g 
from utilities.unit_conversion import LengthUnit, WeightUnit, AngleUnit, convert_length, convert_weight, convert_angle
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
    Z_q: float
    Zdelta_e: float
    Mu: float
    Mw: float
    Mw_dot: float
    M_q: float
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


# @dataclass
# class StabilityDerivatives:
#     long: LongitudinalStabilityDerivatives
#     latdir: LateralDirectionalStabilityDerivatives


class StabilityDerivatives:
    def __init__(self, Aircraft):
        self.aircraft = Aircraft
        
    def calculate_all(self) -> None:
        """
        Calculate all stability derivatives.

        """
        for function in dir(self):
            if function.startswith('calculate_') and not function.endswith('all'):
                value = getattr(self, function)()
                setattr(self, function.replace('calculate_', ''), value)
    
    # Longitudinal derivatives
    def calculate_Xu(self) -> float:
        """
        Calculate Xu derivative.

        :rtype: float
        """
        W = self.aircraft.mass_prop.W
        theta_s = self.aircraft.flight_cond.theta_s
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S
        CT_u = self.aircraft.stab_coeffs.long.CT_u
        CD_u = self.aircraft.stab_coeffs.long.CD_u
        eps_s = self.aircraft.geom.eps_s

        CX_s = W * g * sin(theta_s) / (0.5 * rho * S * u_s ** 2)
        CX_u = CT_u * cos(eps_s) - CD_u

        Xu = rho * S * u_s * (CX_s + 0.5 * CX_u)
        return Xu
    
    def calculate_Xw(self) -> float:
        """
        Calculate Xw derivative.

        :rtype: float
        """
        CL_s = self.aircraft.cond_coeffs.CL_s
        CD_alpha = self.aircraft.stab_coeffs.long.CD_alpha
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S

        CX_alpha = CL_s - CD_alpha

        Xw = 0.5 * rho * S * u_s * CX_alpha
        return Xw
    
    def calculate_Xdelta_e(self) -> float:
        """
        Calculate Xdelta_e derivative.

        :rtype: float
        """
        CD_delta_e = self.aircraft.stab_coeffs.long.CD_delta_e
        CX_delta_e = - CD_delta_e
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S

        Xdelta_e = 0.5 * rho * S * u_s ** 2 * CX_delta_e
        return Xdelta_e
    
    def calculate_Zu(self) -> float:
        """
        Calculate Zu derivative.

        :rtype: float
        """
        W = self.aircraft.mass_prop.W
        theta_s = self.aircraft.flight_cond.theta_s
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S
        CT_u = self.aircraft.stab_coeffs.long.CT_u
        CL_u = self.aircraft.stab_coeffs.long.CL_u
        eps_s = self.aircraft.geom.eps_s

        CZ_s = - W * g * cos(theta_s) / (0.5 * rho * S * u_s ** 2)
        CZ_u = - CT_u * sin(eps_s) - CL_u

        Zu = rho * S * u_s * (CZ_s + 0.5 * CZ_u)
        return Zu
    
    def calculate_Zw(self) -> float:
        """
        Calculate Zw derivative.

        :rtype: float
        """
        CL_alpha = self.aircraft.stab_coeffs.long.CL_alpha
        CD_s =  self.aircraft.cond_coeffs.CD_s
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S

        CZ_alpha = - CL_alpha - CD_s

        Zw = 0.5 * rho * S * u_s * CZ_alpha
        return Zw
    
    def calculate_Zw_dot(self) -> float:
        """
        Calculate Zw_dot derivative.

        :rtype: float
        """
        CL_alpha_dot = self.aircraft.stab_coeffs.long.CL_alpha_dot
        rho = self.aircraft.flight_cond.rho
        S = self.aircraft.geom.S
        c = self.aircraft.geom.c

        CZ_alpha_dot = - CL_alpha_dot

        Zw_dot = 0.25 * rho * S * c * CZ_alpha_dot
        return Zw_dot
    
    def calculate_Zq(self) -> float:
        """
        Calculate Zq derivative.

        :rtype: float
        """
        CL_q = self.aircraft.stab_coeffs.long.CL_q
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S
        c = self.aircraft.geom.c

        CZ_q = - CL_q

        Zq = 0.25 * rho * S * c * u_s * CZ_q
        return Zq
    
    def calculate_Zdelta_e(self) -> float:
        """
        Calculate Zdelta_e derivative.

        :rtype: float
        """
        CL_delta_e = self.aircraft.stab_coeffs.long.CL_delta_e
        CZ_delta_e = - CL_delta_e
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S

        Zdelta_e = 0.5 * rho * S * u_s ** 2 * CZ_delta_e
        return Zdelta_e
    
    def calculate_Mu(self) -> float:
        """
        Calculate Mu derivative.

        :rtype: float
        """
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S
        c = self.aircraft.geom.c
        Cm_u = self.aircraft.stab_coeffs.long.Cm_u

        Mu = 0.5* rho * S * c * u_s * Cm_u
        return Mu
    
    def calculate_Mw(self) -> float:
        """
        Calculate Mw derivative.

        :rtype: float
        """
        Cm_alpha = self.aircraft.stab_coeffs.long.Cm_alpha
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S
        c = self.aircraft.geom.c

        Mw = 0.5 * rho * S * c * u_s * Cm_alpha
        return Mw
    
    def calculate_Mw_dot(self) -> float:
        """
        Calculate Mw_dot derivative.

        :rtype: float
        """
        Cm_alpha_dot = self.aircraft.stab_coeffs.long.Cm_alpha_dot
        rho = self.aircraft.flight_cond.rho
        S = self.aircraft.geom.S
        c = self.aircraft.geom.c

        Mw_dot = 0.25 * rho * S * c ** 2 * Cm_alpha_dot
        return Mw_dot
    
    def calculate_Mq(self) -> float:
        """
        Calculate Mq derivative.

        :rtype: float
        """
        Cm_q = self.aircraft.stab_coeffs.long.Cm_q
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S
        c = self.aircraft.geom.c

        Mq = 0.25 * rho * S * c ** 2 * u_s * Cm_q
        return Mq
    
    def calculate_Mdelta_e(self) -> float:
        """
        Calculate Mdelta_e derivative.

        :rtype: float
        """
        Cm_delta_e = self.aircraft.stab_coeffs.long.Cm_delta_e
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S
        c = self.aircraft.geom.c

        Mdelta_e = 0.5 * rho * S * c * u_s ** 2 * Cm_delta_e
        return Mdelta_e
    
    # Lateral-directional derivatives
    def calculate_Yv(self) -> float:
        """
        Calculate Yv derivative.

        :rtype: float
        """
        CY_beta = self.aircraft.stab_coeffs.latdir.CY_beta
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S

        Yv = 0.5 * rho * S * u_s * CY_beta
        return Yv
    
    def calculate_Yp(self) -> float:
        """
        Calculate Yp derivative.

        :rtype: float
        """
        CY_p = self.aircraft.stab_coeffs.latdir.CY_p
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S
        b = self.aircraft.geom.b

        Yp = 0.25 * rho * S * b * u_s * CY_p
        return Yp
    
    def calculate_Yr(self) -> float:
        """
        Calculate Yr derivative.

        :rtype: float
        """
        CY_r = self.aircraft.stab_coeffs.latdir.CY_r
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S
        b = self.aircraft.geom.b

        Yr = 0.25 * rho * S * b * u_s * CY_r
        return Yr
    
    def calculate_Ydelta_r(self) -> float:
        """
        Calculate Y_delta_r derivative.

        :rtype: float
        """
        CY_delta_r = self.aircraft.stab_coeffs.latdir.CY_delta_r
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S

        Ydelta_r = 0.5 * rho * S * u_s ** 2 * CY_delta_r
        return Ydelta_r
    
    def calculate_Lv(self) -> float:
        """
        Calculate Lv derivative.

        :rtype: float
        """
        Cl_beta = self.aircraft.stab_coeffs.latdir.Cl_beta
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S
        b = self.aircraft.geom.b

        Lv = 0.5 * rho * S * b * u_s * Cl_beta
        return Lv
    
    def calculate_Lp(self) -> float:
        """
        Calculate Lp derivative.

        :rtype: float
        """
        Cl_p = self.aircraft.stab_coeffs.latdir.Cl_p
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S
        b = self.aircraft.geom.b

        Lp = 0.25 * rho * S * b ** 2 * u_s * Cl_p
        return Lp
    
    def calculate_Lr(self) -> float:
        """
        Calculate Lr derivative.

        :rtype: float
        """
        Cl_r = self.aircraft.stab_coeffs.latdir.Cl_r
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S
        b = self.aircraft.geom.b

        Lr = 0.25 * rho * S * b ** 2 * u_s * Cl_r
        return Lr
    
    def calculate_Ldelta_a(self) -> float:
        """
        Calculate Ldelta_a derivative.

        :rtype: float
        """
        Cl_delta_a = self.aircraft.stab_coeffs.latdir.Cl_delta_a
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S
        b = self.aircraft.geom.b

        Ldelta_a = 0.5 * rho * S * b * u_s ** 2 * Cl_delta_a
        return Ldelta_a
    
    def calculate_Ldelta_r(self) -> float:
        """
        Calculate Ldelta_r derivative.

        :rtype: float
        """
        Cl_delta_r = self.aircraft.stab_coeffs.latdir.Cl_delta_r
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S
        b = self.aircraft.geom.b

        Ldelta_r = 0.5 * rho * S * b * u_s ** 2 * Cl_delta_r
        return Ldelta_r
    
    def calculate_Nv(self) -> float:
        """
        Calculate Nv derivative.

        :rtype: float
        """
        Cn_beta = self.aircraft.stab_coeffs.latdir.Cn_beta
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S
        b = self.aircraft.geom.b

        Nv = 0.5 * rho * S * b * u_s * Cn_beta
        return Nv
    
    def calculate_Np(self) -> float:
        """
        Calculate Np derivative.

        :rtype: float
        """
        Cn_p = self.aircraft.stab_coeffs.latdir.Cn_p
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S
        b = self.aircraft.geom.b

        Np = 0.25 * rho * S * b ** 2 * u_s * Cn_p
        return Np
    
    def calculate_Nr(self) -> float:
        """
        Calculate Nr derivative.

        :rtype: float
        """
        Cn_r = self.aircraft.stab_coeffs.latdir.Cn_r
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S
        b = self.aircraft.geom.b

        Nr = 0.25 * rho * S * b ** 2 * u_s * Cn_r
        return Nr
    
    def calculate_Ndelta_a(self) -> float:
        """
        Calculate Ndelta_a derivative.

        :rtype: float
        """
        Cn_delta_a = self.aircraft.stab_coeffs.latdir.Cn_delta_a
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S
        b = self.aircraft.geom.b

        Ndelta_a = 0.5 * rho * S * b * u_s ** 2 * Cn_delta_a
        return Ndelta_a
    
    def calculate_Ndelta_r(self) -> float:
        """
        Calculate Ndelta_r derivative.

        :rtype: float
        """
        Cn_delta_r = self.aircraft.stab_coeffs.latdir.Cn_delta_r
        rho = self.aircraft.flight_cond.rho
        u_s = self.aircraft.flight_cond.u_s
        S = self.aircraft.geom.S
        b = self.aircraft.geom.b

        Ndelta_r = 0.5 * rho * S * b * u_s ** 2 * Cn_delta_r
        return Ndelta_r


class Aircraft:
    def __init__(self, model):
        main_dir = Path.cwd()

        self.model = model

        # Load aircraft data
        data_path = Path(main_dir) / "aircraft" / f"{model}.json"
        with open(data_path) as data:
            aircraft_data = json.load(data)
        aircraft_data = self.convert_to_SI(aircraft_data)
        aircraft_data['flight_cond'].update({'rho': 2 * aircraft_data['flight_cond']['q'] / aircraft_data['flight_cond']['u_s'] ** 2})

        self.geom = Geometry(**aircraft_data['geom'])
        self.mass_prop = MassProperties(**aircraft_data['mass_prop'])
        self.flight_cond = FlightCondition(**aircraft_data['flight_cond'])
        self.cond_coeffs = ConditionCoefficients(**aircraft_data['cond_coeffs'])

        long = LongitudinalCoefficients(**aircraft_data['stab_coeffs']['long'])
        latdir = LateralDirectionalCoefficients(**aircraft_data['stab_coeffs']['latdir'])

        self.stab_coeffs = StabilityCoefficients(long=long, latdir=latdir)

        # Prepare for stability derivatives calculation
        self.stab_der = StabilityDerivatives(self)

        a = 1

    def convert_to_SI(self, aircraft_data) -> dict:
        """
        Convert every parameter to SI units.

        :param aircraft_data: Aircraft data.
        :type aircraft_data: dict

        :rtype: dict
        """
        # Geometry
        aircraft_data['geom']['S'] = convert_length(convert_length(aircraft_data['geom']['S'], from_=LengthUnit.FEET), from_=LengthUnit.FEET)
        aircraft_data['geom']['c'] = convert_length(aircraft_data['geom']['c'], from_=LengthUnit.FEET)
        aircraft_data['geom']['b'] = convert_length(aircraft_data['geom']['b'], from_=LengthUnit.FEET)
        aircraft_data['geom']['eps_s'] = convert_angle(aircraft_data['geom']['eps_s'], from_=AngleUnit.DEGREES)

        # Flight condition
        aircraft_data['flight_cond']['h'] = convert_length(aircraft_data['flight_cond']['h'], from_=LengthUnit.FEET)
        aircraft_data['flight_cond']['u_s'] = convert_length(aircraft_data['flight_cond']['u_s'], from_=LengthUnit.FEET)
        aircraft_data['flight_cond']['q'] = convert_weight(aircraft_data['flight_cond']['q'], from_=WeightUnit.POUNDS) * g
        aircraft_data['flight_cond']['q'] = 1 / convert_length(convert_length(1 / aircraft_data['flight_cond']['q'], from_=LengthUnit.FEET), from_=LengthUnit.FEET)
        aircraft_data['flight_cond']['theta_s'] = convert_angle(aircraft_data['flight_cond']['theta_s'], from_=AngleUnit.DEGREES)

        # Mass properties
        aircraft_data['mass_prop']['W'] = convert_weight(aircraft_data['mass_prop']['W'], from_=WeightUnit.POUNDS)
        aircraft_data['mass_prop']['I_xx'] = convert_weight(aircraft_data['mass_prop']['I_xx'], from_=WeightUnit.SLUG)
        aircraft_data['mass_prop']['I_xx'] = convert_length(convert_length(aircraft_data['mass_prop']['I_xx'], from_=LengthUnit.FEET), from_=LengthUnit.FEET)
        aircraft_data['mass_prop']['I_yy'] = convert_weight(aircraft_data['mass_prop']['I_yy'], from_=WeightUnit.SLUG)
        aircraft_data['mass_prop']['I_yy'] = convert_length(convert_length(aircraft_data['mass_prop']['I_yy'], from_=LengthUnit.FEET), from_=LengthUnit.FEET)
        aircraft_data['mass_prop']['I_zz'] = convert_weight(aircraft_data['mass_prop']['I_zz'], from_=WeightUnit.SLUG)
        aircraft_data['mass_prop']['I_zz'] = convert_length(convert_length(aircraft_data['mass_prop']['I_zz'], from_=LengthUnit.FEET), from_=LengthUnit.FEET)
        aircraft_data['mass_prop']['I_xz'] = convert_weight(aircraft_data['mass_prop']['I_xz'], from_=WeightUnit.SLUG)
        aircraft_data['mass_prop']['I_xz'] = convert_length(convert_length(aircraft_data['mass_prop']['I_xz'], from_=LengthUnit.FEET), from_=LengthUnit.FEET)

        return aircraft_data