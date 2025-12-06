from pathlib import Path
import os
import json
import isacalc as isa
from math import sin, cos
from utilities.constants import GRAVITY as g 

class Aircraft:
    def __init__(self, model):
        main_dir = Path.cwd()

        # Load aircraft data
        with open(os.path.join(main_dir, 'aircraft', f'{model}.json'), 'r') as data:
            aircraft_data = json.load(data)
        self.geom = aircraft_data['geom']
        self.mass_prop = aircraft_data['mass_prop']
        self.cond_coeffs = aircraft_data['cond_coeffs']
        self.stab_coeffs = aircraft_data['stab_coeffs']
        atm = isa.Atmosphere().calculate(h=aircraft_data['flight_cond']['h'])
        aircraft_data['flight_cond'].update({p: val for p, val in zip(['T', 'p', 'rho', 'a', 'mu'], atm)})
        self.stab_der = StabilityDerivatives(self, aircraft_data['flight_cond'])


class StabilityDerivatives:
    def __init__(self, Aircraft, FlightCondition):
        self.Aircraft = Aircraft
        self.FlightCondition = FlightCondition
        
    def calculate_all(self):
        for function in dir(self):
            if function.startswith('calculate_') and not function.endswith('all'):
                value = getattr(self, function)()
                setattr(self, function.replace('calculate_', ''), value)
    
    # Longitudinal derivatives
    def calculate_Xu(self):
        W = self.Aircraft.mass_prop['W']
        theta_s = self.FlightCondition['theta_s']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']
        CT_u = self.Aircraft.stab_coeffs['long']['CT_u']
        CD_u = self.Aircraft.stab_coeffs['long']['CD_u']
        eps_s = self.Aircraft.geom['eps_s']

        CX_s = W * g * sin(theta_s) / (0.5 * rho * S * u_s ** 2)
        CX_u = CT_u * cos(eps_s) - CD_u

        Xu = rho * S * u_s * (CX_s + 0.5 * u_s * CX_u)
        return Xu
    
    def calculate_Xw(self):
        CL_s = self.Aircraft.cond_coeffs['CL_s']
        CD_alpha = self.Aircraft.stab_coeffs['long']['CD_alpha']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']

        CX_alpha = CL_s - CD_alpha

        Xw = 0.5 * rho * S * u_s * CX_alpha
        return Xw
    
    def calculate_Xdelta_e(self):
        CD_delta_e = self.Aircraft.stab_coeffs['long']['CD_delta_e']
        CX_delta_e = - CD_delta_e
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']

        Xdelta_e = 0.5 * rho * S * u_s ** 2 * CX_delta_e
        return Xdelta_e
    
    def calculate_Zu(self):
        W = self.Aircraft.mass_prop['W']
        theta_s = self.FlightCondition['theta_s']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']
        CT_u = self.Aircraft.stab_coeffs['long']['CT_u']
        CL_u = self.Aircraft.stab_coeffs['long']['CD_u']
        eps_s = self.Aircraft.geom['eps_s']

        CZ_s = - W * g * cos(theta_s) / (0.5 * rho * S * u_s ** 2)
        CZ_u = - CT_u * sin(eps_s) - CL_u

        Zu = rho * S * u_s * (CZ_s + 0.5 * u_s * CZ_u)
        return Zu
    
    def calculate_Zw(self):
        CL_alpha = self.Aircraft.stab_coeffs['long']['CL_alpha']
        CD_s =  self.Aircraft.cond_coeffs['CD_s']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']

        CZ_alpha = - CL_alpha - CD_s

        Zw = 0.5 * rho * S * u_s * CZ_alpha
        return Zw
    
    def calculate_Zw_dot(self):
        CL_alpha_dot = self.Aircraft.stab_coeffs['long']['CL_alpha_dot']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']

        CZ_alpha_dot = - CL_alpha_dot

        Zw_dot = 0.5 * rho * S * u_s * CZ_alpha_dot
        return Zw_dot
    
    def calculate_Zq(self):
        CL_q = self.Aircraft.stab_coeffs['long']['CL_q']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']

        CZ_q = - CL_q

        Zq = 0.5 * rho * S * u_s ** 2 * CZ_q
        return Zq
    
    def calculate_Zdelta_e(self):
        CL_delta_e = self.Aircraft.stab_coeffs['long']['CL_delta_e']
        CZ_delta_e = - CL_delta_e
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']

        Zdelta_e = 0.5 * rho * S * u_s ** 2 * CZ_delta_e
        return Zdelta_e
    
    def calculate_Mu(self):
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']
        c = self.Aircraft.geom['c']
        Cm_u = self.Aircraft.stab_coeffs['long']['Cm_u']

        Mu = 0.5* rho * S * c * u_s ** 2 * Cm_u
        return Mu
    
    def calculate_Mw(self):
        Cm_alpha = self.Aircraft.stab_coeffs['long']['Cm_alpha']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']
        c = self.Aircraft.geom['c']

        Mw = 0.5 * rho * S * c * u_s * Cm_alpha
        return Mw
    
    def calculate_Mw_dot(self):
        CL_alpha_dot = self.Aircraft.stab_coeffs['long']['CL_alpha_dot']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']
        c = self.Aircraft.geom['c']

        CZ_alpha_dot = - CL_alpha_dot

        Mw_dot = 0.5 * rho * S * c * u_s * CZ_alpha_dot
        return Mw_dot
    
    def calculate_Mq(self):
        Cm_q = self.Aircraft.stab_coeffs['long']['Cm_q']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']
        c = self.Aircraft.geom['c']

        Mq = 0.5 * rho * S * c * u_s ** 2 * Cm_q
        return Mq
    
    def calculate_Mdelta_e(self):
        Cm_delta_e = self.Aircraft.stab_coeffs['long']['Cm_delta_e']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']
        c = self.Aircraft.geom['c']

        Mdelta_e = 0.5 * rho * S * c * u_s ** 2 * Cm_delta_e
        return Mdelta_e
    
    # Lateral-directional derivatives
    def calculate_Yv(self):
        CY_beta = self.Aircraft.stab_coeffs['latdir']['CY_beta']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']

        Yv = 0.5 * rho * S * u_s * CY_beta
        return Yv
    
    def calculate_Yp(self):
        CY_p = self.Aircraft.stab_coeffs['latdir']['CY_p']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']

        Yp = 0.5 * rho * S * u_s ** 2 * CY_p
        return Yp
    
    def calculate_Yr(self):
        CY_r = self.Aircraft.stab_coeffs['latdir']['CY_r']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']

        Yr = 0.5 * rho * S * u_s ** 2 * CY_r
        return Yr
    
    def calculate_Ydelta_r(self):
        CY_delta_r = self.Aircraft.stab_coeffs['latdir']['CY_delta_r']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']

        Ydelta_r = 0.5 * rho * S * u_s ** 2 * CY_delta_r
        return Ydelta_r
    
    def calculate_Lv(self):
        Cl_beta = self.Aircraft.stab_coeffs['latdir']['Cl_beta']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']

        Lv = 0.5 * rho * S * u_s * Cl_beta
        return Lv
    
    def calculate_Lp(self):
        Cl_p = self.Aircraft.stab_coeffs['latdir']['Cl_p']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']
        b = self.Aircraft.geom['b']

        Lp = 0.5 * rho * S * b * u_s ** 2 * Cl_p
        return Lp
    
    def calculate_Lr(self):
        Cl_r = self.Aircraft.stab_coeffs['latdir']['Cl_r']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']
        b = self.Aircraft.geom['b']

        Lr = 0.5 * rho * S * b * u_s ** 2 * Cl_r
        return Lr
    
    def calculate_Ldelta_a(self):
        Cl_delta_a = self.Aircraft.stab_coeffs['latdir']['Cl_delta_a']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']
        b = self.Aircraft.geom['b']

        Ldelta_a = 0.5 * rho * S * b * u_s ** 2 * Cl_delta_a
        return Ldelta_a
    
    def calculate_Ldelta_r(self):
        Cl_delta_r = self.Aircraft.stab_coeffs['latdir']['Cl_delta_r']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']
        b = self.Aircraft.geom['b']

        Ldelta_r = 0.5 * rho * S * b * u_s ** 2 * Cl_delta_r
        return Ldelta_r
    
    def calculate_Nv(self):
        Cn_beta = self.Aircraft.stab_coeffs['latdir']['Cn_beta']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']

        Nv = 0.5 * rho * S * u_s * Cn_beta
        return Nv
    
    def calculate_Np(self):
        Cn_p = self.Aircraft.stab_coeffs['latdir']['Cn_p']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']
        b = self.Aircraft.geom['b']

        Np = 0.5 * rho * S * b * u_s ** 2 * Cn_p
        return Np
    
    def calculate_Nr(self):
        Cn_r = self.Aircraft.stab_coeffs['latdir']['Cn_r']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']
        b = self.Aircraft.geom['b']

        Nr = 0.5 * rho * S * b * u_s ** 2 * Cn_r
        return Nr
    
    def calculate_Ndelta_a(self):
        Cn_delta_a = self.Aircraft.stab_coeffs['latdir']['Cn_delta_a']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']
        b = self.Aircraft.geom['b']

        Ndelta_a = 0.5 * rho * S * b * u_s ** 2 * Cn_delta_a
        return Ndelta_a
    
    def calculate_Ndelta_r(self):
        Cn_delta_r = self.Aircraft.stab_coeffs['latdir']['Cn_delta_r']
        rho = self.FlightCondition['rho']
        u_s = self.FlightCondition['u_s']
        S = self.Aircraft.geom['S']
        b = self.Aircraft.geom['b']

        Ndelta_r = 0.5 * rho * S * b * u_s ** 2 * Cn_delta_r
        return Ndelta_r
    
    