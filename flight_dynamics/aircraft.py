# -*- coding: utf-8 -*-
"""
Aircraft
=========

Aircraft parameters and calculations.

"""

from pathlib import Path
import json
from utilities.constants import g 
from utilities.unit_conversion import LengthUnit, WeightUnit, AngleUnit, convert_length, convert_weight, convert_angle
from .data_classes import (Geometry, MassProperties, FlightCondition, ConditionCoefficients, LongitudinalCoefficients, 
                           LateralDirectionalCoefficients, StabilityCoefficients, LinearizationParameters, StabilityDerivatives)
from .stability_derivatives import StabilityDerivativesCalculator


class Aircraft:
    def __init__(self, model):
        main_dir = Path.cwd()

        self.model = model

        # Load aircraft data
        data_path = Path(main_dir) / "flight_dynamics" / f"{model}.json"
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

        # Stability derivatives calculation
        params = LinearizationParameters(
            geom=self.geom,
            mass_prop=self.mass_prop,
            flight_cond=self.flight_cond,
            cond_coeffs=self.cond_coeffs,
            stab_coeffs=self.stab_coeffs
        )

        long = StabilityDerivativesCalculator.calculate_longitudinal(params)
        latdir = StabilityDerivativesCalculator.calculate_lateral_directional(params)
        self.stab_der = StabilityDerivatives(long=long, latdir=latdir)


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