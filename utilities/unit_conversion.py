# -*- coding: utf-8 -*-
"""
Unit Conversion
-------------

Unit conversion module.

"""

from enum import Enum, auto
from math import atan, pi, tan

from utilities import constants

# Length conversion factors
FEET_TO_METERS = 0.3048
KILOMETERS_TO_METERS = 1000.0
NAUTICAL_MILES_TO_METERS = 1852.0

# Mass conversion factors
POUNDS_TO_KILOGRAMS = 0.45359237
SLUG_TO_KILOGRAMS = 14.59390

# Time conversion factors
MINUTES_TO_SECONDS = 60.0
HOURS_TO_SECONDS = 3600.0

# Speed conversion factors
KNOTS_TO_METERS_PER_SECOND = 0.5144444444

# Temperature conversion constants
CELSIUS_TO_KELVIN = 273.15
FAHRENHEIT_TO_KELVIN_C = 459.67
FAHRENHEIT_TO_KELVIN_F = 1.0 / 1.8


class WeightUnit(Enum):
    """
    Available weight units.
    """
    KILOGRAMS = auto()
    POUNDS = auto()
    SLUG = auto()


class LengthUnit(Enum):
    """
    Available length units.
    """
    METERS = auto()
    FEET = auto()
    NAUTICAL_MILES = auto()
    KILOMETERS = auto()


class TemperatureUnit(Enum):
    """
    Available temperature units.
    """
    KELVIN = auto()
    CELSIUS = auto()
    FAHRENHEIT = auto()


class SpeedUnit(Enum):
    """
    Available speed units.
    """
    METER_PER_SECOND = auto()
    KILOMETERS_PER_HOUR = auto()
    KNOTS = auto()
    FEET_PER_MINUTE = auto()


class ForceUnit(Enum):
    """
    Available force units.
    """
    NEWTONS = auto()
    KILOGRAMS = auto()
    POUNDS = auto()


class AngleUnit(Enum):
    """
    Available angle units.
    """
    RADIANS = auto()
    DEGREES = auto()
    GRADIENT = auto()


_WEIGHT_TO_KILOGRAMS = {
    WeightUnit.KILOGRAMS: 1.0,
    WeightUnit.POUNDS: POUNDS_TO_KILOGRAMS,
    WeightUnit.SLUG: SLUG_TO_KILOGRAMS,
}

_LENGTH_TO_METERS = {
    LengthUnit.METERS: 1.0,
    LengthUnit.FEET: FEET_TO_METERS,
    LengthUnit.NAUTICAL_MILES: NAUTICAL_MILES_TO_METERS,
    LengthUnit.KILOMETERS: KILOMETERS_TO_METERS,
}

_SPEED_TO_METERS_PER_SECOND = {
    SpeedUnit.METER_PER_SECOND: 1.0,
    SpeedUnit.KILOMETERS_PER_HOUR: KILOMETERS_TO_METERS / HOURS_TO_SECONDS,
    SpeedUnit.KNOTS: KNOTS_TO_METERS_PER_SECOND,
    SpeedUnit.FEET_PER_MINUTE: FEET_TO_METERS / MINUTES_TO_SECONDS,
}

_FORCE_TO_NEWTONS = {
    ForceUnit.NEWTONS: 1.0,
    ForceUnit.KILOGRAMS: constants.g,
    ForceUnit.POUNDS: POUNDS_TO_KILOGRAMS * constants.g,
}


def _convert_linear(
    value: float, 
    from_unit: WeightUnit | LengthUnit | TemperatureUnit | SpeedUnit | ForceUnit | AngleUnit, 
    to_unit: WeightUnit | LengthUnit | TemperatureUnit | SpeedUnit | ForceUnit | AngleUnit, 
    factors: dict,
) -> float:
    """
    Convert units where every unit scales linearly from a base unit.

    Parameters
    ----------
    value :float
        Value to convert.
    from_unit : WeightUnit | LengthUnit | TemperatureUnit | SpeedUnit | ForceUnit | AngleUnit
        Input unit.
    to_unit : WeightUnit | LengthUnit | TemperatureUnit | SpeedUnit | ForceUnit | AngleUnit
        Output unit.
    factors : dict
        Factor dictionary for the selected unit type.

    Returns
    -------
    float
        Converted value.
    """

    return value * factors[from_unit] / factors[to_unit]


def convert_weight(
    value, 
    from_: WeightUnit=WeightUnit.KILOGRAMS, 
    to: WeightUnit=WeightUnit.KILOGRAMS,
) -> float:
    """
    Convert weight to the desired unit.

    Parameters
    ----------
    value :float
        Weight value to convert.
    from_ : WeightUnit
        Input weight unit.
    to : WeightUnit
        Output weight unit.

    Returns
    -------
    float
        Converted value.
    """

    return _convert_linear(value, from_, to, _WEIGHT_TO_KILOGRAMS)


def convert_length(
    value, 
    from_: LengthUnit=LengthUnit.METERS, 
    to: LengthUnit=LengthUnit.METERS,
) -> float:
    """
    Convert length to the desired unit.

    Parameters
    ----------
    value :float
        Length value to convert.
    from_ : LengthUnit
        Input length unit.
    to : LengthUnit
        Output length unit.

    Returns
    -------
    float
        Converted value.
    """

    return _convert_linear(value, from_, to, _LENGTH_TO_METERS)


def convert_temperature(
    value, 
    from_: TemperatureUnit=TemperatureUnit.KELVIN, 
    to: TemperatureUnit=TemperatureUnit.KELVIN
) -> float:
    """
    Convert temperature to the desired unit.

    Parameters
    ----------
    value :float
        Temperature value to convert.
    from_ : TemperatureUnit
        Input temperature unit.
    to : TemperatureUnit
        Output temperature unit.

    Returns
    -------
    float
        Converted value.
    """

    converted = value

    if from_ is TemperatureUnit.CELSIUS:
        converted += CELSIUS_TO_KELVIN
    elif from_ is TemperatureUnit.FAHRENHEIT:
        converted = (converted + FAHRENHEIT_TO_KELVIN_C) * FAHRENHEIT_TO_KELVIN_F

    if to is TemperatureUnit.CELSIUS:
        converted -= CELSIUS_TO_KELVIN
    elif to is TemperatureUnit.FAHRENHEIT:
        converted = (converted / FAHRENHEIT_TO_KELVIN_F) - FAHRENHEIT_TO_KELVIN_C

    return converted


def convert_speed(
    value, 
    from_: SpeedUnit=SpeedUnit.METER_PER_SECOND, 
    to: SpeedUnit=SpeedUnit.METER_PER_SECOND
) -> float:
    """
    Convert speed to the desired unit.

    Parameters
    ----------
    value :float
        Speed value to convert.
    from_ : SpeedUnit
        Input speed unit.
    to : SpeedUnit
        Output speed unit.

    Returns
    -------
    float
        Converted value.
    """

    return _convert_linear(value, from_, to, _SPEED_TO_METERS_PER_SECOND)


def convert_force(
    value, 
    from_: ForceUnit=ForceUnit.NEWTONS, 
    to: ForceUnit=ForceUnit.NEWTONS
) -> float:
    """
    Convert force to the desired unit.

    Parameters
    ----------
    value :float
        Force value to convert.
    from_ : ForceUnit
        Input force unit.
    to : ForceUnit
        Output force unit.

    Returns
    -------
    float
        Converted value.
    """

    return _convert_linear(value, from_, to, _FORCE_TO_NEWTONS)


def convert_angle(
    value, 
    from_: AngleUnit=AngleUnit.RADIANS, 
    to: AngleUnit=AngleUnit.RADIANS
) -> float:
    """
    Convert angle to the desired unit.

    Parameters
    ----------
    value :float
        Angle value to convert.
    from_ : AngleUnit
        Input angle unit.
    to : AngleUnit
        Output angle unit.

    Returns
    -------
    float
        Converted value.
    """

    converted = value

    if from_ is AngleUnit.DEGREES:
        converted *= pi / 180.0
    elif from_ is AngleUnit.GRADIENT:
        converted = atan(converted / 100.0)

    if to is AngleUnit.DEGREES:
        converted *= 180.0 / pi
    elif to is AngleUnit.GRADIENT:
        converted = tan(converted) * 100.0

    return converted
