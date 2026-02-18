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


def _convert_linear(value, from_unit, to_unit, factors):
    """
    Convert units where every unit scales linearly from a base unit.

    """
    return value * factors[from_unit] / factors[to_unit]


def convert_weight(value, from_=WeightUnit.KILOGRAMS, to=WeightUnit.KILOGRAMS):
    """
    Convert weight to the desired unit.

    :param value: Weight value to convert.
    :type value: float
    :param from_: Input weight unit.
    :type from_: :class:`utilities.unit_conversion.WeightUnit`, optional
    :param to: Output weight unit.
    :type to: :class:`utilities.unit_conversion.WeightUnit`, optional
    :return: Converted value
    :rtype: float
    """

    return _convert_linear(value, from_, to, _WEIGHT_TO_KILOGRAMS)


def convert_length(value, from_=LengthUnit.METERS, to=LengthUnit.METERS):
    """
    Convert length to the desired unit.

    :param value: Length value to convert.
    :type value: float
    :param from_: Input length unit.
    :type from_: :class:`utilities.unit_conversion.LengthUnit`, optional
    :param to: Output length unit.
    :type to: :class:`utilities.unit_conversion.LengthUnit`, optional
    :return: Converted value
    :rtype: float
    """

    return _convert_linear(value, from_, to, _LENGTH_TO_METERS)


def convert_temperature(value, from_=TemperatureUnit.KELVIN, to=TemperatureUnit.KELVIN):
    """
    Convert temperature to the desired unit.

    :param value: Temperature value to convert.
    :type value: float
    :param from_: Input temperature unit.
    :type from_: :class:`utilities.unit_conversion.TemperatureUnit`, optional
    :param to: Output temperature unit.
    :type to: :class:`utilities.unit_conversion.TemperatureUnit`, optional
    :return: Converted value
    :rtype: float
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


def convert_speed(value, from_=SpeedUnit.METER_PER_SECOND, to=SpeedUnit.METER_PER_SECOND):
    """
    Convert speed to the desired unit.

    :param value: Speed value to convert.
    :type value: float
    :param from_: Input speed unit.
    :type from_: :class:`utilities.unit_conversion.SpeedUnit`, optional
    :param to: Output speed unit.
    :type to: :class:`utilities.unit_conversion.SpeedUnit`, optional
    :return: Converted value
    :rtype: float
    """

    return _convert_linear(value, from_, to, _SPEED_TO_METERS_PER_SECOND)


def convert_force(value, from_=ForceUnit.NEWTONS, to=ForceUnit.NEWTONS):
    """
    Convert force to the desired unit.

    :param value: Force value to convert.
    :type value: float
    :param from_: Input force unit.
    :type from_: :class:`utilities.unit_conversion.ForceUnit`, optional
    :param to: Output force unit.
    :type to: :class:`utilities.unit_conversion.ForceUnit`, optional
    :return: Converted value
    :rtype: float
    """

    return _convert_linear(value, from_, to, _FORCE_TO_NEWTONS)


def convert_angle(value, from_=AngleUnit.RADIANS, to=AngleUnit.RADIANS):
    """
    Convert angle to the desired unit.

    :param value: Angle value to convert.
    :type value: float
    :param from_: Input angle unit.
    :type from_: :class:`utilities.unit_conversion.AngleUnit`, optional
    :param to: Output angle unit.
    :type to: :class:`utilities.unit_conversion.AngleUnit`, optional
    :return: Converted value
    :rtype: float
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
