# -*- coding: utf-8 -*-
"""
Unit Conversion
-------------

Unit conversion module.

"""

from math import pi, tan, atan
from enum import Enum, auto
from utilities import constants

# ft to m conversion factor
FEET_TO_METERS = 0.3048
# km to m conversion factor
KILOMETERS_TO_METERS = 1000.0
# NM to m conversion factor
NAUTICAL_MILES_TO_METERS = 1852.0
# lb to kg conversion factor
POUNDS_TO_KILOGRAMS = 0.45359237
# Slug to kg conversion factor
SLUG_TO_KILOGRAMS = 14.59390
# min to s conversion factor
MINUTES_TO_SECONDS = 60.0
# h to s conversion factor
HOURS_TO_SECONDS = 3600.0
# kt to m/s conversion factor
KNOTS_TO_METERS_PER_SECOND = 0.5144444444
# degC to K conversion constant
CELSIUS_TO_KELVIN = 273.15
# Fahrenheit to K conversion constant
FAHRENHEIT_TO_KELVIN_C = 459.67
# Fahrenheit to K conversion factor
FAHRENHEIT_TO_KELVIN_F = 1.0/1.8


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


def convert_weight(value, from_ = WeightUnit.KILOGRAMS, to = WeightUnit.KILOGRAMS):
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

    converted = value

    if from_ is WeightUnit.POUNDS:
        converted *= POUNDS_TO_KILOGRAMS
    elif from_ is WeightUnit.SLUG:
        converted *= SLUG_TO_KILOGRAMS

    if to is WeightUnit.POUNDS:
        converted /= POUNDS_TO_KILOGRAMS
    elif to is WeightUnit.SLUG:
        converted /= SLUG_TO_KILOGRAMS

    return converted


def convert_length(value, from_ = LengthUnit.METERS, to = LengthUnit.METERS):
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

    converted = value

    if from_ is LengthUnit.FEET:
        converted *= FEET_TO_METERS
    elif from_ is LengthUnit.NAUTICAL_MILES:
        converted *= NAUTICAL_MILES_TO_METERS
    elif from_ is LengthUnit.KILOMETERS:
        converted *= KILOMETERS_TO_METERS

    if to is LengthUnit.FEET:
        converted /= FEET_TO_METERS
    elif to is LengthUnit.NAUTICAL_MILES:
        converted /= NAUTICAL_MILES_TO_METERS
    elif to is LengthUnit.KILOMETERS:
        converted /= KILOMETERS_TO_METERS

    return converted


def convert_temperature(value, from_ = TemperatureUnit.KELVIN, to = TemperatureUnit.KELVIN):
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


def convert_speed(value, from_ = SpeedUnit.METER_PER_SECOND, to = SpeedUnit.METER_PER_SECOND):
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

    converted = value

    if from_ is SpeedUnit.KILOMETERS_PER_HOUR:
        converted = (converted * KILOMETERS_TO_METERS) / HOURS_TO_SECONDS
    elif from_ is SpeedUnit.KNOTS:
        converted *= KNOTS_TO_METERS_PER_SECOND
    elif from_ is SpeedUnit.FEET_PER_MINUTE:
        converted = (converted * FEET_TO_METERS) / MINUTES_TO_SECONDS

    if to is SpeedUnit.KILOMETERS_PER_HOUR:
        converted = (converted * HOURS_TO_SECONDS) / KILOMETERS_TO_METERS
    elif to is SpeedUnit.KNOTS:
        converted /= KNOTS_TO_METERS_PER_SECOND
    elif to is SpeedUnit.FEET_PER_MINUTE:
        converted = (converted * MINUTES_TO_SECONDS) / FEET_TO_METERS

    return converted


def convert_force(value, from_ = ForceUnit.NEWTONS, to = ForceUnit.NEWTONS):
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

    converted = value

    if from_ is ForceUnit.KILOGRAMS:
        converted *= constants.g
    elif from_ is ForceUnit.POUNDS:
        converted *= (POUNDS_TO_KILOGRAMS * constants.g)

    if to is ForceUnit.KILOGRAMS:
        converted /= constants.g
    elif to is ForceUnit.POUNDS:
        converted /= (POUNDS_TO_KILOGRAMS * constants.g)

    return converted


def convert_angle(value, from_ = AngleUnit.RADIANS, to = AngleUnit.RADIANS):
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
        converted *= (pi / 180.0)
    elif from_ is AngleUnit.GRADIENT:
        converted = atan(converted / 100.0)

    if to is AngleUnit.DEGREES:
        converted *= (180.0 / pi)
    elif to is AngleUnit.GRADIENT:
        converted = tan(converted) * 100.0

    return converted