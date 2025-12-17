# -*- coding: utf-8 -*-
"""
Units
=====

A unit conversion module.

.. rubric:: Examples

Import the units module by just:

.. code-block:: python

    >>> import utilities.units as units

Then units can be converted as follows:


.. code-block:: python

    >>> units.convert_length(1000, units.LengthUnit.METERS, units.LengthUnit.NAUTICAL_MILES)
    0.54
    >>> units.convert_length(100, to = units.LengthUnit.FEET)
    328.084
    >>> units.convert_length(100, from_ = units.LengthUnit.FEET)
    30.48

.. rubric:: Members
"""

from math import pi, tan, atan
from enum import Enum, auto

from utilities import constants

#: Pounds to kilograms conversion factor
POUNDS_TO_KILOGRAMS = 0.453592
#: Celsius to Kelvin conversion constant
CELSIUS_TO_KELVIN = 273.15
#: Fahrenheit to Kelvin conversion constant
FAHRENHEIT_TO_KELVIN_C = 459.67
#: Fahrenheit to Kelvin conversion factor
FAHRENHEIT_TO_KELVIN_F = 1.0/1.8
#: Feet to meters conversion factor
FEET_TO_METERS = 0.3048
#: Nautical miles to meters conversion factor
NAUTICAL_MILES_TO_METERS = 1852.0
#: Kilometers to meters conversion factor
KILOMETERS_TO_METERS = 1000.0
#: Minutes to seconds conversion factor
MINUTES_TO_SECONDS = 60.0
#: Hours to seconds conversion factor
HOURS_TO_SECONDS = 3600.0
#: Knots to meters per second conversion factor
KNOTS_TO_METERS_PER_SECOND = 0.514444
#: Slug to kg conversion factor
SLUG_TO_KILOGRAMS = 14.59390


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
    Convert weight to the selected unit.

    .. code-block:: python

        >>> convert_weight(2000, WeightUnit.KILOGRAMS, WeightUnit.POUNDS)
        4409.2488
        >>> convert_weight(2000, to = WeightUnit.POUNDS)
        4409.2488

    :param value: Weight value to convert.
    :type value: float
    :param from_: Original weight unit.
    :type from_: :class:`utilities.units.WeightUnit`, optional
    :param to: Desired weight unit.
    :type to: :class:`utilities.units.WeightUnit`, optional
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
    Convert length to the selected unit.

    .. code-block:: python

        >>> convert_length(1000, LengthUnit.METERS, LengthUnit.NAUTICAL_MILES)
        0.54
        >>> convert_length(100, to = LengthUnit.FEET)
        328.084
        >>> convert_length(100, from_ = LengthUnit.FEET)
        30.48

    :param value: Length value to convert.
    :type value: float
    :param from_: Original length unit.
    :type from_: :class:`utilities.units.LengthUnit`, optional
    :param to: Desired length unit.
    :type to: :class:`utilities.units.LengthUnit`, optional
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
    Convert temperature to the selected unit.

    .. code-block:: python

        >>> convert_temperature(300.0, TemperatureUnit.KELVIN, TemperatureUnit.CELSIUS)
        26.850
        >>> convert_temperature(300.0, to = TemperatureUnit.FAHRENHEIT)
        80.330
        >>> convert_temperature(15.0, from_ = TemperatureUnit.CELSIUS)
        288.150

    :param value: Temperature value to convert.
    :type value: float
    :param from_: Original temperature unit.
    :type from_: :class:`utilities.units.TemperatureUnit`, optional
    :param to: Desired temperature unit.
    :type to: :class:`utilities.units.TemperatureUnit`, optional
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
    Convert speed to the selected unit.

    .. code-block:: python

        >>> convert_speed(100.0, SpeedUnit.METER_PER_SECOND, SpeedUnit.KNOTS)
        194.385
        >>> convert_speed(100.0, to = SpeedUnit.KILOMETERS_PER_HOUR)
        360.000
        >>> convert_speed(2000.0, from_ = SpeedUnit.FEET_PER_MINUTE)
        10.160

    :param value: Speed value to convert.
    :type value: float
    :param from_: Original speed unit.
    :type from_: :class:`utilities.units.SpeedUnit`, optional
    :param to: Desired speed unit.
    :type to: :class:`utilities.units.SpeedUnit`, optional
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
    Convert force to the selected unit.

    .. code-block:: python

        >>> convert_force(1000.0, ForceUnit.NEWTONS, ForceUnit.KILOGRAMS)
        101.972
        >>> convert_force(1000.0, to = ForceUnit.POUNDS)
        224.809
        >>> convert_force(100.0, from_ = ForceUnit.KILOGRAMS)
        980.665

    :param value: Force value to convert.
    :type value: float
    :param from_: Original force unit.
    :type from_: :class:`utilities.units.ForceUnit`, optional
    :param to: Desired force unit.
    :type to: :class:`utilities.units.ForceUnit`, optional
    :return: Converted value
    :rtype: float
    """

    converted = value

    if from_ is ForceUnit.KILOGRAMS:
        converted *= constants.GRAVITY
    elif from_ is ForceUnit.POUNDS:
        converted *= (POUNDS_TO_KILOGRAMS * constants.GRAVITY)

    if to is ForceUnit.KILOGRAMS:
        converted /= constants.GRAVITY
    elif to is ForceUnit.POUNDS:
        converted /= (POUNDS_TO_KILOGRAMS * constants.GRAVITY)

    return converted


def convert_angle(value, from_ = AngleUnit.RADIANS, to = AngleUnit.RADIANS):
    """
    Convert angles to the selected unit.

    .. code-block:: python

        >>> convert_angle(0.5, AngleUnit.RADIANS, AngleUnit.DEGREES)
        28.648
        >>> convert_angle(0.05, to = AngleUnit.GRADIENT)
        5.004
        >>> convert_angle(2.4, from_ = AngleUnit.GRADIENT)
        0.023995

    :param value: Angle value to convert.
    :type value: float
    :param from_: Original angle unit.
    :type from_: :class:`utilities.units.AngleUnit`, optional
    :param to: Desired angle unit.
    :type to: :class:`utilities.units.AngleUnit`, optional
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