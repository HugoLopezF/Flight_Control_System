# -*- coding: utf-8 -*-
"""
Speed
=====

A speed conversion module.

.. rubric:: Examples

Import the speed module by just:

.. code-block:: python

    >>> import utilities.speed as speed

Then speeds can be converted as follows:

.. code-block:: python

    >>> speed.convert_speed_type(100.0, SpeedType.CAS, SpeedType.TAS, 1000.0, 10.0)
    106.67953374075054

.. rubric:: Members
"""

from enum import IntEnum, auto
from math import copysign
from scipy.optimize import fsolve

from utilities import constants
from utilities import atmosphere

#: Initial Mach to consider when solving the Mach
CAS_MACH_INI = 1.0
#: Initial CAS to consider when solving the CAS
MACH_CAS_INI = 340.0


class SpeedType(IntEnum):
    """
    Available speed type.
    """
    TAS = auto()
    EAS = auto()
    CAS = auto()
    MACH = auto()


def calculate_stall_speed(pressure_altitude, delta_isa, weight, wing_surface, cl_max):
    """
    Calculate stall speed.

    .. code-block:: python

        >>> calculate_stall_speed(1000.0, 0.0, 18000, 60, 1.61)
        57.33774210607001

    :param pressure_altitude: Pressure altitude (m)
    :type pressure_altitude: float
    :param delta_isa: ISA deviation
    :type delta_isa: float
    :param weight: Aircraft weight (kg)
    :type weight: float
    :param wing_surface: Wing surface (m2)
    :type wing_surface: float
    :param cl_max: CL maximum (-)
    :type cl_max: float
    :return: Stall speed (m/s)
    :rtype: float
    """

    density = atmosphere.calculate_density(pressure_altitude, delta_isa)

    return ((constants.GRAVITY * weight) / (0.5 * density * wing_surface * cl_max)) ** 0.5


def convert_speed_type(speed, from_, to, pressure_altitude, delta_isa):
    """
    Convert speed to the selected type.

    .. code-block:: python

        >>> convert_speed_type(100.0, SpeedType.TAS, SpeedType.CAS, 5000.0, 10.0)
        76.45558623411038

    :param speed: Speed value to convert.
    :type speed: float
    :param from_: Original speed type.
    :type from_: :class:`utilities.speed.SpeedType`
    :param to: Desired speed type.
    :type to: :class:`utilities.speed.SpeedType`
    :param pressure_altitude: Pressure altitude (m)
    :type pressure_altitude: float
    :param delta_isa: ISA deviation
    :type delta_isa: float
    :return: Converted speed value.
    :rtype: float
    """

    converted = speed

    if from_ is SpeedType.CAS:
        converted = cas_to_mach(converted, pressure_altitude)
        converted = mach_to_tas(converted, pressure_altitude, delta_isa)
    elif from_ is SpeedType.MACH:
        converted = mach_to_tas(converted, pressure_altitude, delta_isa)
    elif from_ is SpeedType.EAS:
        converted = eas_to_tas(converted, pressure_altitude, delta_isa)

    if to is SpeedType.CAS:
        converted = tas_to_mach(converted, pressure_altitude, delta_isa)
        converted = mach_to_cas(converted, pressure_altitude)
    elif to is SpeedType.MACH:
        converted = tas_to_mach(converted, pressure_altitude, delta_isa)
    elif to is SpeedType.EAS:
        converted = tas_to_eas(converted, pressure_altitude, delta_isa)

    return converted


def tas_to_eas(tas, pressure_altitude, delta_isa):
    """
    Convert TAS to EAS.

    .. code-block:: python

        >>> tas_to_eas(100.0, 1000.0, 0.0)
        95.26086961336307

    :param tas: TAS (m/s)
    :type tas: float
    :param pressure_altitude: Pressure altitude (m)
    :type pressure_altitude: float
    :param delta_isa: ISA deviation
    :type delta_isa: float
    :return: EAS (m/s)
    :rtype: float
    """

    sigma = atmosphere.calculate_sigma(pressure_altitude, delta_isa)

    return tas * sigma ** 0.5


def eas_to_tas(eas, pressure_altitude, delta_isa):
    """
    Convert TAS to EAS.

    .. code-block:: python

        >>> eas_to_tas(100.0, 1000.0, 0.0)
        104.97489725411045

    :param eas: EAS (m/s)
    :type eas: float
    :param pressure_altitude: Pressure altitude (m)
    :type pressure_altitude: float
    :param delta_isa: ISA deviation
    :type delta_isa: float
    :return: TAS (m/s)
    :rtype: float
    """

    sigma = atmosphere.calculate_sigma(pressure_altitude, delta_isa)

    return eas / sigma ** 0.5


def tas_to_mach(tas, pressure_altitude, delta_isa):
    """
    Convert TAS to Mach.

    .. code-block:: python

        >>> tas_to_mach(100.0, 1000.0, 0.0)
        0.2972350405685903

    :param tas: TAS (m/s)
    :type tas: float
    :param pressure_altitude: Pressure altitude (m)
    :type pressure_altitude: float
    :param delta_isa: ISA deviation
    :type delta_isa: float
    :return: Mach (-)
    :rtype: float
    """

    # Calculate speed of temperature
    t = atmosphere.calculate_temperature(pressure_altitude, delta_isa)

    # Calculate speed of sound
    a = atmosphere.calculate_speed_of_sound(t)

    return tas / a


def mach_to_tas(mach, pressure_altitude, delta_isa):
    """
    Convert Mach to TAS.

    .. code-block:: python

        >>> mach_to_tas(0.25, 1000.0, 0.0)
        84.1085221721393

    :param mach: Mach (-)
    :type mach: float
    :param pressure_altitude: Pressure altitude (m)
    :type pressure_altitude: float
    :param delta_isa: ISA deviation
    :type delta_isa: float
    :return: TAS (m/s)
    :rtype: float
    """

    # Calculate speed of temperature
    t = atmosphere.calculate_temperature(pressure_altitude, delta_isa)

    # Calculate speed of sound
    a = atmosphere.calculate_speed_of_sound(t)

    return mach * a


def cas_to_mach(cas, pressure_altitude):
    """
    Convert CAS to Mach.

    .. code-block:: python

        >>> cas_to_mach(100.0, 1000.0)
        0.31160541969551825

    :param cas: CAS (m/s)
    :type cas: float
    :param pressure_altitude: Pressure altitude (m)
    :type pressure_altitude: float
    :return: Mach (-)
    :rtype: float
    """

    mach = cas_to_mach_subsonic(cas, pressure_altitude)

    if mach > 1.0:
        mach = cas_to_mach_supersonic(cas, pressure_altitude)

    return mach


def cas_to_mach_subsonic(cas, pressure_altitude):
    """
    Convert CAS to Mach for subsonic regime using Saint Venant equation.

    .. code-block:: python

        >>> cas_to_mach_subsonic(100.0, 1000.0)
        0.31160541969551825

    :param cas: CAS (m/s)
    :type cas: float
    :param pressure_altitude: Pressure altitude (m)
    :type pressure_altitude: float
    :return: Subsonic Mach (-)
    :rtype: float
    """

    # Calculate mu
    mu = (atmosphere.GAMMA - 1.0) / atmosphere.GAMMA

    # Calculate speed of sound
    a0 = atmosphere.calculate_speed_of_sound(atmosphere.T0)

    # Calculate pressure
    p = atmosphere.calculate_pressure(pressure_altitude)

    # Calculate term to term
    mach = (atmosphere.GAMMA - 1.0) / 2.0 * (cas / a0) ** 2.0
    mach = ((1.0 + mach) ** (1.0 / mu)) - 1.0
    mach = ((1.0 + (atmosphere.P0 / p * mach)) ** mu) - 1.0
    mach = (2.0 / (atmosphere.GAMMA - 1.0) * mach) ** 0.5

    # Return the calculated value
    return copysign(mach, cas)


def cas_to_mach_supersonic(cas, pressure_altitude):
    """
    Convert CAS to Mach for supersonic regime using Rayleigh equation.

    .. code-block:: python

        >>> cas_to_mach_supersonic(600.0, 1000.0)
        1.859753835602196

    :param cas: CAS (m/s)
    :type cas: float
    :param pressure_altitude: Pressure altitude (m)
    :type pressure_altitude: float
    :return: Subsonic Mach (-)
    :rtype: float
    """

    # Solve the Mach speed
    mach = fsolve(cas_mach_supersonic_equation, CAS_MACH_INI, args=(cas, pressure_altitude))

    # Return the calculated value
    return copysign(mach[0], cas)


def cas_mach_supersonic_equation(mach, cas, pressure_altitude):
    """
    CAS-Mach equation obtained from CAS definition and Rayleigh equation.

    :param mach: Mach (-)
    :type mach: float
    :param cas: CAS (m/s)
    :type cas: float
    :param pressure_altitude: Pressure altitude (m)
    :type pressure_altitude: float
    :return: Error margin (-)
    :rtype: float
    """

    # Calculate mu
    mu = (atmosphere.GAMMA - 1.0) / atmosphere.GAMMA

    # Calculate speed of sound
    a0 = atmosphere.calculate_speed_of_sound(atmosphere.T0)

    # Calculate pressure
    p = atmosphere.calculate_pressure(pressure_altitude)

    # Calculate left term (in CAS)
    if cas <= a0:
        left_term = (cas / a0) ** 2.0
        left_term *= (atmosphere.GAMMA - 1.0) / 2.0
        left_term = atmosphere.P0 * ((1.0 + left_term) ** (1.0 / mu) - 1.0)
    else:
        left_term_num = (cas / a0) ** 2.0
        left_term_num *= (atmosphere.GAMMA + 1.0) / 2.0
        left_term_num **= 1.0 / mu

        left_term_den = (cas / a0) ** 2.0 - 1.0
        left_term_den *= 2.0 * atmosphere.GAMMA / (atmosphere.GAMMA + 1.0)
        left_term_den += 1.0
        left_term_den **= 1.0 / (atmosphere.GAMMA - 1.0)

        left_term = atmosphere.P0 * (left_term_num / left_term_den - 1.0)

    # Calculate right term (in Mach)
    right_term_num = (atmosphere.GAMMA + 1.0) / 2.0 * mach ** 2.0
    right_term_num **= 1.0 / mu

    right_term_den = 2.0 * atmosphere.GAMMA / (atmosphere.GAMMA + 1.0) * mach ** 2.0
    right_term_den -= (atmosphere.GAMMA - 1.0) / (atmosphere.GAMMA + 1.0)
    right_term_den **= 1.0/(atmosphere.GAMMA - 1.0)

    right_term = p * (right_term_num / right_term_den - 1.0)

    # Return the equation: 0 = f(mach)
    return left_term - right_term


def mach_to_cas(mach, pressure_altitude):
    """
    Convert Mach to CAS.

    .. code-block:: python

        >>> mach_to_cas(0.5, 1000.0)
        160.77583596741968

    :param mach: Mach (-)
    :type mach: float
    :param pressure_altitude: Pressure altitude (m)
    :type pressure_altitude: float
    :return: CAS (m/s)
    :rtype: float
    """

    if mach <= 1.0:
        cas = mach_to_cas_subsonic(mach, pressure_altitude)
    else:
        cas = mach_to_cas_supersonic(mach, pressure_altitude)

    return cas


def mach_to_cas_subsonic(mach, pressure_altitude):
    """
    Convert Mach to CAS for subsonic regime.

    .. code-block:: python

        >>> mach_to_cas_subsonic(0.5, 1000.0)
        160.77583596741968

    :param mach: Mach (-)
    :type mach: float
    :param pressure_altitude: Pressure altitude (m)
    :type pressure_altitude: float
    :return: CAS (m/s)
    :rtype: float
    """

    # Calculate mu
    mu = (atmosphere.GAMMA - 1.0) / atmosphere.GAMMA

    # Calculate speed of sound
    a0 = atmosphere.calculate_speed_of_sound(atmosphere.T0)

    # Calculate pressure
    p = atmosphere.calculate_pressure(pressure_altitude)

    # Calculate term to term
    cas = (atmosphere.GAMMA - 1.0) / 2.0 * mach ** 2.0
    cas = (1.0 + cas) ** (1.0 / mu) - 1.0
    cas = (1.0 + p / atmosphere.P0 * cas) ** mu - 1.0
    cas = (2.0 / (atmosphere.GAMMA - 1.0) * cas) ** 0.5

    # Return calculated value
    return copysign(cas * a0, mach)


def mach_to_cas_supersonic(mach, pressure_altitude):
    """
    Convert Mach to CAS for supersonic regime.

    .. code-block:: python

        >>> mach_to_cas_supersonic(1.5, 1000.0)
        485.19927112667864

    :param mach: Mach (-)
    :type mach: float
    :param pressure_altitude: Pressure altitude (m)
    :type pressure_altitude: float
    :return: CAS (m/s)
    :rtype: float
    """

    # Solve the CAS speed
    cas = fsolve(mach_cas_supersonic_equation, MACH_CAS_INI, args=(mach, pressure_altitude))

    # Return the calculated value
    return copysign(cas[0], mach)


def mach_cas_supersonic_equation(cas, mach, pressure_altitude):
    """
    Mach-CAS equation obtained from CAS definition and Rayleigh equation.

    :param cas: CAS (m/s)
    :type cas: float
    :param mach: Mach (-)
    :type mach: float
    :param pressure_altitude: Pressure altitude (m)
    :type pressure_altitude: float
    :return: Error margin (-)
    :rtype: float
    """

    return cas_mach_supersonic_equation(mach, cas, pressure_altitude)