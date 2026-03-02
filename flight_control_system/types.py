from enum import Enum

class Axis(str, Enum):
    """
    Dynamic system axes.

    Attributes
    ----------
    LONG : str
        Longitudinal axis.
    LATDIR : str
        Lateral-directional axis.
    """
    
    LONG = "long"
    LATDIR = "latdir"

class OutputChannel(str, Enum):
    """
    Dynamic system output channels.

    Attributes
    ----------
    U : str
        x-axis aerodynamic speed.
    ALPHA : str
        Aircraft's angle of attack.
    THETA : str
        Pitch angle.
    Q : str
        Pitch rate.
    BETA : str
        Sideslip angle.
    P : str
        Roll rate.
    R : str
        Yaw rate.
    PHI : str
        Roll angle.
    """

    U = "u"
    ALPHA = "alpha"
    THETA = "theta"
    Q = "q"
    BETA = "beta"
    P = "p"
    R = "r"
    PHI = "phi"

class InputChannel(str, Enum):
    """
    Dynamic system input channels.

    Attributes
    ----------
    ELEVATOR : str
        Elevator deflection.
    AILERON : str
        Aileron deflection.
    RUDDER : str
        Rudder deflection.
    """

    ELEVATOR = "elevator"
    AILERON = "aileron"
    RUDDER = "rudder"

class ActuatorOrder(str, Enum):
    """
    Actuator model order.

    Attributes
    ----------
    IDEAL : str
        Ideal actuator.
    FIRST : str
        First order actuator.
    SECOND : str
        Second order actuator.
    """

    IDEAL = "ideal"
    FIRST = "first"
    SECOND = "second"

class SensorModel(str, Enum):
    """
    Sensor model.

    Attributes
    ----------
    IDEAL : str
        Ideal sensor.
    PADE : str
        Pade approximation sensor.
    FIRST_ORDER : str
        First order sensor.
    """

    IDEAL = "ideal"
    PADE = "pade"
    FIRST_ORDER = "first_order"

class FilterType(str, Enum):
    """
    Filter type.

    Attributes
    ----------
    IDEAL : str
        Ideal filter.
    LOWPASS : str
        Lowpass filter.
    HIGHPASS : str
        Highpass filter.
    """

    IDEAL = "ideal"
    LOWPASS = "lowpass"
    HIGHPASS = "highpass"