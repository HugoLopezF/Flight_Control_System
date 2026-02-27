from enum import Enum

class Axis(str, Enum):
    LONG = "long"
    LATDIR = "latdir"

class OutputChannel(str, Enum):
    U = "u"
    ALPHA = "alpha"
    THETA = "theta"
    Q = "q"
    BETA = "beta"
    P = "p"
    R = "r"
    PHI = "phi"

class InputChannel(str, Enum):
    ELEVATOR = "elevator"
    AILERON = "aileron"
    RUDDER = "rudder"

class ActuatorOrder(str, Enum):
    IDEAL = "ideal"
    FIRST = "first"
    SECOND = "second"

class SensorModel(str, Enum):
    IDEAL = "ideal"
    PADE = "pade"
    FIRST_ORDER = "first_order"

class FilterType(str, Enum):
    IDEAL = "ideal"
    LOWPASS = "lowpass"
    HIGHPASS = "highpass"