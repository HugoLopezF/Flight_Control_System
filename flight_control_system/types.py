from enum import Enum

class Axis(str, Enum):
    LONG = "long"
    LATDIR = "latdir"

class InputChannel(str, Enum):
    ELEVATOR = "elevator"
    AILERON = "aileron"
    RUDDER = "rudder"

class ActuatorOrder(str, Enum):
    FIRST = "first"
    SECOND = "second"