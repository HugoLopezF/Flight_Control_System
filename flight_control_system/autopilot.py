from collections.abc import Mapping
import control
import numpy as np

from .sas_design import compute_DL_gain
from .types import Axis, InputChannel, OutputChannel
from .axis_metadata import AXIS_LABELS
from .actuator import Actuator
from .sensor import Sensor
from .filter import Filter
from .state_space import LinearizedSystem
from .sas import SAS


class AP:
    """
    Autopilot class.

    Attributes
    ----------

    Methods
    ----------
    """

    SUPPORTED_AUTOPILOT = {
        Axis.LONG: (OutputChannel.THETA,),
        Axis.LATDIR: (OutputChannel.PHI,),
    }

    def __init__(self, sas_sys: SAS):
        """
        Autopilot class.

        Parameters
        ----------
        sas_sys : SAS
            Stability Augmentation System.
        """

        self.sys = sas_sys