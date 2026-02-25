# -*- coding: utf-8 -*-
"""
SAS
-------------
Stability augmentation system definition, calculations and plots.

"""
from .sas_design import compute_DL_gain
from flight_control_system.types import Axis
from flight_control_system.actuator import Actuator
from flight_control_system.sensor import Sensor
from flight_control_system.state_space import LinearizedSystem
from typing import Iterable, Mapping
from .types import ActuatorOrder
from collections.abc import Iterable
import control
import numpy as np

class SAS:
    SUPPORTED = {
        Axis.LONG: ("Kalpha", "Kq"),
        Axis.LATDIR: ("Kbeta", "Kr"),
    }

    def __init__(self, lin_sys: LinearizedSystem):
        self.sys = lin_sys

    def build_sas(self, axis: Axis, feedback_gains: Mapping[str, Iterable], actuator: Actuator = Actuator(order=ActuatorOrder.IDEAL), sensor: Sensor = Sensor(delay=0.0)):
        DL_gain = compute_DL_gain(lin_sys=self.sys, axis=axis, actuator=actuator)

        ax_sys = self.sys.get_sys(axis)
        ol_sys = control.series(actuator.tf(), ax_sys)

        Ka, Kb = self.SUPPORTED[axis]
        feedback = np.zeros(ax_sys.noutputs, ax_sys.ninputs)
        if axis is Axis.LONG:
            feedback[[1, ax_sys.n_outputs], 0] = []
        elif axis is Axis.LATDIR:
            feedback[[0, 2], ax_sys.n_inputs] = []
        feedback = control.series(sensor.tf(), feedback)
        cl_sys = control.feedback(ol_sys, feedback)