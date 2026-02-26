from collections.abc import Mapping
import control
import numpy as np

from .sas_design import compute_DL_gain
from .types import Axis, InputChannel, OutputChannel, ActuatorOrder
from .axis_metadata import AXIS_LABELS
from .actuator import Actuator
from .sensor import Sensor
from .state_space import LinearizedSystem


class SAS:
    # Which outputs are fed back to which input channel
    SUPPORTED_FEEDBACK = {
        Axis.LONG: {
            InputChannel.ELEVATOR: (OutputChannel.ALPHA, OutputChannel.Q),
        },
        Axis.LATDIR: {
            InputChannel.RUDDER: (OutputChannel.BETA, OutputChannel.R),
        },
    }

    # Channel where direct-link gain is applied
    SUPPORTED_DL_GAIN = {
        Axis.LONG: InputChannel.ELEVATOR,
        Axis.LATDIR: InputChannel.AILERON,
    }

    def __init__(self, lin_sys: LinearizedSystem):
        self.sys = lin_sys

    @staticmethod
    def _out_idx(axis: Axis, channel: OutputChannel) -> int:
        return AXIS_LABELS[axis].state_channels.index(channel)

    @staticmethod
    def _inp_idx(axis: Axis, channel: InputChannel) -> int:
        return AXIS_LABELS[axis].input_channels.index(channel)

    def _actuator_block(self, axis: Axis, actuators: Mapping[InputChannel, Actuator]):
        blocks = []
        for ch in AXIS_LABELS[axis].input_channels:
            act = actuators.get(ch, Actuator(order=ActuatorOrder.IDEAL))
            blocks.append(act.tf())
        return control.append(*blocks)

    def _sensor_block(self, axis: Axis, sensors: Mapping[OutputChannel, Sensor]):
        blocks = []
        for ch in AXIS_LABELS[axis].state_channels:
            sens = sensors.get(ch, Sensor(delay=0.0))
            blocks.append(sens.tf())
        return control.append(*blocks)

    def _feedback_gain_matrix(
        self,
        axis: Axis,
        feedback_gains: Mapping[OutputChannel, float],
    ) -> np.ndarray:
        ny = len(AXIS_LABELS[axis].state_channels)
        nu = len(AXIS_LABELS[axis].input_channels)
        K = np.zeros((nu, ny), dtype=float)

        for in_ch, out_list in self.SUPPORTED_FEEDBACK[axis].items():
            i = self._inp_idx(axis, in_ch)
            for out_ch in out_list:
                j = self._out_idx(axis, out_ch)
                K[i, j] = float(feedback_gains.get(out_ch, 0.0))
        return K

    def build_sas(
        self,
        axis: Axis,
        feedback_gains: Mapping[OutputChannel, float],
        actuators: Mapping[InputChannel, Actuator],
        sensors: Mapping[OutputChannel, Sensor],
        desired_out: float = 1.0,
    ):
        # Open loop
        G = self.sys.get_sys(axis)
        A = self._actuator_block(axis, actuators)
        P = control.series(A, G)

        # Feedback sensors and gains
        S = self._sensor_block(axis, sensors)
        K = self._feedback_gain_matrix(axis, feedback_gains)
        K_sys = control.ss([], [], [], K)
        H = control.series(S, K_sys)

        # Closed-loop
        cl_nom = control.feedback(P, H, sign=-1)

        # Direct-link only on selected channel
        dl_in = self.SUPPORTED_DL_GAIN[axis]
        dl_act = actuators.get(dl_in, Actuator(order=ActuatorOrder.IDEAL))
        dl_gain = compute_DL_gain(
            lin_sys=self.sys,
            axis=axis,
            actuator=dl_act,
            desired_out=desired_out,
        )

        nu = len(AXIS_LABELS[axis].input_channels)
        F = np.eye(nu)
        i_dl = self._inp_idx(axis, dl_in)
        F[i_dl, i_dl] = dl_gain
        F_sys = control.ss([], [], [], F)

        sas_sys = control.series(F_sys, cl_nom)
        return sas_sys, dl_gain, K