from collections.abc import Mapping
import control
import numpy as np
from enum import Enum
from dataclasses import dataclass

from .types import Axis, InputChannel, OutputChannel
from .axis_metadata import AXIS_LABELS
from .sensor import Sensor

class AutopilotMode(str, Enum):
    THETA_HOLD = "theta_hold"
    PHI_HOLD = "phi_hold"
    # ALTITUDE_HOLD = "altitude_hold"

@dataclass(frozen=True)
class ControlledVarSpec:
    # z = W*y, optional dynamics for derived vars (e.g. altitude estimator)
    weights: Mapping[OutputChannel, float]
    dynamics: control.TransferFunction | control.StateSpace | None = None

@dataclass(frozen=True)
class PIDGains:
    kp: float = 0.0
    ki: float = 0.0
    kd: float = 0.0
    n: float = 20.0  # derivative filter pole

@dataclass(frozen=True)
class AutopilotSpec:
    axis: Axis
    controlled: ControlledVarSpec
    feedback: Mapping[InputChannel, tuple[OutputChannel, ...]]
    command_input: InputChannel


class AP:
    """
    Autopilot class.

    Attributes
    ----------

    Methods
    ----------
    """

    SUPPORTED_AUTOPILOT: dict[AutopilotMode, AutopilotSpec] = {
        AutopilotMode.THETA_HOLD: AutopilotSpec(
            axis=Axis.LONG,
            controlled=ControlledVarSpec(weights={OutputChannel.THETA: 1.0}),
            feedback={InputChannel.ELEVATOR: (OutputChannel.THETA)},
            command_input=InputChannel.ELEVATOR,
        ),
        AutopilotMode.PHI_HOLD: AutopilotSpec(
            axis=Axis.LATDIR,
            controlled=ControlledVarSpec(weights={OutputChannel.PHI: 1.0}),
            feedback={InputChannel.AILERON: (OutputChannel.PHI)},
            command_input=InputChannel.AILERON,
        ),
    }

    def __init__(self, sas_sys: control.StateSpace):
        """
        Autopilot class.

        Parameters
        ----------
        sas_sys : control.StateSpace
            Stability Augmentation System.
        """

        self.sys = sas_sys

    def _feedback_channels(self, spec: AutopilotSpec) -> tuple[OutputChannel, ...]:
        seen: list[OutputChannel] = []
        for outs in spec.feedback.values():
            for ch in outs:
                if ch not in seen:
                    seen.append(ch)
        return tuple(seen)

    def _sensor_block(
        self,
        channels: tuple[OutputChannel, ...],
        sensors: Mapping[OutputChannel, Sensor],
    ) -> control.StateSpace:
        blocks = []
        for ch in channels:  # only feedback channels
            blocks.append(control.ss(sensors.get(ch, Sensor.ideal()).tf()))
        return control.append(*blocks)

    def _output_selector(self, axis: Axis, channels: tuple[OutputChannel, ...]) -> control.StateSpace:
        ny = len(AXIS_LABELS[axis].state_channels)
        D = np.zeros((len(channels), ny), dtype=float)
        for r, ch in enumerate(channels):
            c = AXIS_LABELS[axis].state_channels.index(ch)
            D[r, c] = 1.0
        return control.ss([], [], [], D, 0)

    def _controlled_var_block(self, axis: Axis, cv: ControlledVarSpec) -> control.StateSpace:
        ny = len(AXIS_LABELS[axis].state_channels)
        W = np.zeros((1, ny), dtype=float)
        for ch, w in cv.weights.items():
            j = AXIS_LABELS[axis].state_channels.index(ch)
            W[0, j] = float(w)
        W_sys = control.ss([], [], [], W, 0)
        return W_sys if cv.dynamics is None else control.series(control.ss(cv.dynamics), W_sys)
        
    def _pid_tf(self, g: PIDGains) -> control.TransferFunction:
        P = control.tf([g.kp], [1], 0)
        I = control.tf([g.ki], [1, 0], 0) if g.ki else control.tf([0], [1], 0)
        D = control.tf([g.kd * g.n, 0], [1, g.n], 0) if g.kd else control.tf([0], [1], 0)
        return P + I + D

    def _controller_matrix(
        self,
        axis: Axis,
        spec: AutopilotSpec,
        fb_channels: tuple[OutputChannel, ...],
        pid_gains: Mapping[tuple[InputChannel, OutputChannel], PIDGains],
    ):
        in_channels = AXIS_LABELS[axis].input_channels
        zero = control.tf([0], [1], 0)
        C = [[zero for _ in fb_channels] for _ in in_channels]

        for in_ch, outs in spec.feedback.items():
            i = in_channels.index(in_ch)
            for out_ch in outs:
                j = fb_channels.index(out_ch)
                C[i][j] = self._pid_tf(pid_gains.get((in_ch, out_ch), PIDGains()))
        return control.tf(C)

    def build_ap(
        self,
        axis: Axis,
        mode: AutopilotMode,
        sensors: Mapping[OutputChannel, Sensor] | None = None,
        pid_gains: Mapping[tuple[InputChannel, OutputChannel], PIDGains] | None = None,
    ):
        sensors = {} if sensors is None else sensors
        pid_gains = {} if pid_gains is None else pid_gains

        spec = self.SUPPORTED_AUTOPILOT[mode]
        if spec.axis is not axis:
            raise ValueError(f"{mode.value} is not valid for axis {axis.value}")

        G = control.ss(self.sys)
        fb_channels = self._feedback_channels(spec)

        Sel_fb = self._output_selector(axis, fb_channels)
        S_fb = self._sensor_block(fb_channels, sensors)
        H = control.series(S_fb, Sel_fb)

        C = self._controller_matrix(axis, spec, fb_channels, pid_gains)
        L = control.series(C, H)
        cl_y = control.feedback(G, L, sign=-1)

        # optional scalar command channel
        nu = len(AXIS_LABELS[axis].input_channels)
        E = np.zeros((nu, 1))
        cmd_idx = AXIS_LABELS[axis].input_channels.index(spec.command_input)
        E[cmd_idx, 0] = 1.0
        cl_y_cmd = control.series(cl_y, control.ss([], [], [], E, 0))

        W = self._controlled_var_block(axis, spec.controlled)
        cl_z = control.series(W, cl_y_cmd)

        return cl_z, cl_y_cmd, C