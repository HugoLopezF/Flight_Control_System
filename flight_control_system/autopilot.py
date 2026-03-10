from collections.abc import Mapping
import control
import numpy as np
from enum import Enum
from dataclasses import dataclass

from .types import Axis, InputChannel, OutputChannel
from .axis_metadata import AXIS_LABELS
from .sensor import Sensor

class AutopilotMode(str, Enum):
    """
    Autopilot modes.

    Attributes
    ----------
    THETA_HOLD : str
        Hold pitch angle.
    PHI_HOLD : str
        Hold roll angle.
    """

    THETA_HOLD = "theta_hold"
    PHI_HOLD = "phi_hold"
    # ALTITUDE_HOLD = "altitude_hold"

@dataclass(frozen=True)
class ControlledVarSpec:
    """
    Controlled variable definition as weighted sum of outputs or dynamics.

    Attributes
    ----------
    weights: Mapping[OutputChannel, float]
        Weights of each output variable.
    dynamics: control.TransferFunction | control.StateSpace | None
        Dynamics of output variables (default is None).
    """

    # z = W*y, optional dynamics for derived vars (e.g. altitude estimator)
    weights: Mapping[OutputChannel, float]
    dynamics: control.TransferFunction | control.StateSpace | None = None

@dataclass(frozen=True)
class PIDGains:
    """
    PID controller gains.

    Attributes
    ----------
    kp: float
        Proportional gain (default is 1.0).
    ki: float
        Integral gain (default is 0.0).
    kd: float
        Derivative gain (default is 0.0).
    n: float
        Derivative pole (default is 0.0).
    """

    kp: float = 1.0
    ki: float = 0.0
    kd: float = 0.0
    n: float = 0.0

@dataclass(frozen=True)
class AutopilotSpec:
    """
    Autopilot specifications.

    Attributes
    ----------
    axis: Axis
        Axis to study.
    controlled: ControlledVarSpec
        Controlled variable definition as weighted sum of outputs or dynamics.
    feedback: Mapping[InputChannel, tuple[OutputChannel, ...]]
        Feedback channels for each input channel.
    command_input: InputChannel
        Input channel controlled by autopilot.
    """

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
        """
        Obtain feedback channels to use in Autopilot.

        Parameters
        ----------
        spec: AutopilotSpec
            Autopilot specifications.

        Returns
        ----------
        tuple[OutputChannel, ...]
            Feedback channels.
        """

        seen: list[OutputChannel] = []
        for outs in spec.feedback.values():
            if isinstance(outs, OutputChannel):
                outs = (outs,)
            elif isinstance(outs, str):
                raise TypeError("feedback values must be OutputChannel or tuple[OutputChannel]")
            for ch in outs:
                if ch not in seen:
                    seen.append(ch)
        return tuple(seen)

    def _sensor_block(
        self,
        channels: tuple[OutputChannel, ...],
        sensors: Mapping[OutputChannel, Sensor],
    ) -> control.StateSpace:
        """
        Obtain sensors block to use in Autopilot.

        Parameters
        ----------
        axis: Axis
            Axis to study.
        sensors: Mapping[OutputChannel, Sensor]
            Sensors for each output channel.

        Returns
        ----------
        control.StateSpace
            Sensors block.
        """

        blocks = []
        for ch in channels:  # only feedback channels
            blocks.append(control.ss(sensors.get(ch, Sensor.ideal()).tf()))
        return control.append(*blocks)

    def _output_selector(
        self,
        axis: Axis,
        channels: tuple[OutputChannel, ...],
    ) -> control.StateSpace:
        """
        Create matrix to select outputs for feedback in Autopilot.

        Parameters
        ----------
        axis: Axis
            Axis to study.
        channels: tuple[OutputChannel, ...]
            Feedback channels.

        Returns
        ----------
        control.StateSpace
            Matrix to select feedback outputs.
        """

        ny = len(AXIS_LABELS[axis].state_channels)
        D = np.zeros((len(channels), ny), dtype=float)
        for r, ch in enumerate(channels):
            c = AXIS_LABELS[axis].state_channels.index(ch)
            D[r, c] = 1.0
        return control.ss([], [], [], D, 0)
    
    def _command_allocation(
        self,
        axis: Axis,
        command_input: InputChannel
    ) -> control.StateSpace:
        """
        Allocate commanded input channel in Autopilot.

        Parameters
        ----------
        axis: Axis
            Axis to study.
        command_input: InputChannel
            Command input channel.

        Returns
        ----------
        control.StateSpace
            Matrix to select commanded input.
        """
                
        in_channels = AXIS_LABELS[axis].input_channels
        D = np.zeros((len(in_channels), 1), dtype=float)
        D[in_channels.index(command_input), 0] = 1.0
        return control.ss([], [], [], D, 0)

    def _controlled_var_block(
        self,
        channels: tuple[OutputChannel, ...],
        cv: ControlledVarSpec,
    ) -> control.StateSpace:
        """
        Obtain controlled variable from weighted outputs or dynamics in Autopilot.

        Parameters
        ----------
        channels: tuple[OutputChannel, ...]
            Axis to study.
        cv: ControlledVarSpec
            Weights or dynamics or controlled variable.

        Returns
        ----------
        control.StateSpace
            Matrix with weights/dynamics to obtain controlled variable from outputs.
        """

        D = np.zeros((1, len(channels)), dtype=float)
        for ch, w in cv.weights.items():
            if ch not in channels:
                raise ValueError(f"{ch.value} not in controlled/feedback channels")
            D[0, channels.index(ch)] = float(w)
        W = control.ss([], [], [], D, 0)
        return W if cv.dynamics is None else control.series(control.ss(cv.dynamics), W)
        
    def _pid_tf(self, g: PIDGains) -> control.TransferFunction:
        """
        Obtain PID controller from gains.

        Parameters
        ----------
        g: PIDGains
            PID controller gains.

        Returns
        ----------
        control.TransferFunction
            PID controller transfer function.
        """

        P = control.tf([g.kp], [1], 0)
        I = control.tf([g.ki], [1, 0], 0) if g.ki else control.tf([0], [1], 0)
        D = control.tf([g.kd * g.n, 0], [1, g.n], 0) if g.kd else control.tf([0], [1], 0)
        return P + I + D

    def build_ap(
        self,
        axis: Axis,
        mode: AutopilotMode,
        sensors: Mapping[OutputChannel, Sensor] | None = None,
        pid_gains: PIDGains | None = None,
    ):
        """
        Construct Autopilot.

        Parameters
        ----------
        axis: Axis
            Axis to study.
        mode : AutopilotMode
            Autopilot mode.
        sensors: Mapping[OutputChannel, Sensor]
            Sensors for each output channel.
        pid_gains : PIDGains
            PID controller gains.

        Returns
        ----------
        tuple[control.StateSpace, control.StateSpace, np.ndarray]
            Autopilot, Controller input and PID matrix.
        """

        sensors = {} if sensors is None else sensors
        pid_gains = PIDGains() if pid_gains is None else pid_gains

        spec = self.SUPPORTED_AUTOPILOT[mode]
        if spec.axis is not axis:
            raise ValueError(f"{mode.value} is not valid for axis {axis.value}")

        G = control.ss(self.sys)  # ny x nu

        # y -> selected feedback channels -> sensor outputs
        fb_channels = self._feedback_channels(spec)
        Sel_fb = self._output_selector(axis, fb_channels)   # nfb x ny
        S_fb = self._sensor_block(fb_channels, sensors)     # nfb x nfb
        H_ch = control.series(S_fb, Sel_fb)                 # nfb x ny

        # collapse feedback channels into controlled measured variable z_m
        W_meas = self._controlled_var_block(fb_channels, spec.controlled)  # 1 x nfb
        H_ap = control.series(W_meas, H_ch)                                 # 1 x ny

        # e -> PID -> allocate to commanded input channel
        C_pid = control.ss(self._pid_tf(pid_gains))                          # 1 x 1
        B_cmd = self._command_allocation(axis, spec.command_input)          # nu x 1
        C_ap = control.series(B_cmd, C_pid)                                 # nu x 1

        # forward path from scalar command r to y
        P_ap = control.series(C_ap, G)                                     # ny x 1

        # Closed loop: e = r - H_ap*y
        ap_sys = control.feedback(P_ap, H_ap, sign=-1)                    # ny x 1

        # Controlled variable from true outputs
        full_channels = AXIS_LABELS[axis].state_channels
        W_out = self._controlled_var_block(full_channels, spec.controlled)  # 1 x ny
        cl_z = control.series(W_out, ap_sys)                              # 1 x 1

        return ap_sys, cl_z, C_ap