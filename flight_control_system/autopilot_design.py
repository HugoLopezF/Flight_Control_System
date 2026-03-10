from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass
from itertools import product
import control
import numpy as np

from .types import Axis, OutputChannel
from .sensor import Sensor
from .autopilot import AP, PIDGains, AutopilotMode


@dataclass(frozen=True)
class APDesignPoint:
    """
    Autopilot design point in a PID controller gain sweep.

    Attributes
    ----------
    axis : Axis
        Axis to study.
    mode: AutopilotMode
        Autopilot mode.
    pid_gains: PIDGains
        PID controller gains.
    ap_sys: control.StateSpace
        Autopilot.
    cl_z: control.StateSpace
        Controller input.
    C: control.StateSpace 
        PID controller matrix.
    poles : np.ndarray
        System poles.
    """

    axis: Axis
    mode: AutopilotMode
    pid_gains: PIDGains
    ap_sys: control.StateSpace    # r -> full outputs
    cl_z: control.StateSpace      # r -> controlled variable
    C: control.StateSpace         # error -> plant inputs

    @property
    def poles(self) -> np.ndarray:
        return np.asarray(control.poles(self.ap_sys), dtype=complex)


def iter_pid_gain_sets(
    gain_values: Mapping[str, Iterable[float]],
    base_pid_gains: PIDGains | None = None,
) -> Iterator[PIDGains]:
    """
    Obtain PID gain combinations.

    If no base_pid_gains are provided, no PID is set (PID = 1).

    Parameters
    ----------
    gain_values: Mapping[str, Iterable[float]]
        PID gain values.
    base_pid_gains: PIDGains | None, optional
        Baseline PID gains (default is None).

    Returns
    -------
    Iterator[PIDGains]
        PID gains combinations.
    """

    base = PIDGains() if base_pid_gains is None else base_pid_gains

    kp_vals = [float(v) for v in gain_values.get("kp", [])] or [base.kp]
    ki_vals = [float(v) for v in gain_values.get("ki", [])] or [base.ki]
    kd_vals = [float(v) for v in gain_values.get("kd", [])] or [base.kd]

    for kp, ki, kd in product(kp_vals, ki_vals, kd_vals):
        yield PIDGains(kp=kp, ki=ki, kd=kd)


def sweep_pid_gains(
    sas_sys: control.StateSpace,
    axis: Axis,
    mode: AutopilotMode,
    gain_values: Mapping[str, Iterable[float]],
    *,
    base_pid_gains: PIDGains | None = None,
    sensors: Mapping[OutputChannel, Sensor] | None = None,
) -> list[APDesignPoint]:
    """
    Sweep PID gain combinations.

    If a component is not filled in, assume it is ideal type.

    Parameters
    ----------
    axis : Axis
        Axis to study.
    mode: AutopilotMode
        Autopilot mode.
    gain_values : Mapping[OutputChannel, Iterable[float]]
        PID gain values.
    base_pid_gains: PIDGains | None = None, optional
        Baseline PID gains (default is None).
    sensors : Mapping[OutputChannel, Sensor] | None, optional
        Autopilot sensors (default is None).

    Returns
    -------
    list[APDesignPoint]
        Autopilot design points for this
        PID gain sweep.
    """

    ap = AP(sas_sys)
    sens_cfg = dict(sensors or {})

    points: list[APDesignPoint] = []
    for pid in iter_pid_gain_sets(gain_values, base_pid_gains):
        ap_sys, cl_z, C = ap.build_ap(
            axis=axis,
            mode=mode,
            sensors=sens_cfg,
            pid_gains=pid,
        )
        points.append(
            APDesignPoint(
                axis=axis,
                mode=mode,
                pid_gains=pid,
                ap_sys=ap_sys,
                cl_z=cl_z,
                C=C,
            )
        )
    return points