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
    axis: Axis
    mode: AutopilotMode
    pid_gains: PIDGains
    cl_z: control.StateSpace      # r -> controlled variable
    cl_y_cmd: control.StateSpace  # r -> full outputs
    c_ap: control.StateSpace      # error -> plant inputs

    @property
    def poles(self) -> np.ndarray:
        return np.asarray(control.poles(self.cl_y_cmd), dtype=complex)


def iter_pid_gain_sets(
    gain_values: Mapping[str, Iterable[float]],
    base_pid_gains: PIDGains | None = None,
) -> Iterator[PIDGains]:
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
    ap = AP(sas_sys)
    sens_cfg = dict(sensors or {})

    points: list[APDesignPoint] = []
    for pid in iter_pid_gain_sets(gain_values, base_pid_gains):
        cl_z, cl_y_cmd, c_ap = ap.build_ap(
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
                cl_z=cl_z,
                cl_y_cmd=cl_y_cmd,
                c_ap=c_ap,
            )
        )
    return points