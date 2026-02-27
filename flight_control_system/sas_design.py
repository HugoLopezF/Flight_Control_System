from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from itertools import product
import control
import numpy as np

from .state_space import LinearizedSystem
from .types import Axis, OutputChannel, InputChannel
from .actuator import Actuator
from .sensor import Sensor
from .filter import Filter

def build_1dof_roll(lin_sys: LinearizedSystem):
    ac = lin_sys.aircraft
    rho = ac.flight_cond.rho
    u_s = ac.flight_cond.u_s
    S = ac.geom.S
    b = ac.geom.b
    I_xx = ac.mass_prop.I_xx
    Cl_delta_a = ac.stab_coeffs.latdir.Cl_delta_a
    Cl_p = ac.stab_coeffs.latdir.Cl_p

    K = - 2 * u_s * Cl_delta_a / (b * Cl_p)
    tau = - 4 * I_xx / (rho * S * b**2 * u_s * Cl_p)

    return control.tf([K], [tau, 1])

def compute_DL_gain(
    lin_sys: LinearizedSystem,
    axis: Axis,
    actuator: Actuator = Actuator(),
    desired_out: float = 1.0,
) -> float:
    if axis is Axis.LONG:
        plant = lin_sys.get_sys(axis)[2, 0]  # theta / delta_e channel from full longitudinal model
    elif axis is Axis.LATDIR:
        plant = build_1dof_roll(lin_sys) # p / delta_a approximation
    else:
        raise ValueError(f"Unsupported axis: {axis}")

    ol_sys = control.series(actuator.tf(), plant)
    dcgain = float(np.real_if_close(control.dcgain(ol_sys)))

    if not np.isfinite(dcgain) or np.isclose(dcgain, 0.0):
        raise ValueError(f"Invalid dcgain for DL gain computation: {dcgain}")

    return float(desired_out / abs(dcgain))

@dataclass(frozen=True)
class SASDesignPoint:
    axis: Axis
    feedback_gains: dict[OutputChannel, float]
    sas_sys: control.StateSpace
    dl_gain: float
    k_matrix: np.ndarray

    @property
    def poles(self) -> np.ndarray:
        return np.asarray(control.poles(self.sas_sys), dtype=complex)


def _feedback_channels(axis: Axis) -> tuple[OutputChannel, ...]:
    # Local import avoids circular import (sas.py imports compute_DL_gain from here)
    from .sas import SAS
    seen: list[OutputChannel] = []
    for outputs in SAS.SUPPORTED_FEEDBACK[axis].values():
        for ch in outputs:
            if ch not in seen:
                seen.append(ch)
    return tuple(seen)


def iter_feedback_gain_sets(
    axis: Axis,
    gain_values: Mapping[OutputChannel, Iterable[float]],
    base_feedback_gains: Mapping[OutputChannel, float] | None = None,
):
    base = dict(base_feedback_gains or {})
    chans = _feedback_channels(axis)

    val_lists: list[list[float]] = []
    for ch in chans:
        vals = [float(v) for v in gain_values.get(ch, [])]
        # If not swept, keep fixed at base gain (default 0.0)
        val_lists.append(vals if vals else [float(base.get(ch, 0.0))])

    for combo in product(*val_lists):
        g = dict(base)
        g.update({ch: v for ch, v in zip(chans, combo)})
        yield g


def sweep_feedback_gains(
    lin_sys: LinearizedSystem,
    axis: Axis,
    gain_values: Mapping[OutputChannel, Iterable[float]],
    *,
    base_feedback_gains: Mapping[OutputChannel, float] | None = None,
    actuators: Mapping[InputChannel, Actuator] | None = None,
    sensors: Mapping[OutputChannel, Sensor] | None = None,
    filters: Mapping[OutputChannel, Filter] | None = None,
    desired_out: float = 1.0,
) -> list[SASDesignPoint]:
    from .sas import SAS

    sas_builder = SAS(lin_sys)

    act_cfg = dict(actuators or {})
    sens_cfg = dict(sensors or {})
    filt_cfg = dict(filters or {})

    points: list[SASDesignPoint] = []
    for gains in iter_feedback_gain_sets(axis, gain_values, base_feedback_gains):
        sas_sys, dl_gain, k_matrix = sas_builder.build_sas(
            axis=axis,
            feedback_gains=gains,
            actuators=act_cfg,
            sensors=sens_cfg,
            filters=filt_cfg,
            desired_out=desired_out,
        )
        points.append(
            SASDesignPoint(
                axis=axis,
                feedback_gains=dict(gains),
                sas_sys=sas_sys,
                dl_gain=dl_gain,
                k_matrix=k_matrix,
            )
        )
    return points