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

@dataclass(frozen=True)
class SASDesignPoint:
    """
    Stability Augmentation System design point in a feedback gain sweep.

    Attributes
    ----------
    axis : Axis
        Axis to study.
    feedback_gains : dict[OutputChannel, float]
        Feedback gains to sweep for each feedback channel.
    sas_sys : control.StateSpace
        Stability Augmentation System.
    dl_gain : float
        Direct-Link gain.
    k_matrix : np.ndarray
        Feedback gain matrix.
    poles : np.ndarray
        System poles.
    """

    axis: Axis
    feedback_gains: dict[OutputChannel, float]
    sas_sys: control.StateSpace
    dl_gain: float
    k_matrix: np.ndarray

    @property
    def poles(self) -> np.ndarray:
        """
        Obtain system poles.

        Returns
        -------
        np.ndarray
            System poles.
        """

        return np.asarray(control.poles(self.sas_sys), dtype=complex)

def build_1dof_roll(lin_sys: LinearizedSystem) -> control.TransferFunction:
    """
    Builds 1 degree of freedom roll model.

    Parameters
    ----------
    lin_sys : LinearizedSystem
        Aircraft linearized system to study.

    Returns
    -------
    control.TransferFunction
        1DOF roll model transfer function.
    """

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
    """
    Compute Direct-Link gain for an aircraft linearized system.

    Parameters
    ----------
    lin_sys : LinearizedSystem
        Aircraft linearized system to study.
    axis : Axis
        Axis to study.
    actuator : Actuator, optional
        Actuator for this input channel (default is Actuator() i.e. ideal).
    desired_out : float, optional
        Desired steady-state value of the output variable for a
        step input of 1deg (default is 1.0).

    Returns
    -------
    float
        Direct-Link gain.
    """

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

def _feedback_channels(axis: Axis) -> tuple[OutputChannel, ...]:
    """
    Obtain feedback channels for the selected axis.

    Parameters
    ----------
    axis : Axis
        Axis to study.

    Returns
    -------
    tuple[OutputChannel, ...]
        Feedback channels.
    """

    # Local import avoids circular import
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
) -> dict[dict[Mapping[OutputChannel, float]]]:
    """
    Obtain feedback gain combinations.

    If no base_feedback_gains are provided, feedback is set to zero.

    Parameters
    ----------
    axis : Axis
        Axis to study.
    gain_values : Mapping[OutputChannel, Iterable[float]]
        Feedback gain values of each feedback channel.
    base_feedback_gains : Mapping[OutputChannel, float] | None, optional
        Baseline feedback gains (default is None).

    Returns
    -------
    dict[dict[Mapping[OutputChannel, float]]]
        Feedback gain combinations.
    """

    base = dict(base_feedback_gains or {})
    chans = _feedback_channels(axis)

    val_lists: list[list[float]] = []
    for ch in chans:
        vals = [float(v) for v in gain_values.get(ch, [])]
        # If not swept, keep fixed at base gain (default 0.0)
        val_lists.append(vals if vals else [float(base.get(ch, 0.0))])

    for comb in product(*val_lists):
        g = dict(base)
        g.update({ch: v for ch, v in zip(chans, comb)})
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
    """
    Sweep feedback gain combinations.

    If a component is not filled in, assume it is ideal type.

    Parameters
    ----------
    axis : Axis
        Axis to study.
    gain_values : Mapping[OutputChannel, Iterable[float]]
        Feedback gain values of each feedback channel.
    base_feedback_gains : Mapping[OutputChannel, float] | None, optional
        Baseline feedback gains (default is None).
    actuators : Mapping[InputChannel, Actuator] | None = None, optional
        Stability Augmentation System actuators (default is None).
    sensors : Mapping[OutputChannel, Sensor] | None, optional
        Stability Augmentation System sensors (default is None).
    filters : Mapping[OutputChannel, Filter] | None, optional
        Stability Augmentation System filters (default is None).
    desired_out : float, optional
        Desired steady-state value of the output variable for a
        step input of 1deg (default is 1.0).

    Returns
    -------
    list[SASDesignPoint]
        Stability Augmentation System design points for this
        feedback gain sweep.
    """

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