from .state_space import LinearizedSystem
from .types import Axis
from .actuator import Actuator
import control
import numpy as np

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