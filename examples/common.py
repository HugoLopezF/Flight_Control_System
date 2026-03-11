import numpy as np

from flight_dynamics.aircraft import Aircraft
from flight_control_system.state_space import LinearizedSystem
from flight_control_system.types import Axis, InputChannel, OutputChannel
from flight_control_system.actuator import Actuator
from flight_control_system.sensor import Sensor
from flight_control_system.filter import Filter
from flight_control_system.autopilot import PIDGains

MODEL = "Learjet_24"


def initialize(model: str = MODEL) -> tuple[Aircraft, LinearizedSystem]:
    aircraft = Aircraft(model)
    lin_sys = LinearizedSystem(aircraft)
    return aircraft, lin_sys


# Aircraft stability coefficient sweep factors
long_factors = {
    "Cm_alpha": [-0.5, -0.2, 0.0, 0.5, 1.0, 2.0, 3.0, 4.0],
    "Cm_q": [-0.5, -0.2, 0.0, 0.5, 1.0, 2.0, 3.0, 4.0],
}

latdir_factors = {
    "Cn_beta": [0.0, 0.5, 1.0, 2.0, 3.0],
    "Cn_r": [0.0, 0.5, 1.0, 2.0, 3.0, 5.0],
}


# SAS components (defaults)
long_sas_actuators = {
    InputChannel.ELEVATOR: Actuator.first_order(omega_b=10.0),
}
long_sas_sensors = {
    OutputChannel.ALPHA: Sensor(),
    OutputChannel.Q: Sensor(),
}

latdir_sas_actuators = {
    InputChannel.AILERON: Actuator.first_order(omega_b=10.0),
    InputChannel.RUDDER: Actuator.first_order(omega_b=10.0),
}
latdir_sas_sensors = {
    OutputChannel.BETA: Sensor(),
    OutputChannel.R: Sensor(),
}
latdir_sas_filters = {
    OutputChannel.R: Filter.highpass(omega_b=2.0),  # Washout filter
}


# SAS fixed gains (example values)
long_feedback_gains = {
    OutputChannel.ALPHA: -0.75,
    OutputChannel.Q: -0.2,
}
latdir_feedback_gains = {
    OutputChannel.BETA: -2.0,
    OutputChannel.R: 0.55,
}


# Autopilot PID sweep defaults
long_pid_values = {
    "kp": np.linspace(-0.5, 0.1, 4),
    "ki": np.linspace(-0.5, 0.1, 4),
}

latdir_pid_values = {
    "kp": np.linspace(0.1, 1, 8),
}

# Autopilot fixed PID gains (from main.py)
long_ap_pid_gains = PIDGains(kp=-0.1, ki=-0.1)
latdir_ap_pid_gains = PIDGains(kp=1.0)


def compute_long_sas_gain_sweep(aircraft: Aircraft, factors: dict[str, list[float]]) -> dict[OutputChannel, list[float]]:
    Cm_alpha = aircraft.stab_coeffs.long.Cm_alpha
    Cm_q = aircraft.stab_coeffs.long.Cm_q
    Cm_delta_e = aircraft.stab_coeffs.long.Cm_delta_e
    c = aircraft.geom.c
    u_s = aircraft.flight_cond.u_s

    return {
        OutputChannel.ALPHA: [-(F_alpha - 1.0) * Cm_alpha / Cm_delta_e for F_alpha in factors["Cm_alpha"]],
        OutputChannel.Q: [-(F_q - 1.0) * Cm_q / Cm_delta_e * c / (2.0 * u_s) for F_q in factors["Cm_q"]],
    }


def compute_latdir_sas_gain_sweep(
    aircraft: Aircraft,
    factors: dict[str, list[float]],
) -> dict[OutputChannel, list[float]]:
    Cn_beta = aircraft.stab_coeffs.latdir.Cn_beta
    Cn_r = aircraft.stab_coeffs.latdir.Cn_r
    Cn_delta_r = aircraft.stab_coeffs.latdir.Cn_delta_r
    b = aircraft.geom.b
    u_s = aircraft.flight_cond.u_s

    return {
        OutputChannel.BETA: [-(F_beta - 1.0) * Cn_beta / Cn_delta_r for F_beta in factors["Cn_beta"]],
        OutputChannel.R: [-(F_r - 1.0) * Cn_r / Cn_delta_r * b / (2.0 * u_s) for F_r in factors["Cn_r"]],
    }
