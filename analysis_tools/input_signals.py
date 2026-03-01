from enum import Enum
import numpy as np


class InputSignal(str, Enum):
    STEP = "step"
    SAT_RAMP = "sat_ramp"
    PULSE = "pulse"


def step(t: np.ndarray, amp: float = 1.0) -> np.ndarray:
    """
    Define step input signal.

    :param t: Time vector.
    :type t: np.linspace
    :param amp: Step amplitude (nondimensional).
    :type amp: float

    :rtype: np.ndarray
    """
    return amp * np.ones_like(t)


def pulse(t: np.ndarray, amp: float = 1.0, t_end: float = 10.0) -> np.ndarray:
    """
    Define pulse input signal.

    u(t) = amp for t < t_end
    u(t) = 0   for t >= t_end

    :param t: Time vector.
    :type t: np.linspace
    :param amp: Step amplitude (nondimensional).
    :type amp: float
    :param t_end: Time instant at which the step stops.
    :type t_end: float

    :rtype: np.ndarray
    """
    if t_end <= 0:
        raise ValueError("t_end must be > 0")
    return np.where(t < t_end, amp, 0.0)


def sat_ramp(t: np.ndarray, amp: float = 1.0, t_end: float = 10.0) -> np.ndarray:
    """
    Define saturated ramp input signal.

    :param t: Time vector.
    :type t: np.linspace
    :param amp: Ramp amplitude (nondimensional).
    :type amp: float
    :param t_end: Time instant at which the ramp stops.
    :type t_end: float

    :rtype: np.ndarray
    """
    if t_end <= 0:
        raise ValueError("t_end must be > 0")
    slope = amp / t_end
    return np.minimum(slope * t, amp)


def build_input(signal: InputSignal, t: np.ndarray, amp: float = 1.0, t_end: float = 10.0) -> np.ndarray:
    """
    Calculate input signal.

    :param signal: Signal type.
    :type signal: InputSignal
    :param t: Time vector.
    :type t: np.linspace
    :param amp: Ramp amplitude (nondimensional).
    :type amp: float
    :param t_end: Time instant at which the ramp stops.
    :type t_end: float

    :rtype: np.ndarray
    """
    if signal is InputSignal.STEP:
        return step(t, amp)
    if signal is InputSignal.SAT_RAMP:
        return sat_ramp(t, amp, t_end)
    if signal is InputSignal.PULSE:
        return pulse(t, amp, t_end)
    raise ValueError(f"Unsupported signal type: {signal}")