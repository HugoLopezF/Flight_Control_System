from enum import Enum
import numpy as np


class InputSignal(str, Enum):
    """
    Input signal type.

    Attributes
    ----------
    STEP : str
        Step input.
    SAT_RAMP : str
        Saturated ramp input.
    PULSE : str
        Pulse input.
    """

    STEP = "step"
    SAT_RAMP = "sat_ramp"
    PULSE = "pulse"


def step(t: np.ndarray, amp: float = 1.0) -> np.ndarray:
    """
    Define step input signal.

    Parameters
    ----------
    t : np.ndarray
        Time vector in seconds.
    amp : float, optional
        Step amplitude (default is 1.0).

    Returns
    -------
    np.ndarray
        Step input signal.
    """

    return amp * np.ones_like(t)


def pulse(t: np.ndarray, amp: float = 1.0, t_i: float = 0.0, t_end: float = 10.0) -> np.ndarray:
    """
    Define pulse input signal.

    u(t) = 0   for t < t_i
    u(t) = amp for t_i <= t <= t_end
    u(t) = 0   for t > t_end

    Parameters
    ----------
    t : np.ndarray
        Time vector in seconds.
    amp : float, optional
        Pulse amplitude (default is 1.0).
    t_i : float, optional
        Pulse start time in seconds (default is 0.0).
    t_end : float, optional
        Pulse end time in seconds (default is 10.0).

    Returns
    -------
    np.ndarray
        Pulse input signal.
    """

    if t_i < 0:
        raise ValueError("t_i must be >= 0")
    if t_end <= t_i:
        raise ValueError("t_end must be > t_i")

    u = np.zeros_like(t, dtype=float)
    u[(t >= t_i) & (t <= t_end)] = amp
    return u


def sat_ramp(t: np.ndarray, amp: float = 1.0, t_end: float = 10.0) -> np.ndarray:
    """
    Define saturated ramp input signal.

    u(t) = 0   for t = 0
    u(t) = amp / t_end * t for t < t_end
    u(t) = amp for t >= t_end

    Parameters
    ----------
    t : np.ndarray
        Time vector in seconds.
    amp : float, optional
        Saturated ramp amplitude (default is 1.0).
    t_end : float, optional
        Saturated ramp saturation time in seconds (default is 10.0).

    Returns
    -------
    np.ndarray
        Saturated ramp input signal.
    """
    
    if t_end <= 0:
        raise ValueError("t_end must be > 0")
    slope = amp / t_end
    return np.minimum(slope * t, amp)


def build_input(signal: InputSignal, t: np.ndarray, amp: float = 1.0, t_i: float = 0.0, t_end: float = 10.0) -> np.ndarray:
    """
    Calculate input signal.

    Parameters
    ----------
    signal : InputSignal
        Signal type.
    t : np.ndarray
        Time vector in seconds.
    amp : float, optional
        Input signal amplitude (default is 1.0).
    t_i : float, optional
        Pulse start time in seconds (default is 0.0).
    t_end : float, optional
        Input signal end time in seconds (default is 10.0).

    Returns
    -------
    np.ndarray
        Input signal.
    """

    if signal is InputSignal.STEP:
        return step(t, amp)
    if signal is InputSignal.SAT_RAMP:
        return sat_ramp(t, amp, t_end)
    if signal is InputSignal.PULSE:
        return pulse(t, amp=amp, t_i=t_i, t_end=t_end)
    raise ValueError(f"Unsupported signal type: {signal}")