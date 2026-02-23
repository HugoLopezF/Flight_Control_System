import numpy as np
import control
from enum import Enum


class Sensor:
    def __init__(
        self, 
        delay: float = 0.01, 
        gain: float = 1.0,
        num_order: int = 1,
        den_order: int = 1,
    ):
        """
        Padé delay sensor model.

        Parameters
        ----------
        delay : float
            Sensor delay in seconds.
        gain : float
            Sensor static gain.
        den_order : int
            Denominator order of Padé approximation (n).
        num_degree : int | None
            Numerator degree (numdeg). If None, uses same order as denominator.
        """
        if delay < 0:
            raise ValueError("delay must be >= 0")
        if den_order < 0:
            raise ValueError("den_order must be >= 0")
        if num_order < 0:
            raise ValueError("num_degree must be >= 0")
        if num_order > den_order:
            raise ValueError("num_degree must be <= den_order")

        self.delay = delay
        self.gain = gain
        self.num_order = num_order
        self.den_order = den_order

    def tf(self) -> control.TransferFunction:
        """
        Define sensor transfer function.
        
        """
        if self.delay == 0:
            return control.tf([self.gain], [1])

        num, den = control.pade(self.delay, self.den_order, self.num_order)
        return self.gain * control.tf(num, den) 