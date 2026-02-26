import control
from .types import ActuatorOrder

class Actuator:
    def __init__(
        self, 
        omega_b: float = 25.0, 
        damp: float = 0.0, 
        gain: float = 1.0, 
        order: ActuatorOrder | str = ActuatorOrder.IDEAL,
    ):
        """
        Actuator first and second order model.

        Parameters
        ----------
        omega_b : float
            Actuator breakout frequency in rad/s.
        damp : float
            Actuator damping.
        gain : float
            Actuator static gain.
        order : ActuatorOrder | str
            Actuator order.
        """
        if omega_b <= 0:
            raise ValueError("omega_b must be > 0")
        if damp < 0:
            raise ValueError("damp must be >= 0")
        try:
            self.order = ActuatorOrder(order)
        except ValueError as exc:
            raise ValueError("order must be 'ideal', 'first' or 'second'") from exc

        self.omega_b = omega_b
        self.damp = damp
        self.gain = gain


    def tf(self) -> control.TransferFunction:
        """
        Define actuator transfer function
        
        """
        if self.order is ActuatorOrder.IDEAL:
            return control.tf(
                [self.gain], 
                [1],
            )
        if self.order is ActuatorOrder.FIRST:
            return control.tf(
                [self.gain], 
                [1 / self.omega_b, 1],
            )
        if self.order is ActuatorOrder.SECOND:
            return control.tf(
                [self.gain * self.omega_b ** 2], 
                [1, 2 * self.damp * self.omega_b, self.omega_b ** 2],
            )
        raise RuntimeError("Unsupported actuator order")