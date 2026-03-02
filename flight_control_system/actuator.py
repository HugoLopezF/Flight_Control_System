import control
from .types import ActuatorOrder

class Actuator:
    """
    Actuator class.

    If no model is passed in, an ideal actuator model is selected.

    Parameters
    ----------
    order : ActuatorOrder | str
        Actuator order.
    omega_b : float
        Actuator breakout frequency in rad/s.
    damp : float
        Actuator damping.
    gain : float
        Actuator static gain.

    Methods
    ----------
    ideal()
        Construct ideal actuator instance.
    first_order()
        Construct first order actuator instance.
    second_order()
        Construct second order actuator instance.
    tf()
        Construct actuator transfer function.
    """
        
    def __init__(
        self,
        order: ActuatorOrder | str = ActuatorOrder.IDEAL,
        omega_b: float = 25.0, 
        damp: float = 0.0, 
        gain: float = 1.0, 
    ):
        """
        Actuator class.

        Parameters
        ----------
        order : ActuatorOrder | str
            Actuator order (default is ActuatorOrder.IDEAL).
        omega_b : float
            Actuator breakout frequency in rad/s (default is 25.0).
        damp : float
            Actuator damping (default is 0.0).
        gain : float
            Actuator static gain (default is 1.0).
        """

        try:
            self.order = ActuatorOrder(order)
        except ValueError as exc:
            allowed = ", ".join(repr(o.value) for o in ActuatorOrder)
            raise ValueError(f"order must be one of: {allowed}") from exc

        if gain <= 0:
            raise ValueError("gain must be > 0")
        if self.order is not ActuatorOrder.IDEAL and omega_b <= 0:
            raise ValueError("omega_b must be > 0")
        if damp < 0:
            raise ValueError("damp must be >= 0")

        self.omega_b = omega_b
        self.damp = damp
        self.gain = gain


    @classmethod
    def ideal(cls, gain: float = 1.0) -> "Actuator":
        """
        Construct ideal actuator instance.

        Parameters
        ----------
        gain : float
            Actuator static gain (default is 1.0).

        Returns
        -------
        "Actuator"
            Ideal actuator instance.
        """
        
        return cls(order=ActuatorOrder.IDEAL, gain=gain)

    @classmethod
    def first_order(cls, omega_b: float = 25.0, gain: float = 1.0) -> "Actuator":
        """
        Construct first order actuator instance.

        Parameters
        ----------
        omega_b : float
            Actuator breakout frequency in rad/s (default is 25.0).
        gain : float
            Actuator static gain (default is 1.0).

        Returns
        -------
        "Actuator"
            First order actuator instance.
        """
                
        return cls(order=ActuatorOrder.FIRST, omega_b=omega_b, gain=gain)

    @classmethod
    def second_order(
        cls,
        omega_b: float = 25.0,
        damp: float = 0.7,
        gain: float = 1.0,
    ) -> "Actuator":
        """
        Construct second order actuator instance.

        Parameters
        ----------
        omega_b : float
            Actuator breakout frequency in rad/s (default is 25.0).
        damp : float
            Actuator damping (default is 0.7).
        gain : float
            Actuator static gain (default is 1.0).

        Returns
        -------
        "Actuator"
            Second order actuator instance.
        """
                
        return cls(order=ActuatorOrder.SECOND, omega_b=omega_b, damp=damp, gain=gain)

    def tf(self) -> control.TransferFunction:
        """
        Construct actuator transfer function.

        Returns
        -------
        control.TransferFunction
            Actuator transfer function.

        Raises
        ------
        RuntimeError
            If incorrect actuator order passed in as a parameter.
        """
                
        if self.order is ActuatorOrder.IDEAL:
            return control.tf([self.gain], [1], 0)
        if self.order is ActuatorOrder.FIRST:
            return control.tf([self.gain], [1 / self.omega_b, 1], 0)
        if self.order is ActuatorOrder.SECOND:
            return control.tf(
                [self.gain * self.omega_b**2],
                [1, 2 * self.damp * self.omega_b, self.omega_b**2],
                0,
            )
        raise RuntimeError("Unsupported actuator order")