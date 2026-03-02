import control
from .types import SensorModel


class Sensor:
    """
    Sensor class.

    If no model is passed in, an ideal sensor model is selected.

    Attributes
    ----------
    model : SensorModel | str
        Sensor model.
    delay : float
        Sensor delay in seconds.
    gain : float
        Sensor static gain.
    num_order : int
        Numerator degree (numdeg). If None, uses same order as denominator.
    den_order : int
        Denominator order of Pade approximation (n).
    tau : float
        Sensor time constant if SensorModel is first order.

    Methods
    ----------
    pade()
        Construct pade sensor instance.
    first_order()
        Construct first order sensor instance.
    ideal()
        Construct ideal sensor instance.
    tf()
        Construct sensor transfer function.
    """

    def __init__(
        self,
        model: SensorModel | str = SensorModel.IDEAL,
        delay: float = 0.01, 
        gain: float = 1.0,
        num_order: int = 1,
        den_order: int = 1,
        tau: float = 0.0,
    ):
        """
        Sensor class.

        If no model is passed in, an ideal sensor model is selected.

        Parameters
        ----------
        model : SensorModel | str
            Sensor model (default is SensorModel.IDEAL).
        delay : float, optional
            Sensor delay in seconds (default is 0.01).
        gain : float, optional
            Sensor static gain (default is 1.0).
        num_order : int, optional
            Numerator degree (numdeg). If None, uses same order as denominator (default is 1).
        den_order : int, optional
            Denominator order of Pade approximation (n) (default is 1).
        tau : float
            Sensor time constant if SensorModel is first order.
        """

        try:
            self.model = SensorModel(model)
        except ValueError as exc:
            allowed = ", ".join(repr(m.value) for m in SensorModel)
            raise ValueError(f"model must be one of: {allowed}") from exc

        if gain <= 0:
            raise ValueError("gain must be > 0")

        self.gain = gain
        self.delay = delay
        self.num_order = num_order
        self.den_order = den_order
        self.tau = tau

        if self.model is SensorModel.PADE:
            if delay < 0:
                raise ValueError("delay must be >= 0 for PADE")
            if den_order < 0 or num_order < 0:
                raise ValueError("num_order and den_order must be >= 0 for PADE")
            if num_order > den_order:
                raise ValueError("num_order must be <= den_order for PADE")
        elif self.model is SensorModel.FIRST_ORDER:
            if tau < 0:
                raise ValueError("tau must be >= 0 for FIRST_ORDER")
            
    @classmethod
    def pade(
        cls,
        delay: float = 0.01,
        gain: float = 1.0,
        num_order: int = 1,
        den_order: int = 1,
    ) -> "Sensor":
        """
        Construct pade sensor instance.

        Parameters
        ----------
        delay : float, optional
            Sensor delay in seconds (default is 0.01).
        gain : float, optional
            Sensor static gain (default is 1.0).
        num_order : int, optional
            Numerator degree (numdeg). If None, uses same order as denominator (default is 1).
        den_order : int, optional
            Denominator order of Pade approximation (n) (default is 1).

        Returns
        -------
        "Sensor"
            Pade sensor instance.
        """

        return cls(
            model=SensorModel.PADE,
            gain=gain,
            delay=delay,
            num_order=num_order,
            den_order=den_order,
        )

    @classmethod
    def first_order(cls, tau: float = 0.0, gain: float = 1.0) -> "Sensor":
        """
        Construct first order sensor instance.

        Parameters
        ----------
        tau : float, optional
            Sensor time constant in seconds (default is 0.0).
        gain : float, optional
            Sensor static gain (default is 1.0).

        Returns
        -------
        "Sensor"
            First order sensor instance.
        """

        return cls(model=SensorModel.FIRST_ORDER, gain=gain, tau=tau)
    
    @classmethod
    def ideal(cls, gain: float = 1.0) -> "Sensor":
        """
        Construct ideal sensor instance.

        Parameters
        ----------
        gain : float, optional
            Sensor static gain (default is 1.0).

        Returns
        -------
        "Sensor"
            Ideal sensor instance.
        """

        return cls(model=SensorModel.IDEAL, gain=gain)

    def tf(self) -> control.TransferFunction:
        """
        Construct sensor transfer function.

        Returns
        -------
        control.TransferFunction
            Sensor transfer function.

        Raises
        ------
        RuntimeError
            If incorrect sensor model passed in as a parameter.
        """

        if self.model is SensorModel.PADE:
            if self.delay == 0:
                return control.tf([self.gain], [1])
            num, den = control.pade(self.delay, self.den_order, self.num_order)
            return self.gain * control.tf(num, den, 0)

        if self.model is SensorModel.FIRST_ORDER:
            if self.tau == 0:
                return control.tf([self.gain], [1], 0)
            return control.tf([self.gain], [self.tau, 1], 0)
        
        if self.model is SensorModel.IDEAL:
            return control.tf([self.gain], [1], 0)

        raise RuntimeError("Unsupported sensor model")