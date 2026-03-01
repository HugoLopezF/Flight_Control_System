import control
from .types import SensorModel


class Sensor:
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
        Sensor model.

        Parameters
        ----------
        model : float
            Sensor model.
        delay : float
            Sensor delay in seconds.
        gain : float
            Sensor static gain.
        den_order : int
            Denominator order of Pade approximation (n).
        num_order : int | None
            Numerator degree (numdeg). If None, uses same order as denominator.
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
        return cls(
            model=SensorModel.PADE,
            gain=gain,
            delay=delay,
            num_order=num_order,
            den_order=den_order,
        )

    @classmethod
    def first_order(cls, tau: float = 0.0, gain: float = 1.0) -> "Sensor":
        return cls(model=SensorModel.FIRST_ORDER, gain=gain, tau=tau)
    
    @classmethod
    def ideal(cls, gain: float = 1.0) -> "Sensor":
        return cls(model=SensorModel.IDEAL, gain=gain)

    def tf(self) -> control.TransferFunction:
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