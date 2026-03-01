import control
from .types import FilterType

class Filter():
    def __init__(
        self, 
        filter_type: FilterType | str = FilterType.IDEAL,
        omega_b: float = 0.5, 
        gain: float = 1.0, 
    ):
        """
        Filter model.

        Parameters
        ----------
        filter_type : FilterType | str = FilterType.IDEAL
            Filter type.
        omega_b : float
            Actuator breakout frequency in rad/s.
        gain : float
            Actuator static gain.
        """
        try:
            self.filter_type = FilterType(filter_type)
        except ValueError as exc:
            allowed = ", ".join(repr(o.value) for o in FilterType)
            raise ValueError(f"filter_type must be one of: {allowed}") from exc

        if gain <= 0:
            raise ValueError("gain must be > 0")
        if self.filter_type is not FilterType.IDEAL and omega_b <= 0:
            raise ValueError("omega_b must be > 0")

        self.omega_b = omega_b
        self.gain = gain


    @classmethod
    def ideal(cls, gain: float = 1.0) -> "Filter":
        return cls(filter_type=FilterType.IDEAL, gain=gain)

    @classmethod
    def lowpass(cls, omega_b: float = 0.5, gain: float = 1.0) -> "Filter":
        return cls(filter_type=FilterType.LOWPASS, omega_b=omega_b, gain=gain)

    @classmethod
    def highpass(cls, omega_b: float = 0.5, gain: float = 1.0) -> "Filter":
        return cls(filter_type=FilterType.HIGHPASS, omega_b=omega_b, gain=gain)

    def tf(self) -> control.TransferFunction:
        if self.filter_type is FilterType.IDEAL:
            return control.tf([self.gain], [1], 0)
        if self.filter_type is FilterType.LOWPASS:
            return control.tf([self.gain], [1 / self.omega_b, 1], 0)
        if self.filter_type is FilterType.HIGHPASS:
            return control.tf([self.gain / self.omega_b, 0], [1 / self.omega_b, 1], 0)
        raise RuntimeError("Unsupported filter type")