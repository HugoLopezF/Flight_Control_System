import control
from .types import FilterType

class Filter():
    """
    Filter class.

    If no model is passed in, an ideal filter model is selected.

    Attributes
    ----------
    filter_type : FilterType | str
        Filter model.
    omega_b : float
        Filter breakout frequency in rad/s.
    gain : float
        Filter static gain.

    Methods
    ----------
    ideal()
        Construct ideal filter instance.
    lowpass()
        Construct lowpass filter instance.
    highpass()
        Construct highpass filter instance.
    tf()
        Construct filter transfer function.
    """
        
    def __init__(
        self, 
        filter_type: FilterType | str = FilterType.IDEAL,
        omega_b: float = 0.5, 
        gain: float = 1.0, 
    ):
        """
        Filter class.

        Parameters
        ----------
        filter_type : FilterType | str, optional
            Filter type (default is FilterType.IDEAL).
        omega_b : float, optional
            Filter breakout frequency in rad/s (default is 0.5).
        gain : float, optional
            Filter static gain (default is 1.0).
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
        """
        Construct ideal filter instance.

        Parameters
        ----------
        gain : float, optional
            Filter static gain (default is 1.0).

        Returns
        -------
        "Filter"
            Ideal filter instance.
        """

        return cls(filter_type=FilterType.IDEAL, gain=gain)

    @classmethod
    def lowpass(cls, omega_b: float = 0.5, gain: float = 1.0) -> "Filter":
        """
        Construct lowpass filter instance.

        Parameters
        ----------
        omega_b : float, optional
            Filter breakout frequency in rad/s (default is 0.5).
        gain : float, optional
            Filter static gain (default is 1.0).

        Returns
        -------
        "Filter"
            Lowpass filter instance.
        """

        return cls(filter_type=FilterType.LOWPASS, omega_b=omega_b, gain=gain)

    @classmethod
    def highpass(cls, omega_b: float = 0.5, gain: float = 1.0) -> "Filter":
        """
        Construct highpass filter instance.

        Parameters
        ----------
        omega_b : float, optional
            Filter breakout frequency in rad/s (default is 0.5).
        gain : float, optional
            Filter static gain (default is 1.0).

        Returns
        -------
        "Filter"
            Highpass filter instance.
        """
                
        return cls(filter_type=FilterType.HIGHPASS, omega_b=omega_b, gain=gain)

    def tf(self) -> control.TransferFunction:
        """
        Construct filter transfer function.

        Returns
        -------
        control.TransferFunction
            Filter transfer function.

        Raises
        ------
        RuntimeError
            If incorrect filter model passed in as a parameter.
        """
                
        if self.filter_type is FilterType.IDEAL:
            return control.tf([self.gain], [1], 0)
        if self.filter_type is FilterType.LOWPASS:
            return control.tf([self.gain], [1 / self.omega_b, 1], 0)
        if self.filter_type is FilterType.HIGHPASS:
            return control.tf([self.gain / self.omega_b, 0], [1 / self.omega_b, 1], 0)
        raise RuntimeError("Unsupported filter type")