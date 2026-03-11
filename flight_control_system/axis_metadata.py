from dataclasses import dataclass
from .types import Axis, InputChannel, OutputChannel


@dataclass(frozen=True)
class AxisLabels:
    """
    Axis input and output (state) channels labels.

    Attributes
    ----------
    states : tuple[str, ...]
        Output (state) channels LaTeX-like labels.
    state_channels : tuple[OutputChannel, ...]
        Output (state) channels variables.
    state_units : tuple[str, ...]
        Output (state) channels units.
    inputs : tuple[str, ...]
        Input channels LaTeX-like labels.
    input_channels : tuple[InputChannel, ...]
        Input channels variables.
    input_units : tuple[str, ...]
        Input channels units.
    """

    states: tuple[str, ...]
    state_channels: tuple[OutputChannel, ...]
    state_units: tuple[str, ...]
    inputs: tuple[str, ...]
    input_channels: tuple[InputChannel, ...]
    input_units: tuple[str, ...]


AXIS_LABELS = {
    Axis.LONG: AxisLabels(
        states=(r"\Delta u", r"\Delta\alpha", r"\Delta\theta", r"\Delta q",),
        state_channels=(OutputChannel.U, OutputChannel.ALPHA, OutputChannel.THETA, OutputChannel.Q,),
        state_units=("m/s", "$^\circ$", "$^\circ$", "$^\circ$/s",),
        inputs=(r"\Delta\delta_e",),
        input_channels=(InputChannel.ELEVATOR,),
        input_units=("$^\circ$",),
    ),
    Axis.LATDIR: AxisLabels(
        states=(r"\Delta\beta", r"\Delta p", r"\Delta r", r"\Delta\phi",),
        state_channels=(OutputChannel.BETA, OutputChannel.P, OutputChannel.R, OutputChannel.PHI,),
        state_units=("$^\circ$", "$^\circ$/s", "$^\circ$/s", "$^\circ$",),
        inputs=(r"\Delta\delta_a", r"\Delta\delta_r",),
        input_channels=(InputChannel.AILERON, InputChannel.RUDDER,),
        input_units=("$^\circ$", "$^\circ$",),
    ),
}