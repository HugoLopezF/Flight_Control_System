from dataclasses import dataclass
from .types import Axis, InputChannel, OutputChannel


@dataclass(frozen=True)
class AxisLabels:
    states: tuple[str, ...]
    state_channels: tuple[OutputChannel, ...]
    inputs: tuple[str, ...]
    input_channels: tuple[InputChannel, ...]


AXIS_LABELS = {
    Axis.LONG: AxisLabels(
        states=(r"\Delta u", r"\Delta\alpha", r"\Delta\theta", r"\Delta q",),
        state_channels=(OutputChannel.U, OutputChannel.ALPHA, OutputChannel.THETA, OutputChannel.Q,),
        inputs=(r"\Delta\delta_e",),
        input_channels=(InputChannel.ELEVATOR,),
    ),
    Axis.LATDIR: AxisLabels(
        states=(r"\Delta\beta", r"\Delta p", r"\Delta r", r"\Delta\phi",),
        state_channels=(OutputChannel.BETA, OutputChannel.P, OutputChannel.R, OutputChannel.PHI,),
        inputs=(r"\Delta\delta_a", r"\Delta\delta_r",),
        input_channels=(InputChannel.AILERON, InputChannel.RUDDER,),
    ),
}