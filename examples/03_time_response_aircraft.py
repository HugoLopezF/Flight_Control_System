import os
import sys
import numpy as np

if __package__ is None or __package__ == "":
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    from examples.common import initialize, MODEL
else:
    from .common import initialize, MODEL
from analysis_tools.time_response import TimeResponseAnalyzer
from analysis_tools.input_signals import InputSignal
from flight_control_system.types import Axis, InputChannel


def main() -> None:
    _, lin_sys = initialize(MODEL)
    tra = TimeResponseAnalyzer(lin_sys)

    t = np.linspace(0.0, 100.0, 5000)
    tra.plot_aircraft_response(
        t=t,
        axis=Axis.LONG,
        channel=InputChannel.ELEVATOR,
        input_type=InputSignal.STEP,
        amp=1.0,
        showfig=True,
    )

    # Full response for all inputs/axes (can be slow)
    # tra.plot_full_response(t=t, input_type=InputSignal.STEP, amp=1.0, showfig=True)


if __name__ == "__main__":
    main()
