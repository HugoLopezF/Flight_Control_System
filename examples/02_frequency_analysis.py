import os
import sys

if __package__ is None or __package__ == "":
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    from examples.common import initialize, MODEL
else:
    from .common import initialize, MODEL
from analysis_tools.freq_analysis import FrequencyAnalyzer
from flight_control_system.types import Axis


def main() -> None:
    _, lin_sys = initialize(MODEL)
    fa = FrequencyAnalyzer(lin_sys)

    fa.plot_bode(axis=Axis.LONG, showfig=True)
    fa.plot_nichols(axis=Axis.LONG, showfig=True)

    metrics = fa.channel_metrics(axis=Axis.LONG, i=0, j=0)
    print("Channel metrics (LONG, i=0, j=0):")
    for k, v in metrics.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
