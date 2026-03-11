import os
import sys

if __package__ is None or __package__ == "":
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    from examples.common import MODEL, long_factors, latdir_factors
else:
    from .common import MODEL, long_factors, latdir_factors

from analysis_tools.sensitivity_analisis import SensitivityAnalyzer
from flight_control_system.types import Axis


def main() -> None:
    sa = SensitivityAnalyzer(MODEL)
    sa.plot_pzmap(Axis.LONG, long_factors, showfig=True)
    sa.plot_pzmap(Axis.LATDIR, latdir_factors, showfig=True)


if __name__ == "__main__":
    main()
