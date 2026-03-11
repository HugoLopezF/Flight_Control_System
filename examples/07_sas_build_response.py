import os
import sys
import numpy as np

if __package__ is None or __package__ == "":
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    from examples.common import (
        initialize,
        MODEL,
        long_sas_actuators,
        long_sas_sensors,
        latdir_sas_actuators,
        latdir_sas_sensors,
        latdir_sas_filters,
        long_feedback_gains,
        latdir_feedback_gains,
    )
else:
    from .common import (
        initialize,
        MODEL,
        long_sas_actuators,
        long_sas_sensors,
        latdir_sas_actuators,
        latdir_sas_sensors,
        latdir_sas_filters,
        long_feedback_gains,
        latdir_feedback_gains,
    )

from analysis_tools.time_response import TimeResponseAnalyzer
from analysis_tools.input_signals import InputSignal
from flight_control_system.sas import SAS
from flight_control_system.types import Axis, InputChannel


def main() -> None:
    _, lin_sys = initialize(MODEL)
    tra = TimeResponseAnalyzer(lin_sys)
    sas_builder = SAS(lin_sys)

    # Longitudinal SAS
    sas_long, _, _ = sas_builder.build_sas(
        axis=Axis.LONG,
        feedback_gains=long_feedback_gains,
        actuators=long_sas_actuators,
        sensors=long_sas_sensors,
    )
    tra.plot_sas_response(
        sas_sys=sas_long,
        t=np.linspace(0.0, 100.0, 10000),
        axis=Axis.LONG,
        channel=InputChannel.ELEVATOR,
        input_type=InputSignal.SAT_RAMP,
        amp=1.0,
        t_end=0.2,
        showfig=True,
    )

    # Lateral-directional SAS
    sas_latdir, _, _ = sas_builder.build_sas(
        axis=Axis.LATDIR,
        feedback_gains=latdir_feedback_gains,
        actuators=latdir_sas_actuators,
        sensors=latdir_sas_sensors,
        filters=latdir_sas_filters,
    )
    tra.plot_sas_response(
        sas_sys=sas_latdir,
        t=np.linspace(0.0, 100.0, 10000),
        axis=Axis.LATDIR,
        channel=InputChannel.AILERON,
        input_type=InputSignal.PULSE,
        amp=1.0,
        t_i=5.0,
        t_end=7.0,
        showfig=True,
    )


if __name__ == "__main__":
    main()
