import os
import sys
import numpy as np

if __package__ is None or __package__ == "":
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    from examples.common import (
        initialize,
        MODEL,
        long_feedback_gains,
        latdir_feedback_gains,
        long_sas_actuators,
        long_sas_sensors,
        latdir_sas_actuators,
        latdir_sas_sensors,
        long_ap_pid_gains,
        latdir_ap_pid_gains,
    )
else:
    from .common import (
        initialize,
        MODEL,
        long_feedback_gains,
        latdir_feedback_gains,
        long_sas_actuators,
        long_sas_sensors,
        latdir_sas_actuators,
        latdir_sas_sensors,
        long_ap_pid_gains,
        latdir_ap_pid_gains,
    )

from analysis_tools.time_response import TimeResponseAnalyzer
from analysis_tools.input_signals import InputSignal
from flight_control_system.autopilot import AP, AutopilotMode
from flight_control_system.sas import SAS
from flight_control_system.types import Axis


def main() -> None:
    _, lin_sys = initialize(MODEL)
    tra = TimeResponseAnalyzer(lin_sys)

    # Build SAS first
    sas_builder = SAS(lin_sys)
    sas_long, _, _ = sas_builder.build_sas(
        axis=Axis.LONG,
        feedback_gains=long_feedback_gains,
        actuators=long_sas_actuators,
        sensors=long_sas_sensors,
    )
    sas_latdir, _, _ = sas_builder.build_sas(
        axis=Axis.LATDIR,
        feedback_gains=latdir_feedback_gains,
        actuators=latdir_sas_actuators,
        sensors=latdir_sas_sensors,
    )

    # Longitudinal AP
    ap_long = AP(sas_long)
    ap_sys_long, _, _ = ap_long.build_ap(
        axis=Axis.LONG,
        mode=AutopilotMode.THETA_HOLD,
        sensors={},
        pid_gains=long_ap_pid_gains,
    )
    tra.plot_ap_response(
        ap_sys=ap_sys_long,
        t=np.linspace(0.0, 200.0, 10000),
        axis=Axis.LONG,
        input_type=InputSignal.SAT_RAMP,
        amp=1.0,
        t_end=0.2,
        showfig=True,
    )

    # Lateral-directional AP
    ap_latdir = AP(sas_latdir)
    ap_sys_latdir, _, _ = ap_latdir.build_ap(
        axis=Axis.LATDIR,
        mode=AutopilotMode.PHI_HOLD,
        sensors={},
        pid_gains=latdir_ap_pid_gains,
    )
    tra.plot_ap_response(
        ap_sys=ap_sys_latdir,
        t=np.linspace(0.0, 200.0, 10000),
        axis=Axis.LATDIR,
        input_type=InputSignal.SAT_RAMP,
        amp=15.0,
        t_end=15.0,
        showfig=True,
    )


if __name__ == "__main__":
    main()
