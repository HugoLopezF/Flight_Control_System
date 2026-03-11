import os
import sys

if __package__ is None or __package__ == "":
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    from examples.common import (
        initialize,
        MODEL,
        long_pid_values,
        latdir_pid_values,
        long_feedback_gains,
        latdir_feedback_gains,
        long_sas_actuators,
        long_sas_sensors,
        latdir_sas_actuators,
        latdir_sas_sensors,
    )
else:
    from .common import (
        initialize,
        MODEL,
        long_pid_values,
        latdir_pid_values,
        long_feedback_gains,
        latdir_feedback_gains,
        long_sas_actuators,
        long_sas_sensors,
        latdir_sas_actuators,
        latdir_sas_sensors,
    )

from analysis_tools.sensitivity_analisis import SensitivityAnalyzer
from analysis_tools.freq_analysis import FrequencyAnalyzer
from flight_control_system.autopilot import AutopilotMode
from flight_control_system.autopilot_design import sweep_pid_gains
from flight_control_system.types import Axis, OutputChannel, InputChannel
from flight_control_system.sensor import Sensor
from flight_control_system.sas import SAS


def main() -> None:
    _, lin_sys = initialize(MODEL)
    sa = SensitivityAnalyzer(MODEL)

    # Build SAS first (autopilot is built on top of SAS)
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

    # Longitudinal AP sweep
    long_ap_points = sweep_pid_gains(
        sas_sys=sas_long,
        axis=Axis.LONG,
        mode=AutopilotMode.THETA_HOLD,
        gain_values=long_pid_values,
        sensors={OutputChannel.THETA: Sensor()},
    )
    # Requires plot_ap_pzmap implemented for PID gains (kp/ki/kd)
    sa.plot_ap_pzmap(long_ap_points, Axis.LONG, showfig=True)

    sys_map_long = {
        f"$K_p$={p.pid_gains.kp:.1f}, $K_i$={p.pid_gains.ki:.1f}": p.ap_sys
        for p in long_ap_points
    }
    FrequencyAnalyzer.compare_ap_nichols(
        sys_map=sys_map_long,
        axis=Axis.LONG,
        out_ch=OutputChannel.ALPHA,
        in_ch=InputChannel.ELEVATOR,
        showfig=True,
    )

    # Lateral-directional AP sweep
    latdir_ap_points = sweep_pid_gains(
        sas_sys=sas_latdir,
        axis=Axis.LATDIR,
        mode=AutopilotMode.PHI_HOLD,
        gain_values=latdir_pid_values,
        sensors={OutputChannel.PHI: Sensor()},
    )
    # Requires plot_ap_pzmap implemented for PID gains (kp/ki/kd)
    sa.plot_ap_pzmap(latdir_ap_points, Axis.LATDIR, showfig=True)

    sys_map_latdir = {
        f"$K_p$={p.pid_gains.kp:.1f}": p.ap_sys
        for p in latdir_ap_points
    }
    FrequencyAnalyzer.compare_ap_nichols(
        sys_map=sys_map_latdir,
        axis=Axis.LATDIR,
        out_ch=OutputChannel.BETA,
        in_ch=InputChannel.AILERON,
        showfig=True,
    )


if __name__ == "__main__":
    main()
