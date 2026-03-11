import os
import sys

if __package__ is None or __package__ == "":
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    from examples.common import (
        initialize,
        MODEL,
        long_factors,
        latdir_factors,
        long_sas_actuators,
        long_sas_sensors,
        latdir_sas_actuators,
        latdir_sas_sensors,
        compute_long_sas_gain_sweep,
        compute_latdir_sas_gain_sweep,
    )
else:
    from .common import (
        initialize,
        MODEL,
        long_factors,
        latdir_factors,
        long_sas_actuators,
        long_sas_sensors,
        latdir_sas_actuators,
        latdir_sas_sensors,
        compute_long_sas_gain_sweep,
        compute_latdir_sas_gain_sweep,
    )

from analysis_tools.sensitivity_analisis import SensitivityAnalyzer
from analysis_tools.freq_analysis import FrequencyAnalyzer
from flight_control_system.sas_design import sweep_feedback_gains
from flight_control_system.types import Axis, OutputChannel, InputChannel


def main() -> None:
    aircraft, lin_sys = initialize(MODEL)
    sa = SensitivityAnalyzer(MODEL)

    # Longitudinal SAS sweep
    long_sas_gains = compute_long_sas_gain_sweep(aircraft, long_factors)
    long_sas_points = sweep_feedback_gains(
        lin_sys=lin_sys,
        axis=Axis.LONG,
        gain_values=long_sas_gains,
        actuators=long_sas_actuators,
        sensors=long_sas_sensors,
    )
    sa.plot_sas_pzmap(long_sas_points, Axis.LONG, OutputChannel.ALPHA, OutputChannel.Q, showfig=True)

    # Optional Nichols comparison
    sys_map_long = {
        rf"$K_{{\alpha}}$={p.feedback_gains[OutputChannel.ALPHA]:.1f}, $K_{{q}}$={p.feedback_gains[OutputChannel.Q]:.1f}": p.sas_sys
        for p in long_sas_points
    }
    # FrequencyAnalyzer.compare_sas_nichols(
    #     sys_map=sys_map_long,
    #     axis=Axis.LONG,
    #     out_ch=OutputChannel.ALPHA,
    #     in_ch=InputChannel.ELEVATOR,
    #     showfig=True,
    # )

    # Lateral-directional SAS sweep
    latdir_sas_gains = compute_latdir_sas_gain_sweep(aircraft, latdir_factors)
    latdir_sas_points = sweep_feedback_gains(
        lin_sys=lin_sys,
        axis=Axis.LATDIR,
        gain_values=latdir_sas_gains,
        actuators=latdir_sas_actuators,
        sensors=latdir_sas_sensors,
    )
    sa.plot_sas_pzmap(latdir_sas_points, Axis.LATDIR, OutputChannel.BETA, OutputChannel.R, showfig=True)

    # Optional Nichols comparison
    sys_map_latdir = {
        rf"$K_{{\beta}}$={p.feedback_gains[OutputChannel.BETA]:.1f}, $K_{{r}}$={p.feedback_gains[OutputChannel.R]:.1f}": p.sas_sys
        for p in latdir_sas_points
    }
    # FrequencyAnalyzer.compare_sas_nichols(
    #     sys_map=sys_map_latdir,
    #     axis=Axis.LATDIR,
    #     out_ch=OutputChannel.BETA,
    #     in_ch=InputChannel.AILERON,
    #     showfig=True,
    # )


if __name__ == "__main__":
    main()
