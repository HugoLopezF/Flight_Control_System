import numpy as np
from flight_dynamics.aircraft import Aircraft
from flight_control_system.sas import SAS
from flight_control_system.state_space import LinearizedSystem
from analysis_tools.freq_analysis import FrequencyAnalyzer
from analysis_tools.time_response import TimeResponseAnalyzer
from analysis_tools.input_signals import InputSignal
from analysis_tools.sensitivity_analisis import SensitivityAnalyzer
from flight_control_system.types import Axis, OutputChannel, InputChannel
from flight_control_system.actuator import Actuator
from flight_control_system.sensor import Sensor
from flight_control_system.filter import Filter, FilterType
from flight_control_system.sas_design import sweep_feedback_gains
from flight_control_system.autopilot import AP, AutopilotMode, PIDGains
from flight_control_system.autopilot_design import sweep_pid_gains

def main():
    AXES = (Axis.LONG, Axis.LATDIR)

    # Initialize aircraft to study
    myaircraft = Aircraft("Learjet_24")

    # Get matrices
    mylin_sys = LinearizedSystem(myaircraft)

    # # Frequency analysis
    # myfreq_analyzer = FrequencyAnalyzer(mylin_sys)
    # myfreq_analyzer.plot_bode(axis=Axis.LONG, showfig=True)
    # myfreq_analyzer.plot_bode(axis=Axis.LATDIR, showfig=True)
    # myfreq_analyzer.plot_nichols(axis=Axis.LONG)
    # freq_metrics = myfreq_analyzer.channel_metrics(axis=Axis.LONG, i=0, j=0)

    # # Response analysis
    # myresponse_analyzer = TimeResponseAnalyzer(mylin_sys)
    # t = np.linspace(0, 1000, 10000)
    # myresponse_analyzer.plot_full_response(t=t, input_type=InputSignal.STEP, amp=1, showfig=True)

    # # Compare actuators
    # act_map = {
    #     r'1st order lag $\omega_{B} = 1rad/s$': Actuator.first_order(omega_b=1.0).tf(),
    #     r'1st order lag $\omega_{B} = 25rad/s$': Actuator.first_order(omega_b=25.0).tf(),
    #     r'2nd order lag $\omega_{B} = 25rad/s$, $\xi = 0.3$': Actuator.second_order(omega_b=25.0, damp=0.3).tf(),
    #     r'2nd order lag $\omega_{B} = 25rad/s$, $\xi = 0.8$': Actuator.second_order(omega_b=25.0, damp=0.8).tf()
    # }
    # FrequencyAnalyzer.compare_components(tf_map=act_map, showfig=True)

    # # Compare sensors
    # sens_map = {
    #     r'Pade 1st order $\tau = 1s$': Sensor.pade(delay=1, num_order=1, den_order=1).tf(),
    #     r'Pade 2nd order $\tau = 1s$': Sensor.pade(delay=1, num_order=2, den_order=2).tf()
    # }
    # FrequencyAnalyzer.compare_components(tf_map=sens_map, omega_limits=(1e-2, 1e1), showfig=True)
    # t = np.linspace(0, 3, 10000)
    # TimeResponseAnalyzer.compare_components(tf_map=sens_map, t=t, input_type=InputSignal.SAT_RAMP, amp=1, t_end=1, showfig=True)

    # Stability coefficients sweep
    sa = SensitivityAnalyzer("Learjet_24")
    long_factors = {
        "Cm_alpha": [-0.5, -0.2, 0, 0.5, 1, 2, 3, 4], 
        "Cm_q": [-0.5, -0.2, 0, 0.5, 1, 2, 3, 4],
    }
    sa.plot_pzmap(Axis.LONG, long_factors, showfig=True)    # Cm_alpha + Cm_q in one grid

    latdir_factors = {
        "Cn_beta": [0, 0.5, 1, 2, 3],
        "Cn_r": [0, 0.5, 1, 2, 3, 5], 
    }
    sa.plot_pzmap(Axis.LATDIR, latdir_factors, showfig=True)  # Cn_beta + Cn_r in one grid

    # Feedback gain sensitivity analysis
    ## Longitudinal poles
    Cm_alpha = myaircraft.stab_coeffs.long.Cm_alpha
    Cm_q = myaircraft.stab_coeffs.long.Cm_q
    Cm_delta_e = myaircraft.stab_coeffs.long.Cm_delta_e
    c = myaircraft.geom.c
    u_s = myaircraft.flight_cond.u_s

    long_gain_values = {
        OutputChannel.ALPHA: [-(F_alpha - 1.0) * Cm_alpha / Cm_delta_e for F_alpha in long_factors['Cm_alpha']],
        OutputChannel.Q: [-(F_q - 1.0) * Cm_q / Cm_delta_e * c / (2.0 * u_s) for F_q in long_factors['Cm_q']],
    }

    actuators = {
        InputChannel.ELEVATOR: Actuator.first_order(omega_b=10),
    }
    sensors = {
        OutputChannel.ALPHA: Sensor(), 
        OutputChannel.Q: Sensor(),
    }

    long_points = sweep_feedback_gains(
        lin_sys=mylin_sys,
        axis=Axis.LONG,
        gain_values=long_gain_values,
        actuators=actuators,
        sensors=sensors,
    )

    sa = SensitivityAnalyzer("Learjet_24")
    sa.plot_sas_pzmap(long_points, Axis.LONG, OutputChannel.ALPHA, OutputChannel.Q, showfig=True)

    ## Nichols
    sys_map = {
        rf"$K_{{\alpha}}$={p.feedback_gains[OutputChannel.ALPHA]:.1f}, $K_{{q}}$={p.feedback_gains[OutputChannel.Q]:.1f}": p.sas_sys
        for p in long_points
    }

    # FrequencyAnalyzer.compare_sas_nichols(
    #     sys_map=sys_map,
    #     axis=Axis.LONG,
    #     out_ch=OutputChannel.ALPHA,
    #     in_ch=InputChannel.ELEVATOR,
    #     showfig=True,
    # )

    # FrequencyAnalyzer.compare_sas_nichols(
    #     sys_map=sys_map,
    #     axis=Axis.LONG,
    #     out_ch=OutputChannel.Q,
    #     in_ch=InputChannel.ELEVATOR,
    #     showfig=True,
    # )

    ## Lateral-directional poles
    Cn_beta = myaircraft.stab_coeffs.latdir.Cn_beta
    Cn_r = myaircraft.stab_coeffs.latdir.Cn_r
    Cn_delta_r = myaircraft.stab_coeffs.latdir.Cn_delta_r
    b = myaircraft.geom.b
    u_s = myaircraft.flight_cond.u_s

    # latdir_gain_values = {
    #     OutputChannel.BETA: [-(F_beta - 1.0) * Cn_beta / Cn_delta_r for F_beta in latdir_factors['Cn_beta']],
    #     OutputChannel.R: [-(F_r - 1.0) * Cn_r / Cn_delta_r * b / (2.0 * u_s) for F_r in latdir_factors['Cn_r']],
    # }
    latdir_gain_values = {
        OutputChannel.BETA: [-1, -2, -3, -4],
        OutputChannel.R: [0.25, 0.5, 0.75, 1],
    }

    actuators = {
        InputChannel.AILERON: Actuator.first_order(omega_b=10),
        InputChannel.RUDDER: Actuator.first_order(omega_b=10),
    }
    sensors = {
        OutputChannel.BETA: Sensor(),
        OutputChannel.R: Sensor(),
    }

    latdir_points = sweep_feedback_gains(
        lin_sys=mylin_sys,
        axis=Axis.LATDIR,
        gain_values=latdir_gain_values,
        actuators=actuators,
        sensors=sensors,
    )

    sa.plot_sas_pzmap(latdir_points, Axis.LATDIR, OutputChannel.BETA, OutputChannel.R, showfig=True)

    ## Nichols
    sys_map = {
        rf"$K_{{\beta}}$={p.feedback_gains[OutputChannel.BETA]:.1f}, $K_{{r}}$={p.feedback_gains[OutputChannel.R]:.1f}": p.sas_sys
        for p in latdir_points
    }

    # FrequencyAnalyzer.compare_sas_nichols(
    #     sys_map=sys_map,
    #     axis=Axis.LATDIR,
    #     out_ch=OutputChannel.BETA,
    #     in_ch=InputChannel.AILERON,
    #     showfig=True,
    # )

    # FrequencyAnalyzer.compare_sas_nichols(
    #     sys_map=sys_map,
    #     axis=Axis.LATDIR,
    #     out_ch=OutputChannel.R,
    #     in_ch=InputChannel.AILERON,
    #     showfig=True,
    # )

    # Select gains after sensitivity study and build your SAS
    mySAS = SAS(mylin_sys)

    ## Longitudinal SAS
    feedback_gains = {
        OutputChannel.ALPHA: -0.75,
        OutputChannel.Q: -0.2,
    }
    actuators = {
        InputChannel.ELEVATOR: Actuator.first_order(omega_b=10),
    }
    sensors = {
        OutputChannel.ALPHA: Sensor(),
        OutputChannel.Q: Sensor(),
    }
    sas_long, dl_long, K_long = mySAS.build_sas(
        axis=Axis.LONG,
        feedback_gains=feedback_gains,
        actuators=actuators,
        sensors=sensors,
    )

    tra = TimeResponseAnalyzer(mylin_sys)
    tra.plot_sas_response(
        sas_sys=sas_long,
        t=np.linspace(0, 100, 10000),
        axis=Axis.LONG,
        channel=InputChannel.ELEVATOR,
        input_type=InputSignal.SAT_RAMP,
        amp=1.0,
        t_end = 0.2,
        showfig=True,
    )

    ## Lateral-directional SAS
    feedback_gains = {
        OutputChannel.BETA: -2,
        OutputChannel.R: 0.55,
    }
    actuators = {
        InputChannel.AILERON: Actuator.first_order(omega_b=10),
        InputChannel.RUDDER: Actuator.first_order(omega_b=10),
    }
    sensors = {
        OutputChannel.BETA: Sensor(),
        OutputChannel.R: Sensor(),
    }
    filters = {
        OutputChannel.R: Filter.highpass(omega_b=2), # Washout filter
    }
    sas_latdir, dl_latdir, K_latdir = mySAS.build_sas(
        axis=Axis.LATDIR,
        feedback_gains=feedback_gains,
        actuators=actuators,
        sensors=sensors,
        filters=filters,
    )

    tra.plot_sas_response(
        sas_sys=sas_latdir,
        t=np.linspace(0, 100, 10000),
        axis=Axis.LATDIR,
        channel=InputChannel.AILERON,
        input_type=InputSignal.PULSE,
        amp=1.0,
        t_i=5.0,
        t_end=7.0,
        showfig=True,
    )

    # Autopilot
    ## Longitudinal AP
    long_pid_values = {
        "kp": np.linspace(0.1, 3.0, 3),
        "ki": np.linspace(0.0, 0.8, 3),
    }

    long_ap_points = sweep_pid_gains(
        sas_sys=sas_long,
        axis=Axis.LONG,
        mode=AutopilotMode.THETA_HOLD,
        gain_values=long_pid_values,
        sensors={OutputChannel.THETA: Sensor()},
    )

    sa.plot_ap_pzmap(long_ap_points, Axis.LONG, showfig=True)

    ## Nichols
    sys_map = {
        f"$K_{{p}}$={p.pid_gains.kp:.1f}, $K_{{i}}$={p.pid_gains.ki:.1f}": p.ap_sys
        for p in long_ap_points
    }

    FrequencyAnalyzer.compare_ap_nichols(
        sys_map=sys_map,
        axis=Axis.LONG,
        out_ch=OutputChannel.ALPHA,
        in_ch=InputChannel.ELEVATOR,
        showfig=True,
    )

    FrequencyAnalyzer.compare_ap_nichols(
        sys_map=sys_map,
        axis=Axis.LONG,
        out_ch=OutputChannel.Q,
        in_ch=InputChannel.ELEVATOR,
        showfig=True,
    )

    ## Lateral-directional AP
    latdir_pid_values = {
        "kp": np.linspace(0.1, 3.0, 8),
    }

    latdir_ap_points = sweep_pid_gains(
        sas_sys=sas_latdir,
        axis=Axis.LATDIR,
        mode=AutopilotMode.PHI_HOLD,
        gain_values=latdir_pid_values,
        sensors={OutputChannel.THETA: Sensor()},
    )

    sa.plot_ap_pzmap(latdir_ap_points, Axis.LATDIR, showfig=True)

    ## Nichols
    sys_map = {
        f"$K_{{p}}$={p.pid_gains.kp:.1f}": p.ap_sys
        for p in latdir_ap_points
    }

    FrequencyAnalyzer.compare_ap_nichols(
        sys_map=sys_map,
        axis=Axis.LATDIR,
        out_ch=OutputChannel.BETA,
        in_ch=InputChannel.AILERON,
        showfig=True,
    )

    FrequencyAnalyzer.compare_ap_nichols(
        sys_map=sys_map,
        axis=Axis.LATDIR,
        out_ch=OutputChannel.R,
        in_ch=InputChannel.AILERON,
        showfig=True,
    )

    # Select PID gains after sensitivity study and build your AP
    my_long_AP = AP(sas_long)

    ## Longitudinal AP
    PID_gains = PIDGains(
        kp=-1.0, 
        ki=0.8,
    )
    sensors = {
        OutputChannel.THETA: Sensor(),
    }
    ap_long, cl_z_long, C_long = my_long_AP.build_ap(
        axis=Axis.LONG,
        mode=AutopilotMode.THETA_HOLD,
        sensors=sensors,
        pid_gains=PID_gains,
    )

    tra.plot_ap_response(
        ap_sys=ap_long,
        t=np.linspace(0, 200, 10000),
        axis=Axis.LONG,
        input_type=InputSignal.SAT_RAMP,
        amp=1.0,
        t_end = 0.2,
        showfig=True,
    )

    ## Lateral-directional AP
    my_latdir_AP = AP(sas_latdir)
    PID_gains = PIDGains(
        kp=1.0, 
    )
    sensors = {
        OutputChannel.PHI: Sensor(),
    }
    ap_latdir, cl_z_long, C_long = my_latdir_AP.build_ap(
        axis=Axis.LATDIR,
        mode=AutopilotMode.PHI_HOLD,
        sensors=sensors,
        pid_gains=PID_gains,
    )

    tra.plot_ap_response(
        ap_sys=ap_latdir,
        t=np.linspace(0, 200, 10000),
        axis=Axis.LATDIR,
        input_type=InputSignal.SAT_RAMP,
        amp=15.0,
        t_end = 15.0,
        showfig=True,
    )

if __name__ == "__main__":
    main()
