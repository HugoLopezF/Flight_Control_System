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
from flight_control_system.sas_design import sweep_feedback_gains

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
    # myresponse_analyzer.plot_full_response(t=t, input_type=InputSignal.STEP, amp=1/4.785714285714261, showfig=True)

    # # Compare actuators
    # act_map = {
    #     r'1st order lag $\omega_{B} = 1rad/s$': Actuator(omega_b=1.0, order='first').tf(),
    #     r'1st order lag $\omega_{B} = 25rad/s$': Actuator(omega_b=25.0, order='first').tf(),
    #     r'2nd order lag $\omega_{B} = 25rad/s$, $\xi = 0.3$': Actuator(omega_b=25.0, damp=0.3, order='second').tf(),
    #     r'2nd order lag $\omega_{B} = 25rad/s$, $\xi = 0.8$': Actuator(omega_b=25.0, damp=0.8, order='second').tf()
    # }
    # FrequencyAnalyzer.compare_components(tf_map=act_map, showfig=True)

    # sens_map = {
    #     r'Pade 1st order $\tau = 1s$': Sensor.pade(delay=1, num_order=1, den_order=1).tf(),
    #     r'Pade 2nd order $\tau = 1s$': Sensor.pade(delay=1, num_order=2, den_order=2).tf()
    # }
    # FrequencyAnalyzer.compare_components(tf_map=sens_map, omega_limits=(1e-2, 1e1), showfig=True)
    # t = np.linspace(0, 3, 10000)
    # TimeResponseAnalyzer.compare_components(tf_map=sens_map, t=t, input_type=InputSignal.SAT_RAMP, amp=1, t_end=1, showfig=True)

    # # Stability coefficients sweep
    # sa = SensitivityAnalyzer("Learjet_24")
    long_factors = {
        "Cm_alpha": np.linspace(0.5, 20.0, 5),
        "Cm_q": np.linspace(0.5, 20.0, 5),
    }
    # sa.plot_pzmap(Axis.LONG, long_factors, showfig=True)    # Cm_alpha + Cm_q in one grid

    latdir_factors = {
        "Cn_beta": np.linspace(0.05, 1.25, 15),
        "Cn_r": np.linspace(0.05, 1.25, 15),
    }
    # sa.plot_pzmap(Axis.LATDIR, latdir_factors, showfig=True)  # Cn_beta + Cn_r in one grid

    # Feedback gain sensitivity analysis
    ## Longitudinal poles
    Cm_alpha = myaircraft.stab_coeffs.long.Cm_alpha
    Cm_q = myaircraft.stab_coeffs.long.Cm_q
    Cm_delta_e = myaircraft.stab_coeffs.long.Cm_delta_e
    c = myaircraft.geom.c
    u_s = myaircraft.flight_cond.u_s

    gain_values = {
        OutputChannel.ALPHA: [-(F_alpha - 1.0) * Cm_alpha / Cm_delta_e for F_alpha in long_factors['Cm_alpha']],
        OutputChannel.Q: [-(F_q - 1.0) * Cm_q / Cm_delta_e * c / (2.0 * u_s) for F_q in long_factors['Cm_q']],
    }

    actuators = {
        InputChannel.ELEVATOR: Actuator.second_order(omega_b=25.0, damp=0.35),
    }
    sensors = {
        OutputChannel.ALPHA: Sensor.first_order(tau=0.1), 
        OutputChannel.Q: Sensor.first_order(tau=0.1),
    }

    long_points = sweep_feedback_gains(
        lin_sys=mylin_sys,
        axis=Axis.LONG,
        gain_values=gain_values,
        actuators=actuators,
        sensors=sensors,
    )

    sa = SensitivityAnalyzer("Learjet_24")
    sa.plot_sas_pzmap(long_points, Axis.LONG, OutputChannel.ALPHA, OutputChannel.Q, showfig=True)

    ## Lateral-directional poles
    Cn_beta = myaircraft.stab_coeffs.latdir.Cn_beta
    Cn_r = myaircraft.stab_coeffs.latdir.Cn_r
    Cn_delta_r = myaircraft.stab_coeffs.latdir.Cn_delta_r
    b = myaircraft.geom.b
    u_s = myaircraft.flight_cond.u_s

    gain_values = {
        OutputChannel.BETA: [-(F_beta - 1.0) * Cn_beta / Cn_delta_r for F_beta in latdir_factors['Cn_beta']],
        OutputChannel.R: [-(F_r - 1.0) * Cn_r / Cn_delta_r * b / (2.0 * u_s) for F_r in latdir_factors['Cn_r']],
    }

    actuators = {
        InputChannel.AILERON: Actuator.second_order(omega_b=25.0, damp=0.35),
        InputChannel.RUDDER: Actuator.second_order(omega_b=25.0, damp=0.35),
    }
    sensors = {
        OutputChannel.BETA: Sensor.first_order(tau=0.1),
        OutputChannel.R: Sensor.first_order(tau=0.1),
    }

    latdir_points = sweep_feedback_gains(
        lin_sys=mylin_sys,
        axis=Axis.LATDIR,
        gain_values=gain_values,
        actuators=actuators,
        sensors=sensors,
    )

    sa.plot_sas_pzmap(latdir_points, Axis.LATDIR, OutputChannel.BETA, OutputChannel.R, showfig=True)

    # ## Nichols
    # sys_map = {
    #     f"Kalpha={p.feedback_gains[OutputChannel.ALPHA]:.2f}, Kq={p.feedback_gains[OutputChannel.Q]:.2f}": p.sas_sys
    #     for p in long_points
    # }
    # FrequencyAnalyzer.compare_sas_nichols(
    #     sys_map=sys_map,
    #     axis=Axis.LONG,
    #     out_ch=OutputChannel.Q,
    #     in_ch=InputChannel.ELEVATOR,
    #     showfig=True,
    # )

    # Initialize SAS
    mySAS = SAS(mylin_sys)

    # LONG
    feedback_gains = {
        OutputChannel.ALPHA: -8.0,
        OutputChannel.Q: -3.0,
    }
    actuators = {
        InputChannel.ELEVATOR: Actuator(),
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

    # LATDIR
    feedback_gains = {
        OutputChannel.BETA: 0.0, # 0.6,
        OutputChannel.R: 0.0, # 1.5,
    }
    actuators = {
        InputChannel.AILERON: Actuator(),
        InputChannel.RUDDER: Actuator(),
    }
    sensors = {
        OutputChannel.BETA: Sensor(),
        OutputChannel.R: Sensor(),
    }
    sas_latdir, dl_latdir, K_latdir = mySAS.build_sas(
        axis=Axis.LATDIR,
        feedback_gains=feedback_gains,
        actuators=actuators,
        sensors=sensors,
    )

    tra = TimeResponseAnalyzer(mylin_sys)
    tra.plot_sas_response(
        sas_sys=sas_long,
        t=np.linspace(0, 1000, 10000),
        axis=Axis.LONG,
        channel=InputChannel.ELEVATOR,
        input_type=InputSignal.SAT_RAMP,
        amp=1.0,
        t_end = 0.2,
        showfig=True,
    )

    tra.plot_sas_response(
        sas_sys=sas_latdir,
        t=np.linspace(0, 1000, 10000),
        axis=Axis.LATDIR,
        channel=InputChannel.AILERON,
        input_type=InputSignal.PULSE,
        amp=1.0,
        t_end = 22.0,
        showfig=True,
    )

if __name__ == "__main__":
    main()
