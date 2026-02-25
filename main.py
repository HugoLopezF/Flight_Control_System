from flight_dynamics.aircraft import Aircraft
from flight_control_system.sas import SAS
from flight_control_system.state_space import LinearizedSystem
from analysis_tools.freq_analysis import FrequencyAnalyzer
from analysis_tools.time_response import TimeResponseAnalyzer
from analysis_tools.input_signals import InputSignal
from analysis_tools.sensitivity_analisis import SensitivityAnalyzer
from flight_control_system.types import Axis 
from flight_control_system.actuator import Actuator
from flight_control_system.sensor import Sensor
import numpy as np

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
    #     r'Pade 1st order $\tau = 1s$': Sensor(delay=1, num_order=1, den_order=1).tf(),
    #     r'Pade 2nd order $\tau = 1s$': Sensor(delay=1, num_order=2, den_order=2).tf()
    # }
    # FrequencyAnalyzer.compare_components(tf_map=sens_map, omega_limits=(1e-2, 1e1), showfig=True)
    # t = np.linspace(0, 3, 10000)
    # TimeResponseAnalyzer.compare_components(tf_map=sens_map, t=t, input_type=InputSignal.SAT_RAMP, amp=1, t_end=1, showfig=True)

    # # Stability coefficients sweep
    # sa = SensitivityAnalyzer("Learjet_24")
    # vals = {
    #     "Cm_alpha": np.linspace(0.5, 20.0, 11),
    #     "Cm_q": np.linspace(0.5, 20.0, 6),
    # }
    # sa.plot_pzmap(Axis.LONG, vals, showfig=True)    # Cm_alpha + Cm_q in one grid

    # vals = {
    #     "Cn_beta": np.linspace(0.5, 12.0, 6),
    #     "Cn_r": np.linspace(0.5, 20.0, 11),
    # }
    # sa.plot_pzmap(Axis.LATDIR, vals, showfig=True)  # Cn_beta + Cn_r in one grid

    # Initialize SAS
    mySAS = SAS(mylin_sys)
    feedback_gains = {
        'Kalpha': 1.0,
        'Kq': 1.0,
    }
    mySAS.build_sas(axis=Axis.LONG, feedback_gains=feedback_gains)

if __name__ == "__main__":
    main()