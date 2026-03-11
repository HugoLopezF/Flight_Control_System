import os
import sys
import numpy as np

if __package__ is None or __package__ == "":
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from analysis_tools.freq_analysis import FrequencyAnalyzer
from analysis_tools.time_response import TimeResponseAnalyzer
from analysis_tools.input_signals import InputSignal
from flight_control_system.actuator import Actuator
from flight_control_system.sensor import Sensor


def main() -> None:
    # Compare actuators (frequency domain)
    act_map = {
        r"1st order lag $\omega_{B} = 1rad/s$": Actuator.first_order(omega_b=1.0).tf(),
        r"1st order lag $\omega_{B} = 25rad/s$": Actuator.first_order(omega_b=25.0).tf(),
        r"2nd order lag $\omega_{B} = 25rad/s$, $\xi = 0.3$": Actuator.second_order(omega_b=25.0, damp=0.3).tf(),
        r"2nd order lag $\omega_{B} = 25rad/s$, $\xi = 0.8$": Actuator.second_order(omega_b=25.0, damp=0.8).tf(),
    }
    FrequencyAnalyzer.compare_components(tf_map=act_map, showfig=True)

    # Compare sensors (frequency)
    sens_map = {
        r"Pade 1st order $\tau = 1s$": Sensor.pade(delay=1, num_order=1, den_order=1).tf(),
        r"Pade 2nd order $\tau = 1s$": Sensor.pade(delay=1, num_order=2, den_order=2).tf(),
    }
    t = np.linspace(0.0, 3.0, 10000)
    TimeResponseAnalyzer.compare_components(
        tf_map=sens_map,
        t=t,
        input_type=InputSignal.SAT_RAMP,
        amp=1.0,
        t_end=1.0,
        showfig=True,
    )


if __name__ == "__main__":
    main()
