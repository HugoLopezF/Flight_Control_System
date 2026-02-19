from flight_dynamics.aircraft import Aircraft
from flight_control_system.sas import SAS
from flight_control_system.state_space import LinearizedSystem
from analysis_tools.freq_analysis import FrequencyAnalyzer
from analysis_tools.time_response import TimeResponseAnalyzer
from analysis_tools.input_signals import InputSignal
from flight_control_system.types import Axis 
import numpy as np

def main():
    AXES = (Axis.LONG, Axis.LATDIR)

    # Initialize aircraft to study
    myaircraft = Aircraft("Learjet24")

    # Get matrices
    mylin_sys = LinearizedSystem(myaircraft)
    mylin_sys.compute_all_axes()

    # Frequency analysis
    myfreq_analyzer = FrequencyAnalyzer(mylin_sys)
    # myfreq_analyzer.plot_bode(axis=Axis.LONG)
    # myfreq_analyzer.plot_nichols(axis=Axis.LONG)
    freq_metrics = myfreq_analyzer.channel_metrics(axis=Axis.LONG, i=0, j=0)

    # Response analysis
    myresponse_analyzer = TimeResponseAnalyzer(mylin_sys)
    t = np.linspace(0, 500, 10000)
    myresponse_analyzer.plot_cbc_response(t=t, input_type=InputSignal.STEP, showfig=True)

    # Initialize SAS
    mySAS = SAS(myaircraft)


if __name__ == "__main__":
    main()