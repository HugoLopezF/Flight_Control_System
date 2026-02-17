from flight_dynamics.aircraft import Aircraft
from flight_control_system.sas import SAS
import numpy as np

def main():
    # Initialize aircraft to study
    myaircraft = Aircraft("Learjet_24")

    # Initialize SAS
    mySAS = SAS(myaircraft)

    # Get matrices
    sys = mySAS.get_sys('long')

    # Plot step response applied
    t = np.linspace(0, 500, 10000)
    mySAS.plot_cbc_response(t, input_type="sat_ramp", amp=1)

    # Plot ramp response applied during 10s
    t = np.linspace(0, 500, 10000)
    mySAS.plot_cbc_response(t, input_type="sat_ramp", amp=1, t_end=10)


if __name__ == "__main__":
    main()