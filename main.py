from aircraft.aircraft import Aircraft
from flight_control_system.sas import SAS
import numpy as np

def main():
    # Declare aircraft to study
    myaircraft = Aircraft("Learjet_24")

    # Calculate stability derivatives
    myaircraft.stab_der.calculate_all()

    # Declare SAS
    mySAS = SAS(myaircraft)
    t = np.linspace(0, 150, 10000)
    mySAS.plot_cbc_response(t, input_type="sat_ramp", amp=1, t_end=20)
    a = 1


if __name__ == "__main__":
    main()