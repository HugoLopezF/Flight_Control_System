from aircraft.aircraft import Aircraft
from flight_control_system.sas import SAS

def main():
    myaircraft = Aircraft("Learjet_24")
    myaircraft.stab_der.calculate_all()
    mySAS = SAS(myaircraft)
    mySAS.plot_cbc_response(input_type="sat_ramp", amp=1, t_end=20)
    a = 1


if __name__ == "__main__":
    main()