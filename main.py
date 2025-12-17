from aircraft.aircraft import Aircraft
from flight_control_system.sas import SAS

def main():
    myaircraft = Aircraft("Learjet_24")
    myaircraft.stab_der.calculate_all()
    mySAS = SAS(myaircraft)
    mySAS.plot_bode()
    a = 1


if __name__ == "__main__":
    main()