from aircraft.aircraft import Aircraft

def main():
    myaircraft = Aircraft("S211")
    myaircraft.stab_der.calculate_all()


if __name__ == "__main__":
    main()