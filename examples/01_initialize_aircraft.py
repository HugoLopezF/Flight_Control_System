import os
import sys

if __package__ is None or __package__ == "":
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    from examples.common import initialize, MODEL
else:
    from .common import initialize, MODEL

from flight_control_system.types import Axis


def main() -> None:
    aircraft, lin_sys = initialize(MODEL)
    sys_long = lin_sys.get_sys(Axis.LONG)
    sys_latdir = lin_sys.get_sys(Axis.LATDIR)

    print(f"Model: {aircraft.model}")
    print(f"Longitudinal system shape: {sys_long.shape}")
    print(f"Lateral-directional system shape: {sys_latdir.shape}")


if __name__ == "__main__":
    main()
