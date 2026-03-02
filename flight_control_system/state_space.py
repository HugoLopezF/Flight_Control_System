from utilities.constants import g
import numpy as np
import sympy as sp
import control
from math import sin, cos, tan
from dataclasses import dataclass
from .types import Axis
from flight_dynamics.aircraft import Aircraft


@dataclass(frozen=True)
class LinearizedMatrices:
    E: np.ndarray
    A_prime: np.ndarray
    B_prime: np.ndarray


@dataclass(frozen=True)
class StandardMatrices:
    A: np.ndarray
    B: np.ndarray


@dataclass(frozen=True)
class ChannelTF:
    tf: control.TransferFunction
    zeros: np.ndarray
    poles: np.ndarray
    gain: float
    dc_gain: float
    wn: np.ndarray
    zeta: np.ndarray
    factored: str


class LinearizedSystem:
    """
    Aircraft's linearized dynamic system.

    Attributes
    ----------
    AXES : tuple
        Dynamic system axes.
    aircraft : Aircraft
        Aircraft to analyze.
    _raw : dict[Axis, LinearizedMatrices]
        Raw form matrices of the linearized dynamic system as in E * x_dot(t) = A' * x(t) + B' * u(t).
    _std : dict[Axis, StandardMatrices]
        Standard form matrices of the linearized dynamic system as in x_dot(t) = A * x(t) + B * u(t).
    _tf : dict[Axis, sp.Matrix]
        Transfer functions of the linearized dynamic system.

    Methods
    ----------
    get_TF()
        
    get_sys()
        Calculate dynamic space state system.
    compute_all_axes()
        Calculate all standard form space state matrices.
    standard_matrices()
        Calculate standard form space state matrices for the selected axis.
    linearized_matrices()
        Calculate linearized system matrices for the selected axis.
    get_long_matrices()
        Calculate longitudinal space state matrices.
    get_latdir_matrices()
        Calculate lateral-directional space state matrices.
    """

    AXES = (Axis.LONG, Axis.LATDIR)

    def __init__(self, aircraft: Aircraft):
        """
        Aircraft's linearized dynamic system.

        Attributes
        ----------
        AXES : tuple
            Dynamic system axes.
        aircraft : Aircraft
            Aircraft to analyze.
        _raw : dict[Axis, LinearizedMatrices]
            Raw form matrices of the linearized dynamic system as in E * x_dot(t) = A' * x(t) + B' * u(t).
        _std : dict[Axis, StandardMatrices]
            Standard form matrices of the linearized dynamic system as in x_dot(t) = A * x(t) + B * u(t).
        _tf : dict[Axis, sp.Matrix]
            Transfer functions of the linearized dynamic system.
        """

        self.aircraft = aircraft
        self._raw: dict[Axis, LinearizedMatrices] = {}
        self._std: dict[Axis, StandardMatrices] = {}
        self._tf: dict[Axis, sp.Matrix] = {}


    def get_TF(self): # TODO: Retrieve transfer functions efficiently
        pass

    def get_sys(self, axis: Axis) -> control.StateSpace:
        """
        Calculate dynamic space state system.

        Attributes
        ----------
        axis : Axis
            Axis to obtain matrices for.

        Returns
        ----------
        control.StateSpace
            Dynamic system state space.
        """

        # Calculate system matrices
        std = self.standard_matrices(axis)
        A, B = std.A, std.B
        m = A.shape[0]
        n = B.shape[1] if len(B.shape) > 1 else 1
        C = np.eye(A.shape[0])
        D = np.zeros((m, n))

        return control.ss(A, B, C, D)


    def compute_all_axes(self) -> None:
        """
        Calculate all standard form space state matrices.

        """

        for axis in self.AXES:
            self.standard_matrices(axis)


    def standard_matrices(self, axis: Axis) -> StandardMatrices:
        """
        Calculate standard form space state matrices for the selected axis.

        Attributes
        ----------
        axis : Axis
            Axis to obtain matrices for.

        Returns
        ----------
        StandardMatrices
            Standard form matrices of the linearized dynamic system 
            as in x_dot(t) = A * x(t) + B * u(t).
        """
        
        if axis not in self._std:
            raw = self.linearized_matrices(axis)
            A = np.linalg.solve(raw.E, raw.A_prime)
            B = np.linalg.solve(raw.E, raw.B_prime)
            self._std[axis] = StandardMatrices(A=A, B=B)
        return self._std[axis]


    def linearized_matrices(self, axis: Axis) -> LinearizedMatrices:
        """
        Calculate linearized system matrices for the selected axis.

        Attributes
        ----------
        axis : Axis
            Axis to obtain matrices for.

        Returns
        ----------
        LinearizedMatrices
            Raw form matrices of the linearized dynamic system
            as in E * x_dot(t) = A' * x(t) + B' * u(t).
        """

        if axis not in self._raw:
            if axis is Axis.LONG:
                E, A_prime, B_prime = self.get_long_matrices()
            elif axis is Axis.LATDIR:
                E, A_prime, B_prime = self.get_latdir_matrices()
            else:
                raise ValueError(f"Invalid axis: {axis}")
            self._raw[axis] = LinearizedMatrices(E=E, A_prime=A_prime, B_prime=B_prime)
        return self._raw[axis]


    def get_long_matrices(self) -> tuple[np.array, np.array, np.array]:
        """
        Calculate longitudinal space state matrices.

        Attributes
        ----------
        axis : Axis
            Axis to obtain matrices for.

        Returns
        ----------
        tuple[np.array, np.array, np.array]
            Raw form matrices of the linearized dynamic system
            as in E * x_dot(t) = A' * x(t) + B' * u(t).
        """

        W = self.aircraft.mass_prop.W
        u_s = self.aircraft.flight_cond.u_s
        theta_s = self.aircraft.flight_cond.theta_s
        Zw_dot = self.aircraft.stab_der.long.Zw_dot
        Zq = self.aircraft.stab_der.long.Zq
        Mw_dot = self.aircraft.stab_der.long.Mw_dot
        Mq = self.aircraft.stab_der.long.Mq
        I_yy = self.aircraft.mass_prop.I_yy
        Xu = self.aircraft.stab_der.long.Xu
        Xw = self.aircraft.stab_der.long.Xw
        Zu = self.aircraft.stab_der.long.Zu
        Zw = self.aircraft.stab_der.long.Zw
        Mu = self.aircraft.stab_der.long.Mu
        Mw = self.aircraft.stab_der.long.Mw
        Xdelta_e = self.aircraft.stab_der.long.Xdelta_e
        Zdelta_e = self.aircraft.stab_der.long.Zdelta_e
        Mdelta_e = self.aircraft.stab_der.long.Mdelta_e

        E = np.array(
            [     
                [W, 0.0, 0.0, 0.0],
                [0.0, u_s * (W - Zw_dot), 0.0, 0.0],
                [0.0, -u_s * Mw_dot, 0.0, I_yy],
                [0.0, 0.0, 1, 0.0]
            ],
            dtype=float,
        )
        A_prime = np.array(
            [
                [Xu, u_s * Xw, -W * g * cos(theta_s), 0.0],
                [Zu, u_s * Zw, -W * g * sin(theta_s), W * u_s + Zq],
                [Mu, u_s * Mw, 0.0, Mq],
                [0.0, 0.0, 0.0, 1.0]
            ],
            dtype=float,
        )
        B_prime = np.array([[Xdelta_e], [Zdelta_e], [Mdelta_e], [0.0]], dtype=float)

        return E, A_prime, B_prime
    
    
    def get_latdir_matrices(self) -> tuple[np.array, np.array, np.array]:
        """
        Calculate lateral-directional space state matrices.

        Attributes
        ----------
        axis : Axis
            Axis to obtain matrices for.

        Returns
        ----------
        tuple[np.array, np.array, np.array]
            Raw form matrices of the linearized dynamic system
            as in E * x_dot(t) = A' * x(t) + B' * u(t).
        """

        W = self.aircraft.mass_prop.W
        u_s = self.aircraft.flight_cond.u_s
        theta_s = self.aircraft.flight_cond.theta_s
        I_xx = self.aircraft.mass_prop.I_xx
        I_zz = self.aircraft.mass_prop.I_zz
        I_xz = self.aircraft.mass_prop.I_xz
        Yv = self.aircraft.stab_der.latdir.Yv
        Yp = self.aircraft.stab_der.latdir.Yp
        Yr = self.aircraft.stab_der.latdir.Yr
        Lv = self.aircraft.stab_der.latdir.Lv
        Lp = self.aircraft.stab_der.latdir.Lp
        Lr = self.aircraft.stab_der.latdir.Lr
        Nv = self.aircraft.stab_der.latdir.Nv
        Np = self.aircraft.stab_der.latdir.Np
        Nr = self.aircraft.stab_der.latdir.Nr
        Ydelta_r = self.aircraft.stab_der.latdir.Ydelta_r
        Ldelta_a = self.aircraft.stab_der.latdir.Ldelta_a
        Ldelta_r = self.aircraft.stab_der.latdir.Ldelta_r
        Ndelta_a = self.aircraft.stab_der.latdir.Ndelta_a
        Ndelta_r = self.aircraft.stab_der.latdir.Ndelta_r

        E = np.array(
            [
                [u_s * W, 0.0, 0.0, 0.0],
                [0.0, I_xx, -I_xz, 0.0],
                [0.0, -I_xz, I_zz, 0.0],
                [0.0, 0.0, 0.0, 1.0]
            ],
            dtype=float,
        )

        A_prime = np.array(
            [
                [u_s * Yv, Yp, Yr - W * u_s, W * g * cos(theta_s)],
                [u_s * Lv, Lp, Lr, 0.0],
                [u_s * Nv, Np, Nr, 0.0],
                [0.0, 1.0, tan(theta_s), 0.0]
            ],
            dtype=float,
        )

        B_prime = np.array(
            [
                [0.0, Ydelta_r],
                [Ldelta_a, Ldelta_r],
                [Ndelta_a, Ndelta_r],
                [0.0, 0.0]
            ],
            dtype=float,
        )
        
        return E, A_prime, B_prime