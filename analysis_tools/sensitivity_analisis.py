import numpy as np
import matplotlib.pyplot as plt
import control
from pathlib import Path
from typing import Iterable, Mapping
from dataclasses import dataclass
from collections.abc import Iterable
from flight_dynamics.aircraft import Aircraft
from flight_dynamics.data_classes import LinearizationParameters
from flight_dynamics.stability_derivatives import StabilityDerivativesCalculator
from flight_control_system.state_space import LinearizedSystem
from flight_control_system.types import Axis


@dataclass(frozen=True)
class SweepPoint:
    value: float
    axis: Axis
    poles: np.ndarray
    zeros: np.ndarray
    sys: control.StateSpace


class SensitivityAnalyzer:
    LONG_COEFFS = ("Cm_alpha", "Cm_q")
    LATDIR_COEFFS = ("Cn_beta", "Cn_r")

    SUPPORTED = {
        "Cm_alpha": (Axis.LONG, "long", "Cm_alpha"),
        "Cm_q": (Axis.LONG, "long", "Cm_q"),
        "Cn_r": (Axis.LATDIR, "latdir", "Cn_r"),
        "Cn_beta": (Axis.LATDIR, "latdir", "Cn_beta"),
    }

    def __init__(self, model: str, fig_root: str | Path = "sens_study"):
        self.model = model
        self.fig_root = Path(fig_root)

    def _figure_dir(self) -> Path:
        """
        Construct figure saving directory.
        """
        figdir = self.fig_root / self.model
        figdir.mkdir(parents=True, exist_ok=True)
        return figdir

    @staticmethod
    def _recompute_derivatives(ac: Aircraft) -> None:
        params = LinearizationParameters(
            geom=ac.geom,
            mass_prop=ac.mass_prop,
            flight_cond=ac.flight_cond,
            cond_coeffs=ac.cond_coeffs,
            stab_coeffs=ac.stab_coeffs,
        )
        ac.stab_der.long = StabilityDerivativesCalculator.calculate_longitudinal(params)
        ac.stab_der.latdir = StabilityDerivativesCalculator.calculate_lateral_directional(params)

    def _group_coeffs(self, group: Axis | str) -> tuple[str, str]:
        if group in (Axis.LONG, "long"):
            return self.LONG_COEFFS
        if group in (Axis.LATDIR, "latdir"):
            return self.LATDIR_COEFFS
        raise ValueError("group must be 'long' or 'latdir'")

    def build_systems(self, coeff_name: str, values: Iterable[float]) -> list[tuple[float, Axis, LinearizedSystem]]:
        if coeff_name not in self.SUPPORTED:
            allowed = ", ".join(self.SUPPORTED.keys())
            raise ValueError(f"Unsupported coefficient '{coeff_name}'. Allowed: {allowed}")

        axis, group, attr = self.SUPPORTED[coeff_name]
        systems: list[tuple[float, Axis, LinearizedSystem]] = []

        for v in values:
            v = float(v)
            ac = Aircraft(self.model)
            setattr(getattr(ac.stab_coeffs, group), attr, v)
            self._recompute_derivatives(ac)
            lin_sys = LinearizedSystem(ac)
            systems.append((v, axis, lin_sys))

        return systems

    def sweep(self, coeff_name: str, values: Iterable[float]) -> list[SweepPoint]:
        points: list[SweepPoint] = []
        for v, axis, lin_sys in self.build_systems(coeff_name, values):
            sys = lin_sys.get_sys(axis)
            poles = np.asarray(control.poles(sys), dtype=complex)
            zeros = np.asarray(control.zeros(sys), dtype=complex)
            points.append(SweepPoint(value=v, axis=axis, poles=poles, zeros=zeros, sys=sys))
        return points

    def plot_pzmap(
        self,
        group: Axis | str,
        coeff_values: Mapping[str, Iterable[float]],
        show_zeros: bool = True,
        savefig: bool = False,
        showfig: bool = False,
    ):
        coeffs = self._group_coeffs(group)
        sweeps = {c: self.sweep(c, coeff_values[c]) for c in coeffs}

        fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharex=True, sharey=True)
        axes = np.atleast_1d(axes)

        for ax, coeff in zip(axes, coeffs):
            points = sweeps[coeff]
            vals = np.array([p.value for p in points], dtype=float)
            vmin, vmax = vals.min(), vals.max()
            if np.isclose(vmin, vmax):
                vmax = vmin + 1.0

            norm = plt.Normalize(vmin, vmax)
            cmap = plt.get_cmap("viridis")

            for p in points:
                color = cmap(norm(p.value))
                ax.plot(p.poles.real, p.poles.imag, "x", color=color, ms=7)
                if show_zeros and p.zeros.size:
                    ax.plot(p.zeros.real, p.zeros.imag, "o", mfc="none", mec=color, ms=6)

            sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
            sm.set_array([])
            fig.colorbar(sm, ax=ax, pad=0.02, label=coeff)

            ax.axhline(0.0, color="k", lw=0.8, alpha=0.6)
            ax.axvline(0.0, color="k", lw=0.8, alpha=0.6)
            ax.grid(True, which="both", linestyle=":", alpha=0.5)
            ax.set_title(coeff)
            ax.set_xlabel("Real [1/s]")
            ax.set_ylabel("Imag [rad/s]")

        group_name = "long" if group in (Axis.LONG, "long") else "latdir"
        fig.suptitle(f"PZ sensitivity grid ({group_name})")
        fig.tight_layout()

        if savefig:
            fig.savefig(self._figure_dir() / f"sensitivity_{group_name}_grid.png", dpi=180)
        if showfig:
            plt.show()
        else:
            plt.close(fig)