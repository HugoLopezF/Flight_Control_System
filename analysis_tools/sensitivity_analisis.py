import numpy as np
import matplotlib.pyplot as plt
import control
from pathlib import Path
from dataclasses import dataclass
from collections.abc import Iterable, Mapping
from flight_dynamics.aircraft import Aircraft
from flight_dynamics.data_classes import LinearizationParameters
from flight_dynamics.stability_derivatives import StabilityDerivativesCalculator
from flight_control_system.state_space import LinearizedSystem
from flight_control_system.types import Axis, OutputChannel
from itertools import product
from matplotlib.colors import hsv_to_rgb


@dataclass(frozen=True)
class SweepPoint:
    """
    Aircraft stability coefficients sweep point.

    Attributes
    ----------
    factors : dict[str, float]
        Factor to multiply stability coefficients by.
    axis : Axis
        Axis to analyze.
    poles : np.ndarray
        Aircraft poles.
    """
        
    factors: dict[str, float]
    axis: Axis
    poles: np.ndarray


class SensitivityAnalyzer:
    """
    Sensitivity analyzer class.

    Creates pole-zero map for aircraft, SAS and AP when sweeping parameters.

    Attributes
    ----------
    SUPPORTED : dict[Axis, tuple(str, str)]
        Supported stability coefficients to modify.
    model : tuple
        Aircraft model to analyze.
    fig_root : Aircraft, optional
        Figure root saving path (default is "flight_dynamics").

    Methods
    ----------
    _figure_dir()
        Construct figure saving directory.
    _recompute_derivatives()
        Recompute stability derivatives for updated stability coefficients.
    _add_omega_grid()
        Adds iso-omega grid to pole-zero map plot.
    _add_damping_grid()
        Adds iso-damping grid to pole-zero map plot.
    build_systems()
        Construct aircraft for all factor combinations in the sweep.
    sweep()
        Obtain sweep points.
    plot_pzmap()
        Plot aircraft pole-zero map for each factor combination.
    plot_sas_pzmap()
        Plot Stability Augmentation System pole-zero map for each gain combination.
    """
        
    SUPPORTED = {
        Axis.LONG: ("Cm_alpha", "Cm_q"),
        Axis.LATDIR: ("Cn_beta", "Cn_r"),
    }

    def __init__(self, model: str, fig_root: str | Path = "sens_study"):
        """
        Frequency analyzer class.

        Creates Bode and Nichols plots for aircraft, components, SAS and AP.

        Attributes
        ----------
        SUPPORTED : dict[Axis, tuple(str, str)]
            Supported stability coefficients to modify.
        model : tuple
            Aircraft model to analyze.
        fig_root : Aircraft, optional
            Figure root saving path (default is "sens_study").
        """
            
        self.model = model
        self.fig_root = Path(fig_root)

    def _figure_dir(self) -> Path:
        """
        Construct figure saving directory.

        Returns
        ----------
        Path
            Figure saving directory.
        """

        figdir = self.fig_root / self.model
        figdir.mkdir(parents=True, exist_ok=True)
        return figdir

    @staticmethod
    def _recompute_derivatives(ac: Aircraft) -> None:
        """
        Recompute stability derivatives for updated stability coefficients.

        Attributes
        ----------
        ac : Aircraft
            Aircraft class to analyze.
        """
                
        params = LinearizationParameters(
            geom=ac.geom,
            mass_prop=ac.mass_prop,
            flight_cond=ac.flight_cond,
            cond_coeffs=ac.cond_coeffs,
            stab_coeffs=ac.stab_coeffs,
        )
        ac.stab_der.long = StabilityDerivativesCalculator.calculate_longitudinal(params)
        ac.stab_der.latdir = StabilityDerivativesCalculator.calculate_lateral_directional(params)

    @staticmethod
    def _add_omega_grid(
        ax: plt.Axes,
        omegas: tuple[float, ...] = (0.5, 1.0, 2.0, 5.0, 10.0),
        color: str = "0.35",
        alpha: float = 0.45,
        lw: float = 0.8,
    ) -> None:
        """
        Adds iso-omega grid to pole-zero map plot.

        Attributes
        ----------
        ax : plt.Axes
            Figure axes.
        omegas : tuple[float, ...], optional
            Omegas to plot iso lines for in rad/s (default is (0.5, 1.0, 2.0, 5.0, 10.0)).
        color : str, optional
            Line color (default is "0.35").
        alpha : float, optional
            Line transparency (default is 0.45).
        lw : float, optional
            Linewidth (default is 0.8).
        """
                
        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()

        # Left-half-plane semicircles: sigma^2 + omega_d^2 = omega_n^2
        th = np.linspace(-np.pi / 2, np.pi / 2, 600)

        for wn in omegas:
            if wn <= 0:
                continue

            x = -wn * np.cos(th)
            y =  wn * np.sin(th)

            mask = (x >= x_min) & (x <= x_max) & (y >= y_min) & (y <= y_max)
            if not np.any(mask):
                continue

            x_vis = x[mask]
            y_vis = y[mask]
            ax.plot(x_vis, y_vis, ":", color=color, alpha=alpha, lw=lw, zorder=0)

            # Label from visible segment (not fixed angle)
            idx = int(0.75 * (x_vis.size - 1))  # upper branch
            x_lbl = x_vis[idx]
            y_lbl = y_vis[idx]
            ax.text(
                x_lbl, y_lbl, fr"$\omega_n={wn:g}$",
                fontsize=7, color=color, ha="left", va="bottom",
                bbox=dict(fc="white", ec="none", alpha=0.6, pad=0.15),
                clip_on=True,
            )

    @staticmethod
    def _add_damping_grid(
        ax: plt.Axes,
        zetas: tuple[float, ...] = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
        color: str = "0.35",
        alpha: float = 0.45,
        lw: float = 0.8,
    ) -> None:
        """
        Adds iso-damping grid to pole-zero map plot.

        Attributes
        ----------
        ax : plt.Axes
            Figure axes.
        zetas : tuple[float, ...], optional
            Damping to plot iso lines for 
            (default is (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)).
        color : str, optional
            Line color (default is "0.35").
        alpha : float, optional
            Line transparency (default is 0.45).
        lw : float, optional
            Linewidth (default is 0.8).
        """
                
        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()

        # Left half-plane only
        x_left = np.linspace(min(x_min, -1e-6), 0.0, 400)
        y_abs_max = max(abs(y_min), abs(y_max))

        for zeta in zetas:
            if not (0.0 < zeta < 1.0):
                continue

            # tan(theta), theta = arccos(zeta)
            m = np.sqrt(1.0 - zeta**2) / zeta
            y = m * (-x_left)

            mask = y <= 1.05 * y_abs_max
            if not np.any(mask):
                continue

            ax.plot(x_left[mask],  y[mask], ":", color=color, alpha=alpha, lw=lw, zorder=0)
            ax.plot(x_left[mask], -y[mask], ":", color=color, alpha=alpha, lw=lw, zorder=0)

            x_vis = x_left[mask]
            y_vis = y[mask]

            if x_vis.size > 0:
                # pick a point in the visible segment (closer to upper side)
                idx = int(0.25 * (x_vis.size - 1))
                x_lbl = x_vis[idx]
                y_lbl = y_vis[idx]

                ax.text(
                    x_lbl, y_lbl, fr"$\zeta={zeta:.1f}$",
                    fontsize=7, color=color, ha="left", va="bottom",
                    bbox=dict(fc="white", ec="none", alpha=0.6, pad=0.15),
                    clip_on=True,
                )

    def build_systems(
        self, 
        group: Axis, 
        coeff_values: Mapping[str, Iterable[float]]
    ) -> list[tuple[dict[str, float], Axis, LinearizedSystem]]:
        """
        Construct aircraft for all factor combinations in the sweep.

        Attributes
        ----------
        group : Axis
            Axis to analyze.
        coeff_values : Mapping[str, Iterable[float]]
            Factor values for each coefficient.

        Returns
        ----------
        list[tuple[dict[str, float], Axis, LinearizedSystem]]
            Factors, axis and aircraft for each case.
        """
                
        coeff_a, coeff_b = self.SUPPORTED[group]
        vals_a = [float(v) for v in coeff_values.get(coeff_a, [])]
        vals_b = [float(v) for v in coeff_values.get(coeff_b, [])]

        # No artificial full [1.0] arrays:
        # both provided -> full grid
        # only one provided -> sweep only that one
        # none provided -> single point
        if vals_a and vals_b:
            comb = ((a, b) for a, b in product(vals_a, vals_b))
        elif vals_a:
            comb = ((a, 1.0) for a in vals_a)
        elif vals_b:
            comb = ((1.0, b) for b in vals_b)
        else:
            comb = ((1.0, 1.0),)

        systems: list[tuple[dict[str, float], Axis, LinearizedSystem]] = []
        coeff_block_name = group.value

        for fa, fb in comb:
            factors = {coeff_a: fa, coeff_b: fb}
            ac = Aircraft(self.model)
            coeff_block = getattr(ac.stab_coeffs, coeff_block_name)

            for name, factor in factors.items():
                base = getattr(coeff_block, name)
                setattr(coeff_block, name, base * factor)

            self._recompute_derivatives(ac)
            systems.append((factors, group, LinearizedSystem(ac)))

        return systems

    def sweep(
        self, 
        group: Axis, 
        coeff_values: Mapping[str, Iterable[float]]
    ) -> list[SweepPoint]:
        """
        Obtain sweep points.

        Attributes
        ----------
        group : Axis
            Axis to analyze.
        coeff_values : Mapping[str, Iterable[float]]
            Factor values for each coefficient.

        Returns
        ----------
        list[SweepPoint]
            Sweep points.
        """

        points: list[SweepPoint] = []
        for factors, axis, lin_sys in self.build_systems(group, coeff_values):
            sys = lin_sys.get_sys(axis)
            poles = np.asarray(control.poles(sys), dtype=complex)
            points.append(SweepPoint(factors=factors, axis=axis, poles=poles))
        return points

    def plot_pzmap(
        self,
        group: Axis,
        coeff_values: Mapping[str, Iterable[float]],
        savefig: bool = False,
        showfig: bool = False,
    ) -> None:
        """
        Plot aircraft pole-zero map for each factor combination.

        Attributes
        ----------
        group : Axis
            Axis to analyze.
        coeff_values : Mapping[str, Iterable[float]]
            Factor values for each coefficient.
        savefig : bool, optional
            Option to save figure (default is False).
        showfig : bool, optional
            Option to show figure (default is False).
        """
                
        coeff_a, coeff_b = self.SUPPORTED[group]
        vals_a = [float(v) for v in coeff_values.get(coeff_a, [])]
        vals_b = [float(v) for v in coeff_values.get(coeff_b, [])]
        points = self.sweep(group, {coeff_a: vals_a, coeff_b: vals_b})

        varied = [name for name, vals in ((coeff_a, vals_a), (coeff_b, vals_b)) if vals]

        fig, ax = plt.subplots(figsize=(9, 6))
        if len(varied) <= 1:
            # One-parameter color map (or nominal single point)
            name = varied[0] if varied else coeff_a
            cvals = np.array([p.factors[name] for p in points], dtype=float)
            vmin, vmax = cvals.min(), cvals.max()
            if np.isclose(vmin, vmax):
                vmax = vmin + 1.0

            norm = plt.Normalize(vmin, vmax)
            cmap = plt.get_cmap("viridis")

            for p in points:
                cv = p.factors[name]
                ax.plot(p.poles.real, p.poles.imag, "x", color=cmap(norm(cv)), ms=7, mew=3)

            sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
            sm.set_array([])
            fig.colorbar(sm, ax=ax, pad=0.02, label=f"{name} factor")

        else:
            # Two-parameter bivariate color:
            # hue <- coeff_a, value/brightness <- coeff_b
            a_vals = np.array([p.factors[coeff_a] for p in points], dtype=float)
            b_vals = np.array([p.factors[coeff_b] for p in points], dtype=float)

            n1 = plt.Normalize(a_vals.min(), a_vals.max() if not np.isclose(a_vals.min(), a_vals.max()) else a_vals.min() + 1.0)
            n2 = plt.Normalize(b_vals.min(), b_vals.max() if not np.isclose(b_vals.min(), b_vals.max()) else b_vals.min() + 1.0)

            for p, va, vb in zip(points, a_vals, b_vals):
                h = n1(va)
                v = 0.30 + 0.70 * n2(vb)
                color = hsv_to_rgb((h, 0.90, v))
                ax.plot(p.poles.real, p.poles.imag, "x", color=color, ms=7, mew=3)

            # 2D legend inset
            h = np.linspace(0, 1, 200)
            v = np.linspace(0, 1, 200)
            H, V = np.meshgrid(h, v)
            legend_rgb = hsv_to_rgb(np.dstack((H, np.full_like(H, 0.90), 0.30 + 0.70 * V)))

            iax = ax.inset_axes([0.62, 0.06, 0.33, 0.33])
            iax.imshow(
                legend_rgb,
                origin="lower",
                aspect="auto",
                extent=[a_vals.min(), a_vals.max(), b_vals.min(), b_vals.max()],
            )
            iax.set_xlabel(f"{coeff_a} factor", fontsize=8)
            iax.set_ylabel(f"{coeff_b} factor", fontsize=8)
            iax.tick_params(labelsize=7)

        ax.axhline(0.0, color="k", lw=0.8, alpha=0.6)
        ax.axvline(0.0, color="k", lw=0.8, alpha=0.6)
        self._add_damping_grid(ax)
        self._add_omega_grid(ax)
        ax.grid(True, which="both", linestyle=":", alpha=0.5)
        ax.set_xlabel("Re [1/s]")
        ax.set_ylabel("Im [rad/s]")
        ax.set_title(f"Pole sensitivity ({group.value})")

        fig.tight_layout()
        if savefig:
            fig.savefig(self._figure_dir() / f"sensitivity_{group.value}_poles.png", dpi=180)
        if showfig:
            plt.show()
        else:
            plt.close(fig)

    def plot_sas_pzmap(
        self,
        points: list,  # list[SASDesignPoint]
        axis: Axis,
        gain_a: OutputChannel,
        gain_b: OutputChannel,
        savefig: bool = False,
        showfig: bool = False,
    ) -> None:
        """
        Plot Stability Augmentation System pole-zero map for each gain combination.

        Attributes
        ----------
        points : list
            Stability Augmentation System sweep point.
        axis : Axis
            Axis to analyze.
        gain_a : OutputChannel
            First feedback gain channel.
        gain_b : OutputChannel
            Second feedback gain channel.
        savefig : bool, optional
            Option to save figure (default is False).
        showfig : bool, optional
            Option to show figure (default is False).
        """

        fig, ax = plt.subplots(figsize=(9, 6))

        vals_a = np.array([p.feedback_gains[gain_a] for p in points], dtype=float)
        vals_b = np.array([p.feedback_gains[gain_b] for p in points], dtype=float)

        varied_a = len(np.unique(vals_a)) > 1
        varied_b = len(np.unique(vals_b)) > 1

        if varied_a and varied_b:
            n1 = plt.Normalize(vals_a.min(), vals_a.max() if not np.isclose(vals_a.min(), vals_a.max()) else vals_a.min() + 1.0)
            n2 = plt.Normalize(vals_b.min(), vals_b.max() if not np.isclose(vals_b.min(), vals_b.max()) else vals_b.min() + 1.0)
            for p, va, vb in zip(points, vals_a, vals_b):
                color = hsv_to_rgb((n1(va), 0.90, 0.30 + 0.70 * n2(vb)))
                ax.plot(p.poles.real, p.poles.imag, "x", color=color, ms=7, mew=3)
            # 2D legend inset
            h = np.linspace(0, 1, 200)
            v = np.linspace(0, 1, 200)
            H, V = np.meshgrid(h, v)
            legend_rgb = hsv_to_rgb(np.dstack((H, np.full_like(H, 0.90), 0.30 + 0.70 * V)))

            iax = ax.inset_axes([0.62, 0.06, 0.33, 0.33])
            iax.imshow(
                legend_rgb,
                origin="lower",
                aspect="auto",
                extent=[vals_a.min(), vals_a.max(), vals_a.min(), vals_a.max()],
            )
            iax.set_xlabel(f"{gain_a.value} gain", fontsize=8)
            iax.set_ylabel(f"{gain_b.value} gain", fontsize=8)
            iax.tick_params(labelsize=7)
        else:
            name = gain_a if varied_a else gain_b
            cvals = np.array([p.feedback_gains[name] for p in points], dtype=float)
            vmin, vmax = cvals.min(), cvals.max()
            if np.isclose(vmin, vmax):
                vmax = vmin + 1.0
            norm = plt.Normalize(vmin, vmax)
            cmap = plt.get_cmap("viridis")
            for p in points:
                ax.plot(p.poles.real, p.poles.imag, "x", color=cmap(norm(p.feedback_gains[name])), ms=7, mew=3)
            sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
            sm.set_array([])
            fig.colorbar(sm, ax=ax, pad=0.02, label=f"{name.value} gain")

        ax.axhline(0.0, color="k", lw=0.8, alpha=0.6)
        ax.axvline(0.0, color="k", lw=0.8, alpha=0.6)
        self._add_damping_grid(ax)
        self._add_omega_grid(ax)
        ax.grid(True, which="both", linestyle=":", alpha=0.5)
        ax.set_xlabel("Re [1/s]")
        ax.set_ylabel("Im [rad/s]")
        ax.set_title(f"SAS pole sensitivity ({axis.value})")

        fig.tight_layout()
        if savefig:
            fig.savefig(self._figure_dir() / f"sas_sensitivity_{axis.value}_poles.png", dpi=180)
        if showfig:
            plt.show()
        else:
            plt.close(fig)