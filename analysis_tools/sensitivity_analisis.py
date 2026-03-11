import numpy as np
import matplotlib as mpl
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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from itertools import product
from matplotlib.colors import rgb_to_hsv, hsv_to_rgb

mpl.rcParams["savefig.dpi"] = 300

# Plot styling (keep consistent across analysis scripts)
mpl.rcParams.update(
    {
        # LaTeX-like look without requiring a TeX install
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "DejaVu Serif", "Times New Roman"],
        "mathtext.fontset": "stix",
        "axes.labelsize": 12,
        "axes.titlesize": 14,
        "axes.titleweight": "normal",
        "figure.titlesize": 14,
        "figure.titleweight": "normal",
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
    }
)

# LaTeX-like labels for coefficients, gains, and channels
LATEX_LABELS = {
    # Stability coefficients
    "Cm_alpha": r"$C_{m_{\alpha}}$",
    "Cm_q": r"$C_{m_{q}}$",
    "Cm_delta_e": r"$C_{m_{\delta_e}}$",
    "Cn_beta": r"$C_{n_{\beta}}$",
    "Cn_r": r"$C_{n_{r}}$",
    "Cn_delta_r": r"$C_{n_{\delta_r}}$",
    "Cl_beta": r"$C_{l_{\beta}}$",
    "Cl_p": r"$C_{l_{p}}$",
    "Cl_r": r"$C_{l_{r}}$",
    "Cl_delta_a": r"$C_{l_{\delta_a}}$",
    "Cl_delta_r": r"$C_{l_{\delta_r}}$",
    # Output channels
    "alpha": r"$\alpha$",
    "beta": r"$\beta$",
    "theta": r"$\theta$",
    "phi": r"$\phi$",
    "p": r"$p$",
    "q": r"$q$",
    "r": r"$r$",
    "u": r"$u$",
    # PID gains
    "kp": r"$K_p$",
    "ki": r"$K_i$",
    "kd": r"$K_d$",
}


def _latex_label(name: str) -> str:
    return LATEX_LABELS.get(name, f"${name}$")


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
    model : str
        Aircraft model to analyze.
    fig_root : str | Path, optional
        Figure root saving path (default is "sens_study").

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
        Sensitivity analyzer class.

        Creates pole-zero map for aircraft, SAS and AP when sweeping parameters.

        Parameters
        ----------
        model : str
            Aircraft model to analyze.
        fig_root : str | Path, optional
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

        Parameters
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

        Parameters
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
        th = np.linspace(-np.pi / 2, np.pi / 2, 1000)

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
            ax.plot(x_vis, y_vis, "--", color=color, alpha=alpha, lw=lw, zorder=0)

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

        Parameters
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
        x_left = np.linspace(min(x_min, -1e-6), 0.0, 1000)
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

        Parameters
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

        Parameters
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

        Parameters
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
            cmap = plt.get_cmap("jet")

            for p in points:
                cv = p.factors[name]
                ax.plot(p.poles.real, p.poles.imag, "x", color=cmap(norm(cv)), ms=7, mew=3)

            # colorbar outside
            sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
            sm.set_array([])
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="3.5%", pad=0.12)
            cbar = fig.colorbar(sm, cax=cax, label=f"{_latex_label(name)} factor")
            cbar.set_label(f"{_latex_label(name)} factor", fontsize=12)
            cbar.ax.tick_params(labelsize=12)

        else:
            # Two-parameter bivariate color:
            # hue <- coeff_a, value/brightness <- coeff_b
            fig.subplots_adjust(right=0.78)
            a_vals = np.array([p.factors[coeff_a] for p in points], dtype=float)
            b_vals = np.array([p.factors[coeff_b] for p in points], dtype=float)

            cmap = plt.get_cmap("jet")
            n1 = plt.Normalize(a_vals.min(), a_vals.max())
            n2 = plt.Normalize(b_vals.min(), b_vals.max())

            # --- point colors (jet hue + brightness from b) ---
            for p, va, vb in zip(points, a_vals, b_vals):
                rgb = np.array(cmap(n1(va))[:3])
                h, s, _ = rgb_to_hsv(rgb.reshape(1, -1))[0]
                v = 0.30 + 0.70 * n2(vb)
                color = hsv_to_rgb([h, s, v])
                ax.plot(p.poles.real, p.poles.imag, "x", color=color, ms=7, mew=3)

            # --- legend outside: jet hue on x, brightness on y ---
            h = np.linspace(0, 1, 200)
            v = np.linspace(0.30, 1.0, 200)
            H, V = np.meshgrid(h, v)

            rgb = cmap(H)[..., :3]                    # jet colors along H
            hsv = rgb_to_hsv(rgb)
            hsv[..., 2] = V                           # replace value with brightness
            legend_rgb = hsv_to_rgb(hsv)

            fig.subplots_adjust(right=0.78)           # leave space
            bbox = ax.get_position()
            iax = fig.add_axes([bbox.x1 + 0.02, bbox.y0 + 0.05, 0.18, 0.28])
            iax.imshow(
                legend_rgb,
                origin="lower",
                aspect="auto",
                extent=[a_vals.min(), a_vals.max(), b_vals.min(), b_vals.max()],
            )
            iax.set_xlabel(f"{_latex_label(coeff_a)} factor", fontsize=12)
            iax.set_ylabel(f"{_latex_label(coeff_b)} factor", fontsize=12)
            iax.tick_params(labelsize=12)

        ax.axhline(0.0, color="k", lw=0.8, alpha=0.6)
        ax.axvline(0.0, color="k", lw=0.8, alpha=0.6)
        self._add_damping_grid(ax)
        self._add_omega_grid(ax)
        ax.grid(True, which="both", linestyle=":", alpha=0.5)
        ax.set_xlabel("Re [1/s]")
        ax.set_ylabel("Im [rad/s]")
        ax.set_title(f"Aircraft pole sensitivity ({group.value})")

        fig.tight_layout()
        if savefig:
            fig.savefig(self._figure_dir() / f"sensitivity_{group.value}_poles.png", dpi=300)
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

        Parameters
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
            cmap = plt.get_cmap("jet")
            n1 = plt.Normalize(vals_a.min(), vals_a.max())
            n2 = plt.Normalize(vals_b.min(), vals_b.max())

            # --- point colors (jet hue + brightness from b) ---
            for p, va, vb in zip(points, vals_a, vals_b):
                rgb = np.array(cmap(n1(va))[:3])
                h, s, _ = rgb_to_hsv(rgb.reshape(1, -1))[0]
                v = 0.30 + 0.70 * n2(vb)
                color = hsv_to_rgb([h, s, v])
                ax.plot(p.poles.real, p.poles.imag, "x", color=color, ms=7, mew=3)

            # --- legend outside: jet hue on x, brightness on y ---
            h = np.linspace(0, 1, 200)
            v = np.linspace(0.30, 1.0, 200)
            H, V = np.meshgrid(h, v)

            rgb = cmap(H)[..., :3]                    # jet colors along H
            hsv = rgb_to_hsv(rgb)
            hsv[..., 2] = V                           # replace value with brightness
            legend_rgb = hsv_to_rgb(hsv)

            fig.subplots_adjust(right=0.78)           # leave space
            bbox = ax.get_position()
            iax = fig.add_axes([bbox.x1 + 0.02, bbox.y0 + 0.05, 0.18, 0.28])
            iax.imshow(
                legend_rgb,
                origin="lower",
                aspect="auto",
                extent=[vals_a.min(), vals_a.max(), vals_b.min(), vals_b.max()],
            )
            iax.set_xlabel(f"{_latex_label(gain_a.value)} gain", fontsize=12)
            iax.set_ylabel(f"{_latex_label(gain_b.value)} gain", fontsize=12)
            iax.tick_params(labelsize=12)
        else:
            name = gain_a if varied_a else gain_b
            cvals = np.array([p.feedback_gains[name] for p in points], dtype=float)
            vmin, vmax = cvals.min(), cvals.max()
            if np.isclose(vmin, vmax):
                vmax = vmin + 1.0
            norm = plt.Normalize(vmin, vmax)
            cmap = plt.get_cmap("jet")
            for p in points:
                ax.plot(p.poles.real, p.poles.imag, "x", color=cmap(norm(p.feedback_gains[name])), ms=7, mew=3)
            sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, pad=0.02, label=f"{_latex_label(name.value)} gain")
            cbar.set_label(f"{_latex_label(name.value)} gain", fontsize=12)
            cbar.ax.tick_params(labelsize=12)

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
            fig.savefig(self._figure_dir() / f"sas_sensitivity_{axis.value}_poles.png", dpi=300)
        if showfig:
            plt.show()
        else:
            plt.close(fig)

    def plot_ap_pzmap(
        self,
        points: list,  # list[APDesignPoint]
        axis: Axis,
        gain_a: str = "kp",
        gain_b: str = "ki",
        gain_c: str = "kd",
        savefig: bool = False,
        showfig: bool = False,
    ) -> None:
        """
        Plot Autopilot pole-zero map for each PID gain combination.

        Parameters
        ----------
        points : list
            Autopilot sweep point.
        axis : Axis
            Axis to analyze.
        gain_a : str
            Proportional gain.
        gain_b : str
            Integral gain.
        gain_c : str
            Derivative gain.
        savefig : bool, optional
            Option to save figure (default is False).
        showfig : bool, optional
            Option to show figure (default is False).
        """
                
        allowed = {"kp", "ki", "kd"}
        if gain_a not in allowed or gain_b not in allowed or gain_c not in allowed:
            raise ValueError(f"gain names must be in {allowed}")
        if len({gain_a, gain_b, gain_c}) != 3:
            raise ValueError("gain_a, gain_b and gain_c must be different")
        if not points:
            raise ValueError("points is empty")

        vals_a = np.array([getattr(p.pid_gains, gain_a) for p in points], dtype=float)
        vals_b = np.array([getattr(p.pid_gains, gain_b) for p in points], dtype=float)
        vals_c = np.array([getattr(p.pid_gains, gain_c) for p in points], dtype=float)

        def _norm(vals: np.ndarray) -> plt.Normalize:
            vmin = float(vals.min())
            vmax = float(vals.max())
            if np.isclose(vmin, vmax):
                vmax = vmin + 1.0
            return plt.Normalize(vmin, vmax)

        var_a = len(np.unique(vals_a)) > 1
        var_b = len(np.unique(vals_b)) > 1
        var_c = len(np.unique(vals_c)) > 1
        varied = [(gain_a, vals_a, var_a), (gain_b, vals_b, var_b), (gain_c, vals_c, var_c)]
        varied = [v for v in varied if v[2]]

        # If only one gain is swept, show a 1D colorbar instead of a 2D colormap
        if len(varied) == 1:
            name, vals, _ = varied[0]
            fig, ax = plt.subplots(figsize=(9, 6))
            cmap = plt.get_cmap("jet")
            norm = _norm(vals)

            for p, v in zip(points, vals):
                ax.plot(p.poles.real, p.poles.imag, "x", color=cmap(norm(v)), ms=7, mew=3)

            sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
            sm.set_array([])
            cbar = fig.colorbar(sm, ax=ax, pad=0.02, label=f"{_latex_label(name)} gain")
            cbar.set_label(f"{_latex_label(name)} gain", fontsize=12)
            cbar.ax.tick_params(labelsize=12)

            ax.axhline(0.0, color="k", lw=0.8, alpha=0.6)
            ax.axvline(0.0, color="k", lw=0.8, alpha=0.6)
            self._add_damping_grid(ax)
            self._add_omega_grid(ax)
            ax.grid(True, which="both", linestyle=":", alpha=0.5)
            ax.set_xlabel("Re [1/s]")
            ax.set_ylabel("Im [rad/s]")
            ax.set_title(f"AP PID pole sensitivity ({axis.value})")

            fig.tight_layout()
            if savefig:
                fig.savefig(self._figure_dir() / f"ap_pid_sensitivity_{axis.value}_poles.png", dpi=300)
            if showfig:
                plt.show()
            else:
                plt.close(fig)
            return

        n1 = _norm(vals_a)  # hue <- gain_a
        # value/brightness <- gain_b (force mid if constant)
        if np.isclose(vals_b.min(), vals_b.max()):
            n2 = lambda _: 0.5
        else:
            n2 = _norm(vals_b)

        c_unique = np.unique(vals_c)
        n_pan = len(c_unique)
        ncols = min(4, n_pan)
        nrows = int(np.ceil(n_pan / ncols))

        fig, axs = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            figsize=(4.8 * ncols, 4.2 * nrows),
            sharex=True,
            sharey=True,
        )
        axs = np.atleast_1d(axs).ravel()

        for k, cval in enumerate(c_unique):
            ax = axs[k]
            mask = np.isclose(vals_c, cval)

            for p, va, vb, m in zip(points, vals_a, vals_b, mask):
                if not m:
                    continue
                cmap = plt.get_cmap("jet")

                rgb = np.array(cmap(n1(va))[:3])
                h, s, _ = rgb_to_hsv(rgb.reshape(1, -1))[0]
                v = 0.30 + 0.70 * n2(vb)
                color = hsv_to_rgb([h, s, v])

                ax.plot(p.poles.real, p.poles.imag, "x", color=color, ms=7, mew=3)

            ax.axhline(0.0, color="k", lw=0.8, alpha=0.6)
            ax.axvline(0.0, color="k", lw=0.8, alpha=0.6)
            self._add_damping_grid(ax)
            self._add_omega_grid(ax)
            ax.grid(True, which="both", linestyle=":", alpha=0.5)
            ax.set_title(f"{gain_c} = {cval:g}")
            ax.set_xlabel("Re [1/s]")
            ax.set_ylabel("Im [rad/s]")

        # Hide unused subplot slots
        for j in range(n_pan, len(axs)):
            axs[j].set_visible(False)

        # 2D legend for gain_a/gain_b on first active axis
        a_min, a_max = float(vals_a.min()), float(vals_a.max())
        b_min, b_max = float(vals_b.min()), float(vals_b.max())
        if np.isclose(a_min, a_max):
            a_max = a_min + 1.0
        if np.isclose(b_min, b_max):
            b_max = b_min + 1.0

        h = np.linspace(0, 1, 200)
        v = np.linspace(0.30, 1.0, 200)
        H, V = np.meshgrid(h, v)

        rgb = cmap(H)[..., :3]      # jet colors along H
        hsv = rgb_to_hsv(rgb)
        hsv[..., 2] = V             # replace value with brightness
        legend_rgb = hsv_to_rgb(hsv)

        # reserve space on the right
        fig.subplots_adjust(right=0.78)

        bbox = axs[0].get_position()
        iax = fig.add_axes([bbox.x1 + 0.02, bbox.y0 + 0.05, 0.18, 0.28])
        iax.imshow(
            legend_rgb,
            origin="lower",
            aspect="auto",
            extent=[a_min, a_max, b_min, b_max],
        )
        iax.set_xlabel(_latex_label(gain_a), fontsize=12)
        iax.set_ylabel(_latex_label(gain_b), fontsize=12)
        iax.tick_params(labelsize=12)

        fig.suptitle(f"AP PID pole sensitivity ({axis.value})")
        fig.tight_layout()

        if savefig:
            fig.savefig(self._figure_dir() / f"ap_pid_sensitivity_{axis.value}_poles.png", dpi=300)
        if showfig:
            plt.show()
        else:
            plt.close(fig)
