from pathlib import Path
import numpy as np
import control
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from flight_control_system.state_space import LinearizedSystem
from flight_control_system.axis_metadata import AXIS_LABELS, AxisLabels
from flight_control_system.types import Axis, InputChannel, OutputChannel
from collections.abc import Mapping

class FrequencyAnalyzer:
    """
    Frequency analyzer class.

    Creates Bode and Nichols plots for aircraft, components, SAS and AP.

    Attributes
    ----------
    AXES : tuple
        Dynamic system axes.
    sys : LinearizedSystem
        Aircraft linearized system.
    fig_root : str | Path
        Figure root saving path (default is "flight_dynamics").
    labels : dict[Axis, AxisLabels]
        Axis input and output (state) channels labels.
    bode_lims : dict[Axis, list]
        Bode plot x-axis limits in rad/s.

    Methods
    ----------
    _labels()
        Obtain axis labels.
    _clean_label()
        Clean labels for figure naming.
    _figure_dir()
        Construct figure saving directory.
    channel_metrics()
        Calculate frequency-domain metrics for channel y_i / u_j.
    plot_all()
        Plot system's Bode and Nichols plots for all axes.
    plot_freq_analysis()
        Plot system's Bode and Nichols plots for the selected axis.
    plot_bode()
        Plot system's Bode plots for the selected axis.
    plot_nichols()
        Plot system's Nichols plots for the selected axis.
    compare_components()
        Compare components frequency response in a Bode plot.
    compare_sas_nichols()
        Compare Stability Augmentation System Nichols plot when varying feedback gains.
    """
        
    AXES = (Axis.LONG, Axis.LATDIR)

    def __init__(self, lin_sys: LinearizedSystem, fig_root: str | Path = "flight_dynamics"):
        """
        Frequency analyzer class.

        Creates Bode and Nichols plots for aircraft, components, SAS and AP.

        Parameters
        ----------
        lin_sys : LinearizedSystem
            Dynamic system axes.
        fig_root : str | Path
            Figure root saving path (default is "flight_dynamics").
        """
            
        self.sys = lin_sys
        self.fig_root = Path(fig_root)
        self.labels = AXIS_LABELS
        self.bode_lims = {
            Axis.LONG: [1e-2, 1e2],
            Axis.LATDIR: [1e-6, 1e2],
        }

    def _labels(self, axis: Axis) -> AxisLabels:
        """
        Obtain axis labels.

        Parameters
        ----------
        axis : Axis
            Axis to get labels for.

        Returns
        ----------
        AxisLabels
            Axis labels.
        """

        return self.labels[axis]
    
    @staticmethod
    def _clean_label(label: str) -> str:
        """
        Clean labels for figure naming.

        Parameters
        ----------
        label : str
            Label to clean.

        Returns
        ----------
        str
            Cleaned label.
        """

        return label.replace("\\Delta", "").replace("\\", "").replace(" ", "")
    
    def _figure_dir(self) -> Path:
        """
        Construct figure saving directory.

        Returns
        ----------
        Path
            Figure saving directory.
        """

        figdir = self.fig_root / self.sys.aircraft.model
        figdir.mkdir(parents=True, exist_ok=True)
        return figdir
    
    def channel_metrics(self, axis: Axis, i: int, j: int) -> dict:
        """
        Calculate frequency-domain metrics for channel y_i / u_j.

        Parameters
        ----------
        axis : Axis
            Axis to analyze.
        i : int
            Output/State index to analyze.
        j : int
            Input index to analyze.

        Returns
        ----------
        dict
            Frequency-domain metrics for channel y_i / u_j.
        """

        gij = self.sys.get_sys(axis)[i, j]

        poles = np.asarray(control.poles(gij), dtype=complex)
        zeros = np.asarray(control.zeros(gij), dtype=complex)

        gm, pm, wgc, wpc = control.margin(gij)

        try:
            dc_gain = float(np.real_if_close(control.dcgain(gij)))
        except Exception:
            dc_gain = np.nan

        try:
            bandwidth = float(np.real_if_close(control.bandwidth(gij)))
        except Exception:
            bandwidth = np.nan

        return {
            "axis": axis,
            "output_index": i,
            "input_index": j,
            "poles": poles,
            "zeros": zeros,
            "gain_margin": gm,
            "phase_margin_deg": pm,
            "gain_crossover_rad_s": wgc,
            "phase_crossover_rad_s": wpc,
            "dc_gain": dc_gain,
            "bandwidth_rad_s": bandwidth,
        }
    
    def plot_all(self, savefig: bool = False, showfig: bool = False) -> None:
        """
        Plot system's Bode and Nichols plots for all axes.

        Parameters
        ----------
        savefig : bool, optional
            Option to save figure (default is False).
        showfig : bool, optional
            Option to show figure (default is False).
        """

        for axis in self.AXES:
            self.plot_freq_analysis(axis=axis, savefig=savefig, showfig=showfig)

    def plot_freq_analysis(self, axis: Axis, savefig: bool = False, showfig: bool = False) -> None:
        """
        Plot system's Bode and Nichols plots for the selected axis.

        Parameters
        ----------
        axis : Axis
            Axis to analyze.
        savefig : bool, optional
            Option to save figure (default is False).
        showfig : bool, optional
            Option to show figure (default is False).
        """

        # Define dynamic system
        sys = self.sys.get_sys(axis)
        labels = self._labels(axis)
        # Loop through every input and output channel
        for i, state_l in enumerate(labels.states):
            for j, input_l in enumerate(labels.inputs):
                # Figure layout
                fig = plt.figure(figsize=(14, 8))
                gs = gridspec.GridSpec(2, 2)

                ax_mag = fig.add_subplot(gs[0, 0])
                ax_phase = fig.add_subplot(gs[1, 0], sharex=ax_mag)
                ax_nichols = fig.add_subplot(gs[:, 1])

                # Bode plot
                control.bode_plot(sys[i, j], omega_limits=self.bode_lims[axis], dB=True, ax=[ax_mag, ax_phase])

                # Nichols plot
                control.nichols_plot(sys[i, j], omega=self.bode_lims[axis], ax=ax_nichols)
                ax_nichols.set_ylim(ax_mag.get_ylim())
                # Replot to convert phase (y axis) to positive if necessary
                if ax_phase.get_ylim()[1] < -170:
                    plt.close(fig)
                    fig = plt.figure(figsize=(14, 8))
                    gs = gridspec.GridSpec(2, 2)
                    ax_mag = fig.add_subplot(gs[0, 0])
                    ax_phase = fig.add_subplot(gs[1, 0], sharex=ax_mag)
                    ax_nichols = fig.add_subplot(gs[:, 1])

                    mag, phase, omega = control.frequency_response(sys[i, j], omega_limits=self.bode_lims[axis])
                    mag = np.squeeze(mag)
                    phase = np.squeeze(phase)

                    ax_mag.semilogx(omega, 20 * np.log10(mag))
                    ax_mag.set_ylabel("Magnitude [dB]")
                    ax_mag.grid(True, which="both")
                    ax_mag.set_xlim(self.bode_lims[axis])

                    ax_phase.semilogx(omega, np.unwrap(np.degrees(phase), discont=8*np.pi))
                    ax_phase.set_ylabel("Phase [deg]")
                    ax_phase.set_xlabel("Frequency [rad/s]")
                    ax_phase.grid(True, which="both")

                    control.nichols_plot(sys[i, j], omega=self.bode_lims[axis], ax=ax_nichols)
                    ax_nichols.set_ylim(ax_mag.get_ylim())
                    min_phase = ax_phase.get_ylim()[0]
                    max_phase = ax_phase.get_ylim()[1]
                    ax_nichols.set_xlim((min_phase, max_phase))
                else:
                    ax_nichols.set_xlim(ax_phase.get_ylim())
                ax_mag.set_title("Bode")
                ax_nichols.set_title("Nichols")
                plt.tight_layout()

                # Save/show figure
                if savefig:
                    input_var = self._clean_label(input_l)
                    state_var = self._clean_label(state_l)
                    figname = self.sys.aircraft.model + f'_{state_var}_{input_var}_Bode_Nichols.png'
                    fig.savefig(self._figure_dir() / figname)
                if showfig:
                    plt.show()
                else:
                    plt.close(fig)

    def plot_bode(self, axis: Axis, savefig: bool = False, showfig: bool = False) -> None:
        """
        Plot system's Bode plots for the selected axis.

        Parameters
        ----------
        axis : Axis
            Axis to analyze.
        savefig : bool, optional
            Option to save figure (default is False).
        showfig : bool, optional
            Option to show figure (default is False).
        """

        # Define dynamic system
        sys = self.sys.get_sys(axis)
        labels = self._labels(axis)
        # Loop through every input and output channel
        for i, state_l in enumerate(labels.states):
            for j, input_l in enumerate(labels.inputs):
                # Figure layout
                fig = plt.figure()
                control.bode_plot(sys[i, j], omega_limits=self.bode_lims[axis], dB=True)
                plt.suptitle(fr'Bode diagram of $[G(i\omega)]_{{{state_l+input_l}}} = \frac{{{state_l}}}{{{input_l}}}$')
                plt.tight_layout()

                # Save/show figure
                if savefig:
                    input_var = self._clean_label(input_l)
                    state_var = self._clean_label(state_l)
                    figname = self.sys.aircraft.model + f'_{state_var}_{input_var}_Bode.png'
                    fig.savefig(self._figure_dir() / figname)
                if showfig:
                    plt.show()
                else:
                    plt.close(fig)
    
    def plot_nichols(self, axis: Axis, savefig: bool = False, showfig: bool = False) -> None:
        """
        Plot system's Nichols plots for the selected axis.

        Parameters
        ----------
        axis : Axis
            Axis to analyze.
        savefig : bool, optional
            Option to save figure (default is False).
        showfig : bool, optional
            Option to show figure (default is False).
        """

        # Define dynamic system
        sys = self.sys.get_sys(axis)
        labels = self._labels(axis)
        # Loop through every input and output channel
        for i, state_l in enumerate(labels.states):
            for j, input_l in enumerate(labels.inputs):
                # Figure layout
                fig = plt.figure()
                control.nichols_plot(sys[i, j], omega=self.bode_lims[axis])
                plt.suptitle(fr'Nichols diagram of $[G(i\omega)]_{{{state_l+input_l}}} = \frac{{{state_l}}}{{{input_l}}}$')
                plt.tight_layout()

                # Save/show figure
                if savefig:
                    input_var = self._clean_label(input_l)
                    state_var = self._clean_label(state_l)
                    figname = self.sys.aircraft.model + f'_{state_var}_{input_var}_Nichols.png'
                    fig.savefig(self._figure_dir() / figname)
                if showfig:
                    plt.show()
                else:
                    plt.close(fig)

    @staticmethod
    def compare_components(
        tf_map: Mapping[str, control.TransferFunction],
        omega_limits: tuple[float, float] = (1e-1, 1e3),
        title: str = "Actuator Bode Comparison", # TODO: Change with component type
        showfig: bool = False,
        savefig: bool = False,
        save_path: str | Path | None = None,
    ) -> None:
        """
        Compare components frequency response in a Bode plot.

        Parameters
        ----------
        tf_map : Mapping[str, control.TransferFunction]
            Transfer functions of components.
        omega_limits : tuple[float, float], optional
            Bode plot x-axis limits in rad/s (default is (1e-1, 1e3)).
        title : str, optional
            Figure title (default is "Actuator Bode Comparison").
        savefig : bool, optional
            Option to save figure (default is False).
        showfig : bool, optional
            Option to show figure (default is False).
        save_path : str | Path | None, optional
            Figure saving path (default is None).
        """
                
        if not tf_map:
            raise ValueError("tf_map is empty.")
        else:
            systems = list(tf_map.values())
            labels = list(tf_map.keys())

        fig, (ax_mag, ax_phase) = plt.subplots(2, 1, sharex=True)
        control.bode_plot(systems, omega_limits=omega_limits, dB=True, deg=True, ax=[ax_mag, ax_phase])

        # Attach labels to plotted lines
        for line, lbl in zip(ax_mag.get_lines(), labels):
            line.set_label(lbl)
        for line, lbl in zip(ax_phase.get_lines(), labels):
            line.set_label(lbl)

        fig.suptitle(title)
        ax_mag.legend(loc="best")
        ax_phase.legend(loc="best")
        plt.tight_layout()

        if savefig and save_path is not None:
            figname = f'{title.replace(" ", "_")}.png'
            save_path = Path(save_path)
            save_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path / figname)
        if showfig:
            plt.show()
        else:
            plt.close(fig)

    @staticmethod
    def compare_sas_nichols(
        sys_map: Mapping[str, control.StateSpace],
        axis: Axis,
        out_ch: OutputChannel,
        in_ch: InputChannel,
        omega_limits: tuple[float, float] = (1e-1, 1e3),
        title: str = "SAS Nichols plot sensitivity analysis",
        showfig: bool = False,
        savefig: bool = False,
        save_path: str | Path | None = None,
    ) -> None:
        """
        Compare Stability Augmentation System Nichols plot when varying feedback gains.

        Parameters
        ----------
        sys_map : Mapping[str, control.StateSpace]
            Stability augmentation systems to compare.
        axis : Axis
            Axis to analyze.
        out_ch : OutputChannel
            Output/State channel to analyze.
        in_ch : InputChannel
            Input channel to analyze.
        omega_limits : tuple[float, float], optional
            Bode plot x-axis limits in rad/s (default is (1e-1, 1e3)).
        title : str, optional
            Figure title (default is "SAS Nichols plot sensitivity analysis").
        savefig : bool, optional
            Option to save figure (default is False).
        showfig : bool, optional
            Option to show figure (default is False).
        save_path : str | Path | None, optional
            Figure saving path (default is None).
        """
                
        if not sys_map:
            raise ValueError("sys_map is empty.")

        i = AXIS_LABELS[axis].state_channels.index(out_ch)
        j = AXIS_LABELS[axis].input_channels.index(in_ch)
        omega = np.logspace(np.log10(omega_limits[0]), np.log10(omega_limits[1]), 10000)

        fig, ax = plt.subplots(figsize=(8, 6))
        control.nichols_grid()  # Draw grid first

        for name, sys in sys_map.items():
            n0 = len(ax.get_lines())
            control.nichols_plot(sys[i, j], omega=omega, ax=ax)
            new_lines = ax.get_lines()[n0:]

            if not new_lines:
                continue

            # One legend entry per system
            new_lines[0].set_label(name)
            for ln in new_lines[1:]:
                ln.set_label("_nolegend_")

        fig.suptitle(title)
        ax.legend(loc="best")
        ax.set_title(f"Nichols comparison: {out_ch.value}/{in_ch.value}")
        ax.grid(True, which="both", linestyle=":", alpha=0.5)

        if savefig and save_path is not None:
            figname = f'{title.replace(" ", "_")}.png'
            save_path = Path(save_path)
            save_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path / figname)
        if showfig:
            plt.show()
        else:
            plt.close(fig)

    @staticmethod
    def compare_ap_nichols(
        sys_map: Mapping[str, control.StateSpace],
        axis: Axis,
        out_ch: OutputChannel,
        in_ch: InputChannel,
        omega_limits: tuple[float, float] = (1e-1, 1e3),
        title: str = "AP Nichols plot sensitivity analysis",
        showfig: bool = False,
        savefig: bool = False,
        save_path: str | Path | None = None,
    ) -> None:
        """
        Compare Autopilot Nichols plot when varying PID gains.

        Parameters
        ----------
        sys_map : Mapping[str, control.StateSpace]
            Autopilots to compare.
        axis : Axis
            Axis to analyze.
        out_ch : OutputChannel
            Output/State channel to analyze.
        in_ch : InputChannel
            Input channel to analyze.
        omega_limits : tuple[float, float], optional
            Bode plot x-axis limits in rad/s (default is (1e-1, 1e3)).
        title : str, optional
            Figure title (default is "AP Nichols plot sensitivity analysis").
        savefig : bool, optional
            Option to save figure (default is False).
        showfig : bool, optional
            Option to show figure (default is False).
        save_path : str | Path | None, optional
            Figure saving path (default is None).
        """
                
        if not sys_map:
            raise ValueError("sys_map is empty.")

        i = AXIS_LABELS[axis].state_channels.index(out_ch)
        j = AXIS_LABELS[axis].input_channels.index(in_ch)
        omega = np.logspace(np.log10(omega_limits[0]), np.log10(omega_limits[1]), 10000)

        fig, ax = plt.subplots(figsize=(8, 6))
        control.nichols_grid()  # Draw grid first

        for name, sys in sys_map.items():
            n0 = len(ax.get_lines())
            control.nichols_plot(sys[i, j], omega=omega, ax=ax)
            new_lines = ax.get_lines()[n0:]

            if not new_lines:
                continue

            # One legend entry per system
            new_lines[0].set_label(name)
            for ln in new_lines[1:]:
                ln.set_label("_nolegend_")

        fig.suptitle(title)
        ax.legend(loc="best")
        ax.set_title(f"Nichols comparison: {out_ch.value}/{in_ch.value}")
        ax.grid(True, which="both", linestyle=":", alpha=0.5)

        if savefig and save_path is not None:
            figname = f'{title.replace(" ", "_")}.png'
            save_path = Path(save_path)
            save_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path / figname)
        if showfig:
            plt.show()
        else:
            plt.close(fig)