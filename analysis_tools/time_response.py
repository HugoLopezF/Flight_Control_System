import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import control
from .input_signals import InputSignal, build_input
from flight_control_system.types import Axis, InputChannel
from flight_control_system.axis_metadata import AXIS_LABELS
from flight_control_system.state_space import LinearizedSystem
from collections.abc import Mapping


class TimeResponseAnalyzer:
    AXES = (Axis.LONG, Axis.LATDIR)

    def __init__(self, lin_sys: LinearizedSystem, fig_root: str | Path = "flight_dynamics"):
        self.sys = lin_sys
        self.fig_root = Path(fig_root)
        self.labels = AXIS_LABELS

    def _labels(self, axis: Axis):
        """
        Obtain axis labels.

        :param axis: Axis to get labels for.
        :type axis: Axis
        """
        return self.labels[axis]
    
    def _figure_dir(self) -> Path:
        """
        Construct figure saving directory.
        """
        figdir = self.fig_root / self.sys.aircraft.model
        figdir.mkdir(parents=True, exist_ok=True)
        return figdir
    
    @staticmethod
    def get_component_response(tf: control.TransferFunction | control.StateSpace, t: np.ndarray, input_type: InputSignal, amp: float = 1.0, t_end: float = 10.0) -> tuple[control.StateSpace, np.ndarray, np.ndarray, np.ndarray]:
        """
        Obtain component's response for selected input.

        :param tf: Transfer function to obtain response from
        :type tf: control.TransferFunction | control.StateSpace
        :param t: Time vector.
        :type t: np.ndarray
        :param input_type: Input signal type.
        :type input_type: InputSignal
        :param amp: Input signal amplitude (nondimensional).
        :type amp: float
        :param t_end: Time instant at which the ramp stops (if ramp).
        :type t_end: float

        :rtype: tuple[control.StateSpace, np.ndarray, np.ndarray, np.ndarray]
        """
        # Fill input channels
        u = np.zeros((1, len(t)))
        amp_rad = np.deg2rad(amp)
        u = build_input(signal=input_type, t=t, amp=amp_rad, t_end=t_end)

        # Calculate response
        t_out, y_out = control.forced_response(tf, t, u)
        # Convert response back to degrees except for u channel (m/(s*rad))
        y_out *= 180 / np.pi
        u *= 180 / np.pi

        return tf, t_out, y_out, u
    
    @staticmethod
    def compare_components(
        tf_map: Mapping[str, control.TransferFunction],
        t: np.ndarray,
        input_type: InputSignal, 
        amp: float = 1.0, 
        t_end: float = 10.0,
        title: str = "Actuator response Comparison",
        showfig: bool = False,
        savefig: bool = False,
        save_path: str | Path | None = None,
    ):
        """

        """

        if not tf_map:
            raise ValueError("tf_map is empty.")

        fig = plt.figure()
        fig.suptitle(f'Component\'s response  to a {amp}deg {input_type.value}', fontsize=13)
        for label, tf in tf_map.items():
            tf, t_out, y_out, u = TimeResponseAnalyzer.get_component_response(tf=tf, t=t, input_type=input_type, amp=amp, t_end=t_end)
            plt.plot(t_out, y_out, label=label)
        plt.plot(t_out, u, label='Input')
        plt.grid(which='both')
        plt.ylabel('$y(t)$ [deg]')
        plt.xlabel('$t$ [s]')
        plt.legend(loc="best")
        plt.tight_layout()

        # Save/show plot
        if savefig and save_path is not None:
            figname = f'{title.replace(" ", "_")}.png'
            save_path = Path(save_path)
            save_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path / figname)
        if showfig:
            plt.show()
        else:
            plt.close(fig)

    def get_aircraft_response(self, t: np.ndarray, axis: Axis, channel: InputChannel, input_type: InputSignal, amp: float = 1.0, t_end: float = 10.0) -> tuple[control.StateSpace, np.ndarray, np.ndarray, np.ndarray]:
        """
        Obtain aircraft's response for selected input and axis.

        :param t: Time vector.
        :type t: np.ndarray
        :param axis: Axis to get response from.
        :type axis: Axis
        :param channel: Input channel for this axis.
        :type channel: InputChannel
        :param input_type: Input signal type.
        :type input_type: InputSignal
        :param amp: Input signal amplitude (nondimensional).
        :type amp: float
        :param t_end: Time instant at which the ramp stops (if ramp).
        :type t_end: float

        :rtype: tuple[control.StateSpace, np.ndarray, np.ndarray, np.ndarray]
        """
        # Get dynamic system from matrices
        sys = self.sys.get_sys(axis)
        labels = self._labels(axis)

        # Check channel and convert input to rad
        if channel not in labels.input_channels:
            allowed = ", ".join(ch.value for ch in labels.input_channels)
            raise ValueError(f"Invalid channel '{channel.value}' for axis '{axis.value}'. Allowed: {allowed}")
        
        # Fill input channels
        u = np.zeros((sys.ninputs, len(t)))
        amp_rad = np.deg2rad(amp)
        idx = labels.input_channels.index(channel)
        u[idx, :] = build_input(signal=input_type, t=t, amp=amp_rad, t_end=t_end)

        # Calculate response
        t_out, y_out = control.forced_response(sys, t, u)
        # Convert response back to degrees except for u channel (m/(s*rad))
        if axis is Axis.LONG:
            y_out[1:, :] *= 180 / np.pi
        else:
            y_out *= 180 / np.pi

        return sys, t_out, y_out, u
        
    def plot_aircraft_response(self, t: np.ndarray, axis: Axis, channel: InputChannel, input_type: InputSignal, amp: float = 1.0, t_end: float = 10.0, savefig: bool = False, showfig: bool = False) -> None:
        """
        Plot aircraft's response for selected input type and axis.

        :param t: Time vector.
        :type t: np.ndarray
        :param axis: Axis to get response from.
        :type axis: Axis
        :param channel: Input channel for this axis.
        :type channel: InputChannel
        :param input_type: Input signal type.
        :type input_type: InputSignal
        :param amp: Input signal amplitude (nondimensional).
        :type amp: float
        :param t_end: Time instant at which the ramp stops (if ramp).
        :type t_end: float
        :param savefig: Option to save figure.
        :type savefig: bool
        :param showfig: Option to show figure.
        :type showfig: bool
        """
        # Get system response
        sys, t_out, y_out, u = self.get_aircraft_response(t=t, axis=axis, amp=amp, t_end=t_end, channel=channel, input_type=input_type)
        # Reformat input/output labels for plotting
        labels = self._labels(axis)
        sys.input_labels = [f'${var}$' for var in labels.inputs]
        sys.output_labels = [f'${var}$' for var in labels.states]

        # Plot response
        fig, axes = plt.subplots(len(sys.output_labels), 1, sharex=True)
        fig.suptitle(f'System response for {axis.value} axis to a {amp}deg {input_type.value} in {channel.value}', fontsize=13)
        for i, ax in enumerate(axes):
            y = y_out[i, :]
            ax.plot(t_out, y)

            # Steady-state line (mean of last 10% of samples)
            n_tail = max(1, int(0.1 * len(y)))
            y_ss = float(np.mean(y[-n_tail:]))
            ax.axhline(y_ss, color="red", linestyle="--", linewidth=1.2, label=f"ss = {y_ss:.3g}")

            ax.set_ylabel(sys.output_labels[i])

            ax.minorticks_on()
            ax.grid(True, which="major", linestyle="-", alpha=0.5)
            ax.grid(True, which="minor", linestyle=":", alpha=0.35)

            ax.legend(loc="best")
        axes[-1].set_xlabel('$t$ [s]')

        # Save/show plot
        if savefig:
            figname = self.sys.aircraft.model + f'{amp}_deg_{input_type.value}_{channel.value}_response.png'
            fig.savefig(self._figure_dir() / figname)
        if showfig:
            plt.show()
        plt.close(fig)

    def plot_full_response(self, t: np.ndarray, input_type: InputSignal, amp: float = 1.0, t_end: float = 10.0, savefig: bool = False, showfig: bool = False) -> None:
        """
        Plot aircraft's response for selected input type for all axes.


        :param t: Time vector.
        :type t: np.ndarray
        :param input_type: Input signal type.
        :type input_type: InputSignal
        :param amp: Input signal amplitude (nondimensional).
        :type amp: float
        :param t_end: Time instant at which the ramp stops (if ramp).
        :type t_end: float
        :param savefig: Option to save figure.
        :type savefig: bool
        :param showfig: Option to show figure.
        :type showfig: bool
        """
        # Loop through all axes and input channels
        for axis in self.AXES:        
            labels = self._labels(axis)
            for ch in labels.input_channels:
                self.plot_aircraft_response(t, axis=axis, input_type=input_type, amp=amp, t_end=t_end, channel=ch, savefig=savefig, showfig=showfig)