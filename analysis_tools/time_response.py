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
    """
    Time response analyzer class.

    Creates time response plots for aircraft, components, SAS and AP.

    Attributes
    ----------
    AXES : dict[Axis, tuple(str, str)]
        Dynamic system axes.
    sys : LinearizedSystem
        Linearized system to analyze.
    fig_root : str | Path
        Figure root saving path (default is "flight_dynamics").
    labels : dict[Axis, AxisLabels]
        Axis input and output (state) channels labels.

    Methods
    ----------
    _labels()
        Construct figure saving directory.
    _figure_dir()
        Recompute stability derivatives for updated stability coefficients.
    get_component_response()
        Obtain component's response for selected input.
    compare_components()
        Plots component's response for selected input.
    _get_mimo_response()
        Wrapper to obtain response from a system whether it is an aircraft, SAS or AP.
    get_aircraft_response()
        Obtain aircraft's response.
    get_sas_response()
        Obtain Stability Augmentation System response.
    plot_aircraft_response()
        Plot aircraft's response for selected input type and axis.
    plot_sas_response()
        Plot Stability Augmentation System response for selected input type and axis.
    plot_full_response()
        Plot aircraft's response for selected input type for all axes.
    """
        
    AXES = (Axis.LONG, Axis.LATDIR)

    def __init__(self, lin_sys: LinearizedSystem, fig_root: str | Path = "flight_dynamics"):
        """
        Time response analyzer class.

        Creates time response plots for aircraft, components, SAS and AP.

        Attributes
        ----------
        AXES : dict[Axis, tuple(str, str)]
            Dynamic system axes.
        sys : LinearizedSystem
            Linearized system to analyze.
        fig_root : str | Path
            Figure root saving path (default is "flight_dynamics").
        labels : dict[Axis, AxisLabels]
            Axis input and output (state) channels labels.
        """

        self.sys = lin_sys
        self.fig_root = Path(fig_root)
        self.labels = AXIS_LABELS

    def _labels(self, axis: Axis):
        """
        Obtain axis labels.

        Attributes
        ----------
        axis : Axis
            Axis to get labels for.

        Returns
        ----------
        AxisLabels
            Axis labels.
        """

        return self.labels[axis]
    
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
    
    @staticmethod
    def get_component_response(
        tf: control.TransferFunction | control.StateSpace,
        t: np.ndarray, 
        input_type: InputSignal,
        amp: float = 1.0, 
        t_end: float = 10.0,
    ) -> tuple[control.StateSpace | control.StateSpace, np.ndarray, np.ndarray, np.ndarray]:
        """
        Obtain component's response for selected input.

        Attributes
        ----------
        tf : control.TransferFunction | control.StateSpace
            Transfer function to obtain response from.
        t : np.ndarray
            Time vector in seconds.
        input_type : InputSignal
            Input signal type.
        amp : float, optional
            Input signal amplitude (default is 1.0).
        t_end : float, optional
            Input signal end time in seconds (default is 10.0).

        Returns
        ----------
        tuple[control.StateSpace, np.ndarray, np.ndarray, np.ndarray]
            Component transfer function and response.   
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
    ) -> None:
        """
        Plots component's response for selected input.

        Attributes
        ----------
        tf_map : Mapping[str, control.TransferFunction]
            Transfer functions of components.
        t : np.ndarray
            Time vector in seconds.
        input_type : InputSignal
            Input signal type.
        amp : float, optional
            Input signal amplitude (default is 1.0).
        t_end : float, optional
            Input signal end time in seconds (default is 10.0).
        title : str, optional
            Figure title (default is "Actuator response Comparison").
        savefig : bool, optional
            Option to save figure (default is False).
        showfig : bool, optional
            Option to show figure (default is False).
        save_path : str | Path | None, optional
            Figure saving path (default is None).  
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

    def _get_mimo_response(
        self,
        sys: control.StateSpace,
        t: np.ndarray,
        axis: Axis,
        channel: InputChannel,
        input_type: InputSignal,
        amp: float = 1.0,
        t_end: float = 10.0,
    ) -> tuple[control.StateSpace, np.ndarray, np.ndarray, np.ndarray]:
        """
        Wrapper to obtain response from a system whether it is an aircraft, SAS or AP.

        Attributes
        ----------
        sys : control.StateSpace
            System to analyze.
        t : np.ndarray
            Time vector in seconds.
        axis : Axis
            Axis to analyze.
        channel : InputChannel
            Input channel.
        input_type : InputSignal
            Input signal type.
        amp : float, optional
            Input signal amplitude (default is 1.0).
        t_end : float, optional
            Input signal end time in seconds (default is 10.0).

        Returns
        ----------
        tuple[control.StateSpace, np.ndarray, np.ndarray, np.ndarray]
            System's response.  
        """
                
        labels = self._labels(axis)

        if channel not in labels.input_channels:
            allowed = ", ".join(ch.value for ch in labels.input_channels)
            raise ValueError(f"Invalid channel '{channel.value}' for axis '{axis.value}'. Allowed: {allowed}")

        if sys.ninputs != len(labels.input_channels):
            raise ValueError(
                f"System input count ({sys.ninputs}) does not match axis '{axis.value}' "
                f"input labels ({len(labels.input_channels)})."
            )

        u = np.zeros((sys.ninputs, len(t)))
        idx = labels.input_channels.index(channel)
        u[idx, :] = build_input(signal=input_type, t=t, amp=np.deg2rad(amp), t_end=t_end)

        t_out, y_out = control.forced_response(sys, t, u)
        y_out = np.atleast_2d(y_out)

        if axis is Axis.LONG:
            y_out[1:, :] *= 180 / np.pi
        else:
            y_out *= 180 / np.pi

        return sys, t_out, y_out, u

    def get_aircraft_response(
        self,
        t: np.ndarray,
        axis: Axis,
        channel: InputChannel,
        input_type: InputSignal,
        amp: float = 1.0,
        t_end: float = 10.0,
    ) -> tuple[control.StateSpace, np.ndarray, np.ndarray, np.ndarray]:
        """
        Obtain aircraft's response.

        Attributes
        ----------
        t : np.ndarray
            Time vector in seconds.
        axis : Axis
            Axis to analyze.
        channel : InputChannel
            Input channel.
        input_type : InputSignal
            Input signal type.
        amp : float, optional
            Input signal amplitude (default is 1.0).
        t_end : float, optional
            Input signal end time in seconds (default is 10.0).

        Returns
        ----------
        tuple[control.StateSpace, np.ndarray, np.ndarray, np.ndarray]
            Aircraft's response.  
        """
                
        sys = self.sys.get_sys(axis)
        return self._get_mimo_response(sys, t, axis, channel, input_type, amp, t_end)


    def get_sas_response(
        self,
        sas_sys: control.StateSpace,
        t: np.ndarray,
        axis: Axis,
        channel: InputChannel,
        input_type: InputSignal,
        amp: float = 1.0,
        t_end: float = 10.0,
    ) -> tuple[control.StateSpace, np.ndarray, np.ndarray, np.ndarray]:
        """
        Obtain Stability Augmentation System response.

        Attributes
        ----------
        sas_sys : control.StateSpace
            Stability Augmentation System system.
        t : np.ndarray
            Time vector in seconds.
        axis : Axis
            Axis to analyze.
        channel : InputChannel
            Input channel.
        input_type : InputSignal
            Input signal type.
        amp : float, optional
            Input signal amplitude (default is 1.0).
        t_end : float, optional
            Input signal end time in seconds (default is 10.0).

        Returns
        ----------
        tuple[control.StateSpace, np.ndarray, np.ndarray, np.ndarray]
            Stability Augmentation System response.  
        """
            
        return self._get_mimo_response(sas_sys, t, axis, channel, input_type, amp, t_end)
        
    def plot_aircraft_response(
        self, 
        t: np.ndarray, 
        axis: Axis, 
        channel: InputChannel, 
        input_type: InputSignal, 
        amp: float = 1.0, 
        t_end: float = 10.0, 
        savefig: bool = False, 
        showfig: bool = False
    ) -> None:
        """
        Plot aircraft's response for selected input type and axis.

        Attributes
        ----------
        t : np.ndarray
            Time vector in seconds.
        axis : Axis
            Axis to analyze.
        channel : InputChannel
            Input channel.
        input_type : InputSignal
            Input signal type.
        amp : float, optional
            Input signal amplitude (default is 1.0).
        t_end : float, optional
            Input signal end time in seconds (default is 10.0).
        savefig : bool, optional
            Option to save figure (default is False).
        showfig : bool, optional
            Option to show figure (default is False).
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

    def plot_sas_response(
        self,
        sas_sys: control.StateSpace,
        t: np.ndarray,
        axis: Axis,
        channel: InputChannel,
        input_type: InputSignal,
        amp: float = 1.0,
        t_end: float = 10.0,
        savefig: bool = False,
        showfig: bool = False,
    ):
        """
        Plot Stability Augmentation System response for selected input type and axis.

        Attributes
        ----------
        sas_sys : control.StateSpace
            Stability Augmentation System system.
        t : np.ndarray
            Time vector in seconds.
        axis : Axis
            Axis to analyze.
        channel : InputChannel
            Input channel.
        input_type : InputSignal
            Input signal type.
        amp : float, optional
            Input signal amplitude (default is 1.0).
        t_end : float, optional
            Input signal end time in seconds (default is 10.0).
        savefig : bool, optional
            Option to save figure (default is False).
        showfig : bool, optional
            Option to show figure (default is False).
        """

        sys, t_out, y_out, u = self.get_sas_response(
            sas_sys=sas_sys, t=t, axis=axis, channel=channel,
            input_type=input_type, amp=amp, t_end=t_end
        )

        labels = self._labels(axis)
        sys.input_labels = [f'${v}$' for v in labels.inputs]
        sys.output_labels = [f'${v}$' for v in labels.states]

        fig, axes = plt.subplots(len(sys.output_labels), 1, sharex=True)
        fig.suptitle(f'SAS response for {axis.value} axis to {amp}deg {input_type.value} in {channel.value}')

        for i, ax in enumerate(np.atleast_1d(axes)):
            y = y_out[i, :]
            ax.plot(t_out, y)
            n_tail = max(1, int(0.1 * len(y)))
            y_ss = float(np.mean(y[-n_tail:]))
            ax.axhline(y_ss, color="red", linestyle="--", linewidth=1.2, label=f"ss = {y_ss:.3g}")
            ax.set_ylabel(sys.output_labels[i])
            ax.minorticks_on()
            ax.grid(True, which="major", linestyle="-", alpha=0.5)
            ax.grid(True, which="minor", linestyle=":", alpha=0.35)
            ax.legend(loc="best")

        np.atleast_1d(axes)[-1].set_xlabel("$t$ [s]")

        if savefig:
            fig.savefig(self._figure_dir() / f"SAS_{axis.value}_{channel.value}_{input_type.value}_{amp}deg.png")
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