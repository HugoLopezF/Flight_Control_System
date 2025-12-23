from aircraft import aircraft
from utilities.constants import GRAVITY as g 
from math import sin, cos, tan
import os
import numpy as np
from sympy import Matrix, Symbol, Number, simplify
from sympy.physics.control.lti import TransferFunction
import control
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

class SAS:
    def __init__(self, Aircraft):
        self.aircraft = Aircraft
        self.matrices = {}
        self.std_matrices = {}
        self.TF = {}
        self.output_labels = {
            'long': [r'\Delta u', r'\Delta\alpha', r'\Delta\theta', r'\Delta q'],
            'latdir': [r'\Delta\beta', r'\Delta p', r'\Delta r', r'\Delta\phi']
        }
        self.input_labels = {
            'long': [r'\Delta\delta_e'],
            'latdir': [r'\Delta\delta_a', r'\Delta\delta_r']
        }
        self.bode_lims = {
            'long': [10e-3, 10e1],
            'latdir': [10e-6, 10e1]
        }

    def round_expr(self, expr, num_digits):
        return expr.xreplace({n : round(n, num_digits) for n in expr.atoms(Number)})
    
    def plot_bode_nichols(self, savefig=False, showfig=False):
        self.get_std_matrices()
        for ax in ['long', 'latdir']:
            A = self.std_matrices[ax]['A']
            B = self.std_matrices[ax]['B']
            m = A.shape[0]
            n = B.shape[1] if len(B.shape) > 1 else 1
            C = np.eye(A.shape[0])
            D = np.zeros((m, n))
            sys = control.ss(A, B, C, D)
            sys.input_labels = self.input_labels[ax]
            sys.output_labels = self.output_labels[ax]
            for i, inp in enumerate(sys.output_labels):
                for j, out in enumerate(sys.input_labels):
                    # Figure layout
                    fig = plt.figure(figsize=(14, 8))
                    gs = gridspec.GridSpec(2, 2)

                    ax_mag = fig.add_subplot(gs[0, 0])
                    ax_phase = fig.add_subplot(gs[1, 0], sharex=ax_mag)
                    ax_nichols = fig.add_subplot(gs[:, 1])

                    # Bode plot
                    control.bode_plot(sys[i, j], omega_limits=self.bode_lims[ax], dB=True, ax=[ax_mag, ax_phase])

                    # Nichols plot
                    control.nichols_plot(sys[i, j], omega=self.bode_lims[ax], ax=ax_nichols)
                    ax_nichols.set_ylim(ax_mag.get_ylim())
                    if ax_phase.get_ylim()[1] < -170:
                        plt.close()
                        fig = plt.figure(figsize=(14, 8))
                        gs = gridspec.GridSpec(2, 2)

                        ax_mag = fig.add_subplot(gs[0, 0])
                        ax_phase = fig.add_subplot(gs[1, 0], sharex=ax_mag)
                        ax_nichols = fig.add_subplot(gs[:, 1])
                        mag, phase, omega = control.frequency_response(sys[i, j], omega_limits=self.bode_lims[ax])

                        ax_mag.semilogx(omega, 20 * np.log10(mag))
                        ax_mag.set_ylabel("Magnitude [dB]")
                        ax_mag.grid(True, which="both")
                        ax_mag.set_xlim(self.bode_lims[ax])

                        ax_phase.semilogx(omega, np.unwrap(np.degrees(phase), discont=8 * np.pi))
                        ax_phase.set_ylabel("Phase [deg]")
                        ax_phase.set_xlabel("Frequency [rad/s]")
                        ax_phase.grid(True, which="both")

                        control.nichols_plot(sys[i, j], omega=self.bode_lims[ax], ax=ax_nichols)
                        ax_nichols.set_ylim(ax_mag.get_ylim())
                        min_phase = ax_phase.get_ylim()[0]
                        max_phase = ax_phase.get_ylim()[1]
                        ax_nichols.set_xlim((min_phase, max_phase))
                    else:
                        ax_nichols.set_xlim(ax_phase.get_ylim())
                    
                    ax_mag.set_title("Bode")
                    ax_nichols.set_title("Nichols")
                    plt.tight_layout()
                    if savefig:
                        inp_var = inp.replace('\\Delta', '').replace('\\', '').replace(' ', '')
                        out_var = out.replace('\\Delta', '').replace('\\', '')
                        figname = self.aircraft.model + f'_{inp_var}_{out_var}_Bode_Nichols.png'
                        figdir = os.path.join(os.getcwd(), 'aircraft', self.aircraft.model)
                        os.makedirs(figdir, exist_ok=True)
                        plt.savefig(os.path.join(figdir, figname))
                    if showfig:
                        plt.show()
                    plt.close()
    
    def plot_bode(self, savefig=False, showfig=False):
        self.get_std_matrices()
        for ax in ['long', 'latdir']:
            A = self.std_matrices[ax]['A']
            B = self.std_matrices[ax]['B']
            m = A.shape[0]
            n = B.shape[1] if len(B.shape) > 1 else 1
            C = np.eye(A.shape[0])
            D = np.zeros((m, n))
            sys = control.ss(A, B, C, D)
            sys.input_labels = self.input_labels[ax]
            sys.output_labels = self.output_labels[ax]
            for i, inp in enumerate(sys.output_labels):
                for j, out in enumerate(sys.input_labels):
                    plt.figure()
                    control.bode_plot(sys[i, j], omega_limits=self.bode_lims[ax], dB=True)
                    plt.suptitle(fr'Bode diagram of $[G(i\omega)]_{{{inp+out}}} = \frac{{{inp}}}{{{out}}}$')
                    plt.tight_layout()
                    if savefig:
                        inp_var = inp.replace('\\Delta', '').replace('\\', '').replace(' ', '')
                        out_var = out.replace('\\Delta', '').replace('\\', '')
                        figname = self.aircraft.model + f'_{inp_var}_{out_var}_Bode.png'
                        figdir = os.path.join(os.getcwd(), 'aircraft', self.aircraft.model)
                        os.makedirs(figdir, exist_ok=True)
                        plt.savefig(os.path.join(figdir, figname))
                    if showfig:
                        plt.show()
                    plt.close()
    
    def plot_nichols(self, savefig=False, showfig=False):
        self.get_std_matrices()
        for ax in ['long', 'latdir']:
            A = self.std_matrices[ax]['A']
            B = self.std_matrices[ax]['B']
            m = A.shape[0]
            n = B.shape[1] if len(B.shape) > 1 else 1
            C = np.eye(A.shape[0])
            D = np.zeros((m, n))
            sys = control.ss(A, B, C, D)
            sys.input_labels = self.input_labels[ax]
            sys.output_labels = self.output_labels[ax]
            for i, inp in enumerate(sys.output_labels):
                for j, out in enumerate(sys.input_labels):
                    plt.figure()
                    control.nichols_plot(sys[i, j], omega=self.bode_lims[ax])
                    plt.suptitle(fr'Nichols diagram of $[G(i\omega)]_{{{inp+out}}} = \frac{{{inp}}}{{{out}}}$')
                    plt.tight_layout()
                    if savefig:
                        inp_var = inp.replace('\\Delta', '').replace('\\', '').replace(' ', '')
                        out_var = out.replace('\\Delta', '').replace('\\', '')
                        figname = self.aircraft.model + f'_{inp_var}_{out_var}_Nichols.png'
                        figdir = os.path.join(os.getcwd(), 'aircraft', self.aircraft.model)
                        os.makedirs(figdir, exist_ok=True)
                        plt.savefig(os.path.join(figdir, figname))
                    if showfig:
                        plt.show()
                    plt.close()

    def print_TF(self):
        self.get_TF()
        for ax in ['long', 'latdir']:
            TF = self.round_expr(self.TF[ax], 3)
            print(f'Printing {ax} transfer functions...\n')
            for row in range(0, self.TF[ax].shape[0]):
                print('\t' + str(TF[row]) + '\n')

    def get_TF(self):
        self.get_std_matrices()
        s = Symbol('s')
        for ax in ['long', 'latdir']:
            A = self.std_matrices[ax]['A']
            B = self.std_matrices[ax]['B']
            I = np.eye(A.shape[0])
            self.TF[ax] = Matrix(Matrix(s * I - A).inv().applyfunc(simplify) @ B).applyfunc(simplify)
 
    def get_std_matrices(self):
        self.get_all_matrices()
        for ax in ['long', 'latdir']:
            self.std_matrices[ax] = {
                'A': np.linalg.solve(self.matrices[ax]['E'], self.matrices[ax]['A_prime']),
                'B': np.linalg.solve(self.matrices[ax]['E'], self.matrices[ax]['B_prime'])
            }

    def get_all_matrices(self):
        E, A, B = self.get_long_matrices()
        self.matrices['long'] = {
            'E': E, 
            'A_prime': A,
            'B_prime': B
        }

        E, A, B = self.get_latdir_matrices()
        self.matrices['latdir'] = {
            'E': E, 
            'A_prime': A, 
            'B_prime': B
        }

    def get_long_matrices(self):
        W = self.aircraft.mass_prop['W']
        u_s = self.aircraft.stab_der.FlightCondition['u_s']
        theta_s = self.aircraft.stab_der.FlightCondition['theta_s']
        Zw_dot = self.aircraft.stab_der.Zw_dot
        Zq = self.aircraft.stab_der.Zq
        Mw_dot = self.aircraft.stab_der.Mw_dot
        Mq = self.aircraft.stab_der.Mq
        I_yy = self.aircraft.mass_prop['I_yy']
        Xu = self.aircraft.stab_der.Xu
        Xw = self.aircraft.stab_der.Xw
        Zu = self.aircraft.stab_der.Zu
        Zw = self.aircraft.stab_der.Zw
        Mu = self.aircraft.stab_der.Mu
        Mw = self.aircraft.stab_der.Mw
        Xdelta_e = self.aircraft.stab_der.Xdelta_e
        Zdelta_e = self.aircraft.stab_der.Zdelta_e
        Mdelta_e = self.aircraft.stab_der.Mdelta_e

        E = np.array([     
                [W, 0, 0, 0],
                [0, u_s * (W - Zw_dot), 0, 0],
                [0, -u_s * Mw_dot, 0, I_yy],
                [0, 0, 1, 0]
        ])
        A_prime = np.array([
            [Xu, u_s * Xw, -W * g * cos(theta_s), 0],
            [Zu, u_s * Zw, -W * g * sin(theta_s), W * u_s + Zq],
            [Mu, u_s * Mw, 0, Mq],
            [0, 0, 0, 1]
        ])
        B_prime = np.array([Xdelta_e, Zdelta_e, Mdelta_e, 0])

        return E, A_prime, B_prime
    
    def get_latdir_matrices(self):
        W = self.aircraft.mass_prop['W']
        u_s = self.aircraft.stab_der.FlightCondition['u_s']
        theta_s = self.aircraft.stab_der.FlightCondition['theta_s']
        I_xx = self.aircraft.mass_prop['I_xx']
        I_zz = self.aircraft.mass_prop['I_zz']
        I_xz = self.aircraft.mass_prop['I_xz']
        Yv = self.aircraft.stab_der.Yv
        Yp = self.aircraft.stab_der.Yp
        Yr = self.aircraft.stab_der.Yr
        Lv = self.aircraft.stab_der.Lv
        Lp = self.aircraft.stab_der.Lp
        Lr = self.aircraft.stab_der.Lr
        Nv = self.aircraft.stab_der.Nv
        Np = self.aircraft.stab_der.Np
        Nr = self.aircraft.stab_der.Nr
        Ydelta_r = self.aircraft.stab_der.Ydelta_r
        Ldelta_a = self.aircraft.stab_der.Ldelta_a
        Ldelta_r = self.aircraft.stab_der.Ldelta_r
        Ndelta_a = self.aircraft.stab_der.Ndelta_a
        Ndelta_r = self.aircraft.stab_der.Ndelta_r

        E = np.array([
            [u_s * W, 0, 0, 0],
            [0, I_xx, -I_xz, 0],
            [0, -I_xz, I_zz, 0],
            [0, 0, 0, 1]
        ])
        A_prime = np.array([
            [u_s * Yv, Yp, Yr - W * u_s, W * g * cos(theta_s)],
            [u_s * Lv, Lp, Lr, 0],
            [u_s * Nv, Np, Nr, 0],
            [0, 1, tan(theta_s), 0]
        ])
        B_prime = np.array([
            [0, Ydelta_r],
            [Ldelta_a, Ldelta_r],
            [Ndelta_a, Ndelta_r],
            [0, 0]
        ])
        
        return E, A_prime, B_prime