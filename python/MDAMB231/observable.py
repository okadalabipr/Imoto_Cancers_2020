import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from scipy.integrate import ode

from .name2idx import C, V
from .set_model import DifferentialEquation
from load_csv import mul_ratio2mcf7


observables = [
    'Phosphorylated_Akt',
    'Phosphorylated_ERK',
    'Phosphorylated_cFos',
]


class NumericalSimulation(DifferentialEquation):
    """ Simulate a model using scipy.integrate.ode

    Attributes
    ----------
    normalization : bool
        if True, simulation results in each observable are divided by their
        maximum values.
    """
    def __init__(self):
        super().__init__(perturbation={})
        self.normalization = True

    t = range(7201)  # 0, 1, 2, ..., 7200 (Unit: sec.)

    # Experimental conditions
    conditions = ['EGF', 'HRG']

    simulations = np.empty((len(observables), len(t), len(conditions)))

    def simulate(self, x, y0, _perturbation={}):
        if _perturbation:
            self.perturbation = _perturbation
        # (x, y0) = mul_ratio2mcf7(x, y0, 'MDA-MB-231')
        # get steady state
        y0[V.E] = 0.0
        y0[V.H] = 0.0
        y0 = self._get_steady_state(self.diffeq, y0, tuple(x))
        if not y0:
            return False
        # add ligand
        for i, condition in enumerate(self.conditions):
            if condition == 'EGF':
                y0[V.E] = 10.0
                y0[V.H] = 0.0
            elif condition == 'HRG':
                y0[V.E] = 0.0
                y0[V.H] = 10.0

            (T, Y) = self._solveode(self.diffeq, y0, self.t, tuple(x))

            if T[-1] < self.t[-1]:
                return False
            else:
                self.simulations[observables.index('Phosphorylated_Akt'), :, i] = (
                    Y[:, V.Aktstar]
                )
                self.simulations[observables.index('Phosphorylated_ERK'), :, i] = (
                    (Y[:, V.pERKn] + Y[:, V.ppERKn]) * (x[C.Vn]/x[C.Vc])
                    + Y[:, V.pERKc] + Y[:, V.ppERKc]
                )
                self.simulations[observables.index('Phosphorylated_cFos'), :, i] = (
                    Y[:, V.pcFOSn]*(x[C.Vn]/x[C.Vc]) + Y[:, V.pcFOSc]
                )

    def _solveode(self, diffeq, y0, tspan, args):
        """
        Solve a system of ordinary differential equations.

        Parameters
        ----------
        diffeq : callable f(y, t, f_args)
            Right-hand side of the differential equation.
        y0 : array
            Initial condition on y (can be a vector).

        tspan : array
            A sequence of time points for which to solve for y.

        args : tuple
            Model parameters.

        Returns
        -------
        T, Y : tuple
            T : array, shape (len(t))
                Evaluation points.
            Y : array, shape (len(t), len(y0))
                Array containing the value of y for each desired time in t,
                with the initial value y0 in the first row.
        """
        dt = (self.t[-1] - self.t[0]) / (len(self.t) - 1)
        sol = ode(diffeq)
        sol.set_integrator(
            'vode', method='bdf', with_jacobian=True,
            atol=1e-9, rtol=1e-9, min_step=1e-8
        )
        sol.set_initial_value(y0, tspan[0])
        sol.set_f_params(args)

        T = [tspan[0]]
        Y = [y0]

        while sol.successful() and sol.t < tspan[-1]:
            sol.integrate(sol.t+dt)
            T.append(sol.t)
            Y.append(sol.y)

        return np.array(T), np.array(Y)

    def _get_steady_state(self, diffeq, y0, args, eps=1e-6):
        """
        Find the steady state for the untreated condition.

        Parameters
        ----------
        diffeq : callable f(y, t, f_args)
            Right-hand side of the differential equation.
        y0 : array
            Initial condition on y (can be a vector).

        args : tuple
            Model parameters.

        eps : float (default: 1e-6)
            Run until a time t for which the maximal absolutevalue of the
            regularized relative derivative was smaller than eps.

        Returns
        -------
        y0 : array
            Steady state concentrations of all species.
        """
        while True:
            (T, Y) = self._solveode(diffeq, y0, range(2), args)
            if T[-1] < 1 or \
                    np.max(
                        np.abs((Y[-1, :] - y0) / (np.array(y0) + eps))
                    ) < eps:
                break
            else:
                y0 = Y[-1, :].tolist()

        return [] if T[-1] < 1 else Y[-1, :].tolist()

class ExperimentalData(object):
    """
    Set experimental data.

    Attributes
    ----------
    experiments : list of dict
        Time series data.

    error_bars : list of dict
        Error bars to show in figures.

    """
    def __init__(self):
        self.experiments = [None] * len(observables)
        self.error_bars = [None] * len(observables)

    @staticmethod
    def _norm01(egf1, hrg1, egf2, hrg2, egf3, hrg3):
        """Normalize data from 0 to 1 and Get standard error of the mean

        Parameters
        ----------
        egf%d : list
            EGF stimulation #%d
        hrg%d : list
            HRG stimulation #%d

        Returns
        -------
        egf_ave, hrg_ave, egf_sem, hrg_sem : tuple
            Averaged values and their standard error of the mean
        """
        data1 = np.array([egf1, hrg1])
        data2 = np.array([egf2, hrg2])
        data3 = np.array([egf3, hrg3])

        # max -> 1 to compare biological replicates
        data1 = data1 / np.max(data1)
        data2 = data2 / np.max(data2)
        data3 = data3 / np.max(data3)

        egf_ave = np.mean(np.stack([data1[0], data2[0], data3[0]]), axis=0)
        hrg_ave = np.mean(np.stack([data1[1], data2[1], data3[1]]), axis=0)

        ave_vec = np.stack([egf_ave, hrg_ave])
        ave_min = np.min(ave_vec)
        ave_max = np.max(ave_vec)

        # To normalize max -> 1, min -> 0
        data1 = (data1 - ave_min) / (ave_max - ave_min)
        data2 = (data2 - ave_min) / (ave_max - ave_min)
        data3 = (data3 - ave_min) / (ave_max - ave_min)

        egf_ave = np.mean(np.stack([data1[0], data2[0], data3[0]]), axis=0)
        hrg_ave = np.mean(np.stack([data1[1], data2[1], data3[1]]), axis=0)
        egf_sem = np.std(np.stack([data1[0], data2[0], data3[0]]), axis=0, ddof=1)/(3**0.5)
        hrg_sem = np.std(np.stack([data1[1], data2[1], data3[1]]), axis=0, ddof=1)/(3**0.5)

        return egf_ave, hrg_ave, egf_sem, hrg_sem

    def set_data(self):
        data_pAkt = self._norm01(
            egf1=[86491, 826975, 575400, 354446, 143728, 107326, 88082, 108892],
            hrg1=[86491, 85990, 114224, 94636, 66139, 92359, 89540, 126099],
            egf2=[56004, 816102, 296684, 165282, 92711, 83090, 55002, 50605],
            hrg2=[56004, 50640, 64646, 55814, 70627, 66657, 68994, 46922],
            egf3=[36061, 729961, 218315, 166109, 68746, 53331, 45191, 21852],
            hrg3=[36061, 75029, 174050, 89546, 88262, 99161, 78673, 58412],
        )
        self.experiments[observables.index('Phosphorylated_Akt')] = {
            'EGF': data_pAkt[0],
            'HRG': data_pAkt[1],
        }
        self.error_bars[observables.index('Phosphorylated_Akt')] = {
            'EGF': data_pAkt[2],
            'HRG': data_pAkt[3],
        }

        data_pERK = self._norm01(
            egf1=[314472, 504819, 607786, 618492, 475195, 376035, 293988, 324600],
            hrg1=[314472, 156705, 183456, 277862, 141450, 199719, 253331, 407923],
            egf2=[458693, 1001334, 875594, 834259, 782639, 815888, 629576, 539187],
            hrg2=[458693, 322542, 403985, 331734, 263578, 262142, 276371, 313541],
            egf3=[365691, 747932, 937413, 945635, 870059, 706306, 510590, 451927],
            hrg3=[365691, 340876, 428510, 303543, 269653, 195660, 215972, 388050],
        )
        self.experiments[observables.index('Phosphorylated_ERK')] = {
            'EGF': data_pERK[0],
            'HRG': data_pERK[1],
        }
        self.error_bars[observables.index('Phosphorylated_ERK')] = {
            'EGF': data_pERK[2],
            'HRG': data_pERK[3],
        }

        data_pcFos = self._norm01(
            egf1=[124104, 175391, 230831, 210449, 337825, 317246, 375420, 311936],
            hrg1=[124104, 124400, 102863, 125643, 124228, 144544, 156551, 163124],
            egf2=[60197, 45929, 63355, 77438, 115070, 122509, 166049, 129227],
            hrg2=[60197, 56254, 57765, 67667, 60445, 66214, 33948, 32766],
            egf3=[47449, 63064, 67735, 183991, 275103, 183211, 552227, 362380],
            hrg3=[47449, 88082, 93539, 73279, 44637, 55030, 56284, 35506],
        )
        self.experiments[observables.index('Phosphorylated_cFos')] = {
            'EGF': data_pcFos[0],
            'HRG': data_pcFos[1],
        }
        self.error_bars[observables.index('Phosphorylated_cFos')] = {
            'EGF': data_pcFos[2],
            'HRG': data_pcFos[3],
        }

    @staticmethod
    def get_timepoint(obs_name):
        if obs_name in observables:
            return [i*60 for i in [0, 5, 15, 30, 45, 60, 90, 120]]