from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from .observable import observables, NumericalSimulation


class Visualization(NumericalSimulation):
    """
    Plotting parameters for customizing figure properties

    Attributes
    ----------
    cm : matplotlib.colors.ListedColormap (default: plt.cm.get_cmap('tab20'))
        Choosing colormaps for 'cmap'.
        If num of conditions is less than 11, use dark color only ('tab20').

    timecourse_options : list of dict
        Plotting options for figure/simulation/<viz_type>/<each_observable>.

            Keys
            ----
            * 'divided_by' : int or float (default: 1)
                Convert time unit. (e.g. sec -> min).

            * 'xlim' : tuple
                Set the x limits of the current axes.

            * 'xticks' : list (default: None)
                Set the current tick locations of the x-axis.

            * 'xlabel' : str (default: 'Time')
                Set the label for the x-axis.

            * 'ylim' : tuple
                Set the y limits of the current axes.

            * 'yticks' : list (default: None)
                Set the current tick locations of the y-axis.

            * 'ylabel' : str (default: observables[i].replace('__', '\n').replace('_', ' '))
                Set the label for the y-axis.

            * 'exp_data' : bool (default: True)
                if False, experimental data will not be shown.

            * 'legend_loc' : Location String (default: None)
                Set the location of the legend. If 'legend_loc' is None,
                pyplot.legend will not be shown.

            * 'cmap' : list or tuple
                Set colormap.

            * 'shape' : list or tuple of strings (default: Line2D.filled_markers)
                Set markers.

            * 'dont_show' : list of strings
                Set conditions you don't want to plot.

    multiplot_options : dict
        Plotting options for figure/simulation/<viz_type>/multiplot_observables.

    sensitivity_options : dict
        Plotting options for figure/sensitivity.

    """
    def __init__(self):
        self.cm = plt.cm.get_cmap('tab20')

        self.timecourse_options = [
            {
                'divided_by' : 1,
                'xlim' : (),
                'xticks' : None,
                'xlabel': 'Time',
                'ylim' : (),
                'yticks' : None,
                'ylabel': observables[i].replace('__', '\n').replace('_', ' '),
                'exp_data' : True,
                'legend_loc' : None,
                'cmap' : [self.cm.colors[j] for j in (range(20)
                            if len(self.conditions) > 10
                            else range(0, 20, 2))],
                'shape' : Line2D.filled_markers,
                'dont_show' : [],
            } for i, _ in enumerate(observables)]

        self.multiplot_options = {
            'fig_name' : 'multiplot_observables',
            'observables' : [],
            'condition' : None,
            'xlim' : (),
            'xticks' : None,
            'xlabel': 'Time',
            'ylim' : (),
            'yticks' : None,
            'ylabel': '',
            'cmap' : [self.cm.colors[j] for j in (range(20)
                        if len(self.conditions) > 10
                        else range(0, 20, 2))],
            'shape' : Line2D.filled_markers,
        }

        self.sensitivity_options = {
            'figsize' : (12, 5),
            'width' : 0.3,
            'legend_loc' : 'upper left',
            'cmap' : ['darkslateblue', 'orangered'],
        }

    def get_timecourse_options(self):
        for i, _ in enumerate(observables):
            self.timecourse_options[i]['divided_by'] = 60  # sec. -> min.
            self.timecourse_options[i]['xlim'] = (-5, 125)
            self.timecourse_options[i]['xticks'] = [0, 30, 60, 90, 120]
            self.timecourse_options[i]['xlabel'] = 'Time (min)'
            self.timecourse_options[i]['ylim'] = (-0.1, 1.3)
            self.timecourse_options[i]['yticks'] = [0.0, 0.3, 0.6, 0.9, 1.2]
            self.timecourse_options[i]['cmap'] = ['darkslateblue', 'orangered']
            self.timecourse_options[i]['shape'] = ['D', 's']
            self.timecourse_options[i]['dont_show'] = []

        self.timecourse_options[
            observables.index('Phosphorylated_cFos')
        ]['ylabel'] = 'Phosphorylated c-Fos'

        return self.timecourse_options

    def multiplot_observables(self):
        return self.multiplot_options

    @staticmethod
    def set_timecourse_rcParams():
        """ figure/simulation
        """
        plt.rcParams['font.size'] = 20
        plt.rcParams['axes.linewidth'] = 1.5
        plt.rcParams['xtick.major.width'] = 1.5
        plt.rcParams['ytick.major.width'] = 1.5
        plt.rcParams['lines.linewidth'] = 1.8
        plt.rcParams['lines.markersize'] = 12
        plt.rcParams['font.family'] = 'Arial'
        # plt.rcParams['mathtext.fontset'] = 'custom'
        # plt.rcParams['mathtext.it'] = 'Arial:italic'

    @staticmethod
    def set_param_range_rcParams():
        """ figure/param_range
        """
        plt.rcParams['font.size'] = 12
        plt.rcParams['axes.linewidth'] = 1.2
        plt.rcParams['xtick.major.width'] = 1.2
        plt.rcParams['ytick.major.width'] = 1.2
        plt.rcParams['font.family'] = 'Arial'

    @staticmethod
    def set_sensitivity_rcParams():
        """ figure/sensitivity
        """
        plt.rcParams['font.size'] = 12
        plt.rcParams['axes.linewidth'] = 1.2
        plt.rcParams['xtick.major.width'] = 1.2
        plt.rcParams['ytick.major.width'] = 1.2
        plt.rcParams['font.family'] = 'Arial'

    @staticmethod
    def convert_species_name(name):
        """ figure/sensitivity/initial_condition
        - Sensitivity for species with nonzero initial conditions
        """
        if name == 'E1':
            return 'ErbB1'
        elif name == 'E2':
            return 'ErbB2'
        elif name == 'E3':
            return 'ErbB3'
        elif name == 'E4':
            return 'ErbB4'
        elif name == 'G':
            return 'Grb2'
        elif name == 'S':
            return 'Shc'
        elif name == 'A':
            return 'Gab1'
        elif name == 'I':
            return 'PI3K'
        elif name == 'O':
            return 'SOS'
        elif name == 'P2':
            return 'PIP2'
        elif name == 'R':
            return 'RasGAP'
        elif name == 'RsD':
            return 'RasGDP'
        elif name == 'T':
            return 'PTP1B'
        elif name == 'ERKc':
            return 'ERK'
        elif name == 'RSKc':
            return 'RSK'
        elif name == 'CREBn':
            return 'CREB'
        elif name == 'Elk1n':
            return 'Elk1'
        else:
            return name
