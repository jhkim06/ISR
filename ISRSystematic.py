from ISRAnalyzer import ISRAnalyzer
from CMSData import CMSData
import numpy as np


# calculate for each period and channel about the default condition
class ISRSystematic:
    def __init__(self, period, channel, lepton_selection):

        self.period = period
        self.channel = channel
        self.lepton_selection = lepton_selection

        self.mass_bins = [
            (55.0, 64.0),
            (64.0, 81.0),
            (81.0, 101.0),
            (101.0, 200.0),
            (200.0, 1000.0),
        ]
        self.pt_bins = (0.0, 100.0)

        self.sample_base_dir = '/Users/junhokim/Work/cms_snu/data/Ultralegacy/default/'
        # Default data, analyzer
        self.cms_data = CMSData(self.sample_base_dir)
        self.data, self.signal, self.bg = self.cms_data.get_isr_samples(self.period,
                                                         self.channel,
                                                         self.lepton_selection)

        self.default_analyzer = ISRAnalyzer(self.data, self.signal, self.bg,
                                            self.mass_bins, self.pt_bins)

        self.default_mean_pt = self.default_analyzer.get_unfolded_mean_pt_2d(do_acceptance_correction=True,
                                                                             correct_binned_mean=True, as_df=True)
        self.default_mean_mass = self.default_analyzer.get_unfolded_mean_mass(do_acceptance_correction=True, as_df=True)

        # sys_name: [variation]
        self.systematic_analyzers = {}  # sys_name: [analyzer]
        self.mean_pt_systematics = {}  # sys_name: [mean pt] (dataframe)
        self.mean_mass_systematics = {}

        self.mean_pt_summary = self.default_mean_pt.copy(deep=True)  # dataframe
        self.mean_mass_summary = self.default_mean_mass.copy(deep=True)

    # TODO enable all systematic cases
    def do_experiment(self, sys_name, **kwargs):
        data = self.data
        signal = self.signal
        bg = self.bg

        if 'use_minnlo' in kwargs:
            data, signal, bg = self.cms_data.get_isr_samples(self.period,
                                                             self.channel,
                                                             self.lepton_selection,
                                                             use_minnlo=kwargs['use_minnlo'])

        if sys_name not in self.systematic_analyzers:
            self.systematic_analyzers[sys_name] = []
            self.mean_pt_systematics[sys_name] = []
            self.mean_mass_systematics[sys_name] = []

        systematic_analyzer = ISRAnalyzer(data, signal, bg,
                                          self.mass_bins, self.pt_bins,
                                          acceptance=self.signal)

        bg_scale = 1.0
        if 'bg_scale' in kwargs:
            bg_scale = kwargs['bg_scale']
            print(f'bg_scale: {bg_scale}')

        self.systematic_analyzers[sys_name].append(systematic_analyzer)
        self.mean_pt_systematics[sys_name].append(
            systematic_analyzer.get_unfolded_mean_pt_2d(do_acceptance_correction=True,
                                                        correct_binned_mean=True, as_df=True,
                                                        bg_scale=bg_scale))
        self.mean_mass_systematics[sys_name].append(
            systematic_analyzer.get_unfolded_mean_mass(do_acceptance_correction=True, as_df=True,
                                                       bg_scale=bg_scale))

    def merge_results(self):
        if 'total error' in self.mean_pt_summary.columns:
            del self.mean_pt_summary['total error']
        if 'total error' in self.mean_mass_summary.columns:
            del self.mean_mass_summary['total error']

        default_mean_pt = self.default_mean_pt["mean"]
        default_mean_mass = self.default_mean_mass["mean"]

        for sys_name, _ in self.systematic_analyzers.items():
            squared_sum_pt = np.zeros(len(self.default_mean_pt))
            squared_sum_mass = np.zeros(len(self.default_mean_pt))

            for index in range(len(self.systematic_analyzers[sys_name])):
                sys_mean_pt = self.mean_pt_systematics[sys_name][index]["mean"]
                sys_mean_mass = self.mean_mass_systematics[sys_name][index]["mean"]

                diff_squared_pt = (default_mean_pt-sys_mean_pt)**2
                diff_squared_mass = (default_mean_mass-sys_mean_mass)**2

                squared_sum_pt += diff_squared_pt
                squared_sum_mass += diff_squared_mass

            self.mean_pt_summary[f'{sys_name} error'] = np.sqrt(squared_sum_pt)
            self.mean_mass_summary[f'{sys_name} error'] = np.sqrt(squared_sum_mass)

        self.attach_total_error(r".*error$")

    def attach_total_error(self, reg_expression, new_col_name="total error"):

        selected_cols = self.mean_pt_summary.filter(regex=reg_expression)
        squared_sum = (selected_cols ** 2).sum(axis=1)
        root_of_squared_sum = np.sqrt(squared_sum)
        self.mean_pt_summary[new_col_name] = root_of_squared_sum

        selected_cols = self.mean_mass_summary.filter(regex=reg_expression)
        squared_sum = (selected_cols ** 2).sum(axis=1)
        root_of_squared_sum = np.sqrt(squared_sum)
        self.mean_mass_summary[new_col_name] = root_of_squared_sum
