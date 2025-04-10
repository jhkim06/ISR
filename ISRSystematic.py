from ISRAnalyzer import ISRAnalyzer
from CMSData import CMSData
from ISRCombiner import ISRCombiner

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
                                                                             correct_binned_mean=True)
        self.default_mean_mass = self.default_analyzer.get_unfolded_mean_mass(do_acceptance_correction=True)

        self.experiments = {}
        self.mean_pt_measurement = {} # dataframe
        self.mean_mass_measurement = {}
        # combiner

    def do_experiment(self, sys_name, **kwargs):
        data = self.data
        signal = self.signal
        bg = self.bg

        # {"use_minnlo": True,
        if 'use_minnlo' in kwargs:
            data, signal, bg = self.cms_data.get_isr_samples(self.period,
                                                             self.channel,
                                                             self.lepton_selection,
                                                             use_minnlo=kwargs['use_minnlo'])

        self.experiments[sys_name] = ISRAnalyzer(data, signal, bg,
                                                 self.mass_bins, self.pt_bins,
                                                 acceptance=self.signal)

        self.mean_pt_measurement[sys_name] = (
            self.experiments[sys_name].get_unfolded_mean_pt_2d(do_acceptance_correction=True,
                                                               correct_binned_mean=True))
        self.mean_mass_measurement[sys_name] = (
            self.experiments[sys_name].get_unfolded_mean_mass(do_acceptance_correction=True))