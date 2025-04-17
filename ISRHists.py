from Hist import Hist
import copy
from ISR2DHist import ISR2DHist

'''
class ISRMatrix(Hist):
    def __init__(self, hist,  # TH2
                 tunfold_folded_bin,
                 tunfold_unfolded_bin,
                 hist_name='',
                 experiment='CMS',
                 label='',  # to be used in legend of histogram
                 channel='',
                 year='',
                 is_measurement=True,
                 is_mc_signal=False,):

        super(ISR2DHist, self).__init__(hist,
                                        hist_name=hist_name,
                                        experiment=experiment,
                                        label=label,
                                        channel=channel,
                                        year=year,
                                        is_measurement=is_measurement,
                                        is_mc_signal=is_mc_signal,)
'''

# group of hists for ISR
class ISRHistSet:
    def __init__(self, measurement_hist, signal_hist, signal_fake_hist, background_hists):
        # detector hist
        self.measurement_hist = measurement_hist
        self.signal_hist = signal_hist
        self.signal_fake_hist = signal_fake_hist
        self.background_hists = background_hists  # Note this is dictionary of hists

        self.response_matrix = None

        self.unfolded_measurement_hist = None
        self.acceptance_corrected_measurement_hist = None


class ISRHists:
    def __init__(self,
                 mass_bins,  # mass windows
                 pt_bins,
                 is_2d, is_pt,
                 folded_tunfold_bin=None, unfolded_tunfold_bin=None):  # pt cut

        self.mass_bins = mass_bins
        self.pt_bins = pt_bins
        self.is_2d = is_2d
        self.is_pt = is_pt

        # detector level hists
        self.isr_hists = []

        self.folded_tunfold_bin = folded_tunfold_bin
        self.unfolded_tunfold_bin = unfolded_tunfold_bin

        #
        self.mean_values = None

    # Note set *detector* level hists
    def set_isr_hists(self, measurement_hist, signal_hist, signal_fake_hist, background_hists):

        if self.is_pt and len(self.mass_bins) == len(self.isr_hists):
            print('Check the number of mass bins and the number of ISR hists...')
            return
        else:
            if self.is_2d:
                new_background_hists = {}
                for key, value in background_hists.items():
                    new_background_hists[key] = ISR2DHist(value, self.folded_tunfold_bin,)

                self.isr_hists.append(
                    ISRHistSet(
                        ISR2DHist(measurement_hist, self.folded_tunfold_bin),
                        ISR2DHist(signal_hist, self.folded_tunfold_bin),
                        ISR2DHist(signal_fake_hist, self.folded_tunfold_bin),
                        new_background_hists
                    )
                )
            else:
                self.isr_hists.append(
                    ISRHistSet(measurement_hist, signal_hist, signal_fake_hist, background_hists)
                )

    def get(self, hist_type='signal', mass_window_index=-1):
        # get hist object
        if self.is_pt:
            if self.is_2d:
                isr_hist_to_draw = self.isr_hists[0]
            else:
                isr_hist_to_draw = self.isr_hists[mass_window_index]
        else:
            # for mass, only one isr_hist will be used
            isr_hist_to_draw = self.isr_hists[0]

        attr_name = hist_type + "_hist"
        target_hist = getattr(isr_hist_to_draw, attr_name)
        if target_hist is None:
            raise ValueError(f"No attribute names '{attr_name}' found")

        # Note for backgrounds are dictionary of hists
        if hist_type == 'background':
            pass
        else:
            if self.is_2d:
                if -1 < mass_window_index < len(self.mass_bins):
                    return target_hist.extract_1d_hist(mass_window_index)
                else:
                    return target_hist
            else:
                return target_hist

    def draw_detector_level(self, mass_window_index=-1):

        measurement_hist = self.get(hist_type='measurement', mass_window_index=mass_window_index)
        signal_hist = self.get(hist_type='signal', mass_window_index=mass_window_index)

        plotter = measurement_hist.plotter

        plotter.init_plotter()
        plotter.add_hist(measurement_hist, **{'histtype':'errorbar'})
        plotter.add_hist(signal_hist, as_stack=True, as_denominator=True)
        plotter.draw_hist()

        plotter.save_and_reset_plotter(measurement_hist.hist_name)

    def draw_unfolded_level(self):
        pass

    def draw_acceptance_corrected_level(self):
        pass
