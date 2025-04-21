from Hist import Hist, change_to_greek
from Analyzer import labels, colors, get_hist_kwargs
from ISR2DHist import ISR2DHist
import pandas as pd
import numpy as np


# group of hists for ISR
class ISRHistSet:
    def __init__(self, measurement_hist, signal_hist, signal_fake_hist, background_hists,
                 response_matrix):
        # detector hist
        self.measurement_hist = measurement_hist
        self.signal_hist = signal_hist
        self.signal_fake_hist = signal_fake_hist
        self.background_hist = background_hists  # Note this is dictionary of hists

        self.response_matrix = response_matrix

        self.tunfolder = None
        self.unfolded_measurement_hist = None
        self.unfolded_signal_hist = None  # request to tunfolder

        self.acceptance_corrected_measurement_hist = None
        self.acceptance_corrected_signal_hist = None  # request to acceptance_corrector


# ISRPtHists ISRMassHists
class ISRHists:
    def __init__(self,
                 mass_bins,  # mass windows
                 pt_bins,
                 is_2d, is_pt,
                 year='', channel='',
                 folded_tunfold_bin=None, unfolded_tunfold_bin=None):  # pt cut

        self.year = year
        self.channel = channel

        self.mass_bins = mass_bins
        self.pt_bins = pt_bins

        self.is_2d = is_2d
        self.is_pt = is_pt

        # detector level hists
        self.isr_hists = []

        self.folded_tunfold_bin = folded_tunfold_bin
        self.unfolded_tunfold_bin = unfolded_tunfold_bin

        #
        self.measurement_mean_values = None
        self.measurement_bg_subtracted_mean_values = None
        self.signal_mean_values = None
        self.background_mean_values = None

        self.unfolded_measurement_mean_values = None
        self.unfolded_signal_mean_values = None

        self.acceptance_corrected_measurement_mean_values = []
        self.acceptance_corrected_signal_mean_values = []

        self.binned_mean_correction_factors = None

        # common plot cosmetics
        if is_pt:
            self.x_axis_label = r"$p_{T}^{"+change_to_greek(self.channel)+"}$"
        else:
            self.x_axis_label = r"$m^{"+change_to_greek(self.channel)+"}$"

    def update_mean_values(self, other, sys_name):
        for index, mean_value in enumerate(self.acceptance_corrected_measurement_mean_values):
            sys_mean = other.acceptance_corrected_measurement_mean_values[index]['mean']
            mean = mean_value['mean']
            insert_idx = mean_value.columns.get_loc('total_error')
            mean_value.insert(insert_idx, sys_name, abs(sys_mean-mean))

            error_columns = mean_value.columns.difference(['mean'])
            mean_value['total_error'] = np.sqrt((mean_value[error_columns] ** 2).sum(axis=1))

    # Note set *detector* level hists
    def set_isr_hists(self, measurement_hist, signal_hist, signal_fake_hist, background_hists, matrix):

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
                        new_background_hists,
                        matrix
                    )
                )
            else:
                self.isr_hists.append(
                    ISRHistSet(measurement_hist, signal_hist, signal_fake_hist, background_hists, matrix)
                )

    def set_acceptance_corrected_hist(self, measurement_hist, signal_hist, mass_window_index=0,):

        self.isr_hists[mass_window_index].acceptance_corrected_measurement_hist = measurement_hist
        self.isr_hists[mass_window_index].acceptance_corrected_signal_hist = signal_hist

        # set mean values
        if self.is_2d or self.is_pt==False:
            # loop over mass bins
            if self.is_pt:
                for index in range(len(self.mass_bins)):
                    extracted_hist = measurement_hist.extract_1d_hist(index=index)
                    mean = extracted_hist.get_mean_df()
                    self.acceptance_corrected_measurement_mean_values.append(mean)
            else:
                for low, high in self.mass_bins:
                    mean = measurement_hist.get_mean_df(range_min=low, range_max=high)
                    self.acceptance_corrected_measurement_mean_values.append(mean)
        else:
            # set mean for mass_window_index
            # TODO 1D case self.is_pt == True, self.is_2d == False
            if self.is_pt:
                mean = measurement_hist.get_mean_df()
                self.acceptance_corrected_measurement_mean_values.append(mean)

    def get_isr_hists(self, mass_window_index=-1):
        if self.is_2d or self.is_pt==False:
            isr_hist = self.isr_hists[0]
        else:
            isr_hist = self.isr_hists[mass_window_index]
        return (isr_hist.measurement_hist,
                isr_hist.signal_fake_hist,
                isr_hist.background_hist,
                isr_hist.response_matrix)

    def set_isr_response_matrices(self, response_matrix, mass_window_index=-1):
        if self.is_pt:
            if self.is_2d:
                self.isr_hists[0].response_matrix = response_matrix
            else:
                pass
        else:
            self.isr_hists[0].response_matrix = response_matrix

    def get(self, hist_type='signal', mass_window_index=-1, bin_width_norm=False):
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

        if self.is_2d:
            # decide whether to extract 1d hist or not
            if -1 < mass_window_index < len(self.mass_bins):
                if hist_type == 'background':
                    extracted_background_hists = {}
                    for bg_name, bg_hist in target_hist.items():
                        extracted_background_hists[bg_name] = bg_hist.extract_1d_hist(mass_window_index)
                    return extracted_background_hists
                else:
                    return target_hist.extract_1d_hist(mass_window_index)  # bin_width_norm applied in extract_1d_hist
            else:
                return target_hist
        else:
            return target_hist
        # hist.bin_width_norm()

    def get_additional_text_on_plot(self, mass_window_index=-1):
        text = ''
        if self.is_pt:
            if self.is_2d:
                text = (str(int(self.mass_bins[mass_window_index][0]))+
                        "$<m^{"+change_to_greek(self.channel)+"}<$"+
                        str(int(self.mass_bins[mass_window_index][1])) + " GeV")
        else:
            text = str(r"$p_{T}^{" + change_to_greek(self.channel) + "}<$" + str(int(self.pt_bins[1])) + " (GeV)")
        return text

    def draw_isr_plot(self, other, save_and_reset_plotter=True):
        if self.is_pt and other.is_pt == False:
            isr_pt = self
            isr_mass = other
        elif self.is_pt == False and other.is_pt:
            isr_pt = other
            isr_mass = self
        else:
            raise ValueError("Cannot draw ISR plot between two ISR histograms with different dimensionality")

        measurement_hist = self.get(hist_type='acceptance_corrected_measurement', mass_window_index=0)
        plotter = measurement_hist.plotter

        plotter.init_plotter(figsize=(10,8), rows=1, cols=1)
        plotter.set_experiment_label(**{'year': measurement_hist.year})

        mass=pd.concat(isr_mass.acceptance_corrected_measurement_mean_values, ignore_index=True)
        # pt=pd.concat(isr_pt.acceptance_corrected_measurement_mean_values, ignore_index=True)
        pt_scaled=pd.concat(isr_pt.acceptance_corrected_measurement_mean_values, ignore_index=True)
        pt_scaled['mean'] = pt_scaled['mean'] * self.binned_mean_correction_factors

        # plotter.add_errorbar((mass, pt), color='black', marker='.')
        plotter.add_errorbar((mass, pt_scaled), color='black', marker='.')
        plotter.draw_errorbar()

        plotter.set_isr_plot_cosmetics(channel=change_to_greek(self.channel),)
        text = isr_mass.get_additional_text_on_plot()
        plotter.add_text(text=text, location=(0, 0), do_magic=False, **{"frameon": False, "loc": "upper left", })

        if save_and_reset_plotter:
            plotter.save_and_reset_plotter("isr_test"+"_"+self.channel+self.year)

    def draw_detector_level(self, mass_window_index=-1):

        measurement_hist = self.get(hist_type='measurement', mass_window_index=mass_window_index)
        signal_hist = self.get(hist_type='signal', mass_window_index=mass_window_index)
        background_hists = self.get(hist_type='background', mass_window_index=mass_window_index)

        plotter = measurement_hist.plotter

        plotter.init_plotter(rows=2, cols=1)
        plotter.set_experiment_label(**{'year': measurement_hist.year})  # to avoid clipping of labels
        for bg_name, bg_hist in background_hists.items():
            plotter.add_hist(bg_hist, as_stack=True, as_denominator=True, **{'label': labels[bg_hist.get_label()]})
        plotter.add_hist(signal_hist, as_stack=True, as_denominator=True, **get_hist_kwargs(signal_hist.get_label()))
        plotter.add_hist(measurement_hist, **get_hist_kwargs(measurement_hist.get_label()))
        plotter.add_ratio_hists(location=(1, 0))

        plotter.draw_hist()
        # TODO define common cosmetics in class definition
        x_log_scale = True
        y_log_scale = True
        if self.is_pt:
            x_log_scale = False
            y_log_scale = False
        plotter.set_common_comparison_plot_cosmetics(self.x_axis_label, x_log_scale=x_log_scale, y_log_scale=y_log_scale)
        text = self.get_additional_text_on_plot(mass_window_index)
        plotter.add_text(text=text, location=(0, 0),  **{"frameon": False, "loc": "upper left", })

        plotter.save_and_reset_plotter(measurement_hist.hist_name)

    def draw_unfolded_level(self, mass_window_index=-1):

        measurement_hist = self.get(hist_type='unfolded_measurement', mass_window_index=mass_window_index)
        signal_hist = self.get(hist_type='unfolded_signal', mass_window_index=mass_window_index)

        plotter = measurement_hist.plotter

        plotter.init_plotter(rows=2, cols=1)
        plotter.set_experiment_label(**{'year': measurement_hist.year})
        plotter.add_hist(measurement_hist, **get_hist_kwargs(measurement_hist.get_label()))
        plotter.add_hist(signal_hist, as_stack=True, as_denominator=True,  **{'label': 'Drell-Yan'})
        plotter.add_ratio_hists(location=(1, 0))

        plotter.draw_hist()
        x_log_scale = True
        y_log_scale = True
        if self.is_pt:
            x_log_scale = False
            y_log_scale = False
        plotter.set_common_comparison_plot_cosmetics(self.x_axis_label, x_log_scale=x_log_scale, y_log_scale=y_log_scale)
        text = self.get_additional_text_on_plot(mass_window_index)
        plotter.add_text(text=text, location=(0, 0), **{"frameon": False, "loc": "upper left", })
        plotter.save_and_reset_plotter(measurement_hist.hist_name)

    def draw_acceptance_corrected_level(self, mass_window_index=-1,):
        measurement_hist = self.get(hist_type='acceptance_corrected_measurement', mass_window_index=mass_window_index)
        signal_hist = self.get(hist_type='acceptance_corrected_signal', mass_window_index=mass_window_index)

        plotter = measurement_hist.plotter

        plotter.init_plotter(rows=2, cols=1)
        plotter.set_experiment_label(**{'year': measurement_hist.year})
        plotter.add_hist(measurement_hist, **get_hist_kwargs(measurement_hist.get_label()))
        plotter.add_hist(signal_hist, as_stack=True, as_denominator=True,  **{'label': 'Drell-Yan'})
        plotter.add_ratio_hists(location=(1, 0))

        plotter.draw_hist()
        x_log_scale = True
        y_log_scale = True
        if self.is_pt:
            x_log_scale = False
            y_log_scale = False
        plotter.set_common_comparison_plot_cosmetics(self.x_axis_label, x_log_scale=x_log_scale, y_log_scale=y_log_scale)
        text = self.get_additional_text_on_plot(mass_window_index)
        plotter.add_text(text=text, location=(0, 0), **{"frameon": False, "loc": "upper left", })
        plotter.save_and_reset_plotter(measurement_hist.hist_name)
