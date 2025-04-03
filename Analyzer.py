from Plotter import Plotter
from TUnFolder import TUnFolder


labels = {
    "gg": r"$\gamma\gamma$",
    "tau": r'$\tau\tau$',
    "ttbar": r'$t\bar{t}$',
    "VV": "$VV$"
}


class Analyzer:
    # plotter + sample group?
    def __init__(self, data, signal, background, experiment='cms', year='2016', channel='ee'):
        self.experiment = experiment

        self.year = year
        self.channel = channel

        # sample group
        self.data = data  # (Label, SampleGroup)
        self.signal = signal
        self.background = background  # [(Label, SampleGroup)]

    def set_base_hist_path(self, hist_path):
        self.data[1].set_hist_path_prefix(hist_path)
        self.signal[1].set_hist_path_prefix(hist_path)

    # FIXME
    def get_unfold_bin_maps(self, unfolded_bin_name, folded_bin_name):
        return self.signal[1].get_tobject(unfolded_bin_name), self.signal[1].get_tobject(folded_bin_name)

    def do_unfold(self, input_hist_name, matrix_name, fake_hist_name, bg_hist_name,
                  unfolded_bin_name=None, folded_bin_name=None, variable_name=''):

        data_hist = self.get_data_hist(input_hist_name)
        response_matrix = self.get_signal_hist(matrix_name)
        fake_hist = self.get_signal_hist(fake_hist_name)
        backgrounds = self.get_background_hist(bg_hist_name)

        unfolded_bin = None
        folded_bin = None
        if unfolded_bin_name and folded_bin_name:
            unfolded_bin, folded_bin = self.get_unfold_bin_maps(unfolded_bin_name, folded_bin_name)

        # draw folded histograms
        unfold = TUnFolder(response_matrix,
                           data_hist,
                           fake_hist,
                           bg_hists=backgrounds, unfolded_bin=unfolded_bin, folded_bin=folded_bin,
                           year=self.year, variable_name=variable_name)
        unfold.unfold()

        return unfold

    def draw_measurement_signal_comparison_plot(self, hist_name,
                                                figsize=(8,8),
                                                text='',
                                                bin_width_norm=False,
                                                x_variable_name='',
                                                y_log_scale=False, x_log_scale=False,
                                                additional_hist={},
                                                custom_x_labels=None,
                                                custom_x_locates=None,
                                                vlines=None):
        data_bg_subtracted = self.get_bg_subtracted_data(hist_name, bin_width_norm=bin_width_norm)
        signal_hist = self.get_signal_hist(hist_name, bin_width_norm=bin_width_norm)

        # comparison plot
        plotter = Plotter('CMS', './Plots')  # FIXME use self.plotter
        plotter.create_subplots(2, 1, figsize=figsize,
                                left=0.15, right=0.95, hspace=0.0, bottom=0.15, height_ratios=[1, 0.3])
        # expectations
        # measurement
        if additional_hist:
            plotter.add_hist(additional_hist['hist'], **{"color": 'cyan', 'yerr': False,
                                                         "label": 'Fake', "zorder": 999})
        plotter.add_hist(signal_hist, as_stack=True, **{ "color": 'red',
                                                         "label": 'Drell-Yan'})
        plotter.add_hist(data_bg_subtracted, **{"histtype": 'errorbar', "color": 'black',
                                                'label': 'Data'})
        if additional_hist:
            plotter.add_ratio_hist(nominator_index=0,
                                   denominator_index=1,  # add all simulation except data
                                   location=(1, 0), color='cyan')

        plotter.add_ratio_hist(nominator_index=2,
                               denominator_index=1,  # add all simulation except data
                               location=(1, 0), color='black', histtype='errorbar')
        plotter.draw_hist()
        plotter.set_experiment_label(**{"year": self.year})
        plotter.comparison_plot_cosmetics(x_variable_name, y_log_scale, x_log_scale, bin_width_norm)
        plotter.adjust_y_scale()
        if vlines:
            plotter.draw_vlines(vlines=vlines, location=(0, 0))
            plotter.draw_vlines(vlines=vlines, location=(1, 0))

        if custom_x_labels:
            plotter.add_custom_axis_tick_labels(custom_x_locates, custom_x_labels, location=(1, 0))
        plotter.add_text(text=text, location=(0,0), **{"frameon": False, "loc": "upper left",})  # Note: required to after setting legend
        plotter.save_fig(hist_name + "_bg_subtracted" + self.year)

    def draw_measurement_expectation_comparison_plot(self, hist_name,
                                                     text='',
                                                     bin_width_norm=False,
                                                     x_variable_name='',
                                                     y_log_scale=False, x_log_scale=False):

        data_hist = self.get_data_hist(hist_name, bin_width_norm=bin_width_norm)
        signal_hist = self.get_signal_hist(hist_name, bin_width_norm=bin_width_norm)
        background_hists = self.get_background_hist(hist_name, bin_width_norm=bin_width_norm)
        # total_expectation = self.get_total_expectation_hist(hist_name)

        plotter = Plotter('CMS', './Plots')  # FIXME use self.plotter
        plotter.create_subplots(2, 1,
                                left=0.15, right=0.95, hspace=0.0, bottom=0.15, height_ratios=[1, 0.3])
        plotter.set_experiment_label(**{"year": self.year})

        # expectations as stack
        # seems labels of stack written later
        for bg_label in background_hists:
            plotter.add_hist(background_hists[bg_label], as_stack=True, **{"label": labels[bg_label]})
        plotter.add_hist(signal_hist, as_stack=True, **{"color": 'red', "linestyle": '--',
                                                        "label": 'Drell-Yan'})
        # measurement
        plotter.add_hist(data_hist, **{"histtype": 'errorbar', "color": 'black',
                                       'label': 'Data'})
        plotter.add_ratio_hist(nominator_index=5,
                               denominator_index=[0,1,2,3,4],  # add all simulation except data
                               location=(1, 0), color='black', histtype='errorbar')
        # just to show x bin widths
        plotter.draw_hist()
        plotter.comparison_plot_cosmetics(x_variable_name, y_log_scale, x_log_scale, bin_width_norm)
        # plotter.get_axis(location=(0, 0)).set_ylim(ymin=1e-2)
        plotter.adjust_y_scale()
        plotter.add_text(text=text, location=(0,0), **{"frameon":False, "loc": "upper left",})  # Note: required to after setting legend
        plotter.save_fig(hist_name + "_" + self.year)

    def get_data_hist(self, hist_name, hist_path='', bin_width_norm=False):
        return self.data[1].get_combined_root_hists(hist_name, bin_width_norm=bin_width_norm)

    def get_signal_hist(self, hist_name, hist_path='', bin_width_norm=False):
        return self.signal[1].get_combined_root_hists(hist_name, bin_width_norm=bin_width_norm)

    def get_background_hist(self, hist_name, hist_path='', bin_width_norm=False):
        # return dictionary of root hists
        temp_dict = {}
        for bg in self.background:
            temp_dict[bg[0]] = bg[1].get_combined_root_hists(hist_name, bin_width_norm=bin_width_norm)
        return temp_dict

    def get_total_bg_hist(self, hist_name, hist_path='', bin_width_norm=False):
        total_bg = self.background[0][1].get_combined_root_hists(hist_name, bin_width_norm=bin_width_norm)
        for bg in self.background[1:]:
            total_bg.Add(bg[1].get_combined_root_hists(hist_name, bin_width_norm=bin_width_norm))
        return total_bg

    def get_total_expectation_hist(self, hist_name, hist_path='', bin_width_norm=False):
        total_expectation_hist = self.signal[1].get_combined_root_hists(hist_name, bin_width_norm=bin_width_norm)
        total_expectation_hist.Add(self.get_total_bg_hist(hist_name, hist_path, bin_width_norm=bin_width_norm))
        return total_expectation_hist

    def get_bg_subtracted_data(self, hist_name, hist_path='', bin_width_norm=False):
        raw_data = self.data[1].get_combined_root_hists(hist_name, bin_width_norm=bin_width_norm)
        total_bg = self.get_total_bg_hist(hist_name, hist_path, bin_width_norm=bin_width_norm)
        raw_data.Add(total_bg, -1)
        return raw_data
