from Plotter import Plotter
from TUnFolder import TUnFolder

# simulation legends used in this analysis
labels = {
    # group_name of ROOTFileGroup: "legend"
    "DY": r"$Drell\!\!-\!\!Yan$",
    "gg": r"$\gamma\gamma$",
    "tau": r'$\tau\tau$',
    "ttbar": r'$t\bar{t}$',
    "vv": "$VV$"
}


class Analyzer:
    def __init__(self, data, signal, background, experiment='cms', year='2016', channel='ee'):
        self.experiment = experiment

        self.year = year
        self.channel = channel

        # sample group
        self.data = data  # ROOTFileGroup
        self.signal = signal
        self.background = background  # [ROOTFileGroup]

        # TODO
        # self.plotter

    def set_base_hist_path(self, hist_path):
        self.data.set_hist_path_prefix(hist_path)
        self.signal.set_hist_path_prefix(hist_path)

    # FIXME
    def get_unfold_bin_maps(self, unfolded_bin_name, folded_bin_name):
        return self.signal.get_tobject(unfolded_bin_name), self.signal.get_tobject(folded_bin_name)

    def do_unfold(self, input_hist_name, matrix_name, fake_hist_name, bg_hist_name,
                  unfolded_bin_name=None, folded_bin_name=None, variable_name=''):

        data_hist = self.get_data_hist(input_hist_name)
        response_matrix = self.get_signal_hist(matrix_name)
        fake_hist = self.get_signal_hist(fake_hist_name)
        backgrounds = self.get_background_hist(bg_hist_name, raw_hists=True)

        unfolded_bin = None
        folded_bin = None
        if unfolded_bin_name and folded_bin_name:
            unfolded_bin, folded_bin = self.get_unfold_bin_maps(unfolded_bin_name, folded_bin_name)

        # draw folded histograms
        unfold = TUnFolder(response_matrix.get_raw_hist(),
                           data_hist.get_raw_hist(),
                           fake_hist.get_raw_hist(),
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
        # data_bg_subtracted.hist, data_bg_subtracted.label
        # data.hist, data.label
        # data.year, data.lumi

        # signal.hist, signal.label
        data_bg_subtracted = self.get_bg_subtracted_data_hist(hist_name, bin_width_norm=bin_width_norm)
        signal_hist = self.get_signal_hist(hist_name, bin_width_norm=bin_width_norm)

        # comparison plot
        plotter = Plotter('CMS', './Plots')  # FIXME use self.plotter
        plotter.create_subplots(2, 1, figsize=figsize,
                                left=0.15, right=0.95, hspace=0.0, bottom=0.15, height_ratios=[1, 0.3])
        # expectations
        if additional_hist:
            plotter.add_hist(additional_hist['hist'], **{"color": 'cyan', 'yerr': False,
                                                         "label": 'Fake', "zorder": 999})
        signal_index = plotter.add_hist(signal_hist, as_stack=True, **{ "color": 'red',
                                                         "label": f'{labels[signal_hist.get_label()]}'})
        # measurement
        data_index = plotter.add_hist(data_bg_subtracted, **{"histtype": 'errorbar', "color": 'black',
                                                'label': f'{data_bg_subtracted.get_label()} (bkg. subtracted)'})
        if additional_hist:
            plotter.add_ratio_hist(nominator_index=0,
                                   denominator_index=signal_index,
                                   location=(1, 0), color='cyan')

        plotter.add_ratio_hist(nominator_index=data_index,
                               denominator_index=signal_index,
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
        denominator_index_list = []
        for _, bg in background_hists.items():
            index = plotter.add_hist(bg, as_stack=True, **{"label": labels[bg.get_label()]})
            denominator_index_list.append(index)
        index = plotter.add_hist(signal_hist, as_stack=True, **{"color": 'red', "linestyle": '--',
                                                                "label": f'{labels[signal_hist.get_label()]}'})
        denominator_index_list.append(index)
        # measurement
        nominator_index = plotter.add_hist(data_hist, **{"histtype": 'errorbar', "color": 'black',
                                       'label': f'{data_hist.get_label()}'})
        plotter.add_ratio_hist(nominator_index=nominator_index,
                               denominator_index=denominator_index_list,  # add all simulation except data
                               location=(1, 0), color='black', histtype='errorbar')
        # just to show x bin widths
        plotter.draw_hist()
        plotter.comparison_plot_cosmetics(x_variable_name, y_log_scale, x_log_scale, bin_width_norm)
        # plotter.get_axis(location=(0, 0)).set_ylim(ymin=1e-2)
        plotter.adjust_y_scale()
        plotter.add_text(text=text, location=(0,0), **{"frameon":False, "loc": "upper left",})  # Note: required to after setting legend
        plotter.save_fig(hist_name + "_" + self.year)

    def get_data_hist(self, hist_name, hist_path='', bin_width_norm=False):
        return self.data.get_combined_root_hists(hist_name, bin_width_norm=bin_width_norm)

    def get_signal_hist(self, hist_name, hist_path='', bin_width_norm=False):
        return self.signal.get_combined_root_hists(hist_name, bin_width_norm=bin_width_norm)

    def get_background_hist(self, hist_name, hist_path='', bin_width_norm=False, raw_hists=False):
        # return dictionary of root hists
        temp_dict = {}
        for bg in self.background:
            if raw_hists:
                temp_dict[bg.get_name()] = (
                    bg.get_combined_root_hists(hist_name, bin_width_norm=bin_width_norm).get_raw_hist())
            else:
                temp_dict[bg.get_name()] = bg.get_combined_root_hists(hist_name, bin_width_norm=bin_width_norm)
        return temp_dict

    def get_total_bg_hist(self, hist_name, hist_path='', bin_width_norm=False):
        # total_bg = self.background[0].get_combined_root_hists(hist_name, bin_width_norm=bin_width_norm)
        total_bg = None
        for bg in self.background:
            # total_bg.Add(bg.get_combined_root_hists(hist_name, bin_width_norm=bin_width_norm))
            total_bg = bg.get_combined_root_hists(hist_name, bin_width_norm=bin_width_norm) + total_bg
        return total_bg

    def get_total_expectation_hist(self, hist_name, hist_path='', bin_width_norm=False):
        total_expectation_hist = self.signal.get_combined_root_hists(hist_name, bin_width_norm=bin_width_norm)
        total_expectation_hist = (self.get_total_bg_hist(hist_name, hist_path, bin_width_norm=bin_width_norm) +
                                  total_expectation_hist)
        return total_expectation_hist

    def get_bg_subtracted_data_hist(self, hist_name, hist_path='', bin_width_norm=False):
        raw_data = self.data.get_combined_root_hists(hist_name, bin_width_norm=bin_width_norm)
        total_bg = self.get_total_bg_hist(hist_name, hist_path, bin_width_norm=bin_width_norm)
        raw_data = raw_data - total_bg
        return raw_data
