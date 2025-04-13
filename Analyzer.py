from CMSData import CMSData
from Plotter import Plotter
from TUnFolder import TUnFolder
from HistSystematic import HistSystematic


# simulation legends used in this analysis
labels = {
    # group_name of ROOTFileGroup: "legend"
    "Data": "Data",
    "DYJetsToEE_MiNNLO": r"Drell-Yan",
    "gg": r"$\gamma\gamma$",
    "DYJetsToTauTau_MiNNLO": r'$\tau\tau$',
    "TTLL": r'$t\bar{t}$',
    "vv": "$VV$"
}


colors = {
    "DYJetsToEE_MiNNLO": "red",
    "Data": "black"
}


def get_hist_kwargs(label):
    kwargs = {
        "color": f'{colors[label]}',
        "label": f'{labels[label]}'
    }
    if label == 'Data':
        kwargs.update({'histtype': 'errorbar'})

    return kwargs

def with_data_hist_info(func):
    def wrapper(self, *args, **kwargs):
        result = func(self, *args, **kwargs)
        self.reset_data_info()  # call reset after method finishes
        return result
    return wrapper


class Analyzer:
    def __init__(self, sample_base_dir, signal='DY',
                 backgrounds=['DYJetsToTauTau_MiNNLO', 'TTLL'],
                 analysis_name=''):

        self.data = CMSData(sample_base_dir)

        self.experiment = self.data.experiment
        self.year = ""
        self.channel = ""
        self.event_selection = ""

        # TODO use sample group!
        # self.cms_data
        # self.data = None
        self.signal = signal # process name ex) DY
        self.backgrounds = backgrounds

        self.analysis_name = analysis_name
        # self.year, self.channel,
        # TODO root file wrapper root_file().draw()
        self.plotter = Plotter(self.experiment,
                               '/Users/junhokim/Work/cms_snu/ISR/Plots')

    def set_data_info(self, year, channel, event_selection):
        self.year = year
        self.channel = channel
        self.event_selection = event_selection

    def reset_data_info(self):
        self.year = ""
        self.channel = ""
        self.event_selection = ""

    def get_data(self):
        return self.data.get_data(self.year, self.channel, self.event_selection)

    def get_mc(self, process_name, label=''):
        return self.data.get_mc(process_name, self.year, self.channel, self.event_selection, label)

    def plot_name_postfix(self, postfix=""):
        postfix = "_" + postfix + "_" if postfix else '_'
        return postfix + self.channel + "_" + self.year

    def get_unfold_bin_maps(self, unfolded_bin_name, folded_bin_name):
        file_group = self.get_mc(self.signal)
        return file_group.get_tobject(unfolded_bin_name), file_group.get_tobject(folded_bin_name)

    def do_unfold(self,
                  input_hist_name,
                  matrix_name,
                  fake_hist_name,
                  bg_hist_name,
                  unfolded_bin_name=None, folded_bin_name=None, variable_name='',
                  bg_scale=1.0):

        #
        data_hist = self.get_data_hist(input_hist_name)
        response_matrix = self.get_mc_hist(self.signal, matrix_name)
        fake_hist = self.get_mc_hist(self.signal, fake_hist_name)
        backgrounds = self.get_background_hist(bg_hist_name, bg_scale=bg_scale)

        unfolded_bin = None
        folded_bin = None
        if unfolded_bin_name and folded_bin_name:
            unfolded_bin, folded_bin = self.get_unfold_bin_maps(unfolded_bin_name, folded_bin_name)

        # draw folded histograms
        unfold = TUnFolder(response_matrix,
                           data_hist,
                           fake_hist,
                           bg_hists=backgrounds, unfolded_bin=unfolded_bin, folded_bin=folded_bin,
                           variable_name=variable_name)
        unfold.unfold()

        return unfold

    def draw_measurement_signal_comparison_plot(self, hist_name,
                                                figsize=(8,8),
                                                text='',
                                                bin_width_norm=False,
                                                x_variable_name='',
                                                y_log_scale=False,
                                                x_log_scale=False,
                                                additional_hist={},
                                                custom_x_labels=None,
                                                custom_x_locates=None,
                                                vlines=None):

        data_bg_subtracted = self.get_bg_subtracted_data_hist(hist_name, bin_width_norm=bin_width_norm)
        signal_hist = self.get_signal_hist(hist_name, bin_width_norm=bin_width_norm)

        self.init_plotter(figsize=figsize)

        # expectations
        if additional_hist:
            self.plotter.add_hist(additional_hist['hist'], **{"color": 'cyan', 'yerr': False,
                                                              "label": 'Fake', "zorder": 999})
        signal_index = self.plotter.add_hist(signal_hist, as_stack=True,
                                             **get_hist_kwargs(signal_hist.get_label()))
        # measurement
        data_index = self.plotter.add_hist(data_bg_subtracted,
                                           **get_hist_kwargs(data_bg_subtracted.get_label()))
        if additional_hist:
            self.plotter.add_ratio_hist(nominator_index=0,
                                        denominator_index=signal_index,
                                        location=(1, 0))

        self.plotter.add_ratio_hist(nominator_index=data_index,
                                    denominator_index=signal_index,
                                    location=(1, 0))
        self.plotter.draw_hist()
        self.plotter.comparison_plot_cosmetics(x_variable_name, y_log_scale, x_log_scale, bin_width_norm)
        self.plotter.adjust_y_scale()

        if vlines:
            self.plotter.draw_vlines(vlines=vlines, location=(0, 0))
            self.plotter.draw_vlines(vlines=vlines, location=(1, 0))
        if custom_x_labels:
            self.plotter.add_custom_axis_tick_labels(custom_x_locates, custom_x_labels, location=(1, 0))

        self.plotter.add_text(text=text, location=(0,0), **{"frameon": False, "loc": "upper left",})  # Note: required to after setting legend
        self.plotter.save_fig(self.plot_name_postfix(hist_name + "_bg_subtracted"))
        self.plotter.reset()

    @with_data_hist_info
    def draw_measurement_expectation_comparison_plot(self,
                                                     year, channel, event_selection,
                                                     hist_name,
                                                     text='',
                                                     bin_width_norm=False,
                                                     x_variable_name='',
                                                     y_log_scale=False, x_log_scale=False):

        self.set_data_info(year, channel, event_selection)
        data_hist = self.get_data_hist(hist_name, bin_width_norm=bin_width_norm)
        signal_hist = self.get_mc_hist(self.signal, hist_name, bin_width_norm=bin_width_norm)
        background_hists = self.get_background_hist(hist_name, bin_width_norm=bin_width_norm)

        # total_expectation = self.get_total_expectation_hist(hist_name)
        self.init_plotter(figsize=(8,8))

        # expectations as stack
        denominator_index_list = []
        for _, bg in background_hists.items():
            index = self.plotter.add_hist(bg, as_stack=True, **{"label": labels[bg.get_label()]})
            denominator_index_list.append(index)
        # self.add_hist_to_plotter(signal_hist, as_stack=True)
        index = self.plotter.add_hist(signal_hist, as_stack=True, **get_hist_kwargs(signal_hist.get_label()))
        denominator_index_list.append(index)
        # measurement
        nominator_index = self.plotter.add_hist(data_hist, **get_hist_kwargs(data_hist.get_label()))
        self.plotter.add_ratio_hist(nominator_index=nominator_index,
                                    denominator_index=denominator_index_list,
                                    location=(1, 0))
        # just to show x bin widths
        # self.draw_plot()
        self.plotter.draw_hist()
        self.plotter.comparison_plot_cosmetics(x_variable_name, y_log_scale, x_log_scale, bin_width_norm)
        # plotter.get_axis(location=(0, 0)).set_ylim(ymin=1e-2)
        self.plotter.adjust_y_scale()

        self.plotter.add_text(text=text, location=(0,0), **{"frameon":False, "loc": "upper left",})  # Note: required to after setting legend
        self.plotter.save_fig(self.plot_name_postfix(hist_name))
        self.plotter.reset()

    # methods for Plotter
    def init_plotter(self, figsize):
        self.plotter.create_subplots(2, 1, figsize=figsize,
                                     left=0.15, right=0.95, hspace=0.0, bottom=0.15, height_ratios=[1, 0.3])
        self.plotter.set_experiment_label(**{"year": self.year})

    def get_data_hist(self, hist_name, hist_path='', bin_width_norm=False):
        file_group = self.get_data()
        return file_group.get_combined_root_hists(hist_name, bin_width_norm=bin_width_norm)

    def get_mc_hist(self, process_name, hist_name, hist_path='', bin_width_norm=False, scale=1.0):
        file_group = self.get_mc(process_name)
        return file_group.get_combined_root_hists(hist_name, bin_width_norm=bin_width_norm, scale=scale)

    def get_background_hist(self, hist_name, hist_path='', bin_width_norm=False, bg_scale=1.0):
        # return dictionary of root hists
        temp_dict = {}
        for bg in self.backgrounds:
            temp_dict[bg] = self.get_mc_hist(bg, hist_name, bin_width_norm=bin_width_norm,
                                             scale=bg_scale)
        return temp_dict

    def get_total_bg_hist(self, hist_name, hist_path='', bin_width_norm=False):
        total_bg = None
        for bg in self.backgrounds:
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
