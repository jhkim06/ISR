from CMSData import CMSData
from Plotter import Plotter
from TUnFolder import TUnFolder


# simulation legends used in this analysis
labels = {
    # group_name of ROOTFileGroup: "legend"
    "Data": "Data",
    "DYJetsToEE_MiNNLO": r"Drell-Yan",
    "GGLL": r"$\gamma\gamma$",
    "DYJetsToTauTau_MiNNLO": r'$\tau\tau$',
    "TTLL": r'$t\bar{t}$',
    "vv": "$VV$",
    "ZZ": "$ZZ$",
    "WZ": "$WZ$",
    "WW": "$WW$",
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
                 backgrounds=['TTLL', 'GGLL', 'ZZ', 'WZ', 'WW', 'DYJetsToTauTau_MiNNLO'],
                 analysis_name=''):

        self.data = CMSData(sample_base_dir)

        self.experiment = self.data.experiment
        # Note set year channel event_selection before do anything
        self.year = ""
        self.channel = ""
        self.event_selection = ""

        self.signal_name = signal # process name ex) DY
        self.background_names = backgrounds

        self.analysis_name = analysis_name
        # TODO root file wrapper root_file().draw()
        self.plotter = Plotter(self.experiment,
                               '/Users/junhokim/Work/cms_snu/ISR/Plots')

        self.systematics = {
            # apply only to background Hist
            "bg_normalization:background": {"up": ("", 1.05), "down": ("", 0.95)},
        }
        #self.systematics = {}

    def set_data_info(self, year, channel, event_selection):
        self.year = year
        self.channel = channel
        self.event_selection = event_selection

    def reset_data_info(self):
        self.year = ""
        self.channel = ""
        self.event_selection = ""

    def plot_name_postfix(self, postfix=""):
        postfix = "_" + postfix + "_" if postfix else '_'
        return postfix + self.channel + "_" + self.year

    def get_measurement(self):
        return self.data.get_measurement(self.year, self.channel, self.event_selection)

    def get_mc(self, process_name, label=''):
        return self.data.get_mc(process_name, self.year, self.channel, self.event_selection, label)

    def set_systematics_on_hist(self, hist, file_group, hist_name, bin_width_norm=False):
        is_measurement = hist.is_measurement
        is_signal = hist.is_mc_signal

        for key, value in self.systematics.items():
            sys_name, apply_to = key.split(":")

            if apply_to == "all":
                use_sys_config = True
            elif apply_to == "measurement":
                use_sys_config = is_measurement
            elif apply_to == "signal":
                use_sys_config = not is_measurement and is_signal
            elif apply_to == "background":
                use_sys_config = not is_measurement and not is_signal
            else:
                use_sys_config = False  # fallback if unknown apply_to

            for variation_name, sys_config in value.items():
                if use_sys_config:
                    hist_name = hist_name+"_"+sys_config[0] if sys_config[0] else hist_name
                    sys_hist = file_group.get_combined_root_hists(hist_name,
                                                                  scale=sys_config[1],
                                                                  bin_width_norm=bin_width_norm)
                else:
                    sys_hist = file_group.get_combined_root_hists(hist_name, bin_width_norm=bin_width_norm)
                hist.set_systematic_hist(sys_name, variation_name, sys_hist.get_raw_hist())

    # Methods to get Hist
    def get_measurement_hist(self, hist_name, hist_path='', bin_width_norm=False, norm=False):
        file_group = self.get_measurement()

        hist = file_group.get_combined_root_hists(hist_name, bin_width_norm=bin_width_norm, norm=norm)
        self.set_systematics_on_hist(hist, file_group, hist_name, bin_width_norm=bin_width_norm)  # Analyzer knows systematics
        return hist

    def get_mc_hist(self, process_name, hist_name, hist_path='', bin_width_norm=False, scale=1.0, norm=False):
        file_group = self.get_mc(process_name)
        is_signal = False
        if process_name == self.signal_name:
            is_signal = True

        hist = file_group.get_combined_root_hists(hist_name, bin_width_norm=bin_width_norm, norm=norm,
                                                  is_signal=is_signal,
                                                  scale=scale)
        self.set_systematics_on_hist(hist, file_group, hist_name, bin_width_norm=bin_width_norm)
        return hist

    def get_background_hists(self, hist_name, hist_path='', bin_width_norm=False, bg_scale=1.0, norm=False):
        # return dictionary of root hists
        temp_dict = {}
        for bg in self.background_names:
            temp_dict[bg] = self.get_mc_hist(bg, hist_name, bin_width_norm=bin_width_norm,
                                             scale=bg_scale, norm=norm)
        return temp_dict

    def get_total_bg_hist(self, hist_name, hist_path='', bin_width_norm=False, norm=False):
        total_bg = None
        for name in self.background_names:
            bg = self.get_mc_hist(name, hist_name, bin_width_norm=bin_width_norm)
            total_bg = bg + total_bg

        if norm:
            total_bg.normalize()
        return total_bg

    def get_total_expectation_hist(self, hist_name, hist_path='', bin_width_norm=False, norm=False):
        file_group = self.get_measurement(self.signal_name)
        total_expectation_hist = file_group.get_combined_root_hists(hist_name, bin_width_norm=bin_width_norm)
        total_expectation_hist = (self.get_total_bg_hist(hist_name, hist_path, bin_width_norm=bin_width_norm) +
                                  total_expectation_hist)

        if norm:
            total_expectation_hist.normalize()
        return total_expectation_hist

    def get_bg_subtracted_measurement_hist(self, hist_name, hist_path='', bin_width_norm=False, norm=False):
        raw_data = self.get_measurement_hist(hist_name, hist_path, bin_width_norm=bin_width_norm)
        total_bg = self.get_total_bg_hist(hist_name, hist_path, bin_width_norm=bin_width_norm)
        raw_data = raw_data-total_bg

        if norm:
            raw_data.normalize()
        return raw_data

    # get TUnfoldBinning from root file
    def get_unfold_bin_maps(self, unfolded_bin_name, folded_bin_name):
        file_group = self.get_mc(self.signal_name)
        return file_group.get_tobject(unfolded_bin_name), file_group.get_tobject(folded_bin_name)

    # loop over Systematics in Hist
    def do_unfold_sys(self):
        pass

    def do_unfold(self,
                  input_hist_name,
                  matrix_name,
                  fake_hist_name,
                  bg_hist_name,
                  unfolded_bin_name=None, folded_bin_name=None, variable_name='',
                  bg_scale=1.0):

        # Systematic in Hist?
        data_hist = self.get_measurement_hist(input_hist_name)
        response_matrix = self.get_mc_hist(self.signal_name, matrix_name)  # FIXME properly handle TH2
        fake_hist = self.get_mc_hist(self.signal_name, fake_hist_name)
        backgrounds = self.get_background_hists(bg_hist_name, bg_scale=bg_scale)

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
                                                bin_width_norm=False,
                                                figsize=(8,8),
                                                x_variable_name='',
                                                y_log_scale=False,
                                                x_log_scale=False,
                                                save_and_reset=True):  # ?

        self.plotter.init_plotter(figsize=figsize, rows=2, cols=1)

        self.add_data_hist_to_plotter(hist_name, bin_width_norm=bin_width_norm, subtract_bg=True)
        self.add_signal_hist_to_plotter(hist_name, bin_width_norm=bin_width_norm, as_stack=True)

        self.plotter.add_ratio_hists(location=(1, 0))
        self.plotter.draw_hist()

        if save_and_reset:
            # Note if Hist added to plotter, legend might not be included
            self.plotter.set_common_comparison_plot_cosmetics(x_variable_name, y_log_scale, x_log_scale, bin_width_norm)
            self.plotter.adjust_y_scale()
            self.plotter.save_and_reset_plotter(hist_name, self.plot_name_postfix("bg_subtracted"))

    def draw_measurement_expectation_comparison_plot(self,
                                                     hist_name,
                                                     bin_width_norm=False,
                                                     x_axis_label='',
                                                     y_log_scale=False,
                                                     x_log_scale=False,
                                                     save_and_reset=True,):

        # check if Plotter is ready?
        self.plotter.init_plotter(figsize=(8,8), rows=2, cols=1)
        # add hists
        self.add_measurement_and_expectation_hists_to_plotter(hist_name,
                                                              bin_width_norm=bin_width_norm,)
        # TODO add option to reverse label order
        self.plotter.add_ratio_hists(location=(1, 0))

        # draw hists
        self.plotter.draw_hist()
        self.plotter.set_common_comparison_plot_cosmetics(x_axis_label, y_log_scale, x_log_scale, bin_width_norm)
        if save_and_reset:
            self.plotter.save_and_reset_plotter(hist_name)

    def draw_measurement_comparison_plot(self,
                                         *setups,
                                         hist_name,
                                         bin_width_norm=False,
                                         norm=False,
                                         x_axis_label='',
                                         x_log_scale=False,
                                         save_and_reset=True,):
        if len(setups) < 2:
            raise ValueError("At least two setups are required for comparison.")

        # Use first setup as reference (denominator in ratio)
        year_ref, channel_ref, event_sel_ref = setups[0]
        # self.plotter.init_plotter(figsize=(8, 8), rows=1, cols=1, year=" ".join([s[0] for s in setups]))
        self.plotter.init_plotter(figsize=(8, 8), rows=1, cols=1)

        # Add reference histogram
        self.set_data_info(year_ref, channel_ref, event_sel_ref)
        self.add_data_hist_to_plotter(hist_name, bin_width_norm=bin_width_norm,
                                      subtract_bg=True,
                                      location=-999,  # do not draw
                                      as_denominator=True,
                                      norm=norm,
                                      label=year_ref)
        self.reset_data_info()

        # Add comparison histograms
        for i, (year, channel, event_sel) in enumerate(setups[1:], start=1):
            self.set_data_info(year, channel, event_sel)
            self.add_data_hist_to_plotter(hist_name, bin_width_norm=bin_width_norm,
                                          subtract_bg=True,
                                          location=-999,  # do not draw
                                          as_denominator=False,
                                          norm=norm,
                                          label=year,
                                          color=f"C{i}")
            self.reset_data_info()

        self.plotter.add_ratio_hists(location=(0, 0))
        self.plotter.draw_hist()

        ratio_name = f"{','.join(s[0] for s in setups[1:])}/{year_ref}"
        self.plotter.set_common_ratio_plot_cosmetics(x_axis_label, x_log_scale=x_log_scale, y_axis_name=ratio_name,
                                                    y_min=0.5, y_max=1.5)
        if save_and_reset:
            combined_years = ''.join(s[0] for s in setups)
            self.plotter.save_and_reset_plotter(hist_name, channel_ref, self.plot_name_postfix(combined_years))

    # methods for Plotter
    def add_measurement_and_expectation_hists_to_plotter(self, hist_name,
                                                         bin_width_norm=False):
        self.add_background_hists_to_plotter(hist_name, bin_width_norm=bin_width_norm,
                                             as_stack=True, as_denominator=True)
        self.add_signal_hist_to_plotter(hist_name,
                                        bin_width_norm=bin_width_norm,
                                        as_stack=True, as_denominator=True)
        self.add_data_hist_to_plotter(hist_name,
                                      bin_width_norm=bin_width_norm, as_denominator=False)

    def add_data_hist_to_plotter(self, hist_name, bin_width_norm=False, subtract_bg=False,
                                 location=(0,0),
                                 as_denominator=False, norm=False, **kwargs):
        if subtract_bg:
            data_hist = self.get_bg_subtracted_measurement_hist(hist_name, bin_width_norm=bin_width_norm, norm=norm)
        else:
            data_hist = self.get_measurement_hist(hist_name, bin_width_norm=bin_width_norm, norm=norm)
        if not kwargs:
            kwargs = get_hist_kwargs(data_hist.get_label())
        else:
            kwargs =  get_hist_kwargs(data_hist.get_label()) | kwargs
        index = self.plotter.add_hist(data_hist, as_denominator=as_denominator, location=location,
                                      **kwargs)

        return index

    def add_signal_hist_to_plotter(self, hist_name, bin_width_norm=False, as_stack=False, as_denominator=True,
                                   norm=False):
        signal_hist = self.get_mc_hist(self.signal_name, hist_name, bin_width_norm=bin_width_norm, norm=norm)
        index = self.plotter.add_hist(signal_hist, as_stack=as_stack, as_denominator=as_denominator, 
                                      **get_hist_kwargs(signal_hist.get_label()))
        return index

    def add_background_hists_to_plotter(self, hist_name, bin_width_norm=False, as_stack=False, as_denominator=True,
                                        norm=False):
        background_hists = self.get_background_hists(hist_name, bin_width_norm=bin_width_norm, norm=norm)
        index_list = []
        for _, bg in background_hists.items():
            index = self.plotter.add_hist(bg, as_stack=as_stack, as_denominator=as_denominator,
                                          **{"label": labels[bg.get_label()]})  # mark as denominator index
            index_list.append(index)
        return index_list

