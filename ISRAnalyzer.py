from Analyzer import Analyzer
from Hist import Hist
from Acceptance import Acceptance
import pandas as pd
import numpy as np


# TODO check if this is needed here
def calculate_squared_root_sum(raw_df, reg_expression, new_col_name="total error"):

    selected_cols = raw_df.filter(regex=reg_expression)
    squared_sum = (selected_cols ** 2).sum(axis=1)
    root_of_squared_sum = np.sqrt(squared_sum)

    raw_df[new_col_name] = root_of_squared_sum


def change_to_greek(raw_string):
    if raw_string == 'mm':
        return "\mu\mu"
    else:
        return raw_string


class ISRAnalyzer(Analyzer):
    def __init__(self,
                 sample_base_dir,
                 mass_bins,
                 pt_bins,
                 acceptance = None,  # DY for acceptance correction
                 folded_bin_name='fine', unfolded_bin_name='coarse',
                 unfolded_space_name='dressed', acceptance_space_name='dressed',
                 mass_folded_bin_name='fine', mass_unfolded_bin_name='coarse', mass_unfolded_space_name='dressed',):

        super(ISRAnalyzer, self).__init__(sample_base_dir,
                                          signal='DY', )

        if acceptance is None:
            self.acceptance = self.signal_name
        else:
            self.acceptance = acceptance
        self.mass_bins = mass_bins
        self.pt_bins = pt_bins

        # self.axis_steering = 'dipt[O];dimass[UOC' + str(index) + ']'
        # dipt_prefix, dimass_prefix

        self.dipt_label = r'$p_{T}$^{' + change_to_greek(self.channel) + '}'
        self.dimass_label = r'$M$^{' + change_to_greek(self.channel) + '}'

        self.folded_bin_name = folded_bin_name
        self.unfolded_bin_name = unfolded_bin_name
        self.unfolded_space_name = unfolded_space_name
        self.acceptance_space_name = acceptance_space_name

        self.mass_folded_bin_name = mass_folded_bin_name
        self.mass_unfolded_bin_name = mass_unfolded_bin_name

        # set histogram name of "pt and mass"
        self._set_isr_1d_hist_names()
        self._set_isr_2d_hist_names()

    def get_mass_bins(self):
        return self.mass_bins

    def _set_isr_1d_hist_names(self):

        # 1D case: set pt and mass histogram prefix
        self.pt_hist_name_prefix = 'dipt_[reco__' + self.folded_bin_name + ']_dimass'
        self.pt_fake_hist_name_prefix = ('dipt_[reco_gen_' + self.unfolded_space_name
                                         + '_fake__' + self.folded_bin_name + ']_dimass')
        self.pt_hist_full_phase_name_prefix = ('dipt_[gen_' + self.acceptance_space_name +
                                               '_acceptance__' + self.unfolded_bin_name + ']_dimass')

        self.pt_matrix_name_prefix = ('dipt_[reco__' + self.folded_bin_name + ']_[gen_' +
                                      self.unfolded_space_name + '__' + self.unfolded_bin_name + ']_dimass')

        self.mass_hist_name_prefix = 'dimass_[reco__' + self.mass_folded_bin_name + ']_dipt'
        self.mass_fake_hist_name_prefix = ('dimass_[reco_gen_' + self.unfolded_space_name
                                           + '_fake__' + self.mass_folded_bin_name + ']_dipt')
        self.mass_hist_full_phase_name_prefix = ('dimass_[gen_' + self.acceptance_space_name +
                                                 '_acceptance__' + self.mass_unfolded_bin_name + ']_dipt')

        self.mass_matrix_name_prefix = ('dimass_[reco__' + self.mass_folded_bin_name + ']_[gen_' +
                                        self.unfolded_space_name + '__' + self.mass_unfolded_bin_name + ']_dipt')

    def _set_isr_2d_hist_names(self):

        type_name = 'bin'
        var_name = 'dipt-dimass'
        tunfold_prefix = f"[tunfold-{type_name}]_[{var_name}]"
        bin_definition = 'O-window_v1_UO'
        # 2D case: set pt and mass histogram name
        self.pt_mass_detector_bin_postfix = "[reco__" + self.folded_bin_name + "_" + bin_definition + "]"
        self.pt_mass_unfolded_bin_postfix = (
                "[gen_" + self.unfolded_space_name + "__" + self.unfolded_bin_name + "_" + bin_definition + "]")
        self.pt_mass_detector_bin_name = tunfold_prefix + "_[" + self.folded_bin_name + "_" + bin_definition + "]"
        self.pt_mass_unfolded_bin_name = tunfold_prefix + "_[" + self.unfolded_bin_name + "_" + bin_definition + "]"

        # even if the same bin definition used different contents can be stored
        type_name = 'hist'
        tunfold_prefix = f"[tunfold-{type_name}]_[{var_name}]"
        self.pt_mass_hist_name = tunfold_prefix + "_" + self.pt_mass_detector_bin_postfix
        self.pt_mass_fake_hist_name = (tunfold_prefix + "_[reco_gen_" + self.unfolded_space_name
                                       + "_fake__" + self.folded_bin_name + "_" + bin_definition + "]")
        self.pt_mass_hist_full_phase_name = (tunfold_prefix + "_[gen_" + self.acceptance_space_name
                                             + "_acceptance__" + self.unfolded_bin_name + "_" + bin_definition + "]")

        type_name = 'matrix'
        tunfold_prefix = f"[tunfold-{type_name}]_[{var_name}]"
        self.pt_mass_matrix_name = (tunfold_prefix + "_" + self.pt_mass_detector_bin_postfix
                                    + "_" + self.pt_mass_unfolded_bin_postfix)

    def get_hist_names_for_1d_dipt(self, bin_postfix):

        # input, matrix, bg, fake
        input_hist_name = self.pt_hist_name_prefix + bin_postfix
        matrix_name = self.pt_matrix_name_prefix + bin_postfix
        fake_hist_name = self.pt_fake_hist_name_prefix + bin_postfix
        bg_hist_name = self.pt_hist_name_prefix + bin_postfix

        return input_hist_name, matrix_name, fake_hist_name, bg_hist_name

    def get_hist_names_for_1d_dimass(self, postfix):

        # input, matrix, bg, fake
        input_hist_name = self.mass_hist_name_prefix + postfix
        matrix_name = self.mass_matrix_name_prefix + postfix
        fake_hist_name = self.mass_fake_hist_name_prefix + postfix
        bg_hist_name = self.mass_hist_name_prefix + postfix

        return input_hist_name, matrix_name, fake_hist_name, bg_hist_name

    def get_bin_labels_locates_mass_window_boundaries(self, hist_name):

        custom_x_labels=None
        custom_x_locates=None
        vlines=None

        if '[tunfold-hist]' in hist_name:  # check if it is 2D
            _, folded_bin = self.get_unfold_bin_maps(self.pt_mass_unfolded_bin_name,
                                                     self.pt_mass_detector_bin_name)
            n_bins = folded_bin.GetDistributionNumberOfBins()
            custom_x_labels = []
            custom_x_locates = []
            vlines = []
            for i in range(1, n_bins+1):
                custom_x_locates.append(i)
                label = ':'.join(str(folded_bin.GetBinName(i)).split(':')[1:])
                custom_x_labels.append(label)
                if 'dipt[0,2]' in label:  # mass window edge
                    vlines.append(i-0.5)

        return custom_x_labels, custom_x_locates, vlines

    # basic plots
    def draw_isr_measurement_signal_plot(self,
                                         hist_name_prefix, bin_postfix='', text="", bin_width_norm=False,
                                         x_variable_name='',
                                         figsize=(8,8),
                                         y_log_scale=False,
                                         x_log_scale=False,
                                         fake_hist=None):

        hist_name = hist_name_prefix + bin_postfix

        self.draw_measurement_signal_comparison_plot(
            hist_name,
            figsize=figsize,
            bin_width_norm=bin_width_norm,
            x_variable_name=x_variable_name,
            y_log_scale=y_log_scale,
            x_log_scale=x_log_scale,
            save_and_reset=False)

        if fake_hist:
            fake_hist['hist'] = self.get_mc_hist(self.signal_name, fake_hist['hist_name'] + bin_postfix,
                                                 bin_width_norm=bin_width_norm,)
            self.plotter.add_hist(fake_hist['hist'], **{"color": 'cyan', 'yerr': False,
                                                        "label": 'Fake', "zorder": 999})
            self.plotter.draw_hist()
            self.plotter.set_common_comparison_plot_cosmetics(x_variable_name, y_log_scale, x_log_scale, bin_width_norm)
            self.plotter.adjust_y_scale()
            # self.plotter.add_ratio_hists(location=(1, 0))

        custom_x_labels, custom_x_locates, vlines = self.get_bin_labels_locates_mass_window_boundaries(hist_name)
        if vlines:
            self.plotter.draw_vlines(vlines=vlines, location=(0, 0))
            self.plotter.draw_vlines(vlines=vlines, location=(1, 0))
        if custom_x_labels:
            self.plotter.add_custom_axis_tick_labels(custom_x_locates, custom_x_labels, location=(1, 0))
        self.plotter.add_text(text=text, location=(0,0), **{"frameon": False, "loc": "upper left",})

        self.save_and_reset_plotter(hist_name, "bg_subtracted")

    def draw_isr_measurement_expectation_plot_dimass(self):
        bin_postfix = '_' + str(self.pt_bins[0]) + 'to' + str(self.pt_bins[1])
        hist_name = self.mass_hist_name_prefix + bin_postfix

        x_axis_label = r"$m^{" + change_to_greek(self.channel) +"}$"
        y_log_scale = True
        x_log_scale = True
        bin_width_norm = True
        text = str(r"$p_{T}^{"+change_to_greek(self.channel)+"}<$"+str(int(self.pt_bins[1]))+" (GeV)")

        self.draw_measurement_expectation_comparison_plot(
                                                          hist_name,
                                                          bin_width_norm=bin_width_norm,
                                                          x_axis_label=x_axis_label,
                                                          y_log_scale=y_log_scale,
                                                          x_log_scale=x_log_scale,
                                                          save_and_reset=False)

        self.plotter.add_text(text=text, location=(0, 0), **{"frameon": False, "loc": "upper left", })
        self.save_and_reset_plotter(hist_name)

    def draw_isr_measurement_expectation_plot_dipt(self,):
        for mass_bin in self.mass_bins:

            mass_bin_postfix = '_' + str(mass_bin[0]) + 'to' + str(mass_bin[1])
            hist_name = self.pt_hist_name_prefix + mass_bin_postfix
            x_axis_label = r"$p_{T}^{" + change_to_greek(self.channel) + "}$"
            bin_width_norm = True
            text = str(int(mass_bin[0]))+"$<m^{"+change_to_greek(self.channel)+"}<$"+str(int(mass_bin[1])) + " GeV"

            self.draw_measurement_expectation_comparison_plot(
                                                              hist_name,
                                                              bin_width_norm=bin_width_norm,
                                                              x_axis_label=x_axis_label,
                                                              save_and_reset=False)
            self.plotter.add_text(text=text, location=(0, 0), **{"frameon": False, "loc": "upper left", })
            self.save_and_reset_plotter(hist_name)

    def draw_isr_measurement_comparison_plot_dimass(self, *setups):

        bin_postfix = '_' + str(self.pt_bins[0]) + 'to' + str(self.pt_bins[1])
        hist_name = self.mass_hist_name_prefix + bin_postfix

        x_log_scale = True
        bin_width_norm = True
        norm = True
        channels = [ch for _, ch, _ in setups]
        if len(set(channels)) != 1:
            x_axis_label = r"$m^{ll}$"
            text = str(r"$p_{T}^{ll}<$"+str(int(self.pt_bins[1]))+" (GeV)")
            plot_postfix_channel = 'll'
        else:
            x_axis_label = r"$m^{" + change_to_greek(channels[0]) + "}$"
            text = str(r"$p_{T}^{" + change_to_greek(channels[0]) + "}<$" + str(int(self.pt_bins[1])) + " (GeV)")
            plot_postfix_channel = channels[0]

        # TODO show systematic band
        # error band for reference hist (denominator)
        self.draw_measurement_comparison_plot(*setups,
                                              hist_name=hist_name,
                                              bin_width_norm=bin_width_norm,
                                              norm=norm,
                                              x_axis_label=x_axis_label,
                                              x_log_scale=x_log_scale,
                                              save_and_reset=False)
        self.plotter.add_text(text=text, location=(0, 0), **{"frameon": False, "loc": "lower left", },
                              do_magic=False)
        years = [year for year, _, _ in setups]
        if len(set(years)) > 2:
            plot_postfix_year = 'years'
        else:
            plot_postfix_year = "_".join(years)
        self.save_and_reset_plotter(hist_name, plot_postfix_channel+"_"+plot_postfix_year)


    def draw_isr_measurement_signal_plot_dimass(self,):
        postfix = '_' + str(self.pt_bins[0]) + 'to' + str(self.pt_bins[1])

        # dimass_[reco_gen_dressed_fake__fine]_dipt
        fake_hist = {"hist_name": self.mass_fake_hist_name_prefix}
        self.draw_isr_measurement_signal_plot(
                                              self.mass_hist_name_prefix, postfix,
                                              bin_width_norm=True,
                                              text=str(r"$p_{T}<$" + str(int(self.pt_bins[1])) + " GeV"),
                                              x_variable_name=r"$m^{" + change_to_greek(self.channel) + "}$",
                                              y_log_scale=True, x_log_scale=True,
                                              fake_hist=fake_hist)

    def draw_isr_measurement_signal_plot_dipt(self):
        x_axis_label = r"$p_{T}^{" + change_to_greek(self.channel) + "}$"
        for mass_bin in self.mass_bins:
            mass_bin_postfix = '_' + str(mass_bin[0]) + 'to' + str(mass_bin[1])

            # show fake histogram also
            fake_hist = {"hist_name": self.pt_fake_hist_name_prefix}
            self.draw_isr_measurement_signal_plot(
                                                  self.pt_hist_name_prefix, mass_bin_postfix,
                                                  bin_width_norm=True,
                                                  text=str(int(mass_bin[0]))+"$<m<$"+str(int(mass_bin[1])) + " GeV",
                                                  x_variable_name=x_axis_label,
                                                  fake_hist=fake_hist,)

    def draw_isr_measurement_signal_plot_2d_dipt(self,):
        fake_hist = {"hist_name": self.pt_mass_fake_hist_name}
        self.draw_isr_measurement_signal_plot(
                                              self.pt_mass_hist_name,
                                              figsize=(15,8),
                                              bin_width_norm=False,
                                              text= "2D",
                                              x_variable_name=
                                              r'$p_{T}$, m bin index',
                                              fake_hist=fake_hist)

    def do_isr_unfold_sys(self):
        # loop over systematics in Hist
        # manual systematic for example use different simulation for response matrix, acceptance variation
        pass

    def do_isr_unfold(self, input_hist_name, matrix_name, fake_hist_name, bg_hist_name,
                      unfolded_bin_name=None, folded_bin_name=None,
                      do_acceptance_correction=False, hist_full_phase_name='',
                      variable_label='', bg_scale=1.0):
        #
        unfold_result = self.do_unfold(
            input_hist_name, matrix_name, fake_hist_name, bg_hist_name,
            unfolded_bin_name, folded_bin_name, variable_label, bg_scale=bg_scale)

        # TODO handle 2D "dipt[O];dimass[UOC1]"
        use_axis_binning = True
        if unfolded_bin_name is not None and folded_bin_name is not None:
            # for 2D, need to loop over each mass bins
            # FIXME change the output plot names according to bin definitions
            for index, _ in enumerate(self.mass_bins):
                axis_steering = 'dipt[O];dimass[UOC' + str(index) + ']'
                #
                unfold_result.bottom_line_test(draw_plot=True, out_name=input_hist_name+'_'+axis_steering,
                                               use_axis_binning=use_axis_binning, projection_mode=axis_steering)
            use_axis_binning = False
        else:
            unfold_result.bottom_line_test(draw_plot=True,
                                           out_name=input_hist_name,
                                           use_axis_binning=use_axis_binning)
        unfold_result.draw_response_matrix(out_name=input_hist_name)

        if do_acceptance_correction:
            mc_hist_full_phase = self.get_acceptance_hist(hist_full_phase_name)
            raw_response_matrix = self.get_acceptance_hist(matrix_name)
            mc_hist_acceptance = (
                unfold_result.extract_truth_from_2d_hist_using_bin_definition(raw_response_matrix.raw_root_hist))
            acceptance = Acceptance(mc_hist_full_phase, Hist(mc_hist_acceptance))
            acceptance_corrected = acceptance.do_correction(
                Hist(unfold_result.get_unfolded_hist(use_axis_binning=use_axis_binning)))
            result = acceptance_corrected
        else:
            result = Hist(unfold_result.get_unfolded_hist(use_axis_binning=use_axis_binning))

        # FIXME return TUnFolder object?
        return result

    def get_isr_measurement(self, values, as_df):
        if as_df:
            dict_list = []

            for ith_mass in range(len(self.mass_bins)):
                low_mass_edge = self.mass_bins[ith_mass][0]
                high_mass_edge = self.mass_bins[ith_mass][1]

                mass_window_str = str(low_mass_edge) + ':' + str(high_mass_edge)
                dict_temp = dict(mass_window=mass_window_str, mean=values[ith_mass][0])
                dict_temp['stat error'] = values[ith_mass][1]
                dict_list.append(dict_temp)

            df = pd.DataFrame(dict_list)

            reg_expression = r".*error"
            calculate_squared_root_sum(df, reg_expression)

            return df
        else:
            return values

    # TODO use postfix for systematic later
    def get_unfolded_mean_pt_1d(self, do_acceptance_correction=False, correct_binned_mean=False):
        # Use mass_bin as postfix
        pt_data = []
        for mass_bin in self.mass_bins:
            mass_bin_postfix = '_' + str(mass_bin[0]) + 'to' + str(mass_bin[1])
            # contains histograms for the analysis
            input_hist_name, matrix_name, fake_hist_name, bg_hist_name = (
                self.get_hist_names_for_1d_dipt(mass_bin_postfix))
            result = self.do_isr_unfold(input_hist_name, matrix_name, fake_hist_name, bg_hist_name,
                                        do_acceptance_correction=do_acceptance_correction,
                                        hist_full_phase_name=self.pt_hist_full_phase_name_prefix + mass_bin_postfix,
                                        variable_label=self.dipt_label)
            # get unbinned mean and binned mean
            pt_data.append(result.get_mean())
        if correct_binned_mean:
            self.correction_to_unbinned_prefsr_mean(pt_data, self.unfolded_space_name, self.unfolded_bin_name)
        return pt_data

    def get_sim_prefsr_mean_pt_1d(self, unfolded_space_name='', unfolded_bin_name=''):
        if unfolded_space_name == '':
            unfolded_space_name = self.unfolded_space_name
        if unfolded_bin_name == '':
            unfolded_bin_name = self.unfolded_bin_name
        pt_hist_full_phase_name_prefix = ('dipt_[gen_' + unfolded_space_name +
                                          '_acceptance__' + unfolded_bin_name + ']_dimass')
        pt_data = []
        for mass_bin in self.mass_bins:
            mass_bin_postfix = '_' + str(mass_bin[0]) + 'to' + str(mass_bin[1])

            mc_hist_full_phase = self.get_signal_hist(pt_hist_full_phase_name_prefix + mass_bin_postfix)
            pt_data.append(mc_hist_full_phase.get_mean(binned_mean=False))
        return pt_data

    def get_sim_prefsr_mean_pt_2d(self, unfolded_space_name='', unfolded_bin_name='',
                                  correct_binned_mean=False):
        if unfolded_space_name == '':
            unfolded_space_name = self.unfolded_space_name
        if unfolded_bin_name == '':
            unfolded_bin_name = self.unfolded_bin_name
        pt_mass_unfolded_bin_name = "[tunfold-bin]_[dipt-dimass]_[" + unfolded_bin_name + "_O-window_v1_UO]"
        pt_mass_hist_full_phase_name = ("[tunfold-hist]_[dipt-dimass]_[gen_" + unfolded_space_name +
                                        "_acceptance__" + unfolded_bin_name + "_O-window_v1_UO]")
        mc_hist_full_phase = self.get_signal_hist(pt_mass_hist_full_phase_name)

        unfolded_bin, _ = self.get_unfold_bin_maps(pt_mass_unfolded_bin_name, self.pt_mass_detector_bin_name)

        pt_data = self.extract_mean_pt_from_2d_hist(mc_hist_full_phase, unfolded_bin)
        if correct_binned_mean:
            self.correction_to_unbinned_prefsr_mean(pt_data, unfolded_space_name, unfolded_bin_name)
        return pt_data

    def correction_to_unbinned_prefsr_mean(self, pt_data, space_name, bin_name):
        pt_hist_full_phase_name_prefix = ('dipt_[gen_' + space_name +
                                          '_acceptance__' + bin_name + ']_dimass')

        for index, mass_bin in enumerate(self.mass_bins):
            mass_bin_postfix = '_' + str(mass_bin[0]) + 'to' + str(mass_bin[1])
            mc_hist_full_phase = self.get_mc_hist(self.signal_name, pt_hist_full_phase_name_prefix + mass_bin_postfix)
            correction_factor = (
                    mc_hist_full_phase.get_mean(binned_mean=False)[0]/mc_hist_full_phase.get_mean(binned_mean=True)[0])
            pt_data[index] = (pt_data[index][0] * correction_factor, pt_data[index][1])
        return pt_data

    def extract_mean_pt_from_2d_hist(self, pt_2d_hist, unfold_bin):
        pt_data = []
        for index, _ in enumerate(self.mass_bins):
            axis_steering = 'dipt[O];dimass[UOC' + str(index) + ']'
            temp_result = Hist(unfold_bin.ExtractHistogram("", pt_2d_hist.get_raw_hist(), 0,
                                                           True,
                                                           axis_steering))
            pt_data.append(temp_result.get_mean())
        return pt_data


    def get_unfolded_mean_pt_2d_systematic(self):
        pass

    def get_unfolded_mean_pt_2d(self,
                                do_acceptance_correction=False,
                                correct_binned_mean=False,
                                as_df=False, bg_scale=1.0, **kwargs):

        # set year, channel, event?
        unfolded_bin, _ = self.get_unfold_bin_maps(
            self.pt_mass_unfolded_bin_name,
            self.pt_mass_detector_bin_name)

        # use default histograms
        # TODO enable to use different hist name
        result = self.do_isr_unfold(self.pt_mass_hist_name, self.pt_mass_matrix_name,
                                    self.pt_mass_fake_hist_name, self.pt_mass_hist_name,
                                    bg_scale=bg_scale,
                                    unfolded_bin_name=self.pt_mass_unfolded_bin_name,
                                    folded_bin_name=self.pt_mass_detector_bin_name,
                                    do_acceptance_correction=do_acceptance_correction,
                                    hist_full_phase_name=self.pt_mass_hist_full_phase_name)

        pt_data = self.extract_mean_pt_from_2d_hist(result, unfolded_bin)
        if correct_binned_mean:
            self.correction_to_unbinned_prefsr_mean(pt_data, self.unfolded_space_name, self.unfolded_bin_name)

        return self.get_isr_measurement(pt_data, as_df)

    def get_unfolded_mean_mass(self,
                               do_acceptance_correction=False, as_df=False, bg_scale=1.0):

        postfix = '_' + str(self.pt_bins[0]) + 'to' + str(self.pt_bins[1])
        input_hist_name, matrix_name, fake_hist_name, bg_hist_name = self.get_hist_names_for_1d_dimass(postfix)
        result = self.do_isr_unfold(input_hist_name, matrix_name, fake_hist_name, bg_hist_name,
                                    bg_scale=bg_scale,
                                    do_acceptance_correction=do_acceptance_correction,
                                    hist_full_phase_name=self.mass_hist_full_phase_name_prefix + postfix,
                                    variable_label=r"$m^{" + change_to_greek(self.channel) + "}$")

        mass_data = []
        for mass_bin in self.mass_bins:
            mass_data.append(result.get_mean(range_min=mass_bin[0], range_max=mass_bin[1]))

        return self.get_isr_measurement(mass_data, as_df)

    def get_acceptance_hist(self, hist_name, hist_path='', bin_width_norm=False):
        return self.get_mc_hist(self.signal_name, hist_name, bin_width_norm=bin_width_norm)
