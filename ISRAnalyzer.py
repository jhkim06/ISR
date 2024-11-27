from Analyzer import Analyzer
from Hist import Hist
from Acceptance import Acceptance


def change_to_greek(raw_string):
    if raw_string == 'mm':
        return "\mu\mu"
    else:
        return raw_string

class ISRAnalyzer(Analyzer):
    def __init__(self, data, signal, background,
                 mass_bins, pt_bins,
                 experiment='cms', year='2016', channel='ee'):
        super(ISRAnalyzer, self).__init__(data, signal, background,
                                          experiment, year, channel)
        self.mass_bins = mass_bins
        self.pt_bins = pt_bins
        # dipt_prefix, dimass_prefix
        self._set_isr_hist_names()
        self.systematics = None

    def _set_isr_hist_names(self):

        # 1D
        self.pt_hist_name_prefix = 'dipt_[reco__fine]_dimass'
        self.pt_matrix_name_prefix = 'dipt_[reco__fine]_[gen_dressed__coarse]_dimass'
        self.pt_fake_hist_name_prefix = 'dipt_[reco_gen_dressed_fake__fine]_dimass'
        self.pt_hist_full_phase_name_prefix = 'dipt_[gen_dressed_acceptance__coarse]_dimass'

        self.mass_hist_name_prefix = 'dimass_[reco__fine]_dipt'
        self.mass_matrix_name_prefix = 'dimass_[reco__fine]_[gen_dressed__coarse]_dipt'
        self.mass_fake_hist_name_prefix = 'dimass_[reco_gen_dressed_fake__fine]_dipt'
        self.mass_hist_full_phase_name_prefix = 'dimass_[gen_dressed_acceptance__coarse]_dipt'

        # 2D

    def get_hists_for_unfolding(self, hist_name_prefix, matrix_name_prefix, fake_hist_name_prefix,
                                bin_postfix=''):
        # input, matrix, bg, fake
        hist_name = hist_name_prefix + bin_postfix
        matrix_name = matrix_name_prefix + bin_postfix
        fake_hist_name = fake_hist_name_prefix + bin_postfix

        data_hist = self.get_data_hist(hist_name)
        response_matrix = self.get_signal_hist(matrix_name)
        fake_hist = self.get_signal_hist(fake_hist_name)
        background_hists = self.get_background_hist(hist_name)

        return data_hist, response_matrix, fake_hist, background_hists

    # basic plots
    def draw_isr_measurement_expectation_plot(self, hist_name_prefix, bin_postfix='', text="",
                                              bin_width_norm=False,
                                              x_variable_name='',
                                              y_log_scale=False, x_log_scale=False):
        hist_name = hist_name_prefix + bin_postfix
        self.draw_measurement_expectation_comparison_plot(hist_name, text=text,
                                                          bin_width_norm=bin_width_norm,
                                                          x_variable_name=x_variable_name,
                                                          y_log_scale=y_log_scale, x_log_scale=x_log_scale)

    def draw_isr_measurement_signal_plot(self, hist_name_prefix, bin_postfix='', text="",
                                         bin_width_norm=False,
                                         x_variable_name='',
                                         figsize=(8,8),
                                         y_log_scale=False, x_log_scale=False,
                                         fake_hist=None):
        hist_name = hist_name_prefix + bin_postfix
        if fake_hist:
            print(fake_hist['hist_name'] + bin_postfix)
            fake_hist['hist'] = self.get_signal_hist(fake_hist['hist_name'] + bin_postfix,
                                                     bin_width_norm=bin_width_norm,)
            del fake_hist['hist_name']

        self.draw_measurement_signal_comparison_plot(hist_name, text=text,
                                                     figsize=figsize,
                                                          bin_width_norm=bin_width_norm,
                                                          x_variable_name=x_variable_name,
                                                          y_log_scale=y_log_scale, x_log_scale=x_log_scale,
                                                     additional_hist=fake_hist)

    def draw_isr_measurement_expectation_plot_dipt(self, hist_name_prefix):
        for mass_bin in self.mass_bins:
            mass_bin_postfix = '_' + str(mass_bin[0]) + 'to' + str(mass_bin[1])

            self.draw_isr_measurement_expectation_plot(hist_name_prefix, mass_bin_postfix,
                                                       bin_width_norm=True,
                                                       text=str(int(mass_bin[0])) + "<m<" +
                                                            str(int(mass_bin[1])) + " GeV",
                                                       x_variable_name=
                                                       r'$p_{T}^{' + change_to_greek(self.channel) + '}$')

    def draw_isr_measurement_expectation_plot_dimass(self, hist_name_prefix):
        bin_postfix = '_' + str(self.pt_bins[0]) + 'to' + str(self.pt_bins[1])
        self.draw_isr_measurement_expectation_plot(hist_name_prefix, bin_postfix,
                                                   bin_width_norm=True,
                                                   text=str(r"$p_{T}$<" + str(int(self.pt_bins[1])) + " GeV"),
                                                   x_variable_name=r'$M^{' + change_to_greek(self.channel) + '}$',
                                                   y_log_scale=True, x_log_scale=True)

    def draw_isr_measurement_signal_plot_dimass(self, hist_name_prefix):
        bin_postfix = '_' + str(self.pt_bins[0]) + 'to' + str(self.pt_bins[1])

        # dimass_[reco_gen_dressed_fake__fine]_dipt
        fake_hist = {"hist_name": "dimass_[reco_gen_dressed_fake__fine]_dipt"}
        self.draw_isr_measurement_signal_plot(hist_name_prefix, bin_postfix,
                                                  bin_width_norm=True,
                                                  text=str(r"$p_{T}$<" + str(int(self.pt_bins[1])) + " GeV"),
                                                  x_variable_name=r'$M^{' + change_to_greek(self.channel) + '}$',
                                                  y_log_scale=True, x_log_scale=True,
                                                  fake_hist=fake_hist
                                                  )

    def draw_isr_measurement_signal_plot_dipt(self, hist_name_prefix):
        for mass_bin in self.mass_bins:
            mass_bin_postfix = '_' + str(mass_bin[0]) + 'to' + str(mass_bin[1])

            fake_hist = {"hist_name": "dipt_[reco_gen_dressed_fake__fine]_dimass"}
            self.draw_isr_measurement_signal_plot(hist_name_prefix, mass_bin_postfix,
                                                  bin_width_norm=True,
                                                  text=str(int(mass_bin[0])) + "<m<" +
                                                  str(int(mass_bin[1])) + " GeV",
                                                  x_variable_name=
                                                  r'$p_{T}^{' + change_to_greek(self.channel) + '}$',
                                                  fake_hist=fake_hist)

    def draw_isr_measurement_signal_plot_2d_dipt(self, hist_name):

        fake_hist = {"hist_name": "[tunfold-hist]_[dipt-dimass]_[reco_gen_dressed_fake__fine_O-window_v1_UO]"}
        self.draw_isr_measurement_signal_plot(hist_name,
                                              figsize=(15,8),
                                                  bin_width_norm=False,
                                                  text= "2D",
                                                  x_variable_name=
                                                  r'$p_{T},M bin index$',
                                                  fake_hist=fake_hist)

    def get_unfold_bin_maps(self, unfolded_bin_name, folded_bin_name):
        return self.signal[1].get_tobject(unfolded_bin_name), self.signal[1].get_tobject(folded_bin_name)

    def do_isr_unfold(self, data_hist, response_matrix, fake_hist, backgrounds,
                      unfolded_bin=None, folded_bin=None,
                      do_acceptance_correction=False, hist_full_phase_name=''):

        unfold_result = self.do_unfold(data_hist, response_matrix, fake_hist, backgrounds,
                                       unfolded_bin, folded_bin)

        if do_acceptance_correction:
            mc_hist_full_phase = self.get_signal_hist(hist_full_phase_name)
            mc_hist_acceptance = unfold_result.get_mc_truth_from_response_matrix()
            acceptance = Acceptance(mc_hist_full_phase, mc_hist_acceptance)

            acceptance_corrected = acceptance.do_correction(unfold_result.get_unfolded_hist())
            result = Hist(acceptance_corrected)
        else:
            result = Hist(unfold_result.get_unfolded_hist())

        # FIXME return TUnFolder object
        return result

    def get_unfolded_mean_pt_1d(self, pt_hist_name_prefix, pt_matrix_name_prefix,
                                pt_fake_hist_name_prefix='',
                                do_acceptance_correction=False, pt_hist_full_phase_name_prefix=''):
        pt_data = []
        for mass_bin in self.mass_bins:
            mass_bin_postfix = '_' + str(mass_bin[0]) + 'to' + str(mass_bin[1])
            # FIXME ISRAnalyzer knows the name of histogram to unfold, but Analyzer don't. But Analyzer
            # contains histograms for the analysis
            data_hist, response_matrix, fake_hist, backgrounds = self.get_hists_for_unfolding(pt_hist_name_prefix,
                                                                                              pt_matrix_name_prefix,
                                                                                              pt_fake_hist_name_prefix,
                                                                                              mass_bin_postfix)
            result = self.do_isr_unfold(data_hist, response_matrix, fake_hist, backgrounds,
                                        do_acceptance_correction=do_acceptance_correction,
                                        hist_full_phase_name=pt_hist_full_phase_name_prefix+mass_bin_postfix)

            pt_data.append(result.get_mean())
        return pt_data

    def get_unfolded_mean_pt_2d(self, pt_hist_name, pt_matrix_name,
                                folded_bin_name, unfolded_bin_name,
                                pt_fake_hist_name='',
                                do_acceptance_correction=False, pt_hist_full_phase_name=''):

        data_hist, response_matrix, fake_hist, backgrounds = self.get_hists_for_unfolding(pt_hist_name,
                                                                                          pt_matrix_name,
                                                                                          pt_fake_hist_name)
        # get bin map
        unfolded_bin, folded_bin = self.get_unfold_bin_maps(unfolded_bin_name, folded_bin_name)
        result = self.do_isr_unfold(data_hist, response_matrix, fake_hist, backgrounds,
                                    unfolded_bin=unfolded_bin, folded_bin=folded_bin,
                                    do_acceptance_correction=do_acceptance_correction,
                                    hist_full_phase_name=pt_hist_full_phase_name)

        pt_data = []
        for index, _ in enumerate(self.mass_bins):
            axis_steering = 'dipt[O];dimass[UOC' + str(index) + ']'
            temp_result = Hist(unfolded_bin.ExtractHistogram("",
                                                             result.raw_root_hist, 0, True, axis_steering))
            pt_data.append(temp_result.get_mean())
        return pt_data

    def get_unfolded_mean_mass(self, mass_hist_name_prefix, mass_matrix_name_prefix,
                               mass_fake_hist_name_prefix='',
                               do_acceptance_correction=False, mass_hist_full_phase_name_prefix=''):

        bin_postfix = '_' + str(self.pt_bins[0]) + 'to' + str(self.pt_bins[1])
        data_hist, response_matrix, fake_hist, backgrounds = self.get_hists_for_unfolding(mass_hist_name_prefix,
                                                                                          mass_matrix_name_prefix,
                                                                                          mass_fake_hist_name_prefix,
                                                                                          bin_postfix)
        result = self.do_isr_unfold(data_hist, response_matrix, fake_hist, backgrounds,
                                    do_acceptance_correction=do_acceptance_correction,
                                    hist_full_phase_name=mass_hist_full_phase_name_prefix+bin_postfix)

        mass_data = []
        for mass_bin in self.mass_bins:
            mass_data.append(result.get_mean(range_min=mass_bin[0], range_max=mass_bin[1]))

        return mass_data
