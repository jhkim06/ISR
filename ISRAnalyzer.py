from Analyzer import Analyzer
from Acceptance import Acceptance
from Hist import Hist, change_to_greek
from ISRHists import ISRHists
from ISR2DHist import ISR2DHist
from TUnFolder import TUnFolder


class ISRAnalyzer(Analyzer):
    def __init__(self,
                 sample_base_dir,
                 mass_bins,
                 pt_bins,
                 signal = "DY",
                 acceptance = None,  # DY for acceptance correction

                 pt_folded_bin_name='fine_O', pt_unfolded_bin_name='coarse_O',
                 unfolded_space_name='dressed', acceptance_space_name='dressed',

                 mass_folded_bin_name='fine_O', mass_unfolded_bin_name='coarse_O',
                 mass_unfolded_space_name='dressed', ):

        super(ISRAnalyzer, self).__init__(sample_base_dir,
                                          signal=signal, )

        if acceptance is None:
            self.acceptance_name = self.signal_name
        else:
            self.acceptance_name = acceptance
        self.mass_bins = mass_bins
        self.pt_bins = pt_bins

        # self.axis_steering = 'dipt[O];dimass[UOC' + str(index) + ']'
        # dipt_prefix, dimass_prefix
        self.dipt_label = r'$p_{T}$^{'+change_to_greek(self.channel)+'}'
        self.dimass_label = r'$M$^{'+change_to_greek(self.channel)+'}'

        # bin names for unfolding
        self.pt_folded_bin_name = pt_folded_bin_name
        self.pt_unfolded_bin_name = pt_unfolded_bin_name
        self.unfolded_space_name = unfolded_space_name
        self.acceptance_space_name = acceptance_space_name

        self.mass_folded_bin_name = mass_folded_bin_name
        self.mass_unfolded_bin_name = mass_unfolded_bin_name

        # set histogram name of "pt and mass"
        self._set_isr_1d_hist_names()
        self._set_isr_2d_hist_names()

        self.isr_pt = None
        self.isr_mass = None

    def get_mass_bins(self):
        return self.mass_bins

    def get_isr_results(self):
        return self.isr_pt, self.isr_mass

    def _set_isr_1d_hist_names(self):

        # 1D case: set pt and mass histogram prefix
        self.pt_hist_name_prefix = 'dipt_[reco__' + self.pt_folded_bin_name + ']_dimass'
        self.pt_fake_hist_name_prefix = ('dipt_[reco_gen_' + self.unfolded_space_name
                                         + '_fake__' + self.pt_folded_bin_name + ']_dimass')
        self.pt_hist_full_phase_name_prefix = ('dipt_[gen_' + self.acceptance_space_name +
                                               '_acceptance__' + self.pt_unfolded_bin_name + ']_dimass')
        self.pt_matrix_name_prefix = ('dipt_[reco__' + self.pt_folded_bin_name + ']_[gen_' +
                                      self.unfolded_space_name + '__' + self.pt_unfolded_bin_name + ']_dimass')

        self.pt_detector_bin_name = '[tunfold-bin]_dipt_[' + self.pt_folded_bin_name + ']'
        self.pt_truth_bin_name = '[tunfold-bin]_dipt_[' + self.pt_unfolded_bin_name + ']'

        self.mass_hist_name_prefix = 'dimass_[reco__' + self.mass_folded_bin_name + ']_dipt'
        self.mass_fake_hist_name_prefix = ('dimass_[reco_gen_' + self.unfolded_space_name
                                           + '_fake__' + self.mass_folded_bin_name + ']_dipt')
        self.mass_hist_full_phase_name_prefix = ('dimass_[gen_' + self.acceptance_space_name +
                                                 '_acceptance__' + self.mass_unfolded_bin_name + ']_dipt')
        self.mass_matrix_name_prefix = ('dimass_[reco__' + self.mass_folded_bin_name + ']_[gen_' +
                                        self.unfolded_space_name + '__' + self.mass_unfolded_bin_name + ']_dipt')

        self.mass_detector_bin_name = '[tunfold-bin]_dimass_['+self.mass_folded_bin_name+']'
        self.mass_truth_bin_name = '[tunfold-bin]_dimass_[' + self.mass_unfolded_bin_name + ']'

    def _set_isr_2d_hist_names(self):

        type_name = 'bin'
        var_name = 'dipt-dimass'
        tunfold_prefix = f"[tunfold-{type_name}]_[{var_name}]"
        bin_definition = '-window_v1_UO'
        # 2D case: set pt and mass histogram name
        self.pt_mass_detector_bin_postfix = "[reco__" + self.pt_folded_bin_name + bin_definition + "]"
        self.pt_mass_unfolded_bin_postfix = (
                "[gen_" + self.unfolded_space_name + "__" + self.pt_unfolded_bin_name + bin_definition + "]")
        self.pt_mass_detector_bin_name = tunfold_prefix + "_[" + self.pt_folded_bin_name + bin_definition + "]"
        self.pt_mass_truth_bin_name = tunfold_prefix + "_[" + self.pt_unfolded_bin_name + bin_definition + "]"

        # even if the same bin definition is used, different contents can be stored
        type_name = 'hist'
        tunfold_prefix = f"[tunfold-{type_name}]_[{var_name}]"
        self.pt_mass_hist_name = tunfold_prefix + "_" + self.pt_mass_detector_bin_postfix
        self.pt_mass_fake_hist_name = (tunfold_prefix + "_[reco_gen_" + self.unfolded_space_name
                                       + "_fake__" + self.pt_folded_bin_name + bin_definition + "]")
        self.pt_mass_hist_full_phase_name = (tunfold_prefix + "_[gen_" + self.acceptance_space_name
                                             + "_acceptance__" + self.pt_unfolded_bin_name + bin_definition + "]")

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
            _, folded_bin = self.get_unfold_bin_maps(self.pt_mass_truth_bin_name,
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

    def setup_isr_detector_hists(self, year, channel, event_selection, use_2d_pt=True,
                                 hist_prefix=''):
        self.set_data_info(year, channel, event_selection)
        # if use_2d_pt, no need to loop over mass windows
        if use_2d_pt:  # FIXME use_bin_map
            unfolded_bin, folded_bin = self.get_unfold_bin_maps(
                self.pt_mass_truth_bin_name,
                self.pt_mass_detector_bin_name
            )
            self.isr_pt = ISRHists(self.mass_bins, self.pt_bins, is_2d=True, is_pt=True,
                                   year=year, channel=channel,
                                   folded_tunfold_bin=folded_bin, unfolded_tunfold_bin=unfolded_bin)
            measurement, signal, signal_fake, bgs, matrix = self.get_isr_hist(hist_prefix=hist_prefix)
            self.isr_pt.set_isr_hists(measurement, signal, signal_fake, bgs, matrix)
        else:
            unfolded_bin, folded_bin = self.get_unfold_bin_maps(
                self.pt_truth_bin_name,
                self.pt_detector_bin_name
            )
            self.isr_pt = ISRHists(self.mass_bins, self.pt_bins, is_2d=False, is_pt=True,
                                   year=year, channel=channel,
                                   folded_tunfold_bin=folded_bin, unfolded_tunfold_bin=unfolded_bin)
            for index, _ in enumerate(self.mass_bins):
                measurement, signal, signal_fake, bgs, matrix = self.get_isr_hist(is_pt=True, is_2d=False,
                                                                                  mass_window_index=index,
                                                                                  hist_prefix=hist_prefix)
                self.isr_pt.set_isr_hists(measurement, signal, signal_fake, bgs, matrix)

        unfolded_bin, folded_bin = self.get_unfold_bin_maps(
            self.mass_truth_bin_name,
            self.mass_detector_bin_name
        )
        # for mass, always 1D unfolding
        self.isr_mass = ISRHists(self.mass_bins, self.pt_bins, is_2d=False, is_pt=False,
                                 year=year, channel=channel,
                                 folded_tunfold_bin=folded_bin, unfolded_tunfold_bin=unfolded_bin)
        measurement, signal, signal_fake, bgs, matrix = self.get_isr_hist(is_pt=False, hist_prefix=hist_prefix)
        self.isr_mass.set_isr_hists(measurement, signal, signal_fake, bgs, matrix)

    def setup_isr_acceptance_hists(self, year, channel, event_selection, is_2d=True, mass_bins=None,
                                   binned_mean=True, range_min=None, range_max=None,):
        self.set_data_info(year, channel, event_selection)

        if mass_bins is None:
            mass_bins = self.mass_bins

        # pt
        if is_2d:
            pass
        else:
            self.isr_pt = ISRHists(mass_bins, self.pt_bins, is_2d=False, is_pt=True,
                                   year=year, channel=channel, )
            for index, mass_bin in enumerate(mass_bins):
                mass_bin_postfix = '_' + str(mass_bin[0]) + 'to' + str(mass_bin[1])

                mc_hist_full_phase = self.get_acceptance_hist(self.pt_hist_full_phase_name_prefix + mass_bin_postfix)
                self.isr_pt.set_isr_hists(acceptance_corrected_signal_hist=mc_hist_full_phase)
                self.isr_pt.set_acceptance_corrected_mean_values(mass_window_index=index, key='simulation',
                                                                 binned_mean=binned_mean,
                                                                 range_min=range_min, range_max=range_max)

            self.isr_pt.binned_mean_correction_factors = self.get_correction_factors(
                self.unfolded_space_name, self.pt_unfolded_bin_name
            )

        # mass
        postfix = '_' + str(self.pt_bins[0]) + 'to' + str(self.pt_bins[1])  # FIXME get this info from self.isr_mass
        hist_full_phase_name = self.mass_hist_full_phase_name_prefix + postfix
        mc_hist_full_phase = self.get_acceptance_hist(hist_full_phase_name)

        self.isr_mass = ISRHists(mass_bins, self.pt_bins, is_2d=False, is_pt=False,
                                 year=year, channel=channel)
        self.isr_mass.set_isr_hists(acceptance_corrected_signal_hist=mc_hist_full_phase)
        self.isr_mass.set_acceptance_corrected_mean_values(key='simulation')

    def get_isr_hist(self, is_pt=True, is_2d=True, mass_window_index=0,
                     hist_prefix=''):
        if is_pt:
            if is_2d:
                hist_name = self.pt_mass_hist_name
                fake_hist_name = self.pt_mass_fake_hist_name
                matrix_name = self.pt_mass_matrix_name
            else:
                mass_bin_postfix = '_' + str(self.mass_bins[mass_window_index][0]) + 'to' + str(self.mass_bins[mass_window_index][1])
                hist_name = self.pt_hist_name_prefix+mass_bin_postfix
                fake_hist_name = self.pt_fake_hist_name_prefix+mass_bin_postfix
                matrix_name = self.pt_matrix_name_prefix+mass_bin_postfix

        else:
            postfix = '_' + str(self.pt_bins[0]) + 'to' + str(self.pt_bins[1])
            hist_name, matrix_name, fake_hist_name, _ = self.get_hist_names_for_1d_dimass(postfix)

        measurement = self.get_measurement_hist(hist_name, hist_name_prefix=hist_prefix)
        signal = self.get_mc_hist(self.signal_name, hist_name, hist_name_prefix=hist_prefix)
        bgs = self.get_background_hists(hist_name, hist_name_prefix=hist_prefix)

        signal_fake = self.get_mc_hist(self.signal_name, fake_hist_name)
        matrix = self.get_mc_hist(self.signal_name, matrix_name)
        return measurement, signal, signal_fake, bgs, matrix

    def isr_unfolds(self):
        self.pt_isr_unfold()
        self.mass_isr_unfold()

    def pt_isr_unfold(self):
        def run_unfold(mass_window_index, is2d=False):
            # Fetch ISR inputs
            input_hist, signal_hist, signal_fake_hist, background_hists, matrix = (
                self.isr_pt.get_isr_hists(mass_window_index=mass_window_index)
            )

            unfolded_bin = self.isr_pt.unfolded_tunfold_bin
            folded_bin = self.isr_pt.folded_tunfold_bin

            # Perform unfolding
            unfold = TUnFolder(matrix, input_hist, signal_fake_hist,
                               bg_hists=background_hists,
                               unfolded_bin=unfolded_bin,
                               folded_bin=folded_bin,
                               variable_name='')
            unfold.unfold()
            # closure (simple)
            closure = TUnFolder(matrix, signal_hist, signal_fake_hist,
                               unfolded_bin=unfolded_bin,
                               folded_bin=folded_bin,
                               variable_name='')
            closure.unfold()

            self.isr_pt.isr_hists[mass_window_index].tunfolder = unfold

            # create unfold input
            unfold_input_hist = input_hist.create(
                hist=unfold.get_input_hist(),
                label=input_hist.label + '(unfold input)'
            )
            # Create unfolded measurement
            unfolded_hist = input_hist.create(
                hist=unfold.get_unfolded_hist(),
                label=input_hist.label + '(unfolded)'
            )
            unfolded_hist.systematic_raw_root_hists = unfold.sys_unfold()
            unfolded_hist.compute_systematic_rss_per_sysname()
            # Create truth signal
            truth_signal_hist = signal_fake_hist.create(
                hist=unfold.get_mc_truth_from_response_matrix(use_axis_binning=False),
                hist_name="_", label='Truth DY',
            )
            unfolded_signal_hist = input_hist.create(
                hist=closure.get_unfolded_hist(),
                hist_name="_", label='Unfolded DY',
            )

            self.isr_pt.isr_hists[mass_window_index].unfold_input_hist = (
                ISR2DHist(unfold_input_hist, folded_bin))
            self.isr_pt.isr_hists[mass_window_index].unfolded_measurement_hist = (
                ISR2DHist(unfolded_hist, unfolded_bin))
            self.isr_pt.isr_hists[mass_window_index].truth_signal_hist = (
                ISR2DHist(truth_signal_hist, unfolded_bin))
            self.isr_pt.isr_hists[mass_window_index].unfolded_signal_hist = (
                ISR2DHist(unfolded_signal_hist, unfolded_bin))

        if self.isr_pt.is_2d:
            run_unfold(0, is2d=True)
        else:
            for index, _ in enumerate(self.mass_bins):
                run_unfold(index)

    def mass_isr_unfold(self, tau=0.0):
        input_hist, signal_hist, signal_fake_hist, background_hists, matrix = self.isr_mass.get_isr_hists()

        unfolded_bin = self.isr_mass.unfolded_tunfold_bin
        folded_bin = self.isr_mass.folded_tunfold_bin

        # TODO set bin map
        unfold = TUnFolder(matrix,
                           input_hist,
                           signal_fake_hist,
                           bg_hists=background_hists, variable_name='',
                           folded_bin=folded_bin, unfolded_bin=unfolded_bin)

        unfold.unfold(tau=tau)
        self.isr_mass.isr_hists[0].tunfolder = unfold

        # closure (simple)
        closure = TUnFolder(matrix, signal_hist, signal_fake_hist,
                            variable_name='',
                            folded_bin=folded_bin, unfolded_bin=unfolded_bin)
        closure.unfold()

        unfold_input_hist = input_hist.create(
            hist=unfold.get_input_hist(),
            label=input_hist.label + '(unfold input)'
        )

        # Check difference input with 1) bin_width_norm=True 2) False
        # seems unfolding input should be bin_width_norm=False
        unfolded_hist = input_hist.create(hist=unfold.get_unfolded_hist(),
                                          label=input_hist.label + '(unfolded)')
        sys_unfolded_raw_hist = unfold.sys_unfold()
        unfolded_hist.systematic_raw_root_hists = sys_unfolded_raw_hist
        unfolded_hist.compute_systematic_rss_per_sysname()

        unfolded_signal_hist = signal_fake_hist.create(
            hist=closure.get_unfolded_hist(),
            hist_name="_", label='Unfolded DY',
        )

        raw_hist = unfold.get_mc_truth_from_response_matrix(use_axis_binning=False)
        # raw_hist.Scale(1, "width")
        truth_signal_hist = signal_fake_hist.create(
            hist=raw_hist,
            hist_name='_',
            label='Truth DY'
        )
        self.isr_mass.isr_hists[0].unfold_input_hist =  ISR2DHist(unfold_input_hist, folded_bin)
        self.isr_mass.isr_hists[0].unfolded_measurement_hist = ISR2DHist(unfolded_hist, unfolded_bin)
        self.isr_mass.isr_hists[0].truth_signal_hist = ISR2DHist(truth_signal_hist, unfolded_bin)
        self.isr_mass.isr_hists[0].unfolded_signal_hist = ISR2DHist(unfolded_signal_hist, unfolded_bin)

    def isr_acceptance_corrections(self):
        self.pt_isr_acceptance_correction()
        self.mass_isr_acceptance_correction()

    def pt_isr_acceptance_correction(self):
        def run_acceptance_correction(index=0, is2d=False):
            if is2d:
                mc_hist_full_phase = self.get_acceptance_hist(self.pt_mass_hist_full_phase_name)
                matrix_name = self.pt_mass_matrix_name
            else:
                mass_bin_postfix = f'_{self.mass_bins[index][0]}to{self.mass_bins[index][1]}'
                mc_hist_full_phase = self.get_acceptance_hist(self.pt_hist_full_phase_name_prefix + mass_bin_postfix)

                mass_bin_postfix = '_' + str(self.mass_bins[index][0]) + 'to' + str(self.mass_bins[index][1])
                matrix_name = self.pt_matrix_name_prefix+mass_bin_postfix

            raw_response_matrix = self.get_acceptance_hist(matrix_name)
            mc_acceptance_hist = Hist(
                self.isr_pt.isr_hists[index].tunfolder.extract_truth_from_2d_hist_using_bin_definition(
                    raw_response_matrix.raw_root_hist))

            unfolded_hist = self.isr_pt.isr_hists[index].unfolded_measurement_hist

            acceptance_corr = Acceptance(mc_hist_full_phase, mc_acceptance_hist)
            acceptance_corrected = acceptance_corr.do_correction(unfolded_hist)

            unfolded_bin = self.isr_pt.unfolded_tunfold_bin
            self.isr_pt.set_acceptance_corrected_hist(
                hist=ISR2DHist(acceptance_corrected, unfolded_bin), mass_window_index=index)
            self.isr_pt.set_acceptance_corrected_hist(
                hist=ISR2DHist(mc_hist_full_phase, unfolded_bin), mass_window_index=index, key='simulation')

        if self.isr_pt.is_2d:
            run_acceptance_correction(index=0, is2d=True)
        else:
            for i in range(len(self.mass_bins)):
                run_acceptance_correction(index=i, is2d=False)

        # This is common in both 1D and 2D cases
        self.isr_pt.binned_mean_correction_factors = self.get_correction_factors(
            self.unfolded_space_name, self.pt_unfolded_bin_name
        )

    def mass_isr_acceptance_correction(self):
        postfix = '_' + str(self.pt_bins[0]) + 'to' + str(self.pt_bins[1])  # FIXME get this info from self.isr_mass
        hist_full_phase_name = self.mass_hist_full_phase_name_prefix + postfix
        mc_hist_full_phase = self.get_acceptance_hist(hist_full_phase_name)

        unfolded_hist = self.isr_mass.isr_hists[0].unfolded_measurement_hist
        mc_acceptance_hist = self.isr_mass.isr_hists[0].truth_signal_hist

        acceptance_corr = Acceptance(mc_hist_full_phase, mc_acceptance_hist)
        acceptance_corrected = acceptance_corr.do_correction(unfolded_hist)

        unfolded_bin = self.isr_mass.unfolded_tunfold_bin
        self.isr_mass.set_acceptance_corrected_hist(hist=ISR2DHist(acceptance_corrected, unfolded_bin))
        self.isr_mass.set_acceptance_corrected_hist(hist=ISR2DHist(mc_hist_full_phase, unfolded_bin), key='simulation')

    def get_correction_factors(self, space_name, bin_name):
        pt_hist_full_phase_name_prefix = ('dipt_gen_' + space_name +
                                          '_acceptance_dipt_0.0to100.0_dimass')
        correction_factors = []
        for index, mass_bin in enumerate(self.mass_bins):
            mass_bin_postfix = '_' + str(mass_bin[0]) + 'to' + str(mass_bin[1])
            mc_hist_full_phase = self.get_mc_hist(self.acceptance_name, pt_hist_full_phase_name_prefix + mass_bin_postfix)
            correction_factor = mc_hist_full_phase.get_mean(binned_mean=False)[0]
            correction_factors.append(correction_factor)

        return correction_factors

    def get_sim_prefsr_mean_pt_1d(self, unfolded_space_name='', unfolded_bin_name='', mass_bins=None):
        if unfolded_space_name == '':
            unfolded_space_name = self.unfolded_space_name  # dressed,
        if unfolded_bin_name == '':
            unfolded_bin_name = self.pt_unfolded_bin_name
        pt_hist_full_phase_name_prefix = ('dipt_[gen_' + unfolded_space_name +
                                          '_acceptance__' + unfolded_bin_name + ']_dimass')
        pt_data = []
        for mass_bin in self.mass_bins:
            mass_bin_postfix = '_' + str(mass_bin[0]) + 'to' + str(mass_bin[1])

            mc_hist_full_phase = self.get_mc_hist(self.signal_name, pt_hist_full_phase_name_prefix + mass_bin_postfix)
            pt_data.append(mc_hist_full_phase.get_mean(binned_mean=False))
        return pt_data

    def get_sim_prefsr_mean_pt_2d(self, unfolded_space_name='', unfolded_bin_name='',
                                  correct_binned_mean=False):
        if unfolded_space_name == '':
            unfolded_space_name = self.unfolded_space_name
        if unfolded_bin_name == '':
            unfolded_bin_name = self.pt_unfolded_bin_name
        pt_mass_unfolded_bin_name = "[tunfold-bin]_[dipt-dimass]_[" + unfolded_bin_name + "_O-window_v1_UO]"
        pt_mass_hist_full_phase_name = ("[tunfold-hist]_[dipt-dimass]_[gen_" + unfolded_space_name +
                                        "_acceptance__" + unfolded_bin_name + "_O-window_v1_UO]")
        mc_hist_full_phase = self.get_signal_hist(self.signal_name, pt_mass_hist_full_phase_name)

        unfolded_bin, _ = self.get_unfold_bin_maps(pt_mass_unfolded_bin_name, self.pt_mass_detector_bin_name)

        pt_data = self.extract_mean_pt_from_2d_hist(mc_hist_full_phase, unfolded_bin)
        if correct_binned_mean:
            self.correction_to_unbinned_prefsr_mean(pt_data, unfolded_space_name, unfolded_bin_name)
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

    def get_acceptance_hist(self, hist_name, hist_path='', bin_width_norm=False):
        return self.get_mc_hist(self.acceptance_name, hist_name, bin_width_norm=bin_width_norm)
