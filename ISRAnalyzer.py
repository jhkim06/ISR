from Analyzer import Analyzer
from Acceptance import Acceptance
from Hist import Hist, change_to_greek
from ISRHists import ISRHists
from HistTUnfoldBin import HistTUnfoldBin
from TUnFolder import TUnFolder


class ISRAnalyzer(Analyzer):
    def __init__(self,
                 sample_base_dir,
                 mass_bins,
                 pt_bins,
                 signal = "DY",
                 sys_on = True,
                 acceptance = None,  # DY for acceptance correction

                 pt_folded_bin_name='fine_O', pt_unfolded_bin_name='coarse_O',
                 folded_space_name='reco',
                 unfolded_space_name='gen_dressed', acceptance_space_name='dressed',

                 mass_folded_bin_name='fine_O', mass_unfolded_bin_name='coarse_O',
                 mass_folded_space_name='reco',
                 mass_unfolded_space_name='gen_dressed', mass_acceptance_space_name='dressed'):

        super(ISRAnalyzer, self).__init__(sample_base_dir,
                                          signal=signal, sys_on=sys_on,)

        if acceptance is None:
            self.acceptance_name = self.signal_name
        else:
            self.acceptance_name = acceptance
        self.mass_bins = mass_bins
        self.pt_bins = pt_bins

        self.dipt_label = r'$p_{T}$^{'+change_to_greek(self.channel)+'}'
        self.dimass_label = r'$M$^{'+change_to_greek(self.channel)+'}'

        # bin names for unfolding
        self.pt_folded_bin_name = pt_folded_bin_name
        self.pt_unfolded_bin_name = pt_unfolded_bin_name
        self.folded_space_name = folded_space_name
        self.unfolded_space_name = unfolded_space_name
        self.acceptance_space_name = acceptance_space_name

        self.mass_folded_bin_name = mass_folded_bin_name
        self.mass_unfolded_bin_name = mass_unfolded_bin_name
        self.mass_folded_space_name = mass_folded_space_name
        self.mass_unfolded_space_name = mass_unfolded_space_name
        self.mass_acceptance_space_name = mass_acceptance_space_name

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

        # dimass
        self.mass_hist_name_prefix = 'dimass_['+self.mass_folded_space_name+'__' + self.mass_folded_bin_name + ']_dipt'
        self.mass_fake_hist_name_prefix = ('dimass_['+self.mass_folded_space_name+'_' + self.unfolded_space_name
                                           + '_fake__' + self.mass_folded_bin_name + ']_dipt')
        self.mass_hist_full_phase_name_prefix = ('dimass_[gen_' + self.mass_acceptance_space_name +
                                                 '_acceptance__' + self.mass_unfolded_bin_name + ']_dipt')
        self.mass_hist_efficiency_name_prefix = ('dimass_[gen_' + self.mass_acceptance_space_name +
                                                 '_efficiency__' + self.mass_unfolded_bin_name + ']_dipt')
        self.mass_matrix_name_prefix = ('dimass_['+self.mass_folded_space_name+'__' + self.mass_folded_bin_name + ']_[' +
                                        self.mass_unfolded_space_name + '__' + self.mass_unfolded_bin_name + ']_dipt')

        self.mass_detector_bin_name = '[tunfold-bin]_dimass_['+self.mass_folded_bin_name+']'
        self.mass_truth_bin_name = '[tunfold-bin]_dimass_[' + self.mass_unfolded_bin_name + ']'

    def _set_isr_2d_hist_names(self):

        type_name = 'bin'
        var_name = 'dipt-dimass'
        tunfold_prefix = f"[tunfold-{type_name}]_[{var_name}]"
        bin_definition = '-window_v1_UO'
        # 2D case: set pt and mass histogram name
        self.pt_mass_detector_bin_postfix = "[" + self.folded_space_name + "__" + self.pt_folded_bin_name + bin_definition + "]"
        self.pt_mass_unfolded_bin_postfix = (
                "[" + self.unfolded_space_name + "__" + self.pt_unfolded_bin_name + bin_definition + "]")
        self.pt_mass_detector_bin_name = tunfold_prefix + "_[" + self.pt_folded_bin_name + bin_definition + "]"
        self.pt_mass_truth_bin_name = tunfold_prefix + "_[" + self.pt_unfolded_bin_name + bin_definition + "]"

        # even if the same bin definition is used, different contents can be stored
        type_name = 'hist'
        tunfold_prefix = f"[tunfold-{type_name}]_[{var_name}]"
        self.pt_mass_hist_name = tunfold_prefix + "_" + self.pt_mass_detector_bin_postfix
        self.pt_mass_fake_hist_name = (tunfold_prefix + "_["+self.folded_space_name+"_" + self.unfolded_space_name
                                       + "_fake__" + self.pt_folded_bin_name + bin_definition + "]")
        self.pt_mass_hist_full_phase_name = (tunfold_prefix + "_[gen_" + self.acceptance_space_name
                                             + "_acceptance__" + self.pt_unfolded_bin_name + bin_definition + "]")
        self.pt_mass_hist_efficiency_name = (tunfold_prefix + "_[gen_" + self.acceptance_space_name
                                             + "_efficiency__" + self.pt_unfolded_bin_name + bin_definition + "]")

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

    # setup using other results
    def setup_isr_hists(self, year, channel, event_selection, isr_pt, isr_mass):
        # dipt
        self.set_data_info(year, channel, event_selection)
        unfolded_bin, folded_bin = self.get_unfold_bin_maps(
            self.pt_mass_truth_bin_name,
            self.pt_mass_detector_bin_name
        )
        self.isr_pt = ISRHists(self.mass_bins, self.pt_bins, is_2d=True, is_pt=True,
                               year=year, channel=channel,
                               folded_tunfold_bin=folded_bin, unfolded_tunfold_bin=unfolded_bin)

        measurement = isr_pt.isr_hists[0].acceptance_corrected_hist['measurement']  # previous unfolded result
        signal = self.get_mc_hist(self.signal_name, self.pt_mass_hist_name)
        matrix = self.get_mc_hist(self.signal_name, self.pt_mass_matrix_name)
        self.isr_pt.set_isr_hist(measurement, signal, matrix)

        # dimass
        unfolded_bin, folded_bin = self.get_unfold_bin_maps(
            self.mass_truth_bin_name,
            self.mass_detector_bin_name
        )
        # for mass, always 1D unfolding
        self.isr_mass = ISRHists(self.mass_bins, self.pt_bins, is_2d=False, is_pt=False,
                                 year=year, channel=channel,
                                 folded_tunfold_bin=folded_bin, unfolded_tunfold_bin=unfolded_bin)
        measurement = isr_mass.isr_hists[0].acceptance_corrected_hist['measurement']
        postfix = '_' + str(self.pt_bins[0]) + 'to' + str(self.pt_bins[1])
        matrix_name = self.mass_matrix_name_prefix + postfix
        matrix = self.get_mc_hist(self.signal_name, matrix_name)
        signal = self.get_mc_hist(self.signal_name, self.mass_hist_name_prefix + postfix)
        self.isr_mass.set_isr_hist(measurement, signal, matrix)

    def setup_isr_acceptance_hists(self, year, channel, event_selection, is_2d=True, mass_bins=None,
                                   binned_mean=True, range_min=None, range_max=None,):
        self.set_data_info(year, channel, event_selection, only_theory_sys=True)
        unfolded_bin, folded_bin = self.get_unfold_bin_maps(
            self.pt_mass_truth_bin_name,
            self.pt_mass_detector_bin_name
        )

        if mass_bins is None:
            mass_bins = self.mass_bins

        # pt
        if is_2d:
            self.isr_pt = ISRHists(mass_bins, self.pt_bins, is_2d=True, is_pt=True,
                                   year=year, channel=channel,
                                   folded_tunfold_bin=folded_bin, unfolded_tunfold_bin=unfolded_bin)
            mc_hist_full_phase = self.get_acceptance_hist(self.pt_mass_hist_full_phase_name)
            mc_hist_efficiency = self.get_acceptance_hist(self.pt_mass_hist_efficiency_name)
            # add efficiency hist
            self.isr_pt.set_isr_hists(acceptance_corrected_signal_hist=mc_hist_full_phase,
                                      efficiency_signal_hist=mc_hist_efficiency,)
            self.isr_pt.set_ISRHistSet_per_mass_window()
            self.isr_pt.set_acceptance_corrected_mean_values_from_isr_hists_per_window(key='simulation')

            self.isr_pt.binned_mean_correction_factors = self.get_correction_factors(self.acceptance_space_name)
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
        hist_efficiency_name = self.mass_hist_efficiency_name_prefix + postfix
        mc_hist_full_phase = self.get_acceptance_hist(hist_full_phase_name)
        mc_efficiency_phase = self.get_acceptance_hist(hist_efficiency_name)

        unfolded_bin, folded_bin = self.get_unfold_bin_maps(
            self.mass_truth_bin_name,
            self.mass_detector_bin_name
        )
        self.isr_mass = ISRHists(mass_bins, self.pt_bins, is_2d=False, is_pt=False,
                                 year=year, channel=channel,
                                 folded_tunfold_bin=folded_bin, unfolded_tunfold_bin=unfolded_bin)
        self.isr_mass.set_isr_hists(acceptance_corrected_signal_hist=mc_hist_full_phase,
                                    efficiency_signal_hist=mc_efficiency_phase,)
        self.isr_mass.set_ISRHistSet_per_mass_window()
        self.isr_mass.set_acceptance_corrected_mean_values_from_isr_hists_per_window(key='simulation')
        self.isr_mass.binned_mean_correction_factors = self.get_correction_factors(self.acceptance_space_name,
                                                                                   is_pt=False, force_sys_off=True)

    def get_isr_hist(self, is_pt=True, is_2d=True, mass_window_index=0,
                     hist_prefix=''):
        if is_pt:
            if is_2d:
                hist_name = self.pt_mass_hist_name
                fake_hist_name = self.pt_mass_fake_hist_name
                matrix_name = self.pt_mass_matrix_name
            else:
                mass_bin_postfix = ('_' + str(self.mass_bins[mass_window_index][0]) +
                                    'to' + str(self.mass_bins[mass_window_index][1]))
                hist_name = self.pt_hist_name_prefix+mass_bin_postfix
                fake_hist_name = self.pt_fake_hist_name_prefix+mass_bin_postfix
                matrix_name = self.pt_matrix_name_prefix+mass_bin_postfix
        else:
            postfix = '_' + str(self.pt_bins[0]) + 'to' + str(self.pt_bins[1])
            hist_name, matrix_name, fake_hist_name, _ = self.get_hist_names_for_1d_dimass(postfix)

        measurement = self.get_measurement_hist(hist_name, hist_name_prefix=hist_prefix)
        signal = self.get_mc_hist(self.signal_name, hist_name, hist_name_prefix=hist_prefix)
        bgs = self.get_background_hists(hist_name, hist_name_prefix=hist_prefix)

        signal_fake = self.get_mc_hist(self.signal_name, fake_hist_name,)
        matrix = self.get_mc_hist(self.signal_name, matrix_name,)  #
        return measurement, signal, signal_fake, bgs, matrix

    def isr_unfolds(self):
        self.pt_isr_unfold()
        self.mass_isr_unfold()

    def _run_isr_unfold(
            self,
            isr_obj,
            mass_window_index: int = 0,
            tau: float = 0.0,
            max_iter: int = 4,
            reg_mode: str = 'None',
            tau_scan_method=None,
            do_iterative: bool = False,
            apply_custom_regularization: bool = False,
            custom_regularization_fn=None,
            do_iterative_for_last_window: bool = False,
            return_input_sys = False
    ) -> tuple:
        """
        Generic helper for ISR unfolding (used for both pt and mass, 1D or 2D).
        """
        # Prepare input
        input_hist, signal_hist, signal_fake_hist, background_hists, matrix = isr_obj.get_isr_hists(
            mass_window_index=mass_window_index)
        unfolded_bin = isr_obj.unfolded_tunfold_bin
        folded_bin = isr_obj.folded_tunfold_bin

        # Optionally turn on iterative for last window
        if do_iterative_for_last_window and hasattr(isr_obj, "is_2d") and not isr_obj.is_2d:
            if mass_window_index == len(self.mass_bins) - 1:
                do_iterative = True

        # Setup TUnFolder for measurement
        unfold = TUnFolder(
            matrix,
            input_hist,
            signal_fake_hist,
            bg_hists=background_hists,
            unfolded_bin=unfolded_bin,
            folded_bin=folded_bin,
            variable_name='',
            iterative=do_iterative,
            reg_mode=reg_mode,
            tau_scan_method=tau_scan_method,
        )

        # Apply regularization if needed
        if (not do_iterative and tau_scan_method is not None and custom_regularization_fn is not None):
            getattr(unfold, custom_regularization_fn)()
        if apply_custom_regularization and custom_regularization_fn is not None:
            getattr(unfold, custom_regularization_fn)()

        tau, max_iter = unfold.unfold(tau=tau, max_iter=max_iter)

        # Closure test (truth MC input, same tau)
        closure = TUnFolder(
            matrix,
            signal_hist,
            signal_fake_hist,
            unfolded_bin=unfolded_bin,
            folded_bin=folded_bin,
            variable_name='',
            iterative=do_iterative,
            reg_mode=reg_mode,
            tau_scan_method=None,
        )
        closure.unfold(tau=tau, max_iter=max_iter)

        # Create output hists
        unfold_input_hist = input_hist.create(hist=unfold.get_input_hist(), label=input_hist.label + '(unfold input)')
        unfolded_hist = input_hist.create(hist=unfold.get_unfolded_hist(), label=input_hist.label + '(unfolded)')
        # no reason to repeat unfolding for scale alpha_s and pdf since,
        # modelling uncertainty for unfold done with reweighted DY MC
        sys_names_to_skip = ['scale', 'alpha_s', 'pdf']
        unfolded_hist.systematic_raw_root_hists, sys_input_hist = unfold.sys_unfold(sys_names_to_skip=sys_names_to_skip,
                                                                                    return_input_sys=return_input_sys)
        if sys_input_hist is not None:
            unfold_input_hist.systematic_raw_root_hists=sys_input_hist
            unfold_input_hist.compute_systematic_rss_per_sysname()
        unfolded_hist.systematic_raw_root_hists.update({"matrix_stat": unfold.get_matrix_stat()})
        unfolded_hist.compute_systematic_rss_per_sysname()

        truth_signal_hist = unfold.get_mc_truth_from_response_matrix(sys_on=True)
        raw_hist = unfold.get_mc_reco_from_response_matrix()
        reco_signal_hist = signal_hist.create(hist=raw_hist, hist_name='_', label=signal_hist.label)
        unfolded_signal_hist = input_hist.create(hist=closure.get_unfolded_hist(), hist_name="_", label='Unfolded DY')

        # Store results in isr_obj
        isr_obj.isr_hists[mass_window_index].tunfolder = unfold
        isr_obj.isr_hists[mass_window_index].unfold_input_hist = HistTUnfoldBin(unfold_input_hist, folded_bin)
        isr_obj.isr_hists[mass_window_index].unfolded_measurement_hist = HistTUnfoldBin(unfolded_hist, unfolded_bin)
        isr_obj.isr_hists[mass_window_index].truth_signal_hist = HistTUnfoldBin(truth_signal_hist, unfolded_bin)
        isr_obj.isr_hists[mass_window_index].reco_signal_hist = HistTUnfoldBin(reco_signal_hist, folded_bin)
        isr_obj.isr_hists[mass_window_index].unfolded_signal_hist = HistTUnfoldBin(unfolded_signal_hist, unfolded_bin)
        return tau, max_iter

    def pt_isr_unfold(
            self, tau=0.0, max_iter=4, reg_mode='None', tau_scan_method=None,
            do_iterative=False, do_iterative_for_last_window=False,
            apply_custom_regularization_for_pt=False,
            return_input_sys=False
    ):
        """
        Unfold the pt distribution, handling 2D and 1D cases.
        """
        custom_reg_fn = "apply_custom_regularization_for_pt"
        if getattr(self.isr_pt, "is_2d", False):
            return self._run_isr_unfold(
                self.isr_pt, 0, tau, max_iter, reg_mode, tau_scan_method,
                do_iterative, apply_custom_regularization_for_pt,
                custom_regularization_fn=custom_reg_fn,
                return_input_sys=return_input_sys,
            )
        else:
            # Unlikely used for this time...
            for idx, _ in enumerate(self.mass_bins):
                self._run_isr_unfold(
                    self.isr_pt, idx, tau, max_iter, reg_mode, tau_scan_method,
                    do_iterative, apply_custom_regularization_for_pt,
                    custom_regularization_fn=custom_reg_fn,
                    do_iterative_for_last_window=do_iterative_for_last_window
                )

    def mass_isr_unfold(
            self, tau=0.0, max_iter=4, reg_mode='None', tau_scan_method=None,
            do_iterative=False, apply_custom_regularization_for_mass=False,
            return_input_sys=False
    ):
        """
        Unfold the mass distribution (always 1D).
        """
        custom_reg_fn = "apply_custom_regularization_for_mass"
        return self._run_isr_unfold(
            self.isr_mass, 0, tau, max_iter, reg_mode, tau_scan_method,
            do_iterative, apply_custom_regularization_for_mass,
            custom_regularization_fn=custom_reg_fn,
            return_input_sys=return_input_sys,
        )

    def isr_acceptance_corrections(self):
        self.pt_isr_acceptance_correction()
        self.mass_isr_acceptance_correction()

    def _run_acceptance_correction(
            self,
            isr_obj,
            mass_bins,
            acceptance_space_name: str,
            is_pt: bool,
            is_2d: bool,
            pt_bins=None,
            get_acceptance_hist=None,
            get_correction_factors=None,
            hist_names=None
    ):
        """
        Generic helper for acceptance correction (handles both pt and mass, 1D and 2D).
        """
        # Set default prefix accessors if not provided
        hist_names = hist_names or {}

        # (is_pt and is_2d) or (not is_pt)
        #
        n_windows = 1 if (is_pt and is_2d) or (not is_pt) else len(mass_bins)
        for index in range(n_windows):
            if (is_2d and is_pt) or (not is_pt):
                # 2D pt: all mass windows together
                mc_hist_full_phase = get_acceptance_hist(hist_names['full_phase'])
                mc_hist_efficiency = get_acceptance_hist(hist_names['efficiency'])
                matrix_name = hist_names['matrix']
            else:
                # This is 1D pt case, will not be used for now
                mass_bin = mass_bins[index]
                mass_bin_postfix = f'_{mass_bin[0]}to{mass_bin[1]}'
                mc_hist_full_phase = get_acceptance_hist(hist_names['full_phase'] + mass_bin_postfix)
                mc_hist_efficiency = None  # Can add if needed for 1D pt
                matrix_name = hist_names['matrix'] + mass_bin_postfix

            mc_acceptance_hist = get_acceptance_hist(matrix_name, is_matrix=True)
            unfolded_hist = isr_obj.isr_hists[index].unfolded_measurement_hist

            acceptance_corr = Acceptance(mc_hist_full_phase, mc_acceptance_hist)
            acceptance_corrected = acceptance_corr.do_correction(unfolded_hist)
            acceptance_corrected.systematic_raw_root_hists.update(
                {"accept_stat": acceptance_corr.get_accept_stat(unfolded_hist)}
            )
            acceptance_corrected.compute_systematic_rss_per_sysname()

            unfolded_bin = isr_obj.unfolded_tunfold_bin
            isr_obj.set_acceptance_corrected_hist(
                hist=HistTUnfoldBin(acceptance_corrected, unfolded_bin), mass_window_index=index
            )
            isr_obj.set_acceptance_corrected_hist(
                hist=HistTUnfoldBin(mc_hist_full_phase, unfolded_bin), mass_window_index=index, key='simulation'
            )
            # Only set efficiency hist if it exists (i.e., for 2D pt or for mass with efficiency available)
            if mc_hist_efficiency is not None:
                isr_obj.isr_hists[index].efficiency_signal_hist = HistTUnfoldBin(mc_hist_efficiency, unfolded_bin)

        # Set binned mean correction factors
        isr_obj.binned_mean_correction_factors = get_correction_factors(
            acceptance_space_name, is_pt=is_pt
        )

    def pt_isr_acceptance_correction(self):
        """
        Apply acceptance correction for pt (handles both 1D and 2D).
        """
        # Prepare prefixes
        hist_names = {
            'full_phase': self.pt_hist_full_phase_name_prefix if not self.isr_pt.is_2d else self.pt_mass_hist_full_phase_name,
            'efficiency': self.pt_mass_hist_efficiency_name if self.isr_pt.is_2d else None,
            'matrix': self.pt_matrix_name_prefix if not self.isr_pt.is_2d else self.pt_mass_matrix_name
        }
        self._run_acceptance_correction(
            isr_obj=self.isr_pt,
            mass_bins=self.mass_bins,
            acceptance_space_name=self.acceptance_space_name,
            is_pt=True,
            is_2d=self.isr_pt.is_2d,
            pt_bins=self.pt_bins,
            get_acceptance_hist=self.get_acceptance_hist,
            get_correction_factors=self.get_correction_factors,
            hist_names=hist_names
        )

    def mass_isr_acceptance_correction(self):
        """
        Apply acceptance correction for mass (always 1D, single mass window).
        """
        postfix = '_' + str(self.pt_bins[0]) + 'to' + str(self.pt_bins[1])
        hist_names = {
            'full_phase': self.mass_hist_full_phase_name_prefix + postfix,
            'efficiency': self.mass_hist_efficiency_name_prefix + postfix,
            'matrix': self.mass_matrix_name_prefix + postfix
        }
        self._run_acceptance_correction(
            isr_obj=self.isr_mass,
            mass_bins=[self.mass_bins[0]],  # mass is always a single window in this context
            acceptance_space_name=self.acceptance_space_name,
            is_pt=False,
            is_2d=False,
            pt_bins=self.pt_bins,
            get_acceptance_hist=self.get_acceptance_hist,
            get_correction_factors=lambda *args, **kwargs: self.get_correction_factors(
                self.acceptance_space_name, is_pt=False, force_sys_off=True
            ),
            hist_names=hist_names
        )

    # get binned from isr_pt.per_mass_ or isr_mass
    # check variations of correction values according to systematics
    def get_binned_mean_correction_values(self, sys_name='default', mass_index=0):
        unbinned_full_phase_name_prefix = ('dipt_gen_' + self.acceptance_space_name +
                                           '_acceptance_dipt_0.0to100.0_dimass')
        mass_bin_postfix = '_' + str(self.mass_bins[mass_index][0]) + 'to' + str(self.mass_bins[mass_index][1])
        if sys_name == 'default':
            mc_hist_full_phase = self.get_mc_hist(self.acceptance_name,
                                                  unbinned_full_phase_name_prefix + mass_bin_postfix,
                                                  force_sys_off=True)  # TODO test systematics at generator level
            unbinned_mean = mc_hist_full_phase.get_mean(binned_mean=False)[0]
            binned_mean = self.isr_pt.isr_hists_per_mass_window[mass_index].acceptance_corrected_hist['simulation'].get_mean(binned_mean=False)[0]
            # print('correction factor', binned_mean/unbinned_mean)
            return unbinned_mean, binned_mean
        else:
            mc_hist_full_phase = self.get_mc_hist(self.acceptance_name,
                                                  unbinned_full_phase_name_prefix + mass_bin_postfix,
                                                  force_sys_off=False)  # TODO test systematics at generator level
            unbinned_mean = mc_hist_full_phase.get_sys_mean_dfs(sys_name, binned_mean=False)
            binned_mean = \
            self.isr_pt.isr_hists_per_mass_window[mass_index].acceptance_corrected_hist['simulation'].get_sys_mean_dfs(sys_name,
                                                                                                                       binned_mean=False)
            #print('binned_mean', binned_mean)
            #print('unbinned_mean', unbinned_mean)
            return unbinned_mean, binned_mean

    def get_correction_factors(self, space_name, is_pt=True, force_sys_off=False):
        # dipt
        if is_pt:
            unbinned_full_phase_name_prefix = ('dipt_gen_' + space_name +
                                              '_acceptance_dipt_0.0to100.0_dimass')
        else:
            unbinned_full_phase_name_prefix = ('dimass_gen_' + space_name +
                                               '_acceptance_dipt_0.0to100.0_dimass')
        correction_factors = []
        for index, mass_bin in enumerate(self.mass_bins):
            mass_bin_postfix = '_' + str(mass_bin[0]) + 'to' + str(mass_bin[1])
            # FIXME
            mc_hist_full_phase = self.get_mc_hist(self.acceptance_name,
                                                  unbinned_full_phase_name_prefix + mass_bin_postfix,
                                                  force_sys_off=force_sys_off)
            correction_factor = mc_hist_full_phase.get_mean(binned_mean=False)[0]
            correction_factors.append(correction_factor)
        # correction_factor["nominal"]["nominal"]
        return correction_factors

    def get_acceptance_hist(self, hist_name, hist_path='', bin_width_norm=False,
                            is_matrix=False):
        if is_matrix:
            matrix = self.get_mc_hist(self.acceptance_name, hist_name, sys_names_to_skip=['matrix_model'])
            # get projected hists
            truth_hist = Hist(matrix.get_raw_hist().ProjectionX("projected_truth",  0, -1, "e"),
                              hist_name=hist_name + "_projected_truth",
                              year=self.year, channel=self.channel, label=matrix.label)

            for sys_name, variations in matrix.systematic_raw_root_hists.items():
                for var_name, hist in variations.items():
                    temp_hist = hist.ProjectionX("projected_truth"+sys_name+var_name,  0, -1, "e")
                    truth_hist.set_systematic_hist(sys_name, var_name, temp_hist)
            truth_hist.compute_systematic_rss_per_sysname()
            return truth_hist
        else:
            return self.get_mc_hist(self.acceptance_name, hist_name, bin_width_norm=bin_width_norm,
                                    force_sys_off=False, sys_names_to_skip=['matrix_model'])

    def extract_mean_pt_from_2d_hist(self, pt_2d_hist, unfold_bin):
        pt_data = []
        for index, _ in enumerate(self.mass_bins):
            axis_steering = 'dipt[O];dimass[UOC' + str(index) + ']'
            temp_result = Hist(unfold_bin.ExtractHistogram("", pt_2d_hist.get_raw_hist(), 0,
                                                           True,
                                                           axis_steering))
            pt_data.append(temp_result.get_mean())
        return pt_data

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
