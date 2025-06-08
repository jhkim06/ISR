from traceback import print_tb

from Hist import Hist, change_to_greek
from Analyzer import labels, colors, get_hist_kwargs
from HistTUnfoldBin import HistTUnfoldBin
import pandas as pd
import numpy as np
from ISRLinearFitter import ISRLinearFitter


def normalize(hist, bin_width_norm=True, scale=1.0):
    return hist.bin_width_norm(scale) if bin_width_norm else hist


# group of hists for ISR
class ISRHistSet:
    def __init__(self,
                 measurement_hist=None,
                 signal_hist=None,
                 signal_fake_hist=None,
                 background_hists=None,
                 response_matrix=None,
                 truth_signal_hist=None,
                 acceptance_corrected_signal_hist=None):
        # detector hist
        self.measurement_hist = measurement_hist
        self.signal_hist = signal_hist
        self.signal_fake_hist = signal_fake_hist
        self.background_hist = background_hists  # Note this is a dictionary of hists

        self.response_matrix = response_matrix

        # unfolded hist
        self.tunfolder = None
        self.unfold_input_hist = None  #
        self.reco_signal_hist = None
        self.unfolded_measurement_hist = None
        self.unfolded_signal_hist = None  # simple closure
        self.truth_signal_hist = truth_signal_hist  # request to tunfolder

        # acceptance corrected hist
        self.acceptance_corrected_hist = {"measurement": None,
                                          "simulation": acceptance_corrected_signal_hist}

    def get_extracted_1d_hist_set(self, index):
        measurement_hist_extracted = self.measurement_hist.extract_1d_hist(index)
        signal_hist_extracted = self.signal_hist.extract_1d_hist(index)
        if self.signal_fake_hist is not None:
            signal_fake_hist_extracted = self.signal_fake_hist.extract_1d_hist(index)
        else:
            signal_fake_hist_extracted = None
        background_hist_extracted = {}
        if self.background_hist is not None:
            for key, value in self.background_hist.items():
                background_hist_extracted[key] = self.background_hist[key].extract_1d_hist(index)
        else:
            background_hist_extracted = None

        unfold_input_hist_extracted = self.unfold_input_hist.extract_1d_hist(index) if self.unfold_input_hist else None
        reco_signal_hist_extracted = self.reco_signal_hist.extract_1d_hist(index) if self.reco_signal_hist else None
        unfolded_measurement_hist_extracted = self.unfolded_measurement_hist.extract_1d_hist(index) if self.unfolded_measurement_hist else None
        unfolded_signal_hist_extracted = self.unfolded_signal_hist.extract_1d_hist(index) if self.unfolded_signal_hist else None
        truth_signal_hist_extracted = self.truth_signal_hist.extract_1d_hist(index) if self.truth_signal_hist else None

        acceptance_corrected_hist_extracted = {
            "measurement": self.acceptance_corrected_hist["measurement"].extract_1d_hist(index)
            if self.acceptance_corrected_hist["measurement"] else None,
            "simulation": self.acceptance_corrected_hist["simulation"].extract_1d_hist(index)
            if self.acceptance_corrected_hist["simulation"] else None,
        }

        extracted_set = ISRHistSet()
        extracted_set.measurement_hist = measurement_hist_extracted
        extracted_set.signal_hist = signal_hist_extracted
        extracted_set.signal_fake_hist = signal_fake_hist_extracted
        extracted_set.background_hist = background_hist_extracted
        extracted_set.unfold_input_hist = unfold_input_hist_extracted
        extracted_set.reco_signal_hist = reco_signal_hist_extracted
        extracted_set.unfolded_measurement_hist = unfolded_measurement_hist_extracted
        extracted_set.unfolded_signal_hist = unfolded_signal_hist_extracted
        extracted_set.response_matrix = None
        extracted_set.truth_signal_hist = truth_signal_hist_extracted
        extracted_set.acceptance_corrected_hist = acceptance_corrected_hist_extracted

        return extracted_set


    def get_extracted_hist_set(self,):
        measurement_hist_extracted = self.measurement_hist.extract_hist()
        signal_hist_extracted = self.signal_hist.extract_hist()
        if self.signal_fake_hist is not None:
            signal_fake_hist_extracted = self.signal_fake_hist.extract_hist()
        else:
            signal_fake_hist_extracted = None
        background_hist_extracted = {}
        if self.background_hist is not None:
            for key, value in self.background_hist.items():
                background_hist_extracted[key] = self.background_hist[key].extract_hist()
        else:
            background_hist_extracted = None

        unfold_input_hist_extracted = self.unfold_input_hist.extract_hist() if self.unfold_input_hist else None
        reco_signal_hist_extracted = self.reco_signal_hist.extract_hist() if self.reco_signal_hist else None
        unfolded_measurement_hist_extracted = self.unfolded_measurement_hist.extract_hist() if self.unfolded_measurement_hist else None
        unfolded_signal_hist_extracted = self.unfolded_signal_hist.extract_hist() if self.unfolded_signal_hist else None
        truth_signal_hist_extracted = self.truth_signal_hist.extract_hist() if self.truth_signal_hist else None

        acceptance_corrected_hist_extracted = {
            "measurement": self.acceptance_corrected_hist["measurement"].extract_hist() if self.acceptance_corrected_hist["measurement"] else None,
            "simulation": self.acceptance_corrected_hist["simulation"].extract_hist() if self.acceptance_corrected_hist["simulation"] else None,
        }

        extracted_set = ISRHistSet()
        extracted_set.measurement_hist = measurement_hist_extracted
        extracted_set.signal_hist = signal_hist_extracted
        extracted_set.signal_fake_hist = signal_fake_hist_extracted
        extracted_set.background_hist = background_hist_extracted
        extracted_set.unfold_input_hist = unfold_input_hist_extracted
        extracted_set.reco_signal_hist = reco_signal_hist_extracted
        extracted_set.unfolded_measurement_hist = unfolded_measurement_hist_extracted
        extracted_set.unfolded_signal_hist = unfolded_signal_hist_extracted
        extracted_set.response_matrix = None
        extracted_set.truth_signal_hist = truth_signal_hist_extracted
        extracted_set.acceptance_corrected_hist = acceptance_corrected_hist_extracted

        return extracted_set

# ISRPtHists ISRMassHists
class ISRHists:
    def __init__(self,
                 mass_bins,  # mass windows
                 pt_bins,

                 is_2d,
                 is_pt,

                 year='', channel='',
                 folded_tunfold_bin=None, unfolded_tunfold_bin=None):  # pt cut

        self.year = year
        self.channel = channel

        self.mass_bins = mass_bins
        self.pt_bins = pt_bins

        self.is_2d = is_2d
        '''
        what doest self.is_2d mean?
        1) histograms ue TUnfoldBinning 
        2) it is 2 dimensional  
        '''
        self.is_pt = is_pt

        # detector level hists
        self.isr_hists = []  # depend on input format
        self.isr_hists_per_mass_window = []

        self.use_tunfoldbinning = False
        self.folded_tunfold_bin = folded_tunfold_bin
        self.unfolded_tunfold_bin = unfolded_tunfold_bin
        if folded_tunfold_bin and unfolded_tunfold_bin:
            self.use_tunfoldbinning = True
        #
        self.measurement_mean_values = None
        self.measurement_bg_subtracted_mean_values = None
        self.signal_mean_values = None
        self.background_mean_values = None

        self.unfolded_measurement_mean_values = None
        self.unfolded_signal_mean_values = None

        # N mean values
        self.acceptance_corrected_mean_values = {"measurement": [],
                                                 "simulation": []}

        self.binned_mean_correction_factors = None
        # common plot cosmetics
        if is_pt:
            self.x_axis_label = r"$p_{T}^{"+change_to_greek(self.channel)+"}$ [GeV]"
        else:
            self.x_axis_label = r"$m^{"+change_to_greek(self.channel)+"}$ [GeV]"

    def set_ISRHistSet_per_mass_window(self):
        if self.is_pt:
            if self.is_2d:
                # extract 1D and then copy to self.isr_hists_per_mass_window
                for index, _ in enumerate(self.mass_bins):
                    self.isr_hists_per_mass_window.append(self.isr_hists[0].get_extracted_1d_hist_set(index))
            else:
                # just copy it to self.isr_hists_per_mass_window
                for isr_hist in self.isr_hists:
                    self.isr_hists_per_mass_window.append(isr_hist.get_extracted_hist_set())
        else:
            # just copy it to self.isr_hists_per_mass_window
            self.isr_hists_per_mass_window.append(self.isr_hists[0].get_extracted_hist_set())

    # from other.isr_hists_per_mass_window
    def add_external_hist_as_sys_hist(self, other, sys_name, var_name=None):
        if var_name is None:
            var_name = sys_name
        for index, isr_hist in enumerate(self.isr_hists_per_mass_window):
            other_isr_hist = other.isr_hists_per_mass_window[index]

            # since FSR use different unfolding phase space, do not include
            if sys_name != "FSR":
                hist = other_isr_hist.unfolded_measurement_hist.get_raw_hist()
                isr_hist.unfolded_measurement_hist.set_systematic_hist(sys_name, var_name, hist)

            hist = other_isr_hist.acceptance_corrected_hist["measurement"].get_raw_hist()
            isr_hist.acceptance_corrected_hist["measurement"].set_systematic_hist(sys_name, var_name, hist)

            if sys_name != "FSR":
                hist = other_isr_hist.acceptance_corrected_hist["simulation"].get_raw_hist()
                isr_hist.acceptance_corrected_hist["simulation"].set_systematic_hist(sys_name, var_name, hist)

    def update_systematics(self):
        for index, isr_hist in enumerate(self.isr_hists_per_mass_window):
            # separate
            isr_hist.unfolded_measurement_hist.compute_systematic_rss_per_sysname()
            isr_hist.acceptance_corrected_hist["measurement"].compute_systematic_rss_per_sysname()
            isr_hist.acceptance_corrected_hist["simulation"].compute_systematic_rss_per_sysname()

            # update mean value
            # set_acceptance_corrected_mean_values
            self.set_acceptance_corrected_mean_values_(mass_window_index=index, key='measurement')
            self.set_acceptance_corrected_mean_values_(mass_window_index=index, key='simulation')

    # TODO update hist also
    # add systematic from other ISRHist object
    def update_mean_values(self, other, sys_name):
        key = 'measurement'
        for index, mean_value in enumerate(self.acceptance_corrected_mean_values[key]):
            sys_mean = other.acceptance_corrected_mean_values[key][index]['mean']
            mean = mean_value['mean']
            insert_idx = mean_value.columns.get_loc('total_error')
            mean_value.insert(insert_idx, sys_name, abs(sys_mean-mean))

            error_columns = mean_value.columns.difference(['mean'])
            mean_value['total_error'] = np.sqrt((mean_value[error_columns] ** 2).sum(axis=1))

    def set_isr_hist(self, measurement, signal, matrix):
        self.isr_hists.append(ISRHistSet(measurement_hist=measurement,
                                         signal_hist=HistTUnfoldBin(signal, self.folded_tunfold_bin),
                                         response_matrix=matrix))

    # Note set *detector* level hists
    def set_isr_hists(self, measurement_hist=None,
                      signal_hist=None, signal_fake_hist=None,
                      background_hists=None, matrix=None,
                      acceptance_corrected_signal_hist=None):
        if self.is_pt and len(self.mass_bins) == len(self.isr_hists):
            print('Check the number of mass bins and the number of ISR hists...')
            return
        else:
            if self.use_tunfoldbinning:
                new_background_hists = {}
                if background_hists:
                    for key, value in background_hists.items():
                        new_background_hists[key] = HistTUnfoldBin(value, self.folded_tunfold_bin, )
                else:
                    new_background_hists = None


                self.isr_hists.append(
                    ISRHistSet(
                        measurement_hist=HistTUnfoldBin(measurement_hist, self.folded_tunfold_bin) if measurement_hist else None,
                        signal_hist=HistTUnfoldBin(signal_hist, self.folded_tunfold_bin) if signal_hist else None,
                        signal_fake_hist=HistTUnfoldBin(signal_fake_hist, self.folded_tunfold_bin) if signal_fake_hist else None,
                        background_hists=new_background_hists,
                        response_matrix=matrix if matrix else None,
                        acceptance_corrected_signal_hist=acceptance_corrected_signal_hist if acceptance_corrected_signal_hist else None,
                    )
                )
            else:
                self.isr_hists.append(
                    ISRHistSet(measurement_hist=measurement_hist,
                               signal_hist=signal_hist,
                               signal_fake_hist=signal_fake_hist, background_hists=background_hists,
                               response_matrix=matrix,
                               acceptance_corrected_signal_hist=acceptance_corrected_signal_hist,)
                )

    def set_acceptance_corrected_hist(self, hist,
                                      mass_window_index=0, key='measurement'):
        self.isr_hists[mass_window_index].acceptance_corrected_hist[key] = hist
        self.set_acceptance_corrected_mean_values(mass_window_index, key)

    # TODO set detector/unfolded mean values also
    def set_acceptance_corrected_mean_values_(self, mass_window_index, key='measurement',
                                              binned_mean=True, range_min=None, range_max=None,):
        self.acceptance_corrected_mean_values[key].clear()
        if self.is_pt:
            for index in range(len(self.mass_bins)):
                mean = self.isr_hists_per_mass_window[index].acceptance_corrected_hist[key].get_mean_df(
                    binned_mean=binned_mean,)
                self.acceptance_corrected_mean_values[key].append(mean)
        else:
            for low, high in self.mass_bins:
                mean = self.isr_hists_per_mass_window[0].acceptance_corrected_hist[key].get_mean_df(range_min=low,
                                                                                                    range_max=high)
                self.acceptance_corrected_mean_values[key].append(mean)

    def set_acceptance_corrected_mean_values(self, mass_window_index=0, key='measurement',
                                             binned_mean=True, range_min=None, range_max=None,):
        hist = self.isr_hists[mass_window_index].acceptance_corrected_hist[key]

        # is_pt
        if self.is_2d or self.is_pt==False:
            # loop over mass bins
            if self.is_pt:
                for index in range(len(self.mass_bins)):
                    extracted_hist = hist.extract_1d_hist(index=index)
                    mean = extracted_hist.get_mean_df(binned_mean=binned_mean)
                    self.acceptance_corrected_mean_values[key].append(mean)
            else:
                # mass
                extracted_hist = hist.extract_hist()
                for low, high in self.mass_bins:
                    mean = extracted_hist.get_mean_df(range_min=low, range_max=high)
                    self.acceptance_corrected_mean_values[key].append(mean)
        else:
            # set mean for mass_window_index
            # TODO 1D case self.is_pt == True, self.is_2d == False
            if self.is_pt:
                extracted_hist = hist.extract_hist()
                mean = extracted_hist.get_mean_df(binned_mean=binned_mean, range_min=range_min, range_max=range_max)
                self.acceptance_corrected_mean_values[key].append(mean)

    def get_isr_hists(self, mass_window_index=-1):
        if self.is_2d or self.is_pt==False:
            isr_hist = self.isr_hists[0]
        else:
            isr_hist = self.isr_hists[mass_window_index]
        return (isr_hist.measurement_hist,
                isr_hist.signal_hist,
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

    def set_1d_isr_hists(self):
        for index, _ in enumerate(self.mass_bins):
            pass

    def get_hist_in_mass_window(self, hist_type, index, bin_width_norm=True, scale=1,
                                key='measurement'):
        # TODO check if index is valid
        isr_hist_set = self.isr_hists_per_mass_window[index]

        # Fetch the target hist or dict of hists
        attr_name = hist_type + "_hist"
        target_hist = getattr(isr_hist_set, attr_name)
        if hist_type == 'acceptance_corrected':
            target_hist = target_hist[key]
        if target_hist is None:
            raise ValueError(f"No attribute named '{attr_name}' found")

        if hist_type == 'background':
            return {bg_name: normalize(bg_hist, bin_width_norm=bin_width_norm, scale=scale)
                    for bg_name, bg_hist in target_hist.items()}
        else:
            return normalize(target_hist, bin_width_norm=bin_width_norm, scale=scale)

    # get from self.isr_hists which could have different formats
    def get(self, hist_type='signal', mass_window_index=-1,
            bin_width_norm=False, scale=1, key='measurement'):
        # Choose the relevant ISR hist container
        if self.is_pt and not self.is_2d:
            isr_hist_to_draw = self.isr_hists[mass_window_index]
        else:
            isr_hist_to_draw = self.isr_hists[0]

        # Fetch the target hist or dict of hists
        attr_name = hist_type + "_hist"
        target_hist = getattr(isr_hist_to_draw, attr_name)
        if hist_type == 'acceptance_corrected':
            target_hist = target_hist[key]
        if target_hist is None:
            raise ValueError(f"No attribute named '{attr_name}' found")

        # Handle 2D case
        if self.is_pt:  #
            if self.is_2d:
                if 0 <= mass_window_index < len(self.mass_bins):  # -1 for 2D
                    if hist_type == 'background':
                        return {k: normalize(v.extract_1d_hist(mass_window_index), bin_width_norm=bin_width_norm)
                                for k, v in target_hist.items()}
                    return normalize(target_hist.extract_1d_hist(mass_window_index), bin_width_norm=bin_width_norm)
                return normalize(target_hist, bin_width_norm=bin_width_norm)
            else:
                if hist_type == 'background':
                    return {k: normalize(v.extract_hist(), bin_width_norm=bin_width_norm)
                            for k, v in target_hist.items()}
                return normalize(target_hist.extract_hist(), bin_width_norm=bin_width_norm)
        else:
            if hist_type == 'background':
                return {k: normalize(v.extract_hist(), bin_width_norm=bin_width_norm)
                        for k, v in target_hist.items()}
            return normalize(target_hist.extract_hist(), bin_width_norm=bin_width_norm)

    def get_additional_text_on_plot(self, mass_window_index=-1):
        text = ''
        if self.is_pt:
            if mass_window_index == -1:
                return ""
            else:
                text = (str(int(self.mass_bins[mass_window_index][0]))+
                        "$<m^{"+change_to_greek(self.channel)+"}<$"+
                        str(int(self.mass_bins[mass_window_index][1])) + " GeV")
        else:
            text = str(r"$p_{T}^{" + change_to_greek(self.channel) + "}<$" + str(int(self.pt_bins[1])) + " (GeV)")
        return text

    # TODO systematics on binned_mean_correction?
    def get_mean_df(self, key='measurement', other=None, binned_mean_correction=True, ):
        # Choose whether to pull from self or another instance
        source = other if other is not None else self
        # Build the base dataframe
        df = pd.concat(source.acceptance_corrected_mean_values[key], ignore_index=True)
        # If PT case, compute and apply correction
        if source.is_pt:
            bmc = source.binned_mean_correction_factors  # unbinned mean values
            df_sim = pd.concat(source.acceptance_corrected_mean_values['simulation'], ignore_index=True)
            bmc = bmc / df_sim['mean']
            if binned_mean_correction:
                df['mean'] *= bmc
        return df

    def get_sys_mean_dfs(self, pt, mass, sys_name=None, key='measurement', other=None, binned_mean_correction=True,):
        if sys_name is None:
            # get all systematic names
            sys_names = mass.isr_hists_per_mass_window[0].acceptance_corrected_hist[key].systematic_raw_root_hists.keys()
        else:
            sys_names = [sys_name]

        sys_mean_pt_dfs = {}
        sys_mean_mass_dfs = {}
        for sys_name in sys_names:
            # loop over variations and plot the means!
            sys_mean_pt = []  # first mass window sys, second mass window sys
            sys_mean_mass = []

            sys_mean_pt_dfs[sys_name] = {}
            sys_mean_mass_dfs[sys_name] = {}
            for index, bin_edge in enumerate(self.mass_bins):
                sys_mean_pt.append(pt.isr_hists_per_mass_window[index].acceptance_corrected_hist[key].get_sys_mean_dfs(
                    sys_name))
                sys_mean_mass.append(mass.isr_hists_per_mass_window[0].acceptance_corrected_hist[key].get_sys_mean_dfs(
                    sys_name, range_min=bin_edge[0], range_max=bin_edge[1]))

            for sys_index in range(len(sys_mean_pt[0])):
                sys_list_pt = []
                sys_list_mass = []
                for mass_window_index, bin_edge in enumerate(self.mass_bins):
                    sys_list_pt.append(sys_mean_pt[mass_window_index][sys_index])
                    sys_list_mass.append(sys_mean_mass[mass_window_index][sys_index])

                mass_df = pd.concat(sys_list_mass, ignore_index=True)
                pt_df = pd.concat(sys_list_pt, ignore_index=True)

                bmc = pt.binned_mean_correction_factors  # unbinned mean values
                df_sim = pd.concat(pt.acceptance_corrected_mean_values['simulation'], ignore_index=True)
                bmc = bmc / df_sim['mean']
                pt_df['mean'] *= bmc

                sys_mean_pt_dfs[sys_name][sys_index] = pt_df
                sys_mean_mass_dfs[sys_name][sys_index] = mass_df
        return sys_mean_pt_dfs, sys_mean_mass_dfs

    def add_isr_plot(self, plotter, mass, pt, key='measurement', do_fit=False,
                     sys_mean_mass=None, sys_mean_pt=None,
                     **kwargs):
        if isinstance(mass, pd.DataFrame):
            mass_mean = mass
        else:
            mass_mean = self.get_mean_df(key=key, other=mass,)

        if isinstance(pt, pd.DataFrame):
            pt_mean = pt
        else:
            pt_mean = self.get_mean_df(key=key, other=pt,)

        print(pt_mean, mass_mean)
        plotter.add_errorbar((mass_mean, pt_mean), **kwargs)  # show total uncertainty

        if do_fit:
            fitter = ISRLinearFitter(mass_mean, pt_mean)
            slope, slope_err, intercept, intercept_err = fitter.do_fit()

            # Need to extract mean of systematics
            if sys_mean_pt is None and sys_mean_mass is None:
                sys_mean_pt, sys_mean_mass = self.get_sys_mean_dfs(pt, mass, sys_name=None)
            #
            fit_slope = []
            fit_intercept = []
            for sys_name in sys_mean_pt.keys():
                for sys_index in range(len(sys_mean_pt[sys_name])):
                    fitter_sys = ISRLinearFitter(sys_mean_mass[sys_name][sys_index], sys_mean_pt[sys_name][sys_index])
                    slope_sys, _, intercept_sys, _ = fitter_sys.do_fit()
                    fit_slope.append(slope_sys-slope)
                    fit_intercept.append(intercept_sys-intercept)
            corr = np.corrcoef(fit_slope, fit_intercept)
            cov_matrix = np.cov(fit_slope, fit_intercept, ddof=0)
            print(corr)
            print(cov_matrix)
            slope_sys = np.sqrt(cov_matrix[0, 0])
            intercept_sys = np.sqrt(cov_matrix[1, 1])

            # draw fit result
            x = np.linspace(50, 400, 350)
            y = 2.0 * slope * np.log(x) + intercept
            plotter.current_axis.plot(x, y, color='black', linewidth=0.7,
                                      label=f'Fit $a={slope:.2f}\pm{slope_sys:.2f}\pm{slope_err:.2f},\ '
                                                  f'b={intercept:.2f}\mp{intercept_sys:.2f}\mp{intercept_err:.2f}$\n'
                                                  r'$y = b + 2 \cdot a  \cdot ln(x)$')
            # update labels
            plotter.update_legend((0,0))
            # add fit result
        return pt_mean, mass_mean

    def draw_isr_plot(self, other, save_and_reset_plotter=True, postfix='', key='measurement',
                      do_fit=True, save_as_csv=False, show_this_sys=None, y_min=13, y_max=29,
                      **kwargs):

        if self.is_pt and other.is_pt == False:
            isr_pt = self
            isr_mass = other
        elif self.is_pt == False and other.is_pt:
            isr_pt = other
            isr_mass = self
        else:
            raise ValueError("Cannot draw ISR plot between two ISR histograms with different dimensionality")

        measurement_hist = self.get(hist_type='acceptance_corrected', mass_window_index=0, key=key)
        plotter = measurement_hist.plotter

        plotter.init_plotter(figsize=(10,8), rows=1, cols=1)
        pt_mean, mass_mean = self.add_isr_plot(plotter, isr_mass, isr_pt, key=key, do_fit=do_fit, **kwargs,)
        # show simulation
        self.add_isr_plot(plotter, isr_mass, isr_pt, key='simulation', do_fit=False,
                          color='red', linestyle='--', linewidth=0.7, marker='o', markersize=4, capsize=3, label='MiNNLO')

        if show_this_sys:
            # loop over variations and plot the means!
            sys_mean_pt = []  # first mass window sys, second mass window sys
            sys_mean_mass = []
            sys_mean_pt, sys_mean_mass = self.get_sys_mean_dfs(isr_pt, isr_mass, sys_name=show_this_sys)

            for sys_index in range(len(sys_mean_pt[show_this_sys])):
                plotter.add_errorbar((sys_mean_mass[show_this_sys][sys_index],
                                      sys_mean_pt[show_this_sys][sys_index]), linestyle='--', linewidth=0.5,)

        if save_as_csv:
            pt_mean.to_csv(f"/Users/junhokim/Work/cms_snu/ISR/results/pt_{self.channel}{self.year}.csv")
            mass_mean.to_csv(f"/Users/junhokim/Work/cms_snu/ISR/results/mass_{self.channel}{self.year}.csv")
        plotter.set_isr_plot_cosmetics(channel=change_to_greek(self.channel), y_min=y_min, y_max=y_max)
        text = isr_mass.get_additional_text_on_plot()
        plotter.add_text(text=text, location=(0, 0), do_magic=False, **{"frameon": False, "loc": "lower right", })

        if save_and_reset_plotter:
            if key == 'simulation':
                plotter.set_experiment_label(**{'year': measurement_hist.year, 'label': 'Simulation'})
            else:
                plotter.set_experiment_label(**{'year': measurement_hist.year, })
            plotter.draw_errorbar()
            plotter.show_legend(location=(0, 0), **{"loc": "upper left"})
            plotter.save_and_reset_plotter("isr_test"+"_"+self.channel+self.year+postfix)
            return None
        else:
            return plotter

    def draw_unfold_inputs(self, mass_window_index=-1, bin_width_norm=False):
        unfold_input_hist = self.get_hist_in_mass_window('unfold_input', mass_window_index,
                                                         bin_width_norm=bin_width_norm)
        measurement_hist = self.get_hist_in_mass_window('measurement', mass_window_index,
                                                       bin_width_norm=bin_width_norm)
        signal_fake_hist = self.get_hist_in_mass_window('signal_fake', mass_window_index,
                                                        bin_width_norm=bin_width_norm)
        background_hists = self.get_hist_in_mass_window('background', mass_window_index,
                                                        bin_width_norm=bin_width_norm)

        dummy_hist = measurement_hist.create(reset_hist=True)
        total_bg = sum(background_hists.values(), dummy_hist)
        measurement_hist = measurement_hist-total_bg-signal_fake_hist

        suffix = '_unfold_input'
        draw_as_2d = False
        if self.is_2d:
            if mass_window_index == -1:
                draw_as_2d = True
                suffix = '_unfold_input_pt_mass'
            else:
                suffix = '_unfold_input_'+str(mass_window_index)
        plotter = measurement_hist.plotter
        plotter.init_plotter(rows=2, cols=1)
        plotter.set_experiment_label(**{'year': measurement_hist.year})

        kwargs = get_hist_kwargs(measurement_hist.get_label()) | {"color": "black"}
        plotter.add_hist(unfold_input_hist, **kwargs)
        kwargs = get_hist_kwargs(measurement_hist.get_label()) | {"mfc": "none", "color": "gray"}
        plotter.add_hist(measurement_hist, as_denominator=True, **kwargs)

        plotter.draw_hist()
        plotter.draw_ratio_hists(location=(1, 0))

        plotter.save_and_reset_plotter(measurement_hist.hist_name + suffix + "_" + self.channel + self.year)
        del measurement_hist
        del background_hists
        del unfold_input_hist
        del signal_fake_hist

    def draw_fake_hists(self, mass_window_index=-1, bin_width_norm=False):
        signal_fake_hist = self.get_hist_in_mass_window('signal_fake', mass_window_index,
                                                        bin_width_norm=bin_width_norm)

        plotter = signal_fake_hist.plotter
        plotter.init_plotter(rows=1, cols=1)
        plotter.set_experiment_label(**{'year': signal_fake_hist.year})

        plotter.add_hist(signal_fake_hist, as_stack=True, as_denominator=True,
                         **get_hist_kwargs(signal_fake_hist.get_label()))

        plotter.draw_hist()
        suffix = '_fake_DY'
        if self.is_2d:
            suffix = '_fake_DY_'+str(mass_window_index)
        plotter.save_and_reset_plotter(signal_fake_hist.hist_name + suffix + "_" + self.channel + self.year)
        del signal_fake_hist

    def draw_background_fractions(self, mass_window_index=-1):
        signal_hist = self.get_hist_in_mass_window("signal", mass_window_index,)
        background_hists = self.get_hist_in_mass_window("background", mass_window_index,)

        total_mc_hist = None

        for bg_name, bg_hist in background_hists.items():
            total_mc_hist = bg_hist + total_mc_hist
        total_mc_hist = total_mc_hist + signal_hist

        plotter = signal_hist.plotter
        plotter.init_plotter(rows=1, cols=1)
        plotter.set_experiment_label(label='Simulation', **{'year': signal_hist.year})

        # set fractions
        for bg_name, bg_hist in background_hists.items():
            background_hists[bg_name] = bg_hist.divide(total_mc_hist)
            plotter.add_hist(background_hists[bg_name], as_stack=False, as_denominator=True,
                             **{'label': labels.get(bg_hist.get_label(), bg_name)})

        suffix = '_bg_fraction_detector'
        if self.is_2d:
            suffix = '_bg_fraction_detector_'+str(mass_window_index)

        plotter.draw_hist()
        text = self.get_additional_text_on_plot(mass_window_index)
        plotter.add_text(text=text, location=(0, 0), do_magic=False, **{"frameon": False, "loc": "upper left"})

        plotter.set_common_ratio_plot_cosmetics(self.x_axis_label, y_axis_name='Fractions',
                                                y_min=0, y_max=0.5)
        plotter.save_and_reset_plotter(signal_hist.hist_name + suffix + "_" + self.channel + self.year)
        del signal_hist
        del background_hists

    def _draw_comparison_plot(self, measurement_hist, signal_hist, background_hists=None, text=None, suffix='',
                              draw_as_2d=False, mc_denominator=True, save_and_reset=True, signal_as_stack=True):
        plotter = measurement_hist.plotter
        plotter.init_plotter(rows=2, cols=1)
        plotter.set_experiment_label(**{'year': measurement_hist.year})

        # Add backgrounds if any
        if background_hists:
            for bg_name, bg_hist in background_hists.items():
                plotter.add_hist(bg_hist, as_stack=True, as_denominator=True,
                                 **get_hist_kwargs(bg_hist.get_label()))
                                 #**{'label': labels.get(bg_hist.get_label(), bg_name)})

        # Add signal and measurement
        kwargs =get_hist_kwargs(signal_hist.get_label())
        #if not signal_as_stack:
        #    kwargs['histtype'] = 'errorbar'
        #    kwargs['mfc'] = 'none'
        plotter.add_hist(signal_hist, as_stack=signal_as_stack, as_denominator=mc_denominator,
                         **kwargs)
        plotter.add_hist(measurement_hist, as_denominator=not mc_denominator,
                         **get_hist_kwargs(measurement_hist.get_label()))


        plotter.draw_hist()
        plotter.draw_ratio_hists(location=(1, 0))

        x_log_scale = y_log_scale = True
        if self.is_pt:
            x_log_scale = y_log_scale = False
            y_log_scale = True

        x_axis_label = self.x_axis_label
        if self.is_2d and draw_as_2d:
            x_axis_label = 'Bin index'

        ratio_name="Data/MC"
        if not mc_denominator:
            ratio_name = "MC/Data"
        plotter.set_common_comparison_plot_cosmetics(x_axis_label, x_log_scale=x_log_scale,
                                                     y_log_scale=y_log_scale, ratio_name=ratio_name)

        if text:
            plotter.add_text(text=text, location=(0, 0), do_magic=True, **{"frameon": False, "loc": "upper left"})

        if save_and_reset:
            plotter.save_and_reset_plotter(measurement_hist.hist_name + suffix + "_" + self.channel + self.year)
            return None
        else:
            return plotter

    def draw_memory_test(self):
        measurement_hist = self.get_hist_in_mass_window("measurement", 2,
                                                        bin_width_norm=True)

        plotter = measurement_hist.plotter
        plotter.init_plotter(rows=1, cols=1)
        plotter.save_and_reset_plotter(measurement_hist.hist_name + "memory_" + self.channel + self.year)

    def draw_detector_level(self, mass_window_index=-1, bin_width_norm=False):
        measurement_hist = self.get_hist_in_mass_window("measurement", mass_window_index,
                                                        bin_width_norm=bin_width_norm)
        signal_hist = self.get_hist_in_mass_window("signal", mass_window_index,
                                                   bin_width_norm=bin_width_norm)
        background_hists = self.get_hist_in_mass_window("background", mass_window_index,
                                                        bin_width_norm=bin_width_norm)
        text = self.get_additional_text_on_plot(mass_window_index)

        suffix = '_detector'
        draw_as_2d = False
        if self.is_2d:
            if mass_window_index == -1:
                draw_as_2d = True
                suffix = '_detector_pt_mass'
            else:
                suffix = '_detector_'+str(mass_window_index)
        self._draw_comparison_plot(measurement_hist, signal_hist, background_hists, text, suffix=suffix,
                                   draw_as_2d=draw_as_2d)
        del measurement_hist
        del signal_hist
        del background_hists

    def draw_unfolded_level(self, mass_window_index=-1, bin_width_norm=False, mc_denominator=True):
        measurement_hist = self.get_hist_in_mass_window("unfolded_measurement", mass_window_index,
                                                        bin_width_norm=bin_width_norm)
        signal_hist = self.get_hist_in_mass_window("truth_signal", mass_window_index,
                                                   bin_width_norm=bin_width_norm,)

        text = self.get_additional_text_on_plot(mass_window_index)

        suffix = '_unfolded'
        if self.is_2d:
            suffix = '_unfolded_'+str(mass_window_index)
        # draw unfold input and its comparison to simulation
        plotter = self._draw_comparison_plot(measurement_hist, signal_hist, text=text, suffix=suffix,
                                             mc_denominator=mc_denominator, signal_as_stack=False, save_and_reset=False)
        # folded level
        input_hist = self.get_hist_in_mass_window("unfold_input", mass_window_index,
                                                  bin_width_norm=bin_width_norm)
        sim_input_hist = self.get_hist_in_mass_window("reco_signal", mass_window_index,
                                                      bin_width_norm=bin_width_norm)

        plotter.add_hist(input_hist, as_denominator=not mc_denominator, label='Data (reco)',
                         histtype='errorbar', color='gray', mfc='none',)
        plotter.add_hist(sim_input_hist, as_denominator=mc_denominator, label='Reco DY')

        plotter.draw_hist()
        plotter.draw_ratio_hists(location=(1, 0), show_error_band=False,
                                 show_normalized_error_band=False)

        plotter.show_legend()
        plotter.save_and_reset_plotter(measurement_hist.hist_name + suffix + "_" + self.channel + self.year)

        del measurement_hist
        del signal_hist
        del input_hist
        del sim_input_hist

    def draw_unfold_closure(self, mass_window_index=-1, bin_width_norm=False):
        unfolded_signal_hist = self.get_hist_in_mass_window("unfolded_signal", mass_window_index,
                                                            bin_width_norm=bin_width_norm)
        truth_signal_hist = self.get_hist_in_mass_window("truth_signal", mass_window_index,
                                                         bin_width_norm=bin_width_norm)

        text = self.get_additional_text_on_plot(mass_window_index)

        suffix = '_closure'
        if self.is_2d:
            suffix = '_closure_'+str(mass_window_index)
        self._draw_comparison_plot(unfolded_signal_hist, truth_signal_hist, text=text, suffix=suffix)
        del unfolded_signal_hist
        del truth_signal_hist

    def draw_acceptance(self, mass_window_index=-1, bin_width_norm=False, y_min=0.0, y_max=0.55):
        full_phase_hist = self.get_hist_in_mass_window("acceptance_corrected", mass_window_index,
                                                       bin_width_norm=bin_width_norm, key='simulation')
        fiducial_phase_hist = self.get_hist_in_mass_window("truth_signal", mass_window_index,
                                                          bin_width_norm=bin_width_norm)

        plotter = full_phase_hist.plotter
        plotter.init_plotter(rows=1, cols=1)
        plotter.set_experiment_label(label='Simulation', **{'year': full_phase_hist.year})

        plotter.add_hist(full_phase_hist, as_denominator=True, not_to_draw=True, yerr=False)
        plotter.add_hist(fiducial_phase_hist, as_denominator=False, not_to_draw=True, yerr=False,
                         histtype='errorbar')
        plotter.draw_ratio_hists(location=(0, 0), show_normalized_error_band=False)

        if not self.is_pt:
            plotter.current_axis.set_xscale("log")
        plotter.current_axis.set_xlabel(self.x_axis_label)
        plotter.current_axis.set_ylabel("Acceptance")
        plotter.get_axis(location=(0, 0)).set_ylim(y_min, y_max)
        text = self.get_additional_text_on_plot(mass_window_index)
        plotter.add_text(text=text, location=(0, 0), do_magic=False, **{"frameon": False, "loc": "lower left"})

        plotter.save_and_reset_plotter(full_phase_hist.hist_name + "_acceptance" + "_" +
                                       str(mass_window_index) + "_" + self.channel + self.year)
        del full_phase_hist
        del fiducial_phase_hist

    def draw_acceptance_corrected_level(self, mass_window_index=-1, bin_width_norm=False,
                                        mc_denominator=True, others=None):
        measurement_hist = self.get_hist_in_mass_window("acceptance_corrected", mass_window_index,
                                                        bin_width_norm=bin_width_norm)
        signal_hist = self.get_hist_in_mass_window("acceptance_corrected", mass_window_index,
                                                   bin_width_norm=bin_width_norm, key='simulation')
        text = self.get_additional_text_on_plot(mass_window_index)

        suffix = '_acceptance_corrected'
        if self.is_2d:
            suffix = '_acceptance_corrected_'+str(mass_window_index)
        save_and_reset = True
        if others:
            save_and_reset = False
        plotter = self._draw_comparison_plot(measurement_hist, signal_hist, text=text, suffix=suffix,
                                             mc_denominator=mc_denominator, save_and_reset=save_and_reset,
                                             signal_as_stack=False)
        if others:
            for index, other in enumerate(others):
                other_hist = other.get('acceptance_corrected', mass_window_index, bin_width_norm=bin_width_norm,
                                       key='simulation')
                plotter.add_hist(other_hist, as_denominator=False, **get_hist_kwargs(other_hist.get_label()))
                plotter.add_hist(measurement_hist, as_denominator=True, not_to_draw=True,
                                 **get_hist_kwargs(measurement_hist.get_label()))
                plotter.draw_hist()
                plotter.draw_ratio_hists(location=(1, 0))

            plotter.show_legend()
            plotter.save_and_reset_plotter(measurement_hist.hist_name + suffix + "_" + self.channel + self.year)
        del measurement_hist
        del signal_hist

    def draw_systematic_summary(self, mass_window_index=0,):
        measurement_hist = self.isr_hists_per_mass_window[mass_window_index].acceptance_corrected_hist['measurement']
        relative_systematic_hist = measurement_hist.create_normalized_error_band()

        plotter = measurement_hist.plotter
        plotter.init_plotter(rows=1, cols=1)
        plotter.set_experiment_label(**{'year': measurement_hist.year})

        _, bins, errors = relative_systematic_hist.to_numpy()
        _, _, stat = relative_systematic_hist.to_numpy(stat=True)
        
        plotter.add_hist((errors, bins, None), as_denominator=False, yerr=False, show_err_band=False, color='black',
                         label='Total', linewidth=2)
        plotter.add_hist((stat, bins, None), as_denominator=False, yerr=False, show_err_band=False, color='black',
                         linestyle='--', label='Stat')
        # merge systematics
        systematics = relative_systematic_hist.get_systematics(merge_statistics=True)
        colors = plotter.get_n_colors(len(systematics))
        for i, sys_name_ in enumerate(systematics):
            sys_error = systematics[sys_name_]
            kwargs = {}
            kwargs['color'] = colors[i]
            if sys_name_ == 'matrix_stat':
                kwargs['linestyle'] = '--'

            plotter.add_hist((sys_error, bins, None), as_denominator=False, yerr=False, show_err_band=False,
                             label=sys_name_, **kwargs)
        plotter.draw_hist()
        if not self.is_pt:
            plotter.get_axis(location=(0, 0)).set_xscale("log")
        plotter.show_legend(**{"loc": "upper left"})
        text = self.get_additional_text_on_plot(mass_window_index)
        plotter.add_text(text=text, location=(0, 0), do_magic=True, **{"frameon": False, "loc": "upper right"})
        plotter.get_axis(location=(0, 0)).set_ylabel("Relative uncertainty")
        plotter.get_axis(location=(0, 0)).set_xlabel(self.x_axis_label)
        plotter.save_and_reset_plotter(measurement_hist.hist_name + "_sys_summary_" +
                                       str(mass_window_index) + "_" +self.channel + self.year)
        del measurement_hist
        del relative_systematic_hist

    def draw_systematic_hists(self, sys_name=None, mass_window_index=-1, bin_width_norm=False,
                              hist_type='unfolded_measurement', key='measurement'):
        measurement_hist = self.get_hist_in_mass_window(hist_type, mass_window_index,
                                                        bin_width_norm=bin_width_norm, key=key)

        if sys_name is None:
            sys_names = measurement_hist.systematic_raw_root_hists.keys()
        else:
            sys_names = [sys_name]

        for sys_name in sys_names:
            systematic_hists = measurement_hist.systematic_raw_root_hists[sys_name]

            plotter = measurement_hist.plotter
            plotter.init_plotter(rows=2, cols=1)
            plotter.set_experiment_label(**{'year': measurement_hist.year})

            plotter.add_hist(measurement_hist, as_denominator=True,
                             label='Default', histtype='errorbar', color='black', yerr=False)

            colors = ['red', 'green', 'blue', 'cyan', 'magenta', 'yellow', 'pink']
            index = 0
            for sys_name_, sys_hist in systematic_hists.items():
                # suppress y error for systematic hists
                plotter.add_hist(Hist(sys_hist), as_stack=False, as_denominator=False, yerr=False, show_err_band=False,
                                 color=colors[index%7], label=sys_name_)
                index += 1

            plotter.draw_hist()
            if index > 10:
                plotter.show_legend(ncol=5, fontsize=7)
            else:
                plotter.show_legend()
            plotter.draw_ratio_hists(location=(1, 0), show_y_error=False, show_error_band=False,)
            plotter.get_axis(location=(0,0)).set_ylabel("Events/bin")
            plotter.get_axis(location=(1,0)).set_ylabel("/Default")
            plotter.get_axis(location=(1,0)).set_xlabel(self.x_axis_label)
            plotter.get_axis(location=(1, 0)).set_ylim(0.90, 1.1)

            x_log_scale = y_log_scale = True
            if self.is_pt:
                x_log_scale = y_log_scale = False
                y_log_scale = True

            if y_log_scale:
                plotter.get_axis(location=(0, 0)).set_yscale("log")

            if x_log_scale:
                plotter.get_axis(location=(0, 0)).set_xscale("log")
                plotter.get_axis(location=(1, 0)).set_xscale("log")

            plotter.save_and_reset_plotter(measurement_hist.hist_name + "_sys_" + sys_name + "_" + str(mass_window_index) + "_" +
                                           self.channel + self.year)
        del measurement_hist

    def draw_response_matrix(self, mass_window_index=-1, out_name_postfix='',
                             x_axis_label_prefix="", y_axis_label_prefix="",
                             label_postfix="", show_number=False,
                             **kwargs_hist2dplot):
        # Choose the relevant ISR hist container
        if self.is_pt and not self.is_2d:
            isr_hist_to_draw = self.isr_hists[mass_window_index]
        else:
            isr_hist_to_draw = self.isr_hists[0]
        # Fetch the target hist or dict of hists
        attr_name = 'tunfolder'
        tunfolder = getattr(isr_hist_to_draw, attr_name)
        useAxisBinning = True
        if self.is_pt:
            useAxisBinning = False
        raw_2d, condition_number = tunfolder.get_response_matrix(useAxisBinning)  # directly from TUnfold
        plotter = raw_2d.plotter

        plotter.init_plotter(rows=1, cols=1, left=0.12, right=0.9)
        label='Simulation'
        if label_postfix:
            label = 'Simulation ' + label_postfix
        plotter.set_experiment_label(label=label, **{'year': tunfolder.year})
        if self.is_pt:
            # pre-FSR bin index
            if not x_axis_label_prefix:
                x_axis_label_prefix = "Truth bin index"
            if not y_axis_label_prefix:
                y_axis_label_prefix = "Reconstructed bin index"
            x_axis_label = x_axis_label_prefix + " $(p_{T},m)^{"+change_to_greek(self.channel)+"}$"
            y_axis_label = y_axis_label_prefix + " $(p_{T},m)^{"+change_to_greek(self.channel)+"}$"
        else:
            if not x_axis_label_prefix:
                x_axis_label_prefix = "Truth"
            if not y_axis_label_prefix:
                y_axis_label_prefix = "Reconstructed"
            x_axis_label = x_axis_label_prefix + " $m^{"+ change_to_greek(self.channel) +"}$ [GeV]"
            y_axis_label = y_axis_label_prefix + " $m^{"+ change_to_greek(self.channel) +"}$ [GeV]"
            # self.x_axis_label = r"$p_{T}^{" + change_to_greek(self.channel) + "}$ [GeV]"
        plotter.draw_matrix(raw_2d.to_numpy_2d(), x_axis_label=x_axis_label, y_axis_label=y_axis_label,
                            show_number=show_number,
                            **kwargs_hist2dplot)
        plotter.add_text("condition number: " + str(round(condition_number, 2)),
                         location=(0,0), do_magic=False,
                         **{"loc": "upper left",})
        if out_name_postfix:
            out_name_postfix = "_" + out_name_postfix
        plotter.save_and_reset_plotter(tunfolder.response_matrix.hist_name + "_RM_" + self.channel + self.year
                                       + out_name_postfix)

    def draw_correlations(self, mass_window_index=-1, **kwargs_hist2dplot):
        # Choose the relevant ISR hist container
        if self.is_pt and not self.is_2d:
            isr_hist_to_draw = self.isr_hists[mass_window_index]
        else:
            isr_hist_to_draw = self.isr_hists[0]

        # Fetch the target hist or dict of hists
        attr_name = 'tunfolder'
        tunfolder = getattr(isr_hist_to_draw, attr_name)
        useAxisBinning = True
        if self.is_pt:
            useAxisBinning = False
        raw_2d = Hist(tunfolder.get_correlation_matrix(useAxisBinning))  # directly from TUnfold
        plotter = raw_2d.plotter

        plotter.init_plotter(rows=1, cols=1, left=0.12, right=0.9)
        plotter.set_experiment_label(**{'year': tunfolder.year})
        if self.is_pt:
            x_axis_label = "Unfolded bin index $(p_{T},m)^{"+change_to_greek(self.channel)+"}$"
            y_axis_label = "Unfolded bin index $(p_{T},m)^{"+change_to_greek(self.channel)+"}$"
        else:
            x_axis_label = "Unfolded $m^{"+ change_to_greek(self.channel) +"}$ [GeV]"
            y_axis_label = "Unfolded $m^{"+ change_to_greek(self.channel) +"}$ [GeV]"
            # self.x_axis_label = r"$p_{T}^{" + change_to_greek(self.channel) + "}$ [GeV]"
        plotter.draw_matrix(raw_2d.to_numpy_2d(), x_axis_label=x_axis_label, y_axis_label=y_axis_label,
                            **kwargs_hist2dplot)
        plotter.save_and_reset_plotter(tunfolder.response_matrix.hist_name + "_correlation_" + self.channel + self.year)

    def draw_bin_efficiency(self, mass_window_index=-1,):
        if self.is_pt and not self.is_2d:
            isr_hist_to_draw = self.isr_hists[mass_window_index]
        else:
            isr_hist_to_draw = self.isr_hists[0]

        attr_name = 'tunfolder'
        tunfolder = getattr(isr_hist_to_draw, attr_name)

        x_log = True
        if self.is_pt:
            x_log = False
        plotter = tunfolder.draw_bin_efficiency(x_log=x_log, save_and_reset_plotter=False)
        plotter.show_legend()
        plotter.current_axis.set_xlabel(self.x_axis_label)
        plotter.current_axis.set_ylabel("Fraction")

        plotter.save_and_reset_plotter(tunfolder.response_matrix.hist_name + "_eff_purity_" + self.channel + self.year)

    def draw_pt_comparisons(self, *others, key='acceptance_corrected', index=2,
                            bin_width_norm=False, scale=1):
        reference_hist = self.get_hist_in_mass_window(key, index,
                                                      bin_width_norm=bin_width_norm, scale=scale)

        plotter = reference_hist.plotter
        plotter.init_plotter(rows=1, cols=1)
        plotter.set_experiment_label(**{'year': reference_hist.year})

        plotter.add_hist(reference_hist, as_denominator=True, not_to_draw=True,)
        colors = ['red', 'green', 'blue', 'cyan', 'magenta', 'yellow', 'pink']
        for i, other in enumerate(others):
            hist = other.get_hist_in_mass_window(key, index,
                                                bin_width_norm=bin_width_norm, scale=scale)
            plotter.add_hist(hist, as_denominator=False, not_to_draw=True, label=other.year,
                             color=colors[i], histtype='errorbar')

        plotter.draw_hist()
        plotter.draw_ratio_hists(location=(0, 0))

        plotter.get_axis(location=(0, 0)).set_ylim(0.5, 1.5)
        plotter.set_common_ratio_plot_cosmetics(self.x_axis_label)
        text = self.get_additional_text_on_plot(index)
        plotter.add_text(text=text, location=(0, 0), do_magic=True, **{"frameon": False, "loc": "upper left"})
        plotter.get_axis(location=(0, 0)).set_ylabel("/"+self.year)

        plotter.save_and_reset_plotter(reference_hist.hist_name + "_comparison_" + key + "_" + str(index) + "_" +
                                       self.channel + self.year)

    # comparisons between different mass windows
    def draw_pt_comparison(self):
        # use Z peak mass reason as denominator
        reference_hist = self.get('acceptance_corrected', 2, bin_width_norm=True, scale=-1)
        hist0 = self.get('acceptance_corrected', 0, bin_width_norm=True, scale=-1)
        hist1 = self.get('acceptance_corrected', 1, bin_width_norm=True, scale=-1)
        hist3 = self.get('acceptance_corrected', 3, bin_width_norm=True, scale=-1)
        hist4 = self.get('acceptance_corrected', 4, bin_width_norm=True, scale=-1)

        plotter = reference_hist.plotter
        plotter.init_plotter(rows=1, cols=1)
        plotter.set_experiment_label(**{'year': reference_hist.year})

        plotter.add_hist(reference_hist, as_denominator=True, location=-999)
        plotter.add_hist(hist0, as_denominator=False, location=-999)
        plotter.add_hist(hist1, as_denominator=False, location=-999)
        plotter.add_hist(hist3, as_denominator=False, location=-999)
        plotter.add_hist(hist4, as_denominator=False, location=-999)
        plotter.draw_ratio_hists(location=(0, 0))
        plotter.draw_hist()

        plotter.get_axis(location=(0, 0)).set_ylim(0, 3)
        plotter.save_and_reset_plotter("test_" + self.channel + self.year)
