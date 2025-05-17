from Hist import Hist
import copy


# ISRTUnfoldBinningHist
# for additional info on tunfold
class HistTUnfoldBin(Hist):
    def __init__(self, hist,  # 1D hist
                 tunfold_bin):

        super(HistTUnfoldBin, self).__init__(hist.get_raw_hist(),
                                             hist_name=hist.hist_name,
                                             experiment=hist.experiment,
                                             label=hist.label,
                                             channel=hist.channel,
                                             year=hist.year,
                                             is_measurement=hist.is_measurement,
                                             is_mc_signal=hist.is_mc_signal, )
        self.systematic_raw_root_hists = copy.deepcopy(hist.systematic_raw_root_hists)
        self.compute_systematic_rss_per_sysname()

        self.tunfold_bin = tunfold_bin

    def extract_1d_raw_hist(self, raw_hist, index):
        projection_mode = 'dipt[O];dimass[UOC'+str(index)+']'  # this if for dipt-dimass 'dipt[O]'
        use_axis_binning = True

        # central raw hist
        extracted_raw_hist = self.extract_raw_hist(raw_hist, projection_mode, use_axis_binning)
        return extracted_raw_hist
        #extracted_raw_hist = self.tunfold_bin.ExtractHistogram("unfold_input_extracted",
        #                                                      raw_hist,
        #                                                      0, # error matrix
        #                                                      use_axis_binning,
        #                                                      projection_mode)

    def extract_raw_hist(self, raw_hist, projection_mode, use_axis_binning=True):

        # central raw hist
        extracted_raw_hist = self.tunfold_bin.ExtractHistogram("unfold_input_extracted",
                                                              raw_hist,
                                                              0,  # error matrix
                                                              use_axis_binning,
                                                              projection_mode)
        return extracted_raw_hist

    def extract_hist(self, bin_width_norm=False):
        projection_mode = 'dimass[UO]'
        extracted_raw_hist = self.extract_raw_hist(self.raw_root_hist, projection_mode, use_axis_binning=True)
        # extracted_raw_hist.Scale(1, "width")

        extracted_hist = Hist(
            extracted_raw_hist,
            hist_name=self.hist_name,
            label=self.label,
            channel=self.channel,
            year=self.year,
            is_measurement=self.is_measurement,
            is_mc_signal=self.is_mc_signal
        )

        # for systematics
        extracted_systematic_raw_hists = {}
        for sys_name, variations in self.systematic_raw_root_hists.items():
            extracted_systematic_raw_hists[sys_name] = {}
            for var_name, hist in variations.items():
                extracted_systematic_raw_hists[sys_name][var_name] = self.extract_raw_hist(hist, projection_mode,
                                                                                           use_axis_binning=True)
        extracted_hist.systematic_raw_root_hists = copy.deepcopy(extracted_systematic_raw_hists)
        extracted_hist.compute_systematic_rss_per_sysname()

        return extracted_hist.bin_width_norm(1.) if bin_width_norm else extracted_hist

    # CAUTION bin_width_norm will change binned mean values a lot, so here just stick to False
    def extract_1d_hist(self, index, bin_width_norm=False):
        extracted_raw_hist = self.extract_1d_raw_hist(self.raw_root_hist, index)
        # extracted_raw_hist.Scale(1, "width")

        extracted_hist = Hist(
            extracted_raw_hist,
            hist_name=self.hist_name,
            label=self.label,
            channel=self.channel,
            year=self.year,
            is_measurement=self.is_measurement,
            is_mc_signal=self.is_mc_signal
        )

        # for systematics
        extracted_systematic_raw_hists = {}
        for sys_name, variations in self.systematic_raw_root_hists.items():
            extracted_systematic_raw_hists[sys_name] = {}
            for var_name, hist in variations.items():
                extracted_systematic_raw_hists[sys_name][var_name] = self.extract_1d_raw_hist(hist, index)
                # extracted_systematic_raw_hists[sys_name][var_name].Scale(1, "width")

        extracted_hist.systematic_raw_root_hists = copy.deepcopy(extracted_systematic_raw_hists)
        extracted_hist.compute_systematic_rss_per_sysname()

        return extracted_hist.bin_width_norm(1.) if bin_width_norm else extracted_hist