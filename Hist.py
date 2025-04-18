import numpy as np
import ROOT
import copy


def to_numpy(raw_root_hist):
    values = []
    bins = []
    errors = []

    # check if bins have labels
    n_bins_x = raw_root_hist.GetNbinsX()
    for i_bin in range(n_bins_x):
        value = raw_root_hist.GetBinContent(i_bin + 1)
        error = raw_root_hist.GetBinError(i_bin + 1)
        values.append(value)
        errors.append(error)
        bins.append(raw_root_hist.GetXaxis().GetBinLowEdge(i_bin + 1))

    bins.append(bins[-1] + raw_root_hist.GetBinWidth(n_bins_x))

    values = np.array(values)
    bins = np.array(bins)
    errors = np.array(errors)  # stat

    return values, bins, errors


# ROOT histogram to
class Hist(object):
    def __init__(self, hist,
                 hist_name='',
                 experiment='CMS',
                 label='',  # to be used in legend of histogram
                 channel='',
                 year='',
                 is_measurement=True,
                 is_mc_signal=False,):

        self.raw_root_hist = hist
        # TODO how to handle systematic of TH2
        self.is_TH2 = isinstance(self.raw_root_hist, ROOT.TH2)

        self.hist_name = hist_name
        self.experiment = experiment
        self.label = label
        self.year = year
        self.channel = channel
        self.is_measurement = is_measurement
        self.is_mc_signal = is_mc_signal
        # dictionary for HistSystematic
        self.systematic_raw_root_hists = {}  # sys_name, sys_variation_name, histogram
        self.systematics = {}  # sys_name, up/down variation calculated in Hist

        # only for TH1D, init the total_sym_sys using statistical error
        if not self.is_TH2:
            # add statistical error to total_sym_sys
            self.total_sym_sys = self.create_sys_np_array(to_numpy(self.raw_root_hist)[2])

        # to avoid ImportError
        from Plotter import Plotter
        # TODO how to set base_output_dir
        self.plotter = Plotter(self.experiment,
                               '/Users/junhokim/Work/cms_snu/ISR/Plots')

    def create(self,
               hist=None,
               hist_name='',
               label='',  # to be used in legend of histogram
               year='',
               ):
        if hist is None:
            hist = self.raw_root_hist
        if hist_name == '':
            hist_name = self.hist_name
        if label == '':
            label = self.label
        if year == '':
            year = self.year

        new_hist = Hist(hist,
                        hist_name=hist_name,
                        label=label,
                        year=year,
                        channel=self.channel,
                        experiment=self.experiment,
                        is_measurement=self.is_measurement,
                        is_mc_signal=self.is_mc_signal,
                        )
        return new_hist



    def draw(self, plotter=None, name_postfix='', **kwargs):
        if plotter:
            # add hist to the given plotter
            plotter.add_hist(self)

        else:
            # inti plotter
            self.plotter.init_plotter()
            self.plotter.add_hist(self, **kwargs)
            self.plotter.draw_hist()
            self.plotter.get_axis(location=(0, 0)).set_yscale("log")
            self.plotter.get_axis(location=(0, 0)).set_xscale("log")
            self.plotter.save_and_reset_plotter(self.hist_name, name_postfix)

    def create_sys_np_array(self, error):
        # [[up errors], [down errors]]
        error_ = np.expand_dims(error, axis=0)
        error_ = np.append(error_, error_, axis=0)
        return error_

    def total_error(self):
        if self.is_TH2: return

        total_squared = to_numpy(self.raw_root_hist)[2] ** 2  # stat
        # self.total_sym_sys = self.create_sys_np_array(total_squared)
        for sys_name, sys_array in self.systematics.items():
            total_squared += sys_array ** 2
        total_error = np.sqrt(total_squared)
        return total_error

    def update_sys_np_array(self):
        if self.is_TH2: return

        total_error = self.total_error()
        self.total_sym_sys = self.create_sys_np_array(total_error)

    def compute_systematic_rss_per_sysname(self):
        if self.is_TH2: return

        central_values = to_numpy(self.raw_root_hist)[0]

        for sys_name, variations in self.systematic_raw_root_hists.items():
            squared_sum = np.zeros_like(central_values)
            for var_name, hist in variations.items():
                variation_values = to_numpy(hist)[0]
                squared_diff = (variation_values - central_values) ** 2
                squared_sum += squared_diff

            self.systematics[sys_name] = np.sqrt(squared_sum)  # RSS for this sys_name
        self.update_sys_np_array()

    def set_systematic_hist(self, sys_name, sys_variation_name, raw_hist):

        if self.is_TH2: return

        central_values = to_numpy(self.raw_root_hist)[0]
        variation_values = to_numpy(raw_hist)[0]
        squared_diff = (variation_values - central_values) ** 2

        # Add to systematic container
        if sys_name not in self.systematic_raw_root_hists:
            self.systematic_raw_root_hists[sys_name] = {}
            self.systematics[sys_name] = np.sqrt(squared_diff)
        else:
            existing_squared = self.systematics[sys_name] ** 2
            self.systematics[sys_name] = np.sqrt(existing_squared + squared_diff)

        self.systematic_raw_root_hists[sys_name][sys_variation_name] = raw_hist
        self.update_sys_np_array()

    def normalize_systematic_raw_hists_by_central(self, raw_hist):
        for sys_name, variations in self.systematic_raw_root_hists.items():
            for var_name, hist in variations.items():
                normalized = hist.Clone(f"{hist.GetName()}_norm")
                normalized.Divide(raw_hist)
                self.systematic_raw_root_hists[sys_name][var_name] = normalized

    def create_ratio_error_band_hist(self):
        norm_default = self.raw_root_hist.Clone("norm_default")
        norm_default.Divide(norm_default)
        ratio_error_band = self.create(hist=norm_default,)
        ratio_error_band.systematic_raw_root_hists = copy.deepcopy(self.systematic_raw_root_hists)
        ratio_error_band.normalize_systematic_raw_hists_by_central(self.raw_root_hist)
        ratio_error_band.compute_systematic_rss_per_sysname()

        return ratio_error_band

    def _apply_systematics_operation(self, other, operation, postfix, scale=1):
        result = {}
        if other is not None:
            self_keys = set(self.systematic_raw_root_hists.keys())
            other_keys = set(other.systematic_raw_root_hists.keys())
        else:
            return copy.deepcopy(self.systematic_raw_root_hists)

        if self_keys != other_keys:
            print(f"⚠️ Warning: Systematic keys differ:\n  self = {self_keys}\n  other = {other_keys}")

        for sys_name in self_keys & other_keys:
            self_vars = set(self.systematic_raw_root_hists[sys_name].keys())
            other_vars = set(other.systematic_raw_root_hists[sys_name].keys())
            if self_vars != other_vars:
                print(
                    f"⚠️ Warning: Variation keys differ in '{sys_name}':\n  self = {self_vars}\n  other = {other_vars}")

        for sys_name, variations in self.systematic_raw_root_hists.items():
            result[sys_name] = {}
            for var_name, hist in variations.items():
                new_hist = hist.Clone(f"{sys_name}_{var_name}_{postfix}")
                if other and sys_name in other.systematic_raw_root_hists and var_name in \
                        other.systematic_raw_root_hists[sys_name]:
                    if operation == "Add":
                        new_hist.__getattribute__(operation)(other.systematic_raw_root_hists[sys_name][var_name], scale)
                    else:
                        new_hist.__getattribute__(operation)(other.systematic_raw_root_hists[sys_name][var_name])
                result[sys_name][var_name] = new_hist

        return result

    def __add__(self, other=None, c1=1):
        added_hist = self.raw_root_hist.Clone("added")

        if other is not None:
            added_hist.Add(other.raw_root_hist, c1)

        result = self.create(hist=added_hist)
        # Initialize systematics
        result.systematic_raw_root_hists = self._apply_systematics_operation(other, "Add", "added",
                                                                             scale=c1)
        result.compute_systematic_rss_per_sysname()
        return result

    def __sub__(self, other=None):
        return self.__add__(other, -1)

    def divide(self, other=None):
        divided_hist = self.raw_root_hist.Clone("divided")

        if other is not None:
            divided_hist.Divide(other.raw_root_hist)

        result = self.create(hist=divided_hist)

        result.systematic_raw_root_hists = self._apply_systematics_operation(other, "Divide", "divided")
        result.compute_systematic_rss_per_sysname()
        return result

    def multiply(self, other=None):
        multiplied_hist = self.raw_root_hist.Clone("multiplied")

        if other is not None:
            multiplied_hist.Multiply(other.raw_root_hist)

        result = self.create(hist=multiplied_hist)
        result.systematic_raw_root_hists = self._apply_systematics_operation(other, "Multiply", "multiplied")
        result.compute_systematic_rss_per_sysname()
        return result

    def get_label(self):
        return self.label

    def get_raw_hist(self):
        return self.raw_root_hist

    def normalize(self):
        self.raw_root_hist.Scale(1. / self.raw_root_hist.Integral())

        for sys_name, variations in self.systematic_raw_root_hists.items():
            for var_name, hist in variations.items():
                hist.Scale(1. / hist.Integral())
                self.systematic_raw_root_hists[sys_name][var_name] = hist


    def get_mean(self, binned_mean=True, range_min=None, range_max=None):
        if binned_mean:
            if range_min is None and range_max is None:
                range_min = self.raw_root_hist.GetXaxis().GetXmin()
                range_max = self.raw_root_hist.GetXaxis().GetXmax()
                self.raw_root_hist.GetXaxis().SetRangeUser(range_min, range_max)
                mean, mean_error = self.raw_root_hist.GetMean(), self.raw_root_hist.GetMeanError()
                self.raw_root_hist.GetYaxis().SetRangeUser(0, 0)
            else:
                self.raw_root_hist.GetXaxis().SetRangeUser(range_min, range_max)
                mean, mean_error = self.raw_root_hist.GetMean(), self.raw_root_hist.GetMeanError()
                self.raw_root_hist.GetYaxis().SetRangeUser(0, 0)
            return mean, mean_error
        else:
            return self.raw_root_hist.GetMean(), self.raw_root_hist.GetMeanError()

    def to_numpy(self, stat=False):
        values, bins, stat_error = to_numpy(self.raw_root_hist)
        if stat:
            error = stat_error
        else:
            error = self.total_error()
        return values, bins, error

    def to_numpy_2d(self):
        content_list = []
        for i_x in range(self.raw_root_hist.GetNbinsX()):
            reco_list = []
            for i_y in range(self.raw_root_hist.GetNbinsY()):
                reco_list.append(self.raw_root_hist.GetBinContent(i_x + 1, i_y + 1))
            content_list.append(reco_list)
        content_np = np.array(content_list)

        x_bin_edges = []
        for i_x in range(self.raw_root_hist.GetNbinsX() + 1):
            x_bin_edges.append(self.raw_root_hist.GetXaxis().GetBinLowEdge(i_x + 1))
        x_bin_edges_np = np.array(x_bin_edges)

        y_bin_edges = []
        for i_y in range(self.raw_root_hist.GetNbinsY() + 1):
            y_bin_edges.append(self.raw_root_hist.GetYaxis().GetBinLowEdge(i_y + 1))
        y_bin_edges_np = np.array(y_bin_edges)

        return content_np, x_bin_edges_np, y_bin_edges_np
