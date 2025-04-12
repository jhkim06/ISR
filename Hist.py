import numpy as np
import ROOT

# ROOT histogram to
class Hist(object):
    def __init__(self, hist, label='', channel='', year=''):
        self.raw_root_hist = hist

        self.label = label
        self.year = year
        self.channel = channel

    def __add__(self, other=None, c1=1):
        added_hist = self.raw_root_hist.Clone("added")
        if other is None:
            return Hist(added_hist, self.label)
        else:
            added_hist.Add(other.raw_root_hist, c1)
            return Hist(added_hist, self.label)

    def __sub__(self, other=None):
        return self.__add__(other, -1)

    def divide(self, other=None):
        divided_hist = self.raw_root_hist.Clone("divided")
        if other is None:
            return Hist(divided_hist, self.label)
        else:
            divided_hist.Divide(other.raw_root_hist)
            return Hist(divided_hist, self.label)

    def multiply(self, other=None):
        multiplied_hist = self.raw_root_hist.Clone("multiplied")
        if other is None:
            return Hist(multiplied_hist, self.label)
        else:
            multiplied_hist.Multiply(other.raw_root_hist)
            return Hist(multiplied_hist, self.label)

    def get_label(self):
        return self.label

    def get_raw_hist(self):
        return self.raw_root_hist

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

    def to_numpy(self):
        values = []
        bins = []
        errors = []

        # check if bins have labels
        n_bins_x = self.raw_root_hist.GetNbinsX()
        for i_bin in range(n_bins_x):
            value = self.raw_root_hist.GetBinContent(i_bin + 1)
            error = self.raw_root_hist.GetBinError(i_bin + 1)
            values.append(value)
            errors.append(error)
            bins.append(self.raw_root_hist.GetXaxis().GetBinLowEdge(i_bin + 1))

        bins.append(bins[-1] + self.raw_root_hist.GetBinWidth(n_bins_x))

        values = np.array(values)
        bins = np.array(bins)
        errors = np.array(errors)  # stat

        return values, bins, errors

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

    def get_bin_widths(self):
        pass
