import ROOT
from Hist import Hist
#ROOT.gErrorIgnoreLevel = ROOT.kWarning


# Merge root files to handel them simultaneously
class ROOTFileGrouper:
    def __init__(self, group_name, file_dict,
                 experiment_name='', year='', channel_name='', is_measurement=True,
                 hist_path_prefix='', hist_name_prefix=''):

        self.group_name = group_name
        self.file_dict = file_dict  # dictionary of file paths

        self.experiment_name = experiment_name
        self.year = year
        self.channel_name = channel_name
        self.is_measurement = is_measurement

        self.hist_path_prefix = hist_path_prefix
        self.hist_name_prefix = hist_name_prefix

        # need to have information on systematics

    def get_name(self):
        return self.group_name

    def get_experiment_name(self):
        return self.experiment_name

    def get_year(self):
        return self.year

    def get_channel_name(self):
        return self.channel_name

    def list_samples(self):
        for label in self.file_dict.keys():
            print(label)

    def set_hist_path_prefix(self, hist_path):
        self.hist_path_prefix = hist_path

    def get_hist_path_prefix(self):
        return self.hist_path_prefix

    def get_tobject(self, name):

        file_path = self.file_dict[list(self.file_dict.keys())[0]]  # just first item in the dictionary
        file = ROOT.TFile.Open(file_path, 'r')
        tobject = file.Get(name)
        file.Close()

        return tobject

    def get_combined_root_hists(self,
                                hist_name,
                                # FIXME it's better allow abstract options here or hide these options
                                hist_sys_postfix = '',
                                hist_path='',
                                use_local_hist_path=False,
                                bin_width_norm=False,
                                is_signal=False,
                                norm=False,
                                scale=1.0, raw_hist=False):
        if use_local_hist_path:
            hist_path = hist_path
        else:
            hist_path = self.hist_path_prefix
        hist_name = self.hist_name_prefix + hist_name
        hist_total = None
        ROOT.TH1.AddDirectory(False)
        # loop over files and open histogram with hist path and add them
        file_dict = self.file_dict
        # TODO systematic?
        for file_label, file_path in file_dict.items():

            file = ROOT.TFile.Open(file_path, 'r')
            hist = file.Get(hist_path + hist_name)
            file.Close()

            if hist_total is None:
                # TODO update proper error handle
                try:
                    hist_total = hist.Clone(self.group_name + hist_name)
                except:
                    pass
            else:
                hist_total.Add(hist)

        if bin_width_norm:
            hist_total.Scale(1.0, "width")
        if scale != 1.0:
            hist_total.Scale(scale)
        # normalize to
        if norm:
            hist_total.Scale(1./hist_total.Integral())

        if raw_hist:
            return hist_total
        else:
            # Add systematics?
            return Hist(hist_total,
                        hist_name=hist_name,
                        label=self.group_name,
                        channel=self.channel_name,
                        year=self.year,
                        is_measurement=self.is_measurement,
                        is_mc_signal=is_signal)
