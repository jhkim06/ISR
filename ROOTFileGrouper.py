import ROOT


class ROOTFileGrouper:
    def __init__(self, group_name,
                 file_dict, hist_path_prefix='', hist_name_prefix=''):

        self.group_name = group_name
        self.file_dict = file_dict  # dictionary of file paths
        self.hist_path_prefix = hist_path_prefix
        self.hist_name_prefix = hist_name_prefix

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

    def get_combined_root_hists(self, hist_name, hist_path='', use_local_hist_path=False,
                                bin_width_norm=False):
        if use_local_hist_path:
            hist_path = hist_path
        else:
            hist_path = self.hist_path_prefix
        hist_name = self.hist_name_prefix + hist_name
        hist_total = None
        ROOT.TH1.AddDirectory(False)
        # loop over files and open histogram with hist path and add them
        file_dict = self.file_dict
        for file_label, file_path in file_dict.items():

            file = ROOT.TFile.Open(file_path, 'r')
            hist = file.Get(hist_path + hist_name)
            file.Close()

            if hist_total is None:
                hist_total = hist.Clone(self.group_name + hist_name)
            else:
                hist_total.Add(hist)

        if bin_width_norm:
            hist_total.Scale(1, "width")
        return hist_total
