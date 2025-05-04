import ROOT
from Hist import Hist
#ROOT.gErrorIgnoreLevel = ROOT.kWarning


# Merge root files to handel them simultaneously
class ROOTHistGrouper:
    def __init__(self, group_name,
                 pather,
                 experiment_name='', period_name='', channel_name='',
                 is_measurement=True,
                 mc_process_names=None,
                 hist_name_prefix=''):

        self.period_name = period_name
        self.channel_name = channel_name
        self.group_name = group_name
        self.pather = pather
        self.mc_process_names = mc_process_names
        # data_dict = self.data_pather.get_path_dict((channel_name, period_name))

        self.experiment_name = experiment_name
        self.is_measurement = is_measurement

        self.hist_path_prefix = f"{channel_name}{period_name}/"
        self.hist_name_prefix = hist_name_prefix
        # hist_path_prefix = channel_name + period_name + '/' + event_selection + '/'
        # need to have information on systematics

    def get_name(self):
        return self.group_name

    def get_experiment_name(self):
        return self.experiment_name

    def get_year(self):
        return self.period_name

    def get_channel_name(self):
        return self.channel_name

    def get_file_dict(self, sys_dir_name='default'):

        # TODO for mc, need to handle multiple case
        return self.pather.get_path_dict(*self.mc_process_names,
                                         sys_dir_name=sys_dir_name)

    def get_tobject(self, object_name, sys_dir_name='default'):
        file_dict = self.get_file_dict(sys_dir_name=sys_dir_name)

        file_path = file_dict[list(file_dict.keys())[0]]  # get the first item in the dictionary
        file = ROOT.TFile.Open(file_path, 'r')
        tobject = file.Get(object_name)
        file.Close()

        return tobject

    def get_combined_root_hists(self,
                                hist_name,
                                event_selection,
                                sys_dir_name='default',
                                hist_name_prefix='',
                                is_signal=False,
                                bin_width_norm=False,
                                norm=False,
                                scale=1.0, raw_hist=False):

        hist_path = self.hist_path_prefix + event_selection + "/"
        # hist_name_prefix related to the event selection
        if hist_name_prefix != '':
            hist_name = hist_name_prefix + hist_name
        # hist_name_prefix related to the mc sample
        if self.hist_name_prefix != '':
            hist_name = self.hist_name_prefix + hist_name
        hist_total = None
        ROOT.TH1.AddDirectory(False)
        # loop over files and open histogram with hist path and add them
        # get file_dict here?
        file_dict = self.get_file_dict(sys_dir_name)
        # TODO systematic?
        for file_label, file_path in file_dict.items():
            file = ROOT.TFile.Open(file_path, 'r')
            hist = file.Get(hist_path + hist_name)
            hist.GetNbinsX()
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
            return Hist(hist_total,
                        hist_name=hist_name,
                        label=self.group_name,
                        channel=self.channel_name,
                        year=self.period_name,
                        is_measurement=self.is_measurement,
                        is_mc_signal=is_signal)
