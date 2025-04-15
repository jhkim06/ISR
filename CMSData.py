from ROOTFileGrouper import ROOTFileGrouper
from FilePather import FilePather


def get_mc_key(process_name, channel_name, generator_name='MiNNLO'):
    if process_name == "DY":
        if generator_name == "MiNNLO":
            if channel_name == "ee":
                sample_key = "DYJetsToEE_MiNNLO"
            if channel_name == "mm":
                sample_key = "DYJetsToMuMu_MiNNLO"
            if channel_name == "tautau":
                sample_key = "DYJetsToTauTau_MiNNLO"
        elif generator_name == "aMCNLO":
            sample_key = "DYJets"
        else:
            sample_key = "DYJets_MG"
    else:
        sample_key = process_name
    return sample_key
'''
example.
sys_name: bg_norm, setups to be used when get_combined_root_hists() called
'''


class CMSData(object):
    def __init__(self, sample_base_dir):
        self.experiment = 'CMS'
        self.sample_base_dir = sample_base_dir

        self.data_pather = FilePather(sample_base_dir + 'data/', 'data_config.json')
        self.mc_pather = FilePather(sample_base_dir + 'mc/', 'mc_config.json')

        # systematic
        # for measurement, mc process
        # signal: DY, background: WW, WZ,ttbar, etc


    def get_measurement(self, period_name, channel_name, event_selection):
        data_dict = self.data_pather.get_path_dict((channel_name, period_name))
        hist_path_prefix = channel_name + period_name + '/' + event_selection + '/'

        data_file_group = ROOTFileGrouper('Data', data_dict,
                                          hist_path_prefix=hist_path_prefix,
                                          experiment_name=self.experiment,
                                          year=period_name,
                                          channel_name=channel_name, is_measurement=True)
        return data_file_group


    # process_name = ()
    def get_mc(self, process_name, period_name, channel_name, event_selection, label=''):
        sample_key = get_mc_key(process_name, channel_name)
        mc_dict = self.mc_pather.get_path_dict((period_name, sample_key))
        # for tau
        hist_path_prefix = channel_name + period_name + '/' + event_selection + '/'
        hist_name_prefix =''
        if sample_key == 'DYJetsToTauTau_MiNNLO':
            hist_name_prefix = 'tau_'

        label = label or sample_key
        mc_file_group = ROOTFileGrouper(label, mc_dict,
                                        hist_path_prefix=hist_path_prefix,
                                        hist_name_prefix=hist_name_prefix,
                                        experiment_name=self.experiment, year=period_name,
                                        channel_name=channel_name, is_measurement=False)
        return mc_file_group
