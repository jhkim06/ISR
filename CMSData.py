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


class CMSData(object):
    def __init__(self, sample_base_dir):
        self.experiment = 'CMS'
        self.sample_base_dir = sample_base_dir

        self.data_pather = FilePather(sample_base_dir + 'data/', 'data_config.json')
        self.mc_pather = FilePather(sample_base_dir + 'mc/', 'mc_config.json')

        # systematic
    # signal: DY, background: WW, WZ,ttbar, etc

    def get_measurement(self, period_name, channel_name, event_selection):
        data_dict = self.data_pather.get_path_dict((channel_name, period_name))
        hist_path_prefix = channel_name + period_name + '/' + event_selection + '/'

        data_file_group = ROOTFileGrouper('Data', data_dict,
                                          hist_path_prefix=hist_path_prefix,
                                          experiment_name=self.experiment,
                                          year=period_name,
                                          channel_name=channel_name)

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
                                        hist_name_prefix=hist_name_prefix)
        return mc_file_group

    def get_isr_samples(self,
                        period_name,
                        channel_name,  # used to set hist path
                        lepton_selection='TightID_b_veto',  # used to set hist path
                        use_minnlo=True):

        # Data
        data_dict = self.data_pather.get_path_dict((channel_name, period_name))

        # expectation
        # Signal
        # Note that mc contains both channels in general
        if use_minnlo:
            if channel_name == 'ee':
                dy_key = 'DYJetsToEE_MiNNLO'
            else:
                dy_key = 'DYJetsToMuMu_MiNNLO'
            dy_tau_key = 'DYJetsToTauTau_MiNNLO'
        else:
            dy_key = 'DYJets'
            dy_tau_key = 'DYJetsToTauTau_MiNNLO'
        dy_dict = self.mc_pather.get_path_dict((period_name, dy_key))

        # Background
        # check the format of the dictionary
        dy_tau_dict = self.mc_pather.get_path_dict((period_name, dy_tau_key))
        ttbar_dict = self.mc_pather.get_path_dict((period_name, 'TTLL'))
        vv_dict = self.mc_pather.get_path_dict((period_name, 'ZZ'),
                                               (period_name, 'WZ'),
                                               (period_name, 'WW'))
        gg_dict = self.mc_pather.get_path_dict((period_name, 'GGLL'))

        # common hist path prefix
        hist_path_prefix = channel_name + period_name + '/' + lepton_selection + '/'

        # here set experiment name, year, channel
        # data_file_group = ROOTFileGrouper('Data', 'ee', '2016a', data_dict, hist_path_prefix)
        data_file_group = ROOTFileGrouper('Data', data_dict,
                                          hist_path_prefix=hist_path_prefix,
                                          experiment_name=self.experiment,
                                          year=period_name, channel_name=channel_name)
        dy_file_group = ROOTFileGrouper('DY', dy_dict,
                                        hist_path_prefix=hist_path_prefix)
        dy_tau_file_group = ROOTFileGrouper('tau', dy_tau_dict,
                                            hist_path_prefix=hist_path_prefix,
                                            hist_name_prefix='tau_')
        ttbar_file_group = ROOTFileGrouper('ttbar', ttbar_dict,
                                           hist_path_prefix=hist_path_prefix)
        vv_file_group = ROOTFileGrouper('vv', vv_dict,
                                        hist_path_prefix=hist_path_prefix)
        gg_file_group = ROOTFileGrouper('gg', gg_dict,
                                        hist_name_prefix=hist_path_prefix)

        return (
            # FIXME use group_name of ROOTFileGrouper
            data_file_group,  # data
            dy_file_group,  # signal
            [
                # this will be the order to draw later
                gg_file_group,
                ttbar_file_group,
                vv_file_group,
                dy_tau_file_group,  # backgrounds
            ]
        )
