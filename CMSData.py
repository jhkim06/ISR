from ROOTFileGrouper import ROOTFileGrouper
from FilePather import FilePather


class CMSData(object):
    def __init__(self, sample_base_dir):
        self.sample_base_dir = sample_base_dir

        self.data_pather = FilePather(sample_base_dir + 'data/', 'data_config.json')
        self.mc_pather = FilePather(sample_base_dir + 'mc/', 'mc_config.json')

        # systematic

    # TODO use txt sample information
    def get_isr_samples(self, period_name, channel_name, lepton_selection='TightID_b_veto', use_minnlo=True):

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
            dy_tau_key = 'DYJets'
        dy_dict = self.mc_pather.get_path_dict((period_name, dy_key))

        # Background
        dy_tau_dict = self.mc_pather.get_path_dict((period_name, dy_tau_key))
        ttbar_dict = self.mc_pather.get_path_dict((period_name, 'TTLL'))
        vv_dict = self.mc_pather.get_path_dict((period_name, 'ZZ'),
                                               (period_name, 'WZ'),
                                               (period_name, 'WW'))
        gg_dict = self.mc_pather.get_path_dict((period_name, 'GGLL'))

        # common hist path prefix
        hist_path_prefix = channel_name + period_name + '/' + lepton_selection + '/'

        data_file_group = ROOTFileGrouper('Data', data_dict, hist_path_prefix)
        dy_file_group = ROOTFileGrouper('DY', dy_dict, hist_path_prefix)
        dy_tau_file_group = ROOTFileGrouper('tau', dy_tau_dict, hist_path_prefix,
                                            hist_name_prefix='tau_')
        ttbar_file_group = ROOTFileGrouper('ttbar', ttbar_dict, hist_path_prefix)
        vv_file_group = ROOTFileGrouper('vv', vv_dict, hist_path_prefix)
        gg_file_group = ROOTFileGrouper('gg', gg_dict, hist_path_prefix)

        return (
            ('Data', data_file_group),  # data
            ('Drell-Yan', dy_file_group),  # signal
            [('tau', dy_tau_file_group),  # backgrounds
             ('ttbar', ttbar_file_group),
             ('VV', vv_file_group),
             ('gg', gg_file_group)]
        )
