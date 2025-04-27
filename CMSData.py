from ROOTHistGrouper import ROOTHistGrouper
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

        self.data_pather = FilePather(sample_base_dir, 'data', 'data_config.json')
        self.mc_pather = FilePather(sample_base_dir, 'mc', 'mc_config.json')

    def get_measurement(self, period_name, channel_name):

        keys_tuple = [(channel_name, period_name)]
        data_file_group = ROOTHistGrouper('Data',
                                          self.data_pather,
                                          mc_process_names=keys_tuple,
                                          channel_name=channel_name,
                                          period_name=period_name,
                                          experiment_name=self.experiment,
                                          is_measurement=True)
        return data_file_group


    # process_name = ()
    def get_mc(self, process_name, period_name, channel_name, label=''):
        def parse_process(proc):
            parts = proc.split(":", 1)
            name = parts[0]
            gen = parts[1] if len(parts) == 2 else 'MiNNLO'
            return name, gen

        hist_name_prefix = ''
        # Handle multiple processes
        if isinstance(process_name, tuple):
            process_info = [parse_process(p) for p in process_name]
            keys_tuple = [(period_name, get_mc_key(p, channel_name, generator_name=g)) for p, g in process_info]
            sample_key = '+'.join([k[0] for k in process_info])
            label = label or sample_key
        else:
            # Single process case
            sample_key = parse_process(process_name)
            sample_key = get_mc_key(sample_key[0], channel_name, generator_name=sample_key[1])

            keys_tuple = [(period_name, sample_key)]
            hist_name_prefix = 'tau_' if sample_key == 'DYJetsToTauTau_MiNNLO' else ''
            label = label or sample_key

        mc_file_group = ROOTHistGrouper(label,
                                        self.mc_pather,
                                        mc_process_names=keys_tuple,
                                        hist_name_prefix=hist_name_prefix,
                                        experiment_name=self.experiment, period_name=period_name,
                                        channel_name=channel_name, is_measurement=False)
        return mc_file_group