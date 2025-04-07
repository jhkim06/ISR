from Combiner import Combiner

class ISRCombiner(Combiner):
    def __init__(self):

        # get_isr_result_df
        self.names = []  # keys for mass and pt dataframe
        self.results_mass_dfs = {}
        self.results_pt_dfs = {}

        self.combiner = Combiner()

    def get_results_dfs(self, name, isr_result):
        # experiment, channel, period
        self.names.append(name)
        self.results_mass_dfs[name] = isr_result[0]
        self.results_pt_dfs[name] = isr_result[1]

    def combine(self):
        #
        pass
