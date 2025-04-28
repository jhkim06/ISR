from Combiner import Combiner
import pandas as pd
import numpy as np


def make_combiner_input(df, i):
    index_name = df.iloc[i].index[:-1]  # drop total error
    values = df.iloc[i].values[:-1]

    combiner_input = [(index_name[index], values[index]) for index in range(len(index_name))]
    combiner_input = (combiner_input[0], combiner_input[1:])
    return combiner_input


class ISRCombiner:
    def __init__(self, is_same_channel=True):

        self.is_same_channel = is_same_channel
        # get_isr_result_df
        self.names = []  # keys for mass and pt dataframe
        self.results_mass_dfs = {}
        self.results_pt_dfs = {}

    def get_results_dfs(self, name, isr_mass, isr_pt):
        # experiment, channel, period
        self.names.append(name)
        self.results_mass_dfs[name] = isr_mass
        self.results_pt_dfs[name] = isr_pt

    def combine(self):
        first_key = self.names[0]
        mass_df = self.results_mass_dfs[first_key]
        pt_df = self.results_pt_dfs[first_key]

        # create DF
        number_of_rows = mass_df.shape[0]

        combined_mass_df = pd.DataFrame(np.nan, index=mass_df.index, columns=mass_df.columns)
        combined_pt_df = pd.DataFrame(np.nan, index=pt_df.index, columns=pt_df.columns)

        # loop over mass windows
        for i in range(number_of_rows):
            # make an input for Combiner
            combiner_input_mass = []
            combiner_input_pt = []

            # loop over measurements
            for key, _ in self.results_mass_dfs.items():
                combiner_input = make_combiner_input(self.results_mass_dfs[key], i)  # i: mass window
                combiner_input_mass.append(combiner_input)
                combiner_input = make_combiner_input(self.results_pt_dfs[key], i)
                combiner_input_pt.append(combiner_input)

            # combine!
            # get solved combiner
            mass_combiner = Combiner("Mean mass", combiner_input_mass,
                                     use_rho_same_channel=self.is_same_channel, solve=True)
            pt_combiner = Combiner("Mean pt", combiner_input_pt,
                                   use_rho_same_channel=self.is_same_channel, solve=True)

            # update dataframe
            combined_mass_df.loc[i, "mean":] = mass_combiner.get_result_list()
            combined_pt_df.loc[i, "mean":] = pt_combiner.get_result_list()
        return combined_mass_df, combined_pt_df