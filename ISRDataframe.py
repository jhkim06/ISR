import pandas as pd
import numpy as np


def calculate_squared_root_sum(raw_df, reg_expression, new_col_name="total error"):

    selected_cols = raw_df.filter(regex=reg_expression)
    squared_sum = (selected_cols ** 2).sum(axis=1)
    root_of_squared_sum = np.sqrt(squared_sum)

    raw_df[new_col_name] = root_of_squared_sum


class ISRDataframe(object):
    def __init__(self, mass_bins, mean_mass, mean_pt):

        self.mass_bins = mass_bins
        self.mean_mass = mean_mass
        self.mean_pt = mean_pt

    def get_isr_dataframe(self):

        dipt_dict_list = []
        dimass_dict_list = []

        for ith_mass in range(len(self.mass_bins)):

            low_mass_edge = self.mass_bins[ith_mass][0]
            high_mass_edge = self.mass_bins[ith_mass][1]

            mass_window_str = str(low_mass_edge) + '-' + str(high_mass_edge)
            dipt_dict = dict(mass_window=mass_window_str, mean=self.mean_pt[ith_mass][0])
            dimass_dict = dict(mass_window=mass_window_str, mean=self.mean_mass[ith_mass][0])

            dipt_dict['stat error'] = self.mean_pt[ith_mass][1]
            dimass_dict['stat error'] = self.mean_mass[ith_mass][1]

            dipt_dict_list.append(dipt_dict)
            dimass_dict_list.append(dimass_dict)

        dimass_df = pd.DataFrame(dimass_dict_list)
        dipt_df = pd.DataFrame(dipt_dict_list)

        reg_expression = r".*error"
        calculate_squared_root_sum(dipt_df, reg_expression)
        calculate_squared_root_sum(dimass_df, reg_expression)

        return dimass_df, dipt_df




