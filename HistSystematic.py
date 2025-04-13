from Hist import Hist
from Plotter import Plotter
import numpy as np


class HistSystematic:
    def __init__(self, systematic_name,
                 nominal_hist,
                 systematic_hists,
                 norm_to_nominal_hist=False,):
        # systematic_hists is dictionary {variation name: systematic histogram}
        self.systematic_name = systematic_name
        self.nominal_hist = nominal_hist
        self.systematic_hists = systematic_hists
        self.norm_to_nominal_hist = norm_to_nominal_hist

        if norm_to_nominal_hist:
            for _, var_hist in self.systematic_hists.items():
                var_hist.Scale(self.nominal_hist.Integral()/var_hist.Integral())

        self.nominal_values = Hist(self.nominal_hist).to_numpy()[0]
        self.template = self.nominal_hist.Clone()
        self.template.Reset()

        # numpy array
        self.sym_sys = self.create_sys_np_array()
        self.asym_sys = self.create_sys_np_array()
        self.calculate_sym_systematic()

        # root histograms
        self.sym_sys_up_hist = self.template.Clone("sys_up")
        self.sym_sys_down_hist = self.template.Clone("sys_down")
        self.asym_sys_up_hist = self.template.Clone("asym_sys_up")
        self.asym_sys_down_hist = self.template.Clone("asym_sys_down")

    def create_sys_np_array(self):
        # [[up errors], [down errors]]
        temp = np.zeros(self.nominal_values.shape)

        temp = np.expand_dims(temp, axis=0)
        temp = np.append(temp, temp, axis=0)
        return temp

    def calculate_sym_systematic(self):
        temp_sys = np.zeros(self.nominal_values.shape)
        for variation_name, variation_hist in self.systematic_hists.items():
            temp_delta = Hist(variation_hist).to_numpy()[0]-self.nominal_values
            temp_sys = np.sqrt(np.square(temp_delta) + np.square(temp_sys))
        self.sym_sys[0] = temp_sys
        self.sym_sys[1] = temp_sys

    def calculate_asym_systematic(self):
        # set asym_sys
        pass

    def draw_systematic(self):

        plotter = Plotter('CMS', './Plots')  # FIXME use self.plotter
        plotter.create_subplots(2, 1, figsize=(8, 8),
                                left=0.15, right=0.95, hspace=0.0, bottom=0.15, height_ratios=[1, 0.3])

        nominal_hist = Hist(self.nominal_hist.Clone("nominal_hist"))
        # nominal_hist.Scale(1./nominal_hist.Integral(), 'width')
        plotter.set_experiment_label(**{"year": "year"})
        # loop over systematic hist
        for variation_name, variation_hist in self.systematic_hists.items():
            variation_hist = Hist(variation_hist.Clone("variation_hist"))
            # variation_hist.Scale(1./variation_hist.Integral(), 'width')
            plotter.add_comparison_pair(variation_hist, denominator_hist=nominal_hist,
                                        location=(0, 0), ratio_location=(1, 0),
                                        nominator_args={'label': variation_name, 'color': 'red'},
                                        denominator_args={"histtype": 'errorbar', 'marker': "o", "color": 'black',
                                                          'mfc': 'none',
                                                          "label": 'Nominal'})
        plotter.draw_error_boxes(self.nominal_values, Hist(self.nominal_hist).to_numpy()[1], self.sym_sys,
                                 sys_name='FSR', **{"facecolor": 'black', "alpha": 0.2, "fill": True,
                                                    'hatch': None})

        plotter.draw_hist()
        plotter.set_common_comparison_plot_cosmetics("variable", bin_width_norm=True,
                                                     ratio_name='/Nominal')
        plotter.show_legend(location=(0, 0))
        plotter.get_axis(location=(0, 0)).set_xticklabels([])
        # plotter.adjust_y_scale()
        plotter.get_axis(location=(1, 0)).set_ylim(0.8, 1.2)
        plotter.save_fig("_systematic_test")
