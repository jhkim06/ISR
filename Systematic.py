from Hist import Hist
from Plotter import Plotter


class Systematic:
    def __init__(self, systematic_name, nominal_hist, systematic_hists):
        # systematic_hists is dictionary {variation name: systematic histogram}
        self.systematic_name = systematic_name
        self.nominal_hist = nominal_hist
        self.systematic_hists = systematic_hists

    def calculate_systematic(self):
        # calculate systematic in each bin
        # use root-mean-square of deviations
        pass

    def draw_systematic(self):

        plotter = Plotter('CMS', './Plots')  # FIXME use self.plotter
        plotter.create_subplots(2, 1, figsize=(8, 8),
                                left=0.15, right=0.95, hspace=0.0, bottom=0.15, height_ratios=[1, 0.3])

        nominal_hist = self.nominal_hist.Clone("nominal_hist")
        nominal_hist.Scale(1./nominal_hist.Integral(), 'width')
        plotter.set_experiment_label(**{"year": "year"})
        # loop over systematic hist
        for variation_name, variation_hist in self.systematic_hists.items():
            variation_hist = variation_hist.Clone("variation_hist")
            variation_hist.Scale(1./variation_hist.Integral(), 'width')
            plotter.add_comparison_pair(variation_hist, denominator_hist=nominal_hist, location=(0, 0), ratio_location=(1, 0),
                                        nominator_args={'label': variation_name, 'color': 'red'},
                                        denominator_args={"histtype": 'errorbar', 'marker': "o", "color": 'black', 'mfc': 'none',
                                                          "label": 'Nominal'})

        plotter.draw_hist()
        plotter.comparison_plot_cosmetics("variable", bin_width_norm=True,
                                          ratio_name='/Nominal')
        plotter.show_legend(location=(0, 0))
        plotter.get_axis(location=(0, 0)).set_xticklabels([])
        # plotter.adjust_y_scale()
        plotter.get_axis(location=(1, 0)).set_ylim(0.8, 1.2)
        plotter.save_fig("_systematic_test")
