from ISRAnalyzer import ISRAnalyzer
from ISRCombiner import ISRCombiner
import logging
from Plotter import Plotter
import ROOT
import sys
from ISRLinearFitter import ISRLinearFitter
import numpy as np
import pandas as pd
ROOT.gROOT.SetBatch(True)


def draw_isr_plot_from_df(mass, pt, save_and_reset_plotter=True, channel_label='ll', postfix='', **kwargs):

        plotter = Plotter('CMS',
                          '/Users/junhokim/Work/cms_snu/ISR/results' )
        plotter.init_plotter(figsize=(10,8), rows=1, cols=1)
        plotter.set_experiment_label(year='Run 2')

#        plotter.add_errorbar((mass, pt), **kwargs)
        #plotter.add_errorbar((mass_1d.get_df(key=key), pt_1d.get_df(key=key)),
        #                     color=red, label='NNLO', linestyle='dashdot', linewidth=1.)
        #plotter.add_errorbar((mass_nlo_1d.get_df(key=key), pt_nlo_1d.get_df(key=key)),
        #                     color=soft_blue, label='NLO', linestyle='dashdot', linewidth=1.)
        #plotter.add_errorbar((mass_lo_1d.get_df(key=key), pt_lo_1d.get_df(key=key)),
        #                     color=soft_orange, label='LO', linestyle='dashdot', linewidth=1.)
        fitter = ISRLinearFitter(mass, pt)
        slope, slope_err, intercept, intercept_err = fitter.do_fit()
        slope_err = 0.34
        intercept_err = 2.91

        # draw fit result
        x = np.linspace(50, 400, 350)
        y = 2.0 * slope * np.log(x) + intercept
        plotter.current_axis.plot(x, y, label=r'$y = b + 2 \cdot a  \cdot ln(x)$'
                                              f'\nFit (CMS 13TeV) $a={slope:.2f}\pm{slope_err:.2f},\ '
                                              f'b={intercept:.2f}\mp{intercept_err:.2f}$\n')
        plotter.update_legend((0,0))

        #slope_cms_8tev = 3.18
        #intercept_cms_8tev = -12.17
        #y_cms_8tev = 2.0 * slope_cms_8tev * np.log(x) + intercept_cms_8tev
        #plotter.current_axis.plot(x, y_cms_8tev, label=f'Fit (CMS 8TeV) $a={slope_cms_8tev:.2f}\pm{slope_err:.2f},\ '
        #                                               f'b={intercept_cms_8tev:.2f}\pm{intercept_err:.2f}$\n')
        #plotter.update_legend((0,0))

        slope_cdf = 2.15
        slope_cdf_err = 0.09
        intercept_cdf = -7.56
        intercept_cdf_err = 0.83
        y_cdf = 2. * slope_cdf * np.log(x) + intercept_cdf
        plotter.current_axis.plot(x, y_cdf, label=f'Fit (CDF 1.96 TeV) $a={slope_cdf:.2f}\pm{slope_cdf_err:.2f},\ '
                                                  f'b={intercept_cdf:.2f}\mp{intercept_cdf_err:.2f}$\n')
        plotter.update_legend((0,0))

        plotter.set_isr_plot_cosmetics(channel=channel_label,)
        plotter.get_axis(location=(0, 0)).set_ylim(9, 29)
        text = r"$p_{T}^{"+ channel_label+"}<$ 100 GeV"
        plotter.add_text(text=text, location=(0, 0), do_magic=False, **{"frameon": False, "loc": "lower right", })

        if save_and_reset_plotter:
            plotter.draw_errorbar()
            plotter.show_legend(location=(0, 0), loc='upper left')
            plotter.save_and_reset_plotter("isr_test" + postfix)
            return None
        else:
            return plotter

mass_combined_final = pd.read_csv("/Users/junhokim/Work/cms_snu/ISR/results/mass_combined.csv", index_col=0)
pt_combined_final = pd.read_csv("/Users/junhokim/Work/cms_snu/ISR/results/pt_combined.csv", index_col=0)

draw_isr_plot_from_df(mass_combined_final, pt_combined_final,
                      label=r'$ee$ and $\mu\mu$ combined',
                      color='black', marker='o', markersize=4.0, linestyle='none', linewidth=0.7,)
