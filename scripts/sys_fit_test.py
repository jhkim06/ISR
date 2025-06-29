import pandas as pd
import numpy as np
import Plotter as Plotter
import pickle

pt_combined_df = {}
mass_combined_df = {}
pt_combined_df['ee'] = pd.read_csv("/Users/junhokim/Work/cms_snu/ISR/results/pt_ee_combined.csv", index_col=0)
mass_combined_df['ee'] = pd.read_csv("/Users/junhokim/Work/cms_snu/ISR/results/mass_ee_combined.csv", index_col=0)

pt_combined_df['mm'] = pd.read_csv("/Users/junhokim/Work/cms_snu/ISR/results/pt_mm_combined.csv", index_col=0)
mass_combined_df['mm'] = pd.read_csv("/Users/junhokim/Work/cms_snu/ISR/results/mass_mm_combined.csv", index_col=0)

sys_pt_combined_df = {}
sys_mass_combined_df = {}

for channel in ['ee', 'mm']:
    with open(f'/Users/junhokim/Work/cms_snu/ISR/{channel}_combined_pt_sys.pkl', 'rb') as f:
        sys_pt_combined_df[channel] = pickle.load(f)
    with open(f'/Users/junhokim/Work/cms_snu/ISR/{channel}_combined_mass_sys.pkl', 'rb') as f:
        sys_mass_combined_df[channel] = pickle.load(f)

pt_final_df = pd.read_csv("/Users/junhokim/Work/cms_snu/ISR/results/pt_combined.csv", index_col=0)
mass_final_df = pd.read_csv("/Users/junhokim/Work/cms_snu/ISR/results/mass_combined.csv", index_col=0)

with open(f'/Users/junhokim/Work/cms_snu/ISR/combined_pt_sys.pkl', 'rb') as f:
    final_sys_pt_combined_df = pickle.load(f)
with open(f'/Users/junhokim/Work/cms_snu/ISR/combined_mass_sys.pkl', 'rb') as f:
    final_sys_mass_combined_df = pickle.load(f)

pt_df = {}
mass_df = {}

pt_df["ee"] = {}
pt_df["mm"] = {}

mass_df["ee"] = {}
mass_df["mm"] = {}

sys_pt_df = {}
sys_mass_df = {}

sys_pt_df["ee"] = {}
sys_pt_df["mm"] = {}

sys_mass_df["ee"] = {}
sys_mass_df["mm"] = {}

for period in ["2016a", "2016b", "2017", "2018"]:
    for channel in ['ee', 'mm']:
        pt_df[channel][period] = pd.read_csv(f"/Users/junhokim/Work/cms_snu/ISR/results/pt_{channel}{period}.csv", index_col=0)
        mass_df[channel][period] = pd.read_csv(f"/Users/junhokim/Work/cms_snu/ISR/results/mass_{channel}{period}.csv", index_col=0)

        with open(f'/Users/junhokim/Work/cms_snu/ISR/{channel}{period}_pt_sys.pkl', 'rb') as f:
            sys_pt_df[channel][period] = pickle.load(f)
        with open(f'/Users/junhokim/Work/cms_snu/ISR/{channel}{period}_mass_sys.pkl', 'rb') as f:
            sys_mass_df[channel][period] = pickle.load(f)



from Plotter import Plotter
from ISRLinearFitter import ISRLinearFitter

channel = 'ee'
year = '2016b'
#mass_mean = mass_df[channel][year]
#pt_mean = pt_df[channel][year]
#sys_mean_mass = sys_mass_df[channel][year]
#sys_mean_pt = sys_pt_df[channel][year]
#sys_names = sys_pt_df[channel][year].keys()


mass_mean = mass_combined_df[channel]
pt_mean = pt_combined_df[channel]
sys_mean_mass = sys_mass_combined_df[channel]
sys_mean_pt = sys_pt_combined_df[channel]
sys_names = sys_pt_combined_df[channel].keys()

#mass_mean = mass_final_df
#pt_mean = pt_final_df
#sys_mean_mass = final_sys_mass_combined_df
#sys_mean_pt = final_sys_pt_combined_df
#sys_names = final_sys_pt_combined_df.keys()

slope_sys_err_new = 0
intercept_sys_err_new = 0
fitter_default = ISRLinearFitter(mass_mean, pt_mean, sys_mean_mass, sys_mean_pt)
slope, slope_err, intercept, intercept_err = fitter_default.do_fit()
slope_sys_err_new, intercept_sys_err_new = fitter_default.do_sys_fit()  # do sys fit after nominal fit

plotter = Plotter('CMS',
                  '/Users/junhokim/Work/cms_snu/ISR/Plots')

for show_this_sys in sys_names:


    plotter.init_plotter(figsize=(10,8), rows=1, cols=1)
    plotter.add_errorbar((mass_mean, pt_mean),)
    # current version: repeat fitting for each systematic
    # nominal fit

    # draw fit result
    # TODO option to show chi2
    x = np.linspace(50, 400, 350)
    y = 2.0 * slope * np.log(x) + intercept
    chi2 = f'($\chi^{2}$: {fitter_default.chi2:.2f}, NDOF: {fitter_default.ndof})'
    print(chi2)
    label = (
        rf"Fit {chi2}"
        "\n"  # newline
        r"$y = b + 2\,a\,\ln(x)$"
        "\n"
        rf"$a = {slope:.2f}\pm{slope_err:.2f}\pm{slope_sys_err_new:.2f}$"
        "\n"
        rf"$b = {intercept:.2f}\pm{intercept_err:.2f}\mp{intercept_sys_err_new:.2f}$"
    )
    plotter.current_axis.plot(x, y, color='black', linewidth=0.7,
                              label=label)
    plotter.update_legend((0,0))

    if show_this_sys:
        # loop over variations and plot the means!

        for sys_index in range(len(sys_mean_pt[show_this_sys])):
            plotter.add_errorbar((sys_mean_mass[show_this_sys][sys_index],
                                  sys_mean_pt[show_this_sys][sys_index]), linestyle='--', linewidth=0.5,)

    plotter.set_isr_plot_cosmetics(channel='ee', y_min=13, y_max=29)
    plotter.draw_errorbar()

    plotter.show_legend(location=(0, 0), **{"loc": "upper left"})
    plotter.save_and_reset_plotter("sys_fit_test_" + show_this_sys + "_" + year)


