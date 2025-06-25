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
import matplotlib.colors as mcolors
import pickle


logging.basicConfig(level=logging.INFO)


def unfold_and_correct(analyzer, period, channel, event_selection, is_2d=True):
    analyzer.setup_isr_detector_hists(period, channel, event_selection, use_2d_pt=is_2d)
    analyzer.isr_unfolds()
    analyzer.isr_acceptance_corrections()
    return analyzer.get_isr_results()


soft_blue =    (86 / 255, 180 / 255, 233 / 255)
soft_orange =  (230 / 255, 159 / 255, 0 / 255)
soft_green =   (0 / 255, 158 / 255, 115 / 255)
red =  (213 / 255, 94 / 255, 0 / 255)


mass_bins_1d = [(55.0, 64.0),
                (55.0, 68.0),
                #(55.0, 81.0),
                (64.0, 81.0),
                (72.0, 91.0),
                (81.0, 101.0),
                (91.0, 110.0),
                #(101.0, 150.0),
                (101.0, 200.0),
                (106.0, 220.0),
                (110.0, 243.0),
                (115.0, 273.0),
                (120.0, 320.0),
                (126.0, 380.0),
                (133.0, 440.0),
                (141.0, 510.0),
                (150.0, 600.0),
                (160.0, 700.0),
                (171.0, 830.0),
                (200.0, 1000.0),
                ]

def draw_isr_plot_from_df(mass, pt, sys_mean_mass, sys_mean_pt, 
                          save_and_reset_plotter=True, channel_label='ll', postfix='', **kwargs):

        plotter = Plotter('CMS',
                          '/Users/junhokim/Work/cms_snu/ISR/Plots' )
        plotter.init_plotter(figsize=(10,8), rows=1, cols=1)
        plotter.set_experiment_label(year='Run 2')

        plotter.add_errorbar((mass, pt), **kwargs)
        #plotter.add_errorbar((mass_1d.get_df(key=key), pt_1d.get_df(key=key)),
        #                     color=red, label='NNLO', linestyle='dashdot', linewidth=1.)
        #plotter.add_errorbar((mass_nlo_1d.get_df(key=key), pt_nlo_1d.get_df(key=key)),
        #                     color=soft_blue, label='NLO', linestyle='dashdot', linewidth=1.)
        #plotter.add_errorbar((mass_lo_1d.get_df(key=key), pt_lo_1d.get_df(key=key)),
        #                     color=soft_orange, label='LO', linestyle='dashdot', linewidth=1.)

        fitter = ISRLinearFitter(mass, pt, sys_mean_mass, sys_mean_pt)
        slope, slope_err, intercept, intercept_err = fitter.do_fit()
        slope_sys_err_new = 0
        intercept_sys_err_new = 0
        if sys_mean_mass and sys_mean_pt:
            slope_sys_err_new, intercept_sys_err_new = fitter.do_sys_fit()
                    
        # draw fit result
        x = np.linspace(50, 400, 350)
        y = 2.0 * slope * np.log(x) + intercept
        chi2 = f'($\chi^{2}$: {fitter.chi2:.2f}, NDOF: {fitter.ndof})'
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

        plotter.set_isr_plot_cosmetics(channel=channel_label,)
        text = r"$p_{T}^{"+ channel_label+"}<$ 100 GeV"
        plotter.add_text(text=text, location=(0, 0), do_magic=False, **{"frameon": False, "loc": "lower right", })

        if save_and_reset_plotter:
            plotter.draw_errorbar()
            plotter.show_legend(location=(0, 0), loc='upper left')
            plotter.save_and_reset_plotter("isr_test" + postfix)
            return None
        else:
            return plotter


def main():
    sample_base_dir = '/Users/junhokim/Work/cms_snu/data/Ultralegacy/'

    mass_bins = [(55.0, 64.0),
                 (64.0, 81.0),
                 (81.0, 101.0),
                 (101.0, 200.0),
                 (200.0, 1000.0)]
    pt_bins = (0.0, 100.0)

    ## -------------------
    ## Now combine and plot
    ## -------------------
    # Combiner for "ee" channel
    combiner_ee = ISRCombiner() 
    # Combiner for "mm" channel
    combiner_mm = ISRCombiner()

    for channel in ["ee", "mm"]:
        for period in ["2016a", "2016b", "2017", "2018"]:
            pt = pd.read_csv(f"/Users/junhokim/Work/cms_snu/ISR/results/pt_{channel}{period}.csv", index_col=0)
            mass = pd.read_csv(f"/Users/junhokim/Work/cms_snu/ISR/results/mass_{channel}{period}.csv", index_col=0)
            if channel == "ee":
                combiner_ee.get_results_dfs(channel+period, mass, pt)
            elif channel == "mm":
                combiner_mm.get_results_dfs(channel+period, mass, pt)


    sys_mass = {}
    sys_pt = {}
    for period in ["2016a", "2016b", "2017", "2018"]:
        with open(f'/Users/junhokim/Work/cms_snu/ISR/results/ee{period}_pt_sys.pkl', 'rb') as f:
            sys_pt[period] = pickle.load(f)
        with open(f'/Users/junhokim/Work/cms_snu/ISR/results/ee{period}_mass_sys.pkl', 'rb') as f:
            sys_mass[period] = pickle.load(f)
            

    sys_mass_ee_combined_df = {}
    sys_pt_ee_combined_df = {}

    for sys_name in sys_mass['2016a'].keys():
        sys_mass_ee_combined_df[sys_name] = {}
        sys_pt_ee_combined_df[sys_name] = {}

        for index in sys_mass['2016a'][sys_name].keys():
            combiner_temp = ISRCombiner()
            for period in ['2016a', '2016b', '2017', '2018']:
                combiner_temp.get_results_dfs(period + sys_name + str(index),
                                              sys_mass[period][sys_name][index],
                                              sys_pt[period][sys_name][index])

            mass_combined_ee_temp, pt_combined_ee_temp = combiner_temp.combine() # combined sys
            sys_mass_ee_combined_df[sys_name][index] = mass_combined_ee_temp
            sys_pt_ee_combined_df[sys_name][index] = pt_combined_ee_temp
            del combiner_temp

    del sys_pt
    del sys_mass

    with open(f'./results/ee_combined_mass_sys.pkl', 'wb') as f:
        pickle.dump(sys_mass_ee_combined_df, f)

    with open(f'./results/ee_combined_pt_sys.pkl', 'wb') as f:
        pickle.dump(sys_pt_ee_combined_df, f)

    sys_mass = {}
    sys_pt = {}
    for period in ["2016a", "2016b", "2017", "2018"]:
        with open(f'/Users/junhokim/Work/cms_snu/ISR/results/mm{period}_pt_sys.pkl', 'rb') as f:
            sys_pt[period] = pickle.load(f)
        with open(f'/Users/junhokim/Work/cms_snu/ISR/results/mm{period}_mass_sys.pkl', 'rb') as f:
            sys_mass[period] = pickle.load(f)


    sys_mass_mm_combined_df = {}
    sys_pt_mm_combined_df = {}

    for sys_name in sys_mass['2016a'].keys():
        sys_mass_mm_combined_df[sys_name] = {}
        sys_pt_mm_combined_df[sys_name] = {}

        for index in sys_mass['2016a'][sys_name].keys():
            combiner_temp = ISRCombiner()
            for period in ['2016a', '2016b', '2017', '2018']:
                combiner_temp.get_results_dfs(period + sys_name + str(index),
                                              sys_mass[period][sys_name][index],
                                              sys_pt[period][sys_name][index])

            mass_combined_mm_temp, pt_combined_mm_temp = combiner_temp.combine() # combined sys
            sys_mass_mm_combined_df[sys_name][index] = mass_combined_mm_temp
            sys_pt_mm_combined_df[sys_name][index] = pt_combined_mm_temp
            del combiner_temp

    del sys_pt
    del sys_mass

    with open(f'./results/mm_combined_mass_sys.pkl', 'wb') as f:
        pickle.dump(sys_mass_mm_combined_df, f)

    with open(f'./results/mm_combined_pt_sys.pkl', 'wb') as f:
        pickle.dump(sys_pt_mm_combined_df, f)

    # Combine all periods for each channel
    mass_combined_ee, pt_combined_ee = combiner_ee.combine()
    pt_combined_ee.to_csv(f"/Users/junhokim/Work/cms_snu/ISR/results/pt_ee_combined.csv", float_format='%.4f')
    mass_combined_ee.to_csv(f"/Users/junhokim/Work/cms_snu/ISR/results/mass_ee_combined.csv", float_format='%.4f')
    del combiner_ee

    draw_isr_plot_from_df(mass_combined_ee, pt_combined_ee, sys_mass_ee_combined_df, sys_pt_ee_combined_df, postfix='ee', channel_label='ee', 
                          label=r'$ee$ combined data',
                          linestyle='none', marker='o', color='black', ms=5, zorder=1001, capsize=3)

    mass_combined_mm, pt_combined_mm = combiner_mm.combine()
    pt_combined_mm.to_csv(f"/Users/junhokim/Work/cms_snu/ISR/results/pt_mm_combined.csv", float_format='%.4f')
    mass_combined_mm.to_csv(f"/Users/junhokim/Work/cms_snu/ISR/results/mass_mm_combined.csv", float_format='%.4f')
    del combiner_mm

    draw_isr_plot_from_df(mass_combined_mm, pt_combined_mm, sys_mass_mm_combined_df, sys_pt_mm_combined_df, postfix='mm', channel_label='\mu\mu',
                          label=r'$\mu\mu$ combined data',
                          linestyle='none', marker='o', color='black', ms=5, zorder=1001, capsize=3)

    # ===========================================================================================================================================
    sys_mass_combined_df = {}
    sys_pt_combined_df = {}

    muon_sys_names = sys_pt_mm_combined_df.keys()
    electron_sys_names = sys_pt_ee_combined_df.keys()

    total_sys_names = muon_sys_names | electron_sys_names
    #slice_list = ['mean', 'stat', 'sys', 'total_error']

    for sys_name in total_sys_names:
        if sys_name in sys_pt_mm_combined_df and sys_name in sys_pt_ee_combined_df:
            sys_mass_combined_df[sys_name] = {}
            sys_pt_combined_df[sys_name] = {}

            for index in sys_pt_mm_combined_df[sys_name].keys():
                combiner_temp = ISRCombiner()
                # muon
                combiner_temp.get_results_dfs("mm"+sys_name + str(index),
                                              sys_mass_mm_combined_df[sys_name][index],
                                              sys_pt_mm_combined_df[sys_name][index])
                # electron
                combiner_temp.get_results_dfs("ee"+sys_name + str(index),
                                              sys_mass_ee_combined_df[sys_name][index],
                                              sys_pt_ee_combined_df[sys_name][index])

                mass_combined_temp, pt_combined_temp = combiner_temp.combine() # combined sys
                sys_mass_combined_df[sys_name][index] = mass_combined_temp
                sys_pt_combined_df[sys_name][index] = pt_combined_temp
                del combiner_temp

        elif sys_name in sys_pt_mm_combined_df and sys_name not in sys_pt_ee_combined_df:
            sys_mass_combined_df[sys_name] = {}
            sys_pt_combined_df[sys_name] = {}
            for index in sys_pt_mm_combined_df[sys_name].keys():
                combiner_temp = ISRCombiner()
                combiner_temp.get_results_dfs("mm"+sys_name + str(index),
                                              sys_mass_mm_combined_df[sys_name][index],
                                              sys_pt_mm_combined_df[sys_name][index])

                # use ee default
                combiner_temp.get_results_dfs("ee"+sys_name + str(index),
                                              #mass_combined_ee[slice_list],
                                              #pt_combined_ee[slice_list])
                                              mass_combined_ee,
                                              pt_combined_ee)

                mass_combined_temp, pt_combined_temp = combiner_temp.combine() # combined sys
                sys_mass_combined_df[sys_name][index] = mass_combined_temp
                sys_pt_combined_df[sys_name][index] = pt_combined_temp  # check if they are inside total error?
                del combiner_temp
            
        elif sys_name not in sys_pt_mm_combined_df and sys_name in sys_pt_ee_combined_df:
            sys_mass_combined_df[sys_name] = {}
            sys_pt_combined_df[sys_name] = {}
            for index in sys_pt_ee_combined_df[sys_name].keys():
                combiner_temp = ISRCombiner()
                combiner_temp.get_results_dfs("ee"+sys_name + str(index),
                                              sys_mass_ee_combined_df[sys_name][index],
                                              sys_pt_ee_combined_df[sys_name][index])

                # use mm default
                combiner_temp.get_results_dfs("mm"+sys_name + str(index),
                                              #mass_combined_mm[slice_list],
                                              #pt_combined_mm[slice_list])
                                              mass_combined_mm,
                                              pt_combined_mm)

                mass_combined_temp, pt_combined_temp = combiner_temp.combine() # combined sys
                sys_mass_combined_df[sys_name][index] = mass_combined_temp
                sys_pt_combined_df[sys_name][index] = pt_combined_temp
                del combiner_temp

    with open(f'./results/combined_mass_sys.pkl', 'wb') as f:
        pickle.dump(sys_mass_combined_df, f)

    with open(f'./results/combined_pt_sys.pkl', 'wb') as f:
        pickle.dump(sys_pt_combined_df, f)
           

    # Global Combiner (ee + mm)
    cols = ['electronRECOSF','electronIDSF']
    mass_combined_ee['efficiency'] = np.sqrt(mass_combined_ee[cols].pow(2).sum(axis=1))
    pt_combined_ee['efficiency'] = np.sqrt(pt_combined_ee[cols].pow(2).sum(axis=1))
    print(pt_combined_ee)

    mass_combined_ee.drop(columns=['electronRECOSF'], inplace=True)
    mass_combined_ee.drop(columns=['electronIDSF'], inplace=True)
    pt_combined_ee.drop(columns=['electronRECOSF'], inplace=True)
    pt_combined_ee.drop(columns=['electronIDSF'], inplace=True)

    efficiency = pt_combined_ee.pop('efficiency')
    pt_combined_ee.insert(3, 'efficiency', efficiency)

    efficiency = mass_combined_ee.pop('efficiency')
    mass_combined_ee.insert(3, 'efficiency', efficiency)

    pt_combined_ee.to_csv(f"/Users/junhokim/Work/cms_snu/ISR/results/pt_ee_combined_.csv", float_format='%.4f')
    mass_combined_ee.to_csv(f"/Users/junhokim/Work/cms_snu/ISR/results/mass_ee_combined_.csv", float_format='%.4f')

    # MUON
    # combine 'roccor_scale','roccor_stat'
    cols = ['roccor_scale','roccor_stat']
    mass_combined_mm['momentum_scale'] = np.sqrt(mass_combined_mm[cols].pow(2).sum(axis=1))
    pt_combined_mm['momentum_scale'] = np.sqrt(pt_combined_mm[cols].pow(2).sum(axis=1))

    mass_combined_mm.drop(columns=['roccor_scale'], inplace=True)
    mass_combined_mm.drop(columns=['roccor_stat'], inplace=True)
    pt_combined_mm.drop(columns=['roccor_scale'], inplace=True)
    pt_combined_mm.drop(columns=['roccor_stat'], inplace=True)

    momentum_scale = pt_combined_mm.pop('momentum_scale')
    pt_combined_mm.insert(4, 'momentum_scale', momentum_scale)

    momentum_scale = mass_combined_mm.pop('momentum_scale')
    mass_combined_mm.insert(4, 'momentum_scale', momentum_scale)

    # change name for roccor_resolution
    momentum_res = pt_combined_mm.pop('roccor_resolution')
    pt_combined_mm.insert(5, 'momentum_resolution', momentum_res)

    momentum_res = mass_combined_mm.pop('roccor_resolution')
    mass_combined_mm.insert(5, 'momentum_resolution', momentum_res)

    efficiency = pt_combined_mm.pop('muonIDSF')
    pt_combined_mm.insert(6, 'efficiency', efficiency)

    efficiency = mass_combined_mm.pop('muonIDSF')
    mass_combined_mm.insert(6, 'efficiency', efficiency)

    pt_combined_mm.to_csv(f"/Users/junhokim/Work/cms_snu/ISR/results/pt_mm_combined_.csv", float_format='%.4f')
    mass_combined_mm.to_csv(f"/Users/junhokim/Work/cms_snu/ISR/results/mass_mm_combined_.csv", float_format='%.4f')
    # ============================================================================================= #

    pt_combined_ee_ = pd.read_csv(f"/Users/junhokim/Work/cms_snu/ISR/results/pt_ee_combined_.csv", index_col=0)
    mass_combined_ee_ = pd.read_csv(f"/Users/junhokim/Work/cms_snu/ISR/results/mass_ee_combined_.csv", index_col=0)

    pt_combined_mm_ = pd.read_csv(f"/Users/junhokim/Work/cms_snu/ISR/results/pt_mm_combined_.csv", index_col=0)
    mass_combined_mm_ = pd.read_csv(f"/Users/junhokim/Work/cms_snu/ISR/results/mass_mm_combined_.csv", index_col=0)

    combiner_all = ISRCombiner(is_same_channel=False)
    combiner_all.get_results_dfs("combined_ee", mass_combined_ee_, pt_combined_ee_)                          # Note it is important to match columns
    combiner_all.get_results_dfs("combined_mm", mass_combined_mm_[mass_combined_ee_.keys()], pt_combined_mm_[pt_combined_ee_.keys()]) 

    mass_combined_final, pt_combined_final = combiner_all.combine()
    pt_combined_final.to_csv(f"/Users/junhokim/Work/cms_snu/ISR/results/pt_combined.csv", float_format='%.4f')
    mass_combined_final.to_csv(f"/Users/junhokim/Work/cms_snu/ISR/results/mass_combined.csv", float_format='%.4f')

    #pt_combined_final = pd.read_csv("/Users/junhokim/Work/cms_snu/ISR/results/pt_combined.csv", index_col=0)
    #mass_combined_final = pd.read_csv("/Users/junhokim/Work/cms_snu/ISR/results/mass_combined.csv", index_col=0)
    #sys_mass_combined_df = None
    #sys_pt_combined_df = None
    draw_isr_plot_from_df(mass_combined_final, pt_combined_final, sys_mass_combined_df, sys_pt_combined_df,  
                          label=r'$ee$ and $\mu\mu$ combined data',
                          linestyle='none', marker='o', color='black', ms=5, zorder=1001, capsize=3)

    periods = ["2016a", "2016b", "2017", "2018"]

    color_map = {
        "2016a": (86 / 255, 180 / 255, 233 / 255),
        "2016b": (230 / 255, 159 / 255, 0 / 255),
        "2017": (0 / 255, 158 / 255, 115 / 255),
        "2018": (213 / 255, 94 / 255, 0 / 255),}

    marker_map = {
        "2016a": 'o',
        "2016b": 's',
        "2017": '^',
        "2018": 'v',
    }

    def draw_isr_combination_plot(mass_dict, pt_dict,
                                  mass_combined=None, pt_combined=None,
                                  mass_sys_combined=None, pt_sys_combined=None,
                                  channel="", periods=[],
                                  color_map={}, marker_map={},
                                  save_prefix="ISR_Combined"):
        if not periods:
            raise ValueError("Periods list must not be empty!")

        reference_period = f"{periods[0]}_{channel}"
        plotter = pt_dict[reference_period].draw_isr_plot(
            mass_dict[reference_period],
            key='measurement',
            save_and_reset_plotter=False,
            do_fit=False,
            label=reference_period.split("_")[0],
            color=color_map.get(periods[0], 'black'),
            marker=marker_map.get(periods[0], 'o'),
            linestyle='none',
            mfc='none'
        )
        plotter.set_experiment_label(year='Run 2')

        for period in periods[1:]:
            period_key = f"{period}_{channel}"
            pt_dict[period_key].add_isr_plot(
                plotter,
                mass_dict[period_key],
                pt_dict[period_key],
                do_fit=False,
                label=period,
                color=color_map.get(period, 'black'),
                marker=marker_map.get(period, 'o'),
                linestyle='none',
                mfc='none'
            )

        # --- Add Combined if provided ---
        if mass_combined is not None and pt_combined is not None:
            pt_dict[reference_period].add_isr_plot(
                plotter,
                mass_combined,
                pt_combined,
                sys_mean_mass=mass_sys_combined,
                sys_mean_pt=pt_sys_combined,
                do_fit=True,
                label='Combined',
                zorder=10,
                color='black',
                marker='o',
                linestyle='none',
                linewidth=0.7
            )

        plotter.draw_errorbar()
        plotter.show_legend(location=(0, 0))
        plotter.save_and_reset_plotter(f"{save_prefix}_{channel}")

    #draw_isr_combination_plot(mass_dict, pt_dict,
    #                          mass_combined=mass_combined_ee, pt_combined=pt_combined_ee,
    #                          mass_sys_combined=sys_mass_ee_combined_df, pt_sys_combined=sys_pt_ee_combined_df,
    #                          channel="ee", periods=periods,
    #                          color_map=color_map, marker_map=marker_map)

    #draw_isr_combination_plot(mass_dict, pt_dict,
    #                          mass_combined=mass_combined_mm, pt_combined=pt_combined_mm,
    #                          mass_sys_combined=sys_mass_mm_combined_df, pt_sys_combined=sys_pt_mm_combined_df,
    #                          channel="mm", periods=periods,
    #                          color_map=color_map, marker_map=marker_map)





if __name__ == "__main__":
    main()
    sys.exit(0)
