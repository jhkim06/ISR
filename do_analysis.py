from ISRAnalyzer import ISRAnalyzer
from ISRCombiner import ISRCombiner
import logging
from Plotter import Plotter
import ROOT
import sys
ROOT.gROOT.SetBatch(True)


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
'''
sample_base_dir = '/Users/junhokim/Work/cms_snu/data/Ultralegacy/'
pt_bins_ = (0.0, 100.0)  # actually pt cut
sim_1d = ISRAnalyzer(sample_base_dir,
                     mass_bins_1d,
                     pt_bins_,)

period = "2018"
channel = "mm"
event_selection = "TightID_TightIso_b_veto"

sim_1d.setup_isr_acceptance_hists(period, channel, event_selection, is_2d=False)
pt_1d, mass_1d = sim_1d.get_isr_results()
key='simulation'

test_nlo = ISRAnalyzer(sample_base_dir,
                       mass_bins_1d,
                       pt_bins_, signal="DY:aMCNLO")
test_nlo.setup_isr_acceptance_hists(period, channel, event_selection, is_2d=False)
pt_nlo_1d, mass_nlo_1d = test_nlo.get_isr_results()

period = "2016a"
channel = "mm"
event_selection = ""

test_lo = ISRAnalyzer(sample_base_dir,
                       mass_bins_1d,
                       pt_bins_, signal="DY:MG")
test_lo.setup_isr_acceptance_hists(period, channel, event_selection, is_2d=False)
pt_lo_1d, mass_lo_1d = test_lo.get_isr_results()
'''

def draw_isr_plot_from_df(mass, pt, save_and_reset_plotter=True, postfix='', **kwargs):

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

        plotter.set_isr_plot_cosmetics(channel='ll',)
        text = r"$p_{T}^{ll}<$ 100 GeV"
        plotter.add_text(text=text, location=(0, 0), do_magic=False, **{"frameon": False, "loc": "upper right", })

        if save_and_reset_plotter:
            plotter.draw_errorbar()
            plotter.show_legend(location=(0, 0), loc='upper left')
            plotter.save_and_reset_plotter("isr_test")
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

    setups = [
        {"period": "2016a", "channel": "ee", "event_selection": "TightID_b_veto"},
        {"period": "2016a", "channel": "mm", "event_selection": "TightID_TightIso_b_veto"},
        {"period": "2016b", "channel": "ee", "event_selection": "TightID_b_veto"},
        {"period": "2016b", "channel": "mm", "event_selection": "TightID_TightIso_b_veto"},
        {"period": "2017", "channel": "ee", "event_selection": "TightID_b_veto"},
        {"period": "2017", "channel": "mm", "event_selection": "TightID_TightIso_b_veto"},
        {"period": "2018", "channel": "ee", "event_selection": "TightID_b_veto"},
        {"period": "2018", "channel": "mm", "event_selection": "TightID_TightIso_b_veto"},
    ]

    mass_dict = {}
    pt_dict = {}
    use_2d_pt = True

    for setup in setups:
        period = setup["period"]
        channel = setup["channel"]
        event_selection = setup["event_selection"]

        reg_mode = 'None'
        analyzer = ISRAnalyzer(sample_base_dir, mass_bins, pt_bins)
        analyzer.setup_isr_detector_hists(period, channel, event_selection, use_2d_pt=use_2d_pt)
        tau, _ = analyzer.mass_isr_unfold(do_iterative=False, reg_mode=reg_mode, tau_scan_method='scan_lcurve')
        analyzer.pt_isr_unfold(do_iterative=False)
        analyzer.isr_acceptance_corrections()
        pt, mass = analyzer.get_isr_results()
        pt.set_ISRHistSet_per_mass_window()
        mass.set_ISRHistSet_per_mass_window()

        analyzer_1d = ISRAnalyzer(sample_base_dir, mass_bins, pt_bins, sys_on=False)
        analyzer_1d.setup_isr_detector_hists(period, channel, event_selection, use_2d_pt=use_2d_pt)
        analyzer_1d.mass_isr_unfold(do_iterative=False, tau=tau,  # use the same regularization strength
                                    reg_mode=reg_mode, tau_scan_method=None)
        analyzer_1d.pt_isr_unfold(do_iterative=False)
        analyzer_1d.isr_acceptance_corrections()
        pt_1d, mass_1d = analyzer_1d.get_isr_results()
        pt_1d.set_ISRHistSet_per_mass_window()
        mass_1d.set_ISRHistSet_per_mass_window()

        pt.add_external_hist_as_sys_hist(pt_1d, '1d_2d')  # mass use only 1D unfolding
        #pt_others = []
        #test_aMCNLO = ISRAnalyzer(sample_base_dir,
        #                        mass_bins,
        #                        pt_bins, signal="DY:aMCNLO", sys_on=False)
        #test_aMCNLO.setup_isr_acceptance_hists(period, channel, event_selection, is_2d=True)
        #pt_aMCNLO, mass_aMCNLO = test_aMCNLO.get_isr_results()
        #pt_others.append(pt_aMCNLO)

        #if period == "2016a":
        #    test_LO = ISRAnalyzer(sample_base_dir,
        #                          mass_bins,
        #                          pt_bins, signal="DY:MG")
        #    test_LO.setup_isr_acceptance_hists(period, channel, "", is_2d=False)
        #    pt_LO, mass_LO = test_LO.get_isr_results()
        #    pt_others.append(pt_LO)
        pt.draw_isr_plot(mass)

        for index in range(len(mass_bins)):
            pt.draw_detector_level(index, bin_width_norm=True)
            pt.draw_background_fractions(index)
            pt.draw_unfold_inputs(index, bin_width_norm=True)  # check unfold inputs
            pt.draw_fake_hists(index, bin_width_norm=True)
            pt.draw_unfold_closure(index, bin_width_norm=True)
            pt.draw_unfolded_level(index, bin_width_norm=True, mc_denominator=False)
            pt.draw_acceptance_corrected_level(index, bin_width_norm=True, mc_denominator=False)
            pt.draw_acceptance(mass_window_index=index, bin_width_norm=True)
            pt.draw_systematic_summary(mass_window_index=index)

        #pt.draw_unfold_inputs(-1, bin_width_norm=False)
        #pt.draw_detector_level(-1, bin_width_norm=False)
        #pt.draw_fake_hists(-1, bin_width_norm=False)
        #pt.draw_unfold_closure(-1, bin_width_norm=False)
        mass.draw_background_fractions(0)
        mass.draw_detector_level(0, bin_width_norm=True)
        mass.draw_unfolded_level(0, bin_width_norm=True, mc_denominator=False)
        mass.draw_response_matrix(mass_window_index=0, cbarsize='3%', cbarpad=0)
        mass.draw_acceptance(mass_window_index=0, bin_width_norm=True)
        mass.draw_correlations(mass_window_index=0, cbarsize='3%', cbarpad=0)
        mass.draw_bin_efficiency()
        if use_2d_pt:
            pt.draw_response_matrix(mass_window_index=0, cbarsize='3%', cbarpad=0)
            pt.draw_correlations(mass_window_index=0, cbarsize='3%', cbarpad=0)
            pt.draw_bin_efficiency()
        mass.draw_acceptance_corrected_level(0, bin_width_norm=True, mc_denominator=False)
        mass.draw_systematic_summary(mass_window_index=0)

        # Same-sign
        ss_test = ISRAnalyzer(sample_base_dir,
                              mass_bins,
                              pt_bins)
        ss_test.background_names = [('top', 'antitop'),
                                    'TTLL', 'GGLL', ('ZZ', 'WZ', 'WW'), 'DYJetsToTauTau_MiNNLO']

        ss_test.setup_isr_detector_hists(period, channel, event_selection, use_2d_pt=True, hist_prefix='ss_')
        ss_pt, ss_mass = ss_test.get_isr_results()
        ss_pt.set_ISRHistSet_per_mass_window()
        ss_mass.set_ISRHistSet_per_mass_window()

        for index in range(len(mass_bins)):
            ss_pt.draw_detector_level(index, bin_width_norm=True)
        ss_mass.draw_detector_level(0, bin_width_norm=True)

        # Save to dict
        key = f"{period}_{channel}"
        mass_dict[key] = mass
        pt_dict[key] = pt

    ## -------------------
    ## Now combine and plot
    ## -------------------
    # Combiner for "ee" channel
    combiner_ee = ISRCombiner()
    # Combiner for "mm" channel
    combiner_mm = ISRCombiner()

    for key in mass_dict:
        period, channel = key.split("_")
        if channel == "ee":
            combiner_ee.get_results_dfs(key, mass_dict[key].get_mean_df(), pt_dict[key].get_mean_df())
        elif channel == "mm":
            combiner_mm.get_results_dfs(key, mass_dict[key].get_mean_df(), pt_dict[key].get_mean_df())

    # TODO delete 

    # Combine all periods for each channel
    mass_combined_ee, pt_combined_ee = combiner_ee.combine()
    mass_combined_mm, pt_combined_mm = combiner_mm.combine()


    # Global Combiner (ee + mm)
    combiner_all = ISRCombiner()

    combiner_all.get_results_dfs("combined_ee", mass_combined_ee, pt_combined_ee)
    combiner_all.get_results_dfs("combined_mm", mass_combined_mm, pt_combined_mm)

    mass_combined_final, pt_combined_final = combiner_all.combine()

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
                                  channel="", periods=[],
                                  color_map={}, marker_map={},
                                  save_prefix="ISR_Combined"):
        """
        Draw ISR measurement plot for given channel and periods.

        Args:
            mass_dict (dict): Mass objects keyed by "period_channel".
            pt_dict (dict): Pt objects keyed by "period_channel".
            mass_combined (object, optional): Combined mass object (default: None).
            pt_combined (object, optional): Combined pt object (default: None).
            channel (str): e.g., "ee" or "mm"
            periods (list): List of period names.
            color_map (dict): Mapping from period to color.
            marker_map (dict): Mapping from period to marker.
            save_prefix (str): Output file prefix.
        """

        if not periods:
            raise ValueError("Periods list must not be empty!")

        reference_period = f"{periods[0]}_{channel}"
        plotter = pt_dict[reference_period].draw_isr_plot(
            mass_dict[reference_period],
            key='measurement',
            save_and_reset_plotter=False,
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
                label='Combined',
                color='black',
                marker='*',
                linestyle='--',
                mfc='none',
                linewidth=0.7
            )

        plotter.draw_errorbar()
        plotter.show_legend(location=(0, 0))
        plotter.save_and_reset_plotter(f"{save_prefix}_{channel}")

    draw_isr_combination_plot(mass_dict, pt_dict,
                              mass_combined=mass_combined_ee, pt_combined=pt_combined_ee,
                              channel="ee", periods=periods,
                              color_map=color_map, marker_map=marker_map)

    draw_isr_combination_plot(mass_dict, pt_dict,
                              mass_combined=mass_combined_mm, pt_combined=pt_combined_mm,
                              channel="mm", periods=periods,
                              color_map=color_map, marker_map=marker_map)

    draw_isr_plot_from_df(mass_combined_final, pt_combined_final,
                          label=r'$ee$ and $\mu\mu$ combined',
                          color='black', marker='o', markersize=3.0, linestyle='--', linewidth=0.7,)

if __name__ == "__main__":
    main()
    sys.exit(0)
