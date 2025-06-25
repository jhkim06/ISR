from ISRAnalyzer import ISRAnalyzer
from ISRCombiner import ISRCombiner
import logging
from Plotter import Plotter
import ROOT
import sys
from ISRLinearFitter import ISRLinearFitter
import numpy as np
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

    setups = [
        {"period": "2017", "channel": "ee", "event_selection": "TightID_b_veto"},
        {"period": "2018", "channel": "ee", "event_selection": "TightID_b_veto"},
        {"period": "2018", "channel": "mm", "event_selection": "TightID_TightIso_b_veto"},
        {"period": "2016a", "channel": "ee", "event_selection": "TightID_b_veto"},
        {"period": "2016b", "channel": "ee", "event_selection": "TightID_b_veto"},
        {"period": "2016a", "channel": "mm", "event_selection": "TightID_TightIso_b_veto"},
        {"period": "2016b", "channel": "mm", "event_selection": "TightID_TightIso_b_veto"},
        {"period": "2017", "channel": "mm", "event_selection": "TightID_TightIso_b_veto"},
    ]

    mass_dict = {}
    pt_dict = {}
    use_2d_pt = True
    sys_on = True
    draw_all_plots = False

    for setup in setups:
        period = setup["period"]
        channel = setup["channel"]
        event_selection = setup["event_selection"]

        reg_mode = 'None'
        analyzer = ISRAnalyzer(sample_base_dir, mass_bins, pt_bins, sys_on=sys_on)
        analyzer.setup_isr_detector_hists(period, channel, event_selection, use_2d_pt=use_2d_pt)
        tau, _ = analyzer.mass_isr_unfold(do_iterative=False, reg_mode=reg_mode, tau_scan_method='scan_lcurve')
        tau_scan_method_for_pt = None
        # Note significant change in chi2 of fit result?
        #if period == '2017' and channel == 'ee':
        #    # apply regularisation only for ee 2017
        #    tau_scan_method_for_pt = 'scan_lcurve' 

        tau_pt, _ = analyzer.pt_isr_unfold(do_iterative=False, reg_mode=reg_mode, tau_scan_method=tau_scan_method_for_pt)
        analyzer.isr_acceptance_corrections()
        pt, mass = analyzer.get_isr_results()
        pt.set_ISRHistSet_per_mass_window()
        mass.set_ISRHistSet_per_mass_window()

        # For FSR systematic
        analyzer_1d = ISRAnalyzer(sample_base_dir, mass_bins, pt_bins, sys_on=False,
                                  unfolded_space_name='gen_dressedp1', acceptance_space_name='dressedp1',
                                  mass_unfolded_space_name='gen_dressedp1', mass_acceptance_space_name='dressedp1')
        analyzer_1d.setup_isr_detector_hists(period, channel, event_selection, use_2d_pt=True)
        analyzer_1d.mass_isr_unfold(do_iterative=False,
                                    reg_mode=reg_mode, tau_scan_method='scan_lcurve')
        analyzer_1d.pt_isr_unfold(do_iterative=False,
                                  reg_mode=reg_mode, tau_scan_method=tau_scan_method_for_pt)
        analyzer_1d.isr_acceptance_corrections()

        pt_1d, mass_1d = analyzer_1d.get_isr_results()
        pt_1d.set_ISRHistSet_per_mass_window()
        mass_1d.set_ISRHistSet_per_mass_window()

        analyzer_fsr_nominal = ISRAnalyzer(sample_base_dir, mass_bins, pt_bins, sys_on=False,
                               folded_space_name='gen_dressedp1',
                               unfolded_space_name='gen_dressed', acceptance_space_name='dressed',
                               pt_folded_bin_name='coarse_O',

                               mass_folded_space_name='gen_dressedp1',
                               mass_unfolded_space_name='gen_dressed', mass_acceptance_space_name='dressed',
                               mass_folded_bin_name='coarse_O',)

        # Use aMC@NLO DY where Pythia is used for QED FSR (acceptace use nominal DY sample)
        analyzer_fsr_pythia = ISRAnalyzer(sample_base_dir, mass_bins, pt_bins, sys_on=False,
                               signal='DY:aMCNLO', acceptance='DY',
                               folded_space_name='gen_dressedp1',
                               unfolded_space_name='gen_dressed', acceptance_space_name='dressed',
                               pt_folded_bin_name='coarse_O',

                               mass_folded_space_name='gen_dressedp1',
                               mass_unfolded_space_name='gen_dressed', mass_acceptance_space_name='dressed',
                               mass_folded_bin_name='coarse_O',)

        analyzer_fsr_nominal.setup_isr_hists(period, channel, event_selection, pt_1d, mass_1d)
        analyzer_fsr_pythia.setup_isr_hists(period, channel, event_selection, pt_1d, mass_1d)

        analyzer_fsr_nominal.pt_isr_unfold()
        analyzer_fsr_nominal.mass_isr_unfold()
        analyzer_fsr_nominal.isr_acceptance_corrections()
        pt_fsr_nominal, mass_fsr_nominal = analyzer_fsr_nominal.get_isr_results()
        pt_fsr_nominal.draw_response_matrix(mass_window_index=0, out_name_postfix='nominalFSR', label_postfix='(PHOTOS)', show_number=True, 
                                            x_axis_label_prefix='pre-FSR', y_axis_label_prefix='post-FSR',
                                            cbarsize='3%', cbarpad=0, norm=mcolors.LogNorm())
        mass_fsr_nominal.draw_response_matrix(mass_window_index=0, out_name_postfix='nominalFSR', label_postfix='(PHOTOS)', show_number=True, 
                                              x_axis_label_prefix='pre-FSR', y_axis_label_prefix='post-FSR', 
                                              cbarsize='3%', cbarpad=0, norm=mcolors.LogNorm())

        pt_fsr_nominal.set_ISRHistSet_per_mass_window()
        mass_fsr_nominal.set_ISRHistSet_per_mass_window()

        analyzer_fsr_pythia.pt_isr_unfold()
        analyzer_fsr_pythia.mass_isr_unfold()
        analyzer_fsr_pythia.isr_acceptance_corrections()
        pt_fsr_pythia, mass_fsr_pythia = analyzer_fsr_pythia.get_isr_results()
        pt_fsr_pythia.draw_response_matrix(mass_window_index=0, out_name_postfix='pythiaFSR', label_postfix='(Pythia)', show_number=True,
                                           x_axis_label_prefix='pre-FSR', y_axis_label_prefix='post-FSR',
                                           cbarsize='3%', cbarpad=0, norm=mcolors.LogNorm())
        mass_fsr_pythia.draw_response_matrix(mass_window_index=0, out_name_postfix='pythiaFSR', label_postfix='(Pythia)', show_number=True,
                                             x_axis_label_prefix='pre-FSR', y_axis_label_prefix='post-FSR',
                                             cbarsize='3%', cbarpad=0, norm=mcolors.LogNorm())

        pt_fsr_pythia.set_ISRHistSet_per_mass_window()
        mass_fsr_pythia.set_ISRHistSet_per_mass_window()

        pt.add_external_hist_as_sys_hist(pt_fsr_nominal, 'FSR', 'nominal')
        mass.add_external_hist_as_sys_hist(mass_fsr_nominal, 'FSR', 'nominal')
        pt.add_external_hist_as_sys_hist(pt_fsr_pythia, 'FSR', 'pythia')
        mass.add_external_hist_as_sys_hist(mass_fsr_pythia, 'FSR', 'pythia')
        pt.update_systematics()
        mass.update_systematics()

        # compare aMC@NLO at parton level
        pt_others = []
        mass_others = []
        sys_on_for_NLO=True
        test_aMCNLO = ISRAnalyzer(sample_base_dir,
                                  mass_bins,
                                  pt_bins, signal="DY:aMCNLO", sys_on=sys_on_for_NLO)  # TODO include systematic!
        test_aMCNLO.setup_isr_acceptance_hists(period, channel, event_selection, is_2d=True)
        pt_aMCNLO, mass_aMCNLO = test_aMCNLO.get_isr_results()
        pt_others.append(pt_aMCNLO)
        mass_others.append(mass_aMCNLO)

        plotter = pt.draw_isr_plot(mass, save_as_csv=True, save_and_reset_plotter=True, do_fit=True, show_chi2=True,
                                   linestyle='none', marker='o', color='black', label='Data', ms=5, zorder=1001, capsize=3)

        plotter = pt.draw_isr_plot(mass, save_as_csv=True, save_and_reset_plotter=False, show_chi2=True,
                                   linestyle='none', marker='o', color='black', label='Data', ms=5, zorder=1001, capsize=3)
        # include aMC@NLO estimation
        key='simulation'
        plotter.add_errorbar((mass_aMCNLO.get_mean_df(key=key), pt_aMCNLO.get_mean_df(key=key)), marker='o', ms=4,
                              mfc='none', fill_between_only=True,
                              color="blue", label='aMC@NLO', linestyle='--', linewidth=0.7)
        plotter.draw_errorbar()
        plotter.set_experiment_label(**{'year': period})
        plotter.show_legend(location=(0, 0), **{"loc": "upper left"})
        plotter.save_and_reset_plotter("isr_"+channel+period+"test_with_nlo")

        # draw plots
        if draw_all_plots:
            for index in range(len(mass_bins)):
                pt.draw_detector_level(index, bin_width_norm=True)
                pt.draw_background_fractions(index)
                pt.draw_unfold_inputs(index, bin_width_norm=True)  # check unfold inputs
                pt.draw_fake_hists(index, bin_width_norm=True)
                pt.draw_unfold_closure(index, bin_width_norm=True)
                pt.draw_unfolded_level(index, bin_width_norm=True, mc_denominator=True, **{"histtype":"errorbar", "marker": ".",
                                          "markersize": 0})
                pt.draw_acceptance_corrected_level(index, bin_width_norm=True, mc_denominator=True,
                                                   **{"histtype":"errorbar", "marker": ".", "markersize": 0})
                other_kwargs = {'histtype': 'errorbar', 'marker': 'o', 'markersize': 6, 'mfc':'none', "color":"skyblue", "mec":"blue",}
                pt.draw_acceptance_corrected_level(mass_window_index=index, bin_width_norm=True, others=pt_others, add_more_hist=True,
                                                   mc_denominator=False, other_kwargs=other_kwargs,
                                                   **{"histtype":"errorbar", "marker": "o", "mfc": "none", "mec":"red",
                                                      "markersize": 5})
                pt.draw_acceptance(mass_window_index=index, bin_width_norm=True, y_max=0.9)
                pt.draw_acceptance(mass_window_index=index, bin_width_norm=True, y_max=1.05, show_only_acceptance_times_efficiency=False)
                if sys_on:
                    pt.draw_systematic_summary(mass_window_index=index)
                    pt.draw_systematic_hists(None, hist_type='acceptance_corrected', key='measurement', mass_window_index=index,  bin_width_norm=True)

            pt.draw_unfold_inputs(-1, bin_width_norm=False)
            pt.draw_detector_level(-1, bin_width_norm=False)
            pt.draw_fake_hists(-1, bin_width_norm=False)
            pt.draw_unfold_closure(-1, bin_width_norm=False)

            mass.draw_background_fractions(0)
            mass.draw_detector_level(0, bin_width_norm=True, y_min_scale=0.1)
            mass.draw_unfolded_level(0, bin_width_norm=True, mc_denominator=True, **{"histtype":"errorbar", "marker": ".",
                                          "markersize": 0})
            mass.draw_response_matrix(mass_window_index=0, show_number=True, cbarsize='3%', cbarpad=0, norm=mcolors.LogNorm())
            mass.draw_acceptance(mass_window_index=0, bin_width_norm=True, y_max=0.9)
            mass.draw_acceptance(mass_window_index=0, bin_width_norm=True, y_max=1.05, show_only_acceptance_times_efficiency=False)
            mass.draw_correlations(mass_window_index=0, cbarsize='3%', cbarpad=0)
            mass.draw_bin_efficiency()
            if use_2d_pt:
                pt.draw_response_matrix(mass_window_index=0, show_number=True, cbarsize='3%', cbarpad=0, norm=mcolors.LogNorm())
                pt.draw_correlations(mass_window_index=0, cbarsize='3%', cbarpad=0)
                pt.draw_bin_efficiency()
            mass.draw_acceptance_corrected_level(0, bin_width_norm=True, mc_denominator=True, **{"histtype":"errorbar", "marker": ".", 
                                                                                                 "markersize": 0})
            other_kwargs = {'histtype': 'errorbar', 'marker': 'o', 'markersize': 6, 'mfc':'none', "color":"skyblue", "mec":"blue",}
            mass.draw_acceptance_corrected_level(mass_window_index=0, bin_width_norm=True, others=mass_others, add_more_hist=True,
                                                 mc_denominator=False, other_kwargs=other_kwargs,
                                                 **{"histtype":"errorbar", "marker": "o", "mfc": "none", "mec":"red",
                                                  "markersize": 5})
            if sys_on:
                mass.draw_systematic_summary(mass_window_index=0)
                mass.draw_systematic_hists(None, hist_type='acceptance_corrected', key='measurement', mass_window_index=0,  bin_width_norm=True)

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

            del ss_test
            del ss_pt
            del ss_mass

        # Save to dict
        key = f"{period}_{channel}"
        mass_dict[key] = mass
        pt_dict[key] = pt

    if draw_all_plots:
        for i_mass_index in range(len(mass_bins)):
            pt_dict["2016a_ee"].draw_pt_comparisons(pt_dict['2016b_ee'], pt_dict['2017_ee'], pt_dict['2018_ee'],
                                                    key='unfold_input',
                                                    index=i_mass_index, scale=-1, bin_width_norm=True)

            pt_dict["2016a_ee"].draw_pt_comparisons(pt_dict['2016b_ee'], pt_dict['2017_ee'], pt_dict['2018_ee'],
                                    key='unfolded_measurement',
                                    index=i_mass_index, scale=-1, bin_width_norm=True)

            pt_dict["2016a_mm"].draw_pt_comparisons(pt_dict['2016b_mm'], pt_dict['2017_mm'], pt_dict['2018_mm'],
                                    key='unfold_input',
                                    index=i_mass_index, scale=-1, bin_width_norm=True)

            pt_dict["2016a_mm"].draw_pt_comparisons(pt_dict['2016b_mm'], pt_dict['2017_mm'], pt_dict['2018_mm'],
                                            key='unfolded_measurement',
                                            index=i_mass_index, scale=-1, bin_width_norm=True)


	# save each systematics
    sys_mass = {}
    sys_pt = {}
    sys_pt['2016a'], sys_mass['2016a'] = pt_dict["2016a_ee"].get_sys_mean_dfs(pt_dict["2016a_ee"], mass_dict["2016a_ee"])
    sys_pt['2016b'], sys_mass['2016b'] = pt_dict["2016b_ee"].get_sys_mean_dfs(pt_dict["2016b_ee"], mass_dict["2016b_ee"])
    sys_pt['2017'], sys_mass['2017'] =   pt_dict["2017_ee"].get_sys_mean_dfs(pt_dict["2017_ee"], mass_dict["2017_ee"])
    sys_pt['2018'], sys_mass['2018'] =   pt_dict["2018_ee"].get_sys_mean_dfs(pt_dict["2018_ee"], mass_dict["2018_ee"])

    for period in ['2016a', '2016b', '2017', '2018']:	
        with open(f'./results/ee{period}_mass_sys.pkl', 'wb') as f:
            pickle.dump(sys_mass[period], f)
        with open(f'./results/ee{period}_pt_sys.pkl', 'wb') as f:
            pickle.dump(sys_pt[period], f)
		
    del sys_pt
    del sys_mass

    sys_mass = {}
    sys_pt = {}
    sys_pt['2016a'], sys_mass['2016a'] = pt_dict["2016a_mm"].get_sys_mean_dfs(pt_dict["2016a_mm"], mass_dict["2016a_mm"])
    sys_pt['2016b'], sys_mass['2016b'] = pt_dict["2016b_mm"].get_sys_mean_dfs(pt_dict["2016b_mm"], mass_dict["2016b_mm"])
    sys_pt['2017'], sys_mass['2017'] =   pt_dict["2017_mm"].get_sys_mean_dfs(pt_dict["2017_mm"], mass_dict["2017_mm"])
    sys_pt['2018'], sys_mass['2018'] =   pt_dict["2018_mm"].get_sys_mean_dfs(pt_dict["2018_mm"], mass_dict["2018_mm"])

    for period in ['2016a', '2016b', '2017', '2018']:	
        with open(f'./results/mm{period}_mass_sys.pkl', 'wb') as f:
            pickle.dump(sys_mass[period], f)
        with open(f'./results/mm{period}_pt_sys.pkl', 'wb') as f:
            pickle.dump(sys_pt[period], f)

    del sys_pt
    del sys_mass


if __name__ == "__main__":
    main()
    print("FINISHED")
    sys.exit(0)
