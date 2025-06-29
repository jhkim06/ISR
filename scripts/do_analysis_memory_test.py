from ISRAnalyzer import ISRAnalyzer
from ISRCombiner import ISRCombiner
import logging
from Plotter import Plotter
import ROOT
import sys
from ISRLinearFitter import ISRLinearFitter
import numpy as np
ROOT.gROOT.SetBatch(True)
import gc


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
        #{"period": "2016b", "channel": "ee", "event_selection": "TightID_b_veto"},
        #{"period": "2017", "channel": "ee", "event_selection": "TightID_b_veto"},
        #{"period": "2018", "channel": "ee", "event_selection": "TightID_b_veto"},
        #{"period": "2016a", "channel": "mm", "event_selection": "TightID_TightIso_b_veto"},
        #{"period": "2016b", "channel": "mm", "event_selection": "TightID_TightIso_b_veto"},
        #{"period": "2017", "channel": "mm", "event_selection": "TightID_TightIso_b_veto"},
        #{"period": "2018", "channel": "mm", "event_selection": "TightID_TightIso_b_veto"},
    ]

    mass_dict = {}
    pt_dict = {}
    use_2d_pt = True
    sys_on = True

    for setup in setups:
        period = setup["period"]
        channel = setup["channel"]
        event_selection = setup["event_selection"]

        reg_mode = 'None'
        analyzer = ISRAnalyzer(sample_base_dir, mass_bins, pt_bins, sys_on=sys_on)
        analyzer.setup_isr_detector_hists(period, channel, event_selection, use_2d_pt=use_2d_pt)
        tau, _ = analyzer.mass_isr_unfold(do_iterative=False, reg_mode=reg_mode, tau_scan_method='scan_lcurve')
        tau_scan_method_for_pt = None
        if period == '2017' and channel == 'ee':
            # apply regularisation only for ee 2017
            tau_scan_method_for_pt = 'scan_lcurve' 
        analyzer.pt_isr_unfold(do_iterative=False, reg_mode=reg_mode, tau_scan_method=tau_scan_method_for_pt)
        analyzer.isr_acceptance_corrections()
        pt, mass = analyzer.get_isr_results()
        pt.set_ISRHistSet_per_mass_window()
        mass.set_ISRHistSet_per_mass_window()

        analyzer_1d = ISRAnalyzer(sample_base_dir, mass_bins, pt_bins, sys_on=False)
        analyzer_1d.setup_isr_detector_hists(period, channel, event_selection, use_2d_pt=False)
        analyzer_1d.mass_isr_unfold(do_iterative=False, tau=tau,  # use the same regularization strength
                                    reg_mode=reg_mode, tau_scan_method=None)
        do_iterative_for_last_window = False
        if tau_scan_method_for_pt == 'scan_lcurve':
            do_iterative_for_last_window = True
        analyzer_1d.pt_isr_unfold(do_iterative=False, tau=0, max_iter=4,
                                  do_iterative_for_last_window=do_iterative_for_last_window)
        analyzer_1d.isr_acceptance_corrections()
        pt_1d, mass_1d = analyzer_1d.get_isr_results()
        pt_1d.set_ISRHistSet_per_mass_window()
        mass_1d.set_ISRHistSet_per_mass_window()

        pt.add_external_hist_as_sys_hist(pt_1d, 'FSR')  # mass use only 1D unfolding
        
        for i_th in range(100):
            print(i_th, " TESTING................")
            #pt.draw_memory_test()
            pt.draw_detector_level(2, bin_width_norm=True)
            pt.draw_acceptance(mass_window_index=2, bin_width_norm=True)
            gc.collect()


if __name__ == "__main__":
    main()
    print("======================END==========================")
    sys.exit(0)
