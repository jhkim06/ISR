from ISRAnalyzer import ISRAnalyzer
sample_base_dir = '/Users/junhokim/Work/cms_snu/data/Ultralegacy/'
mass_bins = [(55.0, 64.0),
             (64.0, 81.0),
             (81.0, 101.0),
             (101.0, 200.0),
             (200.0, 1000.0),]
# mass bins to measure average of dilepton transverse momentum
# (160.0, 700.0), (171.0, 830.0), (185.0, 1000.0),]
pt_bins = (0.0, 100.0)  # actually pt cut

def unfold_and_correct(analyzer, period, channel, event_selection, is_2d=True):
    analyzer.setup_isr_detector_hists(period, channel, event_selection, use_2d_pt=is_2d)
    analyzer.isr_unfolds()
    analyzer.isr_acceptance_corrections()
    return analyzer.get_isr_results()

test = ISRAnalyzer(sample_base_dir,
                   mass_bins,
                   pt_bins, mass_folded_bin_name='fine_O', sys_on=True)

# analysis step
period = "2017"
channel = "ee"
event_selection = "TightID_b_veto"
#event_selection = "TightID_TightIso_b_veto"
test.setup_isr_detector_hists(period, channel, event_selection, use_2d_pt=True)

reg_mode = 'None'
tau, _ = test.mass_isr_unfold(do_iterative=False,
                             reg_mode=reg_mode, tau_scan_method='scan_lcurve')  #
tau_pt, _ = test.pt_isr_unfold(do_iterative=False,
                               reg_mode=reg_mode, tau_scan_method='scan_lcurve')  # without regularisation

test.isr_acceptance_corrections()
pt, mass = test.get_isr_results()

pt.set_ISRHistSet_per_mass_window()
mass.set_ISRHistSet_per_mass_window()

# FSR systematic
analyzer_1d = ISRAnalyzer(sample_base_dir, mass_bins, pt_bins, sys_on=False)
analyzer_1d.setup_isr_detector_hists(period, channel, event_selection, use_2d_pt=False)
analyzer_1d.mass_isr_unfold(do_iterative=False, tau=tau,  # use the same regularization strength
                            reg_mode=reg_mode, tau_scan_method=None)
analyzer_1d.pt_isr_unfold(do_iterative=False, tau=0,
                          reg_mode=reg_mode, tau_scan_method=None)
analyzer_1d.isr_acceptance_corrections()

pt_1d, mass_1d = analyzer_1d.get_isr_results()

pt_1d.set_ISRHistSet_per_mass_window()
mass_1d.set_ISRHistSet_per_mass_window()

# different generator for acceptance correction
#analyzer_acceptance = ISRAnalyzer(sample_base_dir, mass_bins, pt_bins, acceptance="DY:aMCNLO", sys_on=False)
#analyzer_acceptance.setup_isr_detector_hists(period, channel, event_selection, use_2d_pt=True)
#analyzer_acceptance.mass_isr_unfold(do_iterative=False, tau=tau,  # use the same regularization strength
#                            reg_mode=reg_mode, tau_scan_method=None)
#analyzer_acceptance.pt_isr_unfold(do_iterative=False)
#analyzer_acceptance.isr_acceptance_corrections()
#
#pt_acceptance, mass_acceptance = analyzer_acceptance.get_isr_results()
#
#pt_acceptance.set_ISRHistSet_per_mass_window()
#mass_acceptance.set_ISRHistSet_per_mass_window()
#
pt.add_external_hist_as_sys_hist(pt_1d, '1d_2d')
# mass.add_external_hist_as_sys_hist(mass_1d, '1d_2d') # for mass, only 1D unfolding done

#pt.add_external_hist_as_sys_hist(pt_acceptance, 'acc')
#mass.add_external_hist_as_sys_hist(mass_acceptance, 'acc')

##mass.draw_detector_level(bin_width_norm=True,)
#pt.draw_detector_level(mass_window_index=4, bin_width_norm=True)

#mass.draw_unfolded_level(bin_width_norm=True,)
#pt.draw_unfolded_level(mass_window_index=4, bin_width_norm=True)
#
#pt.draw_acceptance_corrected_level(mass_window_index=4, bin_width_norm=True, mc_denominator=False)
#mass.draw_acceptance_corrected_level(bin_width_norm=True, mc_denominator=False)
#
#mass.draw_systematics("pdf", mass_window_index=0, bin_width_norm=True)
#pt.draw_systematics("alpha_s", mass_window_index=4, bin_width_norm=True)
#pt.draw_systematics("pdf", mass_window_index=0, bin_width_norm=True)
###pt.draw_systematics("scale", mass_window_index=0, bin_width_norm=True)
#
#mass.draw_correlations()
#mass.draw_bin_efficiency()
#pt.draw_bin_efficiency(mass_window_index=4)
#
#pt.draw_isr_plot(mass)
