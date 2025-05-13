import ROOT
import ctypes
from types import MappingProxyType
from Hist import Hist
import numpy as np
from Plotter import Plotter
from array import array


REG_MODE = MappingProxyType(
    {"None": ROOT.TUnfold.kRegModeNone,
     "Size": ROOT.TUnfold.kRegModeSize,
     "Derivative": ROOT.TUnfold.kRegModeDerivative,
     "Curvature": ROOT.TUnfold.kRegModeCurvature,
     "Mixed": ROOT.TUnfold.kRegModeMixed
     }
)

DENSITY_MODE = MappingProxyType(
    {"None": ROOT.TUnfoldDensity.kDensityModeNone,
     "BinWidth": ROOT.TUnfoldDensity.kDensityModeBinWidth,
     "User": ROOT.TUnfoldDensity.kDensityModeUser,
     "BinWidthAndUser": ROOT.TUnfoldDensity.kDensityModeBinWidthAndUser
     }
)

EX_CONSTRAINTS = MappingProxyType(
    {"None": ROOT.TUnfold.kEConstraintNone,
     "Area": ROOT.TUnfold.kEConstraintArea
     }
)


class TUnFolder:
    def __init__(self, response_matrix, input_hist,
                 input_fake_hist=None, bg_hists=None,
                 label='default',
                 unfolded_bin=None, folded_bin=None,
                 # TODO allow dynamic unfold setup
                 reg_mode='None', ex_constraint='None', density_mode='BinWidth',
                 unfold_method=None, n_iterative=100,

                 analysis_name='',
                 iterative=False,
                 efficiency_correction=True,

                 variable_name=''):

        self.year = input_hist.year
        self.channel = input_hist.channel
        self.variable_name = variable_name

        # default setup
        self.response_matrix = response_matrix
        self.input_hist = input_hist

        # optional histograms
        self.input_fake_hist = input_fake_hist
        #
        self.bg_hists = bg_hists  # dictionary of background hists

        self.use_tunfoldbinning = False
        self.unfolded_bin = unfolded_bin
        self.folded_bin = folded_bin
        if self.folded_bin is not None and self.unfolded_bin is None:
            self.use_tunfoldbinning = True

        self.label = label

        self.reg_mode = reg_mode
        self.ex_constraint = ex_constraint
        self.density_mode = density_mode
        self.unfold_method = unfold_method
        self.n_iterative = n_iterative

        self.reg_strength = 0

        self.graphSURE = ROOT.TGraph()
        self.df_deviance = ROOT.TGraph()
        self.lcurve = ROOT.TGraph()
        self.curvature = ROOT.TSpline3()

        self.tunfolder = self._create_tunfolder()
        self._set_tunfolder_input()
        self._subtract_backgrounds()

    # FIXME check if use_bin_map really needed
    def _create_tunfolder(self, sys_name='', var_name=''):
        if sys_name:
            #print(self.response_matrix.systematic_raw_root_hists)
            matrix = self.response_matrix.systematic_raw_root_hists[sys_name][var_name]
        else:
            matrix = self.response_matrix.get_raw_hist()
        if self.use_tunfoldbinning:  # Use bin map
            return ROOT.TUnfoldDensity(matrix,
                                       ROOT.TUnfold.kHistMapOutputHoriz,
                                       REG_MODE[self.reg_mode],
                                       EX_CONSTRAINTS[self.ex_constraint],
                                       DENSITY_MODE[self.density_mode],
                                       self.unfolded_bin, self.folded_bin)
        else:
            return ROOT.TUnfoldDensity(matrix,
                                       ROOT.TUnfold.kHistMapOutputHoriz,
                                       REG_MODE[self.reg_mode],
                                       EX_CONSTRAINTS[self.ex_constraint],
                                       DENSITY_MODE[self.density_mode])

    def _set_tunfolder_input(self, sys_name='', var_name='', sys_tunfolder=None):

        if sys_name:
            sys_input_hist = self.input_hist.systematic_raw_root_hists[sys_name][var_name]
            sys_tunfolder.SetInput(sys_input_hist)

            if self.input_fake_hist:
                sys_input_fake_hist = self.input_fake_hist.systematic_raw_root_hists[sys_name][var_name]
                sys_tunfolder.SubtractBackground(sys_input_fake_hist, 'input_fake'+sys_name+var_name)
        else:
            self.tunfolder.SetInput(self.input_hist.get_raw_hist())

            if self.input_fake_hist is not None:
                self.tunfolder.SubtractBackground(self.input_fake_hist.get_raw_hist(),
                                                  'input_fake')

    def _subtract_backgrounds(self, sys_name='', var_name='', sys_tunfolder=None):
        if self.bg_hists is not None:
            for name, hist in self.bg_hists.items():
                if sys_name:
                    sys_bg = hist.systematic_raw_root_hists[sys_name][var_name]
                    sys_tunfolder.SubtractBackground(sys_bg, name)
                else:
                    self.tunfolder.SubtractBackground(hist.get_raw_hist(), name)

    def unfold(self, unfold_method=None,
               n_iterative=100,
               tau=0,
               sys_tunfolder=None):

        # FIXME
        if sys_tunfolder:
            tunfolder = sys_tunfolder
        else:
            tunfolder = self.tunfolder

        if unfold_method is None and self.unfold_method is None:
            self.unfold_method = unfold_method
        else:
            if self.unfold_method is None:
                self.unfold_method = unfold_method

        self.n_iterative = n_iterative

        if self.unfold_method is None:
            if sys_tunfolder is None:
                # tunfolder.RegularizeCurvature(5, 7, 10)  # t=0.0001 looks best, for 1D mass unfold
                tunfolder.RegularizeCurvature(5, 7, 10)
                tunfolder.DoUnfold(tau)  # TODO understand what it does!
            else:
                tunfolder.DoUnfold(tau)
        elif self.unfold_method == 'scan_sure':
            if sys_tunfolder is None:
                i_best = tunfolder.ScanSURE(n_iterative,
                                            ctypes.c_double(0.0),
                                            ctypes.c_double(0.0),
                                            self.graphSURE,
                                            self.df_deviance,
                                            self.lcurve)
                self.reg_strength = pow(10, self.graphSURE.GetX()[i_best])
                # print(self.reg_strength)
            else:
                tunfolder.DoUnfold(0.0)
        else:
            print(unfold_method, ' is not supported unfolding method.')

    def sys_unfold(self):
        # loop over hist systematics
        sys_unfolded_hist = {}
        sys_hist = self.input_hist.systematic_raw_root_hists
        for sys_name, variations in sys_hist.items():
            # print('systematic unfolding... ', sys_name)
            sys_unfolded_hist[sys_name] = {}
            for var_name, hist in variations.items():
                tunfolder = self._create_tunfolder(sys_name=sys_name, var_name=var_name)  # create tunfolder with the matrix FIXME allow response matrix variation!
                self._set_tunfolder_input(sys_name=sys_name, var_name=var_name, sys_tunfolder=tunfolder)
                self._subtract_backgrounds(sys_name=sys_name, var_name=var_name, sys_tunfolder=tunfolder)
                self.unfold(sys_tunfolder=tunfolder)
                sys_unfolded_hist[sys_name][var_name] = self.get_unfolded_hist(use_axis_binning=False,
                                                                               tunfolder=tunfolder)

        # TODO different response matrix?

        return sys_unfolded_hist

    def get_unfolded_hist(self, projection_mode="*[*]",
                          use_axis_binning=True,
                          tunfolder=None):

        if tunfolder is None:
            tunfolder = self.tunfolder
        else:
            tunfolder = tunfolder

        if self.use_tunfoldbinning:
            # it seems parameter axis steering not passed properly so use ExtractHistogram()
            unfolded_hist = tunfolder.GetOutput("unfolded_hist",
                                                ctypes.c_char_p(0),
                                                ctypes.c_char_p(0),
                                                ctypes.c_char_p(0),
                                                False)

            # for "*[*]", use_axis_binning=True, makes no sense?
            #unfolded_hist = self.unfolded_bin.ExtractHistogram("unfolded_hist_extracted",
            #                                                   unfolded_hist,
            #                                                   0,  # error matrix
            #                                                   False,
            #                                                   projection_mode)
        else:
            unfolded_hist = tunfolder.GetOutput("unfolded_hist",
                                                ctypes.c_char_p(0),
                                                ctypes.c_char_p(0),
                                                ctypes.c_char_p(0),
                                                True)  # for 1D, False/True makes no difference?

        return unfolded_hist

    def get_input_hist(self, projection_mode="*[*]", use_axis_binning=True):
        if self.use_tunfoldbinning:
            # it seems parameter axis steering not passed properly so use ExtractHistogram()
            data_hist = self.tunfolder.GetInput("unfold_input",  # histogram title
                                                ctypes.c_char_p(0),
                                                ctypes.c_char_p(0),
                                                ctypes.c_char_p(0),
                                                False)
            #data_hist = self.folded_bin.ExtractHistogram("unfold_input_extracted",
            #                                             data_hist,
            #                                             0,  # error matrix
            #                                             False,
            #                                             projection_mode)
        else:
            data_hist = self.tunfolder.GetInput("unfold_input",  # histogram title
                                                ctypes.c_char_p(0),
                                                ctypes.c_char_p(0),
                                                ctypes.c_char_p(0),
                                                True)

        return data_hist

    # TODO understand what is required by this chi2 test
    def get_chi2(self, folded=True, projection_mode="*[*]", use_axis_binning=True):
        if folded:
            data_hist = self.get_input_hist(projection_mode, use_axis_binning)
            expectation_hist = self.get_mc_reco_from_response_matrix(projection_mode, use_axis_binning)
        else:
            data_hist = self.get_unfolded_hist(projection_mode, use_axis_binning)
            expectation_hist = self.get_mc_truth_from_response_matrix(projection_mode, use_axis_binning)

        # negligible effect on ch2 value?
        data_hist.Scale(1, "width")
        expectation_hist.Scale(1, "width")

        data_values = Hist(data_hist).to_numpy()[0]
        expectation_values = Hist(expectation_hist).to_numpy()[0]
        data_errors = Hist(data_hist).to_numpy()[2]

        # TODO maybe need to handle non-diagonal errors
        chi2 = np.sum(np.square((data_values - expectation_values) / data_errors))

        return chi2

    def condition_number(self, draw_matrix=False):
        h_prob_matrix = self.tunfolder.GetProbabilityMatrix("hProb")
        n_bin_x = h_prob_matrix.GetNbinsX()
        n_bin_y = h_prob_matrix.GetNbinsY()

        matrix = ROOT.TMatrixD(n_bin_y, n_bin_x)
        for i in range(n_bin_x):
            for j in range(n_bin_y):
                matrix[j][i] = h_prob_matrix.GetBinContent(i + 1, j + 1)

        decomp = ROOT.TDecompSVD(matrix)
        print("Matrix condition: " + str(decomp.Condition()))
        return decomp.Condition()

    def rebin_y_to_match_x(self):
        """
        Return a new TH2D whose X‐ and Y‐axes both use h2’s X‐axis bin edges.
        The contents are summed (and errors propagated in quadrature) along Y.
        """
        # 1) Grab the X‐axis bin edges
        raw_matrix = self.response_matrix.get_raw_hist().Clone("rebin_y_axis")
        xa = raw_matrix.GetXaxis()
        nbx = xa.GetNbins()
        # we'll need nbx+1 edges
        edges = array('d', (xa.GetBinLowEdge(i) for i in range(1, nbx + 2)))

        # 2) Create the new histogram
        name = raw_matrix.GetName() + "_rebin"
        title = raw_matrix.GetTitle()
        hnew = ROOT.TH2D(name, title, nbx, edges, nbx, edges)
        # detach from any file
        ROOT.TH1.AddDirectory(False)

        # 3) Loop over all original bins, re‐assign them to new Y‐bins
        ya = raw_matrix.GetYaxis()
        for ix in range(1, raw_matrix.GetNbinsX() + 1):
            xval = xa.GetBinCenter(ix)
            for iy in range(1, raw_matrix.GetNbinsY() + 1):
                yval = ya.GetBinCenter(iy)
                c = raw_matrix.GetBinContent(ix, iy)
                e = raw_matrix.GetBinError(ix, iy)
                if c == 0 and e == 0:
                    continue

                # find which new Y‐bin this yval falls into
                new_iy = hnew.GetYaxis().FindBin(yval)

                # sum the content and propagate error in quadrature
                old_c = hnew.GetBinContent(ix, new_iy)
                old_e = hnew.GetBinError(ix, new_iy)
                hnew.SetBinContent(ix, new_iy, old_c + c)
                hnew.SetBinError(ix, new_iy, (old_e ** 2 + e ** 2) ** 0.5)

        return hnew

    def get_purity_hist(self):
        rebinned_matrix = self.rebin_y_to_match_x()
        hist_purity = rebinned_matrix.ProjectionX("hist_purity", 0, -1, "")
        hist_purity.Reset()

        for iy in range(1, rebinned_matrix.GetNbinsY() + 1):
            sum_ = 0.
            for ix in range(1, rebinned_matrix.GetNbinsX() + 1):
                sum_ += rebinned_matrix.GetBinContent(ix, iy)
            purity = 0.
            if sum_ > 0.:
                purity = rebinned_matrix.GetBinContent(iy, iy) / sum_
            hist_purity.SetBinContent(iy, purity)

        return Hist(hist_purity, 'hist_purity')

    def draw_bin_efficiency(self, x_log=True):
        ROOT.TH1.SetDefaultSumw2()  # FIXME
        prob_matrix = self.tunfolder.GetProbabilityMatrix("histProb", ";P_T(gen);P_T(rec)")
        bin_efficiency_hist = Hist(prob_matrix.ProjectionX("histEfficiency"), "bin_efficiency")
        # extract 1D
        # FIXME update for 2D unfold
        bin_purity_hist = self.get_purity_hist()
        # bin_purity_hist = Hist(prob_matrix.ProjectionY("histEfficiency"), "bin_purity")

        plotter = bin_efficiency_hist.plotter
        plotter.init_plotter(rows=1, cols=1)
        plotter.set_experiment_label(year=self.input_hist.year)
        plotter.add_hist(bin_efficiency_hist, as_denominator=False)
        plotter.add_hist(bin_purity_hist, as_denominator=False)
        plotter.draw_hist()
        plotter.get_axis(location=(0, 0)).set_ylim(0., 1.05)
        plotter.get_axis(location=(0, 0)).grid(True, which='both', axis='x', linestyle='--', linewidth=0.7)
        plotter.get_axis(location=(0, 0)).grid(True, which='major', axis='y', linestyle='--', linewidth=0.7)
        if x_log:
            plotter.get_axis(location=(0, 0)).set_xscale("log")
        plotter.save_and_reset_plotter("bin_efficiency")

    def draw_response_matrix(self, out_name=''):
        plotter = Plotter('CMS',
                          '/Users/junhokim/Work/cms_snu/ISR/Plots')  # FIXME use self.plotter
        plotter.create_subplots(1, 1, figsize=(8, 8),
                                left=0.15, right=0.9, hspace=0.0, bottom=0.15)
        plotter.set_experiment_label(label="Simulation")

        rm_np = Hist(self.tunfolder.GetProbabilityMatrix(ctypes.c_char_p(0),
                                                         ctypes.c_char_p(0),
                                                         True)).to_numpy_2d()

        plotter.draw_matrix(rm_np, "ll")
        plotter.add_text("condition number: " + str(round(self.condition_number(), 2)),
                         location=(0,0), do_magic=False,
                         **{"loc": "upper left",})
        plotter.save_fig(out_name + "_rm_" +  self.channel + self.year)

    def get_response_matrix(self):
        return Hist(self.tunfolder.GetProbabilityMatrix(ctypes.c_char_p(0),
                                                        ctypes.c_char_p(0),
                                                        True))

    # FIXME
    # general comparison template?
    def bottom_line_test(self, projection_mode="*[*]", use_axis_binning=True, draw_plot=False,
                         out_name=''):

        folded_chi2 = self.get_chi2(folded=True, projection_mode=projection_mode, use_axis_binning=use_axis_binning)
        unfolded_chi2 = self.get_chi2(folded=False, projection_mode=projection_mode, use_axis_binning=use_axis_binning)

        folded_hist = self.get_input_hist(projection_mode, use_axis_binning)
        folded_expectation_hist = self.get_mc_reco_from_response_matrix(projection_mode, use_axis_binning)
        unfolded_hist = self.get_unfolded_hist(projection_mode, use_axis_binning)
        unfolded_expectation_hist = self.get_mc_truth_from_response_matrix(projection_mode, use_axis_binning)

        folded_hist.Scale(1, "width")
        folded_expectation_hist.Scale(1, "width")
        unfolded_hist.Scale(1, "width")
        unfolded_expectation_hist.Scale(1, "width")

        folded_hist = Hist(folded_hist)
        folded_expectation_hist = Hist(folded_expectation_hist)
        unfolded_hist = Hist(unfolded_hist)
        unfolded_expectation_hist = Hist(unfolded_expectation_hist)

        if draw_plot:
            # TODO make a generic function to draw comparison plot
            plotter = Plotter('CMS',
                              '/Users/junhokim/Work/cms_snu/ISR/Plots')  # FIXME use self.plotter
            plotter.create_subplots(2, 1, figsize=(8,8),
                                    left=0.15, right=0.95, hspace=0.0, bottom=0.15, height_ratios=[1, 0.3])

            plotter.set_experiment_label(**{"year": self.year})
            # measurement
            # plotter.add_comparison()

            plotter.add_comparison_pair(folded_hist, folded_expectation_hist, location=(0,0), ratio_location=(1,0),
                                        nominator_args={"histtype": 'errorbar',
                                                        "color": 'black', 'label': 'RECO Data'},
                                        denominator_args={"histtype": 'errorbar', 'marker':"s",
                                                          "color": 'red', 'mfc': 'none',
                                                          "label": 'RECO Sim'})
            plotter.draw_hist()

            plotter.add_comparison_pair(unfolded_hist, unfolded_expectation_hist, location=(0,0), ratio_location=(1,0),
                                        nominator_args={"histtype": 'errorbar', "color": 'gray', 'label': 'Unfolded Data'},
                                        denominator_args={"histtype": 'errorbar', 'marker':"s", "color": 'blue', 'mfc': 'none',
                                                          "label": 'GEN Sim'})

            plotter.draw_hist()

            plotter.set_common_comparison_plot_cosmetics('ll')
            plotter.show_legend(location=(0, 0))
            plotter.get_axis(location=(0, 0)).set_xticklabels([])
            plotter.adjust_y_scale()
            plotter.add_text(text='$\chi_{folded}^{2}:$' + str(round(folded_chi2, 2))+'\n$\chi_{unfolded}^{2}:$' +
                                  str(round(unfolded_chi2, 2)),
                             location=(0,0),
                             **{"loc": "upper left",})  # Note: required to after setting legend
            plotter.get_axis(location=(1, 0)).set_ylim(0.4, 1.6)
            plotter.save_fig(out_name + "_btl_test_" + self.channel + self.year)

        return folded_chi2 > unfolded_chi2

    def get_mc_truth_from_response_matrix(self, projection_mode="*[*]", use_axis_binning=True):
        return self.extract_truth_from_2d_hist_using_bin_definition(self.response_matrix.get_raw_hist())

    def extract_truth_from_2d_hist_using_bin_definition(self, raw_2d_hist,):
        projected_hist = raw_2d_hist.ProjectionX("histMCTruth", 0, -1, "e")

        #if projection_mode != "*[*]":
        #projected_hist = self.unfolded_bin.ExtractHistogram("truth_extracted",
        #                                                    projected_hist,
        #                                                    0,  # error matrix
        #                                                    use_axis_binning,
        #                                                    projection_mode)
        return projected_hist

    def get_mc_reco_from_response_matrix(self, projection_mode="*[*]", use_axis_binning=True):
        projected_hist = self.response_matrix.get_raw_hist().ProjectionY("histMCReco", 0, -1, "e")

        if projection_mode != "*[*]":
            projected_hist = self.folded_bin.ExtractHistogram("reco_extracted",
                                                                projected_hist,
                                                                0,  # error matrix
                                                                use_axis_binning,
                                                                projection_mode)
        return projected_hist
