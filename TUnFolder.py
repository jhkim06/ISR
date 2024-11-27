import ROOT
import ctypes
from types import MappingProxyType
from Hist import Hist
import numpy as np
from Plotter import Plotter


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
                 use_bin_map=False,

                 reg_mode='None', ex_constraint='None', density_mode='None',
                 unfold_method=None, n_iterative=100,

                 iterative=False,
                 efficiency_correction=True):

        # default setup
        self.response_matrix = response_matrix
        self.input_hist = input_hist

        # optional histograms
        self.input_fake_hist = input_fake_hist
        self.bg_hists = bg_hists  # dictionary of background hists

        self.unfolded_bin = unfolded_bin
        self.folded_bin = folded_bin
        self.use_bin_map = use_bin_map

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

        self._create_tunfolder()
        self._set_tunfolder_input()
        self._subtract_backgrounds()

    def _create_tunfolder(self):
        if self.use_bin_map:  # Use bin map
            self.tunfolder = ROOT.TUnfoldDensity(self.response_matrix,
                                                 ROOT.TUnfold.kHistMapOutputHoriz,
                                                 REG_MODE[self.reg_mode],
                                                 EX_CONSTRAINTS[self.ex_constraint],
                                                 DENSITY_MODE[self.density_mode],
                                                 self.unfolded_bin, self.folded_bin)
        else:
            self.tunfolder = ROOT.TUnfoldDensity(self.response_matrix,
                                                 ROOT.TUnfold.kHistMapOutputHoriz,
                                                 REG_MODE[self.reg_mode],
                                                 EX_CONSTRAINTS[self.ex_constraint],
                                                 DENSITY_MODE[self.density_mode])

    def _set_tunfolder_input(self):
        self.tunfolder.SetInput(self.input_hist)

        if self.input_fake_hist is not None:
            self.tunfolder.SubtractBackground(self.input_fake_hist, 'input_fake')

    def _subtract_backgrounds(self):
        if self.bg_hists is not None:
            for name, hist in self.bg_hists.items():
                self.tunfolder.SubtractBackground(hist, name)

    def unfold(self, unfold_method=None, n_iterative=100):
        # FIXME
        if unfold_method is None and self.unfold_method is None:
            self.unfold_method = unfold_method
        else:
            if self.unfold_method is None:
                self.unfold_method = unfold_method

        self.n_iterative = n_iterative

        if self.unfold_method is None:
            self.tunfolder.DoUnfold(self.reg_strength)
            self.reg_strength = self.tunfolder.GetTau()
        elif self.unfold_method == 'scan_sure':
            i_best = self.tunfolder.ScanSURE(n_iterative,
                                             ctypes.c_double(0),
                                             ctypes.c_double(0),
                                             self.graphSURE,
                                             self.df_deviance,
                                             self.lcurve)

            self.reg_strength = pow(10, self.graphSURE.GetX()[i_best])
        else:
            print(unfold_method, ' is not supported unfolding method.')

    def get_unfolded_hist(self, projection_mode="*[*]", use_axis_binning=True):

        unfolded_hist = self.tunfolder.GetOutput("unfolded_hist",
                                                 ctypes.c_char_p(0),
                                                 ctypes.c_char_p(0),
                                                 ctypes.c_char_p(0),
                                                 use_axis_binning)

        if projection_mode != "*[*]":
            unfolded_hist = self.unfolded_bin.ExtractHistogram("unfolded_hist_extracted",
                                                             unfolded_hist,
                                                             0,  # error matrix
                                                             use_axis_binning,
                                                             projection_mode)
        return unfolded_hist

    def get_input_hist(self, projection_mode="*[*]", use_axis_binning=True):
        # seems parameter axis steering not passed properly so use ExtractHistogram()
        data_hist = self.tunfolder.GetInput("unfold_input",  # histogram title
                                            ctypes.c_char_p(0),
                                            ctypes.c_char_p(0),
                                            ctypes.c_char_p(0),
                                            use_axis_binning)

        if projection_mode != "*[*]":
            data_hist = self.folded_bin.ExtractHistogram("unfold_input_extracted",
                                                         data_hist,
                                                         0,  # error matrix
                                                         use_axis_binning,
                                                         projection_mode)
        return data_hist

    def get_chi2(self, folded=True, projection_mode="*[*]", use_axis_binning=True):
        if folded:
            data_hist = self.get_input_hist(projection_mode, use_axis_binning)
            expectation_hist = self.get_mc_reco_from_response_matrix(projection_mode, use_axis_binning)
        else:
            data_hist = self.get_unfolded_hist(projection_mode, use_axis_binning)
            expectation_hist = self.get_mc_truth_from_response_matrix(projection_mode, use_axis_binning)

        # histogram width doesnt matter
        # data_hist.Scale(1, "width")
        # expectation_hist.Scale(1, "width")

        data_values = Hist(data_hist).to_numpy()[0]
        expectation_values = Hist(expectation_hist).to_numpy()[0]
        data_errors = Hist(data_hist).to_numpy()[2]

        chi2 = np.sum(np.square((data_values - expectation_values) / data_errors))  # TODO maybe need to handle non-diagonal errors
        # draw comparison plot
        return chi2

    def bottom_line_test(self, projection_mode="*[*]", use_axis_binning=True):
        # required chi2_folded > chi2_unfolded
        # case1: 1D, case2: 2D
        # projection_mode ex) dipt[O];dimass[OUC1]
        folded_chi2 = self.get_chi2(folded=True, projection_mode=projection_mode, use_axis_binning=use_axis_binning)
        unfolded_chi2 = self.get_chi2(folded=False, projection_mode=projection_mode, use_axis_binning=use_axis_binning)
        return folded_chi2 > unfolded_chi2

    def get_condition_number(self):
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


    def draw_bottom_line_test(self, data, expectation, text):

        plotter = Plotter('CMS', './Plots')  # FIXME use self.plotter
        plotter.create_subplots(2, 1, figsize=(8,8),
                                left=0.15, right=0.95, hspace=0.0, bottom=0.15, height_ratios=[1, 0.3])
        # measurement
        # plotter.add_comparison()
        plotter.add_hist(data, **{"histtype": 'errorbar', "color": 'black', 'label': 'Data'})
        # expectations
        plotter.add_hist(expectation, **{"histtype": 'errorbar', "color": 'red', 'mfc': 'none',
                                      "label": 'Simulation'})

        plotter.add_ratio_hist(nominator_index=0,
                               denominator_index=1,  # add all simulation except data
                               location=(1, 0), color='black', histtype='errorbar')
        # plotter.cosmetic()

        plotter.draw_hist()
        plotter.show_legend(location=(0, 0))
        plotter.get_axis(location=(0, 0)).set_xticklabels([])
        plotter.adjust_y_scale()
        plotter.add_text(text=text, location=(0,0), **{"loc": "upper left",})  # Note: required to after setting legend
        plotter.get_axis(location=(1, 0)).set_ylim(0.4, 1.6)
        plotter.save_fig("comparison_test")

    def get_mc_truth_from_response_matrix(self, projection_mode="*[*]", use_axis_binning=True):
        projected_hist = self.response_matrix.ProjectionX("histMCTruth", 0, -1, "e")

        if projection_mode != "*[*]":
            projected_hist = self.unfolded_bin.ExtractHistogram("truth_extracted",
                                                                projected_hist,
                                                                0,  # error matrix
                                                                use_axis_binning,
                                                                projection_mode)
        return projected_hist

    def get_mc_reco_from_response_matrix(self, projection_mode="*[*]", use_axis_binning=True):
        projected_hist = self.response_matrix.ProjectionY("histMCReco", 0, -1, "e")

        if projection_mode != "*[*]":
            projected_hist = self.folded_bin.ExtractHistogram("reco_extracted",
                                                                projected_hist,
                                                                0,  # error matrix
                                                                use_axis_binning,
                                                                projection_mode)
        return projected_hist
