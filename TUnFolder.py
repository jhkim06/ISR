import ROOT
import ctypes
from types import MappingProxyType


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

        unfolded_hist = self.tunfolder.GetOutput("default",
                                                 ctypes.c_char_p(0),
                                                 ctypes.c_char_p(0),
                                                 projection_mode,
                                                 use_axis_binning)
        return unfolded_hist

    def get_chi2(self):
        pass

    def get_condition_number(self):
        pass

    def get_mc_truth_from_response_matrix(self):
        return self.response_matrix.ProjectionX("histMCTruth", 0, -1, "e")

    def get_mc_reco_from_response_matrix(self):
        return self.response_matrix.ProjectionY("histMCReco", 0, -1, "e")
