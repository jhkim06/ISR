import ROOT
import ctypes
from types import MappingProxyType
from Hist import Hist, to_numpy
import numpy as np
from Plotter import Plotter
from array import array
import re
import math


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


def _force_negative_bin_as_zero(matrix):
    n_x_bins = matrix.GetXaxis().GetNbins()
    n_y_bins = matrix.GetYaxis().GetNbins()

    for ix in range(0, n_x_bins + 2):
        for iy in range(0, n_y_bins + 2):
            if matrix.GetBinContent(ix, iy) < 0:
                print("negative bin here:", ix, iy, matrix.GetBinContent(ix, iy))
                matrix.SetBinContent(ix, iy, 0)


class TUnFolder:
    def __init__(self, response_matrix, input_hist,
                 input_fake_hist=None, bg_hists=None,
                 label='default',
                 unfolded_bin=None, folded_bin=None,
                 # TODO allow dynamic unfold setup
                 reg_mode='None', ex_constraint='None', density_mode='BinWidth',
                 tau_scan_method=None, n_iterative=100,

                 analysis_name='',
                 iterative=False,
                 efficiency_correction=True,

                 variable_name=''):
        # TODO Calculate bottom line test value before/after unfolding

        self.year = input_hist.year
        self.channel = input_hist.channel
        self.variable_name = variable_name
        self.iterative = iterative
        self.input_hist_for_iterEM = None
        self.chi2_folded = 0
        self.chi2_unfolded = 0

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
        if self.folded_bin is not None and self.unfolded_bin is not None:
            self.use_tunfoldbinning = True

        self.label = label

        self.reg_mode = reg_mode
        self.ex_constraint = ex_constraint
        self.density_mode = density_mode
        self.tau_scan_method = tau_scan_method
        self.n_iterative = n_iterative
        self.custom_regularization_for_mass = False
        self.custom_regularization_for_pt = False

        self.reg_strength = 0
        self.iter_best = 0

        self.graphSURE = ROOT.TGraph()
        self.df_deviance = ROOT.TGraph()
        self.lcurve = ROOT.TGraph()
        self.curvature = ROOT.TSpline3()
        self.logTaux = ROOT.TSpline3()
        self.logTauy = ROOT.TSpline3()
        self.logTaucurvature = ROOT.TSpline3()

        self.tunfolder = self._create_tunfolder()
        self._set_tunfolder_input()
        self._subtract_backgrounds()

    # comparable regularization to iterative 4
    def apply_custom_regularization_for_mass(self, tunfolder=None):
        if tunfolder is None:
            self.tunfolder.RegularizeCurvature(4, 5, 6)
            self.tunfolder.RegularizeCurvature(5, 6, 7)
            #self.tunfolder.RegularizeCurvature(5, 7, 10)

        else:
            tunfolder.RegularizeCurvature(4, 5, 6)
            tunfolder.RegularizeCurvature(5, 6, 7)
            #tunfolder.RegularizeCurvature(5, 7, 10)

        self.custom_regularization_for_mass = True

    def apply_custom_regularization_for_pt(self, tunfolder=None):
        ndim = self.folded_bin.GetDistributionDimension()
        if ndim == 2:
            # for 2d pt
            if tunfolder is None:
                #self.tunfolder.RegularizeCurvature(51, 52, 53)
                self.tunfolder.RegularizeDerivative(53, 43)
                self.tunfolder.RegularizeDerivative(51, 41)
                self.tunfolder.RegularizeDerivative(52, 42, -2)
                #self.tunfolder.RegularizeCurvature(52, 53, 54)
            else:
                #tunfolder.RegularizeCurvature(51, 52, 53)
                tunfolder.RegularizeDerivative(53, 43)
                tunfolder.RegularizeDerivative(51, 41)
                tunfolder.RegularizeDerivative(52, 42, -2)
                #tunfolder.RegularizeCurvature(52, 53, 54)
        else:
            self.tunfolder.RegularizeCurvature(1, 2, 3)
            self.tunfolder.RegularizeCurvature(2, 3, 4)

        self.custom_regularization_for_pt = True

    def _create_tunfolder(self, sys_name='', var_name='', other_matrix=None):
        if sys_name:
            #print(self.response_matrix.systematic_raw_root_hists)
            if sys_name in self.response_matrix.systematic_raw_root_hists:
                matrix = self.response_matrix.systematic_raw_root_hists[sys_name][var_name]
            else:
                matrix = self.response_matrix.get_raw_hist()
        else:
            if other_matrix is None:
                matrix = self.response_matrix.get_raw_hist()
            else:
                matrix = other_matrix

        tunfolder = None
        if self.use_tunfoldbinning:  # Use bin map
            if self.iterative:
                _force_negative_bin_as_zero(matrix) # TODO check effect of setting negative bin as zero
                tunfolder = ROOT.TUnfoldIterativeEM(matrix, ROOT.TUnfold.kHistMapOutputHoriz,
                                                    self.unfolded_bin, self.folded_bin)
            else:
                _force_negative_bin_as_zero(matrix)
                tunfolder = ROOT.TUnfoldDensity(matrix,
                                                ROOT.TUnfold.kHistMapOutputHoriz,
                                                REG_MODE[self.reg_mode],
                                                EX_CONSTRAINTS[self.ex_constraint],
                                                DENSITY_MODE[self.density_mode],
                                                self.unfolded_bin, self.folded_bin)
        else:
            if self.iterative:
                pass
            else:
                _force_negative_bin_as_zero(matrix)
                tunfolder = ROOT.TUnfoldDensity(matrix,
                                                ROOT.TUnfold.kHistMapOutputHoriz,
                                                REG_MODE[self.reg_mode],
                                                EX_CONSTRAINTS[self.ex_constraint],
                                                DENSITY_MODE[self.density_mode])
        if sys_name:
            if self.custom_regularization_for_mass:
                self.apply_custom_regularization_for_mass(tunfolder)
            if self.custom_regularization_for_pt:
                self.apply_custom_regularization_for_pt(tunfolder)
        return tunfolder

    def _set_tunfolder_input(self, sys_name='', var_name='', sys_tunfolder=None):

        if sys_name:
            sys_input_hist = self.input_hist.systematic_raw_root_hists[sys_name][var_name].Clone("")
            sys_tunfolder.SetInput(sys_input_hist)
            #sys_tunfolder.SetInput(self.input_hist.get_raw_hist())

            if self.input_fake_hist:
                sys_input_fake_hist = self.input_fake_hist.systematic_raw_root_hists[sys_name][var_name]
                sys_tunfolder.SubtractBackground(sys_input_fake_hist, 'input_fake'+sys_name+var_name)
        else:
            self.tunfolder.SetInput(self.input_hist.get_raw_hist())
            if self.iterative:
                self.input_hist_for_iterEM = self.input_hist.get_raw_hist().Clone('input_hist_for_iterativeEM')

            if self.input_fake_hist is not None:
                self.tunfolder.SubtractBackground(self.input_fake_hist.get_raw_hist(),
                                                  'input_fake')
                if self.iterative:
                    self.input_hist_for_iterEM.Add(self.input_fake_hist.get_raw_hist(), -1)

    def _subtract_backgrounds(self, sys_name='', var_name='', sys_tunfolder=None):
        if self.bg_hists is not None:
            for name, hist in self.bg_hists.items():
                if sys_name:
                    sys_bg = hist.systematic_raw_root_hists[sys_name][var_name]
                    sys_tunfolder.SubtractBackground(sys_bg, name)
                else:
                    self.tunfolder.SubtractBackground(hist.get_raw_hist(), name)
                    if self.iterative:
                        self.input_hist_for_iterEM.Add(hist.get_raw_hist(), -1)

    def unfold(self, unfold_method=None,
               tau=0, max_iter=0,
               max_scan_iter=100,
               sys_tunfolder=None):
        # FIXME
        if sys_tunfolder:
            tunfolder = sys_tunfolder
        else:
            tunfolder = self.tunfolder

        if self.tau_scan_method is None:
            if sys_tunfolder is None:
                #tunfolder.RegularizeCurvature(5, 7, 10)  # for 2017 ee, t=0.0001 looks best, for 1D mass unfold
                if self.iterative:
                    tunfolder.DoUnfold(int(max_iter))
                    self.iter_best = max_iter
                else:
                    tunfolder.DoUnfold(tau)
                    self.reg_strength = tau
            else:
                if self.iterative:
                    tunfolder.DoUnfold(int(self.iter_best))
                else:
                    tunfolder.DoUnfold(self.reg_strength)

        elif self.tau_scan_method == 'scan_sure':
            if sys_tunfolder is None:  #
                if self.iterative:
                    self.iter_best = tunfolder.ScanSURE(int(max_scan_iter), self.graphSURE, self.df_deviance)
                else:
                    i_best = tunfolder.ScanSURE(max_scan_iter,
                                                ctypes.c_double(5e-6),
                                                ctypes.c_double(1e-2),
                                                self.graphSURE,
                                                self.df_deviance,
                                                self.lcurve)
                    print("TAU: ", tunfolder.GetTau())
                    self.reg_strength = tunfolder.GetTau()
                # print(self.reg_strength)
            else:
                if self.iterative:
                    tunfolder.DoUnfold(int(self.iter_best))
                else:
                    tunfolder.DoUnfold(self.reg_strength)
        elif self.tau_scan_method == 'scan_lcurve':
            if sys_tunfolder is None:  # or force scan_lcurve
                if self.iterative:
                    pass
                else:
                    i_best = tunfolder.ScanLcurve(max_scan_iter,
                                                  ctypes.c_double(5e-5),
                                                  ctypes.c_double(5e-3),
                                                  self.lcurve,
                                                  self.logTaux,
                                                  self.logTauy, self.logTaucurvature)

                    self.reg_strength = tunfolder.GetTau()
                    print("TAU: ", tunfolder.GetTau())
            else:
                if self.iterative:
                    tunfolder.DoUnfold(int(self.iter_best))
                else:
                    tunfolder.DoUnfold(self.reg_strength)
        else:
            print(unfold_method, ' is not supported unfolding method.')

        return self.reg_strength, self.iter_best

    def create_stat_matrix(self, n_matrix=100):
        # get std for each unfolded result, use it as statistical from matrix
        nx = self.response_matrix.get_raw_hist().GetNbinsX()
        ny = self.response_matrix.get_raw_hist().GetNbinsY()

        default_matrix = self.response_matrix.get_raw_hist()
        stat_matrix_template = self.response_matrix.get_raw_hist().Clone('stat_matrix')
        stat_matrix_template.Reset()

        matrixs = []
        for i in range(n_matrix):
            temp_matrix = stat_matrix_template.Clone(f"stat{i}")
            for i_x in range(nx+2):
                for i_y in range(ny+2):
                    mean = default_matrix.GetBinContent(i_x, i_y)
                    error = default_matrix.GetBinError(i_x, i_y)
                    temp_matrix.SetBinContent(i_x,i_y, ROOT.gRandom.Gaus(mean, error))
                    temp_matrix.SetBinError(i_x, i_y, error)
            matrixs.append(temp_matrix)
        return matrixs

    def get_matrix_stat(self):
        # helper to build the up/down TH1s given a list of per-bin shifts Δ[i]
        def build_up_down(nominal_hist, deltas):
            up = nominal_hist.Clone("matrixUnc_up")
            down = nominal_hist.Clone("matrixUnc_down")
            up.Reset()
            down.Reset()
            nbin = nominal_hist.GetNbinsX()
            for i in range(1, nbin + 1):
                val = nominal_hist.GetBinContent(i)
                d = deltas[i - 1]
                up.SetBinContent(i, val + d)
                down.SetBinContent(i, val - d)
                up.SetBinError(i, nominal_hist.GetBinError(i))
                down.SetBinError(i, nominal_hist.GetBinError(i))
            return {"up": up, "down": down}

        # get the nominal unfolded result
        nominal = self.get_unfolded_hist()

        if not self.iterative:
            # simply take the diagonal of the uncorrelated matrix
            cov2 = self.tunfolder.GetEmatrixSysUncorr("mCovMatrix")
            # Δ_i = sqrt( cov2(i,i) )
            #deltas = [math.sqrt(cov2.GetBinContent(i, i))
            #         for i in range(1, nominal.GetNbinsX() + 1)]

            deltas = [
                math.sqrt(max(0.0, cov2.GetBinContent(i, i)))
                for i in range(1, nominal.GetNbinsX() + 1)
            ]
        else:
            # build an array of |Δ_unfolded| from N toy‐matrices
            default_vals = np.array(Hist(nominal).to_numpy()[0])
            toy_diffs = []
            for idx, matrix in enumerate(self.create_stat_matrix()):
                tf = self._create_tunfolder(other_matrix=matrix)
                tf.SetInput(self.get_input_hist().Clone(f"matrix_stat{idx}"))
                self.unfold(sys_tunfolder=tf)
                out_vals = np.array(Hist(self.get_unfolded_hist(tunfolder=tf)).to_numpy()[0])
                toy_diffs.append(np.abs(out_vals - default_vals))
            # per‐bin pop-std
            deltas = np.std(toy_diffs, axis=0, ddof=0).tolist()

        # now build & return the two TH1s
        return build_up_down(nominal, deltas)

    def sys_unfold(self, return_input_sys=False, sys_names_to_skip=[]):
        # loop over hist systematics
        sys_unfolded_hist = {}
        sys_unfold_input_hist = {}
        sys_hist = self.input_hist.systematic_raw_root_hists
        for sys_name, variations in sys_hist.items():
            if sys_name in sys_names_to_skip:
                continue
            # print('systematic unfolding... ', sys_name)
            sys_unfolded_hist[sys_name] = {}
            if return_input_sys:
                sys_unfold_input_hist[sys_name] = {}
            for var_name, hist in variations.items():
                tunfolder = self._create_tunfolder(sys_name=sys_name, var_name=var_name)  # create tunfolder with the matrix FIXME allow response matrix variation!
                self._set_tunfolder_input(sys_name=sys_name, var_name=var_name, sys_tunfolder=tunfolder)
                self._subtract_backgrounds(sys_name=sys_name, var_name=var_name, sys_tunfolder=tunfolder)
                self.unfold(sys_tunfolder=tunfolder)
                sys_unfolded_hist[sys_name][var_name] = self.get_unfolded_hist(use_axis_binning=False,
                                                                               tunfolder=tunfolder)
                if return_input_sys and not self.iterative:
                    # FIXME TUnfoldIterative doesn't have GetInput()
                    sys_unfold_input_hist[sys_name][var_name] = tunfolder.GetInput("unfold_input",  # histogram title
                                                                                   ctypes.c_char_p(0),
                                                                                   ctypes.c_char_p(0),
                                                                                   ctypes.c_char_p(0),
                                                                                   False)
        if return_input_sys and not self.iterative:
            return sys_unfolded_hist, sys_unfold_input_hist
        else:
            return sys_unfolded_hist, None

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

            # To avoid RuntimeWarning: Attempt to multiply histograms with different bin limits
            unfolded_hist = self.unfolded_bin.ExtractHistogram("unfolded_hist_extracted",
                                                               unfolded_hist,
                                                               0,  # error matrix
                                                               False,
                                                               projection_mode)
        else:
            unfolded_hist = tunfolder.GetOutput("unfolded_hist",
                                                ctypes.c_char_p(0),
                                                ctypes.c_char_p(0),
                                                ctypes.c_char_p(0),
                                                True)  # for 1D, False/True makes no difference?

        return unfolded_hist

    def get_input_hist(self, projection_mode="*[*]", use_axis_binning=True):
        if self.iterative:
            data_hist = self.input_hist_for_iterEM
        else:
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

    def get_folded_hist(self,):
        return self.tunfolder.GetFoldedOutput("",
                                              ctypes.c_char_p(0),
                                              ctypes.c_char_p(0),
                                              ctypes.c_char_p(0),
                                              False)

    def get_arr_without_uoflow(self, arr, folded_bin=True):

        tunfold_bin = self.folded_bin
        if not folded_bin:
            tunfold_bin = self.unfolded_bin

        arr_no_uoflow = []
        for i in range(len(arr)):
            if 'ufl' in str(tunfold_bin.GetBinName(i + 1)) or 'ofl' in str(tunfold_bin.GetBinName(i + 1)):
                continue
            arr_no_uoflow.append(arr[i])
        arr_no_uoflow = np.array(arr_no_uoflow)
        return arr_no_uoflow

    def get_unfold_output_input_ibin_shifted(self, ibin, up=True):
        shifted_input = self.get_input_hist().Clone("shifted_input")
        stat_err = shifted_input.GetBinError(ibin)
        if not up:
            stat_err = -1. * stat_err
        shifted_input.SetBinContent(ibin, shifted_input.GetBinContent(ibin)+stat_err)

        tunfolder = self._create_tunfolder()
        tunfolder.SetInput(shifted_input)
        self.unfold(sys_tunfolder=tunfolder)

        return self.get_unfolded_hist(tunfolder=tunfolder)

    def get_ematrix_input(self):
        # TUnfoldDensity have its own GetEmatrixInput()
        # but TUnfoldIterative doesn't
        # for TUnfoldIterative, get error matrix by
        error_matrix = ROOT.TUnfoldBinning.CreateHistogramOfMigrations(self.unfolded_bin, self.unfolded_bin, "ematrix")
        unfolded_hist = self.get_unfolded_hist()

        n_bins = unfolded_hist.GetNbinsX()
        for i in range(n_bins+1):
            up = self.get_unfold_output_input_ibin_shifted(i+1)
            down = self.get_unfold_output_input_ibin_shifted(i+1, up=False)
            for a in range(n_bins):
                da = 0.5 * (up.GetBinContent(a+1)-down.GetBinContent(a+1))
                for b in range(n_bins):
                    db = 0.5 * (up.GetBinContent(b+1)-down.GetBinContent(b+1))
                    cov = error_matrix.GetBinContent(a+1, b+1) + da*db
                    error_matrix.SetBinContent(a+1, b+1, cov)
        return error_matrix

    # TODO understand what is required by this chi2 test
    def get_chi2(self, folded=True, projection_mode="*[*]", use_axis_binning=True):
        # FIXME for 1D case
        if folded:
            data_hist = self.get_input_hist(projection_mode, use_axis_binning)
            expectation_hist = self.get_mc_reco_from_response_matrix(projection_mode, use_axis_binning)

            data_values = Hist(data_hist).to_numpy()[0]
            expectation_values = Hist(expectation_hist).to_numpy()[0]
            data_errors = Hist(data_hist).to_numpy()[2]

            data_values = self.get_arr_without_uoflow(data_values)
            expectation_values = self.get_arr_without_uoflow(expectation_values)
            data_errors = self.get_arr_without_uoflow(data_errors)

            # TODO maybe need to handle non-diagonal errors
            chi2 = np.sum(np.square((data_values-expectation_values) / data_errors))
            print("folded chi2: ", chi2)
            return chi2
        else:
            data_hist = self.get_unfolded_hist(projection_mode, use_axis_binning)
            expectation_hist = self.get_mc_truth_from_response_matrix()

            data_values = Hist(data_hist).to_numpy()[0]
            expectation_values = Hist(expectation_hist).to_numpy()[0]
            data_errors = Hist(data_hist).to_numpy()[2]

            data_values = self.get_arr_without_uoflow(data_values, folded_bin=False)
            expectation_values = self.get_arr_without_uoflow(expectation_values, folded_bin=False)
            data_errors = self.get_arr_without_uoflow(data_errors, folded_bin=False)

            # FIXME IterativeEM doesn't have GetEmatrixInput, use custom method to calculate it?
            # remove underflow and overflow bins (to avoid zero)
            if not self.iterative:
                mCovInput = self.tunfolder.GetEmatrixInput("mCovInput")
                input_stat_error_list = []
                for i_x in range(mCovInput.GetNbinsX()):
                    if 'ufl' in str(self.unfolded_bin.GetBinName(i_x+1)) or 'ofl' in str(self.unfolded_bin.GetBinName(i_x+1)):
                        continue
                    reco_list = []
                    for i_y in range(mCovInput.GetNbinsY()):
                        if 'ufl' in str(self.unfolded_bin.GetBinName(i_y+1)) or 'ofl' in str(self.unfolded_bin.GetBinName(i_y+1)):
                            continue
                        reco_list.append(mCovInput.GetBinContent(i_x + 1, i_y + 1))
                    input_stat_error_list.append(reco_list)

                input_stat_error_np = np.array(input_stat_error_list)
                A_inv = np.linalg.inv(input_stat_error_np)

                result1 = A_inv.dot(data_values-expectation_values)
                chi2 = (data_values-expectation_values).dot(result1)
                print("unfolded chi2: ", chi2)
            else:
                chi2 = np.sum(np.square((data_values-expectation_values) / data_errors))
                #print("data-exp ", data_values-expectation_values)
                #print("error ", data_errors)
                print("unfolded chi2: ", chi2)
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

    # FIXME use tunfoldbinning to get global bin index
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
            # xval = xa.GetBinCenter(ix)  # FIXME use TUnfoldBinning
            for iy in range(1, raw_matrix.GetNbinsY() + 1):
                c = raw_matrix.GetBinContent(ix, iy)
                e = raw_matrix.GetBinError(ix, iy)
                if c == 0 and e == 0:
                    continue

                # find which new Y‐bin this yval falls into
                #yval = ya.GetBinCenter(iy)
                #new_iy = hnew.GetYaxis().FindBin(yval)

                # get iy return new_iy
                new_iy = self.get_bin_index_of_x_axis(iy)

                # sum the content and propagate error in quadrature
                old_c = hnew.GetBinContent(ix, new_iy)
                old_e = hnew.GetBinError(ix, new_iy)
                hnew.SetBinContent(ix, new_iy, old_c + c)
                hnew.SetBinError(ix, new_iy, (old_e ** 2 + e ** 2) ** 0.5)
        return hnew

    def get_bin_index_of_x_axis(self, iy):
        def get_bin_center(rng_str, underflow, overflow):
            """Given 'ufl','ofl' or 'low,high', return the correct bin‐center."""
            if rng_str == 'ufl':
                return underflow
            if rng_str == 'ofl':
                return overflow
            left, right = map(float, rng_str.split(',', 1))
            return 0.5 * (left + right)

        # how many axes are we dealing with?
        dim = self.folded_bin.GetDistributionDimension()

        # primary axis edges & under/overflow for axis 0
        edges0 = list(self.folded_bin.GetDistributionBinning(0))
        under0 = edges0[0] - 0.5
        over0 = edges0[-1] + 0.5

        # grab the raw bin-name (e.g. "…:axis0[0,4]" or "…:axis0[0,4]:axis1[5,7]")
        name = str(self.folded_bin.GetBinName(iy))
        tail = name.split(':', 1)[1]  # drop everything before the first ':'
        ranges = re.findall(r'\[([^\]]*)\]', tail)

        if dim == 1:
            # only one axis → take the *first* range
            axis0 = get_bin_center(ranges[0], under0, over0)
            # 1D GetGlobalBinNumber takes just the single axis coordinate
            return self.unfolded_bin.GetGlobalBinNumber(axis0)

        # otherwise dim >= 2, so also do axis1
        edges1 = list(self.folded_bin.GetDistributionBinning(1))
        under1 = edges1[0] - 0.5
        over1 = edges1[-1] + 0.5

        # ranges[0] ↦ axis1 string, ranges[1] ↦ axis0 string
        axis1 = get_bin_center(ranges[0], under1, over1)
        axis0 = get_bin_center(ranges[1], under0, over0)

        # 2D GetGlobalBinNumber takes (x0, x1)
        return self.unfolded_bin.GetGlobalBinNumber(axis0, axis1)

    def get_purity_hist(self):
        rebinned_matrix = self.rebin_y_to_match_x()
        # FIXME
        hist_purity = rebinned_matrix.ProjectionX("hist_purity", 0, -1, "")
        if self.folded_bin.GetDistributionDimension() == 1:
            hist_purity = self.unfolded_bin.ExtractHistogram("hist_purity_",
                                                             hist_purity,
                                                             0,  # error matrix
                                                             True)
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

    def draw_bin_efficiency(self, x_log=True, save_and_reset_plotter=True):
        ROOT.TH1.SetDefaultSumw2()  # FIXME
        prob_matrix = self.tunfolder.GetProbabilityMatrix("histProb", ";P_T(gen);P_T(rec)")
        bin_efficiency_hist = Hist(prob_matrix.ProjectionX("histEfficiency"), "bin_efficiency")
        bin_purity_hist = self.get_purity_hist()
        # bin_purity_hist = Hist(prob_matrix.ProjectionY("histEfficiency"), "bin_purity")

        plotter = bin_efficiency_hist.plotter
        plotter.init_plotter(rows=1, cols=1)
        plotter.set_experiment_label(label='Simulation', year=self.input_hist.year)
        plotter.add_hist(bin_efficiency_hist, as_denominator=False, label='Efficiency', yerr=False, show_err_band=False)
        plotter.add_hist(bin_purity_hist, as_denominator=False, label='Purity', yerr=False, show_err_band=False)
        plotter.draw_hist()

        plotter.get_axis(location=(0, 0)).set_ylim(0., 1.05)
        plotter.get_axis(location=(0, 0)).grid(True, which='both', axis='x', linestyle='--', linewidth=0.7)
        plotter.get_axis(location=(0, 0)).grid(True, which='major', axis='y', linestyle='--', linewidth=0.7)

        # for 2 dimension show bin info
        if self.folded_bin.GetDistributionDimension() == 2:
            nbin = bin_efficiency_hist.raw_root_hist.GetNbinsX()
            vline = []
            for i in range(nbin):
                bin_name = str(self.unfolded_bin.GetBinName(i + 1))
                after_first_colon = bin_name.split(":", 1)[1]
                results = re.findall(r'\[([^\]]*)\]', after_first_colon)
                if results[1] == 'ofl':
                    vline.append(i + 1.5)
                    vline.append(i + 2.5)
            plotter.draw_vlines(vline, color='black', linewidth=0.7)

        if x_log:
            plotter.get_axis(location=(0, 0)).set_xscale("log")

        if save_and_reset_plotter:
            plotter.save_and_reset_plotter("bin_efficiency")
            return None
        else:
            return plotter

    def draw_response_matrix(self, out_name=''):
        plotter = Plotter('CMS',
                          '/Users/junhokim/Work/cms_snu/ISR/Plots')
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

    def get_correlation_matrix(self, useAxisBinning=False):

        return self.tunfolder.GetRhoIJtotal("correlation_matrix",
                                            ctypes.c_char_p(0), ctypes.c_char_p(0), ctypes.c_char_p(0),
                                            useAxisBinning)

    def get_response_matrix(self, useAxisBinning=True):
        return Hist(self.tunfolder.GetProbabilityMatrix(ctypes.c_char_p(0),
                                                        ctypes.c_char_p(0),
                                                        useAxisBinning)), self.condition_number()

    # FIXME
    # general comparison template?
    def bottom_line_test(self, projection_mode="*[*]", use_axis_binning=True, draw_plot=False,
                         out_name=''):

        folded_chi2 = self.get_chi2(folded=True, projection_mode=projection_mode, use_axis_binning=use_axis_binning)
        unfolded_chi2 = self.get_chi2(folded=False, projection_mode=projection_mode, use_axis_binning=use_axis_binning)

        folded_hist = self.get_input_hist(projection_mode, use_axis_binning)
        folded_expectation_hist = self.get_mc_reco_from_response_matrix(projection_mode, use_axis_binning)
        unfolded_hist = self.get_unfolded_hist(projection_mode, use_axis_binning)
        unfolded_expectation_hist = self.get_mc_truth_from_response_matrix()

        folded_hist.Scale(1, "width")
        folded_expectation_hist.Scale(1, "width")
        unfolded_hist.Scale(1, "width")
        unfolded_expectation_hist.Scale(1, "width")

        folded_hist = Hist(folded_hist)
        folded_expectation_hist = Hist(folded_expectation_hist)
        unfolded_hist = Hist(unfolded_hist)
        unfolded_expectation_hist = Hist(unfolded_expectation_hist)

        return folded_chi2 > unfolded_chi2

    def get_mc_truth_from_response_matrix(self, sys_on=False):
        if sys_on:
            # extract all systematic hist
            truth_hist = Hist(self.projection_matrix(self.response_matrix.get_raw_hist()),
                              hist_name=self.response_matrix.hist_name+"_projected_truth",
                              year=self.year, channel=self.channel, label=self.response_matrix.label,)

            for sys_name, variations in self.response_matrix.systematic_raw_root_hists.items():
                if sys_name == 'matrix_model':
                    continue
                for var_name, hist in variations.items():
                    temp_matrix = self.projection_matrix(hist)
                    truth_hist.set_systematic_hist(sys_name, var_name, temp_matrix)
            truth_hist.compute_systematic_rss_per_sysname()
            return truth_hist
        else:
            return self.projection_matrix(self.response_matrix.get_raw_hist())

    def get_mc_reco_from_response_matrix(self):
        return self.projection_matrix(self.response_matrix.get_raw_hist(), project_x=False)

    def projection_matrix(self, raw_2d_hist, project_x=True, sys_on=False):
        if project_x:
            projected_hist = raw_2d_hist.ProjectionX("histMCTruth", 0, -1, "e")
        else:
            projected_hist = raw_2d_hist.ProjectionY("histMCTruth", 0, -1, "e")

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
