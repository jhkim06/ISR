import ROOT
import numpy as np


ROOT.gSystem.Load("/Users/junhokim/Work/cms_snu/2.4.0/libBlue.so")

for_value = "%5.2f"
for_uncertainty = for_value
for_weight = "%6.4f"
for_rho = "%5.2f"
for_pull = "%4.2f"
for_chi = "%5.3f"
for_unit = " GeV"

# Define global pointer variables
ROOT.gInterpreter.Declare("TString* EstNam;")
ROOT.gInterpreter.Declare("TString* UncNam;")
ROOT.gInterpreter.Declare("TString* ObsNam;")


class Combiner:
    def __init__(self,
                 observable_name,
                 estimation_info,
                 use_rho_same_channel=True,
                 solve=False):
        self.rho_same_channel = {
            "stat": 0.0,
            "IDsf": 1.0,
            "bg_normalization": 1.0,
            "alternative matrix": 1.0,
            "qcd": 1.0,
            "1d_2d": 1.0,
            "matrix": 1.0,
            "pdf": 1.0,
            "alpha_s": 1.0,
            "matrix_model": 1.0,
            "matrix_stat": 0.0,
            "momentum_scale": 0.0,  #
            "momentum_resolution": 0.0,
            "roccor_scale": 0.0,
            "roccor_resolution": 0.0,
            "roccor_stat": 1.0,
            "btagSF": 1.0,
            "puWeight": 1.0,
            "prefireweight": 1.0,
            "electronIDSF": 1.0,
            "electronRECOSF": 1.0,
            "muonIDSF": 1.0,
            "triggerSF": 1.0,
            "Unfolding": 1.0,
            "Roccor": 1.0,
            "momentum_correction": 1.0,
            "FSR": 1.0,
            "accept_stat": 0.0,
            "scale": 0.7,  # TODO need further study
            "efficiency": 1.0,
        }
        self.rho_same_period = {
            "stat": 0.0,
            "IDsf": 0.0,
            "bg_normalization": 1.0,
            "alternative matrix": 1.0,
            "qcd": 1.0,
            "1d_2d": 1.0,
            "matrix": 1.0,
            "pdf": 1.0,
            "alpha_s": 1.0,
            "matrix_model": 1.0,
            "matrix_stat": 0.0,
            "momentum_scale": 0.0,
            "momentum_resolution": 0.0,
            "roccor_scale": 1.0,
            "roccor_resolution": 1.0,
            "roccor_stat": 1.0,
            "btagSF": 1.0,
            "puWeight": 1.0,
            "prefireweight": 1.0,
            "electronIDSF": 0.0,
            "electronRECOSF": 0.0,
            "triggerSF": 0.0,
            "muonIDSF": 0.0,
            "Unfolding": 1.0,
            "Roccor": 0.0,
            "momentum_correction": 0.0,
            "FSR": 1.0,
            "accept_stat": 0.0,
            "scale": 0.7,
            "efficiency": 0.0,
        }

        self.n_est = len(estimation_info)
        self.n_unc = len(estimation_info[0][1])
        self.n_obs = 1
        self.use_rho_same_channel = use_rho_same_channel

        self.observable_name = observable_name
        self.estimation_info = estimation_info
        self.estimation_values = []
        self.unc_names = []
        self.make_estimation_value_list()
        self.make_unc_name_list()

        self.rho_sou = ROOT.TMatrixD(self.n_est, self.n_est)
        self.blue = ROOT.Blue(self.n_est, self.n_unc)
        self.blue.SetFormat(
            for_value,
            for_uncertainty,
            for_weight,
            for_rho,
            for_pull,
            for_chi,
            for_unit
        )
        if solve:
            self.set_names()
            self.set_inputs()
            self.solve()

    def make_estimation_value_list(self):
        for estimation in self.estimation_info:
            self.estimation_values.append(estimation[0][1])
            for uncertainty in estimation[1]:
                self.estimation_values.append(uncertainty[1])

    def make_unc_name_list(self):
        for uncertainty in self.estimation_info[0][1]:
            self.unc_names.append(uncertainty[0])

    def set_name_for_estimation(self):
        ROOT.gInterpreter.ProcessLine(f"EstNam = new TString[{self.n_est}];")
        for index, estimation in enumerate(self.estimation_info):
            ROOT.gInterpreter.ProcessLine(f"EstNam[{index}] = \"{estimation}\";")

        # self.blue.FillNamEst(ROOT.p_NamEst)
        self.blue.FillNamEst(ROOT.EstNam)

    def set_name_for_uncertainty(self):
        ROOT.gInterpreter.ProcessLine(f"UncNam = new TString[{self.n_unc}];")
        for index, uncertainty in enumerate(self.estimation_info[0][1]):
            ROOT.gInterpreter.ProcessLine(f"UncNam[{index}] = \"{uncertainty}\";")

        self.blue.FillNamUnc(ROOT.UncNam)

    def set_name_for_observable(self):
        ROOT.gInterpreter.ProcessLine(f"ObsNam = new TString[1];")
        ROOT.gInterpreter.ProcessLine(f"ObsNam[0] = \"{self.observable_name}\";")

        self.blue.FillNamObs(ROOT.ObsNam)

    def set_names(self):
        self.set_name_for_estimation()
        self.set_name_for_uncertainty()
        self.set_name_for_observable()

    def fill_estimation(self):
        estimation_index = 0
        for i in range(self.n_est):
            self.blue.FillEst(i,
                              np.array(self.estimation_values[estimation_index:])
                              )
            # print(self.estimation_values[estimation_index:])
            estimation_index += self.n_unc + 1

    def fill_correlation(self):
        for i in range(self.n_unc):
            self.fill_matrix(self.unc_names[i])  # FIXME use name
            self.blue.FillCor(i, self.rho_sou)

    def fill_matrix(self, unc_name):

        if self.use_rho_same_channel:
            rho = self.rho_same_channel
        else:
            rho = self.rho_same_period

        for i in range(self.n_est):
            for j in range(self.n_est):
                if i == j:
                    self.rho_sou[i][j] = 1
                else:
                    self.rho_sou[i][j] = rho[unc_name]

    def set_inputs(self):
        self.fill_estimation()
        self.fill_correlation()

    def solve(self):
        self.blue.FixInp()
        self.blue.Solve()
        #self.blue.LatexResult("test")

    def print(self):
        self.blue.PrintEst()
        self.blue.PrintResult()

    def get_result_list(self):
        result_matrix = ROOT.TMatrixD(
            self.blue.GetActObs(),
            self.blue.GetActUnc()+1
        )
        self.blue.GetResult(result_matrix)
        result_list = []
        combined_result = result_matrix[0][0]
        result_list.append(combined_result)

        unc_list = []
        for unc_index in range(1, self.blue.GetActUnc()+1):
            unc_list.append(result_matrix[0][unc_index])

        # calculate total systematic and total uncertainty and attach it
        total_sys = np.sqrt(np.sum(np.square(unc_list[1:])))
        total_unc = np.sqrt(np.sum(np.square(unc_list)))

        result_list.extend(unc_list)
        result_list.append(total_sys)
        result_list.append(total_unc)

        return result_list
