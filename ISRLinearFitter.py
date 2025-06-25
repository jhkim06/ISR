from ROOT import TF1, TGraphErrors, TGraphMultiErrors
from array import array
import numpy as np
import pandas as pd


class ISRLinearFitter:
    def __init__(self, mass_mean, pt_mean, sys_mass_mean=None, sys_pt_mean=None):
        self.mass_mean = mass_mean
        self.pt_mean = pt_mean

        self.fit_funct_nominal = "2.* [0] * TMath::Log(x) + [1]"
        self.fit_funct_sys = "2.* [0] * TMath::Log(x) + [1]"
        self.x_min = 50
        self.x_max = 400

        # nominal fit
        self.funct_nominal = TF1("funct_nominal", self.fit_funct_nominal, self.x_min, self.x_max)
        self.slope = 0
        self.intercept = 0
        self.slope_err = 0
        self.intercept_err = 0
        self.chi2 = 0
        self.ndof = 0

        # for systematic fit
        self.funct_sys = TF1("funct_sys", self.fit_funct_sys, self.x_min, self.x_max)
        self.sys_mean_mass = sys_mass_mean
        self.sys_mean_pt = sys_pt_mean
        self.sys_slope = {}  # "sys_name": [slope1, slope2,...]
        self.sys_intercept = {}
        self.slope_sys = 0
        self.intercept_sys = 0

    def set_tgrapherrors(self, mass_mean, pt_mean):
        graph = TGraphErrors(len(mass_mean['mean']))
        for index in range(len(pt_mean['mean'])):
            graph.SetPoint(index, mass_mean['mean'][index], pt_mean['mean'][index])
            graph.SetPointError(index, mass_mean['stat'][index], pt_mean['stat'][index])
        return graph

    def compute_covariance_matrix(self):
        df = self.pt_mean.drop(columns=['mean', 'sys', 'total_error'])

        # stat_cols = [c for c in df.columns if c == 'stat' or '_stat' in c]
        # these selected systematic sources have low sign-correlation, so treat them in diagonal
        stat_cols = [c for c in df.columns if c == 'stat' or '_stat' in c or
                     c in ['triggerSF', 'electronIDSF', 'electronRECOSF', 'muonIDSF', 'pdf', 'FSR',
                           'efficiency']]
        syst_cols = [c for c in df.columns
                     if c not in stat_cols]

        stat_var = df[stat_cols].pow(2).sum(axis=1).to_numpy()
        cov = np.diag(stat_var)

        for k in syst_cols:
            s = df[k].to_numpy()
            cov += np.outer(s, s)

        return pd.DataFrame(cov, index=df.index, columns=df.index)

    def calculate_chi2(self, delta, cov):
        """
        Compute chi2 = delta^T V^{-1} delta.

        Parameters
        ----------
        delta : np.ndarray
            1D array of length N giving y - f for each point.
        cov : pd.DataFrame
            N×N covariance matrix (as returned by compute_covariance_matrix).

        Returns
        -------
        float
            The chi-squared value.
        """
        # Extract the raw numpy array
        V = cov.values

        # For numerical stability, solve V x = delta rather than explicitly invert
        # i.e. chi2 = delta^T @ V^{-1} @ delta = delta^T @ x  where x = V^{-1} @ delta
        x = np.linalg.solve(V, delta)

        chi2 = float(delta @ x)
        return chi2

    def custom_chi2(self):
        cov_matrix = self.compute_covariance_matrix()
        print(cov_matrix)

        diff = []
        for index in range(len(self.mass_mean['mean'])):
            diff.append(self.pt_mean['mean'][index]-self.funct_nominal.Eval(self.mass_mean['mean'][index]))
        delta = np.array(diff)
        print(delta)
        chi2_val = self.calculate_chi2(delta, cov_matrix)
        print("chi2:", chi2_val)
        return chi2_val

    def do_fit(self, ):

        graph = self.set_tgrapherrors(self.mass_mean, self.pt_mean)
        graph.Fit("funct_nominal")

        self.slope = self.funct_nominal.GetParameter(0)
        self.slope_err = self.funct_nominal.GetParError(0)
        self.intercept = self.funct_nominal.GetParameter(1)
        self.intercept_err = self.funct_nominal.GetParError(1)
        self.chi2 = self.custom_chi2()
        self.ndof = self.funct_nominal.GetNDF()

        return self.slope, self.slope_err, self.intercept, self.intercept_err

    def do_sys_fit(self):
        self.slope_sys = 0
        self.intercept_sys = 0
        for sys_name in self.sys_mean_pt.keys():
            slope_diffs = []
            intercept_diffs = []
            if sys_name == "FSR":
                graph = self.set_tgrapherrors(self.sys_mean_mass[sys_name][0], self.sys_mean_pt[sys_name][0])
                graph.Fit("funct_sys")
                slope_nominal, intercept_nominal = self.funct_sys.GetParameter(0), self.funct_sys.GetParameter(1)
                del graph

                graph = self.set_tgrapherrors(self.sys_mean_mass[sys_name][1], self.sys_mean_pt[sys_name][1])
                graph.Fit("funct_sys")
                slope_pythia, intercept_pythia = self.funct_sys.GetParameter(0), self.funct_sys.GetParameter(1)
                del graph

                slope_diffs.append(slope_pythia - slope_nominal)
                intercept_diffs.append(intercept_pythia - intercept_nominal)
            else:
                for sys_index in range(len(self.sys_mean_pt[sys_name])):
                    graph = self.set_tgrapherrors(self.sys_mean_mass[sys_name][sys_index], self.sys_mean_pt[sys_name][sys_index])
                    graph.Fit("funct_sys")
                    slope_sys, intercept_sys = self.funct_sys.GetParameter(0), self.funct_sys.GetParameter(1)
                    del graph

                    # this is envelope...
                    slope_diffs.append(slope_sys - self.slope)
                    intercept_diffs.append(intercept_sys - self.intercept)

            slope_diffs = np.array(slope_diffs)
            intercept_diffs = np.array(intercept_diffs)
            if sys_name == "pdf" or "_stat" in sys_name:
                self.slope_sys += np.mean(slope_diffs ** 2)
                self.intercept_sys += np.mean(intercept_diffs ** 2)
            else:
                if sys_name == "scale":
                    self.slope_sys += np.sum((slope_diffs) ** 2)
                    self.intercept_sys += np.sum((intercept_diffs) ** 2)
                else:
                    self.slope_sys += np.sum((slope_diffs) ** 2)
                    self.intercept_sys += np.sum((intercept_diffs) ** 2)

        self.slope_sys = np.sqrt(self.slope_sys)
        self.intercept_sys = np.sqrt(self.intercept_sys)
        return self.slope_sys, self.intercept_sys

    def do_multi_error_fit(self):
        graph = TGraphMultiErrors("isr_fit", "isr_fit",
                                  len(self.mass_mean['mean']),
                                  self.mass_mean['mean'].values, self.pt_mean['mean'].values,
                                  self.mass_mean['stat'].values, self.mass_mean['stat'].values,
                                  self.pt_mean['stat'].values, self.pt_mean['stat'].values)

        # graph.SetSumErrorsMode(1)
        for index in range(len(self.pt_mean)):
            graph.SetPointError(index,
                                len(self.mass_mean.iloc[0][2:-2]),
                                self.mass_mean.iloc[index]['sys'],
                                self.mass_mean.iloc[index]['sys'],
                                array('d', self.pt_mean.iloc[index][2:-2].values),
                                array('d', self.pt_mean.iloc[index][2:-2].values),)


        graph.Fit("funct_nominal")

        self.slope = self.funct_nominal.GetParameter(0)
        self.slope_err = self.funct_nominal.GetParError(0)
        self.intercept = self.funct_nominal.GetParameter(1)
        self.intercept_err = self.funct_nominal.GetParError(1)
        self.chi2 = self.funct_nominal.GetChisquare()
        self.ndof = self.funct_nominal.GetNDF()

        return self.slope, self.slope_err, self.intercept, self.intercept_err
