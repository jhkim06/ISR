from ROOT import TF1, TGraphErrors, TGraphMultiErrors
from array import array


class ISRLinearFitter:
    def __init__(self, mass_mean, pt_mean):
        self.mass_mean = mass_mean
        self.pt_mean = pt_mean

        self.fit_funct = "2.* [0] * TMath::Log(x) + [1]"
        self.slope = 0
        self.intercept = 0
        self.slope_err = 0
        self.intercept_err = 0

    def do_fit(self):

        graph = TGraphErrors(len(self.mass_mean['mean']))
        for index in range(len(self.pt_mean['mean'])):
            graph.SetPoint(index, self.mass_mean['mean'][index], self.pt_mean['mean'][index])
            graph.SetPointError(index, self.mass_mean['stat'][index], self.pt_mean['stat'][index])

        funct = TF1("funct", self.fit_funct, 50, 400)
        graph.Fit("funct")

        self.slope = funct.GetParameter(0)
        self.slope_err = funct.GetParError(0)
        self.intercept = funct.GetParameter(1)
        self.intercept_err = funct.GetParError(1)

        #print(f"a = {funct.GetParameter(0):.6f} +/- {funct.GetParError(0):.6f}")
        #print(f"b = {funct.GetParameter(1):.6f} +/- {funct.GetParError(1):.6f}")
        return self.slope, self.slope_err, self.intercept, self.intercept_err

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


        funct = TF1("funct", self.fit_funct, 50, 400)
        graph.Fit("funct")

        self.slope = funct.GetParameter(0)
        self.slope_err = funct.GetParError(0)
        self.intercept = funct.GetParameter(1)
        self.intercept_err = funct.GetParError(1)

        # print(f"a = {funct.GetParameter(0):.6f} +/- {funct.GetParError(0):.6f}")
        # print(f"b = {funct.GetParameter(1):.6f} +/- {funct.GetParError(1):.6f}")
        return self.slope, self.slope_err, self.intercept, self.intercept_err
