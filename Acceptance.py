from pygments.lexer import default
import numpy as np
import ROOT
from Hist import Hist


class Acceptance:
    def __init__(self, mc_hist_full_phase, mc_hist_in_acceptance):
        self.mc_hist_full_phase = mc_hist_full_phase
        self.mc_hist_in_acceptance = mc_hist_in_acceptance
        self._set_acceptance()
        # Systematic on acceptance correction?
        # Statistical,

    def _set_acceptance(self):
        self.acceptance_hist = self.mc_hist_full_phase.divide(self.mc_hist_in_acceptance)
        # TODO check acceptance_hist !

    def do_correction(self, input_hist):
        # input_hist from unfolder
        out_hist = input_hist.multiply(self.acceptance_hist)  #
        return out_hist

    def create_stat_accept(self, n_toy=100):
        default_accept = self.acceptance_hist.get_raw_hist()
        stat_accept_template = self.acceptance_hist.get_raw_hist().Clone('stat_accept')
        stat_accept_template.Reset()

        accepts = []
        nx = default_accept.GetNbinsX()
        for i in range(n_toy):
            temp_accept = stat_accept_template.Clone(f"stat{i}")
            for i_x in range(nx+2):
                mean = default_accept.GetBinContent(i_x)
                error = default_accept.GetBinError(i_x)
                temp_accept.SetBinContent(i_x, ROOT.gRandom.Gaus(mean, error))
                temp_accept.SetBinError(i_x, error)
            accepts.append(Hist(temp_accept, f"stat_accept{i}"))
        return accepts

    def get_accept_stat(self, input_hist):
        # repeat acceptance correction over acceptance stat
        default_hist = input_hist.multiply(self.acceptance_hist)
        default_vals = np.array(default_hist.to_numpy()[0])
        toy_diffs = []
        for idx, stat_accept in enumerate(self.create_stat_accept()):
            temp_hist = input_hist.multiply(stat_accept)
            out_vals = np.array(temp_hist.to_numpy()[0])
            toy_diffs.append(np.abs(out_vals - default_vals))
        # per‐bin pop-std
        deltas = np.std(toy_diffs, axis=0, ddof=0).tolist()

        up = default_hist.get_raw_hist().Clone("accept_stat_up")
        down = default_hist.get_raw_hist().Clone("accept_stat_down")
        up.Reset()
        down.Reset()
        nbin = default_hist.get_raw_hist().GetNbinsX()
        for i in range(1, nbin + 1):
            val = default_hist.get_raw_hist().GetBinContent(i)
            d = deltas[i - 1]  # underflow/overflow not considered
            up.SetBinContent(i, val + d)
            down.SetBinContent(i, val - d)
            up.SetBinError(i, default_hist.get_raw_hist().GetBinError(i))
            down.SetBinError(i, default_hist.get_raw_hist().GetBinError(i))
        return {"up": up, "down": down}
