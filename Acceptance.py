

class Acceptance:
    def __init__(self, mc_hist_full_phase, mc_hist_in_acceptance):
        self.mc_hist_full_phase = mc_hist_full_phase
        self.mc_hist_in_acceptance = mc_hist_in_acceptance

        self._set_acceptance()

    def _set_acceptance(self):
        self.acceptance_hist = self.mc_hist_full_phase.divide(self.mc_hist_in_acceptance)

    def do_correction(self, input_hist):
        out_hist = input_hist.multiply(self.acceptance_hist)
        return out_hist

