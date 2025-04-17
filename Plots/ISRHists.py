from Hist import Hist


class ISRHists(Hist):
    def __init__(self,
                 mass_bins,  # mass windows
                 pt_bins,
                 is_pt=True):  # pt cut

        self.is_pt = is_pt

        # detector level hists
        self.measurement_hist = None
        self.signal_hist = None
        self.signal_fake_hist = None
        self.background_hists = None
        # self.draw()

        '''
        def draw_detector_level(): 
        plotter = self.measurement_hist.draw()
        self.sigal_hist.draw(plotter)
        '''

        # unfolding
        self.response_matrix = None
        self.unfolded_measurement_hist = None

        # acceptance correction
        self.acceptance_corrected_measurement_hist = None
        # check is_2d?