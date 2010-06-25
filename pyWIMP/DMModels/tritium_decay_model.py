import ROOT
from base_model import BaseModel
from beta_decay_model import BetaDecayModel
from flat_model import FlatModel
import pyWIMP.WIMPPdfs as pdfs  

class TritiumDecayModel(BaseModel):
    """
    This class returns a 2-D tritium decay spectrum with 
    a flat background.  Inputs are:

    - tritium_exposure: Tritium Exposure time (days)
    - tritium_activation: Tritium activation rate in atoms/kg/day 
    - mass_of_detector: Mass of the detector
    - flat_background_rate: Flat background rate, in counts/kg/keV/day

    It returns an extended model with the relative amounts are set correctly (but allowed to float)
    This means that if you are doing toy model fitting with the returned model before to generate
    the events (RooFit will automatically generate the correct number of events from an extended pdf)
    and SAVE all the variables. 

        # Save the values of the parameters to reset at the end
        var_cache = ROOT.ostringstream() 
        model.getVariables().writeToStream(var_cache, False)

        data = model.generate()
        ...
        ... (Fitting)
        # Reset the variables 
        model.getVariables().readFromStream(ROOT.istringstream(var_cache.str()), False)
    """
    def __init__(self, 
                 basevars,
                 tritium_exposure,
                 tritium_activation,
                 mass_of_detector,
                 flat_background_rate):
        BaseModel.__init__(self, basevars)


        # Set up the tritium decay
        self.beta_decay = beta_decay_model.BetaDecayModel(basevars, 18.6, 12.36)
        self.beta_model_amp = ROOT.RooRealVar("tritium_amplitude", 
                                              "Tritium Amplitude", 
                                              events, 1e-15, 100000)
        self.beta_model = beta_decay.get_model() 
        self.beta_model_extend = ROOT.RooExtendPdf("tritium_extend_model", 
                                                   "Tritium Extended Model", 
                                                   self.beta_model, self.beta_model_amp)

        self.flat_background = FlatModel(basevars)
        self.flat_amp = ROOT.RooRealVar("flat_amplitude", "Flat Background amplitude", 
                                        events_in_flat, 1e-15, 100000)
        self.flat_model = self.flat_background.get_model()
        self.flat_model_extend = ROOT.RooExtendPdf("flat_extend_model", 
                                                   "Flat Extended Model", 
                                                   self.flat_model, self.flat_amp)
        
        self.total_background = ROOT.RooAddPdf("total_tritium_background", 
                                               "Total Background (Tritium Model)", 
                                               ROOT.RooArgList(self.flat_model_extend, 
                                               self.beta_model_extend))



    def get_model(self):
        return self.total_background


