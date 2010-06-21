import ROOT
from base_model import BaseModel
from gamma_line_model import GammaLineFactory
import math

def get_sigma(energy_in_eV):
    noise = 70.4843
    energy_par = 6.49228e-2
    return math.sqrt(noise*noise + 2.96*2.96*energy_par*energy_in_eV)

class GaussianSignalModel(BaseModel):
    def __init__(self, 
                 basevars, 
                 mean_of_signal=20):
        # Normally, we don't want to do this, but this keeps 
        # it from importing this module until the last moment.
        BaseModel.__init__(self, basevars)
        sig = get_sigma(mean_of_signal*1e3)*1e-3
        self.class_model = GammaLineFactory.generate(mean_of_signal, 0, 
                                                     sig, 0, 
                                                     0, basevars)
        self.get_model().SetName("Gauss_Signal_%g" % mean_of_signal)
        self.get_model().SetTitle("Gauss_Signal_%g" % mean_of_signal)

        self.normalization = 1./(math.sqrt(ROOT.TMath.TwoPi())*sig)
        
    
    def get_model(self):     
        return self.class_model.get_model()

    def get_normalization(self):     
        # Prints normalization because roofit doesn't return correct norm
        return self.normalization
