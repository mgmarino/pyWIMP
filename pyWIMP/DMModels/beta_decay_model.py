import ROOT
from base_model import BaseModel
import pyWIMP.WIMPPdfs as pdfs  

class BetaDecayModel(BaseModel):
    """
    This class returns a 2-D beta decay spectrum with 
    a certain life_time (measured in the unit of time)
    """
    mass_of_electron = ROOT.RooRealVar("mass_of_electron", "Mass of Electron", 511, "keV")
    def __init__(self, 
                 basevars,
                 q_value,
                 lifetime = None):
        BaseModel.__init__(self, basevars)

        # Flat pdf
        tag = self.get_tag()
        name = str(self.get_tag()) + "_" + str(q_value)
        if lifetime:
            name += "_lt_"
            name += str(lifetime) 
 
        if not lifetime:
            self.time_pdf = ROOT.RooPolynomial("time_beta_" + name,
                                               "Time Beta " + name,
                                               self.basevars.get_time())
        else:
            self.lifetime = ROOT.RooRealVar("lifetime" + name,
                                            "lifetime" + name,
                                            lifetime, self.basevars.get_time().getUnit())
            self.local_lifetime = ROOT.RooFormulaVar(
                                    "local_lifetime_%s" % name, 
                                    "local_lifetime_%s" % name, 
                                    "-0.693147181/@0", 
                                    ROOT.RooArgList(self.lifetime))
            self.time_pdf = ROOT.RooExponential("time_beta_" + name, 
                                                "Time Beta " + name, 
                                                basevars.get_time(),
                                                self.local_lifetime)


        self.q_value = ROOT.RooRealVar("q_value" + name, 
                                       "q_value" + name, 
                                        q_value)

        self.energy_pdf = pdfs.MGMBetaDecayFunction("energy_beta_" + name, 
                                                    "Energy Beta " + name, 
                                                    self.basevars.get_energy(), 
                                                    self.mass_of_electron, 
                                                    self.q_value)
        self.model_pdf = ROOT.RooProdPdf("beta_time_and_energy_pdf_%s" % name, 
                                         "Beta Time Energy Pdf " + name, 
                                         self.time_pdf, 
                                         self.energy_pdf)

    def get_model(self):
        return self.model_pdf


