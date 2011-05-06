import ROOT
from base_model import BaseModel
import pyWIMP.WIMPPdfs as pdfs
import math

class GammaLineModel(BaseModel):
    def __init__(self, 
                 basevars,
                 mean,
                 sigma,
                 lifetime = None,
                 live_time = None
                 ):
        BaseModel.__init__(self, basevars)



        name = str(self.get_tag()) + "_" + str(mean.getVal())
        if lifetime:
            name += "_lt_"
            name += str(lifetime.getVal()) 
       
	    # Gamma pdf
        self.energy_pdf = ROOT.RooGaussian("gamma_line_%s" % name, 
                                           "GammaPdf_%s" % name, 
                                           basevars.get_energy(),
                                           mean, sigma)
        self.lifetime = lifetime
        if not lifetime:
            self.gamma_pdf = self.energy_pdf

        if lifetime:
            afactor = math.log(2)*365.25
            self.local_lifetime = ROOT.RooFormulaVar(
                                    "local_lifetime_%s" % name, 
                                    "local_lifetime_%s" % name, 
                                    "-%f/@0" % afactor, 
                                    ROOT.RooArgList(lifetime))
            self.time_pdf = pdfs.MGMExponential("time_pdf_%s" % 
                              str(self.local_lifetime.getVal()), 
                              "TimePdf", 
                              basevars.get_time(),
                              self.local_lifetime)
            self.gamma_pdf = ROOT.RooProdPdf("GammaLine%s" % name, 
                                            "Gamma Line %s" % name, 
                                            self.energy_pdf, 
                                            self.time_pdf, 1e-8)
            self.time_pdf.SetRegionsOfValidity( live_time )

    def get_model(self):
        return self.gamma_pdf

    def get_lifetime(self):
        return self.lifetime
        #return self.energy_pdf


class GammaLineFactory:
    created = {}

    @classmethod
    def generate(cls,mean_value, lifetime_value, basevars):
      return cls.generate(mean_value, mean_value*0.05, 0.1*mean_value, 0.05*mean_value,
                   lifetime_value, basevars)
    @classmethod
    def generate(cls,mean_value, mean_error, sigma, sigma_error, lifetime_value, basevars, name=None, live_time = None):
        if mean_value in cls.created.keys(): return cls.created[mean_value][0]
        if not name:
            name = str(mean_value) + "_" + str(lifetime_value)
        mean = ROOT.RooRealVar("mean_%s" % name,"mean_%s" % name, 
                               mean_value, 
                               mean_value-mean_error, 
                               mean_value+mean_error)
        if mean_error == 0:
            mean.setConstant()
        sigma = ROOT.RooRealVar("sigma_%s" % name,"sigma_%s" % name, 
                                sigma, 
                                sigma - sigma_error, 
                                sigma + sigma_error)
        if sigma_error ==0:
            sigma.setConstant()
        try:
            lifetime_value.IsA()
            lifetime = lifetime_value
        except AttributeError:
            if lifetime_value == 0:
                lifetime = None
            else:
                lifetime = ROOT.RooRealVar("lifetime_%s" % name,
                                           "lifetime_%s" % name, 
                                           lifetime_value, 
                                           0.1*lifetime_value, 
                                           10*lifetime_value)
        gamma_line = GammaLineModel(basevars, mean, sigma, lifetime, live_time)
        cls.created[mean_value] = (gamma_line, mean, sigma, lifetime)
        return cls.created[mean_value][0]
        
