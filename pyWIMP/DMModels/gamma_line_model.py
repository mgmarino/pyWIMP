import ROOT
from base_model import BaseModel
import pyWIMP.WIMPPdfs as pdfs
import math

class GammaLineModel(BaseModel):
    def __init__(self, 
                 basevars,
                 gamma_pdf,
                 name,
                 lifetime,
                 live_time
                 ):
        BaseModel.__init__(self, basevars)

	# Gamma pdf
        self.energy_pdf = gamma_pdf
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

    @classmethod
    def with_mean_sigma(cls, basevars, mean, sigma, lifetime = None, live_time = None):
        name = str(mean.getVal())
        if lifetime:
           name += "_lt_"
           name += str(lifetime.getVal()) 
        energy_pdf = ROOT.RooGaussian("gamma_line_%s" % name, 
                                      "GammaPdf_%s" % name, 
                                      basevars.get_energy(),
                                      mean, sigma)
        return cls(basevars, energy_pdf, name, lifetime, live_time)

    @classmethod
    def with_gamma_pdf(cls, basevars, energy_pdf, name, lifetime, live_time):
        return cls(basevars, energy_pdf, name, lifetime, live_time)

    def get_model(self):
        return self.gamma_pdf

    def get_lifetime(self):
        return self.lifetime

    def get_energy_pdf(self):
        return self.energy_pdf


class GammaLineFactory:
    created = []
    save_list = []

    @classmethod
    def generate(cls, name, mean_value, lifetime_value, basevars):
        return cls.generate(basevars, name, mean_value=mean_value, mean_error=mean_value*0.05, sigma=0.1*mean_value, sigma_error=0.05*mean_value,
                   lifetime_value=lifetime_value)
    @classmethod
    def generate(cls, basevars, name, **kw):
        """
        The inputs to this are as follows:
        mean_vaue: expected mean value
        mean_error: expected mean error.  If 0, will be set to constant
        sigma:
        sigma_error:
        lifetime_value: either 
        lifetime_err: either 
        gamma_pdf: the gamma_pdf
        live_time: should be a MGMPiecewiseRegions object
        """
        lifetime_value = kw['lifetime_value']
        lifetime_err   = kw['lifetime_err']
        live_time      = kw['live_time']
        try:
            lifetime_value.IsA()
            lifetime = lifetime_value
        except AttributeError:
            if lifetime_value == 0:
                lifetime = None
            else:
                min_value = lifetime_value - lifetime_err
                if min_value <= 0: min_value = 1e-16
                lifetime = ROOT.RooRealVar("lifetime_%s" % name,
                                           "lifetime_%s" % name, 
                                           lifetime_value, 
                                           min_value, 
                                           lifetime_value + lifetime_err)
                if lifetime_err == 0: lifetime.setConstant()
                cls.save_list.append(lifetime)

        if "gamma_pdf" not in kw.keys():
            mean_value = kw['mean_value']
            mean_error = kw['mean_error']
            mean = ROOT.RooRealVar("mean_%s" % name,"mean_%s" % name, 
                                   mean_value, 
                                   mean_value-mean_error, 
                                   mean_value+mean_error)
            if mean_error == 0: mean.setConstant()
            sigma_val = kw['sigma']
            sigma_error = kw['sigma_error']
            sigma = ROOT.RooRealVar("sigma_%s" % name,"sigma_%s" % name, 
                                    sigma_val, 
                                    sigma_val - sigma_error, 
                                    sigma_val + sigma_error)
            cls.save_list.append((mean, sigma))

            if sigma_error ==0:
                sigma.setConstant()

            gamma_line = GammaLineModel.with_mean_sigma(basevars, mean, sigma, lifetime, live_time)
        else:
            gamma_line = GammaLineModel.with_gamma_pdf(basevars, kw['gamma_pdf'], name, lifetime, live_time)

        cls.created.append(gamma_line)
        return cls.created[-1]
        
