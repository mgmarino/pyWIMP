import ROOT
from base_model import BaseModel
import gamma_line_model 
import math
import pyWIMP.WIMPPdfs as pdfs

def get_sigma(energy_in_eV):
    noise = 70.4843
    energy_par = 6.49228e-2
    return math.sqrt(noise*noise + 2.96*2.96*energy_par*energy_in_eV)

class FittingModelTime(BaseModel):
    def __init__(self, 
                 basevars,
                 max_energy, 
                 livetime):
        self.max_energy = max_energy
        BaseModel.__init__(self, basevars)
        self.initialize(basevars, livetime)

    def initialize(self, basevars, livetime):
                          #Name  :  Energy  :  EnergyErr  :  SigmaErr :  L/K ratio : L-energy : L-energy error : Livetime : max_ampl 
        mean_list = [
                     ("73As", 11.103, 0.04, 0.04, 0.104, 1.414, 0.00, 80.30,   500), 
                     ("68Ge", 10.367, 0.10, 0.04, 0.108, 1.302, 0.00, 270.95,  10000), 
                     ("68Ga", 9.659,  0.10, 0.04, 0.112, 1.196, 0.00, 100.0,   3000), 
                     ("65Zn", 8.989,  0.10, 0.04, 0.117, 1.100, 0.00, 243.66,  10000), 
                     #("56Ni", 7.710,  0.2,  0.0,  0.),                        
                     ("58Co", 7.110,  0.04, 0.0,  0.127, 0.929, 0.00, 70.86,   300),
                     #("55Fe", 6.540,  0.2,  0.0,  0.133, 0.842, 0.00, 365.25*2.737, 3000), 
                     ("54Mn", 5.990,  0.10, 0.0,  0.139, 0.754, 0.00, 312.12,  500),
                     ("51Cr", 5.460,  0.10, 0.0,  0.146, 0.682, 0.00, 27.7025, 300),
                     ("49V",  4.970,  0.2, 0.0,  0.154, 0.604, 0.00,  300,     500),
                    ]
                    #("Unknown", 0.9, 0.05, 0.1, 0.04, 0), 
        
        self.save_list = []
        self.high_pdf_list = ROOT.RooArgList()
        self.high_coefficienct_list = ROOT.RooArgList()
        self.low_pdf_list = ROOT.RooArgList()
        self.low_coefficienct_list = ROOT.RooArgList()

        for name, mean, mean_error, sigma_error, l_k_ratio, l_energy, l_energy_err, decay_time, max_ampl in mean_list:
            agamma = gamma_line_model.GammaLineFactory.generate(mean, 
                        mean_error, get_sigma(mean*1e3)*1e-3, 
                        sigma_error, decay_time, basevars,name, livetime)
            agamma_low = gamma_line_model.GammaLineFactory.generate(l_energy, 
                        l_energy_err, get_sigma(l_energy*1e3)*1e-3, 
                        0.0, agamma.get_lifetime(), basevars,name+"_l_line", livetime)

            gamma = agamma.get_model()
            low_gamma = agamma_low.get_model()
            the_name = "%s_ampl" % gamma.GetName()
            new_var = ROOT.RooRealVar(the_name, the_name, 
                                      1e-15, max_ampl, "Counts")
            low_ampl = ROOT.RooFormulaVar(the_name + "_low", the_name + "_low", "%f*@0" % l_k_ratio, ROOT.RooArgList(new_var))

            self.save_list.append((agamma, agamma_low, new_var, low_ampl))

            self.high_pdf_list.add(gamma)
            self.high_coefficienct_list.add(new_var)
            self.low_pdf_list.add(low_gamma)
            self.low_coefficienct_list.add(low_ampl)
          
        tag = "high"
        flat_coef = ROOT.RooRealVar("flat_coef_%s" % tag,
                                    "flat_coef_%s" % tag,
                                    10, 4000)
        # Flat pdf
        energy_pdf_flat = ROOT.RooPolynomial("energy_pdf_flat_%s" % tag, 
                                             "energy_pdf_flat_%s" % tag, 
                                             basevars.get_energy())
        energy_pdf_flat = pdfs.MGMPiecewiseFunction("energy_pdf_flat_%s" % tag,  "energy_pdf_flat_%s" % tag, self.basevars.get_energy())
        energy_pdf_flat_regions = pdfs.MGMPiecewiseRegions()
        energy_pdf_flat_regions.InsertNewRegion(3.1760041, 14.6620021)
        energy_pdf_flat.SetRegionsOfValidity(energy_pdf_flat_regions)
        asig = ROOT.RooRealVar("asig", "asig", 0.1)
        amean = ROOT.RooRealVar("amean", "amean", 10.367)
        anoffset = ROOT.RooRealVar("anoffset", "anoffset", 0)
        erf_amp = ROOT.RooRealVar("erf_amp", "erf_amp", 1e-15, 6000)
        erf_pdf = pdfs.MGMErfcFunction("erf_func", "erf_func", basevars.get_energy(), amean, asig, anoffset)
        erf_pdf = pdfs.MGMPiecewiseFunction("erf_func",  "erf_func", self.basevars.get_energy())
        erf_pdf_flat_regions = pdfs.MGMPiecewiseRegions()
        erf_pdf_flat_regions.InsertNewRegion(3.1760041, 10.367)
        erf_pdf.SetRegionsOfValidity( erf_pdf_flat_regions)

        self.save_list.extend([flat_coef, energy_pdf_flat, asig, amean, anoffset, erf_amp, erf_pdf, energy_pdf_flat_regions, erf_pdf_flat_regions])

        self.high_pdf_list.add(energy_pdf_flat)
        self.high_coefficienct_list.add(flat_coef)
        self.high_pdf_list.add(erf_pdf)
        self.high_coefficienct_list.add(erf_amp)

        self.high_extend_list = ROOT.RooArgList()
        for i in range(self.high_pdf_list.getSize()): 
            pdf, coeff = self.high_pdf_list.at(i), self.high_coefficienct_list.at(i)
            extend = ROOT.RooExtendPdf("extend" + pdf.GetName(), "extend" + pdf.GetName(), pdf, coeff)
            self.save_list.append(extend)
            self.high_extend_list.add(extend)



        ########################################################################
        # Low pdf
        tag = "low"
        exp_constant_one = ROOT.RooRealVar("expo_const_one%s" % tag,
                                           "expo_const_one%s" % tag,
                                           -3., -10, -0.8)
        exp_constant_one.setError(0.5)
        exp_coef = ROOT.RooRealVar("exp_coef_%s" % tag,
                                   "exp_coef_%s" % tag,
                                   10, 3000)
        low_flat_coef = ROOT.RooRealVar("flat_coef_%s" % tag,
                                        "flat_coef_%s" % tag,
                                        10, 3000)
        # Flat pdf
        low_erf_sig = ROOT.RooRealVar("low_sig", "low_sig", 0)
        low_erf_mean = ROOT.RooRealVar("low_erf_mean", "low_erf_mean", self.max_energy)
        low_erf_offset = ROOT.RooRealVar("low_offset", "low_offset", 0)
        low_erf_pdf = pdfs.MGMErfcFunction("erf_func_%s" % tag, "erf_func_%s" % tag, 
                                            basevars.get_energy(), 
                                            low_erf_mean, 
                                            low_erf_sig, 
                                            low_erf_offset)

        energy_exp_pdf = ROOT.RooExponential("energy_pdf_exp", 
                                             "energy_pdf_exp", 
                                             basevars.get_energy(),
                                             exp_constant_one)

        self.low_pdf_list.add(low_erf_pdf)
        self.low_coefficienct_list.add(low_flat_coef)

        self.low_pdf_list.add(energy_exp_pdf)
        self.low_coefficienct_list.add(exp_coef)

        self.save_list.extend([energy_exp_pdf, low_erf_pdf, low_erf_offset, low_erf_mean, low_erf_sig,
                               low_flat_coef, exp_coef, exp_constant_one])

        self.time_efficiency = pdfs.MGMPiecewiseFunction("time_eff", "time_eff", self.basevars.get_time())
        #self.energy_range = pdfs.MGMPiecewiseFunction("energy_range", "energy_range", self.basevars.get_energy())
        #self.high_energy_range = pdfs.MGMPiecewiseFunction("high_energy_range", "high_energy_range", self.basevars.get_energy())
        
        self.low_extend_list = ROOT.RooArgList()
        for i in range(self.low_pdf_list.getSize()): 
            pdf, coeff = self.low_pdf_list.at(i), self.low_coefficienct_list.at(i)
            extend = ROOT.RooExtendPdf("extend" + pdf.GetName(), "extend" + pdf.GetName(), pdf, coeff)
            self.save_list.append(extend)
            self.low_extend_list.add(extend)



    def get_high_components(self):
        return self.high_pdf_list, self.high_coefficienct_list

    def get_low_components(self):
        return self.low_pdf_list, self.low_coefficienct_list

    def get_extend_components(self):
        return self.low_extend_list, self.high_extend_list

    def get_time_efficiency(self):
        return self.time_efficiency

    def get_energy_ranges(self):
        return self.energy_range, self.high_energy_range




