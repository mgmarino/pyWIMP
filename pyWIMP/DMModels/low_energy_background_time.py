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
                     ("73As", 11.103, 0.04, 0.04, 0.104, 1.414, 0.00, 80.30,        20,  200), 
                     ("68Ge", 10.367, 0.10, 0.04, 0.108, 1.302, 0.00, 270.95,       80, 5000), 
                     #("68GeFast", 10.367, 0.10, 0.04, 0.108, 1.302, 0.00, 11.26,    0,  10000), 
                     ("68Ga", 9.659,  0.10, 0.04, 0.112, 1.196, 0.00, 270.95,       100, 600), 
                     ("65Zn", 8.989,  0.10, 0.04, 0.117, 1.100, 0.00, 243.66,       100, 3000), 
                     #("56Ni", 7.710,  0.2,  0.0,  0.),                        
                     ("58Co", 7.010,  0.1, 0.0,  0.127, 0.929,  0.00, 70.86,        20, 200),
                     ("55Fe", 6.540,  0.2,  0.0,  0.133, 0.842, 0.00, 100*365.25*2.737, 0, 200), 
                     ("54Mn", 5.990,  0.15, 0.0,  0.139, 0.754, 0.00, 100*312.12,   0, 200),
                     ("51Cr", 5.460,  0.10, 0.0,  0.146, 0.682, 0.00, 27.7025,      40, 50),
                     ("49V",  4.970,  0.2, 0.0,  0.154, 0.604,  0.00, 100*300,          0, 200),
                    ]
                    #("Unknown", 0.9, 0.05, 0.1, 0.04, 0), 
        
        self.save_list = []
        self.high_pdf_list = ROOT.RooArgList()
        self.high_coefficienct_list = ROOT.RooArgList()
        self.low_pdf_list = ROOT.RooArgList()
        self.low_coefficienct_list = ROOT.RooArgList()

        ga68_time_pdf = None
        for name, mean, mean_error, sigma_error, l_k_ratio, l_energy, l_energy_err, decay_time, decay_err, max_ampl in mean_list:
            if name == "68Ga": 
                decay_time = ga68_time_pdf
                decay_err = None

            agamma = gamma_line_model.GammaLineFactory.generate(
                        basevars,name, 
                        mean_value=mean, 
                        mean_error=mean_error, 
                        sigma=get_sigma(mean*1e3)*1e-3, 
                        sigma_error=sigma_error, 
                        lifetime_value=decay_time, 
                        lifetime_err=decay_err, 
                        live_time=livetime)

            next_life = agamma.get_lifetime()
            if name == "68Ge": ga68_time_pdf = next_life 
            if name == "68Ga": next_life = ga68_time_pdf

            agamma_low = gamma_line_model.GammaLineFactory.generate(
                        basevars,name+"_l_line", 
                        mean_value=l_energy, 
                        mean_error=l_energy_err, 
                        sigma=get_sigma(l_energy*1e3)*1e-3, 
                        sigma_error=0.0, 
                        lifetime_value=next_life, 
                        lifetime_err=None, 
                        live_time=livetime)

            gamma = agamma.get_model()
            low_gamma = agamma_low.get_model()
            the_name = "%s_ampl" % gamma.GetName()
            new_var = ROOT.RooRealVar(the_name, the_name, 
                                      5, max_ampl, "Counts")
            low_ampl = ROOT.RooFormulaVar(the_name + "_low", the_name + "_low", "%f*@0" % l_k_ratio, ROOT.RooArgList(new_var))

            self.save_list.append((agamma, agamma_low, new_var, low_ampl))

            self.high_pdf_list.add(gamma)
            self.high_coefficienct_list.add(new_var)
            self.low_pdf_list.add(low_gamma)
            self.low_coefficienct_list.add(low_ampl)
          
        tag = "high"
        flat_coef = ROOT.RooRealVar("flat_coef_%s" % tag,
                                    "flat_coef_%s" % tag,
                                    2000, 10, 4000)
        # Flat pdf
        time_pdf_flat = pdfs.MGMPiecewiseFunction("time_pdf_flat_all", "time_pdf_flat_all", self.basevars.get_time())
        time_pdf_flat.SetRegionsOfValidity(livetime) 

        time_exp = ROOT.RooRealVar("time_exp", "time_exp", -0.01, -1, -1e-16)
        time_pdf_exp = pdfs.MGMExponential("time_pdf_exp", "time_pdf_exp", self.basevars.get_time(), time_exp)
        time_pdf_exp.SetRegionsOfValidity(livetime) 

        energy_pdf_flat = pdfs.MGMPiecewiseFunction("energy_pdf_flat_%s" % tag,  "energy_pdf_flat_%s" % tag, self.basevars.get_energy())
        energy_pdf_flat_regions = pdfs.MGMPiecewiseRegions()
        energy_pdf_flat_regions.InsertNewRegion(2.5, 15)
        energy_pdf_flat.SetRegionsOfValidity(energy_pdf_flat_regions)
        energy_flat_whole = ROOT.RooProdPdf("high_flat_whole", "high_flat_whole", energy_pdf_flat, time_pdf_flat)

        erf_amp = ROOT.RooRealVar("erf_amp", "erf_amp", 1000, 10, 6000)
        erf_pdf = pdfs.MGMPiecewiseFunction("erf_func",  "erf_func", self.basevars.get_energy())
        erf_pdf_flat_regions = pdfs.MGMPiecewiseRegions()
        erf_pdf_flat_regions.InsertNewRegion(2.5, 10.7) 
        erf_pdf.SetRegionsOfValidity(erf_pdf_flat_regions)
        #erf_pdf_whole = ROOT.RooProdPdf("high_erf_whole", "high_erf_whole", erf_pdf, time_pdf_flat)
        erf_pdf_whole = ROOT.RooProdPdf("high_erf_whole", "high_erf_whole", erf_pdf, time_pdf_exp)

        self.save_list.extend([flat_coef, energy_pdf_flat, erf_amp, erf_pdf, energy_pdf_flat_regions, erf_pdf_flat_regions])
        self.save_list.extend([time_pdf_flat, energy_flat_whole, erf_pdf_whole])
        self.save_list.extend([time_pdf_exp, time_exp])

        self.high_pdf_list.add(energy_flat_whole)
        self.high_coefficienct_list.add(flat_coef)

        self.high_pdf_list.add(erf_pdf_whole)
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
                                           -3., -5, -0.8)
        exp_constant_one.setError(0.5)
        exp_coef = ROOT.RooRealVar("exp_coef_%s" % tag,
                                   "exp_coef_%s" % tag,
                                   600, 10, 1000)
        low_flat_coef = ROOT.RooRealVar("flat_coef_%s" % tag,
                                        "flat_coef_%s" % tag,
                                        1000, 10, 2000)
        # Flat pdf
        low_erf_sig = ROOT.RooRealVar("low_sig", "low_sig", 0)
        low_erf_mean = ROOT.RooRealVar("low_erf_mean", "low_erf_mean", self.max_energy)
        low_erf_offset = ROOT.RooRealVar("low_offset", "low_offset", 0)
        low_erf_pdf = pdfs.MGMErfcFunction("erf_func_%s" % tag, "erf_func_%s" % tag, 
                                            basevars.get_energy(), 
                                            low_erf_mean, 
                                            low_erf_sig, 
                                            low_erf_offset)

        self.phase    = ROOT.RooRealVar("phase", "phase", 0, -ROOT.TMath.Pi(), ROOT.TMath.Pi()) 
        self.osc_ampl = ROOT.RooRealVar("osc_amplitude", "osc_amplitude", 0, 1.) 
        self.frequency = ROOT.RooRealVar("osc_frequency", "osc_frequency", 1., 0.1, 10.) 
        self.frequency.setConstant(True)
        osc_time_pdf = ROOT.MGMExponentialPlusSinusoid("osc_pdf", "osc_pdf", 
                                                       self.basevars.get_time(), time_exp, 
                                                       self.osc_ampl, 
                                                       self.frequency, 
                                                       self.phase) 
        osc_time_pdf.SetRegionsOfValidity(livetime)
        self.save_list.append(osc_time_pdf)


        energy_exp_pdf = ROOT.RooExponential("energy_pdf_exp", 
                                             "energy_pdf_exp", 
                                             basevars.get_energy(),
                                             exp_constant_one)

        exp_pdf_whole = ROOT.RooProdPdf("exp_pdf_whole", "exp_pdf_whole", energy_exp_pdf, osc_time_pdf)
        low_erf_whole = ROOT.RooProdPdf("low_erf_whole", "low_erf_whole", low_erf_pdf, time_pdf_exp)

        self.low_pdf_list.add(low_erf_whole)
        self.low_coefficienct_list.add(low_flat_coef)

        self.low_pdf_list.add(exp_pdf_whole)
        self.low_coefficienct_list.add(exp_coef)

        self.save_list.extend([energy_exp_pdf, low_erf_pdf, low_erf_offset, low_erf_mean, low_erf_sig,
                               low_flat_coef, exp_coef, exp_constant_one])
        self.save_list.extend([low_erf_whole, exp_pdf_whole])

        self.low_extend_list = ROOT.RooArgList()
        for i in range(self.low_pdf_list.getSize()): 
            pdf, coeff = self.low_pdf_list.at(i), self.low_coefficienct_list.at(i)
            extend = ROOT.RooExtendPdf("extend" + pdf.GetName(), "extend" + pdf.GetName(), pdf, coeff)
            self.save_list.append(extend)
            self.low_extend_list.add(extend)



    def get_extend_components(self):
        return self.low_extend_list, self.high_extend_list

    def get_osc_ampl(self):
        return self.osc_ampl

    def set_frequency(self, val):
        self.frequency.setVal(val)

    def set_osc_ampl_const(self, val):
        self.set_osc_ampl(val)
        self.osc_ampl.setConstant(True)

    def set_osc_ampl(self, val):
        self.osc_ampl.setVal(val)
        if val == 0.0:
            self.osc_ampl.setConstant(True)
            self.phase.setConstant(True)
        else:
            self.osc_ampl.setConstant(False)
            self.phase.setConstant(False)


