import pyWIMP.DMModels.base_model as base_model
import pyWIMP.WIMPPdfs as pdfs
import ROOT
import sys

#ROOT.gROOT.SetBatch()
ROOT.RooRandom.randomGenerator().SetSeed(1)
basevars = base_model.BaseVariables(0, 1, 0, 3)
time = basevars.get_time()
energy = basevars.get_energy()
livetime = pdfs.MGMPiecewiseRegions()
livetime.InsertNewRegion(0, time.getMax())
time_exp = ROOT.RooRealVar("time_exp", "time_exp", -0.8, -3, 0.5)
tag = "low"
exp_constant_one = ROOT.RooRealVar("expo_const_one%s" % tag,
                                   "expo_const_one%s" % tag,
                                   -3., -5, -0.8)
exp_coef = ROOT.RooRealVar("exp_coef_%s" % tag,
                           "exp_coef_%s" % tag,
                           200, 1e-15, 2000)

# Flat pdf
phase    = ROOT.RooRealVar("phase", "phase", 0.1, -0.5, 0.5)#-ROOT.TMath.Pi(), ROOT.TMath.Pi()) 
period = ROOT.RooRealVar("osc_period", "osc_period", 1, 0.1, 5) 
exp_time_pdf = ROOT.MGMPiecewiseFunction("osc_exp", "osc_exp", 
                                    basevars.get_time()) 
osc_time_pdf = ROOT.MGMExponentialPlusSinusoid("osc_pdf", "osc_pdf", 
                                               basevars.get_time(),
                                               period, phase) 
exp_time_pdf.SetRegionsOfValidity(livetime)
osc_time_pdf.SetRegionsOfValidity(livetime)

osc_ampl = ROOT.RooRealVar("osc_amplitude", "osc_amplitude", 100, 1e-15, 1000) 
#time_pdf = ROOT.RooAddPdf("total_time", "total_time", osc_time_pdf, exp_time_pdf, osc_ampl)


energy_exp_pdf = ROOT.RooExponential("energy_pdf_exp", 
                                     "energy_pdf_exp", 
                                     basevars.get_energy(),
                                     exp_constant_one)

exp_pdf_whole = ROOT.RooProdPdf("exp_pdf_whole", "exp_pdf_whole", energy_exp_pdf, osc_time_pdf)
exp_pdf_whole_osc = ROOT.RooProdPdf("exp_pdf_whole_osc", "exp_pdf_whole_osc", energy_exp_pdf, osc_time_pdf)

low_erf_sig = ROOT.RooRealVar("low_sig", "low_sig", 0)
low_erf_mean = ROOT.RooRealVar("low_erf_mean", "low_erf_mean", energy.getMax())
low_erf_offset = ROOT.RooRealVar("low_offset", "low_offset", 0)
low_erf_pdf = pdfs.MGMErfcFunction("erf_func_%s" % tag, "erf_func_%s" % tag, 
                                    basevars.get_energy(), 
                                    low_erf_mean, 
                                    low_erf_sig, 
                                    low_erf_offset)

low_flat_coef = ROOT.RooRealVar("flat_coef_%s" % tag,
                                "flat_coef_%s" % tag,
                                1500, 1e-15, 5000)

extend = ROOT.RooExtendPdf("extend", "extend", exp_pdf_whole, exp_coef)
extendtwo = ROOT.RooExtendPdf("extendtwo", "extendtwo", low_erf_pdf, low_flat_coef)
extendthree = ROOT.RooExtendPdf("extendthree", "extendthree", exp_pdf_whole_osc, osc_ampl)
final = ROOT.RooAddPdf("final", "final", ROOT.RooArgList(extend, extendtwo, extendthree))

set = ROOT.RooArgSet(energy, time)
data = final.generate(set, -2)
data.Print('v')

nll = final.createNLL(data)
minuit = ROOT.RooMinuit(nll)
minuit.migrad()
#minuit.minos(ROOT.RooArgSet(phase, period))
minuit.minos()

c1 = ROOT.TCanvas()
for var in (energy, time):
    frame = var.frame(ROOT.RooFit.Bins(32))
    data.plotOn(frame)
    final.plotOn(frame)
    final.plotOn(frame, ROOT.RooFit.Components("osc_exp"), ROOT.RooFit.LineColor(ROOT.kRed))
    final.plotOn(frame, ROOT.RooFit.Components("osc_pdf"), ROOT.RooFit.LineColor(ROOT.kRed))
    frame.Draw()
    c1.Update()
    raw_input("E")

contour = minuit.contour(period, phase, 1, 2, 3, 0)
if contour:
    contour.Draw()
    c1.Update()
    raw_input("E")

for var in (period, phase):
    pll = nll.createProfile(ROOT.RooArgSet(var))
    frame = var.frame(ROOT.RooFit.Bins(5))#, ROOT.RooFit.Precision(1e-1))
    pll.plotOn(frame)
    frame.Draw()
    c1.Update()
    raw_input("E")
