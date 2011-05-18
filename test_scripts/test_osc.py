import pyWIMP.DMModels.base_model as base_model
import pyWIMP.WIMPPdfs as pdfs
import ROOT
import sys

#ROOT.gROOT.SetBatch()
basevars = base_model.BaseVariables(0, 1, 0, 13)
time = basevars.get_time()
energy = basevars.get_energy()
livetime = pdfs.MGMPiecewiseRegions()
livetime.InsertNewRegion(0, 1)
time_exp = ROOT.RooRealVar("time_exp", "time_exp", -0.20, -1, 0.2)
tag = "low"
exp_constant_one = ROOT.RooRealVar("expo_const_one%s" % tag,
                                   "expo_const_one%s" % tag,
                                   -3., -5, -0.8)
exp_constant_one.setError(0.5)
exp_coef = ROOT.RooRealVar("exp_coef_%s" % tag,
                           "exp_coef_%s" % tag,
                           1000.2, 10, 2000)

# Flat pdf
phase    = ROOT.RooRealVar("phase", "phase", 0, -0.5, 0.5)#-ROOT.TMath.Pi(), ROOT.TMath.Pi()) 
osc_ampl = ROOT.RooRealVar("osc_amplitude", "osc_amplitude", 0.2, 0, 1.) 
period = ROOT.RooRealVar("osc_period", "osc_period", 1, 0.5, 1.5) 
osc_time_pdf = ROOT.MGMExponentialPlusSinusoid("osc_pdf", "osc_pdf", 
                                               basevars.get_time(), time_exp, 
                                               osc_ampl, 
                                               period, 
                                               phase) 
osc_time_pdf.SetRegionsOfValidity(livetime)


energy_exp_pdf = ROOT.RooExponential("energy_pdf_exp", 
                                     "energy_pdf_exp", 
                                     basevars.get_energy(),
                                     exp_constant_one)

exp_pdf_whole = ROOT.RooProdPdf("exp_pdf_whole", "exp_pdf_whole", energy_exp_pdf, osc_time_pdf)


extend = ROOT.RooExtendPdf("extend", "extend", exp_pdf_whole, exp_coef)

set = ROOT.RooArgSet(energy, time)
data = extend.generate(set, -1)
data.Print('v')

nll = extend.createNLL(data)
minuit = ROOT.RooMinuit(nll)
minuit.migrad()
minuit.minos(ROOT.RooArgSet(phase, osc_ampl, period))

c1 = ROOT.TCanvas()
for var in (energy, time):
    frame = var.frame()
    data.plotOn(frame)
    extend.plotOn(frame)
    frame.Draw()
    c1.Update()
    raw_input("E")

for var in (period, phase):
    pll = nll.createProfile(ROOT.RooArgSet(var))
    frame = var.frame(ROOT.RooFit.Bins(5))#, ROOT.RooFit.Precision(1e-1))
    pll.plotOn(frame)
    frame.Draw()
    c1.Update()
    raw_input("E")
