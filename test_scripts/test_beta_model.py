import ROOT
import pyWIMP.DMModels.base_model as base_model
import sys
import pyWIMP.WIMPPdfs as pdfs  
from pyWIMP.utilities.utilities import rescale_frame
import math


exposure_time_in_days = 15
time_in_days = 3*365
kg = 30 
events = (200*exposure_time_in_days*kg*(1 - math.exp(-math.log(2)*(time_in_days)/(12.36*365.25))))

basevars = base_model.BaseVariables(0,1./365., 0, 20)
time = basevars.get_time()
energy = basevars.get_energy()
q_value = ROOT.RooRealVar("q_value", "q_value", 18.6)
mass_of_electron = ROOT.RooRealVar("mass_of_electron", "mass_of_electron", 511)

model_amp = ROOT.RooRealVar("amp", "amp", events, 1e-15, 100000)
model = pdfs.MGMBetaDecayFunction("beta", "beta", energy, mass_of_electron, q_value)

model_extend = ROOT.RooExtendPdf("extend", "extend", model, model_amp)

data = model_extend.generate(ROOT.RooArgSet(energy))

model_extend.fitTo(data)

c1 = ROOT.TCanvas()
frame = energy.frame()
data.plotOn(frame)
model_extend.plotOn(frame)
frame.Draw()
bin_width = frame.getFitRangeBinW()
axis = rescale_frame(c1, frame, 1./(bin_width*time_in_days*kg), "Events/keV/kg/d") 
c1.Update()
raw_input("E")


