import numpy
import ROOT
import sys
import cPickle as pickle
from pyWIMP.DMModels.base_model import BaseVariables
from pyWIMP.DMModels.wimp_model import WIMPModel
from pyWIMP.DMModels.low_energy_background import TestModel
import pyWIMP.Calculation.DataCalcVerification as dcv 


if len(sys.argv) != 2: 
    sys.exit(1)


#####################################
"""
Initialization stuff
"""
#####################################

#ROOT.RooMsgService.instance().setSilentMode(True)
#ROOT.RooMsgService.instance().setGlobalKillBelow(5)
efficiency = ROOT.RooRealVar("efficiency", "efficiency",
                             0, 1)
m = ROOT.RooRealVar("m","m", 100)

zero_offset = ROOT.RooRealVar("zero_offset", "zero_offset", 0)
background = ROOT.RooRealVar("background", "background", 0, 10)
tau = ROOT.RooRealVar("tau", "tau", 3.5)
signal = ROOT.RooRealVar("signal", "signal", 0, 10)

x = ROOT.RooRealVar("x", "x", 0, 10)
y = ROOT.RooRealVar("y", "y", 0, 10)

linear_var = ROOT.RooRealVar("pois_var", "pois_var", signal, efficiency, signal)
bkgd_linear_var = ROOT.RooRealVar("bkgd_var", "bkgd_var", background, tau, zero_offset)

bkg_pois = ROOT.RooPoisson("bkg_pois", "bkg_pois", y, bkgd_linear_var) 
pois = ROOT.RooPoisson("pois", "pois", x, linear_var) 

variables = ROOT.RooArgSet(basevars.get_energy())
# Following is the number of events
# Remember model_normal is in pb, so this is
# useful for setting things correctly later.
# This give us the expected events for a
# model_normal of 1
scaler = model_extend.expectedEvents(variables)
# Add the models together to get a final, extended model
i = 0
extended_models = []
while 1:
    amod = list_of_models.at(i)
    avar = list_of_coefficients.at(i)
    if not amod: break
    i += 1
    extend = ROOT.RooExtendPdf("extend%s" % amod.GetName(),
                               "extend%s" % amod.GetName(),
                               amod, avar)
    extended_models.append(extend)
temp_list = ROOT.RooArgList()
temp_list.add(model_extend)
for amod in extended_models:
    temp_list.add(amod)

fit_model = ROOT.RooAddPdf("b+s",
                           "Background + Signal",
                           temp_list)
test_variable = model_normal


calc_system = dcv.DataCalcVerification()

#####################################
"""
MPI stuff
"""
#####################################


number_of_iter = 10
step_size = float(exponential_total)/number_of_iter
m=[i*step_size for i in range(number_of_iter)]

c1 = ROOT.TCanvas()
calc_system.set_canvas(c1)
calc_system.set_debug(True)
results = []
for v in reversed(m):
    test_variable.setVal(v/scaler)
    exp_coef.setVal(exponential_total-v)
    results = calc_system.scan_confidence_value_space_for_model(
                      fit_model, 
                      test_variable,
                      variables,
                      total_entries,
                      1,
                      0.9)
    results = (v/scaler, exponential_total-v, results)
    print "Finished: ", comm.Get_rank()
                      
recvbuf = comm.gather(results,root=root)

print "Finishing"
afile = open('output_WM_%g.pkl' % wimp_mass, 'wb')
pickle.dump(recvbuf, afile)
afile.close()
