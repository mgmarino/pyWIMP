import numpy
import ROOT
import sys
import cPickle as pickle
from pyWIMP.DMModels.base_model import BaseVariables
from pyWIMP.DMModels.wimp_model import WIMPModel
from pyWIMP.DMModels.low_energy_background import TestModel
import pyWIMP.Calculation.DataCalcVerification as dcv 
import math



#####################################
"""
Initialization stuff
"""
#####################################

ROOT.RooMsgService.instance().setSilentMode(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(5)
efficiency = ROOT.RooRealVar("efficiency", "efficiency",
                             0.85, 0, 1)
efficiency_sigma = ROOT.RooRealVar("efficiency_sigma", "efficiency_sigma",
                             0.075, 0, 1)

background = ROOT.RooRealVar("background", "background", 0, 20)
background_sigma = ROOT.RooRealVar("background_sigma", 
                                   "background_sigma",
                                   0.075, 0, 1)
signal = ROOT.RooRealVar("signal", "signal", -2, 10)

x = ROOT.RooRealVar("x", "x", 0, 20)
x.setBins(int(x.getMax()-x.getMin()))
y = ROOT.RooRealVar("y", "y", 0, 20)
#x.setBins(int(x.getMax()-x.getMin()))
z = ROOT.RooRealVar("z", "z", 0, 1)

linear_var = ROOT.RooLinearVar("pois_var", "pois_var", signal, efficiency, background)

bkg_gaus = ROOT.RooGaussian("bkg_gaus", "bkg_gaus", y, background, background_sigma) 
eff_gaus = ROOT.RooGaussian("eff_gaus", "eff_gaus", z, efficiency, efficiency_sigma) 
pois = ROOT.RooPoisson("pois", "pois", x, linear_var) 

fit_model = ROOT.RooProdPdf("fit_model", "fit_model", 
                            ROOT.RooArgList(pois, bkg_gaus, eff_gaus))
                            

list_of_everything = [ efficiency, efficiency_sigma, 
                       background, background_sigma,
                       signal, x, y, z, linear_var,
                       bkg_gaus, eff_gaus, pois, 
                       fit_model ] 

variables = ROOT.RooArgSet(x, y, z)
# Following is the number of events
# Remember model_normal is in pb, so this is
# useful for setting things correctly later.
# This give us the expected events for a
# model_normal of 1
test_variable = signal

calc_system = dcv.DataCalcVerification()

#####################################
"""
MPI stuff
"""
#####################################


number_of_jobs = 64

number_of_iter = math.sqrt(number_of_jobs)
step_size = float(10)/(number_of_iter-1)
m= numpy.array([[(i*step_size, j*step_size) 
       for i in range(number_of_iter)] 
       for j in range(number_of_iter)])

m.shape=(number_of_jobs, 2)

c1 = ROOT.TCanvas()
#calc_system.set_canvas(c1)
#calc_system.set_debug(True)
results = []
for v, bgd in reversed(m):
    test_variable.setVal(v)
    background.setVal(bgd)
    results = calc_system.scan_confidence_value_space_for_model(
                      fit_model, 
                      test_variable,
                      variables,
                      100,
                      1,
                      0.9)
    #bf, upper, lower, bound, an_array = results[0]
    bf, upper, lower, bound = results[0]
    """
    curve = ROOT.RooCurve()
    [curve.addPoint(x,y) for x, y in an_array]
    curve.Draw("APL")
    hist = curve.GetHistogram()
    hist.SetMaximum(2)
    hist.SetMinimum(-0.5)
    hist.Draw()
    curve.Draw("L")
    c1.Update()
    """
    print bf, upper, lower, bound, v
    raw_input("E")
    #results = (v, 10-v, results)
    #print "Finished: ", comm.Get_rank()
                      
#recvbuf = comm.gather(results,root=root)

#print "Finishing"
#afile = open('output_WM_%g.pkl' % wimp_mass, 'wb')
#pickle.dump(recvbuf, afile)
#afile.close()
