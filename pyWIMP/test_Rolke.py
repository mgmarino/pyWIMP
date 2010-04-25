import numpy
import ROOT
ROOT.RooCurve() # Apparently important, this blows away 
                # some MPI things if it is instantiated later
import sys
import cPickle as pickle
from mpi4py import MPI
import pyWIMP.Calculation.DataCalcVerification as dcv 
import math


job_number = int(sys.argv[1])

#####################################
"""
Initialization stuff
"""
#####################################

total_MC_iterations = 1500

ROOT.RooRandom.randomGenerator().SetSeed(0)
ROOT.RooMsgService.instance().setSilentMode(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(5)

efficiency = ROOT.RooRealVar("efficiency", 
                             "efficiency",
                             0.85, 0, 1)
efficiency_sigma = ROOT.RooRealVar("efficiency_sigma", 
                                   "efficiency_sigma",
                                    0.075, 0, 1)

background = ROOT.RooRealVar("background", 
                             "background", 
                             -2, 20)

background_sigma = ROOT.RooRealVar("background_sigma", 
                                   "background_sigma",
                                   0.075, 0, 1)

# Signal, test variable
signal = ROOT.RooRealVar("signal", 
                         "signal", 
                         -2, 20)

x = ROOT.RooRealVar("x", "x", -2, 20)
x.setBins(int(x.getMax()-x.getMin()))
y = ROOT.RooRealVar("y", "y", -2, 20)
z = ROOT.RooRealVar("z", "z", 0, 1)

linear_var = ROOT.RooLinearVar("pois_var", 
                               "pois_var", 
                               signal, 
                               efficiency, 
                               background)

bkg_gaus = ROOT.RooGaussian("bkg_gaus", 
                            "bkg_gaus", 
                            y, 
                            background, 
                            background_sigma) 

eff_gaus = ROOT.RooGaussian("eff_gaus", 
                            "eff_gaus", 
                            z, 
                            efficiency, 
                            efficiency_sigma) 

pois = ROOT.RooPoisson(    "pois", 
                        "pois", 
                        x, 
                        linear_var) 

fit_model = ROOT.RooProdPdf("fit_model", "fit_model", 
                            ROOT.RooArgList(pois, bkg_gaus, eff_gaus))
                            

# This shouldn't be necessary, but for some reason some items go out of
# scope.  This keeps this from happening.
list_of_everything = [ efficiency, efficiency_sigma, 
                       background, background_sigma,
                       signal, x, y, z, linear_var,
                       bkg_gaus, eff_gaus, pois, 
                       fit_model ] 

variables = ROOT.RooArgSet(x, y, z)

# Test variable is the signal
test_variable = signal

calc_system = dcv.DataCalcVerification()

#####################################
"""
MPI stuff
"""
#####################################
comm = MPI.COMM_WORLD
sendbuf=[]
root=0

if comm.Get_rank()==0:
    print("Using number of nodes: ", comm.Get_size())
    number_of_jobs = comm.Get_size() - 1
    number_of_iter = int(math.sqrt(number_of_jobs))
    step_size = float(10)/(number_of_iter-1)
    m= numpy.array([[(i*step_size, j*step_size) 
           for i in range(number_of_iter)] 
           for j in range(number_of_iter)])
    m.shape=(number_of_jobs, 2)
    test = numpy.zeros((1,2))
    sendbuf = numpy.concatenate((test, m))

v = comm.scatter(sendobj=sendbuf,root=root)
results = []
if comm.Get_rank() != 0:
    test_value, bgd_value = v
    print comm.Get_rank(), ": ", test_value, bgd_value
    test_variable.setVal(test_value)
    background.setVal(bgd_value)

    results = calc_system.scan_confidence_value_space_for_model(
                      fit_model, 
                      test_variable,
                      variables,
                      100,
                      total_MC_iterations,
                      0.9)
    results = (test_value, bgd_value, results)
    print "Finished: ", comm.Get_rank()
                      
recvbuf = comm.gather(results,root=root)

if comm.Get_rank()==0:
    recvbuf = recvbuf[1:]
    print "Finishing"
    afile = open('output_Rolke_Test_%i.pkl' % job_number, 'wb')
    pickle.dump(recvbuf, afile)
    afile.close()
