import numpy
import ROOT
ROOT.RooCurve() # Apparently important, this blows away 
                # some MPI things if it is instantiated later

import sys
import cPickle as pickle
from mpi4py import MPI
from pyWIMP.DMModels.base_model import BaseVariables
from pyWIMP.DMModels.wimp_model import WIMPModel
from pyWIMP.DMModels.low_energy_background import TestModel
import pyWIMP.Calculation.DataCalcVerification as dcv 


if len(sys.argv) != 2: 
    sys.exit(1)

wimp_mass = float(sys.argv[1])

#####################################
"""
Initialization stuff
"""
#####################################

ROOT.RooMsgService.instance().setSilentMode(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(5)
total_mc_entries = 500
total_entries = 400 
exponential_total = 190
basevars = BaseVariables(0, 0.1444,0.5, 3.5) 
basevars.get_time().setConstant(True)

# Set up the WIMP class
wimp_class = WIMPModel(basevars, 
                       mass_of_wimp = wimp_mass,
                       kilograms=0.4,
                       constant_quenching=False)
model = wimp_class.get_model()

# Set up the background class
background = TestModel(basevars)
list_of_models, list_of_coefficients = \
                          background.get_list_components()

exp_coef = list_of_coefficients.at(
             list_of_coefficients.index("exp_coef_"))
flat_coef = list_of_coefficients.at(
              list_of_coefficients.index("flat_coef_"))
exp_coef.setVal(exponential_total)
flat_coef.setVal(180)
flat_coef.setMin(-5)
exp_coef.setMin(-5)
# Now set up the extended model
model_normal = ROOT.RooRealVar("model_normal",
                               "WIMP-nucleus xs",
                               1, -10, total_entries,
                               "pb")

model_extend = ROOT.RooExtendPdf("model_extend",
                                 "model_extend",
                                 model,
                                 model_normal)

variables = ROOT.RooArgSet(basevars.get_energy())
# Following is the number of events
# Remember model_normal is in pb, so this is
# useful for setting things correctly later.
# This give us the expected events for a
# model_normal of 1
scaler = model_extend.expectedEvents(variables)
model_normal.setMax(total_entries/scaler)
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
comm = MPI.COMM_WORLD
sendbuf=[]
root=0


if comm.Get_rank()==0:
    print("Using number of nodes: ", comm.Get_size())
    step_size = float(exponential_total)/(comm.Get_size()-1)
    m=[(i-1)*step_size if i > 0 else 0 for i in range(comm.Get_size())]
    sendbuf=m

v=comm.scatter(sendobj=sendbuf,root=root)

results = []
if comm.Get_rank() != 0:
    test_variable.setVal(v/scaler)
    exp_coef.setVal(exponential_total-v)
    results = calc_system.scan_confidence_value_space_for_model(
                      fit_model, 
                      test_variable,
                      variables,
                      total_entries,
                      total_mc_entries,
                      0.9)
    results = (v, v/scaler, exponential_total-v, results)
    print "Finished: ", comm.Get_rank()
                      
recvbuf = comm.gather(results,root=root)

if comm.Get_rank()==0:
    print "Finishing"
    afile = open('output_WM_%g.pkl' % wimp_mass, 'wb')
    pickle.dump(recvbuf, afile)
    afile.close()
