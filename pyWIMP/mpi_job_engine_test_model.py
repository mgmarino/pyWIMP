import numpy
import ROOT
ROOT.RooCurve() # Apparently important, this blows away 
                # some MPI things if it is instantiated later
ROOT.gROOT.SetBatch()

import sys
import cPickle as pickle
from mpi4py import MPI
from pyWIMP.DMModels.base_model import BaseVariables
from pyWIMP.DMModels.wimp_model import WIMPModel
from pyWIMP.DMModels.low_energy_background import TestModel
import pyWIMP.Calculation.DataCalcVerification as dcv 
from pyWIMP.DMModels.low_energy_background import LowEnergyBackgroundModel


if len(sys.argv) != 2: 
    sys.exit(1)

test_file = ROOT.TFile(sys.argv[1])
test_tree = test_file.Get("sensitivity_tree")
test_tree.GetEntry(0)

wimp_mass = test_tree.wimp_mass 

adict = {}
list_of_branches = test_tree.GetListOfBranches()
for i in range(list_of_branches.GetEntries()):
    branch = list_of_branches.At(i)
    try:
        temp = float(branch.GetName())
    except ValueError:
        continue
    adict[float(branch.GetName())] = (str(getattr(test_tree, branch.GetName()+'vars')),
                                      getattr(test_tree, branch.GetName()).minNll())

keys = adict.keys()
keys.sort()

var_string_list = []
array_list = []
for akey in keys:
    astring, minNLL = adict[akey]
    array_list.append((akey, minNLL))
    var_string_list.append(astring) 

# No get the index of the minimum ll above model_amplitude of 0
ll_array = numpy.array(array_list)
above_zero_entries = ll_array[:,0] >= 0
min_val = numpy.min(ll_array[above_zero_entries])
# Subtract the minimum
ll_array -= [0, min_val]
below_two_entries = ll_array[:,1] <= 2

minimum_entry = ll_array[:,1] == 0
not_minimum_entry = ll_array[:,1] != 0

scratch_array = ll_array[below_two_entries*above_zero_entries*not_minimum_entry]
minimum = ll_array[above_zero_entries*minimum_entry]
string_np = numpy.array(var_string_list)
string_np = numpy.concatenate((string_np[above_zero_entries*minimum_entry],
                              string_np[below_two_entries*above_zero_entries*not_minimum_entry]))

# Now string_np holds all the necessary values, with the first value being the best fit


#####################################
"""
Initialization stuff
"""
#####################################

ROOT.RooMsgService.instance().setSilentMode(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(5)
total_mc_entries = 500
#total_mc_entries = 10
total_entries = 400 

basevars = BaseVariables(0, 0.1444,0.5, 3.5) 
basevars.get_time().setConstant(True)

# Set up the WIMP class
wimp_class = WIMPModel(basevars, 
                       mass_of_wimp = wimp_mass,
                       kilograms=0.4,
                       constant_quenching=False)
model = wimp_class.get_model()
# Set up the background class
low_energy = LowEnergyBackgroundModel(basevars)
low_energy_model = low_energy.get_model()

list_of_models, list_of_coefficients = low_energy.get_list_components()
background_normal = ROOT.RooRealVar("flat_normal", 
                                    "Background event number", 
                                    0,
                                    total_entries)
background_extend = ROOT.RooExtendPdf("background_extend", 
                                                   "background_extend", 
                                                   low_energy_model, 
                                                   background_normal)


#flat_coef.setMin(-5)
#exp_coef.setMin(-5)
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
    print "Using number of nodes: ", comm.Get_size()

    number = comm.Get_size()-1
    step_size = int((len(string_np))/(number))
    # The last slice just prunes the end
    m =  string_np[0::step_size][:number]
    sendbuf = numpy.concatenate(([''], m))

print "Scatter for: ", comm.Get_rank()
var_cache=comm.scatter(sendbuf,root)
print "Go for: ", comm.Get_rank()
results = []
if comm.Get_rank() != 0:
    temp = ROOT.istringstream(var_cache)
    fit_model.getVariables().readFromStream(temp, False)
    results = calc_system.scan_confidence_value_space_for_model(
                      fit_model, 
                      test_variable,
                      variables,
                      total_entries,
                      total_mc_entries,
                      0.9)
results = (v, v/scaler, exponential_total-v, results)

print "Finished: ", comm.Get_rank()
recvbuf = comm.gather(results,root)
print "Sent: ", comm.Get_rank()

if comm.Get_rank()==0:
    print "Finishing"
    afile = open('output_WM_%g.pkl' % wimp_mass, 'wb')
    recvbuf = recvbuf[1:]
    pickle.dump(recvbuf, afile)
    afile.close()
