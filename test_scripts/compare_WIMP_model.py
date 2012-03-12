#!/usr/bin/env python

import pyWIMP.DMModels.wimp_model as wimp_model
from pyWIMP.DMModels.base_model import BaseVariables
import ROOT
import re
#ROOT.gROOT.SetBatch()
c1 = ROOT.TCanvas()
c1.SetLogy()
basevars = BaseVariables(0, 1, 0, 4.5)
time = basevars.get_time()
energy = basevars.get_energy()
time.setVal(0)
time.setConstant()
mass_of_wimp = 15
wimp_models = []
pdf_list = []
#mass_list = range(8,11)
#mass_list.extend([100])
mass_list = [10.25]
for mass_of_wimp in mass_list: 
    wm = wimp_model.WIMPModel(basevars, mass_of_wimp)
    
    #extend = ROOT.RooExtendPdf("extend_%s" % wm.get_model().GetName(),
    #                           "extend_%s" % wm.get_model().GetName(),
                               
    pdf_list.append((wm.get_model().createIntegral(ROOT.RooArgSet(energy)).getVal(), 
                     wm.get_model())) 
    wm.get_model().SetTitle("WIMP Mass: %g GeV " % mass_of_wimp)
    wimp_models.append(wm)
expo_const = ROOT.RooRealVar("expo_const", "expo_const", -3.3)
expo_model = ROOT.RooExponential("expo_model", "Exponential, e^{-3.3E}", energy, expo_const)
flat_model = ROOT.RooPolynomial("flat_model", "Flat", energy)
#pdf_list.append((500, expo_model))
#pdf_list.append((2.5*(energy.getMax() - energy.getMin()), flat_model))


frame = energy.frame()
energy.setVal(0)
line_style = 0
line_list = [1, 2, 3, 9]
for norm, test1 in pdf_list:
    print norm
    print 
    test1.plotOn(frame, ROOT.RooFit.LineStyle(line_list[line_style % len(line_list)]),ROOT.RooFit.FillColor(10),\
                    ROOT.RooFit.Normalization(norm, ROOT.RooAbsReal.Raw))
    line_style += 1
#together_again.plotOn(frame,\
#                    ROOT.RooFit.Normalization(orig_time_ff, ROOT.RooAbsReal.Raw),\
#                    ROOT.RooFit.Precision(1e-10))
frame.SetMaximum(3000)
frame.SetMinimum(0.01)
frame.GetXaxis().SetTitle("Energy (keVee)")
frame.GetXaxis().CenterTitle()
frame.GetYaxis().SetTitle("#frac{dR}{dE} #sigma^{-1} (counts/keV/kg/yr/pb)")
frame.GetYaxis().SetTitleOffset(1.35)
frame.GetYaxis().CenterTitle()
frame.SetTitle("WIMP Mass: %g GeV" % mass_of_wimp)
c1.SetLeftMargin(0.13)
frame.Draw()
alist = c1.GetListOfPrimitives()
legend = ROOT.TLegend(0.4, 0.6, 0.88, 0.88)
for item in alist:
    match = re.match("Projection of (.*)", item.GetTitle())
    if not match: continue
    item.SetTitle(match.group(1))
    legend.AddEntry(item)
c1.SetGrid(1,1)
legend.Draw()
c1.Update()
raw_input("Enter")
c1.Print("WIMPModelCompare.eps")

