#!/usr/bin/env python

from pyWIMP.DMModels.base_model import BaseVariables
import pyWIMP.DMModels.wimp_model as wimp_model
import ROOT
import re
mass_of_wimp = 7
basevars = BaseVariables(0.0, 4.0, 0, 10, True, 0.0)  
time = basevars.get_time()
time.setConstant()
time.setVal(0)
wm = wimp_model.WIMPModel(basevars, mass_of_wimp,  nucl_recoil=True)
model = wm.get_WIMP_model()
expected_events = model.expectedEvents(ROOT.RooArgSet(basevars.get_energy()))
hist = model.createHistogram("hist", basevars.get_energy(), \
                          ROOT.RooFit.Binning(200),\
                          ROOT.RooFit.YVar(basevars.get_time(), \
                          ROOT.RooFit.Binning(4*36)))
hist.SetLineColor(4)
hist.GetXaxis().CenterTitle()
hist.GetXaxis().SetTitle("Energy (keVnr)")
hist.GetXaxis().SetTitleOffset(1.5)
hist.GetYaxis().CenterTitle()
hist.GetYaxis().SetTitle("Time (years)")
hist.GetYaxis().SetTitleOffset(1.5)
hist.GetYaxis().SetNdivisions(507)
hist.GetZaxis().CenterTitle()
hist.GetZaxis().SetTitle("Amplitude (a.u.)")
hist.GetZaxis().SetTitle("#frac{dR}{dE} #sigma^{-1} (counts/keV/kg/yr/pb)")
hist.GetZaxis().SetTitleOffset(1.1)
hist.GetZaxis().SetNdivisions(509)
ysize = hist.GetYaxis().GetBinWidth(1)
xsize = hist.GetXaxis().GetBinWidth(1)
time = basevars.get_time()
time_scale = time.getMax() - time.getMin()
print ysize, xsize
hist.Scale(expected_events/(xsize*ysize))
c1 = ROOT.TCanvas()
c1.SetLogz()
c1.SetTheta(22.1822);
c1.SetPhi(-24.31034);
hist.Draw("SURF4FB")
hist.GetZaxis().SetRangeUser(0.05, 1000)
c1.Update()
rootFile = raw_input("Enter output ROOT file name: ")
f = ROOT.TFile(rootFile,'recreate')
name = raw_input("Enter name for hist: ")
hist.SetName(name)
hist.Write()

c1.Print(name+".eps")
