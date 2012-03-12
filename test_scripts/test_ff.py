import pyWIMP.DMModels.wimp_model as wimp_model
from pyWIMP.DMModels.base_model import BaseVariables
import ROOT
import re
#ROOT.gROOT.SetBatch()
basevars = BaseVariables(0, 1, 0, 3.5)
time = basevars.get_time()
time.setConstant()
time.setVal(0)
mass_of_wimp = 10
wm = wimp_model.WIMPModel(basevars, mass_of_wimp)
c1 = ROOT.TCanvas()
c1.SetLogy()
#test = ROOT.MGMWimpHelmFFSquared("model", "model", q, r, s)
#energy = q
test = wm.get_helm_form_factor()
test1 = wm.get_WIMP_model()
test2 = wm.get_WIMP_model_with_escape_vel()
test3 = wm.get_WIMP_model_with_escape_vel_no_ff()
energy = basevars.get_energy()
frame = energy.frame()
energy.setVal(0)
print test1.getVal(), test2.getVal(), test3.getVal(), test.getVal()
variables = ROOT.RooArgSet(basevars.get_energy())

integral = test1.createIntegral(variables)
orig = integral.getVal()

integral = test.createIntegral(variables)
ff_int = integral.getVal()

integral = test2.createIntegral(variables)
orig_esc_int = integral.getVal()

integral = test3.createIntegral(variables)
simple_int = integral.getVal()

print orig, orig_esc_int, simple_int, ff_int
test1_var =  test1.expectedEvents(variables)
test2_var =  test2.expectedEvents(variables)
test3_var =  test3.expectedEvents(variables)
print test1_var, test2_var, test3_var
test2.SetTitle("v_{esc} = 600 km s^{-1}")
test1.SetTitle("v_{esc} #rightarrow #infty")
test3.SetTitle("v_{esc} = 600 km s^{-1}, FF #rightarrow 1  ")
test2.plotOn(frame, ROOT.RooFit.LineStyle(9),ROOT.RooFit.FillColor(10),\
                    ROOT.RooFit.Normalization(test2_var, ROOT.RooAbsReal.Raw))
test1.plotOn(frame, ROOT.RooFit.LineStyle(3),ROOT.RooFit.FillColor(10),\
                    ROOT.RooFit.Normalization(test1_var, ROOT.RooAbsReal.Raw))
test3.plotOn(frame, ROOT.RooFit.LineStyle(1),ROOT.RooFit.FillColor(10),\
                    ROOT.RooFit.Normalization(test3_var, ROOT.RooAbsReal.Raw))
#together_again.plotOn(frame,\
#                    ROOT.RooFit.Normalization(orig_time_ff, ROOT.RooAbsReal.Raw),\
#                    ROOT.RooFit.Precision(1e-10))
frame.SetMaximum(10000)
frame.SetMinimum(1e-1)
frame.GetXaxis().SetTitle("Energy (keV)")
frame.GetXaxis().CenterTitle()
frame.GetYaxis().SetTitle("#frac{dR}{dE} #sigma^{-1} (counts/keV/kg/yr/pb)")
frame.GetYaxis().SetTitleOffset(1.35)
frame.GetYaxis().CenterTitle()
frame.SetTitle("WIMP Mass: %g GeV" % mass_of_wimp)
ROOT.gStyle.SetOptTitle(1)
c1.SetLeftMargin(0.13)
frame.Draw()
alist = c1.GetListOfPrimitives()
legend = ROOT.TLegend(0.45, 0.6, 0.88, 0.88)
for item in alist:
    match = re.match("Projection of (.*)", item.GetTitle())
    if not match: continue
    item.SetTitle(match.group(1))
    legend.AddEntry(item)
c1.SetGrid(1,1)
legend.Draw()
c1.Update()
raw_input("Enter")
c1.Print("WIMPMass10GeVExample.eps")
ROOT.gStyle.SetOptTitle(0)
energy.setMax(20)
frame = energy.frame() 
test.plotOn(frame, ROOT.RooFit.Normalization(1, ROOT.RooAbsReal.Raw),
            ROOT.RooFit.Precision(1e-10))
frame.SetMinimum(7e-2)
frame.SetMaximum(4)
yaxis = frame.GetYaxis()
yaxis.SetTitle("Helm FF^{2}")
yaxis.SetTitleOffset(1.3)
#yaxis.CenterTitle()
frame.GetXaxis().SetTitle("Ionization energy (keV)")

frame.Draw()
c1.SetGrid(1,1)
c1.Update()
raw_input("Enter")
c1.Print("HelmFormFactor.eps")


#c1.Print("wimpMass%i.root" % mass_of_wimp)
