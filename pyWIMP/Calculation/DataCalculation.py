import ExclusionCalculation
import ROOT
import os
import math
from exceptions import Exception
from ..utilities.utilities import rescale_frame
import numpy
class DataCalculation(ExclusionCalculation.ExclusionCalculation):

    """
    This class handles calculating the exclusion limits given a 
    certain input data set.  We have to overload 
    scan_confidence_value_space_for_model to do this.
    """

    def print_plot(self, model, data, title = "", scaling = 1.):
        var_iter = model.getObservables(data).createIterator()
        while 1:
            var_obj = var_iter.Next()
            if not var_obj: break
            aframe = var_obj.frame().emptyClone('temp')
            ROOT.RooAbsData.plotOn(data, aframe)
            model.plotOn(aframe)
            model.plotOn(aframe, 
                 ROOT.RooFit.Components("WIMPPDF_With_Time_And_Escape_Vel"), 
                 ROOT.RooFit.LineStyle(ROOT.RooFit.kDashed))
            model.plotOn(aframe, 
                 ROOT.RooFit.Components("energy_pdf_*"), 
                 ROOT.RooFit.LineWidth(4),
                 ROOT.RooFit.LineStyle(ROOT.RooFit.kDotted),
                 ROOT.RooFit.LineColor(ROOT.RooFit.kRed))
            model.plotOn(aframe, 
                 ROOT.RooFit.Components("gamma*"), 
                 ROOT.RooFit.LineWidth(4),
                 ROOT.RooFit.LineColor(ROOT.RooFit.kRed))
            aframe.SetTitle("%s %s" % 
                            (self.plot_base_name, title))
            #bin_width = aframe.getFitRangeBinW()
            #axis = rescale_frame(self.c1, aframe, scaling/bin_width, axis_title)
            #axis.CenterTitle()
            aframe.Draw()
            self.c1.Update()
            if self.print_out_plots:
                title = aframe.GetTitle()
                title = title.replace(' ','').replace('(','').replace(')','').replace(',','') 
                self.c1.Print(title + ("%s.eps" % var_obj.GetName()))
            else:
                raw_input("Hit Enter to continue")


    def find_confidence_value_for_model(self, 
                                        model, 
                                        data, 
                                        model_amplitude, 
                                        conf_level, 
                                        mult_factor, 
                                        print_level = -1, 
                                        verbose = False, 
                                        debug = False,
                                        tolerance = 0.0001):
    
        number_of_points = 100
        distance_from_min = 20.
        #pars = model.getParameters(data)
        nll = model.createNLL(data, ROOT.RooFit.Verbose(verbose))
       
        #model_amplitude.setVal(model_amplitude.getMin())
        model_amplitude.setVal(0)
        model_amplitude.setConstant(True)
        minuit = ROOT.RooMinuit(nll)
        minuit.migrad()

        # Now fit with the model_amplitude
        model_amplitude.setVal(0)
        model_amplitude.setConstant(False)
        minuit.migrad()
        min_nll = nll.getVal()
        best_fit_value = model_amplitude.getVal()

        if conf_level == 0: return None 

        min_value = model_amplitude.getVal() - distance_from_min 
        if min_value > 0: min_value = 0
        if min_value < model_amplitude.getMin(): min_val = model_amplitude.getMin()
        
        max_range = model_amplitude.getVal() + 50
        if max_range < 10: max_range = 50


        # Ensure that we are going 

        model_amplitude.setConstant(True)
        model_amplitude.setVal(max_range)
        minuit.migrad()
        while nll.getVal()-min_nll < 2*conf_level: 
            max_range += 100
            model_amplitude.setVal(max_range)
            if model_amplitude.getVal() == model_amplitude.getMax():
                self.logging("Resetting maximum:", model_amplitude.getMax() )
                model_amplitude.setMax(model_amplitude.getVal()*2)
            minuit.migrad()

        model_amplitude.setVal(0)
        model_amplitude.setConstant(False)
        minuit.migrad()

        first_res = minuit.save()

        output_list = []
        pll_curve = ROOT.RooCurve()
        pll_curve.SetName("pll_frac_plot") 

        orig = first_res.minNll() 
        output_dict = {}
        
        model_amplitude.setConstant(True)
        for i in numpy.arange(min_value, max_range + 1., (max_range - min_value)/number_of_points):
            print "Perfoming: ", i 
            model_amplitude.setVal(i)
            #j = 0
            #while 1:
            minuit.migrad()
            res = minuit.save(str(i)) 
            #expo_coef.setConstant(False) 
            #if j == 1: break
            #if math.fabs(res.correlation('exp_coef_', 'expo_const_one')) > 0.98:
            #    expo_coef.setConstant(True) 
            #    print "Set constant."
            #else: break
            #j += 1
                
            if debug:
                self.print_plot(model, data, str(i))
            output_list.append((i, res))
            if res.minNll() < orig:  orig = res.minNll()
            output_dict[str(i)] = res 
            
        [pll_curve.addPoint(i, res.minNll() - orig) for i, res in output_list]
        output_dict['first_fit'] = first_res
        output_dict['pll_curve'] = pll_curve
        first_res.Print('v')
       
        return output_dict
 
    def scan_confidence_value_space_for_model(self, 
                                              model, 
                                              data_model, 
                                              model_amplitude, 
                                              mult_factor,
                                              variables, 
                                              number_of_events,
                                              number_iterations, 
                                              cl):
    
        print_level = -1
        verbose = self.debug
        if self.debug: print_level = 3
        show_plots = self.show_plots or self.print_out_plots

        list_of_values = []
        i = 0
        confidence_value = ROOT.TMath.ChisquareQuantile(cl, 1) 
        if self.debug:
            self.logging( "Process %s: Fitting to data." )
        # Generate the data, use Extended flag
        # because the number_of_events is just
        # an expected number.
        model_amplitude.setVal(0)
        data_set_func = data_model
        get_val = self.find_confidence_value_for_model(
            model, 
            data_set_func, 
            model_amplitude, 
            0/2,
            mult_factor,
            print_level,
            verbose,
            show_plots) 

        # Looking for whether or not we have info 
        scaling = 1.
        axis_title = "Counts/keV"
        if "mass_of_detector" in self.input_variables.keys()\
           and "total_time" in self.input_variables.keys():
            kilos = self.input_variables["mass_of_detector"]
            time_in_years = self.input_variables["total_time"]
            scaling = 1./(kilos*time_in_years*365.25)
            axis_title = "Counts/keV/kg/d"
        if show_plots:
            var_iter = model.getObservables(data_set_func).createIterator()
            while 1:
                var_obj = var_iter.Next()
                if not var_obj: break
                aframe = var_obj.frame()
                ROOT.RooAbsData.plotOn(data_set_func, aframe)
                model.plotOn(aframe)
                model.plotOn(aframe, 
                     ROOT.RooFit.Components("WIMPPDF_With_Time_And_Escape_Vel"), 
                     ROOT.RooFit.LineStyle(ROOT.RooFit.kDashed))
                model.plotOn(aframe, 
                     ROOT.RooFit.Components("energy_pdf_*"), 
                     ROOT.RooFit.LineWidth(4),
                     ROOT.RooFit.LineStyle(ROOT.RooFit.kDotted),
                     ROOT.RooFit.LineColor(ROOT.RooFit.kRed))
                model.plotOn(aframe, 
                     ROOT.RooFit.Components("gamma*"), 
                     ROOT.RooFit.LineWidth(4),
                     ROOT.RooFit.LineColor(ROOT.RooFit.kRed))
                aframe.SetTitle("%s (Initial fit)" % self.plot_base_name)
                bin_width = aframe.getFitRangeBinW()
                axis = rescale_frame(self.c1, aframe, scaling/bin_width, axis_title)
                axis.CenterTitle()
                self.c1.Update()
                if self.print_out_plots:
                    title = aframe.GetTitle()
                    title = title.replace(' ','').replace('(','').replace(')','').replace(',','') 
                    self.c1.Print(title + ("%s.eps" % var_obj.GetName()))
                else:
                    raw_input("Hit Enter to continue")
          

        # Perform the fit and find the limits
        get_val = None
        while 1:
            get_val = self.find_confidence_value_for_model(
                model, 
                data_set_func, 
                model_amplitude, 
                confidence_value/2,
                mult_factor,
                print_level,
                verbose,
                show_plots) 
            
            if not get_val: 
                # There was an error somewhere downstream
                # or an interrupt was signalled
                # Get out
                break
            elif get_val == self.retry_error: 
                # Calling function requested a retry
                continue
            else: break
    
            # Store the results
        list_of_values.append(get_val)
    
        if show_plots:
            var_iter = model.getObservables(data_set_func).createIterator()
            while 1:
                var_obj = var_iter.Next()
                if not var_obj: break
                aframe = var_obj.frame()
                ROOT.RooAbsData.plotOn(data_set_func, aframe)
                model.plotOn(aframe)
                model.plotOn(aframe, 
                     ROOT.RooFit.Components("WIMPPDF_With_Time_And_Escape_Vel"), 
                     ROOT.RooFit.LineStyle(ROOT.RooFit.kDashed))
                model.plotOn(aframe, 
                     ROOT.RooFit.Components("energy_pdf_*"), 
                     ROOT.RooFit.LineWidth(4),
                     ROOT.RooFit.LineStyle(ROOT.RooFit.kDotted),
                     ROOT.RooFit.LineColor(ROOT.RooFit.kRed))
                model.plotOn(aframe, 
                     ROOT.RooFit.Components("gamma*"), 
                     ROOT.RooFit.LineWidth(4),
                     ROOT.RooFit.LineColor(ROOT.RooFit.kRed))
                aframe.SetTitle("%s (Final fit, %g CL, #sigma: %g pb)" % 
                                (self.plot_base_name, cl, 
                                 mult_factor*model_amplitude.getVal()))
                bin_width = aframe.getFitRangeBinW()
                axis = rescale_frame(self.c1, aframe, scaling/bin_width, axis_title)
                axis.CenterTitle()
                self.c1.Update()
                if self.print_out_plots:
                    title = aframe.GetTitle()
                    title = title.replace(' ','').replace('(','').replace(')','').replace(',','') 
                    self.c1.Print(title + ("%s.eps" % var_obj.GetName()))
                else:
                    raw_input("Hit Enter to continue")
                
        return list_of_values

