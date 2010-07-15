import ROOT
import os
import sys
from exceptions import Exception
from ..utilities.utilities import rescale_frame
class BaseCalculation:

    def __init__(self, exit_manager = None):
        self.exit_manager = exit_manager
        self.retry_error = {'again' : True} 
        self.debug = False
        self.show_plots = False
        self.plot_base_name = ""
        self.input_variables = {}
        self.print_out_plots = False

    def is_exit_requested(self):
        if not self.exit_manager: return False
        return self.exit_manager.is_exit_requested()

    def set_canvas(self, canv): self.c1 = canv
    def set_debug(self, set_d = True): self.debug = set_d 
    def set_print_out_plots(self, set_p = True): 
        self.print_out_plots = set_p 
    def set_show_plots(self, set_p = True): 
        self.show_plots = set_p 
    def set_plot_base_name(self, name): 
        self.plot_base_name = name 

    def set_input_variables(self, input):
        self.input_variables = input

    def logging(self, *strings_to_log): 
        header = "PID(%i), Name(%s): " % (os.getpid(), self.plot_base_name) 
        for astr in strings_to_log: header += " %s" % str(astr)
        sys.stdout.write(header)
        sys.stdout.write('\n')
        sys.stdout.flush()

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
                                        tolerance = 0.001):
 
        print "Base Class: BaseCalculation being called."
        return None

    def scan_confidence_value_space_for_model(self, 
                                              model, 
                                              data_model, 
                                              model_amplitude, 
                                              mult_factor,
                                              variables, 
                                              do_bin_data,
                                              number_iterations, 
                                              cl):
    
        print_level = -1
        verbose = self.debug
        if self.debug: print_level = 3

        list_of_values = []
        i = 0
        confidence_value = ROOT.TMath.ChisquareQuantile(cl, 1) 
        #ROOT.RooTrace.active(True)

        # Save the values of the variables to reset after fitting
        var_cache = ROOT.ostringstream()
        data_model.getVariables().writeToStream(var_cache, False)
        # Looking for whether or not we have info 
        scaling = 1.
        axis_title = "Counts/keV"
        if "mass_of_detector" in self.input_variables.keys()\
           and "total_time" in self.input_variables.keys():
            kilos = self.input_variables["mass_of_detector"]
            time_in_years = self.input_variables["total_time"]
            scaling = 1./(kilos*time_in_years*365.25)
            axis_title = "Counts/keV/kg/d"

        while i < number_iterations:
            #ROOT.RooTrace.dump(ROOT.cout, True)
            #ROOT.RooTrace.mark()
            self.logging("Process %s: Iteration (%i) of (%i)" 
                    % (os.getpid(), i+1, number_iterations))
            model_amplitude.setVal(0)

            # Generate the data
            # Reset the variables to the initial values in the cache
            data_model.getVariables().readFromStream(ROOT.istringstream(var_cache.str()), False)
            data_set_func = data_model.generate(variables)
            data_cache = None
            if do_bin_data:
                data_cache = data_set_func
                data_set_func = data_cache.binnedClone()
            
    
            if not data_set_func:
                print "Background entries are much too low, need to estimate with FC or Rolke."
                break
            # Perform the fit and find the limits
            get_val = self.find_confidence_value_for_model(
                model, 
                data_set_func, 
                model_amplitude, 
                confidence_value/2,
                mult_factor,
                print_level,
                verbose,
                self.show_plots) 
    
            if not get_val or self.is_exit_requested(): 
                # There was an error somewhere downstream
                # or an interrupt was signalled
                # Get out
                break
            elif get_val == self.retry_error: 
                # Calling function requested a retry
                continue
    
            # Store the results
            list_of_values.append(get_val)
            i += 1
    
            if self.show_plots:
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
                         ROOT.RooFit.Components("*Gauss_Signal*"), 
                         ROOT.RooFit.LineStyle(ROOT.RooFit.kDashed))
                
                    model.plotOn(aframe, 
                         ROOT.RooFit.Components("flat_*"), 
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
                    axis.SetTitleOffset(1.09)
                    aframe.Draw()
                    axis.Draw()
                    self.c1.Update()
                    title = aframe.GetTitle()
                    title = title.replace(' ','').replace('(','').replace(')','').replace(',','') 
                    self.c1.Print(title + ("%s.eps" % var_obj.GetName()))
                    if not self.print_out_plots:
                        raw_input("Hit Enter to continue")
 
                   
 
            # ROOT doesn't play nicely with python always, 
            # so we have to delete by hand
            data_set_func.IsA().Destructor(data_set_func)
    
        return list_of_values

