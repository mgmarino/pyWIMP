import DataCalculation
import ROOT
import os
import math
from exceptions import Exception
from ..utilities.utilities import rescale_frame
import numpy

class DataCalcVerification(DataCalculation.DataCalculation):

    """
    This class handles calculating the exclusion limits given a 
    certain input data set.  We have to overload 
    scan_confidence_value_space_for_model to do this.
    """


    def find_confidence_value_for_model(self, 
                                        model, 
                                        data, 
                                        model_amplitude, 
                                        conf_level, 
                                        print_level = -1):
    

        """
        This function scans the profile likelihood space of a parameter.
        It returns the profile likelihood in an array like:

        [(x1, pll(x1)), ... (xN, pll(xN))]

        where xN is the value of the model_amplitude and pll(xN) is the value of
        the profile likelihood at value xN.  The range (currently)
        is automatically calculated, but will not go below -100 for the
        model_amplitude.  This will eventually be an input.  

        The input is as follows:

        model: model to fit (B + S), this is of type RooAbsPdf
        data: data to fit, this is of type RooAbsData
        model_amplitude: test parameter of the signal.  A profile likelihood will be
          generated for this variable.

        conf_level: the delta Likelihood value that this function must go by.
        print_level: printing level of RooMinuit

        
        """
        #absolute_min = -20
        number_of_points = 100
        pars = model.getParameters(data)
        expo_coef = pars.find("expo_const_one")

        #model_amplitude.setError(0.5)
        #model_amplitude.setVal(0)
        #model_amplitude.setMin(-10)
        model_amplitude.setConstant(True)

        nll = model.createNLL(data, ROOT.RooFit.Verbose((print_level > 0)))
       
        minuit = ROOT.RooMinuit(nll)
        minuit.setPrintLevel(print_level)
        minuit.migrad()

        # Now fit with the model_amplitude
        model_amplitude.setVal(0)
        model_amplitude.setConstant(False)
        minuit.migrad()

        distance_from_min = 20.

        # The first point is the minimum likelihood, *but* it's possible
        # we could find other, more negative, points
        # It is also possible to have a linear likelihood which 
        # has no minimum.  In this case it's important to choose a bound for 
        # the function which is generally at a model_amplitude value of 0 
        # If we take the value of the nll at this point, it will be close enough.
        #while math.fabs(model_amplitude.getMin() - model_amplitude.getVal()) < distance_from_min: 
        #   model_amplitude.setMin(model_amplitude.getMin() - distance_from_min) 
        #   if model_amplitude.getMin() < absolute_min: break
        #   minuit.migrad()
        min_nll = nll.getVal()

        best_fit_value = model_amplitude.getVal()
        # Determine the range to scan over
        min_value = model_amplitude.getVal() - distance_from_min 
        if min_value > 0: min_value = 0
        if min_value < model_amplitude.getMin(): min_value = model_amplitude.getMin() 
        max_range = model_amplitude.getVal() + 50
        if max_range < 10: max_range = 50



        # Scan to make sure that the likelihood
        model_amplitude.setConstant(True)
        model_amplitude.setVal(max_range)
        minuit.migrad()
        #
        ## Multiply by two here to make sure we get a large range
        while nll.getVal()-min_nll < 2*conf_level: 
            max_range += 100
            model_amplitude.setVal(max_range)
            if model_amplitude.getVal() == model_amplitude.getMax():
                self.logging("Resetting maximum:", model_amplitude.getMax() )
                model_amplitude.setMax(model_amplitude.getVal()*2)
            minuit.migrad()

        # Reset it back to the original position
        model_amplitude.setVal(0)
        model_amplitude.setConstant(False)
        minuit.migrad()


        
        model_amplitude.setConstant(True)
        
        output_list = numpy.zeros((number_of_points, 2)) 
        j = 0
        step_size = (max_range - min_value)/number_of_points
        test_value = min_value
        min_nll = 1e15
        min_point = 0
        for j in range(number_of_points): 
            model_amplitude.setVal(test_value)
            if self.debug:  
                self.logging("Performing: ", model_amplitude.getVal())
            minuit_val = minuit.migrad()
            saved = minuit.save()
            
            if self.debug:  
                self.logging("Saved: ", saved.minNll())
                self.logging("Minuit: ", minuit_val)
                self.logging("Test value: ", test_value)
                self.logging("NLL: ", nll.getVal())
            min_val = saved.minNll()
            output_list[j] = [test_value, min_val]
            if min_val < min_nll:  
                min_nll = min_val 
                min_point = j
            test_value += step_size
            
        output_list -= [0, min_nll]

        # Now find the confidence_level using unbounded and bounded PLL
        
        best_fit = output_list[min_point][0]
        unbounded_list = output_list
        bounded_list = unbounded_list.copy() 
        if best_fit < 0:
            bounded_list = output_list[numpy.where(bounded_list[:,0] >= 0)]
        bounded_list -= [0, bounded_list[0][1]]
        

        unbounded_curve = ROOT.RooCurve()
        [unbounded_curve.addPoint(x,y) for x, y in unbounded_list]
        bounded_curve = ROOT.RooCurve()
        [bounded_curve.addPoint(x,y) for x, y in bounded_list]
        # Now find where each rises to the particular value
        step_size = 0.01
        
        unbounded_start = best_fit 
        while (unbounded_curve.Eval(unbounded_start) < conf_level and 
               unbounded_start <= model_amplitude.getMax()): unbounded_start += step_size
        unbounded_upper_limit = unbounded_start 

        unbounded_start =  best_fit
        while (unbounded_curve.Eval(unbounded_start) < conf_level and unbounded_start >= 0): unbounded_start -= step_size
        unbounded_lower_limit = unbounded_start 

       
        bounded_limit = unbounded_upper_limit
        if best_fit < 0:
            bounded_start = bounded_curve.Eval(0) 
            while (bounded_curve.Eval(bounded_start) < conf_level and 
                   bounded_start <= model_amplitude.getMax()): bounded_start += step_size
            bounded_limit = bounded_start 
        
        # Now the first value in each of these should be the calculated limit

        if self.debug:
            print best_fit
            model_amplitude.setVal(best_fit)
            minuit.migrad()
            self.print_plot(model, data)
            return (best_fit, unbounded_upper_limit, unbounded_lower_limit, bounded_limit, output_list)
        return (best_fit, unbounded_upper_limit, unbounded_lower_limit, bounded_limit)
 
 
    def scan_confidence_value_space_for_model(self, 
                                              model, 
                                              model_amplitude, 
                                              variables, 
                                              number_of_events,
                                              number_iterations, 
                                              cl):
    
        print_level = -1
        if self.debug: print_level = 1

        confidence_value = ROOT.TMath.ChisquareQuantile(cl, 1) 

        # Save the values of the parameters to reset at the end
        var_cache = ROOT.ostringstream() 
        model.getVariables().writeToStream(var_cache, False)
        
        # Generate the data, use Extended flag
        # because the number_of_events is just
        # an expected number.
        # Generate the full set here, because we want to 
        # make sure it's with the correct model distribution

        data_list = [ model.generate(variables,  
                        ROOT.RooFit.NumEvents(number_of_events),
                        ROOT.RooFit.Extended(),
                        ROOT.RooFit.Name("Data_" + str(i))) 
                      for i in range(number_iterations) ]

       

        # Perform the fit and find the limits
        list_of_values = []
        iter = 1
        for data_model in data_list:
            self.logging("Iteration: %i of %i" % (iter, number_iterations))
            iter += 1
            get_val = self.find_confidence_value_for_model(
                model, 
                data_model, 
                model_amplitude, 
                confidence_value/2,
                print_level)
            
            if get_val is None: 
                # There was an error somewhere downstream
                # or an interrupt was signalled
                # Get out
                break
            #elif get_val == self.retry_error: 
                # Calling function requested a retry
            #    continue
    
            # Store the results
            list_of_values.append(get_val)
    
        # Reset the variables
        model.getVariables().readFromStream(ROOT.istringstream(var_cache.str()), False)
        return list_of_values

