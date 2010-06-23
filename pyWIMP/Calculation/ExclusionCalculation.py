import BaseCalculation
import ROOT
import numpy

class ExclusionCalculation(BaseCalculation.BaseCalculation):
    """
    This class is an engine for calculating an exclusion for
    a certain model given a data set.  All variables are RooFit
    type variables.
    """

    def find_confidence_value_for_model(self, 
                                        model, 
                                        data, 
                                        model_amplitude, 
                                        conf_level, 
                                        mult_factor, 
                                        print_level = -1, 
                                        verbose = False, 
                                        debug = False,
                                        tolerance = 0.001):
    
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


        output_list = []

        output_dict = {}
        
        output_list = numpy.zeros((number_of_points, 2)) 
        model_amplitude.setConstant(True)

        orig = 1e15 
        min_point = 0
        j = 0
        step_size = (max_range - min_value)/(number_of_points-1)
        test_points = numpy.arange(min_value, max_range + step_size*0.5, step_size)
        for test_val in test_points: 
            #print "Performing: ", test_val 
            model_amplitude.setVal(test_val)
            minuit.migrad()
            res = minuit.save(str(test_val)) 
            min_val = res.minNll()
            res.IsA().Destructor(res)
                
            if debug:
                self.print_plot(model, data, str(test_val))
            if min_val < orig:  
                orig = min_val 
                min_point = j
            output_list[j] = [test_val, min_val]
            j += 1
            # This is the most dense loop, so checking if we should get out
            if self.is_exit_requested(): return None
            
       
        # Now find the confidence_level using unbounded and bounded PLL
        
        # Grab the best fit value (at the min_point)
        best_fit = output_list[min_point][0]
        bounded_min_nll = min_nll
        unbounded_list = output_list.copy()
        unbounded_list -= [0, min_nll]
        bounded_list = output_list.copy() 
        if best_fit < 0:
            # Means the best fit was less than 0, in which case, we need to 
            # bound the lower limit, taking the point where the model_amplitude
            # is equal to 0 
            bounded_list = output_list[numpy.where(bounded_list[:,0] >= 0)]
            model_amplitude.setVal(0)
            minuit.migrad()
            res = minuit.save()
            bounded_min_nll = res.minNll()
            res.IsA().Destructor(res)
        # Now this list is shifted correctly
        bounded_list -= [0, bounded_min_nll]
        

        # Generate *helper* curves, these allow us to extrapolate
        # between points assuming that the likelihood is reasonably 
        # smooth
        unbounded_curve = ROOT.RooCurve()
        [unbounded_curve.addPoint(x,y) for x, y in unbounded_list]
        bounded_curve = ROOT.RooCurve()
        [bounded_curve.addPoint(x,y) for x, y in bounded_list]

        # Now find where each rises to the particular confidence level value
        step_size = 0.01
        
        # finding unbounded, upper limit
        unbounded_start = best_fit 
        while (unbounded_curve.Eval(unbounded_start) < conf_level and 
               unbounded_start <= model_amplitude.getMax()): unbounded_start += step_size
        unbounded_upper_limit = unbounded_start 

        # finding unbounded, lower limit
        unbounded_start =  best_fit
        while (unbounded_curve.Eval(unbounded_start) < conf_level and 
               unbounded_start >= 0): unbounded_start -= step_size
        unbounded_lower_limit = unbounded_start 

       
        # We only calculate the bounded upper limit if the best fit is below 0,
        # otherwise it's exactly the same as the unbounded limit
        bounded_limit = unbounded_upper_limit
        if best_fit < 0:
            bounded_start = bounded_curve.Eval(0) 
            while (bounded_curve.Eval(bounded_start) < conf_level and 
                   bounded_start <= model_amplitude.getMax()): bounded_start += step_size
            bounded_limit = bounded_start 
        
        # Save these bounds in the output dictionary
        output_dict['unbounded_lower_limit'] = unbounded_lower_limit*mult_factor
        output_dict['unbounded_upper_limit'] = unbounded_upper_limit*mult_factor
        output_dict['bounded_limit'] = bounded_limit*mult_factor
        

        # Return the results
        for obj in [nll, minuit]:
            obj.IsA().Destructor(obj)
        return output_dict
 
