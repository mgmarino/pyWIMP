#!/usr/bin/env python
import os
import time
import glob
import numpy

# First the vitals of the job
queue = ["default", "scavenge"]
sens_home = "/share/home/mgmarino/sensitivity/SensitivityCalculation"
output_dir = "/share/home/mgmarino/DataVerify3"
path_of_results = output_dir 
cores = 8
nodes = 2
hours, minutes, seconds = ("7", "00", "00")

base_dir = "/share/home/mgmarino/FitBegeData/output_data_sets"
all_root_files = glob.glob("%s/output_for_data_verification/results_rt_99/*.root" % base_dir)
all_root_files = numpy.array(all_root_files)
all_root_files = all_root_files[0::2] # split the size in two
#Now the particulars of the calculations 

file_basename = "mgm_wimp_data_verify" 
max_time = 3600*int(hours) + 60*int(minutes)
max_time -= 60 # allow for some extra padding in the wall clock
#wimp_list = [100, 80, 60, 40]
job = 0

if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

for afile in all_root_files:
    execution_string=""
    execution_string+="""
trap : SIGINT
trap "echo TERM" SIGTERM
cd %s
time mpirun -np %i -hostfile $PBS_NODEFILE python %s/pyWIMP/mpi_job_engine_test_model.py %s 
     """ % (output_dir, cores*nodes, sens_home, afile)
    
    print "Submitting job"
    write_handle = os.popen("""
qsub \\
  -N %s_job_%s \\
  -q %s \\
  -l nodes=%i:ppn=%i\\
  -l walltime=%s:%s:%s\\
  -e %s/output_%s_job_%s.err\\
  -o %s/output_%s_job_%s.out\\
        """ % (file_basename, job, queue[job%2], nodes, cores,
               hours, minutes, seconds,
               path_of_results, file_basename, job,
               path_of_results, file_basename, job), 
            'w')
    write_handle.write(execution_string)
    write_handle.close()
    job += 1
