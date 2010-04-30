#!/usr/bin/env python
import os
import time

# First the vitals of the job
queue = ["default", "scavenge"]
sens_home = "/share/home/mgmarino/sensitivity/SensitivityCalculation"
path_of_results = sens_home
cores = 8
nodes = 4
hours, minutes, seconds = ("6", "00", "00")

wimp_list = range(4, 20+1, 2)
wimp_list.extend(range(30, 100+1, 10))


#Now the particulars of the calculations 

file_basename = "mgm_wimp_many_two" 
max_time = 3600*int(hours) + 60*int(minutes)
max_time -= 60 # allow for some extra padding in the wall clock
j = 0
#wimp_list = [100, 80, 60, 40]
for ajob in wimp_list:
    job = str(ajob)
    execution_string=""
    execution_string+="""
trap : SIGINT
trap "echo TERM" SIGTERM
cd %s
time mpirun -np %i -hostfile $PBS_NODEFILE python %s/pyWIMP/mpi_job_engine.py %g 
     """ % (sens_home, cores*nodes, sens_home, ajob)
    
    print "Submitting job"
    write_handle = os.popen("""
qsub \\
  -N %s_job_%s \\
  -q %s \\
  -l nodes=%i:ppn=%i\\
  -l walltime=%s:%s:%s\\
  -e %s/output_%s_job_%s.err\\
  -o %s/output_%s_job_%s.out\\
        """ % (file_basename, job, queue[j%2], nodes, cores,
               hours, minutes, seconds,
               path_of_results, file_basename, job,
               path_of_results, file_basename, job), 
            'w')
    write_handle.write(execution_string)
    write_handle.close()
    j += 1
