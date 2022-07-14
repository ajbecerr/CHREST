import numpy as np
import os
import subprocess

tag = str(np.random.randint(0,100000))

#-----------------------------------------------------------------
x2 = 0.505
x3 = 0.868
x4 = 51.714
x5 = 0.392

#-----------------------------------------------------------------
with open ('ORIGINAL_FULL_pentane.yaml', "r") as myfile:
    inputfile = myfile.readlines()
    inputfile[1900] = '  rate-constant: {A: 1.411e+10, b: 0.935, Ea: '+str(x2)+'}\n'
    inputfile[1923] = '  rate-constant: {A: 2.732e+07, b: 1.813, Ea: '+str(x3)+'}\n'
    inputfile[1975] = '  rate-constant: {A: 2.03e+05, b: 2.745, Ea: '+str(x4)+'}\n'
    inputfile[2629] = '  - {P: 10.0 atm, A: 1.07e+09, b: 1.33, Ea: -'+str(x5)+'}\n'
    np.savetxt('FULL_pentane_'+tag+'.yaml', inputfile, fmt='%s', newline='')

#-----------------------------------------------------------------
command = 'bash one_sample_parallel.sh '+tag
p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
(stdout, err) = p.communicate()
output = stdout.splitlines()

#-----------------------------------------------------------------
with open ('outputs/IgnitionDelayTimeTthreshold_'+tag+'.dat', "r") as myoutfile:
    outfile = myoutfile.readlines()
    QoI = float(outfile[1])
    print(QoI)

#-----------------------------------------------------------------