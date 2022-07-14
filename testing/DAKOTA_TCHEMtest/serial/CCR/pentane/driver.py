import dakota.interfacing as di
import numpy as np
import os
import subprocess


#-----------------------------------------------------------------
params, results = di.read_parameters_file()
# x1=params['x1']*1E+00
x2=params['x2']*1E-01
x3=params['x3']*1E-01
x4=params['x4']*1E+01
x5=params['x5']*1E-01

#-----------------------------------------------------------------

with open ('ORIGINAL_FULL_pentane.yaml', "r") as myfile:
    inputfile = myfile.readlines()
    # inputfile[1853] = ''
    inputfile[1900] = '  rate-constant: {A: 1.411e+10, b: 0.935, Ea: '+str(x2)+'}\n'
    inputfile[1923] = '  rate-constant: {A: 2.732e+07, b: 1.813, Ea: '+str(x3)+'}\n'
    inputfile[1975] = '  rate-constant: {A: 2.03e+05, b: 2.745, Ea: '+str(x4)+'}\n'
    inputfile[2629] = '  - {P: 10.0 atm, A: 1.07e+09, b: 1.33, Ea: -'+str(x5)+'}\n'
    np.savetxt('FULL_pentane.yaml', inputfile, fmt='%s', newline='')

#-----------------------------------------------------------------
command = 'bash one_sample.sh'
p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

(stdout, err) = p.communicate()

output = stdout.splitlines()


#-----------------------------------------------------------------
with open ('outputs/IgnitionDelayTimeTthreshold.dat', "r") as myoutfile:
    outfile = myoutfile.readlines()
    QoI = float(outfile[1])
    print(QoI)

os.remove('FULL_pentane.yaml')
os.remove('outputs/IgnitionDelayTimeTthreshold.dat')


#-----------------------------------------------------------------
for i, r in enumerate(results.responses()):
    r.function = QoI
results.write()