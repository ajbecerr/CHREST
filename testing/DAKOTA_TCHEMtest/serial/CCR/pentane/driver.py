import dakota.interfacing as di
import os
import subprocess


#-----------------------------------------------------------------
params, results = di.read_parameters_file()
# x1=params['x1']*1E+00
x2=params['x2']*1E-01
x3=params['x3']*1E+01
x4=params['x4']*1E+01
x5=params['x5']*1E-01

#-----------------------------------------------------------------
command = 'bash driver_helper.sh '+str(x2)+' '+str(x3)+' '+str(x4)+' '+str(x5)
p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

(stdout, err) = p.communicate()

output = stdout.splitlines()


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
