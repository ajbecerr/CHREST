import dakota.interfacing as di
import numpy as np
import os
import subprocess
# random file tag (without having to use mpi4py)
tag = str(np.random.randint(0,100000))


#-----------------------------------------------------------------
params, results = di.read_parameters_file()
x1=params['x1']*1E+15
x2=params['x2']*1E+9

#-----------------------------------------------------------------

with open ('ORIGINAL_2S_CH4_CM2.mech.dat', "r") as myfile:
    inputfile = myfile.readlines()
    inputfile[19] = 'CH4+1.5O2=>CO+2H2O '+str(x1)+'  0.00   35000.00\n'
    inputfile[22] = 'CO+0.5O2<=>CO2 '+str(x2)+'  0.000   12000.00\n'
    np.savetxt('2S_CH4_CM2_'+tag+'.mech.dat', inputfile, fmt='%s', delimiter='')

#-----------------------------------------------------------------

with open ('input_template.yaml', "r") as myfile:
    inputfile = myfile.readlines()
    inputfile[2] = '  title: _ignitionDelay2S_CH4_CM2_'+tag+'\n'
    inputfile[24] = '          mechFile: 2S_CH4_CM2_'+tag+'.mech.dat\n'
    inputfile[60] = '        name: ignitionDelayTemperature_'+tag+'.txt\n'
    np.savetxt('2S_CH4_CM2_'+tag+'.yaml', inputfile, fmt='%s', newline='')

#-----------------------------------------------------------------
command = 'module use /projects/academic/chrest/modules; module load chrest/release; $ABLATE_DIR/ablate --input 2S_CH4_CM2_'+tag+'.yaml  '

p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

(stdout, err) = p.communicate()

output = stdout.splitlines()


#-----------------------------------------------------------------
targetline = output[0].split(':')
QoI = float(targetline[1])
print(QoI)

os.remove('2S_CH4_CM2_'+tag+'.yaml')
os.remove('2S_CH4_CM2_'+tag+'.mech.dat')
os.remove('_ignitionDelay2S_CH4_CM2_'+tag+'/ignitionDelayTemperature_'+tag+'.txt')


#-----------------------------------------------------------------
for i, r in enumerate(results.responses()):
    r.function = QoI
results.write()
