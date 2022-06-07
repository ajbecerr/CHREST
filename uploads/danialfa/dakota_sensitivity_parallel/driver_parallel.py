import dakota.interfacing as di
import numpy as np
import os
import subprocess



# random file tag (without having to use mpi4py)
tag = str(np.random.randint(0,100000))


#-----------------------------------------------------------------
params, results = di.read_parameters_file()
x1=params['x1']*100000.
x2=params['x2']*10.
#-----------------------------------------------------------------

with open ('input_template.yaml', "r") as myfile:
    inputfile = myfile.readlines()
    # print(inputfiles[66])
    inputfile[66] = '        formula: "1.0, ' +str(x1)+  ', ' +str(x2)+ ', 0.0"\n'
    for i in range(len(inputfile)):
            inputfile[i] = inputfile[i][0:-1]
    np.savetxt('temp'+tag+'.yaml', inputfile, fmt='%s',delimiter='')

#-----------------------------------------------------------------
# os.system('/projects/academic/danialfa/software/chrest_ablet/debug/ablate --input temp.yaml  > "temp.out"  ')

command = '/projects/academic/danialfa/software/chrest_ablet/debug/ablate --input temp'+tag+'.yaml  '

p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

(stdout, err) = p.communicate()

output = stdout.splitlines()


#-----------------------------------------------------------------
for i in range(len(output)-1,-1, -1):
    if 'euler' in output[i]:
        break

targetline = output[i].split(',')

QoI = float(targetline[1])


os.remove('temp'+tag+'.yaml')

#-----------------------------------------------------------------
for i, r in enumerate(results.responses()):
    r.function = QoI
results.write()
