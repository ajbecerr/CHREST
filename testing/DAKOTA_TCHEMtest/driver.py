import dakota.interfacing as di
import numpy as np
import os
import subprocess


#-----------------------------------------------------------------
params, results = di.read_parameters_file()
x1=params['x1']*1E+00
x2=params['x2']*1E+10
x3=params['x3']*1E+07
x4=params['x4']*1E+05
x5=params['x5']*1E+00

#-----------------------------------------------------------------

with open ('ORIGINAL_pentane.mech.dat', "r") as myfile:
    inputfile = myfile.readlines()
    #print(inputfile[19])
    inputfile[21] = 'OH(6)+OH(6)(+M)=H2O2(10)(+M) '+str(x1)+' 0.000     0.000\n'
    inputfile[31] = 'OH(6)+pentane(1)=H2O(8)+C5H11(426) '+str(x2)+' 0.935     0.505\n'
    inputfile[32] = 'OH(6)+pentane(1)=H2O(8)+C5H11(425) '+str(x3)+' 1.813     0.868\n'
    inputfile[33] = 'O2(2)+CH4(17)=HO2(4)+CH3(16) '+str(x4)+' 2.745     51.714\n'
    inputfile[34] = 'OH(6)+CH2O(14)=H2O(8)+HCO(15) '+str(x5)+' 0.000     0.000\n'
    #for i in range(len(inputfile)):
    #        inputfile[i] = inputfile[i][0:-1]
    np.savetxt('pentane.mech.dat', inputfile, fmt='%s',delimiter='')

#-----------------------------------------------------------------
command = 'module use /projects/academic/chrest/modules; module load chrest/release; export TCHEM_INSTALL_PATH=/projects/academic/chrest/lib/tchem/v2.0.0_25-04-2022_6ae59e8; exec=$TCHEM_INSTALL_PATH/example/TChem_IgnitionZeroD.x; this="$exec --chemfile=inputs/pentane.yaml --use-cvode=false --samplefile=inputs/pentane.dat --outputfile=outputs/IgnSolution.dat --atol-newton=1e-18 --rtol-newton=1e-8 --run-constant-pressure=false --max-newton-iterations=20 --tol-time=1e-6 --useYaml=true --dtmax=1e-3 --dtmin=1e-20 --tend=0.1 --time-iterations-per-interval=10 --jacobian-interval=5 --max-time-iterations=5000 --ignition-delay-time-file=outputs/IgnitionDelayTime.dat --ignition-delay-time-w-threshold-temperature-file=outputs/IgnitionDelayTimeTthreshold.dat --threshold-temperature=1500"; echo $this; eval $this'

p = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

(stdout, err) = p.communicate()

output = stdout.splitlines()


#-----------------------------------------------------------------
with open ('outputs/IgnitionDelayTimeTthreshold.dat', "r") as myoutfile:
    outfile = myoutfile.readlines()
    QoI = float(outfile[1])
    print(QoI)

os.remove('pentane.mech.dat')
os.remove('outputs/IgnitionDelayTimeTthreshold.dat')


#-----------------------------------------------------------------
for i, r in enumerate(results.responses()):
    r.function = QoI
results.write()
