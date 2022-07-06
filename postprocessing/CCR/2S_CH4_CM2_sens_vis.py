import matplotlib.pyplot as plt

with open ('../../testing/DAKOTA_ABLATEtest/serial/CCR/2S_CH4_CM2/out_chem_sens-8784140.out', 'r') as myfile:
    outputfile = myfile.readlines()
    x1 = outputfile[6449].split(' x1')[0][-16:]
    x2 = outputfile[6450].split(' x2')[0][-16:]

sobol_indices = [float(x1), float(x2)]

plt.bar(['$E_{1}$', '$E_{2}$'], sobol_indices)
plt.yscale('log')
plt.savefig('2S_CH4_CM2_sens_vis.PNG')