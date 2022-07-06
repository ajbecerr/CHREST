import matplotlib.pyplot as plt

with open ('../../testing/DAKOTA_TCHEMtest/serial/CCR/pentane/out_chem_sens-8828214.out', 'r') as myfile:
    outputfile = myfile.readlines()
    x2 = outputfile[10849].split(' x2')[0][-16:]
    x3 = outputfile[10850].split(' x3')[0][-16:]
    x4 = outputfile[10851].split(' x4')[0][-16:]
    x5 = outputfile[10852].split(' x5')[0][-16:]

sobol_indices = [float(x2), float(x3), float(x4), float(x5)]

plt.bar(['$E_{92}$', '$E_{99}$', '$E_{116}$', '$E_{258}$'], sobol_indices)
plt.yscale('log')
plt.savefig('pentane_sens_vis.PNG')