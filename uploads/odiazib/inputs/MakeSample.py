#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import numpy as np
import sys
import json

savedir = "inputs"
one_atm = 101325

# temperature [K]
Temp = np.array([1350, 1325, 1300, 1275, 1250, 1225, 1200, 1150, 1100, 1075, 1050, 1025, 1012.5, 1000, 987.5,
              975, 962.5, 950, 937.5, 925, 912.5, 900, 875, 850, 825, 800, 775, 750, 725, 700, 675, 650])

# pressure [atm]
Pressure = 10*one_atm 
Nsamples = len(Temp)
values = np.empty([Nsamples, 5])

values[:,1]=Pressure

for i in range(Nsamples):
    values[i,0] = Temp[i]

# use mass fraction for species
values[:,2] = 0.061626  #"pentane"
values[:,3] = 0.21864  # "O2"
values[:,4] = 1. -  values[:,2] - values[:,3]# "N2"

names = ["T", "P","PENTANE(1)", "O2(2)", "N2"]
# make input for example-interfaces
header = "T"
for iname in names[1:]:
    header +=" "+iname
np.savetxt('sample.dat', values, header=header ,comments='')
# make input for tchem-main interface
sample = {}
sample.update({"variable name": names})
sample.update({"state vector":values.tolist()})

with open('sample.json', 'w') as outfile:
    json.dump(sample, outfile)
