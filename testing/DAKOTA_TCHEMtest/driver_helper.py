import numpy as np
import sys

x2 = str(sys.argv[1])
x3 = str(sys.argv[2])
x4 = str(sys.argv[3])
x5 = str(sys.argv[4])

with open ('ORIGINAL_FULL_pentane.yaml', "r") as myfile:
    inputfile = myfile.readlines()
    # inputfile[1853] = ''
    inputfile[1900] = '  rate-constant: {A: 1.411e+10, b: 0.935, Ea: '+x2+'}\n'
    inputfile[1923] = '  rate-constant: {A: 2.732e+07, b: 1.813, Ea: '+x3+'}\n'
    inputfile[1975] = '  rate-constant: {A: 2.03e+05, b: 2.745, Ea: '+x4+'}\n'
    inputfile[2629] = '  - {P: 10.0 atm, A: 1.07e+09, b: 1.33, Ea: -'+x5+'}\n'
    np.savetxt('FULL_pentane.yaml', inputfile, fmt='%s', newline='')