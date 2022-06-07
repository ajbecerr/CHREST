import numpy as np

Tin = 300 # K?
Pin = 10 # atm

with open ('../DAKOTA_TCHEMtest/ORIGINAL_FULL_pentane.yaml', "r") as myfile:
    inputfile = myfile.readlines()
    cheb_data = inputfile[1854:1863]

Tmin = float(cheb_data[0].split(': ')[1][:-1].split(',')[0][1:])
Tmax = float(cheb_data[0].split(': ')[1][:-1].split(',')[1][1:-1])
Pmin = float(cheb_data[1].split(': ')[1][:-1].split(',')[0][1:-4])
Pmax = float(cheb_data[1].split(': ')[1][:-1].split(',')[1][1:-5])

cheb_coeffs = [cheb_data[-6:][i][5:-2].split(', ') for i in range(6)]

def Tmap(Tin, Tmin, Tmax):
    Tout = (2*Tin**(-1) - Tmin**(-1) - Tmax**(-1))/(Tmax**(-1) - Tmin**(-1))
    return(Tout)

def Pmap(Pin, Pmin, Pmax):
    Pout = (2*np.log(Pin) - np.log(Pmin) - np.log(Pmax))/(np.log(Pmax) - np.log(Pmin))
    return(Pout)
    
Tout = Tmap(Tin, Tmin, Tmax)
Pout = Pmap(Pin, Pmin, Pmax)

def log_k(cheb_coeffs, Tout, Pout):
    S = 0
    for i in range(6):
        TCHEB = np.polynomial.chebyshev.Chebyshev([1 if i == k else 0 for k in range(i)])
        for j in range(4):
            PCHEB = np.polynomial.chebyshev.Chebyshev([1 if j == l else 0 for l in range(j)])
            S = S + float(cheb_coeffs[i][j])*TCHEB(Tout)*PCHEB(Pout)
    return(S)
    
print(log_k(cheb_coeffs, Tout, Pout))
print(np.exp(log_k(cheb_coeffs, Tout, Pout)))