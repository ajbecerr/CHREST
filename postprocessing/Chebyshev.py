import matplotlib.pyplot as plt
import numpy as np

Tlist = np.load('../preprocessing/Tlist.npy')
Plist = np.load('../preprocessing/Plist.npy')
klist = np.load('../preprocessing/klist.npy')

Tlist = np.kron(np.array([[1] for i in range(5)]), Tlist).flatten()
Plist = np.transpose(np.kron(np.array([[1] for i in range(50)]), Plist)).flatten()
klist = np.load('../preprocessing/klist.npy').flatten()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(np.array([1000/temp for temp in Tlist]), np.log10(Plist), np.log10(klist))

def Tmap(Tin, Tmin, Tmax):
    Tout = (2*(Tin**(-1)) - Tmin**(-1) - Tmax**(-1))/(Tmax**(-1) - Tmin**(-1))
    return(Tout)

def Pmap(Pin, Pmin, Pmax):
    Pout = (2*np.log(Pin) - np.log(Pmin) - np.log(Pmax))/(np.log(Pmax) - np.log(Pmin))
    return(Pout)
    
def log10_k(Tin, Pin):
    with open ('../testing/DAKOTA_TCHEMtest/ORIGINAL_FULL_pentane.yaml', "r") as myfile:
        inputfile = myfile.readlines()
        cheb_data = inputfile[1854:1863]
    Tmin = float(cheb_data[0].split(': ')[1][:-1].split(',')[0][1:])
    Tmax = float(cheb_data[0].split(': ')[1][:-1].split(',')[1][1:-1])
    Pmin = float(cheb_data[1].split(': ')[1][:-1].split(',')[0][1:-4])
    Pmax = float(cheb_data[1].split(': ')[1][:-1].split(',')[1][1:-5])
    cheb_coeffs = [cheb_data[-6:][i][5:-2].split(', ') for i in range(6)]
    Tout = Tmap(Tin, Tmin, Tmax)
    Pout = Pmap(Pin, Pmin*101325, Pmax*101325)
    S = 0
    for i in range(1, 7):
        TCHEB_coeffs = [0 for k in range(i+1)]
        TCHEB_coeffs[-1] = 1
        TCHEB = np.polynomial.chebyshev.Chebyshev(TCHEB_coeffs)
        for j in range(1, 5):
            PCHEB_coeffs = [0 for l in range(j+1)]
            PCHEB_coeffs[-1] = 1
            PCHEB = np.polynomial.chebyshev.Chebyshev(PCHEB_coeffs)
            S = S + float(cheb_coeffs[i-1][j-1])*TCHEB(Tout)*PCHEB(Pout)
    return(S)

kpred = np.array([log10_k(Tlist[i], Plist[i]) for i in range(250)])

ax.scatter([1000/temp for temp in Tlist], np.log10(Plist), kpred)
fig.savefig('Chebyshev.PNG')


# for i in range(7):
#     coeffs_i = [0 for j in range(i+1)]
#     coeffs_i[-1] = 1
#     plt.scatter(list(np.linspace(-1, 1)), [np.polynomial.chebyshev.Chebyshev(coeffs_i)(x) for x in np.linspace(-1, 1)])
# plt.savefig('test.PNG')