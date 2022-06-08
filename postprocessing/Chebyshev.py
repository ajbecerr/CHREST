import matplotlib.pyplot as plt
import numpy as np

Tlist = np.load('../preprocessing/Tlist.npy')
Plist = np.load('../preprocessing/Plist.npy')
klist = np.load('../preprocessing/klist.npy')

Tlist = np.array([1000/temp for temp in np.kron(np.array([[1] for i in range(5)]), Tlist).flatten()])
Plist = np.log10(np.transpose(np.kron(np.array([[1] for i in range(50)]), Plist)).flatten())
klist = np.log10(np.load('../preprocessing/klist.npy').flatten())

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(Tlist, Plist, klist)
fig.savefig('Chebyshev.PNG')