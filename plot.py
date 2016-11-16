import numpy as np
import matplotlib.pyplot as plt

nside=128
T=np.loadtxt("T_out.txt")

plt.imshow(T)
plt.title("T for nside="+str(nside))
plt.colorbar()
plt.show()
plt.close()

print "Mean="+str(np.mean(T))
