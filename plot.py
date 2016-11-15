import numpy as np
import matplotlib.pyplot as plt

nside=128
T=np.loadtxt("T_out_128.txt")

plt.imshow(T)
plt.title("T for nside="+str(nside))
plt.show()

print "Mean="+str(np.mean(T))
