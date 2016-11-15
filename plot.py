import numpy as np
import matplotlib.pyplot as plt

T=np.loadtxt("T_out_128.txt")

plt.imshow(T)
plt.show()
