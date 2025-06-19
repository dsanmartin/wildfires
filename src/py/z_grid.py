import sys
import numpy as np
import matplotlib.pyplot as plt
# Numpy formatting
np.set_printoptions(precision=2, suppress=True, linewidth=100)

z_min, z_max = 0, float(sys.argv[1])
Nz = int(sys.argv[2])
k = float(sys.argv[3])
z = np.linspace(z_min, z_max, Nz)
r = lambda z: (z_max - z_min) * (np.exp(k * (z - z_min) / (z_max - z_min)) - 1) / (np.exp(k) - 1) + z_min
z_r = r(z)
z_r_fil = z_r[(z_r <= 0.51) == True]
print("z:")
print(z_r)
print()
print("z <= 0.51:")
print(z_r_fil, z_r_fil.shape[0])
print()
dz = z_r[1:] - z_r[:-1]
dz_kp1 = dz[1:]
dz_k = dz[:-1]
quotient = dz_kp1 / dz_k
print("dz_kp1 / dz_k:")
print(quotient)                                      

plt.figure(figsize=(12, 6))
plt.plot(z_r, z_r, 'bo')
plt.grid(True)
plt.show()