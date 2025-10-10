import numpy as np
import matplotlib.pyplot as plt 
import h5py 

data = h5py.File("./pincflow_output.h5")

tke = data["tke"][:, :, 0, :]

print(tke[1, :, :])


cepsilon = 0.871
lturb = 100.
dt = 1000. 

turb = 4 / (cepsilon / lturb * dt + 2 / np.sqrt(0.1)) ** 2.

print(turb)