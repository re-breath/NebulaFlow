import numpy as np
import matplotlib.pyplot as plt
from pylab import *

energy_train = loadtxt('energy_train.out')
force_train = loadtxt('force_train.out')
virial_train = loadtxt('virial_train.out')

rmse_energy = np.sqrt(np.mean((energy_train[:,0]-energy_train[:,1])**2))
force_diff = np.reshape(force_train[:,3:6]-force_train[:,0:3], (force_train.shape[0]*3, 1))
rmse_force = np.sqrt(np.mean(force_diff**2))
virial_train = virial_train[virial_train[:, 1] > -1000, :]
rmse_virial = np.sqrt(np.mean((virial_train[:, 0:5] - virial_train[:, 6:11])**2))

plt.figure(figsize=(14,10))

plt.subplot(2,2,1)
#loglog(loss[:, 1:7])
#loglog(loss[:, 7:10])
xlabel('Generation/100')
ylabel('Loss')
legend(['Total', 'L1-regularization', 'L2-regularization', 'Energy-train', 'Force-train', 'Virial-train', 'Energy-test', 'Force-test', 'Virial-test'])
tight_layout()

plt.subplot(2,2,2)
plot(energy_train[:, 1], energy_train[:, 0], '.')
plot(linspace(-9,-8), linspace(-9,-8), '-')
xlabel('DFT energy (eV/atom)')
ylabel('NEP energy (eV/atom)')
legend(['train'])
plt.title(f'RMSE = {rmse_energy:.5f} eV/atom')
tight_layout()

plt.subplot(2,2,3)
plot(force_train[:, 3:6], force_train[:, 0:3], '.')
plot(linspace(-10,10), linspace(-10,10), '-')
xlabel('DFT force (eV/A)')
ylabel('NEP force (eV/A)')
legend(['train x direction', 'train y direction', 'train z direction'])
plt.title(f'RMSE = {rmse_force:.5f} eV/A')
tight_layout()

plt.subplot(2,2,4)
plot(virial_train[:, 6:11], virial_train[:, 0:5], '.')
plot(linspace(-2,2), linspace(-2,2), '-')
xlabel('DFT virial (eV/atom)')
ylabel('NEP virial (eV/atom)')
legend(['train'])
plt.title(f'RMSE = {rmse_virial:.5f} meV/atom')
tight_layout()

#plt.savefig('nep.png', dpi=150, bbox_inches='tight')
plt.savefig('nep.png')
