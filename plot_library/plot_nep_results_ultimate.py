#该脚本用于绘制nep训练的最终结果
import numpy as np
import matplotlib.pyplot as plt
from pylab import *

aw = 2
fs = 16
font = {'size'   : fs}
matplotlib.rc('font', **font)
matplotlib.rc('axes' , linewidth=aw)

loss = loadtxt('loss.out')
energy_train = loadtxt('energy_train.out')
energy_test = loadtxt('energy_test.out')
force_train = loadtxt('force_train.out')
force_test = loadtxt('force_test.out')
virial_train = loadtxt('virial_train.out')
virial_test = loadtxt('virial_test.out')

rmse_energy = np.sqrt(np.mean((energy_train[:,0]-energy_train[:,1])**2))
force_diff = np.reshape(force_train[:,3:6]-force_train[:,0:3], (force_train.shape[0]*3, 1))
rmse_force = np.sqrt(np.mean(force_diff**2))
virial_train = virial_train[virial_train[:, 1] > -1000, :]
rmse_virial = np.sqrt(np.mean((virial_train[:, 0:5] - virial_train[:, 6:11])**2))

#计算MAE
mae_energy = np.mean(np.abs(energy_train[:,0]-energy_train[:,1]))
mae_force1 = np.mean(np.abs(force_train[:,3:6]-force_train[:,0:3])) #实验无需重新塑形，mean会自动reshape
mae_force = np.mean(np.abs(force_diff))
mae_virial = np.mean(np.abs(virial_train[:, 0:5] - virial_train[:, 6:11]))
print(f"MAE :\nenergy :{mae_energy:.5f}eV/atom\nforce :{mae_force:.5f}eV/A\nvirial:  {mae_virial:.5f}eV/A")

plt.figure(figsize=(12,10))

plt.subplot(2,2,1)
mulcolor = ['#aedefc','#de95ba','#4ecca3','#bfcfff','#ff6464']
loglog(loss[:, 1:2],color = mulcolor[4])
loglog(loss[:, 3:7])
ylim(6e-4,1)
#loglog(loss[:, 7:9])
xlabel('Generation/100')
ylabel('Loss')
legend(['Total', 'L2','RMSE-Energy-train','RMSE-Force-train','RMSE-Virial-trian' ],fontsize = 10)
plt.title('(a)',y=0.05,x=0.9,size=26)
tight_layout()

plt.subplot(2,2,2)
plot(energy_test[:, 1], energy_test[:, 0], marker = 's',color = '#f67280', markerfacecolor='none')
plot(energy_train[:, 1], energy_train[:, 0], marker = 's',color = '#f67280', markerfacecolor='none')
plot(linspace(-7.5,-6.6), linspace(-7.5,-6.6), '-')
xlabel('DFT energy (eV/atom)')
ylabel('NEP energy (eV/atom)')
#legend(['test', 'train'],fontsize = 10)
#plt.title(f'RMSE = {rmse_energy:.5f} eV/atom',fontsize = 10,x=0.7,y=0.02)
plt.text(-7.5,-6.65,f'RMSE = {rmse_energy:.5f} eV/atom\nMAE = {mae_energy:.5f} eV/atom',fontsize = 10)
plt.title('(b)',y=0.05,x=0.9,size=26)
tight_layout()

plt.subplot(2,2,3)
plt.plot(force_test[:, 3:6], force_test[:, 0:3], 'o', mfc='none', color = '#cca8e9',mew=1, ms=6)
plot(force_train[:, 3:6], force_train[:, 0:3], 'o', mfc='none', mew=1,color = '#cca8e9', ms=6)
plot(linspace(-13,15), linspace(-13,15), '-')
xlabel('DFT force (eV/Å)')
ylabel('NEP force (eV/Å)')
#legend(['test x direction', 'test y direction', 'test z direction', 'train x direction', 'train y direction', 'train z direction'],fontsize = 10)
plt.text(-13,13.5,f'RMSE = {rmse_force:.5f} eV/Å\nMAE = {mae_force:.5f} eV/Å',fontsize = 10)
plt.title('(c)',y=0.05,x=0.9,size=26)
tight_layout()

plt.subplot(2,2,4)
plot(virial_test[:, 6:11], virial_test[:, 0:5],'<', mfc='none', color = '#61c0bf',mew=1, ms=10)
plot(virial_train[:, 6:11], virial_train[:, 0:5],'<', mfc='none', color = '#61c0bf',mew=1, ms=10)
plot(linspace(-3,4), linspace(-3,4), '-')
xlabel('DFT virial (eV/atom)')
ylabel('NEP virial (eV/atom)')
#legend(['test', 'train'],fontsize = 10)
plt.text(-2.9,3.5,f'RMSE = {rmse_virial:.5f} eV/atom\nMAE = {mae_virial:.5f} eV/atom',fontsize = 10)
plt.title('(d)',y=0.05,x=0.9,size=26)
tight_layout()

#plt.savefig('Fig2.pdf',dpi=80)
plt.savefig('nep.png')
