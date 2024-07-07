from pylab import *
from ase.build import graphene_nanoribbon
from gpyumd.atoms import GpumdAtoms
from ase.io import write
from gpyumd.load import load_shc, load_kappa
from gpyumd.math import running_ave
from gpyumd.calc import calc_spectral_kappa
aw = 2
fs = 16
font = {'size'   : fs}
matplotlib.rc('font', **font)
matplotlib.rc('axes' , linewidth=aw)

def set_fig_properties(ax_list):
    tl = 8
    tw = 2
    tlm = 4

    for ax in ax_list:
        ax.tick_params(which='major', length=tl, width=tw)
        ax.tick_params(which='minor', length=tlm, width=tw)
        ax.tick_params(which='both', axis='both', direction='in', right=True, top=True)
kappa = load_kappa()
kappa.keys()
t = np.arange(1,kappa['kxi'].shape[0]+1)*0.001 # ns
kappa['kyi_ra'] = running_ave(kappa['kyi'],t)
kappa['kyo_ra'] = running_ave(kappa['kyo'],t)
kappa['kxi_ra'] = running_ave(kappa['kxi'],t)
kappa['kxo_ra'] = running_ave(kappa['kxo'],t)
kappa['kz_ra'] = running_ave(kappa['kz'],t)
figure(figsize=(12,10))
subplot(1,1,1)
set_fig_properties([gca()])
#plot(t, kappa['kz'],color='C7',alpha=0.5)
plot(t, kappa['kz_ra'], linewidth=2)
print("z方向热导率的最后一个值为：",kappa['kz_ra'][-1])
# xlim([0, 10])
# gca().set_xticks(range(0,11,2))
# ylim([-2000, 4000])
# gca().set_yticks(range(-2000,4001,1000))
xlabel('time (ns)')
ylabel(r'$\kappa$ W/m/K')
title('(a)')



tight_layout()
savefig('hnemd.png')
