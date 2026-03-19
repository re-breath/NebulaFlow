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

t = np.arange(1,kappa['kxi'].shape[0]+1)*0.001  # ns
kappa['kyi_ra'] = running_ave(kappa['kyi'],t)
kappa['kyo_ra'] = running_ave(kappa['kyo'],t)
kappa['kxi_ra'] = running_ave(kappa['kxi'],t)
kappa['kxo_ra'] = running_ave(kappa['kxo'],t)
kappa['kz_ra'] = running_ave(kappa['kz'],t)
figure(figsize=(12,10))
subplot(2,2,1)
set_fig_properties([gca()])
plot(t, kappa['kyi'],color='C7',alpha=0.5)
plot(t, kappa['kyi_ra'], linewidth=2)
#xlim([0, 10])
#gca().set_xticks(range(0,11,2))
#ylim([-2000, 4000])
gca().set_yticks(range(-2000,4001,1000))
xlabel('time (ns)')
ylabel(r'$\kappa_{in}$ W/m/K')
title('(a)')

subplot(2,2,2)
set_fig_properties([gca()])
plot(t, kappa['kyo'],color='C7',alpha=0.5)
plot(t, kappa['kyo_ra'], linewidth=2, color='C3')
#xlim([0, 10])
#gca().set_xticks(range(0,11,2))
#ylim([0, 4000])
#gca().set_yticks(range(0,4001,1000))
xlabel('time (ns)')
ylabel(r'$\kappa_{out}$ (W/m/K)')
title('(b)')

subplot(2,2,3)
set_fig_properties([gca()])
plot(t, kappa['kyi_ra'], linewidth=2)
plot(t, kappa['kyo_ra'], linewidth=2, color='C3')
plot(t, kappa['kyi_ra']+kappa['kyo_ra'], linewidth=2, color='k')
#xlim([0, 10])
#gca().set_xticks(range(0,11,2))
#ylim([0, 4000])
#gca().set_yticks(range(0,4001,1000))
xlabel('time (ns)')
ylabel(r'$\kappa$ (W/m/K)')
legend(['in', 'out', 'total'])
title('(c)')


subplot(2,2,4)
set_fig_properties([gca()])
plot(t, kappa['kyi_ra']+kappa['kyo_ra'],color='k', linewidth=2)
plot(t, kappa['kxi_ra']+kappa['kxo_ra'], color='C0', linewidth=2)
plot(t, kappa['kz_ra'], color='C3', linewidth=2)
print(f"最后一个y方向的热导率为：{kappa['kyi_ra'][-1]+kappa['kyo_ra'][-1]}")
#xlim([0, 10])
#gca().set_xticks(range(0,11,2))
#ylim([0, 4000])
#gca().set_yticks(range(-2000,4001,1000))
xlabel('time (ns)')
ylabel(r'$\kappa$ (W/m/K)')
legend(['yy', 'xy', 'zy'])
title('(d)')

tight_layout()
plt.savefig('hnemd.png')