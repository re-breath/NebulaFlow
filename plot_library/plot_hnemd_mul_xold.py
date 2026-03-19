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

def plot_hnemd_x(filename:str,color:str='C7',linewidth:int=1,alpha:float=0.5):
    kappa = load_kappa(filename)
    kappa.keys()
    t = np.arange(1,kappa['kxi'].shape[0]+1)*0.001  # ns
    kappa['kyi_ra'] = running_ave(kappa['kyi'],t)
    kappa['kyo_ra'] = running_ave(kappa['kyo'],t)
    kappa['kxi_ra'] = running_ave(kappa['kxi'],t)
    kappa['kxo_ra'] = running_ave(kappa['kxo'],t)
    kappa['kz_ra'] = running_ave(kappa['kz'],t)

    #figure(figsize=(12,10))

    set_fig_properties([gca()])
    plot(t, kappa['kxi_ra']+kappa['kxo_ra'], color, linewidth,alpha)

    print(f"最后一个x方向的热导率为 ：{kappa['kxi_ra'][-1]+kappa['kxo_ra'][-1]}")

for i in range(1, 6):
    plot_hnemd_x(f"kappa_{i}.out")

#画出总的kappa
plot_hnemd_x("kappa.out",color='C3',linewidth=2,alpha=1)
# xlim([0, 10])
# gca().set_xticks(range(0,11,2))
# ylim([-500,1000])
# gca().set_yticks(range(-500,1001,300))
xlabel('time (ns)')
ylabel(r'$\kappa$ (W/m/K)')
#legend(['yx', 'xx', 'zx'])
title('kappa x direction')

tight_layout()
savefig('hnemd_mul_x.png')


