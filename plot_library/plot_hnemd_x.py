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
plot(t, kappa['kxi'],color='C7',alpha=0.5)
plot(t, kappa['kxi_ra'], linewidth=2)
# xlim([0, 10])
# gca().set_xticks(range(0,11,2))
# ylim([-2000, 4000])
# gca().set_yticks(range(-2000,4001,1000))
xlabel('time (ns)')
ylabel(r'$\kappa_{in}$ W/m/K')
title('(a)')

subplot(2,2,2)
set_fig_properties([gca()])
plot(t, kappa['kxo'],color='C7',alpha=0.5)
plot(t, kappa['kxo_ra'], linewidth=2, color='C3')
# xlim([0, 10])
# gca().set_xticks(range(0,11,2))
# ylim([-400,1400])
# gca().set_yticks(range(-400,1401,200))
xlabel('time (ns)')
ylabel(r'$\kappa_{out}$ (W/m/K)')
title('(b)')

subplot(2,2,3)
set_fig_properties([gca()])
plot(t, kappa['kxi_ra'], linewidth=2)
plot(t, kappa['kxo_ra'], linewidth=2, color='C3')
plot(t, kappa['kxi_ra']+kappa['kxo_ra'], linewidth=2, color='k')
# xlim([0, 10])
# gca().set_xticks(range(0,11,2))
# ylim([-200,1000])
# gca().set_yticks(range(-200,1001,200))
xlabel('time (ns)')
ylabel(r'$\kappa$ (W/m/K)')
legend(['in', 'out', 'total'])
title('(c)')


subplot(2,2,4)
set_fig_properties([gca()])
plot(t, kappa['kyi_ra']+kappa['kyo_ra'],color='k', linewidth=2)
plot(t, kappa['kxi_ra']+kappa['kxo_ra'], color='C0', linewidth=2)
plot(t, kappa['kz_ra'], color='C3', linewidth=2)
print(f"最后一个x方向的热导率为 ：{kappa['kxi_ra'][-1],kappa['kxo_ra'][-1]}")
# xlim([0, 10])
# gca().set_xticks(range(0,11,2))
# ylim([-500,1000])
# gca().set_yticks(range(-500,1001,300))
xlabel('time (ns)')
ylabel(r'$\kappa$ (W/m/K)')
legend(['yx', 'xx', 'zx'])
title('(d)')

tight_layout()
savefig('hnemd.png')

import numpy as np
# 计算特定时间范围内对应y的所有数据的平均值
# 子图(a)
indices_a = np.where((t >= 6) & (t <= 8))
x_values_a = kappa['kxi_ra'][indices_a]
average_kappa_a = np.mean(x_values_a)

# 子图(b)
indices_b = np.where((t >= 6) & (t <= 8))
x_values_b = kappa['kxo_ra'][indices_b]
average_kappa_b = np.mean(x_values_b)

# 子图(c)
indices_c_in = np.where((t >= 6) & (t <= 8))
indices_c_out = np.where((t >= 6) & (t <= 8))
x_values_c_in = kappa['kxi_ra'][indices_c_in]
x_values_c_out = kappa['kxo_ra'][indices_c_out]
average_kappa_c_in = np.mean(x_values_c_in)
average_kappa_c_out = np.mean(x_values_c_out)

# 子图(d)
indices_d = np.where((t >= 6) & (t <= 8))
x_values_d = kappa['kxi_ra']+kappa['kxo_ra']
average_kappa_d = np.mean(x_values_d)

# 输出平均值到文件
output_file = open('average_kappa_values.txt', 'w')
output_file.write("子图(a) 平均热导率值为： " + str(average_kappa_a) + "\n")
output_file.write("子图(b) 平均热导率值为： " + str(average_kappa_b) + "\n")
output_file.write("子图(c) in 平均热导率值为： " + str(average_kappa_c_in) + "\n")
output_file.write("子图(c) out 平均热导率值为： " + str(average_kappa_c_out) + "\n")
output_file.write("子图(d) 平均热导率值为： " + str(average_kappa_d))
output_file.close()


