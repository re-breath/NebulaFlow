import nebula
import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import ase.io as ai

def elect_rely_force(config,minforce,maxforce):
    """检查二维矩阵中每一个值是否在一定范围内"""
    flag = not np.any(config.force < minforce) and not np.any(config.force > maxforce)
    return flag
def filter_force(configs,minforce,maxforce):
    """筛选出符合要求的构型"""
    fit_file = "force_fited.xyz"
    unfit_file = "force_unfited.xyz"
    
    if os.path.exists(fit_file):
        os.remove(fit_file)
    if os.path.exists(unfit_file):
        os.remove(unfit_file)
    fitindex=[]
    fitconfigs=[]
    unfitconfigs=[]
    for i in range(len(configs)):
        if elect_rely_force(configs[i],minforce,maxforce):
            fitconfigs.append(configs[i])  
            fitindex.append(i)
        else:
            unfitconfigs.append(configs[i])
    # print("Fit configs:",fitindex)
    print("Force range:",np.min(force_all),np.max(force_all))
    print("Force fit number:",len(fitconfigs))
    print("Force unfit number:",len(unfitconfigs))
   
    force_fit = np.concatenate([config.force.flatten() for config in fitconfigs])
    plt.hist(force_fit,bins=100,color=my_favorite_colors[2],alpha=0.8,label='Fit-Force')
    # for i in fitindex:
    #     ai.write(fit_file,atoms[i],append=True)
    nebula.write_xyz(fit_file,fitconfigs)
    nebula.write_xyz(unfit_file,unfitconfigs)
    

readfile = sys.argv[1]

minforce = float(sys.argv[2])
maxforce = float(sys.argv[3])
if len(sys.argv) != 4:
    print("Usage: python elect_rely_force.py readfile minforce maxforce")
    sys.exit(1)
# readfile = "train_E_Vir_rational.xyz"
# maxforce = 5
# minforce = -5

my_favorite_colors = ['#1f77b4', '#a56cc1', '#39bdc8', '#d62728', '#9467bd', '#f5587b', '#fcb1b1', '#cabbe9', '#30e3ca', '#00d1ff']

configs = nebula.read_xyz(readfile)
# atoms = ai.read(readfile)
print(f"Number of configs:",len(configs))
# print(len(configs),type(configs),configs)


force_all = np.concatenate([config.force.flatten() for config in configs])

filter_force(configs,minforce,maxforce)
#print(configs[158].index1)

# force_all = np.array([config.force for config in configs])
# force_all = np.reshape(force_all,(-1))



plt.hist(force_all,bins=100,color=my_favorite_colors[4],alpha=0.5,label='All-Force')
plt.title('Force Distribution')
plt.xlabel('Force (eV/Angstrom)')
plt.ylabel('Frequency')
plt.legend()
plt.savefig('force_distribution.png')