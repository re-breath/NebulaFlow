# 该脚本将会探究LJ势能与力在不同距离下的变化，并进行画图
# 该脚本将会指定epslion和sigma

import matplotlib.pyplot as plt
import numpy as np
import os

epsilon = 0.01032  
sigma = 3.405 


def lj_potential(r):
    return 4 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)

def lj_force(r):
    return -24 * epsilon * (2 * (sigma / r) ** 12 - (sigma / r) ** 6) / r


# 生成0.5埃到15埃的距离数据
r = np.linspace(0.5, 15, 1000)

# 计算势能和力
V = lj_potential(r)
F = lj_force(r)

plt.figure(figsize=(12, 6), dpi=150)

plt.subplot(2, 1, 1)
plt.plot(r, V, color="#1B72B0FF", linewidth=2.4)
plt.xlabel("Distance ($\sigma$)")
plt.ylabel("Potential Energy ($\epsilon$)")
plt.title("Lennard-Jones Potential")
plt.xlim(0.5, 15)
plt.ylim(-1.5, 1.5)
plt.xticks(np.arange(0, 16, 2))
plt.yticks(np.arange(-2, 2.1, 0.5))
plt.grid(True, linestyle="--", alpha=0.3)

plt.subplot(2, 1, 2)
plt.plot(r, F, color="#D95F02FF", linewidth=2.4)
plt.xlabel("Distance ($\sigma$)")
plt.ylabel("Force ($\epsilon/\sigma$)")
plt.title("Lennard-Jones Force")
plt.xlim(0.5, 15)
plt.ylim(-2.5, 2.5)
plt.xticks(np.arange(0, 16, 2))
plt.yticks(np.arange(-3, 3.1, 1))
plt.grid(True, linestyle="--", alpha=0.3)

plt.tight_layout()
plt.legend(["$V(r)$", "$F(r)$"], fontsize=12)

plt.savefig('lj.png', dpi=300)