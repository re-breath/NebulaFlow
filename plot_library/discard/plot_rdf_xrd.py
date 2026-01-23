# 处理xrd.exe输出的SFACTOR_L1_L2_na.data文件
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# 读取 RDF
r, g = np.loadtxt("RDF_L1_L2_na.data", unpack=True)
# 读取 SFACTOR
theta, I = np.loadtxt("SFACTOR_L1_L2_na.data", unpack=True)

# 设置统一风格
plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams.update({
    "font.size": 14,
    "axes.labelsize": 16,
    "axes.titlesize": 18,
    "legend.fontsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "figure.dpi": 300
})

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# RDF 图
axes[0].plot(r, g, color="tab:blue", linewidth=2, label="g(r)")
axes[0].set_xlabel("r (Å)")
axes[0].set_ylabel("g(r)")
axes[0].set_title("Radial Distribution Function")
axes[0].legend()
axes[0].grid(True, linestyle="--", alpha=0.6)
axes[0].set_xlim(left=0)

# SFACTOR 图
axes[1].plot(theta, I, color="tab:red", linewidth=2, label="I(2θ)")
axes[1].set_xlabel("2θ (degrees)")
axes[1].set_ylabel("Intensity (a.u.)")
axes[1].set_title("XRD Pattern")
axes[1].legend()
axes[1].grid(True, linestyle="--", alpha=0.6)
axes[1].set_xlim(left=0)
axes[1].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))

plt.tight_layout()
plt.savefig("rdf_xrd.png", dpi=600)
plt.show()