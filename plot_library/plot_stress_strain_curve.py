import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# 全局字体设置（统一放大）
plt.rcParams.update({
    "font.size": 18,          # 默认字体大小
    "axes.labelsize": 20,     # 坐标轴标签字体
    "axes.titlesize": 22,     # 标题字体
    "xtick.labelsize": 16,    # x轴刻度字体
    "ytick.labelsize": 16,    # y轴刻度字体
    "legend.fontsize": 16     # 图例字体
})

# 读取数据
data = np.loadtxt("stress_strain_curve.txt")
strain = data[:, 0]
stress = data[:, 1]

# 计算抗拉强度
tensile_strength = np.max(stress)
tensile_index = np.argmax(stress)
tensile_strain = strain[tensile_index]

# 手动指定线性区间
mask = (strain >= 0.01) & (strain <= 0.06)
linear_strain = strain[mask]
linear_stress = stress[mask]

# 线性拟合
slope, intercept, r_value, p_value, std_err = linregress(linear_strain, linear_stress)
youngs_modulus = slope

# 绘图
fig, ax = plt.subplots(figsize=(10, 7))  # 图尺寸稍微加大

# 背景填充：左蓝右橙
ax.axvspan(strain.min(), tensile_strain, color="#c6dbef", alpha=0.3)
ax.axvspan(tensile_strain, strain.max(), color="#fdd0a2", alpha=0.3)
plt.xlim(strain.min(), strain.max())

# 主曲线
ax.plot(strain, stress, color="#545557FF", linewidth=2.5, label="Stress-Strain Curve")

# 杨氏模量拟合线
ax.plot(linear_strain, intercept + slope*linear_strain,
        label=f"Linear Fit (E={youngs_modulus:.2f} GPa)",
        color="#b30000", linestyle="--", linewidth=2)

# 抗拉强度点
ax.scatter(tensile_strain, tensile_strength,
           color="#FFBFBDFF", s=150, edgecolor="#5C37CAFF", zorder=5)

# 计算x,y轴最大长度
x_max = max(strain.max(), tensile_strain)
y_max = max(stress.max(), tensile_strength)

# 注释字体加大
ax.annotate(f"Tensile Strength\n{tensile_strength:.2f} GPa",
            xy=(tensile_strain, tensile_strength),
            xytext=(tensile_strain+x_max*0.08, tensile_strength-y_max*0.1),
            fontsize=16, color="#333333",
            arrowprops=dict(arrowstyle="->", color="#555555", lw=2))

# 坐标轴标签
ax.set_xlabel("Strain", fontsize=20, weight="bold")
ax.set_ylabel("Stress (GPa)", fontsize=20, weight="bold")

# 图例美化
ax.legend(frameon=True, fancybox=True, loc="best")

# 网格优化
ax.grid(True, linestyle=":", linewidth=1.2, alpha=0.7)

# 去掉顶部和右侧边框
for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

plt.tight_layout()
plt.savefig("stress_strain_curve_bigfont.png", dpi=300)
# plt.show()

print(f"Young's Modulus: {youngs_modulus:.2f} GPa")
print(f"Tensile Strength: {tensile_strength:.2f} GPa")
