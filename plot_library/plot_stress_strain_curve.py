import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# 全局字体设置（统一放大）
plt.rcParams.update({
    "font.size": 18,         
    "axes.labelsize": 20,    
    "axes.titlesize": 22,     
    "xtick.labelsize": 16,   
    "ytick.labelsize": 16,  
    "legend.fontsize": 16    
})

def sample_stress_strain(strain, stress, target_points: int = 1000):
    """
    均匀采样应力应变数据，保证绘图效率,如果数据量过多，该函数将会对数据进行均匀采样
    :param strain: 原始应变数组
    :param stress: 原始应力数组
    :param target_points: 目标采样点数（默认1000，最优）
    :return: 筛选后的 strain, stress
    """
    # 获取数据总长度
    total_points = len(strain)
    
    # 如果数据量 ≤ 目标点数，直接返回全部数据
    if total_points <= target_points:
        return strain.copy(), stress.copy()
    
    # 数据量过大：均匀等间隔采样（保留首尾点，曲线完整）
    indices = np.linspace(0, total_points - 1, target_points, dtype=int)
    # 去重（防止极端情况重复索引）
    indices = np.unique(indices)
    
    return strain[indices], stress[indices]

# 读取数据
data = np.loadtxt("stress_strain_curve.txt")
strain = data[:, 0]
stress = data[:, 1]

strain, stress = sample_stress_strain(strain, stress)

import numpy as np
from scipy.stats import linregress

def best_linear_window(strain, stress,
                       eps_min=0.0, eps_max=0.05,
                       min_points=20,
                       min_span=0.005,
                       r2_threshold=0.995):
    """
    在 [eps_min, eps_max] 里找最佳线性窗口，返回:
    (best_slope, best_intercept, best_r2, (i, j))
    其中窗口是 strain[i:j]（j不含）。
    """
    # 只在小应变段搜索
    mask = (strain >= eps_min) & (strain <= eps_max)
    s = strain[mask]
    sig = stress[mask]

    if s.size < min_points:
        raise ValueError("Not enough points in the search range.")

    best = None  # (score, slope, intercept, r2, i, j)
    n = s.size

    for i in range(0, n - min_points):
        for j in range(i + min_points, n + 1):
            span = s[j-1] - s[i]
            if span < min_span:
                continue

            slope, intercept, r, p, se = linregress(s[i:j], sig[i:j])
            r2 = r * r

            if r2 < r2_threshold:
                continue

            # RMSE 作为第二指标，避免“R2高但窗口太窄/偶然”
            pred = intercept + slope * s[i:j]
            rmse = float(np.sqrt(np.mean((sig[i:j] - pred) ** 2)))

            score = r2 - 0.1 * rmse  # 0.1 是经验权重，可调

            if best is None or score > best[0]:
                best = (score, slope, intercept, r2, i, j, rmse)

    if best is None:
        raise ValueError("No suitable linear window found. Try lowering r2_threshold or widening eps_max.")

    _, slope, intercept, r2, i, j, rmse = best
    return slope, intercept, r2, (i, j), rmse


tensile_strength = np.max(stress)
tensile_index = np.argmax(stress)
tensile_strain = strain[tensile_index]

#该函数能用就用，不能就跳过
try:
    slope, intercept, r2, (i, j), rmse = best_linear_window(strain, stress,
                                                        eps_min=0.0, eps_max=0.05,
                                                        min_points=20, min_span=0.005,
                                                        r2_threshold=0.995)
except ValueError as e:
    print(e)
    print("warning: No suitable linear window found.")
    slope, intercept, r2, (i, j), rmse = 0.0, 0.0, 0.0, (0, 0), 0.0
    

youngs_modulus = slope
linear_strain = strain[(strain >= 0.0) & (strain <= 0.05)][i:j]
linear_stress = stress[(strain >= 0.0) & (strain <= 0.05)][i:j]

print(f"Best linear window: strain[{linear_strain[0]:.4f}, {linear_strain[-1]:.4f}]")
print(f"Young's Modulus: {youngs_modulus:.6f} GPa (R²={r2:.6f}, RMSE={rmse:.6g})")


fig, ax = plt.subplots(figsize=(10, 7))  # 图尺寸稍微加大

ax.axvspan(strain.min(), tensile_strain, color="#c6dbef", alpha=0.3)
ax.axvspan(tensile_strain, strain.max(), color="#fdd0a2", alpha=0.3)
plt.xlim(strain.min(), strain.max())


ax.plot(strain, stress, color="#545557FF", linewidth=2.5, label="Stress-Strain Curve")

ax.plot(linear_strain, intercept + slope*linear_strain,
        label=f"Linear Fit (E={youngs_modulus:.2f} GPa)",
        color="#b30000", linestyle="--", linewidth=2)

# 抗拉强度点
ax.scatter(tensile_strain, tensile_strength,
           color="#FFBFBDFF", s=150, edgecolor="#5C37CAFF", zorder=5)

# 计算x,y轴最大长度
x_max = max(strain.max(), tensile_strain)
y_max = max(stress.max(), tensile_strength)

ax.annotate(f"Tensile Strength\n{tensile_strength:.2f} GPa",
            xy=(tensile_strain, tensile_strength),
            xytext=(tensile_strain+x_max*0.08, tensile_strength-y_max*0.1),
            fontsize=16, color="#333333",
            arrowprops=dict(arrowstyle="->", color="#555555", lw=2))

ax.set_xlabel("Strain", fontsize=20, weight="bold")
ax.set_ylabel("Stress (GPa)", fontsize=20, weight="bold")

ax.legend(frameon=True, fancybox=True, loc="best")
ax.grid(True, linestyle=":", linewidth=1.2, alpha=0.7)

for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

plt.tight_layout()
plt.savefig("stress_strain_curve.png", dpi=300)
# plt.show()

print(f"Young's Modulus: {youngs_modulus:.6f} GPa")
print(f"Tensile Strength: {tensile_strength:.6f} GPa")
