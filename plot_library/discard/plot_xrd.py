# 绘制 XRD 图,输入文件为debyer计算后的xrd数据，第一列为2θ，第二列为intensity
# 使用方法：python plot_xrd.py xrd.txt
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import sys

# -----------------------------
# 读取数据
# -----------------------------

try:
    filename = sys.argv[1]
except IndexError:
    print("请输入文件名")
    exit(1)

data = np.loadtxt(filename, comments='#')
x_init = data[:, 0]   # 2θ
y_init = data[:, 1]   # intensity

mask_cut = (x_init >= 10) & (x_init <= 90)
x = x_init[mask_cut]
y = y_init[mask_cut]

lambda_A = 1.5406  # Cu Kα
K002 = 0.9
K100 = 1.84


def calc_fwhm(x, y, peak_index):
    """
    x: 2θ 数组
    y: intensity 数组
    peak_index: 峰位置的索引
    返回: FWHM, 区间数据 (x_mid, y_mid) —— 半高宽区间的中间一半
    """
    I_max = y[peak_index]
    half = I_max / 2

    # 左侧
    left = peak_index
    while left > 0 and y[left] > half:
        left -= 1

    # 右侧
    right = peak_index
    while right < len(y)-1 and y[right] > half:
        right += 1

    # 线性插值提高精度
    x_left = np.interp(half, [y[left], y[left+1]], [x[left], x[left+1]])
    x_right = np.interp(half, [y[right-1], y[right]], [x[right-1], x[right]])

    fwhm = x_right - x_left

    # 提取半高宽区间的数据
    mask = (x >= x_left) & (x <= x_right)
    x_seg = x[mask]
    y_seg = y[mask]

    # 只取中间一半的数据
    n = len(x_seg)
    start = n // 4
    end = 3 * n // 4
    x_mid = x_seg[start:end]
    y_mid = y_seg[start:end]

    return fwhm, x_mid, y_mid



def fit_peak(x_seg, y_seg, kind="quadratic"):
    """
    对半高宽区间的数据进行拟合
    kind: "quadratic" 二次多项式拟合, "gaussian" 高斯拟合
    返回: 拟合函数 poly 或参数
    """
    if kind == "quadratic":
        coeffs = np.polyfit(x_seg, y_seg, 2)
        poly = np.poly1d(coeffs)
        return poly
    elif kind == "gaussian":
        from scipy.optimize import curve_fit

    def gauss(x, A, x0, sigma, baseline):
        return A * np.exp(-(x - x0)**2 / (2*sigma**2)) + baseline

    A0 = y_seg.max() - y_seg.min()
    x0 = x_seg[np.argmax(y_seg)]
    sigma0 = (x_seg[-1] - x_seg[0]) / 2
    baseline0 = y_seg.min()

    popt, _ = curve_fit(gauss, x_seg, y_seg, p0=[A0, x0, sigma0, baseline0])

    # 返回一个可调用函数
    return lambda xx: gauss(xx, *popt)


def quadratic_peak_x(poly):
    """
    计算二次函数的顶点横坐标
    poly: np.poly1d 对象 或 系数数组 [a, b, c]
    返回: 顶点横坐标 x_peak
    """
    # 如果传入的是 poly1d
    if isinstance(poly, np.poly1d):
        coeffs = poly.coefficients
    else:
        coeffs = np.array(poly)

    a, b = coeffs[0], coeffs[1]
    if a == 0:
        raise ValueError("不是二次函数，无法计算顶点")
    x_peak = -b / (2 * a)
    return x_peak


# -----------------------------
# 1. 找 002 峰（18–32°）
# -----------------------------
mask_002 = (x >= 18) & (x <= 32)
x_002 = x[mask_002]
y_002 = y[mask_002]

peaks_002, _ = find_peaks(y_002, prominence=1)

if len(peaks_002) == 0:
    print('⚠️  没找到(002)峰，请检查数据范围或阈值')
    main_002 = 0
    d002 = 0
    Lc = 0
else:
    main_002 = peaks_002[np.argmax(y_002[peaks_002])]

if main_002 is None:
    print('⚠️  跳过 FWHM 计算：没有主 (002) 峰')
    beta_002 = x_seg_002 = y_seg_002 = None
else:
    beta_002, x_seg_002, y_seg_002 = calc_fwhm(x_002, y_002, main_002)

    poly_002 = fit_peak(x_seg_002, y_seg_002, kind="quadratic")

    y_seg_fit_002 = poly_002(x_seg_002)
    two_theta_002_fit = quadratic_peak_x(poly_002)
    I_002_fit = poly_002(two_theta_002_fit)

    print(f"FWHM_002 = {beta_002:.3f} °")
    beta_002_rad = np.deg2rad(beta_002)

    two_theta_002 = x_002[main_002]
    I_002 = y_002[main_002]

    # 计算 d002  
    theta_002 = np.deg2rad(two_theta_002_fit / 2)
    d002 = lambda_A / (2 * np.sin(theta_002))

    # 计算 Lc（Scherrer）

    Lc = (K002 * lambda_A) / (beta_002_rad * np.cos(theta_002))
    print(f"Peak 002 Fit = {two_theta_002_fit:.3f} °")
    print(f"d₀₀₂ = {d002:.3f} Å")
    print(f"Lc = {Lc:.2f} Å")


# -----------------------------
# 2. 找 100 峰（38–48°）
# -----------------------------
mask_100 = (x >= 38) & (x <= 48)
x_100 = x[mask_100]
y_100 = y[mask_100]

peaks_100, _ = find_peaks(y_100, prominence=1)
if len(peaks_100) == 0:
    print('⚠️  没找到(100)峰，请检查数据范围或阈值')
    main_100 = 0
    La = 0
else:
    main_100 = peaks_100[np.argmax(y_100[peaks_100])]

if main_100 is None:
    print('⚠️  跳过 FWHM 计算：没有主 (100) 峰')
    beta_100 = x_seg_100 = y_seg_100 = None
else:
    beta_100, x_seg_100, y_seg_100 = calc_fwhm(x_100, y_100, main_100) 

    poly_100 = fit_peak(x_seg_100, y_seg_100, kind="quadratic")
    y_seg_fit_100 = poly_100(x_seg_100)
    two_theta_100_fit = quadratic_peak_x(poly_100)
    I_100_fit = poly_100(two_theta_100_fit)

    print(f"FWHM_100 = {beta_100:.3f} °")
    beta_100_rad = np.deg2rad(beta_100)

    two_theta_100 = x_100[main_100]
    I_100 = y_100[main_100]

    # 计算 La（Scherrer）
    theta_100 = np.deg2rad(two_theta_100_fit / 2)

    La = (K100 * lambda_A) / (beta_100_rad * np.cos(theta_100))

    print(f"Peak 100 Fit = {two_theta_100_fit:.3f} °")

    print(f"La = {La:.2f} Å")


# -----------------------------
# 4. 绘图
# -----------------------------

# 2. 画布 ---------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 6))

# 3. 主曲线 + 填充 ------------------------------------------------------------
ax.plot(x, y, color="#1E4AA0FF", lw=2, label="XRD pattern")
ax.fill_between(x, y, color="#35B17FFF", alpha=0.12, zorder=3)

# 4. peak 标记 ----------------------------------------------------------------
# 002
if main_002 is not None:
    ax.scatter(two_theta_002_fit, I_002_fit, color="#CD5F5FFF", s=100,
            edgecolor="white", linewidth=1.5, zorder=6)
    ax.annotate(f"002 : {two_theta_002_fit:.2f}°",
                xy=(two_theta_002_fit+0.5, I_002_fit+0.5),
                xytext=(two_theta_002_fit + 2.2, I_002_fit * 0.95),
                color="#CD5F5FFF", fontweight="bold", fontsize=13,
                arrowprops=dict(arrowstyle="->", color="#CD5F5FFF", lw=1.5))
    ax.plot(x_seg_002, y_seg_fit_002, color="#E66B6BFF", lw=2, label="Fit 002")



# 100
if main_100 is not None:
    ax.scatter(two_theta_100_fit, I_100_fit, color="#6854B7FF", s=100,
           edgecolor="white", linewidth=1.5, zorder=6)
    ax.annotate(f"100 : {two_theta_100_fit:.2f}°",
                xy=(two_theta_100_fit, I_100_fit+0.5),
                xytext=(two_theta_100_fit - 1.5, I_100_fit * 1.28),
                color="#6854B7FF", fontweight="bold", fontsize=13,
                arrowprops=dict(arrowstyle="->", color="#6854B7FF", lw=1.5))
    ax.plot(x_seg_100, y_seg_fit_100, color="#EB6E6EFF", lw=2, label="Fit 100")

# 5. 信息卡片 ------------------------------------------------------------------
outinfo_txt = (f"d002           {d002:.6f} Å\n"  # 单位 Å
               f"Lc             {Lc:.6f} Å\n"
               f"La             {La:.6f} Å\n"
               f"Peak_002       {two_theta_002_fit:.6f} °\n"
               f"Peak_100       {two_theta_100_fit:.6f} °\n"
               f"FWHM_002       {beta_002:.6f} °\n"
               f"FWHM_100       {beta_100:.6f} °"
               )

info_txt = (
    f"$\\mathbf{{d_{{002}}}}$ = {d002:.3f} Å\n"
    f"$\\mathbf{{L_c}}$   = {Lc:.2f} Å\n"
    f"$\\mathbf{{L_a}}$   = {La:.2f} Å"
)

# 输入到文件中
with open("xrd_info.txt", "w") as f:
    f.write(outinfo_txt)


ax.text(0.82, 0.95, info_txt, transform=ax.transAxes,fontweight='bold',
        fontsize=12, verticalalignment='top',color='#5878F8FF',linespacing=1.5, 
        bbox=dict(boxstyle="round,pad=0.4",
                  facecolor='#FAFCBEFF', alpha=0.85,
                  edgecolor='0.5'))

# 6. 坐标轴 & 标题 -------------------------------------------------------------
ax.set_xlabel(r"2$\mathbf{\theta}$ (deg)", fontsize=15,fontweight='bold')
ax.set_ylabel("Intensity (a.u.)", fontsize=15,fontweight='bold')
ax.set_title("XRD Pattern of Carbon Fiber", fontsize=17, fontweight='bold', pad=15)
ax.set_xlim(10, 50)
ax.set_ylim(bottom=0)

# 7. 细节 ----------------------------------------------------------------------
ax.tick_params(axis='both', labelsize=12)
fig.tight_layout()
ax.xaxis.grid(True, which='major', linestyle='-', linewidth=1.5, color='#0B5D56FF',alpha=0.25)
ax.yaxis.grid(False)          # 确保 y 轴网格不出现
fig.savefig("xrd.png", dpi=400, bbox_inches='tight')
