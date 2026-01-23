# 绘制 XRD 图,输入文件为debyer计算后的xrd数据，第一列为2θ，第二列为intensity
# 使用方法：python plot_xrd.py xrd.txt
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_widths
import sys
from pybaselines import Baseline

# -----------------------------
# 读取数据
# -----------------------------

try:
    filename = sys.argv[1]
except IndexError:
    print("请输入文件名")
    exit(1)

data = np.loadtxt(filename, comments='#')
fig, ax = plt.subplots(figsize=(10, 6))

x_init = data[:, 0]   # 2θ
y_init = data[:, 1]   # intensity

mask_cut = (x_init >= 10) & (x_init <= 60)
x_init = x_init[mask_cut]
y_init = y_init[mask_cut]

def find_peak_region(x, y, target_pos, search_window=8):
    """
    自动查找完整峰形范围（基于峰宽）。
    
    参数：
    x : ndarray 横坐标 (2θ)
    y : ndarray 强度
    target_pos : float 目标峰的大致位置 (例如 26.5 对应002峰)
    search_window : float 搜索窗口范围 (默认 ±8°)
    
    返回：
    (start, end) 峰区在原始数组中的索引
    """
    # 1. 在目标位置附近截取数据
    mask = (x >= target_pos - search_window) & (x <= target_pos + search_window)
    x_sub, y_sub = x[mask], y[mask]
    
    # 2. 找峰
    peaks, _ = find_peaks(y_sub, prominence=10)
    if len(peaks) == 0:
        return None
    
    # 3. 找到最接近 target_pos 的峰
    peak_idx_sub = peaks[np.argmin(np.abs(x_sub[peaks] - target_pos))]
    
    # 4. 计算峰宽（半高宽）
    results_half = peak_widths(y_sub, [peak_idx_sub], rel_height=0.98)
    left_idx_sub  = int(results_half[2][0])   # 子数组中的左索引
    right_idx_sub = int(results_half[3][0])   # 子数组中的右索引
    
    # 5. 映射回原数组索引
    origin_inds = np.where(mask)[0]           # mask 为 True 的位置在原数组中的下标
    left_idx  = origin_inds[left_idx_sub]
    right_idx = origin_inds[right_idx_sub]
    
    return left_idx, right_idx


def find_002_peak(x, y):
    """查找002峰完整范围"""
    return find_peak_region(x, y, target_pos=26.5, search_window=8)


def find_100_peak(x, y):
    """查找100峰完整范围"""
    return find_peak_region(x, y, target_pos=42.0, search_window=8)


def find_noise_region(x, y, step=50):
    """
    自动查找噪音区（最平坦的一段）。
    
    参数：
    x : ndarray 横坐标 (2θ)
    y : ndarray 强度
    step : int 区间宽度 (单位：点数)
    
    返回：
    (start, end) 噪音区范围
    """
    noise_std_list = []
    for i in range(0, len(x)-step, step):
        seg = y[i:i+step]
        noise_std_list.append((np.std(seg), i))
    noise_std_list.sort(key=lambda x: x[0])
    noise_start = noise_std_list[0][1]
    return noise_start, noise_start+step-1

# p002 = find_002_peak(x_init, y_init)
# p100 = find_100_peak(x_init, y_init)
# noise = find_noise_region(x_init, y_init)

# plt.plot([x_init[p002[0]], x_init[p002[1]]], [50, 50], color='#D43358', linestyle='--', label='002 Peak')
# plt.plot([x_init[p100[0]], x_init[p100[1]]], [50, 50], color='#61BA2E', linestyle='--', label='100 Peak')
# plt.plot([x_init[noise[0]], x_init[noise[1]]], [50, 50], color='#2F1DA8', linestyle='--', label='Noise Region')

# plt.plot(x_init[p002[0]:p002[1]], y_init[p002[0]:p002[1]], color='#D43358', label='002 Peak')
# plt.plot(x_init[p100[0]:p100[1]], y_init[p100[0]:p100[1]], color='#61BA2E', label='100 Peak')
# plt.plot(x_init[noise[0]:noise[1]], y_init[noise[0]:noise[1]], color='#2F1DA8', label='Noise Region')

# -----------------------------
# 自动查找最优 lam 值
# -----------------------------
def find_optimal_lam(baseline_fitter, y_init, x_init, p=0.002,
                     lam_start=1e3, lam_factor=10, tol=0.1, lam_max=1e9):
    """
    自动查找合适的 lam 值，用峰高保持和 SNR 判据。
    
    参数：
    baseline_fitter : Baseline 对象
    y_init : ndarray 原始信号
    x_init : ndarray 横坐标 (2θ)
    p : float asls 平滑参数
    lam_start : float 初始 lam
    lam_factor : float 每次增加倍数
    tol : float 峰高差距容忍度
    lam_max : float 最大 lam
    peak_region : tuple 峰区范围 (002峰附近)
    noise_region : tuple 噪声区范围
    
    返回：
    lam_opt, y_correct, bkg
    """
    p002 = find_002_peak(x_init, y_init)
    noise = find_noise_region(x_init, y_init)

    peak_region = (x_init[p002[0]], x_init[p002[1]])
    noise_region = (x_init[noise[0]], x_init[noise[1]])

    lam = lam_start
    peak_mask = (x_init >= peak_region[0]) & (x_init <= peak_region[1])
    noise_mask = (x_init >= noise_region[0]) & (x_init <= noise_region[1])
    
    peak_height_init = np.max(y_init[peak_mask]) - np.min(y_init[peak_mask])
    
    best_lam, best_snr = None, -np.inf
    best_ycorr, best_bkg = None, None
    
    while lam <= lam_max:
        bkg, params = baseline_fitter.asls(y_init, lam=lam, p=p)
        y_corr = y_init - bkg
        
        peak_height_corr = np.max(y_corr[peak_mask]) - np.min(y_corr[peak_mask])
        noise_std = np.std(y_corr[noise_mask])
        snr = peak_height_corr / noise_std if noise_std > 0 else 0
        
        diff_ratio = abs(peak_height_corr - peak_height_init) / peak_height_init
        
        print(f"lam: {lam:.0e}, diff_ratio: {diff_ratio:.3f}, snr: {snr:.3f}")

        # 判据：峰高保持在容忍范围内，同时 SNR 最大化
        if diff_ratio <= tol and snr > best_snr:
            best_snr = snr
            best_lam = lam
            best_ycorr = y_corr
            best_bkg = bkg
        
        lam *= lam_factor
    
    return best_lam, best_ycorr, best_bkg


baseline_fitter = Baseline(x_data=x_init)

lam_opt, ycorrect, bkg = find_optimal_lam(baseline_fitter, y_init, x_init)
print(f"optimum fit lambda: {lam_opt:.0e}")


mask_cut = (x_init >= 0) & (x_init <= 180)
x = x_init[mask_cut]
y = ycorrect[mask_cut]




lambda_A = 1.5406  # Cu Kα
K002 = 0.9
K100 = 1.84

def cut_peak(x, y, center_deg, half_width_deg=3, floor=0.02):
    """
    center_deg: 预期峰位
    half_width_deg: 往左右各延伸几度
    floor: 相对高度低于 2 %max 就当作基线
    """
    mask = (x >= center_deg - half_width_deg) & (x <= center_deg + half_width_deg)
    x_, y_ = x[mask], y[mask]
    # 再把两端“已经到基线”的部分各留 0.5° 缓冲
    ymax = y_.max()
    base = y_ < ymax * floor
    # 找最左/最右的 base 点
    left  = np.where(base)[0]
    right = np.where(base)[0]
    if left.size and right.size:
        i1, i2 = max(left[0] - int(0.5/np.diff(x_).mean()), 0), \
                 min(right[-1] + int(0.5/np.diff(x_).mean()), x_.size-1)
        x_, y_ = x_[i1:i2], y_[i1:i2]
    return x_, y_

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


# 3. 主曲线 + 填充 ------------------------------------------------------------
ax.plot(x_init,y_init,color="#2ECC75", lw=2, alpha=0.3, label="XRD init")
ax.plot(x, y, color="#1E4AA0FF", lw=2, label="XRD corrected")
ax.plot(x, bkg, color="#7D8DCF", linestyle='--', lw=2, alpha=0.3, label="baseline")


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

peak_color = ['#FF4070', '#3783F5', '#109B9E']
low_y = -2

p002 = find_002_peak(x, y)
p100 = find_100_peak(x, y)
noise = find_noise_region(x, y)

plt.plot([x[p002[0]], x[p002[1]]], [0, 0], color=peak_color[0], linestyle='--' , alpha=0.5, linewidth=2)
plt.plot([x[p100[0]], x[p100[1]]], [0, 0], color=peak_color[1], linestyle='--' , alpha=0.5, linewidth=2)
plt.plot([x[noise[0]], x[noise[1]]], [0, 0], color=peak_color[2], linestyle='--', alpha=0.5, linewidth=2)

x_mid_002 = (x[p002[0]] + x[p002[1]]) / 2
plt.text(x_mid_002, low_y, "002 Peak", color=peak_color[0],
         ha='center', va='top', fontsize=10)

x_mid_100 = (x[p100[0]] + x[p100[1]]) / 2
plt.text(x_mid_100, low_y, "100 Peak", color=peak_color[1],
         ha='center', va='top', fontsize=10)

x_mid_noise = (x[noise[0]] + x[noise[1]]) / 2
plt.text(x_mid_noise, low_y, "Noise Region", color=peak_color[2],
         ha='center', va='top', fontsize=10)



ax.text(0.82, 0.95, info_txt, transform=ax.transAxes,fontweight='bold',
        fontsize=12, verticalalignment='top',color='#5878F8FF',linespacing=1.5, 
        bbox=dict(boxstyle="round,pad=0.4",
                  facecolor='#FAFCBEFF', alpha=0.85,
                  edgecolor='0.5'))

# 6. 坐标轴 & 标题 -------------------------------------------------------------
ax.set_xlabel(r"2$\mathbf{\theta}$ (deg)", fontsize=15,fontweight='bold')
ax.set_ylabel("Intensity (a.u.)", fontsize=15,fontweight='bold')
ax.set_title("XRD Pattern of Carbon Fiber", fontsize=17, fontweight='bold', pad=15)
ax.set_xlim(10, 60)
#ax.set_ylim(0,100)

# 7. 细节 ----------------------------------------------------------------------
ax.tick_params(axis='both', labelsize=12)
fig.tight_layout()
ax.xaxis.grid(True, which='major', linestyle='-', linewidth=1.5, color='#0B5D56FF',alpha=0.25)
ax.yaxis.grid(False)          # 确保 y 轴网格不出现
ax.legend(fontsize=12, loc='upper left')
fig.savefig("xrd_baseline.png", dpi=400, bbox_inches='tight')
