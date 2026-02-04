# 绘制 XRD 图,输入文件为debyer计算后的xrd数据，第一列为2θ，第二列为intensity
# 该版本为最终版本，使用viogt方法拟合了002，100峰，并且处理了多峰叠加的情况
# 经过该脚本的测试，debyer计算出来的xrd比xrd.exe计算得到的效果更好，更合理
# 使用方法：python plot_xrd.py xrd.txt

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, peak_widths
import sys
from pybaselines import Baseline
from scipy.optimize import curve_fit


lambda_A = 1.5406  # Cu Kα
K002 = 0.9
K100 = 1.84



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

def find_peak_region(x, y, target_pos, search_window=12):
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
    print(f"查找范围: {target_pos - search_window:.2f}° 到 {target_pos + search_window:.2f}°")
    
    # 2. 找峰
    peaks, _ = find_peaks(y_sub, prominence=1)

    if len(peaks) == 0:
        return None
    
    # 3. 找到最接近 target_pos 的峰
    # peak_idx_sub = peaks[np.argmin(np.abs(x_sub[peaks] - target_pos))]
    # 找到峰值最高的峰
    peak_idx_sub = peaks[np.argmax(y_sub[peaks])]
    
    # 4. 计算峰宽（半高宽）
    results_half = peak_widths(y_sub, [peak_idx_sub], rel_height=1)
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
    return find_peak_region(x, y, target_pos=43.0, search_window=8)


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

def extract_dominant_peak(x, y, start_idx, end_idx, min_peak_prominence=0.1):
    """
    从多个叠加峰中提取最大的主峰
    
    参数:
        x, y: 完整数据数组
        start_idx, end_idx: 峰区域的起始和结束索引
        min_peak_prominence: 最小峰突出度，用于寻峰
        
    返回:
        new_start_idx, new_end_idx: 主峰的索引范围
    """
    # 提取峰区域数据
    x_region = x[start_idx:end_idx]
    y_region = y[start_idx:end_idx]
    
    # 1. 找到所有峰
    peaks_idx, properties = find_peaks(
        y_region, 
        prominence=min_peak_prominence,
        distance=5  # 最小峰间距，避免噪声
    )
    
    if len(peaks_idx) == 0:
        print("未检测到峰，返回原始区域")
        return start_idx, end_idx
    
    print(f"检测到 {len(peaks_idx)} 个峰:")
    for i, peak in enumerate(peaks_idx):
        print(f"  峰{i+1}: 位置={x_region[peak]:.2f}°, 强度={y_region[peak]:.2f}")
    
    # 2. 找到最高的峰
    max_peak_idx = peaks_idx[np.argmax(y_region[peaks_idx])]
    max_peak_x = x_region[max_peak_idx]
    max_peak_y = y_region[max_peak_idx]
    print(f"主峰: 位置={max_peak_x:.2f}°, 强度={max_peak_y:.2f}")
    
    # 3. 寻找主峰两侧的谷底（局部最小值）
    # 在主峰左侧找局部最小值
    left_min_idx = max_peak_idx
    for i in range(max_peak_idx - 1, 0, -1):
        if y_region[i] < y_region[left_min_idx]:
            left_min_idx = i
        # 如果开始上升，停止
        elif i > 0 and y_region[i] > y_region[i-1]:
            break
    
    # 在主峰右侧找局部最小值
    right_min_idx = max_peak_idx
    for i in range(max_peak_idx + 1, len(y_region) - 1):
        if y_region[i] < y_region[right_min_idx]:
            right_min_idx = i
        # 如果开始上升，停止
        elif i < len(y_region) - 1 and y_region[i] > y_region[i+1]:
            break
    
    # 4. 扩展边界以确保包含完整峰形
    # 从谷底向外扩展一定的点数，确保包含峰的尾部
    left_extend = max(5, int((right_min_idx - left_min_idx) * 0.01))  # 扩展20%的峰宽
    right_extend = max(5, int((right_min_idx - left_min_idx) * 0.01))
    
    new_start = max(0, left_min_idx - left_extend)
    new_end = min(len(y_region) - 1, right_min_idx + right_extend)
    
    # 转换为原始数据的索引
    new_start_idx = start_idx + new_start
    new_end_idx = start_idx + new_end
    
    print(f"主峰提取范围: {x[new_start_idx]:.2f}° 到 {x[new_end_idx]:.2f}°")
    print(f"数据点数: {new_end_idx - new_start_idx}")
    
    return new_start_idx, new_end_idx

# -----------------------------
# 自动查找最优 lam 值
# -----------------------------
def find_optimal_lam(baseline_fitter, y_init, x_init, p=0.002,
                     lam_start=1e3, lam_factor=10, tol=0.3, lam_max=1e9):
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


print(f"====================>基线拟合...")

baseline_fitter = Baseline(x_data=x_init)

lam_opt, ycorrect, bkg = find_optimal_lam(baseline_fitter, y_init, x_init)
print(f"optimum fit lambda: {lam_opt:.0e}")

x = x_init
y = ycorrect

# -----------------------------
# 1. 找 002 峰（18–32°）
# -----------------------------

def pseudo_voigt(x, A, x0, fwhm, eta, bg):
        sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
        gauss = np.exp(-(x - x0)**2 / (2 * sigma**2))
        lorentz = 1.0 / (1.0 + ((x - x0) / (fwhm / 2.0))**2)
        return A * (eta * lorentz + (1 - eta) * gauss) + bg

def fit_peak_with_initial_guess(x, y, peaks):
    
    # 提取峰区域数据
    x_peak = x[peaks[0]:peaks[1]]
    y_peak = y[peaks[0]:peaks[1]]
    
    print(f"拟合数据范围: x从{x_peak[0]:.2f}到{x_peak[-1]:.2f}")
    print(f"数据点数: {len(x_peak)}")
    
    # 计算更合理的初始猜测值
    # 1. 峰中心: 寻找y的最大值点
    idx_max = np.argmax(y_peak)
    x0_guess = x_peak[idx_max]
    
    # 2. 振幅: 最大值减去估计的背景
    y_max = y_peak[idx_max]
    # 估计背景为数据两端点的平均值
    bg_guess = (y_peak[0] + y_peak[-1]) / 2
    A_guess = y_max - bg_guess
    
    # 3. 半高宽: 估算为数据范围的一部分
    x_range = x_peak[-1] - x_peak[0]
    fwhm_guess = x_range * 0.4  # 假设峰宽约为数据范围的40%
    
    # 4. 混合参数: 设为中间值
    eta_guess = 0.5
    
    # 改进的初始猜测
    initial_guess = [
        A_guess,      # 振幅 - 基于实际数据计算
        x0_guess,     # 峰中心 - 最大值位置
        fwhm_guess,   # 半高宽 - 数据范围的比例
        eta_guess,    # 混合参数
        bg_guess      # 背景
    ]
    
    print(f"自动计算的初始猜测:")
    print(f"  振幅 A = {initial_guess[0]:.2f}")
    print(f"  峰中心 x0 = {initial_guess[1]:.2f}°")
    print(f"  半高宽 fwhm = {initial_guess[2]:.2f}°")
    print(f"  混合参数 eta = {initial_guess[3]:.3f}")
    print(f"  背景 bg = {initial_guess[4]:.2f}")
    
    try:
        # 增加maxfev（最大函数调用次数）并设置边界
        bounds = ([0, x_peak[0], 0.1, 0, min(y_peak)], 
                  [A_guess*10, x_peak[-1], x_range, 1, max(y_peak)])
        
        popt_002, _ = curve_fit(
            pseudo_voigt, 
            x_peak, 
            y_peak,
            p0=initial_guess,
            bounds=bounds,
            maxfev=5000  # 增加最大迭代次数
        )
        
        print("拟合成功!")
        print(f"拟合参数:")
        print(f"  振幅 A = {popt_002[0]:.2f}")
        print(f"  峰中心 x0 = {popt_002[1]:.2f}°")
        print(f"  半高宽 fwhm = {popt_002[2]:.2f}°")
        print(f"  混合参数 eta = {popt_002[3]:.3f}")
        print(f"  背景 bg = {popt_002[4]:.2f}")
        
        # 计算R²值
        residuals = y_peak - pseudo_voigt(x_peak, *popt_002)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y_peak - np.mean(y_peak))**2)
        r_squared = 1 - (ss_res / ss_tot)
        print(f"  R² = {r_squared:.4f}")
        
        return popt_002, True , r_squared
        
    except Exception as e:
        print(f"拟合失败: {e}")
        
        # 尝试更简单的方法：只拟合高斯峰
        print("\n尝试简单的高斯拟合...")
        
        def simple_gaussian(x, A, x0, sigma, bg):
            return A * np.exp(-(x - x0)**2 / (2 * sigma**2)) + bg
        
        initial_gauss = [A_guess, x0_guess, fwhm_guess/2.355, bg_guess]  # 高斯: fwhm = 2.355*sigma
        
        try:
            popt_gauss, _ = curve_fit(
                simple_gaussian,
                x_peak,
                y_peak,
                p0=initial_gauss,
                maxfev=5000
            )
            
            print("高斯拟合成功!")
            print(f"高斯拟合参数:")
            print(f"  振幅 A = {popt_gauss[0]:.2f}")
            print(f"  峰中心 x0 = {popt_gauss[1]:.2f}°")
            print(f"  标准差 sigma = {popt_gauss[2]:.2f}°")
            print(f"  背景 bg = {popt_gauss[3]:.2f}")
            
            # 将高斯参数转换为伪Voigt格式（eta=0表示纯高斯）
            popt_voigt = [popt_gauss[0], popt_gauss[1], popt_gauss[2]*2.355, 0, popt_gauss[3]]
            return popt_voigt, True
            
        except Exception as e2:
            print(f"高斯拟合也失败: {e2}")
            return None, False, 0

print(f"\n====================>拟合 002 峰...")

peaks_002 = find_002_peak(x, y)

peaks_002 = extract_dominant_peak(x, y, peaks_002[0], peaks_002[1])

# 拟合 002 峰
popt_002, success ,r2_002 = fit_peak_with_initial_guess(x, y, peaks_002)

if success:
    #plt.plot(x[peaks_002[0]:peaks_002[1]], y[peaks_002[0]:peaks_002[1]], label='002 Peak init', alpha=0.7, color='#9546CA')
    p002fit_y = pseudo_voigt(x[peaks_002[0]:peaks_002[1]], *popt_002)
    # plt.plot(x[peaks_002[0]:peaks_002[1]], 
    #          p002fit_y, 
    #          linestyle='--', label='002 Peak fit', 
    #          alpha=0.7,color='#9546CA')

    # 计算 d002  
    two_theta_002_fit = popt_002[1]
    two_theta_002_fit_rad = np.deg2rad(popt_002[1])
    theta_002 = two_theta_002_fit_rad / 2

    d002 = lambda_A / (2 * np.sin(theta_002))
    
    beta_002 = popt_002[2]
    beta_002_rad = np.deg2rad(popt_002[2])

    I_002_fit = max(p002fit_y)

    Lc = (K002 * lambda_A) / (beta_002_rad * np.cos(theta_002))

    # 计算 Lc（Scherrer）
    
    print(f"\n002 峰处理结果:")
    print(f"Peak 002 Fit = {two_theta_002_fit:.3f} °")
    print(f"d₀₀₂ = {d002:.3f} Å")
    print(f"Lc = {Lc:.2f} Å")


# -----------------------------
# 2. 找 100 峰（38–48°）
# -----------------------------

print(f"\n====================>拟合 100 峰...")
peaks_100 = find_100_peak(x, y)
print(peaks_100)
peaks_100 = extract_dominant_peak(x, y, peaks_100[0], peaks_100[1])

# 使用提取的主峰进行拟合
popt_100, success ,r2_100 = fit_peak_with_initial_guess(x, y, peaks_100)

if success:
    #plt.plot(x[peaks_100[0]:peaks_100[1]], y[peaks_100[0]:peaks_100[1]], label='100 Peak init')
    p100fit_y = pseudo_voigt(x[peaks_100[0]:peaks_100[1]], *popt_100)
    # plt.plot(x[peaks_100[0]:peaks_100[1]], 
    #          p100fit_y, 
    #          linestyle='--', label='100 Peak fit', 
    #          alpha=0.7, color='#B70000')

    beta_100 = popt_100[2]    
    beta_100_rad = np.deg2rad(popt_100[2])
    two_theta_100_fit = popt_100[1]
    theta_100 = np.deg2rad(two_theta_100_fit) / 2

    I_100_fit = max(p100fit_y)

    # 计算 La（Scherrer）
    La = (K100 * lambda_A) / (beta_100_rad * np.cos(theta_100))

    print(f"\n100 峰处理结果:")
    print(f"Peak 100 Fit = {two_theta_100_fit:.3f} °")
    print(f"La = {La:.2f} Å")


# -----------------------------
# 4. 绘图
# -----------------------------

# 2. 画布 ---------------------------------------------------------------------


# 3. 主曲线 + 填充 ------------------------------------------------------------
ax.plot(x_init,y_init,color="#2ECC75", lw=2, alpha=0.3, label="XRD init")
ax.plot(x, y, color="#1E4AA0FF", lw=2, label="XRD corrected",alpha=0.7)
ax.plot(x, bkg, color="#7D8DCF", linestyle='--', lw=2, alpha=0.3, label="baseline")


ax.fill_between(x, y, color="#35B17FFF", alpha=0.12, zorder=3)

# 4. peak 标记 ----------------------------------------------------------------
# 002
ax.scatter(two_theta_002_fit, I_002_fit, color="#CD5F5FFF", s=100,
        edgecolor="white", linewidth=1.5, zorder=6)
ax.annotate(f"002 : {two_theta_002_fit:.2f}°",
            xy=(two_theta_002_fit+0.5, I_002_fit+0.5),
            xytext=(two_theta_002_fit + 2.2, I_002_fit * 0.95),
            color="#CD5F5FFF", fontweight="bold", fontsize=13,
            arrowprops=dict(arrowstyle="->", color="#CD5F5FFF", lw=1.5))

plt.plot(x[peaks_002[0]:peaks_002[1]], 
             p002fit_y, 
             linestyle='-.', label='002 Peak fit', 
             alpha=0.9,color='#04B7B1')


# 100
ax.scatter(two_theta_100_fit, I_100_fit, color="#6854B7FF", s=100,
        edgecolor="white", linewidth=1.5, zorder=6)
ax.annotate(f"100 : {two_theta_100_fit:.2f}°",
            xy=(two_theta_100_fit, I_100_fit+0.5),
            xytext=(two_theta_100_fit - 1.5, I_100_fit * 1.28),
            color="#6854B7FF", fontweight="bold", fontsize=13,
            arrowprops=dict(arrowstyle="->", color="#6854B7FF", lw=1.5))

plt.plot(x[peaks_100[0]:peaks_100[1]], 
             p100fit_y, 
             linestyle='-.', label='100 Peak fit', 
             alpha=0.9, color='#9C31C3')

# 5. 信息卡片 ------------------------------------------------------------------
outinfo_txt = (f"d002           {d002:.6f} Å\n"  # 单位 Å
               f"Lc             {Lc:.6f} Å\n"
               f"La             {La:.6f} Å\n"
               f"Peak_002       {two_theta_002_fit:.6f} °\n"
               f"Peak_100       {two_theta_100_fit:.6f} °\n"
               f"FWHM_002       {beta_002:.6f} °\n"
               f"FWHM_100       {beta_100:.6f} °\n"
               f"R2_002         {r2_002:.6f}\n"
               f"R2_100         {r2_100:.6f}"
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
low_y = -1.2

p002 = peaks_002
p100 = peaks_100
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

# 计算 bottom 位置
y_bottom = min(plt.ylim()[0],-3)

# 6. 坐标轴 & 标题 -------------------------------------------------------------
ax.set_xlabel(r"2$\mathbf{\theta}$ (deg)", fontsize=15,fontweight='bold')
ax.set_ylabel("Intensity (a.u.)", fontsize=15,fontweight='bold')
ax.set_title("XRD Pattern of Carbon Fiber", fontsize=17, fontweight='bold', pad=15)
ax.set_xlim(10, 60)
ax.set_ylim(bottom=y_bottom)

# 7. 细节 ----------------------------------------------------------------------
ax.tick_params(axis='both', labelsize=12)
fig.tight_layout()
ax.xaxis.grid(True, which='major', linestyle='-', linewidth=1.5, color='#0B5D56FF',alpha=0.25)
ax.yaxis.grid(False)          # 确保 y 轴网格不出现
ax.legend(fontsize=12, loc='upper left')
fig.savefig("xrd.png", dpi=400, bbox_inches='tight')


print(f"\n====================>分析完毕 -.O")