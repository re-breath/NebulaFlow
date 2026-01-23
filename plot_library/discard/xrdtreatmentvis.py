#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
分析 XRD.exe 输出的结构因子文件，计算：
  d002、Lc（沿 c 轴晶粒厚度）
  La（沿 a/b 轴微晶尺寸）
并绘制 (002) 与 (100) 峰拟合结果
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# ---------- 用户参数 ----------
lambda_A      = 1.5406      # Cu Kα 波长 (Å)
K_002         = 0.89        # Scherrer 常数 002
K_100         = 1.84        # Scherrer 常数 100
beta_inst_deg = 0.10        # 仪器 FWHM (°)
data_file     = 'SFACTOR_L1_L2_na.data'  # 输入文件
# --------------------------------

# ---- 读数据 ----
x, y = np.loadtxt(data_file, unpack=True)

# ---- 峰窗口 ----
win_002 = (22.5, 25.0)
win_100 = (40.0, 47.0)

# ---- 函数定义 ----
def baseline_poly(x, a0, a1):
    return a0 + a1 * x

def pseudo_voigt(x, A, x0, fwhm, eta):
    sigma   = fwhm / (2 * np.sqrt(2 * np.log(2)))
    gauss   = np.exp(-(x - x0)**2 / (2 * sigma**2))
    lorentz = 1.0 / (1.0 + ((x - x0) / (fwhm / 2.0))**2)
    return A * (eta * lorentz + (1 - eta) * gauss)

def fit_peak(x, y, xmin, xmax):
    mask = (x >= xmin) & (x <= xmax)
    xf, yf = x[mask], y[mask]

    # 1. 拟合线性背景
    pbl, _ = curve_fit(baseline_poly, xf, yf, p0=[yf.min(), 0.0])
    y_corr = yf - baseline_poly(xf, *pbl)

    # 2. 初值
    A0   = y_corr.max()
    x00  = xf[np.argmax(y_corr)]
    fwhm0 = (xmax - xmin) / 8.0
    eta0 = 0.5

    # 3. 拟合
    popt, _ = curve_fit(pseudo_voigt, xf, y_corr,
                        p0=[A0, x00, fwhm0, eta0], maxfev=20000)
    A, x0, fwhm, eta = popt
    return (x0, fwhm), (xf, yf), pbl, (A, eta)

# ---- 拟合 (002) ----
(x0_002, fwhm_002), (xf002, yf002), pbl002, (A_002, eta_002) = fit_peak(x, y, *win_002)
theta_002 = np.deg2rad(x0_002 / 2.0)
beta_002  = np.deg2rad(np.sqrt(max(fwhm_002**2 - beta_inst_deg**2, 0.0)))
d002 = lambda_A / (2.0 * np.sin(theta_002))
Lc   = K_002 * lambda_A / (beta_002 * np.cos(theta_002))

# ---- 拟合 (100) ----
(x0_100, fwhm_100), (xf100, yf100), pbl100, (A_100, eta_100) = fit_peak(x, y, *win_100)
theta_100 = np.deg2rad(x0_100 / 2.0)
beta_100  = np.deg2rad(np.sqrt(max(fwhm_100**2 - beta_inst_deg**2, 0.0)))
La = K_100 * lambda_A / (beta_100 * np.cos(theta_100))

# ---- 终端输出 ----
print(f'2θ(002) = {x0_002:6.3f}°   d002 = {d002:7.4f} Å   Lc = {Lc:6.2f} Å')
print(f'2θ(100) = {x0_100:6.3f}°   La   = {La:7.2f} Å')

# ---- 可视化 ----
def plot_peak(xf, yf, pbl, A, x0, fwhm, eta, title, savefile):
    plt.figure(figsize=(5, 3), tight_layout=True)
    baseline = baseline_poly(xf, *pbl)
    y_corr   = yf - baseline

    plt.plot(xf, yf, 'k.', label='raw')
    plt.plot(xf, baseline, 'C7--', label='baseline')
    plt.plot(xf, y_corr, 'C0o', mfc='none', label='bg-corrected')
    plt.plot(xf, pseudo_voigt(xf, A, x0, fwhm, eta), 'C1-', label='Pseudo-Voigt')
    plt.xlabel(r'2θ (°)')
    plt.ylabel('Intensity (arb. u.)')
    plt.title(title)
    plt.legend()
    plt.savefig(savefile, dpi=150)
    plt.show()

plot_peak(xf002, yf002, pbl002, A_002, x0_002, fwhm_002, eta_002,
          title=f'(002) peak  2θ={x0_002:.2f}°', savefile='fit_002.png')

plot_peak(xf100, yf100, pbl100, A_100, x0_100, fwhm_100, eta_100,
          title=f'(100) peak  2θ={x0_100:.2f}°', savefile='fit_100.png')
