# 该脚本使用来分析XRD.exe输出的结果

import numpy as np
from scipy.optimize import curve_fit

# --- config ---
lambda_A = 1.5406       # Cu Kα (Å)
K_002 = 0.89
K_100 = 1.84
beta_inst_deg = 0.0     # instrument FWHM in degrees (replace with your value)

# --- load ---
x, y = np.loadtxt('SFACTOR_L1_L2_na.data', unpack=True)

# --- choose peak windows (adjust as needed) ---
win_002 = (22.5, 25.0)
win_100 = (40.0, 47.0)

def baseline_poly(x, a0, a1):
    return a0 + a1*x

def pseudo_voigt(x, A, x0, fwhm, eta):
    sigma = fwhm / (2*np.sqrt(2*np.log(2)))
    gauss = np.exp(-(x-x0)**2 / (2*sigma**2))
    lorentz = 1.0 / (1.0 + ((x-x0)/(fwhm/2.0))**2)
    return A*(eta*lorentz + (1-eta)*gauss)

def fit_peak(x, y, xmin, xmax):
    mask = (x >= xmin) & (x <= xmax)
    xf, yf = x[mask], y[mask]
    # fit baseline
    pbl, _ = curve_fit(baseline_poly, xf, yf, p0=[yf.min(), 0.0])
    y_corr = yf - baseline_poly(xf, *pbl)
    # initial guesses
    A0 = y_corr.max()
    x0 = xf[np.argmax(y_corr)]
    fwhm0 = (xmax - xmin)/8.0
    eta0 = 0.5
    popt, _ = curve_fit(pseudo_voigt, xf, y_corr, p0=[A0, x0, fwhm0, eta0], maxfev=20000)
    A, x0, fwhm, eta = popt
    return (x0, fwhm), (xf, yf), pbl, (A, eta)

# fit (002)
(x0_002, fwhm_002), _, pbl002, _ = fit_peak(x, y, *win_002)
theta_002 = np.deg2rad(x0_002/2.0)
beta_002 = np.deg2rad(np.sqrt(max(fwhm_002**2 - beta_inst_deg**2, 0.0)))

d002 = lambda_A / (2.0*np.sin(theta_002))
Lc = K_002 * lambda_A / (beta_002 * np.cos(theta_002))

# fit (100)
(x0_100, fwhm_100), _, pbl100, _ = fit_peak(x, y, *win_100)
theta_100 = np.deg2rad(x0_100/2.0)
beta_100 = np.deg2rad(np.sqrt(max(fwhm_100**2 - beta_inst_deg**2, 0.0)))
La = K_100 * lambda_A / (beta_100 * np.cos(theta_100))

print(f"2θ(002) = {x0_002:.3f} deg  -> d002 = {d002:.4f} Å,  Lc = {Lc:.2f} Å")
print(f"2θ(100) = {x0_100:.3f} deg  -> La   = {La:.2f} Å")
