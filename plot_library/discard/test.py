import re
import numpy as np
import matplotlib.pyplot as plt

# =========================
# 1. Read processed kappa file
# =========================
def read_kappa_file(filename):
    header_lines = []
    data_rows = []

    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                header_lines.append(line)
            else:
                parts = line.split()
                if len(parts) >= 2:
                    data_rows.append([float(parts[0]), float(parts[1])])

    data = np.array(data_rows, dtype=float)
    time_ns = data[:, 0]
    kappa = data[:, 1]

    info = {
        "label": None,
        "plateau_mean": None,
        "plateau_std": None,
        "plateau_start": None,
        "plateau_end": None
    }

    header_text = " ".join(header_lines)

    m_label = re.search(r'kappa_([xyz])_running_average', header_text, re.IGNORECASE)
    if m_label:
        info["label"] = m_label.group(1).lower()

    m_mean = re.search(r'plateau_mean=([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)', header_text)
    if m_mean:
        info["plateau_mean"] = float(m_mean.group(1))

    m_std = re.search(r'plateau_std=([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)', header_text)
    if m_std:
        info["plateau_std"] = float(m_std.group(1))

    m_range = re.search(
        r'plateau_range_ns=\[\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s*,\s*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)\s*\]',
        header_text
    )
    if m_range:
        info["plateau_start"] = float(m_range.group(1))
        info["plateau_end"] = float(m_range.group(2))

    return time_ns, kappa, info


# =========================
# 2. Load x/y/z data
# =========================
files = {
    "x": "kappa_dealed_x.txt",
    "y": "kappa_dealed_y.txt",
    "z": "kappa_dealed_z.txt"
}

data_dict = {}
for d, fname in files.items():
    t, k, info = read_kappa_file(fname)
    data_dict[d] = {"time": t, "kappa": k, "info": info}


# =========================
# 3. Plot style
# =========================
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 13,
    "axes.labelsize": 13,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "legend.fontsize": 12,
    "axes.linewidth": 1.2,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 5,
    "ytick.major.size": 5,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
    "savefig.dpi": 600
})

colors = {
    "x": "#1F3A5F",   # deep blue
    "y": "#8A5A36",   # warm brown
    "z": "#566F64"    # muted green
}


# =========================
# 4. Create figure with panels (a) and (b)
# =========================
fig, axes = plt.subplots(
    1, 2, figsize=(12.5, 5.6),
    gridspec_kw={"width_ratios": [2.3, 1.0]}
)

ax1, ax2 = axes

# ---------- Panel (a): running average curves ----------
for d in ["x", "y", "z"]:
    t = data_dict[d]["time"]
    k = data_dict[d]["kappa"]
    info = data_dict[d]["info"]

    ax1.plot(t, k, lw=2.2, color=colors[d], label=f"{d.upper()} direction")

    if info["plateau_mean"] is not None:
        ax1.axhline(info["plateau_mean"], ls="--", lw=1.5, color=colors[d], alpha=0.6)

ax1.set_xlabel("Time (ns)")
ax1.set_ylabel("Thermal conductivity (W/m/K)")
ax1.legend(frameon=False, loc="upper right")
ax1.grid(True, linestyle="--", linewidth=0.7, alpha=0.22)
ax1.set_title("(a)", fontsize=15,fontweight="bold",loc="left")


# ---------- Panel (b): plateau statistics ----------
dirs = ["x", "y", "z"]
xpos = np.arange(len(dirs))
means = [data_dict[d]["info"]["plateau_mean"] for d in dirs]
stds = [data_dict[d]["info"]["plateau_std"] for d in dirs]

for i, d in enumerate(dirs):
    ax2.errorbar(
        xpos[i], means[i], yerr=stds[i],
        fmt='o', markersize=8,
        capsize=5, elinewidth=1.6,
        color=colors[d], markerfacecolor=colors[d],
        markeredgecolor="black", markeredgewidth=0.7
    )
    # ax2.text(
    #     xpos[i]+0.1, means[i] + stds[i] + 0.1,
    #     f"{means[i]:.2f}±{stds[i]:.2f}",
    #     ha="center", va="bottom", fontsize=11
    # )

ax2.set_xticks(xpos)
ax2.set_xticklabels(["X", "Y", "Z"])
ax2.set_ylabel("Plateau thermal conductivity (W/m/K)")
ax2.set_xlim(-0.5, 2.5)
ax2.grid(True, axis="y", linestyle="--", linewidth=0.7, alpha=0.22)

ax2.set_title("(b)", fontsize=15,fontweight="bold",loc="left")

# spine style
for ax in axes:
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)

plt.tight_layout()
plt.savefig("thermal_conductivity_ab.png", bbox_inches="tight")
plt.savefig("thermal_conductivity_ab.pdf", bbox_inches="tight")
plt.show()
