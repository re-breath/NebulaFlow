import numpy as np
import matplotlib.pyplot as plt

# -------------------- config --------------------
VIRIAL_FILTER_THRESHOLD = -1000.0
VIRIAL_SCALE = 1.0   # 1.0 -> eV/atom; 1000.0 -> meV/atom
VIRIAL_UNIT = 'eV/atom' if VIRIAL_SCALE == 1.0 else 'meV/atom'

SCATTER_SIZE = 10
SCATTER_ALPHA = 0.35

LOSS_LABELS = [
    'Total',
    'L1-regularization',
    'L2-regularization',
    'Energy-train',
    'Force-train',
    'Energy-test',
    'Force-test',
]
# ------------------------------------------------


def load_data():
    return {
        'loss': np.loadtxt('loss.out'),
        'energy_train': np.loadtxt('energy_train.out'),
        'energy_test': np.loadtxt('energy_test.out'),
        'force_train': np.loadtxt('force_train.out'),
        'force_test': np.loadtxt('force_test.out'),
        'virial_train': np.loadtxt('virial_train.out'),
        'virial_test': np.loadtxt('virial_test.out'),
    }


def calc_metrics(y_true, y_pred):
    y_true = np.asarray(y_true).ravel()
    y_pred = np.asarray(y_pred).ravel()
    err = y_pred - y_true

    rmse = np.sqrt(np.mean(err**2))
    mae = np.mean(np.abs(err))
    bias = np.mean(err)
    maxae = np.max(np.abs(err))

    ss_res = np.sum(err**2)
    ss_tot = np.sum((y_true - np.mean(y_true))**2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan

    return {
        'rmse': rmse,
        'mae': mae,
        'bias': bias,
        'maxae': maxae,
        'r2': r2,
    }


def energy_xy(arr):
    # x: DFT, y: NEP
    return arr[:, 1], arr[:, 0]


def force_xy(arr):
    # x: DFT, y: NEP
    return arr[:, 3:6].reshape(-1), arr[:, 0:3].reshape(-1)


def clean_virial(arr, threshold=VIRIAL_FILTER_THRESHOLD):
    return arr[arr[:, 1] > threshold, :]


def virial_xy(arr, scale=VIRIAL_SCALE):
    # 如果你的 virial 文件实际有 6 个分量，把下面 0:5/6:11 改成 0:6/6:12
    return arr[:, 6:11].reshape(-1) * scale, arr[:, 0:5].reshape(-1) * scale


def format_label(name, metrics, unit):
    return (
        f"{name}\n"
        f"RMSE={metrics['rmse']:.5f} {unit}, "
        f"MAE={metrics['mae']:.5f} {unit}\n"
        f"R²={metrics['r2']:.6f}, "
        f"Bias={metrics['bias']:.5f} {unit}"
    )


def parity_panel(ax, x_train, y_train, x_test, y_test, unit, title):
    m_train = calc_metrics(x_train, y_train)
    m_test = calc_metrics(x_test, y_test)

    all_vals = np.concatenate([x_train, y_train, x_test, y_test])
    vmin, vmax = np.min(all_vals), np.max(all_vals)
    pad = 0.05 * (vmax - vmin) if vmax > vmin else 1.0
    lo, hi = vmin - pad, vmax + pad

    ax.scatter(
        x_train, y_train,
        s=SCATTER_SIZE, alpha=SCATTER_ALPHA,
        label=format_label('Train', m_train, unit)
    )
    ax.scatter(
        x_test, y_test,
        s=SCATTER_SIZE, alpha=SCATTER_ALPHA,
        label=format_label('Test', m_test, unit)
    )
    ax.plot([lo, hi], [lo, hi], '--', lw=1.2, label='y = x')

    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(alpha=0.25)
    ax.set_title(title)
    ax.legend(fontsize=8, frameon=True, loc='best')

    return m_train, m_test


def residual_panel(ax, x_train, y_train, x_test, y_test, unit, title):
    res_train = y_train - x_train
    res_test = y_test - x_test

    lo = min(res_train.min(), res_test.min())
    hi = max(res_train.max(), res_test.max())
    if hi == lo:
        lo -= 1.0
        hi += 1.0

    bins = np.linspace(lo, hi, 60)

    ax.hist(
        res_train, bins=bins, alpha=0.5, density=True,
        label=f"Train: mean={np.mean(res_train):.5f}, std={np.std(res_train):.5f}"
    )
    ax.hist(
        res_test, bins=bins, alpha=0.5, density=True,
        label=f"Test: mean={np.mean(res_test):.5f}, std={np.std(res_test):.5f}"
    )

    ax.axvline(0.0, ls='--', lw=1.2)
    ax.grid(alpha=0.25)
    ax.set_xlabel(f"Residual (NEP - DFT) [{unit}]")
    ax.set_ylabel('Density')
    ax.set_title(title)
    ax.legend(fontsize=8, frameon=True)


def summary_panel(ax, metrics_dict):
    ax.axis('off')

    def line(name, mtr, mte, unit):
        gap = mte['rmse'] - mtr['rmse']
        ratio = mte['rmse'] / mtr['rmse'] if mtr['rmse'] > 0 else np.nan
        return (
            f"{name}\n"
            f"  Train: RMSE={mtr['rmse']:.5f}, MAE={mtr['mae']:.5f}, R²={mtr['r2']:.6f}\n"
            f"  Test : RMSE={mte['rmse']:.5f}, MAE={mte['mae']:.5f}, R²={mte['r2']:.6f}\n"
            f"  Gap(test-train RMSE)={gap:.5f} {unit}, Ratio={ratio:.3f}\n"
        )

    text = (
        "Model summary\n\n"
        + line('Energy', metrics_dict['energy_train'], metrics_dict['energy_test'], 'eV/atom')
        + line('Force', metrics_dict['force_train'], metrics_dict['force_test'], 'eV/A')
        + line('Virial', metrics_dict['virial_train'], metrics_dict['virial_test'], VIRIAL_UNIT)
        + "\nJudgement tips\n"
        + "1. Test RMSE close to Train RMSE -> better generalization\n"
        + "2. Residual histogram centered at 0 and narrow -> smaller systematic error\n"
        + "3. R² closer to 1 -> parity closer to ideal line\n"
        + "4. If test loss rises while train loss keeps falling -> overfitting\n"
    )

    ax.text(
        0.0, 1.0, text,
        va='top', ha='left',
        fontsize=10, family='monospace'
    )


def main():
    data = load_data()

    data['virial_train'] = clean_virial(data['virial_train'])
    data['virial_test'] = clean_virial(data['virial_test'])

    # x: DFT, y: NEP
    e_tr_x, e_tr_y = energy_xy(data['energy_train'])
    e_te_x, e_te_y = energy_xy(data['energy_test'])

    f_tr_x, f_tr_y = force_xy(data['force_train'])
    f_te_x, f_te_y = force_xy(data['force_test'])

    v_tr_x, v_tr_y = virial_xy(data['virial_train'])
    v_te_x, v_te_y = virial_xy(data['virial_test'])

    metrics = {
        'energy_train': calc_metrics(e_tr_x, e_tr_y),
        'energy_test': calc_metrics(e_te_x, e_te_y),
        'force_train': calc_metrics(f_tr_x, f_tr_y),
        'force_test': calc_metrics(f_te_x, f_te_y),
        'virial_train': calc_metrics(v_tr_x, v_tr_y),
        'virial_test': calc_metrics(v_te_x, v_te_y),
    }

    fig, axes = plt.subplots(2, 3, figsize=(18, 10), constrained_layout=True)

    # (1) loss
    ax = axes[0, 0]
    loss = data['loss']
    ax.loglog(loss[:, 1:6])
    ax.loglog(loss[:, 7:9])
    ax.set_xlabel('Generation / 100')
    ax.set_ylabel('Loss')
    ax.set_title('Loss curves')
    ax.grid(alpha=0.25, which='both')
    ax.legend(LOSS_LABELS, fontsize=8, frameon=True)

    # (2) energy parity
    ax = axes[0, 1]
    parity_panel(ax, e_tr_x, e_tr_y, e_te_x, e_te_y, 'eV/atom', 'Energy parity')
    ax.set_xlabel('DFT energy (eV/atom)')
    ax.set_ylabel('NEP energy (eV/atom)')

    # (3) force parity
    ax = axes[0, 2]
    parity_panel(ax, f_tr_x, f_tr_y, f_te_x, f_te_y, 'eV/A', 'Force parity')
    ax.set_xlabel('DFT force (eV/A)')
    ax.set_ylabel('NEP force (eV/A)')

    # (4) virial parity
    ax = axes[1, 0]
    parity_panel(ax, v_tr_x, v_tr_y, v_te_x, v_te_y, VIRIAL_UNIT, 'Virial parity')
    ax.set_xlabel(f'DFT virial ({VIRIAL_UNIT})')
    ax.set_ylabel(f'NEP virial ({VIRIAL_UNIT})')

    # (5) residual distribution
    ax = axes[1, 1]
    residual_panel(ax, f_tr_x, f_tr_y, f_te_x, f_te_y, 'eV/A', 'Force residual distribution')

    # (6) summary
    ax = axes[1, 2]
    summary_panel(ax, metrics)

    plt.savefig('nep_evaluation.png', dpi=300)
    plt.show()


if __name__ == '__main__':
    main()

