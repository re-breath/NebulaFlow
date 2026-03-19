import json
import numpy as np
import matplotlib.pyplot as plt
import argparse

# 从material project中获取声子谱以及声子DOS数据放到json中
# 使用方式：python plot.py mp-2741_phonon_bs_dfpt.json  将会从mp-2741_phonon_bs_dfpt.json中读取数据并绘制声子谱图

parser = argparse.ArgumentParser(description="Plot phonon bandstructure from Material Project JSON.")
parser.add_argument("bs_file", type=str, help="Path to the JSON file containing phonon bandstructure data.")
args = parser.parse_args()

BS_FILE = args.bs_file
OUT_PNG = BS_FILE.replace(".json", "_mpstyle3.png")



# BS_FILE = "mp-2741_phonon_bs_dfpt.json"
# OUT_PNG = "mp-2741_phonon_bs_dfpt_mpstyle3.png"
THZ_TO_CM1 = 33.35641

def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def get_reciprocal_matrix(data):
    rl = data["reciprocal_lattice"]
    if isinstance(rl, dict) and "matrix" in rl:
        return np.array(rl["matrix"], dtype=float)
    return np.array(rl, dtype=float)

def frac_to_cart(q_frac, rec_mat):
    q = np.array(q_frac, dtype=float)
    return q @ rec_mat

def cumulative_k_distance(qpoints_frac, rec_mat):
    k = frac_to_cart(qpoints_frac, rec_mat)   # (nq, 3)
    dk = np.linalg.norm(k[1:] - k[:-1], axis=1)
    x = np.concatenate([[0.0], np.cumsum(dk)])
    return x

def norm_label(s: str) -> str:
    return s.replace("\\Gamma", "Γ").replace("GAMMA", "Γ").replace("Gamma", "Γ")

def find_label_indices(qpoints, labels_dict, tol=1e-6):
    """
    labels_dict: {label: [qx,qy,qz]}  (你现在这种格式)
    返回: {index:int -> set(labels)}
    """
    q = np.array(qpoints, dtype=float)  # (nq,3)
    out = {}

    for lab, coord in labels_dict.items():
        lab = norm_label(str(lab))
        target = np.array(coord, dtype=float)

        # 找所有与 target 接近的点（允许数值误差）
        d = np.linalg.norm(q - target, axis=1)
        idxs = np.where(d < tol)[0]
        if idxs.size == 0:
            # 容忍一些：有时 tol 太严格
            idxs = np.where(d < 1e-4)[0]

        # 有些路径里同一个高对称点会出现多次（拼接），这里全保留
        for idx in idxs.tolist():
            out.setdefault(int(idx), set()).add(lab)

    return out

def build_ticks(x, index_to_labels):
    """
    把 index->labels 转成 (xticks, xlabels)
    - 同一位置多个标签合并 K|U
    - 若同一 label 在多个 index 出现（Γ 重复），保留全部（官网也是重复 Γ 的）
    """
    idxs = sorted(index_to_labels.keys())
    xticks = []
    xlabels = []

    for idx in idxs:
        labs = sorted(index_to_labels[idx])
        lab = "|".join(labs)
        xticks.append(float(x[idx]))
        xlabels.append(lab)

    # 合并“几乎同一位置”的 tick（浮点导致的重合）
    merged_ticks = []
    merged_labels = []
    eps = 1e-8
    for t, lab in sorted(zip(xticks, xlabels), key=lambda z: z[0]):
        if not merged_ticks or abs(t - merged_ticks[-1]) > eps:
            merged_ticks.append(t)
            merged_labels.append(lab)
        else:
            # 同一个位置合并标签
            merged_labels[-1] = "|".join(sorted(set(merged_labels[-1].split("|") + lab.split("|"))))

    return merged_ticks, merged_labels

def main():
    data = load_json(BS_FILE)

    qpoints = data["qpoints"]
    freqs = np.array(data["frequencies"], dtype=float)  # (nq, nbranch) or (nbranch, nq)

    # shape fix
    if freqs.shape[0] != len(qpoints) and freqs.shape[1] == len(qpoints):
        freqs = freqs.T
    if freqs.shape[0] != len(qpoints):
        raise ValueError(f"qpoints={len(qpoints)} but frequencies shape={freqs.shape}")

    # x axis in 1/Å
    rec_mat = get_reciprocal_matrix(data)
    x = cumulative_k_distance(qpoints, rec_mat)

    # y axis in cm^-1 (like MP plot)
    y = freqs * THZ_TO_CM1

    # plot branches
    for b in range(y.shape[1]):
        plt.plot(x, y[:, b], linewidth=1.2)

    # ticks & vertical lines
    labels_dict = data.get("labels_dict", {})
    index_to_labels = find_label_indices(qpoints, labels_dict, tol=1e-6)
    xticks, xlabels = build_ticks(x, index_to_labels)

    if xticks:
        for t in xticks:
            plt.axvline(t, linewidth=0.8)
        plt.xticks(xticks, [norm_label(s) for s in xlabels])

    plt.xlabel("Wave vector")
    plt.ylabel("Frequencies (cm$^{-1}$)")
    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    print("Saved:", OUT_PNG)

if __name__ == "__main__":
    main()