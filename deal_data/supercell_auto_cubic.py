import math
import numpy as np
from ase.io import read, write


def is_orthorhombic(cell, tol=1e-6):
    """判断晶胞是否正交（允许极小数值误差）。"""
    cell = np.array(cell)
    # 正交：非对角项接近 0
    off = cell.copy()
    np.fill_diagonal(off, 0.0)
    return np.all(np.abs(off) < tol)


def cube_score(Lx, Ly, Lz):
    """
    越接近正方体越好：使用 max/min 作为指标（=1 最好）。
    """
    mx = max(Lx, Ly, Lz)
    mn = min(Lx, Ly, Lz)
    return mx / mn


def find_best_repeat(a, b, c, n0, target_natoms,
                     max_factor=30,
                     w_count=1.0,
                     w_cube=50.0):
    """
    搜索整数 (nx, ny, nz)，使:
    - 原子数接近 target_natoms
    - (a*nx, b*ny, c*nz) 接近正方体
    用一个综合损失函数来选最优。
    """
    # 目标总倍数（连续值）
    target_mult = target_natoms / n0
    if target_mult <= 0:
        raise ValueError("target_natoms 必须为正数")

    # 期望“立方体边长”（连续近似）
    # 令 a*nx ~ b*ny ~ c*nz ~ L
    # 且 nx*ny*nz ~ target_mult
    # => (L/a)*(L/b)*(L/c) = L^3/(abc) ~ target_mult
    L = (target_mult * a * b * c) ** (1/3)

    # 连续近似的 nx,ny,nz
    nx0, ny0, nz0 = L / a, L / b, L / c

    # 搜索范围：围绕近似值，限制在 [1, max_factor]
    # 同时给个冗余窗口，避免错过更优整数组合
    def window(x):
        lo = max(1, int(math.floor(x)) - 6)
        hi = min(max_factor, int(math.ceil(x)) + 6)
        return range(lo, hi + 1)

    rx, ry, rz = window(nx0), window(ny0), window(nz0)

    best = None  # (loss, nx, ny, nz, N, cube_ratio)

    for nx in rx:
        for ny in ry:
            for nz in rz:
                mult = nx * ny * nz
                N = n0 * mult

                # 1) 原子数误差（相对误差更稳）
                count_err = abs(N - target_natoms) / target_natoms

                # 2) 立方体程度
                Lx, Ly, Lz = a * nx, b * ny, c * nz
                cube_ratio = cube_score(Lx, Ly, Lz)  # 1 最好
                cube_err = cube_ratio - 1.0

                # 综合损失：你可以调 w_count / w_cube
                loss = w_count * count_err + w_cube * cube_err

                if best is None or loss < best[0]:
                    best = (loss, nx, ny, nz, N, cube_ratio)

    if best is None:
        raise RuntimeError("没有找到可行的扩胞系数（请增大 max_factor）")

    return best  # loss, nx, ny, nz, N, cube_ratio


def supercell_to_near_cube(infile, target_natoms,
                           fmt_out="extxyz",
                           max_factor=30):
    atoms = read(infile)

    # xyz 若没有 cell/pbc，你需要先确保文件是 extxyz 或你手动设置过 cell
    cell = atoms.cell.array
    if np.allclose(cell, 0):
        raise ValueError("读到的结构没有晶胞信息(cell)。请使用 extxyz，或先给 atoms 设置 cell/pbc。")

    if not is_orthorhombic(cell):
        raise ValueError("晶胞不是正交的（orthorhombic）。如果是一般晶胞，需要用不同策略（如 Niggli/形状优化）。")

    # 正交晶胞边长
    a, b, c = atoms.cell.lengths()
    n0 = len(atoms)

    loss, nx, ny, nz, N, cube_ratio = find_best_repeat(
        a, b, c, n0, target_natoms,
        max_factor=max_factor
    )

    # 扩胞
    sc = atoms.repeat((nx, ny, nz))
    outfile = infile.split(".")[0] + f"_{N}.xyz"

    # 输出
    write(outfile, sc, format=fmt_out)

    print(f"Input atoms: {n0}")
    print(f"Target atoms: {target_natoms}")
    print(f"Chosen repeat: (nx, ny, nz)=({nx}, {ny}, {nz}) -> atoms={N}")
    print(f"New lengths: {a*nx:.4f}, {b*ny:.4f}, {c*nz:.4f} Å")
    print(f"Cube ratio (max/min): {cube_ratio:.6f}  (1.0 is perfect)")
    print(f"Wrote: {outfile}")


if __name__ == "__main__":
    import sys
    infile = sys.argv[1]
    target_natoms = int(sys.argv[2])
    supercell_to_near_cube(infile, target_natoms)
    pass
