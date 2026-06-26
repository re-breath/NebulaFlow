#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from collections import Counter

import numpy as np
from ase.io import read


# 1 amu / Å^3 = 1.66053906660 g/cm^3
AMU_PER_A3_TO_G_CM3 = 1.66053906660


def has_valid_3d_cell(atoms, tol=1e-12):
    """
    判断结构是否有完整三维晶胞。
    普通 xyz 文件通常没有晶格；extxyz/cif/vasp 通常有。
    """
    cell = np.asarray(atoms.get_cell())
    return cell.shape == (3, 3) and abs(np.linalg.det(cell)) > tol


def analyze_structure(filename, index=-1, fmt=None):
    """
    使用 ASE 读取并分析结构。

    Parameters
    ----------
    filename : str
        输入结构文件路径
    index : int or str
        读取第几帧。-1 表示最后一帧；0 表示第一帧。
    fmt : str or None
        ASE 格式名。一般不用填，ASE 会自动判断。
    """
    atoms = read(filename, index=index, format=fmt)

    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()
    pbc = atoms.get_pbc()

    result = {
        "filename": filename,
        "num_atoms": len(atoms),
        "elements": Counter(symbols),
        "pbc": pbc.tolist(),
        "has_cell": has_valid_3d_cell(atoms),
        "cell": None,
        "cell_lengths": None,
        "cell_angles": None,
        "volume_A3": None,
        "density_g_cm3": None,
        "coord_min": None,
        "coord_max": None,
        "coord_span": None,
    }

    if len(atoms) > 0:
        coord_min = positions.min(axis=0)
        coord_max = positions.max(axis=0)
        result["coord_min"] = coord_min
        result["coord_max"] = coord_max
        result["coord_span"] = coord_max - coord_min

    if result["has_cell"]:
        cell = atoms.get_cell()
        a, b, c, alpha, beta, gamma = cell.cellpar()
        volume = atoms.get_volume()

        total_mass_amu = atoms.get_masses().sum()
        density = total_mass_amu / volume * AMU_PER_A3_TO_G_CM3

        result["cell"] = np.asarray(cell)
        result["cell_lengths"] = np.array([a, b, c])
        result["cell_angles"] = np.array([alpha, beta, gamma])
        result["volume_A3"] = volume
        result["density_g_cm3"] = density

    return result


def print_report(result):
    print("=" * 60)
    print("Structure Analysis Report")
    print("=" * 60)

    print(f"\nInput file: {result['filename']}")
    print(f"Atoms: {result['num_atoms']}")

    print("\nElements:")
    for elem, count in result["elements"].items():
        ratio = count / result["num_atoms"] * 100
        print(f"  {elem:>3s}: {count:6d} atoms ({ratio:6.2f}%)")

    print("\nPBC:")
    for name, value in zip(["a", "b", "c"], result["pbc"]):
        print(f"  {name}-direction: {'True' if value else 'False'}")

    print("\nCell:")
    if result["has_cell"]:
        print("  Lattice matrix / Å:")
        for row in result["cell"]:
            print(f"    {row[0]:12.6f} {row[1]:12.6f} {row[2]:12.6f}")

        a, b, c = result["cell_lengths"]
        alpha, beta, gamma = result["cell_angles"]

        print("\n  Cell lengths / Å:")
        print(f"    a = {a:.6f}")
        print(f"    b = {b:.6f}")
        print(f"    c = {c:.6f}")

        print("\n  Cell angles / degree:")
        print(f"    alpha = {alpha:.6f}")
        print(f"    beta  = {beta:.6f}")
        print(f"    gamma = {gamma:.6f}")

        print(f"\n  Volume:  {result['volume_A3']:.6f} Å^3")
        print(f"  Density: {result['density_g_cm3']:.6f} g/cm^3")
    else:
        print("  No valid 3D cell found.")
        print("  Density and cell volume are not available.")

    print("\nCoordinate range / Å:")
    if result["coord_min"] is not None:
        for axis, mn, mx, span in zip(
            ["X", "Y", "Z"],
            result["coord_min"],
            result["coord_max"],
            result["coord_span"],
        ):
            print(f"  {axis}: {mn:12.6f} to {mx:12.6f}, span = {span:12.6f}")
    else:
        print("  No atoms found.")

    print("\n" + "=" * 60)


def save_summary(result, output_file):
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(f"file {result['filename']}\n")
        f.write(f"atoms {result['num_atoms']}\n")

        for elem, count in result["elements"].items():
            f.write(f"element_{elem} {count}\n")

        f.write(f"pbc {' '.join(str(x) for x in result['pbc'])}\n")

        if result["has_cell"]:
            a, b, c = result["cell_lengths"]
            alpha, beta, gamma = result["cell_angles"]

            f.write(f"lattice_a {a:.6f} Angstrom\n")
            f.write(f"lattice_b {b:.6f} Angstrom\n")
            f.write(f"lattice_c {c:.6f} Angstrom\n")
            f.write(f"alpha {alpha:.6f} degree\n")
            f.write(f"beta {beta:.6f} degree\n")
            f.write(f"gamma {gamma:.6f} degree\n")
            f.write(f"volume {result['volume_A3']:.6f} Angstrom^3\n")
            f.write(f"density {result['density_g_cm3']:.6f} g/cm^3\n")
        else:
            f.write("volume NA\n")
            f.write("density NA\n")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze atomic structure using ASE."
    )
    parser.add_argument("filename", help="Input structure file, e.g. xyz, extxyz, cif, POSCAR")
    parser.add_argument(
        "-i", "--index",
        default="-1",
        help="Frame index for multi-frame files. Default: -1, meaning the last frame.",
    )
    parser.add_argument(
        "-f", "--format",
        default=None,
        help="ASE file format. Usually unnecessary because ASE can auto-detect.",
    )
    parser.add_argument(
        "-o", "--output",
        default="structure_summary.txt",
        help="Output summary file.",
    )

    args = parser.parse_args()

    # ASE 的 index 可以是 int，也可以是 ':' 这类字符串
    try:
        index = int(args.index)
    except ValueError:
        index = args.index

    result = analyze_structure(args.filename, index=index, fmt=args.format)
    print_report(result)
    save_summary(result, args.output)

    print(f"\nSummary saved to: {args.output}")


if __name__ == "__main__":
    main()
