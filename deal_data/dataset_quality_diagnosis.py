#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dataset_scientific_diagnosis.py

A reference-style pre-training dataset diagnosis tool for atomistic ML-potential datasets.

The script intentionally avoids making a single hard "good/bad" judgment. Instead it
separates diagnostics into dimensions:

  1. composition diversity and whether element-reference fitting is applicable
  2. geometry sanity: distances, volume, cell/PBC, atom counts
  3. label sanity: energies, forces, robust tails, possible unit/format issues
  4. optional descriptor-space coverage: PCA, nearest-neighbor isolation, clusters
  5. optional train/test or selected/unselected overlap checks

Typical use:

  python dataset_scientific_diagnosis.py train.xyz
  python dataset_scientific_diagnosis.py train.xyz --descriptor descriptor.out
  python dataset_scientific_diagnosis.py train.xyz --descriptor descriptor.out --outdir diagnosis_report

Optional comparison:

  python dataset_scientific_diagnosis.py train.xyz --test-xyz test.xyz
  python dataset_scientific_diagnosis.py train.xyz --descriptor descriptor.out --test-descriptor descriptor_test.out --test-xyz test.xyz

Dependencies:
  numpy, matplotlib, ase, scikit-learn
Optional:
  scipy (only used for Spearman correlation if available)
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import ase.io
    from ase.data import atomic_numbers, covalent_radii
except Exception as exc:  # pragma: no cover
    raise SystemExit("This script requires ASE. Install with: pip install ase") from exc

try:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from sklearn.neighbors import NearestNeighbors
    from sklearn.cluster import KMeans
except Exception as exc:  # pragma: no cover
    raise SystemExit("This script requires scikit-learn. Install with: pip install scikit-learn") from exc

try:
    from scipy.stats import spearmanr
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False


ENERGY_KEYS = [
    "energy", "Energy", "ENERGY", "E", "dft_energy", "DFT_energy",
    "free_energy", "Free_energy", "REF_energy", "ref_energy", "TotEnergy",
]
FORCE_KEYS = ["forces", "force", "Forces", "FORCES", "REF_forces", "ref_forces"]
STRESS_KEYS = ["stress", "virial", "virials", "Stress", "Virial", "REF_virial", "ref_virial"]


# -----------------------------------------------------------------------------
# Small utilities
# -----------------------------------------------------------------------------

def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def first_float(value: Any) -> float:
    """Parse the first floating number from a value, allowing strings like '-1.23eV'."""
    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value)
    s = str(value).strip()
    m = re.search(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", s)
    if m is None:
        raise ValueError(f"cannot parse float from {value!r}")
    return float(m.group(0))


def robust_mad_z(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    if not np.isfinite(mad) or mad < 1e-15:
        std = np.nanstd(x)
        if std < 1e-15:
            return np.zeros_like(x)
        return (x - np.nanmean(x)) / std
    return 0.67448975 * (x - med) / mad


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return float("nan")
    x = x[mask]
    y = y[mask]
    if np.std(x) < 1e-15 or np.std(y) < 1e-15:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    if not _HAS_SCIPY:
        return float("nan")
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return float("nan")
    if np.std(x[mask]) < 1e-15 or np.std(y[mask]) < 1e-15:
        return float("nan")
    r, _ = spearmanr(x[mask], y[mask])
    return float(r)


def percentile_dict(x: np.ndarray, prefix: str = "") -> Dict[str, float]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return {}
    ps = [0, 0.1, 1, 5, 10, 50, 90, 95, 99, 99.9, 100]
    return {f"{prefix}p{str(p).replace('.', '_')}": float(np.percentile(x, p)) for p in ps}


def fmt_float(x: Any, nd: int = 6) -> str:
    try:
        x = float(x)
    except Exception:
        return "NA"
    if not np.isfinite(x):
        return "NA"
    return f"{x:.{nd}f}"


def safe_array_loadtxt(path: str | Path) -> np.ndarray:
    arr = np.loadtxt(path)
    if arr.ndim == 1:
        arr = arr.reshape(1, -1)
    return np.asarray(arr, dtype=float)


# -----------------------------------------------------------------------------
# Reading XYZ labels and structure data
# -----------------------------------------------------------------------------

def read_xyz(path: str | Path):
    try:
        atoms_list = ase.io.read(str(path), index=":")
    except Exception as exc:
        raise SystemExit(f"Failed to read {path}: {exc}") from exc
    if not isinstance(atoms_list, list):
        atoms_list = [atoms_list]
    if len(atoms_list) == 0:
        raise SystemExit(f"No configurations found in {path}")
    return atoms_list


def read_energies(atoms_list) -> Tuple[np.ndarray, str]:
    """Read total energies from atoms.info."""
    for key in ENERGY_KEYS:
        vals = []
        ok = True
        for atoms in atoms_list:
            if key not in atoms.info:
                ok = False
                break
            try:
                vals.append(first_float(atoms.info[key]))
            except Exception:
                ok = False
                break
        if ok and len(vals) == len(atoms_list):
            return np.asarray(vals, dtype=float), f"atoms.info['{key}']"

    # Fallback: try ASE calculator single-point energy if available.
    vals = []
    ok = True
    for atoms in atoms_list:
        try:
            vals.append(float(atoms.get_potential_energy()))
        except Exception:
            ok = False
            break
    if ok and len(vals) == len(atoms_list):
        return np.asarray(vals, dtype=float), "atoms.get_potential_energy()"

    raise SystemExit(
        "Could not read energies. Expected one of these atoms.info keys: "
        + ", ".join(ENERGY_KEYS)
    )


def read_forces(atoms_list) -> Tuple[Optional[List[np.ndarray]], str]:
    for key in FORCE_KEYS:
        vals = []
        ok = True
        for atoms in atoms_list:
            if key not in atoms.arrays:
                ok = False
                break
            arr = np.asarray(atoms.arrays[key], dtype=float)
            if arr.shape != (len(atoms), 3):
                ok = False
                break
            vals.append(arr)
        if ok and len(vals) == len(atoms_list):
            return vals, f"atoms.arrays['{key}']"

    vals = []
    ok = True
    for atoms in atoms_list:
        try:
            arr = np.asarray(atoms.get_forces(), dtype=float)
            if arr.shape != (len(atoms), 3):
                ok = False
                break
            vals.append(arr)
        except Exception:
            ok = False
            break
    if ok and len(vals) == len(atoms_list):
        return vals, "atoms.get_forces()"
    return None, "not found"


def cell_volume_per_atom(atoms) -> float:
    try:
        vol = float(atoms.get_volume())
    except Exception:
        return float("nan")
    if vol <= 0 or len(atoms) == 0:
        return float("nan")
    return vol / len(atoms)


# -----------------------------------------------------------------------------
# Composition analysis
# -----------------------------------------------------------------------------

def composition_matrix(atoms_list) -> Tuple[List[str], np.ndarray, np.ndarray]:
    elements = sorted({s for a in atoms_list for s in a.get_chemical_symbols()})
    counts = np.array(
        [[a.get_chemical_symbols().count(e) for e in elements] for a in atoms_list],
        dtype=float,
    )
    n_atoms = counts.sum(axis=1)
    comp = counts / np.maximum(n_atoms[:, None], 1.0)
    return elements, counts, comp


def chemical_system(symbols: Iterable[str]) -> str:
    return "-".join(sorted(set(symbols)))


def reduced_formula_from_counts(elements: List[str], counts: np.ndarray) -> str:
    cnt = np.asarray(counts, dtype=int)
    nz = cnt[cnt > 0]
    if len(nz) == 0:
        return "EMPTY"
    g = int(nz[0])
    for v in nz[1:]:
        g = math.gcd(g, int(v))
    parts = []
    for e, c in zip(elements, cnt):
        if c <= 0:
            continue
        r = c // max(g, 1)
        parts.append(e if r == 1 else f"{e}{r}")
    return "".join(parts)


def analyze_composition(elements: List[str], counts: np.ndarray, comp: np.ndarray) -> Dict[str, Any]:
    n_cfg, n_elem = counts.shape
    rank_counts = int(np.linalg.matrix_rank(counts)) if n_cfg > 0 else 0
    rank_comp = int(np.linalg.matrix_rank(comp)) if n_cfg > 0 else 0

    # Condition number only meaningful if full column rank and enough variation.
    if rank_counts < n_elem:
        cond_counts = float("inf")
    else:
        try:
            cond_counts = float(np.linalg.cond(counts))
        except Exception:
            cond_counts = float("inf")

    comp_std = comp.std(axis=0)
    comp_range = comp.max(axis=0) - comp.min(axis=0)
    mean_comp_range = float(np.mean(comp_range)) if n_elem else 0.0
    max_comp_range = float(np.max(comp_range)) if n_elem else 0.0

    # A conservative type classifier. This does not judge quality; it tells which
    # diagnostics are applicable.
    if n_elem == 1:
        composition_diversity = "single-element"
        dataset_type = "single-composition / structure-dominated"
    elif max_comp_range < 1e-6 or rank_comp <= 1:
        composition_diversity = "very low"
        dataset_type = "fixed-stoichiometry / structure-dominated"
    elif max_comp_range < 0.05 or rank_comp < min(n_elem, 3):
        composition_diversity = "low"
        dataset_type = "near-fixed-stoichiometry / mostly structure-dominated"
    elif max_comp_range < 0.25:
        composition_diversity = "medium"
        dataset_type = "moderate-composition-diversity"
    else:
        composition_diversity = "high"
        dataset_type = "composition-rich / composition-and-structure mixed"

    formulas = [reduced_formula_from_counts(elements, row) for row in counts]
    unique_formulas, formula_counts = np.unique(formulas, return_counts=True)
    top_formula_idx = np.argsort(formula_counts)[::-1][:10]

    return {
        "n_elements": int(n_elem),
        "elements": elements,
        "rank_counts": rank_counts,
        "rank_comp": rank_comp,
        "cond_counts": cond_counts,
        "comp_std": {e: float(s) for e, s in zip(elements, comp_std)},
        "comp_range": {e: float(r) for e, r in zip(elements, comp_range)},
        "mean_comp_range": mean_comp_range,
        "max_comp_range": max_comp_range,
        "composition_diversity": composition_diversity,
        "dataset_type": dataset_type,
        "top_formulas": [
            {"formula": str(unique_formulas[i]), "count": int(formula_counts[i])}
            for i in top_formula_idx
        ],
    }


# -----------------------------------------------------------------------------
# Atomic reference baseline: now informational, not a global quality judge
# -----------------------------------------------------------------------------

def fit_atomic_reference(counts: np.ndarray, energies: np.ndarray, mode: str = "wls") -> Dict[str, Any]:
    n_atoms = counts.sum(axis=1)
    valid = (n_atoms > 0) & np.isfinite(energies)
    N = counts[valid]
    E = energies[valid]
    Nat = n_atoms[valid]
    if mode == "wls":
        w = 1.0 / np.maximum(Nat, 1.0)
        A = N * w[:, None]
        b = E * w
    elif mode == "ols":
        A, b = N, E
    else:
        raise ValueError("mode must be 'wls' or 'ols'")

    mu, *_ = np.linalg.lstsq(A, b, rcond=None)
    E_ref = counts @ mu
    raw_pa = energies / np.maximum(n_atoms, 1.0)
    ref_pa = E_ref / np.maximum(n_atoms, 1.0)
    shifted = raw_pa - ref_pa

    ss_res_total = float(np.sum((energies - E_ref) ** 2))
    ss_tot_total = float(np.sum((energies - np.mean(energies)) ** 2))
    r2_total = 1.0 - ss_res_total / ss_tot_total if ss_tot_total > 1e-20 else float("nan")

    ss_res_pa = float(np.sum((raw_pa - ref_pa) ** 2))
    ss_tot_pa = float(np.sum((raw_pa - np.mean(raw_pa)) ** 2))
    r2_pa = 1.0 - ss_res_pa / ss_tot_pa if ss_tot_pa > 1e-20 else float("nan")

    return {
        "mode": mode,
        "mu": mu,
        "E_ref": E_ref,
        "raw_per_atom": raw_pa,
        "ref_per_atom": ref_pa,
        "shifted_per_atom": shifted,
        "r2_total": r2_total,
        "r2_per_atom": r2_pa,
        "rmse_shifted_eV_atom": float(np.sqrt(np.mean(shifted ** 2))),
        "mae_shifted_eV_atom": float(np.mean(np.abs(shifted))),
        "range_shifted_eV_atom": [float(np.min(shifted)), float(np.max(shifted))],
    }


# -----------------------------------------------------------------------------
# Geometry and labels
# -----------------------------------------------------------------------------

def min_pair_distance_metrics(atoms, max_atoms_exact: int = 1500) -> Tuple[float, float, str]:
    """
    Return (min_distance_A, min_covalent_ratio, pair_label).
    ratio = d_ij / (r_cov_i + r_cov_j). Lower values indicate suspicious overlap.
    """
    n = len(atoms)
    if n < 2:
        return float("nan"), float("nan"), "NA"
    if n > max_atoms_exact:
        # exact all-distance matrix can become expensive; still most MLP datasets are smaller.
        # Return NaN rather than giving a misleading approximate answer.
        return float("nan"), float("nan"), "too-many-atoms"
    try:
        D = atoms.get_all_distances(mic=True)
    except Exception:
        D = atoms.get_all_distances(mic=False)
    D = np.asarray(D, dtype=float)
    np.fill_diagonal(D, np.inf)
    idx = np.unravel_index(np.argmin(D), D.shape)
    dmin = float(D[idx])
    syms = atoms.get_chemical_symbols()
    s1, s2 = syms[idx[0]], syms[idx[1]]
    try:
        rsum = float(covalent_radii[atomic_numbers[s1]] + covalent_radii[atomic_numbers[s2]])
        ratio = dmin / rsum if rsum > 0 else float("nan")
    except Exception:
        ratio = float("nan")
    return dmin, ratio, f"{s1}-{s2}"


def analyze_geometry(atoms_list) -> Dict[str, Any]:
    n_atoms = np.array([len(a) for a in atoms_list], dtype=float)
    volumes_pa = np.array([cell_volume_per_atom(a) for a in atoms_list], dtype=float)
    pbc_patterns = [tuple(bool(x) for x in a.get_pbc()) for a in atoms_list]
    pbc_counts: Dict[str, int] = {}
    for p in pbc_patterns:
        key = "".join("1" if x else "0" for x in p)
        pbc_counts[key] = pbc_counts.get(key, 0) + 1

    min_d = []
    min_ratio = []
    min_pair = []
    for atoms in atoms_list:
        d, r, pair = min_pair_distance_metrics(atoms)
        min_d.append(d)
        min_ratio.append(r)
        min_pair.append(pair)
    min_d = np.asarray(min_d, dtype=float)
    min_ratio = np.asarray(min_ratio, dtype=float)

    severe_overlap = np.where(min_ratio < 0.55)[0].tolist()
    warn_overlap = np.where((min_ratio >= 0.55) & (min_ratio < 0.70))[0].tolist()

    return {
        "atom_count": percentile_dict(n_atoms, "N_"),
        "atom_count_unique": int(len(np.unique(n_atoms))),
        "atom_count_values_top": [int(x) for x in np.unique(n_atoms)[:20]],
        "volume_per_atom": percentile_dict(volumes_pa, "volpa_"),
        "pbc_counts": pbc_counts,
        "min_distance": percentile_dict(min_d, "dmin_"),
        "min_covalent_ratio": percentile_dict(min_ratio, "ratio_"),
        "severe_overlap_indices": severe_overlap[:200],
        "warn_overlap_indices": warn_overlap[:200],
        "severe_overlap_count": int(len(severe_overlap)),
        "warn_overlap_count": int(len(warn_overlap)),
        "min_pair_labels": min_pair,
        "min_distance_array": min_d,
        "min_ratio_array": min_ratio,
        "volume_per_atom_array": volumes_pa,
        "atom_count_array": n_atoms,
    }


def analyze_labels(energies: np.ndarray, n_atoms: np.ndarray, forces: Optional[List[np.ndarray]]) -> Dict[str, Any]:
    E = np.asarray(energies, dtype=float)
    N = np.asarray(n_atoms, dtype=float)
    Epa = E / np.maximum(N, 1.0)
    label = {
        "energy_total": percentile_dict(E, "E_"),
        "energy_per_atom": percentile_dict(Epa, "Epa_"),
        "energy_per_atom_robust_outlier_fraction_z5": float(np.mean(np.abs(robust_mad_z(Epa)) > 5.0)),
        "E_total_vs_N_pearson": pearson(E, N),
        "Epa_vs_N_pearson": pearson(Epa, N),
    }
    if forces is not None:
        all_f = np.concatenate([np.linalg.norm(f, axis=1) for f in forces])
        cfg_mean = np.array([np.mean(np.linalg.norm(f, axis=1)) for f in forces], dtype=float)
        cfg_max = np.array([np.max(np.linalg.norm(f, axis=1)) for f in forces], dtype=float)
        label.update({
            "forces_found": True,
            "force_magnitude": percentile_dict(all_f, "F_"),
            "force_cfg_mean": percentile_dict(cfg_mean, "Fmean_cfg_"),
            "force_cfg_max": percentile_dict(cfg_max, "Fmax_cfg_"),
            "force_atom_fraction_gt_10": float(np.mean(all_f > 10.0)),
            "force_atom_fraction_gt_20": float(np.mean(all_f > 20.0)),
            "force_cfg_fraction_max_gt_10": float(np.mean(cfg_max > 10.0)),
            "force_cfg_fraction_max_gt_20": float(np.mean(cfg_max > 20.0)),
            "force_cfg_max_array": cfg_max,
            "force_cfg_mean_array": cfg_mean,
        })
    else:
        label["forces_found"] = False
    return label


# -----------------------------------------------------------------------------
# Descriptor analysis
# -----------------------------------------------------------------------------

def analyze_descriptor(descriptor: Optional[np.ndarray], n_cfg: int, shifted: Optional[np.ndarray], n_atoms: np.ndarray) -> Dict[str, Any]:
    if descriptor is None:
        return {"available": False, "note": "descriptor.out not provided"}
    X = np.asarray(descriptor, dtype=float)
    if X.ndim != 2:
        return {"available": False, "note": "descriptor is not 2D"}
    if X.shape[0] != n_cfg:
        return {"available": False, "note": f"descriptor rows {X.shape[0]} != number of configurations {n_cfg}"}
    Xs = StandardScaler().fit_transform(X)
    n_comp = min(5, Xs.shape[0], Xs.shape[1])
    pca = PCA(n_components=n_comp)
    Z = pca.fit_transform(Xs)

    # nearest neighbor distances in standardized descriptor space
    if len(Xs) >= 3:
        nn = NearestNeighbors(n_neighbors=2, metric="euclidean").fit(Xs)
        dists, _ = nn.kneighbors(Xs)
        nn1 = dists[:, 1]
    else:
        nn1 = np.full(len(Xs), np.nan)
    iso_z = robust_mad_z(nn1)
    iso_idx = np.where(iso_z > 5.0)[0].tolist()

    # modest cluster summary; not a judgment, just reveals mixture.
    k = min(8, max(2, int(round(np.sqrt(len(Xs) / 50.0))))) if len(Xs) >= 20 else 1
    if k >= 2:
        km = KMeans(n_clusters=k, random_state=7, n_init=10).fit(Z[:, :min(3, Z.shape[1])])
        labels = km.labels_
        cluster_info = []
        for c in range(k):
            mask = labels == c
            item = {"cluster": int(c), "count": int(mask.sum()), "fraction": float(mask.mean())}
            if shifted is not None:
                item["shifted_mean"] = float(np.mean(shifted[mask]))
                item["shifted_std"] = float(np.std(shifted[mask]))
            item["N_mean"] = float(np.mean(n_atoms[mask]))
            cluster_info.append(item)
    else:
        labels = np.zeros(len(Xs), dtype=int)
        cluster_info = [{"cluster": 0, "count": int(len(Xs)), "fraction": 1.0}]

    return {
        "available": True,
        "shape": [int(X.shape[0]), int(X.shape[1])],
        "pca_explained_variance_ratio": [float(v) for v in pca.explained_variance_ratio_],
        "pca_coordinates": Z,
        "nearest_neighbor_distance": percentile_dict(nn1, "desc_nn_"),
        "descriptor_isolated_fraction_z5": float(np.mean(iso_z > 5.0)) if len(iso_z) else float("nan"),
        "descriptor_isolated_indices": iso_idx[:200],
        "cluster_labels": labels,
        "cluster_info": cluster_info,
    }


# -----------------------------------------------------------------------------
# Optional overlap with test set
# -----------------------------------------------------------------------------

def compare_train_test_descriptor(train_desc: np.ndarray, test_desc: np.ndarray) -> Dict[str, Any]:
    if train_desc is None or test_desc is None:
        return {"available": False, "note": "train/test descriptors not both provided"}
    if train_desc.ndim != 2 or test_desc.ndim != 2 or train_desc.shape[1] != test_desc.shape[1]:
        return {"available": False, "note": "descriptor dimensions do not match"}
    scaler = StandardScaler().fit(train_desc)
    Xtr = scaler.transform(train_desc)
    Xte = scaler.transform(test_desc)
    nn = NearestNeighbors(n_neighbors=1).fit(Xtr)
    dist, _ = nn.kneighbors(Xte)
    dist = dist[:, 0]
    tr_nn = NearestNeighbors(n_neighbors=2).fit(Xtr).kneighbors(Xtr)[0][:, 1]
    threshold = float(np.percentile(tr_nn, 95)) if len(tr_nn) else float("nan")
    return {
        "available": True,
        "test_to_train_nn_distance": percentile_dict(dist, "test_train_nn_"),
        "train_internal_nn_p95": threshold,
        "test_fraction_outside_train_p95_nn": float(np.mean(dist > threshold)) if np.isfinite(threshold) else float("nan"),
    }


# -----------------------------------------------------------------------------
# Ratings: separate dimensions, no single quality score
# -----------------------------------------------------------------------------

def grade_geometry(geom: Dict[str, Any]) -> Tuple[str, List[str]]:
    notes = []
    severe = geom["severe_overlap_count"]
    warn = geom["warn_overlap_count"]
    if severe > 0:
        notes.append(f"{severe} configurations have very small pair distances: d/(r_cov_i+r_cov_j) < 0.55")
        return "FAIL", notes
    if warn > 0:
        notes.append(f"{warn} configurations have suspiciously small pair distances: 0.55 <= ratio < 0.70")
        return "WARN", notes
    p0 = geom["min_covalent_ratio"].get("ratio_p0", float("nan"))
    if np.isfinite(p0):
        notes.append(f"minimum covalent-radius-normalized distance ratio = {p0:.3f}")
    return "PASS", notes


def grade_labels(labels: Dict[str, Any]) -> Tuple[str, List[str]]:
    notes = []
    grade = "PASS"
    out_frac = labels.get("energy_per_atom_robust_outlier_fraction_z5", 0.0)
    if out_frac > 0.02:
        grade = "WARN"
        notes.append(f"energy per atom has {100*out_frac:.2f}% robust outliers beyond |z|>5")
    if labels.get("forces_found", False):
        f20 = labels.get("force_atom_fraction_gt_20", 0.0)
        f10 = labels.get("force_atom_fraction_gt_10", 0.0)
        fmax = labels.get("force_magnitude", {}).get("F_p100", float("nan"))
        if f20 > 0.001:
            grade = "WARN"
            notes.append(f"{100*f20:.3f}% atoms have |F| > 20 eV/A; inspect whether these are intended repulsive configurations")
        elif f10 > 0.01:
            grade = "WARN"
            notes.append(f"{100*f10:.2f}% atoms have |F| > 10 eV/A")
        if np.isfinite(fmax):
            notes.append(f"maximum force magnitude = {fmax:.3f} eV/A")
    else:
        grade = "UNKNOWN"
        notes.append("forces were not found; force-label sanity could not be evaluated")
    return grade, notes


def grade_descriptor(desc: Dict[str, Any]) -> Tuple[str, List[str]]:
    if not desc.get("available", False):
        return "NOT_PROVIDED", [desc.get("note", "descriptor not available")]
    notes = []
    iso = desc.get("descriptor_isolated_fraction_z5", 0.0)
    if iso > 0.03:
        grade = "WARN"
        notes.append(f"{100*iso:.2f}% configurations are isolated in descriptor space by robust nearest-neighbor criterion")
    else:
        grade = "PASS"
        notes.append(f"descriptor isolated fraction = {100*iso:.2f}%")
    evr = desc.get("pca_explained_variance_ratio", [])
    if evr:
        notes.append("PCA explained variance ratio = " + ", ".join(f"{v:.3f}" for v in evr[:5]))
    return grade, notes


def grade_composition(comp_info: Dict[str, Any], baseline: Dict[str, Any]) -> Tuple[str, List[str]]:
    notes = []
    div = comp_info["composition_diversity"]
    dtype = comp_info["dataset_type"]
    notes.append(f"dataset type: {dtype}")
    notes.append(f"composition diversity: {div}")
    if div in ("single-element", "very low", "low"):
        notes.append("atomic-reference baseline is mostly a constant shift here; rank deficiency should not be treated as data-quality failure")
        return "LOW_DIVERSITY_NOT_A_FAILURE", notes
    if comp_info["rank_counts"] < comp_info["n_elements"]:
        notes.append("composition-rich data still has rank-deficient element counts; element references are not unique")
        return "WARN", notes
    if np.isfinite(comp_info["cond_counts"]) and comp_info["cond_counts"] > 1e6:
        notes.append(f"composition matrix is ill-conditioned: cond={comp_info['cond_counts']:.2e}")
        return "WARN", notes
    notes.append("composition basis appears usable for element-reference diagnostics")
    return "PASS", notes


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------

def plot_baseline(outdir: Path, elements: List[str], counts: np.ndarray, comp: np.ndarray,
                  comp_info: Dict[str, Any], baseline: Dict[str, Any], n_atoms: np.ndarray):
    raw = baseline["raw_per_atom"]
    ref = baseline["ref_per_atom"]
    shifted = baseline["shifted_per_atom"]

    if comp.shape[1] >= 2 and len(comp) >= 2:
        try:
            comp_pc = PCA(n_components=1).fit_transform(comp).ravel()
        except Exception:
            comp_pc = np.zeros(len(comp))
    else:
        comp_pc = np.zeros(len(comp))

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    ax = axes[0, 0]
    ax.scatter(ref, raw, s=14, alpha=0.65, c=n_atoms)
    lo = min(np.nanmin(ref), np.nanmin(raw))
    hi = max(np.nanmax(ref), np.nanmax(raw))
    pad = 0.05 * max(hi - lo, 1e-12)
    ax.plot([lo-pad, hi+pad], [lo-pad, hi+pad], "--", lw=1.0, color="black", alpha=0.6)
    ax.set_xlabel("composition baseline E_ref/N (eV/atom)")
    ax.set_ylabel("raw E/N (eV/atom)")
    ax.set_title("Composition baseline: informational")

    text = (
        f"mode={baseline['mode']}\n"
        f"dataset type:\n{comp_info['dataset_type']}\n"
        f"rank={comp_info['rank_counts']}/{comp_info['n_elements']}\n"
        f"R2 total={baseline['r2_total']:.4f}\n"
        f"R2/atom={baseline['r2_per_atom']:.4f}"
    )
    ax.text(0.03, 0.97, text, transform=ax.transAxes, va="top", fontsize=8,
            bbox=dict(boxstyle="round", fc="white", ec="0.8", alpha=0.9))

    ax = axes[0, 1]
    ax.scatter(ref, shifted, s=14, alpha=0.65, c=n_atoms)
    ax.axhline(0, ls="--", lw=1.0, color="black", alpha=0.6)
    ax.set_xlabel("composition baseline E_ref/N (eV/atom)")
    ax.set_ylabel("residual E/N (eV/atom)")
    ax.set_title("Residual after composition shift")

    ax = axes[1, 0]
    ax.scatter(n_atoms, shifted, s=14, alpha=0.65)
    ax.axhline(0, ls="--", lw=1.0, color="black", alpha=0.6)
    ax.set_xlabel("number of atoms")
    ax.set_ylabel("residual E/N (eV/atom)")
    ax.set_title("Residual vs system size")

    ax = axes[1, 1]
    ax.scatter(comp_pc, shifted, s=14, alpha=0.65, c=n_atoms)
    ax.axhline(0, ls="--", lw=1.0, color="black", alpha=0.6)
    ax.set_xlabel("composition PC1")
    ax.set_ylabel("residual E/N (eV/atom)")
    ax.set_title("Residual vs composition coordinate")

    fig.suptitle("Atomic-reference baseline diagnosis (not a universal quality score)", fontweight="bold")
    fig.tight_layout()
    fig.savefig(outdir / "01_composition_baseline_informational.png", dpi=220)
    plt.close(fig)


def plot_geometry_labels(outdir: Path, raw_pa: np.ndarray, shifted: np.ndarray, geom: Dict[str, Any], labels: Dict[str, Any]):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    axes[0, 0].hist(raw_pa, bins="auto", alpha=0.8)
    axes[0, 0].set_xlabel("raw E/N (eV/atom)")
    axes[0, 0].set_ylabel("count")
    axes[0, 0].set_title("Energy per atom distribution")

    axes[0, 1].hist(shifted, bins="auto", alpha=0.8)
    axes[0, 1].set_xlabel("composition-shifted residual E/N (eV/atom)")
    axes[0, 1].set_ylabel("count")
    axes[0, 1].set_title("Residual distribution")

    ratio = geom["min_ratio_array"]
    axes[1, 0].hist(ratio[np.isfinite(ratio)], bins="auto", alpha=0.8)
    axes[1, 0].axvline(0.55, ls="--", lw=1.0, color="black")
    axes[1, 0].axvline(0.70, ls="--", lw=1.0, color="black")
    axes[1, 0].set_xlabel("min d / (r_cov_i + r_cov_j)")
    axes[1, 0].set_ylabel("count")
    axes[1, 0].set_title("Geometry overlap check")

    if labels.get("forces_found", False):
        fmax = labels["force_cfg_max_array"]
        axes[1, 1].hist(fmax[np.isfinite(fmax)], bins="auto", alpha=0.8)
        axes[1, 1].axvline(10.0, ls="--", lw=1.0, color="black")
        axes[1, 1].axvline(20.0, ls="--", lw=1.0, color="black")
        axes[1, 1].set_xlabel("max |F| per configuration (eV/A)")
        axes[1, 1].set_title("Force tail check")
    else:
        axes[1, 1].axis("off")
        axes[1, 1].text(0.5, 0.5, "forces not found", ha="center", va="center")
    axes[1, 1].set_ylabel("count")

    fig.suptitle("Geometry and label sanity", fontweight="bold")
    fig.tight_layout()
    fig.savefig(outdir / "02_geometry_and_label_sanity.png", dpi=220)
    plt.close(fig)


def plot_composition_pca(outdir: Path, comp: np.ndarray, shifted: np.ndarray, n_atoms: np.ndarray):
    if comp.shape[1] < 2 or len(comp) < 2:
        return
    ncomp = min(2, comp.shape[1], len(comp))
    Z = PCA(n_components=ncomp).fit_transform(comp)
    if Z.shape[1] == 1:
        Z = np.column_stack([Z[:, 0], np.zeros(len(Z))])
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    sc = axes[0].scatter(Z[:, 0], Z[:, 1], c=shifted, s=14, alpha=0.75)
    axes[0].set_xlabel("composition PC1")
    axes[0].set_ylabel("composition PC2")
    axes[0].set_title("Composition PCA colored by residual")
    fig.colorbar(sc, ax=axes[0], label="residual E/N")
    sc = axes[1].scatter(Z[:, 0], Z[:, 1], c=n_atoms, s=14, alpha=0.75)
    axes[1].set_xlabel("composition PC1")
    axes[1].set_ylabel("composition PC2")
    axes[1].set_title("Composition PCA colored by N")
    fig.colorbar(sc, ax=axes[1], label="N")
    fig.tight_layout()
    fig.savefig(outdir / "03_composition_pca.png", dpi=220)
    plt.close(fig)


def plot_descriptor_pca(outdir: Path, desc: Dict[str, Any], shifted: np.ndarray, n_atoms: np.ndarray):
    if not desc.get("available", False):
        return
    Z = desc["pca_coordinates"]
    if Z.shape[1] == 1:
        Z = np.column_stack([Z[:, 0], np.zeros(len(Z))])
    labels = desc.get("cluster_labels", np.zeros(len(Z), dtype=int))
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    sc = axes[0].scatter(Z[:, 0], Z[:, 1], c=shifted, s=14, alpha=0.75)
    axes[0].set_title("Descriptor PCA: residual")
    axes[0].set_xlabel("descriptor PC1")
    axes[0].set_ylabel("descriptor PC2")
    fig.colorbar(sc, ax=axes[0], label="residual E/N")
    sc = axes[1].scatter(Z[:, 0], Z[:, 1], c=n_atoms, s=14, alpha=0.75)
    axes[1].set_title("Descriptor PCA: N")
    axes[1].set_xlabel("descriptor PC1")
    axes[1].set_ylabel("descriptor PC2")
    fig.colorbar(sc, ax=axes[1], label="N")
    sc = axes[2].scatter(Z[:, 0], Z[:, 1], c=labels, s=14, alpha=0.75)
    axes[2].set_title("Descriptor PCA: clusters")
    axes[2].set_xlabel("descriptor PC1")
    axes[2].set_ylabel("descriptor PC2")
    fig.colorbar(sc, ax=axes[2], label="cluster")
    fig.tight_layout()
    fig.savefig(outdir / "04_descriptor_pca_coverage.png", dpi=220)
    plt.close(fig)


# -----------------------------------------------------------------------------
# Output files
# -----------------------------------------------------------------------------

def write_config_summary(outdir: Path, atoms_list, elements: List[str], counts: np.ndarray, comp: np.ndarray,
                         energies: np.ndarray, baseline: Dict[str, Any], geom: Dict[str, Any], labels: Dict[str, Any]):
    n_atoms = counts.sum(axis=1).astype(int)
    raw_pa = baseline["raw_per_atom"]
    shifted = baseline["shifted_per_atom"]
    force_max = labels.get("force_cfg_max_array", np.full(len(atoms_list), np.nan))
    force_mean = labels.get("force_cfg_mean_array", np.full(len(atoms_list), np.nan))
    path = outdir / "config_summary.csv"
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        header = [
            "index", "N", "chemical_system", "formula", "energy_total", "energy_per_atom",
            "residual_per_atom", "robust_z_residual", "min_distance_A", "min_covalent_ratio",
            "min_pair", "volume_per_atom", "force_mean", "force_max",
        ] + [f"frac_{e}" for e in elements]
        writer.writerow(header)
        rz = robust_mad_z(shifted)
        for i, atoms in enumerate(atoms_list):
            row = [
                i, n_atoms[i], chemical_system(atoms.get_chemical_symbols()),
                reduced_formula_from_counts(elements, counts[i].astype(int)),
                energies[i], raw_pa[i], shifted[i], rz[i],
                geom["min_distance_array"][i], geom["min_ratio_array"][i], geom["min_pair_labels"][i],
                geom["volume_per_atom_array"][i], force_mean[i], force_max[i],
            ] + list(comp[i])
            writer.writerow(row)


def write_outlier_summary(outdir: Path, atoms_list, elements: List[str], counts: np.ndarray,
                          baseline: Dict[str, Any], geom: Dict[str, Any], labels: Dict[str, Any], top_n: int = 100):
    shifted = baseline["shifted_per_atom"]
    rz = robust_mad_z(shifted)
    order = np.argsort(-np.abs(rz))[:top_n]
    force_max = labels.get("force_cfg_max_array", np.full(len(atoms_list), np.nan))
    path = outdir / "outlier_configs_reference.csv"
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "index", "abs_robust_z_residual", "residual_per_atom", "raw_energy_per_atom", "N",
            "formula", "min_distance_A", "min_covalent_ratio", "force_max", "comment",
        ])
        for i in order:
            comments = []
            if geom["min_ratio_array"][i] < 0.70:
                comments.append("short-distance")
            if np.isfinite(force_max[i]) and force_max[i] > 10:
                comments.append("high-force")
            writer.writerow([
                int(i), float(abs(rz[i])), float(shifted[i]), float(baseline["raw_per_atom"][i]),
                int(counts[i].sum()), reduced_formula_from_counts(elements, counts[i].astype(int)),
                float(geom["min_distance_array"][i]), float(geom["min_ratio_array"][i]),
                float(force_max[i]) if np.isfinite(force_max[i]) else "NA",
                ";".join(comments),
            ])


def markdown_report(outdir: Path, input_file: str, energy_source: str, force_source: str,
                    elements: List[str], comp_info: Dict[str, Any], baseline: Dict[str, Any],
                    geom: Dict[str, Any], label_info: Dict[str, Any], desc: Dict[str, Any],
                    ratings: Dict[str, Any], overlap: Optional[Dict[str, Any]] = None):
    lines = []
    lines.append("# Scientific Dataset Diagnosis Report")
    lines.append("")
    lines.append(f"Input file: `{input_file}`")
    lines.append(f"Energy source: `{energy_source}`")
    lines.append(f"Force source: `{force_source}`")
    lines.append(f"Number of elements: {len(elements)}")
    lines.append(f"Elements: {' '.join(elements)}")
    lines.append("")
    lines.append("> This report is a reference diagnosis, not a verdict. It avoids a single quality score because fixed-stoichiometry, composition-rich, and mixed-source datasets require different criteria.")
    lines.append("")

    lines.append("## Dataset type and applicability")
    lines.append("")
    lines.append(f"- dataset type: **{comp_info['dataset_type']}**")
    lines.append(f"- composition diversity: **{comp_info['composition_diversity']}**")
    lines.append(f"- composition count rank: {comp_info['rank_counts']} / {comp_info['n_elements']}")
    cond = comp_info["cond_counts"]
    lines.append(f"- composition condition number: {'inf' if not np.isfinite(cond) else f'{cond:.3e}'}")
    lines.append("- top reduced formulas:")
    for item in comp_info["top_formulas"][:8]:
        lines.append(f"  - {item['formula']}: {item['count']}")
    lines.append("")

    lines.append("## Dimension ratings")
    lines.append("")
    for key, item in ratings.items():
        lines.append(f"### {key}: {item['grade']}")
        for n in item["notes"]:
            lines.append(f"- {n}")
        lines.append("")

    lines.append("## Atomic-reference baseline, informational only")
    lines.append("")
    lines.append(f"- fit mode: {baseline['mode']}")
    for e, mu in zip(elements, baseline["mu"]):
        lines.append(f"- mu({e}) = {mu:.8f} eV/atom")
    lines.append(f"- R2 total = {fmt_float(baseline['r2_total'])}")
    lines.append(f"- R2 per atom = {fmt_float(baseline['r2_per_atom'])}")
    lines.append(f"- residual RMSE = {1000*baseline['rmse_shifted_eV_atom']:.3f} meV/atom")
    lines.append(f"- residual MAE = {1000*baseline['mae_shifted_eV_atom']:.3f} meV/atom")
    lo, hi = baseline["range_shifted_eV_atom"]
    lines.append(f"- residual range = [{lo:.6f}, {hi:.6f}] eV/atom")
    if comp_info["composition_diversity"] in ("single-element", "very low", "low"):
        lines.append("- interpretation: low composition diversity means this baseline is mostly a constant shift; rank deficiency or low R2/atom is not a data-quality failure by itself.")
    lines.append("")

    lines.append("## Geometry sanity")
    lines.append("")
    lines.append(f"- atom count: p50={fmt_float(geom['atom_count'].get('N_p50'), 1)}, range=[{fmt_float(geom['atom_count'].get('N_p0'), 0)}, {fmt_float(geom['atom_count'].get('N_p100'), 0)}]")
    lines.append(f"- PBC patterns: {geom['pbc_counts']}")
    lines.append(f"- min distance: p0={fmt_float(geom['min_distance'].get('dmin_p0'), 4)} A, p1={fmt_float(geom['min_distance'].get('dmin_p1'), 4)} A")
    lines.append(f"- min covalent ratio: p0={fmt_float(geom['min_covalent_ratio'].get('ratio_p0'), 4)}, p1={fmt_float(geom['min_covalent_ratio'].get('ratio_p1'), 4)}")
    lines.append(f"- severe short-distance configs: {geom['severe_overlap_count']}")
    lines.append(f"- warning short-distance configs: {geom['warn_overlap_count']}")
    lines.append("")

    lines.append("## Label sanity")
    lines.append("")
    lines.append(f"- energy per atom p50 = {fmt_float(label_info['energy_per_atom'].get('Epa_p50'))} eV/atom")
    lines.append(f"- energy per atom p0/p100 = {fmt_float(label_info['energy_per_atom'].get('Epa_p0'))} / {fmt_float(label_info['energy_per_atom'].get('Epa_p100'))} eV/atom")
    lines.append(f"- E_total vs N Pearson = {fmt_float(label_info.get('E_total_vs_N_pearson'), 4)}")
    lines.append(f"- E/N vs N Pearson = {fmt_float(label_info.get('Epa_vs_N_pearson'), 4)}")
    if label_info.get("forces_found", False):
        F = label_info["force_magnitude"]
        lines.append(f"- |F| p50/p95/p99/p100 = {fmt_float(F.get('F_p50'), 4)} / {fmt_float(F.get('F_p95'), 4)} / {fmt_float(F.get('F_p99'), 4)} / {fmt_float(F.get('F_p100'), 4)} eV/A")
        lines.append(f"- atom fraction |F| > 10 eV/A = {100*label_info['force_atom_fraction_gt_10']:.4f}%")
        lines.append(f"- atom fraction |F| > 20 eV/A = {100*label_info['force_atom_fraction_gt_20']:.4f}%")
    else:
        lines.append("- forces not found")
    lines.append("")

    lines.append("## Descriptor coverage")
    lines.append("")
    if desc.get("available", False):
        lines.append(f"- descriptor shape: {desc['shape'][0]} x {desc['shape'][1]}")
        lines.append("- PCA explained variance: " + ", ".join(f"{v:.4f}" for v in desc["pca_explained_variance_ratio"][:5]))
        lines.append(f"- descriptor isolated fraction by robust NN z>5: {100*desc['descriptor_isolated_fraction_z5']:.3f}%")
        lines.append("- cluster summary:")
        for c in desc["cluster_info"]:
            extra = ""
            if "shifted_mean" in c:
                extra = f", residual_mean={c['shifted_mean']:.4f}, residual_std={c['shifted_std']:.4f}, N_mean={c['N_mean']:.1f}"
            lines.append(f"  - cluster {c['cluster']}: count={c['count']} ({100*c['fraction']:.1f}%){extra}")
    else:
        lines.append(f"- {desc.get('note', 'descriptor not provided')}")
    lines.append("")

    if overlap:
        lines.append("## Train/test descriptor overlap")
        lines.append("")
        if overlap.get("available", False):
            lines.append(f"- train internal NN p95 = {fmt_float(overlap['train_internal_nn_p95'], 4)}")
            p = overlap["test_to_train_nn_distance"]
            lines.append(f"- test-to-train NN p50/p95/p99 = {fmt_float(p.get('test_train_nn_p50'),4)} / {fmt_float(p.get('test_train_nn_p95'),4)} / {fmt_float(p.get('test_train_nn_p99'),4)}")
            lines.append(f"- test fraction outside train p95 NN radius = {100*overlap['test_fraction_outside_train_p95_nn']:.2f}%")
        else:
            lines.append(f"- {overlap.get('note', 'not available')}")
        lines.append("")

    lines.append("## Suggested use")
    lines.append("")
    lines.append("1. Treat `01_composition_baseline_informational.png` as a composition diagnostic only; do not use it alone to judge fixed-stoichiometry datasets.")
    lines.append("2. Inspect `02_geometry_and_label_sanity.png` and `config_summary.csv` for short-distance, high-force, or formatting/unit problems.")
    lines.append("3. If `descriptor.out` is available, use `04_descriptor_pca_coverage.png` to find isolated structural families and decide whether group-wise sampling is needed.")
    lines.append("4. Use `outlier_configs_reference.csv` as a checklist for visualization, not as an automatic deletion list.")

    (outdir / "scientific_diagnosis_report.md").write_text("\n".join(lines), encoding="utf-8")


def json_sanitize(obj: Any) -> Any:
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, dict):
        return {str(k): json_sanitize(v) for k, v in obj.items() if not str(k).endswith("_array") and k not in ("pca_coordinates", "cluster_labels", "min_pair_labels")}
    if isinstance(obj, list):
        return [json_sanitize(v) for v in obj]
    return obj


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Scientific pre-training diagnosis for atomistic ML-potential datasets."
    )
    parser.add_argument("xyz", help="input extxyz/xyz dataset")
    parser.add_argument("--descriptor", default=None, help="descriptor.out matching the input xyz")
    parser.add_argument("--outdir", default="scientific_dataset_diagnosis", help="output directory")
    parser.add_argument("--fit-mode", choices=["wls", "ols"], default="wls", help="atomic-reference baseline fit mode; informational")
    parser.add_argument("--test-xyz", default=None, help="optional test/validation xyz for overlap reference")
    parser.add_argument("--test-descriptor", default=None, help="optional descriptor.out for test xyz")
    args = parser.parse_args()

    outdir = ensure_dir(args.outdir)
    atoms_list = read_xyz(args.xyz)
    energies, energy_source = read_energies(atoms_list)
    forces, force_source = read_forces(atoms_list)

    elements, counts, comp = composition_matrix(atoms_list)
    n_atoms = counts.sum(axis=1)
    comp_info = analyze_composition(elements, counts, comp)
    baseline = fit_atomic_reference(counts, energies, mode=args.fit_mode)
    geom = analyze_geometry(atoms_list)
    label_info = analyze_labels(energies, n_atoms, forces)

    desc_arr = safe_array_loadtxt(args.descriptor) if args.descriptor else None
    desc_info = analyze_descriptor(desc_arr, len(atoms_list), baseline["shifted_per_atom"], n_atoms)

    overlap_info = None
    if args.test_descriptor and args.descriptor:
        test_desc = safe_array_loadtxt(args.test_descriptor)
        overlap_info = compare_train_test_descriptor(desc_arr, test_desc)
    elif args.test_xyz is not None:
        overlap_info = {"available": False, "note": "test xyz provided, but descriptor overlap requires --descriptor and --test-descriptor"}

    ratings = {}
    g, notes = grade_composition(comp_info, baseline)
    ratings["Composition diagnosis"] = {"grade": g, "notes": notes}
    g, notes = grade_geometry(geom)
    ratings["Geometry sanity"] = {"grade": g, "notes": notes}
    g, notes = grade_labels(label_info)
    ratings["Label sanity"] = {"grade": g, "notes": notes}
    g, notes = grade_descriptor(desc_info)
    ratings["Descriptor coverage"] = {"grade": g, "notes": notes}

    write_config_summary(outdir, atoms_list, elements, counts, comp, energies, baseline, geom, label_info)
    write_outlier_summary(outdir, atoms_list, elements, counts, baseline, geom, label_info)

    plot_baseline(outdir, elements, counts, comp, comp_info, baseline, n_atoms)
    plot_geometry_labels(outdir, baseline["raw_per_atom"], baseline["shifted_per_atom"], geom, label_info)
    plot_composition_pca(outdir, comp, baseline["shifted_per_atom"], n_atoms)
    plot_descriptor_pca(outdir, desc_info, baseline["shifted_per_atom"], n_atoms)

    markdown_report(outdir, args.xyz, energy_source, force_source, elements, comp_info, baseline,
                    geom, label_info, desc_info, ratings, overlap_info)

    summary = {
        "input": args.xyz,
        "energy_source": energy_source,
        "force_source": force_source,
        "composition": comp_info,
        "baseline_informational": baseline,
        "geometry": geom,
        "labels": label_info,
        "descriptor": desc_info,
        "ratings": ratings,
        "overlap": overlap_info,
    }
    (outdir / "summary.json").write_text(json.dumps(json_sanitize(summary), indent=2), encoding="utf-8")

    print(f"[OK] Scientific dataset diagnosis written to: {outdir}")
    print(f"[OK] Report: {outdir / 'scientific_diagnosis_report.md'}")
    print("[NOTE] No single quality score is produced. Use the dimension ratings and plots as references.")


if __name__ == "__main__":
    main()
