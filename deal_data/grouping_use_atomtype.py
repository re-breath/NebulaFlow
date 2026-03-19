#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
按照原子类型将xyz分组
Read any XYZ/extxyz using ASE, group atoms by element species, and write an
EXTXYZ file that includes per-atom integer property "group".

Output:
  - <out>.xyz : extxyz with Properties including group
  - <map_out> : json mapping element -> group_id and counts

Example:
  python ase_write_grouped_extxyz.py model.xyz -o model_grouped.xyz

Notes:
  - group ids are assigned deterministically by sorted element symbols by default
    (use --order appearance to assign by first occurrence).
  - For multi-frame xyz trajectories, select a frame with --frame (default -1 = last).
"""

import argparse
import json
from collections import Counter, OrderedDict

import numpy as np
from ase.io import read, write


def build_mapping(symbols, order: str, start_id: int):
    if order == "sorted":
        elems = sorted(set(symbols))
    elif order == "appearance":
        seen = OrderedDict()
        for s in symbols:
            if s not in seen:
                seen[s] = None
        elems = list(seen.keys())
    else:
        raise ValueError(f"Unknown order: {order}")

    return {elem: start_id + i for i, elem in enumerate(elems)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("xyz", help="Input XYZ / extxyz file path")
    ap.add_argument("-o", "--out", default="model_grouped.xyz",
                    help="Output EXTXYZ file with per-atom property 'group'")
    ap.add_argument("-m", "--map-out", default="group_map.json",
                    help="Output JSON mapping file (element->group_id, counts)")
    ap.add_argument("--order", choices=["sorted", "appearance"], default="sorted",
                    help="Group id assignment order: sorted (alphabetical) or appearance (first seen)")
    ap.add_argument("--start-id", type=int, default=1, help="Starting group id (default: 1)")
    ap.add_argument("--frame", type=int, default=-1,
                    help="If XYZ contains multiple frames, choose which frame to use "
                         "(default -1 means last frame)")
    args = ap.parse_args()

    atoms = read(args.xyz, index=args.frame)
    symbols = atoms.get_chemical_symbols()
    counts = Counter(symbols)

    mapping = build_mapping(symbols, args.order, args.start_id)
    group_ids = np.array([mapping[s] for s in symbols], dtype=np.int32)

    # Attach per-atom property "group" to atoms
    atoms.arrays["group"] = group_ids

    # Write as extxyz so the group property is preserved
    write(args.out, atoms, format="extxyz")

    # Write mapping json for reference
    payload = {
        "input_file": args.xyz,
        "frame_used": args.frame,
        "order": args.order,
        "start_id": args.start_id,
        "element_to_group": mapping,
        "element_counts": dict(counts),
    }
    with open(args.map_out, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)

    print(f"[OK] Read {len(symbols)} atoms from {args.xyz} (frame {args.frame})")
    print("[SUMMARY] element -> group_id (count)")
    for elem in (sorted(mapping) if args.order == "sorted" else mapping.keys()):
        print(f"  {elem:>3s} -> {mapping[elem]:>3d}  (N={counts[elem]})")
    print(f"[OK] Wrote EXTXYZ with group property: {args.out}")
    print(f"[OK] Wrote mapping JSON           : {args.map_out}")


if __name__ == "__main__":
    main()
