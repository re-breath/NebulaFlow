#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 去重 ref_tidy.bib 中的重复条目


import re
from pathlib import Path
import sys

INPUT = sys.argv[1]
OUTPUT = "ref_dedup.bib"

def split_entries(text):
    parts = re.split(r'(?=^[ \t]*@\w+\s*\{)', text, flags=re.M)
    return [p for p in parts if p.strip()]

def get_key(entry):
    m = re.match(r'^[ \t]*@\w+\s*\{\s*([^,]+)\s*,', entry, flags=re.M)
    return m.group(1).strip() if m else None

def main():
    text = Path(INPUT).read_text(encoding="utf-8")
    entries = split_entries(text)

    seen = set()
    kept = []
    removed = []

    for entry in entries:
        key = get_key(entry)
        if key is None:
            kept.append(entry)
            continue
        if key in seen:
            removed.append(key)
        else:
            seen.add(key)
            kept.append(entry)

    Path(OUTPUT).write_text("\n".join(e.strip() for e in kept) + "\n", encoding="utf-8")

    print(f"[OK] 输出文件: {OUTPUT}")
    print(f"[OK] 保留条目数: {len(kept)}")
    print(f"[OK] 删除重复条目数: {len(removed)}")
    if removed:
        print("[INFO] 删除的重复 key：")
        for k in removed:
            print(k)

if __name__ == "__main__":
    main()
