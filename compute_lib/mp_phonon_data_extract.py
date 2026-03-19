# 从material project中获取声子谱以及声子DOS数据放到json中
# 使用方式：python judge.py mp-2741  将会从material project中获取mp-2741的声子谱以及声子DOS数据放到json中

from mp_api.client import MPRester
import json
import os
import gzip
import urllib.request
from pathlib import Path
import sys
import argparse

API_KEY = "5NxP3sMUSrESWwk7cih4G6yjS8Br7ToB"
import gzip
import urllib.request
from pathlib import Path

#mpid = "mp-2741"
parser = argparse.ArgumentParser(description="Download phonon bandstructure and DOS from Material Project.")
parser.add_argument("mpid", type=str, help="Material Project ID (e.g., mp-2741)")
parser.add_argument("method", type=str, default="dfpt", help="Phonon method (e.g., dfpt)")
args = parser.parse_args()

mpid = args.mpid
method = args.method

BASES = [
    "https://materialsproject-parsed.s3.amazonaws.com/",
]

targets = {
    f"{mpid}_phonon_bs_{method}.json": f"ph-bandstructures/{method}/{mpid}.json.gz",
    f"{mpid}_phonon_dos_{method}.json": f"ph-dos/{method}/{mpid}.json.gz",
}

def download(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read()

for outname, rel in targets.items():
    url = BASES[0] + rel
    print("Downloading:", url)
    raw_gz = download(url)
    raw_json = gzip.decompress(raw_gz)
    Path(outname).write_bytes(raw_json)
    print("Saved:", outname)