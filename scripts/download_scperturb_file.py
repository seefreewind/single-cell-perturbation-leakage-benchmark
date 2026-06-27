#!/usr/bin/env python3
"""Download one file from the scPerturb Zenodo record with resume support."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from urllib.request import urlopen

import requests


ZENODO_RECORD = "https://zenodo.org/api/records/13350497"


def get_file_url(filename: str) -> str:
    with urlopen(ZENODO_RECORD, timeout=30) as response:
        record = json.load(response)
    for item in record.get("files", []):
        if item.get("key") == filename:
            return item["links"]["self"]
    available = ", ".join(item.get("key", "") for item in record.get("files", []))
    raise SystemExit(f"File not found: {filename}\nAvailable files:\n{available}")


def download(filename: str, outdir: Path, chunk_size: int = 1024 * 1024, dry_run: bool = False) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    url = get_file_url(filename)
    outfile = outdir / filename
    if dry_run:
        print(f"filename: {filename}")
        print(f"url: {url}")
        print(f"target: {outfile}")
        return
    existing = outfile.stat().st_size if outfile.exists() else 0
    headers = {"Range": f"bytes={existing}-"} if existing else {}
    mode = "ab" if existing else "wb"

    with requests.get(url, headers=headers, stream=True, timeout=60) as response:
        response.raise_for_status()
        total = response.headers.get("content-length")
        total_int = int(total) + existing if total and existing else int(total or 0)
        downloaded = existing
        with open(outfile, mode) as handle:
            for chunk in response.iter_content(chunk_size=chunk_size):
                if not chunk:
                    continue
                handle.write(chunk)
                downloaded += len(chunk)
                if total_int:
                    pct = downloaded / total_int * 100
                    print(f"\r{outfile.name}: {downloaded/1024**2:.1f} MB / {total_int/1024**2:.1f} MB ({pct:.1f}%)", end="")
                else:
                    print(f"\r{outfile.name}: {downloaded/1024**2:.1f} MB", end="")
    print()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    parser.add_argument("--outdir", type=Path, default=Path("data/raw/scperturb"))
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    download(args.filename, args.outdir, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
