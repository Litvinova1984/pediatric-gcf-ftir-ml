#!/usr/bin/env python3
"""Verify checks/PACKAGE_SHA256.csv using only Python 3 standard library."""
from pathlib import Path
import csv
import hashlib

root = Path(__file__).resolve().parents[1]
manifest = root / "checks" / "PACKAGE_SHA256.csv"
with manifest.open(newline="", encoding="utf-8") as fp:
    rows = list(csv.DictReader(fp))
failed = []
expected = set()
for row in rows:
    rel = row["path"]
    expected.add(rel)
    item = root / rel
    if (not item.is_file()
        or hashlib.sha256(item.read_bytes()).hexdigest() != row["sha256"]):
        failed.append(rel)
def is_generated_or_local(rel: str) -> bool:
    parts = rel.split("/")
    return (rel == "checks/PACKAGE_SHA256.csv" or
            ".git" in parts or "__MACOSX" in parts or
            (parts and parts[0] in (".Rproj.user", "private", "restricted")) or
            rel.startswith(("R/results/", "data/private/")) or
            rel.endswith((".DS_Store", ".Rhistory", ".RData")) or
            any(part.startswith("._") for part in parts))

actual = {rel for p in root.rglob("*") if p.is_file()
          for rel in [p.relative_to(root).as_posix()]
          if not is_generated_or_local(rel)}
failed.extend(sorted(actual - expected))
failed.extend(sorted(expected - actual))
if failed:
    raise SystemExit("SHA-256 FAIL: " + ", ".join(failed))
print(f"PASS: {len(rows)} files verified against SHA-256 manifest")
