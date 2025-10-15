#!/usr/bin/env python3
# Prepare IXI T1 data with AGE labels from IXI.xls
# - Extracts the IXI .tar/.tar.gz into out_dir/raw
# - Finds T1 NIfTI files
# - Joins with IXI.xls by IXI_ID to get AGE as y
# - Writes:
#     out_dir/csv/all_labeled.csv   (id,path,y)
#     out_dir/csv/all_unlabeled.csv (id,path)
# - Optional: --sample N  -> out_dir/csv/U_aux.csv with N random rows from labeled set

import argparse, os, tarfile, re, csv, random
from pathlib import Path

import pandas as pd

def extract_tar(tar_path: Path, out_root: Path) -> Path:
    raw_dir = out_root / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    # Skip extraction if already extracted
    if any(raw_dir.iterdir()):
        print(f"[info] {raw_dir} not empty; assuming already extracted.")
        return raw_dir
    mode = "r:gz" if str(tar_path).lower().endswith((".tar.gz",".tgz")) else "r:"
    with tarfile.open(tar_path, mode) as tf:
        tf.extractall(raw_dir)
    return raw_dir

def find_t1_files(root: Path):
    t1s = []
    for p in root.rglob("*"):
        if not p.is_file():
            continue
        s = p.name.lower()
        if (s.endswith(".nii") or s.endswith(".nii.gz")) and ("t1" in s):
            t1s.append(p)
    return sorted(t1s)

def load_age_map(xls_path: Path):
    # pandas needs 'xlrd' for .xls (pip install xlrd)
    try:
        df = pd.read_excel(xls_path)  # engine auto-chosen
    except Exception as e:
        raise SystemExit(
            f"Failed to read {xls_path}. If it's .xls, install 'xlrd' (pip install xlrd). Error: {e}"
        )
    # Normalize column names
    norm = {c: re.sub(r"[\s\-]+", "_", str(c).strip().lower()) for c in df.columns}
    df.columns = [norm[c] for c in df.columns]
    # Find id + age columns robustly
    id_col = next((c for c in df.columns if c in ("ixi_id","ixi","id","subject_id")), None)
    age_col = next((c for c in df.columns if c.startswith("age")), None)
    if id_col is None or age_col is None:
        raise SystemExit(f"Could not find IXI_ID and AGE columns in {xls_path}. Columns={df.columns.tolist()}")
    # Clean + map to dict: key=int IXI_ID (no leading zeros) -> float AGE
    df = df[[id_col, age_col]].dropna()
    # Some IXI_ID are strings; extract digits and cast to int
    df["_id_int"] = (
        df[id_col]
        .astype(str)
        .str.extract(r"(\d+)", expand=False)
        .astype(int)
    )
    df["_age"] = pd.to_numeric(df[age_col], errors="coerce")
    df = df.dropna(subset=["_age"])
    age_map = {int(i): float(a) for i, a in zip(df["_id_int"], df["_age"])}
    return age_map

def ixi_id_from_path(p: Path):
    # Extract IXI number from filename or any parent string, e.g. "IXI002" -> 2, "IXI115" -> 115
    m = re.search(r"IXI(\d+)", str(p), flags=re.IGNORECASE)
    return int(m.group(1)) if m else None

def main():
    ap = argparse.ArgumentParser("Prepare IXI T1 with AGE labels")
    ap.add_argument("--tar", required=True, help="Path to IXI .tar or .tar.gz with T1 images")
    ap.add_argument("--xls", required=True, help="Path to IXI.xls with metadata (AGE)")
    ap.add_argument("--out_dir", required=True, help="Output directory")
    ap.add_argument("--sample", type=int, default=0,
                    help="If >0, also write U_aux.csv with N random labeled rows (default: 0 = disabled)")
    ap.add_argument("--seed", type=int, default=1337)
    args = ap.parse_args()

    random.seed(args.seed)
    tar_path = Path(args.tar)
    xls_path = Path(args.xls)
    out_dir  = Path(args.out_dir)
    csv_dir  = out_dir / "csv"
    csv_dir.mkdir(parents=True, exist_ok=True)

    print(f"[1/4] Extracting (or reusing) archive: {tar_path}")
    raw_root = extract_tar(tar_path, out_dir)
    print(f"      -> raw root: {raw_root}")

    print("[2/4] Loading AGE from Excel …")
    age_map = load_age_map(xls_path)
    print(f"      -> ages for {len(age_map)} IXI IDs")

    print("[3/4] Scanning for T1 images …")
    t1s = find_t1_files(raw_root)
    print(f"      -> found {len(t1s)} T1 files")

    labeled_rows = []
    missing_age, missing_id = 0, 0
    for p in t1s:
        iid = ixi_id_from_path(p)
        if iid is None:
            missing_id += 1
            continue
        age = age_map.get(int(iid))
        if age is None:
            missing_age += 1
            continue
        # Canonical ID like "IXI002"
        id_str = f"IXI{int(iid):03d}"
        labeled_rows.append((id_str, str(p), float(age)))

    if missing_id:
        print(f"[warn] {missing_id} files without a parsable IXI ID (skipped).")
    if missing_age:
        print(f"[warn] {missing_age} files had no AGE in {xls_path} (skipped).")

    if not labeled_rows:
        raise SystemExit("No labeled rows produced. Check that your T1 filenames include 'T1' and IDs like 'IXI###'.")

    # Write all_labeled.csv
    labeled_csv = csv_dir / "all_labeled.csv"
    with open(labeled_csv, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["id","path","y"])  # y = AGE
        w.writerows(labeled_rows)
    print(f"[4/4] Wrote labeled CSV: {labeled_csv}  (rows={len(labeled_rows)})")

    # Write all_unlabeled.csv
    unlabeled_csv = csv_dir / "all_unlabeled.csv"
    with open(unlabeled_csv, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["id","path"])
        for rid, p, _ in labeled_rows:
            w.writerow([rid, p])
    print(f"      Wrote unlabeled CSV: {unlabeled_csv}  (rows={len(labeled_rows)})")

    # Optional: small subset for a quick U_aux smoke test (disabled by default)
    if args.sample and args.sample > 0:
        k = min(args.sample, len(labeled_rows))
        subset = random.sample(labeled_rows, k)
        uaux_csv = csv_dir / "U_aux.csv"
        with open(uaux_csv, "w", newline="") as f:
            w = csv.writer(f); w.writerow(["id","path","y"])
            w.writerows(sorted(subset))
        print(f"      Also wrote: {uaux_csv} (rows={k})")

    print("Done.")

if __name__ == "__main__":
    main()
