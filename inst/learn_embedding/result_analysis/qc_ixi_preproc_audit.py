#!/usr/bin/env python3
import argparse, csv, random, os
import numpy as np
import nibabel as nib
from pathlib import Path
from nibabel.orientations import aff2axcodes

def robust_stats(data):
    d = data[np.isfinite(data)]
    if d.size == 0:
        return dict(min=np.nan, max=np.nan, mean=np.nan, std=np.nan, p01=np.nan, p99=np.nan, nonzero=0.0)
    nz = np.count_nonzero(d) / d.size
    return dict(
        min=float(np.nanmin(d)),
        max=float(np.nanmax(d)),
        mean=float(np.nanmean(d)),
        std =float(np.nanstd(d)),
        p01 =float(np.nanpercentile(d, 1)),
        p99 =float(np.nanpercentile(d, 99)),
        nonzero=float(nz),
    )

def main():
    ap = argparse.ArgumentParser("Audit IXI T1 preprocessing needs")
    ap.add_argument("--csv", required=True, help="all_labeled.csv with columns id,path,y")
    ap.add_argument("--img_col", default="path")
    ap.add_argument("--id_col", default="id")
    ap.add_argument("--n", type=int, default=50, help="How many images to audit (random sample)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out_csv", default="qc_summary.csv")
    args = ap.parse_args()

    # load rows
    rows = []
    with open(args.csv, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            rows.append(row)
    random.seed(args.seed)
    sample = random.sample(rows, min(args.n, len(rows)))

    out_rows = []
    for row in sample:
        p = row[args.img_col]
        iid = row.get(args.id_col, Path(p).stem)
        try:
            img = nib.load(p)
            data = img.get_fdata(dtype=np.float32)
            zooms = img.header.get_zooms()[:3]
            ax = "".join(aff2axcodes(img.affine))
            shp = img.shape[:3]
            st = robust_stats(data)
            out_rows.append({
                "id": iid,
                "path": p,
                "shape": f"{shp[0]}x{shp[1]}x{shp[2]}",
                "zoom_x": zooms[0], "zoom_y": zooms[1], "zoom_z": zooms[2],
                "orientation": ax,
                "min": st["min"], "max": st["max"], "mean": st["mean"], "std": st["std"],
                "p01": st["p01"], "p99": st["p99"], "nonzero_frac": st["nonzero"],
            })
        except Exception as e:
            out_rows.append({"id": iid, "path": p, "error": str(e)})

    # write per-image CSV
    outp = Path(args.out_csv)
    outp.parent.mkdir(parents=True, exist_ok=True)
    cols = ["id","path","shape","orientation","zoom_x","zoom_y","zoom_z","min","max","mean","std","p01","p99","nonzero_frac","error"]
    with open(outp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in out_rows:
            for c in cols:
                r.setdefault(c, "")
            w.writerow(r)

    # quick aggregates
    ori = [r["orientation"] for r in out_rows if r.get("orientation")]
    shapes = [r["shape"] for r in out_rows if r.get("shape")]
    zooms = [(float(r["zoom_x"]), float(r["zoom_y"]), float(r["zoom_z"])) for r in out_rows if r.get("zoom_x")!=""]
    nz = [float(r["nonzero_frac"]) for r in out_rows if r.get("nonzero_frac")!=""]
    p01 = [float(r["p01"]) for r in out_rows if r.get("p01")!=""]
    p99 = [float(r["p99"]) for r in out_rows if r.get("p99")!=""]

    def uniq(seq): return sorted(set(seq))

    print("\n=== QC summary ===")
    print(f"N audited: {len(out_rows)}  | wrote: {outp}")
    print(f"Orientations (unique): {uniq(ori)}")
    if zooms:
        zx, zy, zz = zip(*zooms)
        print(f"Zooms x: min={min(zx):.3f}, max={max(zx):.3f}  "
              f"y: min={min(zy):.3f}, max={max(zy):.3f}  "
              f"z: min={min(zz):.3f}, max={max(zz):.3f}")
    print(f"Shapes (unique up to sample): {uniq(shapes)[:8]}{' …' if len(uniq(shapes))>8 else ''}")
    if nz:
        print(f"Nonzero frac: median={np.median(nz):.3f}, range=({min(nz):.3f},{max(nz):.3f})")
    if p01 and p99:
        print(f"Intensity p01~p99 (median): {np.median(p01):.2f} .. {np.median(p99):.2f}")

if __name__ == "__main__":
    main()
