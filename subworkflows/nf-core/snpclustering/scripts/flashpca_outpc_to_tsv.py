#!/usr/bin/env python3
import argparse, json
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outpc", required=True, help="File prodotto da flashpca2 --outpc")
    ap.add_argument("--n-pcs", type=int, required=True)
    ap.add_argument("--out-pca", required=True, help="TSV output: sample_id + PC1..PCn")
    ap.add_argument("--out-info", required=True, help="JSON output info")
    ap.add_argument("--id-mode", choices=["fid_iid", "iid"], default="fid_iid")
    args = ap.parse_args()

    outpc = Path(args.outpc)
    n_pcs = args.n_pcs

    rows = []
    with outpc.open() as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            parts = s.split()
            # se c'è header (es. "FID IID PC1 ...") lo saltiamo
            if parts[0].upper() in ("FID", "#FID") or parts[1].upper() == "IID":
                continue
            if len(parts) < 2 + n_pcs:
                raise SystemExit(f"ERROR: riga con poche colonne in outpc: {parts[:5]} ...")

            fid, iid = parts[0], parts[1]
            sample_id = iid if args.id_mode == "iid" else f"{fid}:{iid}"
            pcs = parts[2:2+n_pcs]
            rows.append((sample_id, pcs))

    header = ["sample_id"] + [f"PC{i}" for i in range(1, n_pcs+1)]
    with Path(args.out_pca).open("w") as w:
        w.write("\t".join(header) + "\n")
        for sid, pcs in rows:
            w.write(sid + "\t" + "\t".join(pcs) + "\n")

    info = {"tool": "flashpca2", "n_pcs": n_pcs, "n_samples": len(rows)}
    Path(args.out_info).write_text(json.dumps(info, indent=2))

if __name__ == "__main__":
    main()
