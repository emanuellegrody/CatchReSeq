#!/usr/bin/env python3
"""
Export per-cell isoform1/isoform2 counts from a BRIE2 brie_count.h5ad to CSV.

Usage:
  python export_brie_per_cell_EGS024.py <brie_count.h5ad> <out.csv> [event_id ...]

If no event_ids are passed, all events in adata.var_names are exported.

Output columns:
  cellID, event, isoform1, isoform2, total, psi
"""

import sys
import numpy as np
import pandas as pd
import scanpy as sc


def column(layer, idx):
    col = layer[:, idx]
    if hasattr(col, "toarray"):
        return col.toarray().ravel()
    if hasattr(col, "todense"):
        return np.asarray(col.todense()).ravel()
    return np.asarray(col).ravel()


def main():
    if len(sys.argv) < 3:
        print("Usage: python export_brie_per_cell_EGS024.py "
              "<brie_count.h5ad> <out.csv> [event_id ...]")
        sys.exit(1)

    h5ad_path = sys.argv[1]
    out_csv   = sys.argv[2]
    events    = sys.argv[3:] if len(sys.argv) > 3 else None

    adata = sc.read_h5ad(h5ad_path)
    if events is None:
        events = list(adata.var_names)

    rows = []
    for event_id in events:
        if event_id not in adata.var_names:
            print(f"WARNING: {event_id} not in {h5ad_path}, skipping")
            continue

        idx = adata.var_names.get_loc(event_id)
        iso1 = column(adata.layers["isoform1"], idx)
        iso2 = column(adata.layers["isoform2"], idx)
        total = iso1 + iso2
        with np.errstate(divide="ignore", invalid="ignore"):
            psi = np.where(total > 0, iso1 / total, np.nan)

        for i, barcode in enumerate(adata.obs_names):
            rows.append({
                "cellID":   barcode,
                "event":    event_id,
                "isoform1": int(iso1[i]),
                "isoform2": int(iso2[i]),
                "total":    int(total[i]),
                "psi":      psi[i],
            })

    df = pd.DataFrame(rows)
    df.to_csv(out_csv, index=False)
    print(f"Wrote {len(df)} rows ({df['cellID'].nunique()} cells × "
          f"{df['event'].nunique()} events) to {out_csv}")


if __name__ == "__main__":
    main()
