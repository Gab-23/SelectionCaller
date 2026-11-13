#!/usr/bin/env python
"""
BURNER - quali geni driver finiscono in cima, per residuo e per osservato (21/10/2026)
=====================================================================================
Per ogni omica ordina i bin per osservato e per residuo (stessa direzione del controllo
residual_vs_raw_drivers.py) e riporta i geni dei bin driver nel top 1% e 5%, cosa il residuo
guadagna o perde rispetto all'osservato, e una tabella dei primi N bin con i loro geni.

Geni per bin: geni protein_coding del GTF di GENCODE che si sovrappongono al bin (backbone BED).
Bin "driver": distanza 0 da un OG o TSG (colonne di sola validazione del dataset), come nello
script dei driver. Se si passa --intogen, i nomi dei driver sono presi da quella lista (colonna
dei simboli riconosciuta in automatico); altrimenti si riportano tutti i geni del bin driver.

Uso:
  python driver_genes_in_top.py --oof pilot_v3_genome_nopos/oof_predictions.tsv --dataset $DS \
     --backbone resources/backbone/hg38_binned_20000_0.mask.filtered.sorted.bed \
     --gtf resources/reference_data/gencode/gencode.v45.annotation.gtf.gz \
     --intogen resources/IntOGen/IntOGen-Drivers-20240920/Compendium_Cancer_Genes.tsv --tumor_type BRCA \
     --out_dir driver_genes_v3_genome_nopos
"""
import argparse, os, re
import numpy as np, pandas as pd

p = argparse.ArgumentParser()
p.add_argument("--oof", required=True); p.add_argument("--dataset", required=True)
p.add_argument("--backbone", required=True); p.add_argument("--gtf", required=True)
p.add_argument("--intogen", default=None); p.add_argument("--tumor_type", default="BRCA")
p.add_argument("--top_n", type=int, default=30); p.add_argument("--out_dir", default="driver_genes")
a = p.parse_args(); os.makedirs(a.out_dir, exist_ok=True)

DIRECTION = {"cna_ampl": ("OG",), "cna_del": ("TSG",), "mut": ("OG", "TSG"), "atac": ("OG", "TSG")}
ABS = {"atac"}
strip = lambda c: c.str.replace("^chr", "", regex=True)

# --- geni protein_coding dal GTF ---
genes = pd.read_csv(a.gtf, sep="\t", header=None, comment="#", usecols=[0, 2, 3, 4, 8],
                    names=["chr", "type", "start", "end", "attr"], dtype={"chr": str})
genes = genes[(genes["type"] == "gene") & genes["attr"].str.contains('gene_type "protein_coding"', regex=False)]
genes["name"] = genes["attr"].str.extract(r'gene_name "([^"]+)"')[0]
genes["chr"] = strip(genes["chr"]); genes = genes[["chr", "start", "end", "name"]]

# --- backbone e geni per bin ---
bb = pd.read_csv(a.backbone, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chr", "start", "end", "binID"], dtype={"chr": str})
bb["chr"] = strip(bb["chr"]); bb["start"] += 1                       # BED 0-based -> 1-based chiuso
bin_genes = {}
for ch, b in bb.groupby("chr"):
    g = genes[genes["chr"] == ch]
    if g.empty: continue
    b = b.sort_values("start"); bs, be, bid = b["start"].to_numpy(), b["end"].to_numpy(), b["binID"].to_numpy()
    for s, e, n in g[["start", "end", "name"]].itertuples(index=False):
        i0 = max(np.searchsorted(be, s, side="left"), 0); i1 = np.searchsorted(bs, e, side="right")
        for i in range(i0, i1):
            if bs[i] <= e and be[i] >= s: bin_genes.setdefault(bid[i], set()).add(n)

# --- lista driver (opzionale) ---
drivers = None
if a.intogen:
    ig = pd.read_csv(a.intogen, sep="\t")
    sym = next((c for c in ig.columns if c.upper() in ("SYMBOL", "GENE", "GENE_SYMBOL", "HUGO_SYMBOL")), None)
    ct = next((c for c in ig.columns if c.upper() in ("CANCER_TYPE", "CANCER", "TUMOR_TYPE")), None)
    if sym is None: print(f"ATTENZIONE: colonna dei simboli non trovata in {a.intogen} ({list(ig.columns)[:8]}...)")
    else:
        if ct is not None: ig = ig[ig[ct].astype(str).str.upper() == a.tumor_type.upper()]
        drivers = set(ig[sym].astype(str))
        print(f"Driver IntOGen per {a.tumor_type}: {len(drivers)} geni")

D = pd.read_csv(a.dataset, sep="\t", usecols=["binID", "dist.to.closest.OG", "dist.to.closest.TSG"])
D["is_OG"] = D["dist.to.closest.OG"].fillna(np.inf) == 0; D["is_TSG"] = D["dist.to.closest.TSG"].fillna(np.inf) == 0
M = pd.read_csv(a.oof, sep="\t").merge(D[["binID", "is_OG", "is_TSG"]], on="binID", how="left")

def names(bid, only_drivers):
    g = bin_genes.get(bid, set())
    if only_drivers and drivers is not None: g = g & drivers
    return ",".join(sorted(g))

lines = []
for om, g in M.groupby("omic"):
    if om not in DIRECTION: continue
    g = g.copy(); g["is_driver"] = g[[f"is_{c}" for c in DIRECTION[om]]].any(axis=1)
    sc = {"osservato": g["obs"], "residuo": g["resid"]}
    if om in ABS: sc = {k: v.abs() for k, v in sc.items()}
    lines.append(f"\n===== {om} ({'+'.join(DIRECTION[om])}; {len(g)} bin, {int(g['is_driver'].sum())} bin driver) =====")
    sets = {}
    for frac in (0.01, 0.05):
        k = max(1, int(round(frac * len(g))))
        for sname, s in sc.items():
            top = g.loc[s.sort_values(ascending=False).index[:k]]
            dg = sorted({x for bid in top.loc[top["is_driver"], "binID"] for x in names(bid, True).split(",") if x})
            sets[(frac, sname)] = set(dg)
            lines.append(f"top {int(frac*100)}% per {sname:9s} ({k} bin): {int(top['is_driver'].sum())} bin driver -> {', '.join(dg) if dg else '-'}")
        gain = sorted(sets[(frac, 'residuo')] - sets[(frac, 'osservato')]); lost = sorted(sets[(frac, 'osservato')] - sets[(frac, 'residuo')])
        lines.append(f"   residuo guadagna: {', '.join(gain) if gain else '-'} | perde: {', '.join(lost) if lost else '-'}")
    for sname, s in sc.items():
        top = g.loc[s.sort_values(ascending=False).index[:a.top_n]].copy()
        top["genes"] = top["binID"].map(lambda b: names(b, False)); top["rank"] = range(1, len(top) + 1)
        top[["rank", "binID", "obs", "pred", "resid", "is_driver", "genes"]].to_csv(f"{a.out_dir}/{om}_top{a.top_n}_by_{sname}.tsv", sep="\t", index=False)
txt = "\n".join(lines); print(txt); open(f"{a.out_dir}/summary.txt", "w").write(txt + "\n")
print(f"\nTabelle dei primi {a.top_n} bin per omica e punteggio in {a.out_dir}/")