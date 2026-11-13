#!/usr/bin/env python
"""
BURNER - il nostro background come COVARIATA di dNdScv.
=======================================================
Per ogni gene: log del tasso neutrale per sito previsto dal NOSTRO modello di background per le sinonime
    tasso = atteso delle sinonime del gene (gene_scores_syn.tsv, colonna E) / siti sinonimi del gene
I siti sinonimi del gene (non pesati) vengono da gene_bin_weighted_sites.tsv di context_weighted_sites.py (colonna sites).
Uso le sinonime perche' dNdScv allena la regressione sulle sinonime: la covariata deve descrivere il tasso neutrale.
Il nostro modello delle sinonime e' allenato sulle sole sinonime (neutrali), quindi non porta selezione nella covariata.
Geni senza atteso (copertura insufficiente) o senza siti: lasciati vuoti qui; run_dndscv.R li riempie con la mediana
(dNdScv vuole una riga per ogni gene del suo RefCDS) e lo dichiara nel log.
Uscita: TSV gene_name, our_log_rate_syn.
"""
import argparse, sys
import numpy as np, pandas as pd
p = argparse.ArgumentParser()
p.add_argument("--gene_scores_syn", required=True, help="drivers_gene_level/gene_scores_syn.tsv del nostro run")
p.add_argument("--gene_sites", default=None,
               help="(facoltativo) gene_bin_weighted_sites.tsv di context_weighted_sites.py; senza: colonna 'sites' di gene_scores_syn.tsv "
                    "(scritta da gene_level_drivers_classes.py dal 28/10/2026)")
p.add_argument("--out", required=True)
a = p.parse_args()
G = pd.read_csv(a.gene_scores_syn, sep="\t", index_col=0)
if a.gene_sites is not None:
    S = pd.read_csv(a.gene_sites, sep="\t"); sites = S[S["cls"] == "syn"].groupby("gene")["sites"].sum()
elif "sites" in G.columns:
    sites = G["sites"].astype(float)
else:
    sys.exit("gene_scores_syn.tsv senza colonna 'sites' (versione vecchia di gene_level_drivers_classes.py): passare --gene_sites")
E = G["E"].astype(float)
rate = E / sites.reindex(E.index)
ok = np.isfinite(rate) & (rate > 0)
out = pd.DataFrame({"gene_name": E.index, "our_log_rate_syn": np.where(ok, np.log(rate.where(ok, 1.0)), np.nan)})
out.to_csv(a.out, sep="\t", index=False)
print(f"geni: {len(out)}; con covariata: {int(ok.sum())}; mancanti (verranno riempiti con la mediana): {int((~ok).sum())}", file=sys.stderr)
print(f"log tasso per sito: mediana {np.nanmedian(out.our_log_rate_syn):.2f}, 5-95%: "
      f"{np.nanpercentile(out.our_log_rate_syn, 5):.2f} .. {np.nanpercentile(out.our_log_rate_syn, 95):.2f}", file=sys.stderr)