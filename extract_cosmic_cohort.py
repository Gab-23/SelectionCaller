#!/usr/bin/env python
"""
BURNER - estrae da COSMIC le mutazioni di una coorte TCGA (stessa logica delle regole filter_classification_COSMIC,
filter_samples_COSMIC, get_final_samples_COSMIC e filter_mutations della pipeline), per lanciare dNdScv su altre coorti.
Perche' non la pipeline: il wildcard COSMIC_parameters si divide in sito e istologia sul carattere "_", quindi non funziona
con valori COSMIC che contengono "_" (es. "malignant_melanoma", "urinary_tract").
Regole (come la pipeline): campioni con SAMPLE_NAME che inizia con --sample_name_pattern (default TCGA), sesso --gender
(F, M o ALL), fenotipi con PRIMARY_SITE e PRIMARY_HISTOLOGY dati (piu', facoltativo, HISTOLOGY_SUBTYPE_1).
Differenza dalla pipeline: il filtro dei campioni nelle mutazioni e' sulla COLONNA COSMIC_SAMPLE_ID (match esatto),
non con grep -Fw sull'intera riga.
--list: stampa le combinazioni sito / istologia / sottotipo disponibili (con il numero di campioni TCGA) ed esce.
"""
import argparse, sys, gzip
import pandas as pd
p = argparse.ArgumentParser()
p.add_argument("--classification", required=True); p.add_argument("--samples", required=True); p.add_argument("--mutations")
p.add_argument("--site"); p.add_argument("--histology"); p.add_argument("--subtype", default=None)
p.add_argument("--gender", default="ALL", choices=["F", "M", "ALL"]); p.add_argument("--sample_name_pattern", default="TCGA")
p.add_argument("--list", action="store_true"); p.add_argument("--out")
a = p.parse_args()
C = pd.read_csv(a.classification, sep="\t", dtype=str)
S = pd.read_csv(a.samples, sep="\t", dtype=str)
S = S[S["SAMPLE_NAME"].fillna("").str.startswith(a.sample_name_pattern)]
if a.gender != "ALL": S = S[S["GENDER"] == a.gender]
if a.list:
    M = S.merge(C, on="COSMIC_PHENOTYPE_ID", how="inner")
    t = M.groupby(["PRIMARY_SITE", "PRIMARY_HISTOLOGY", "HISTOLOGY_SUBTYPE_1"], dropna=False)["COSMIC_SAMPLE_ID"].nunique()
    print(t.sort_values(ascending=False).head(60).to_string()); sys.exit(0)
if not (a.site and a.histology and a.mutations and a.out): sys.exit("servono --site, --histology, --mutations, --out (oppure --list)")
c = C[(C["PRIMARY_SITE"] == a.site) & (C["PRIMARY_HISTOLOGY"] == a.histology)]
if a.subtype: c = c[c["HISTOLOGY_SUBTYPE_1"] == a.subtype]
pheno = set(c["COSMIC_PHENOTYPE_ID"])
samples = set(S.loc[S["COSMIC_PHENOTYPE_ID"].isin(pheno), "COSMIC_SAMPLE_ID"])
print(f"fenotipi: {len(pheno)}; campioni ({a.sample_name_pattern}, sesso {a.gender}): {len(samples)}", file=sys.stderr)
if not samples: sys.exit("nessun campione: controlla sito/istologia/sottotipo con --list")
op = gzip.open if a.mutations.endswith(".gz") else open
n_in = n_out = 0
with op(a.mutations, "rt") as f, open(a.out, "w") as g:
    header = f.readline(); g.write(header); col = header.rstrip("\n").split("\t").index("COSMIC_SAMPLE_ID")
    for line in f:
        n_in += 1
        if line.split("\t", col + 1)[col] in samples: g.write(line); n_out += 1
print(f"righe lette {n_in}, scritte {n_out}; campioni con almeno una riga: "
      f"{pd.read_csv(a.out, sep=chr(9), usecols=['COSMIC_SAMPLE_ID'], dtype=str)['COSMIC_SAMPLE_ID'].nunique()}", file=sys.stderr)