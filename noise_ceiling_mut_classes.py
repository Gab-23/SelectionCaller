#!/usr/bin/env python
"""
BURNER - tetto di rumore per i target delle mutazioni PER CLASSE (sinonime, missense, troncanti)
=================================================================================================
Stessa idea di noise_ceiling.py (che resta invariato per CNA e per il target vecchio), applicata ai nuovi
target del pilota: EVENTI (campione x mutazione) per classe, rapportati ai SITI della classe del bin.

Tutto dentro il pilota, NESSUNA dipendenza da regole nuove della pipeline: la classificazione degli eventi e'
riscritta qui in Python con le STESSE regole di compute_coding_mutations_counts.R (NOTE 28/10/2026: con
--events_classified le classi vengono invece dalla riannotazione sui siti di gene_level_drivers_classes.py, come la
pipeline con event_annotation = reference; senza, dalle righe COSMIC come prima):
  - sostituzioni singole (GENOME_START == GENOME_STOP, alleli di una base ACGT);
  - classe dal trascritto Ensembl_canonical (--transcript_mode canonical) o la piu' grave tra tutti (all);
    gravita': troncante > missense > altra codificante > sinonima > non codificante;
  - un evento = (campione, mutazione), assegnato al bin di GENOME_START;
  - esclusi i campioni con piu' di --hypermutator_max eventi codificanti (sostituzioni e indel).
CONTROLLO DI COERENZA (stampato e salvato): i totali per bin e classe calcolati qui devono coincidere con le
colonne mut_n_syn / mut_n_mis / mut_n_trunc del dataset. Se non coincidono, la classificazione qui e' diversa
da quella della pipeline e il tetto NON va usato.

Tetto: affidabilita' tra pazienti della DENSITA' per sito (eventi / siti della classe), sui bin del pilota
(cds_bp >= --min_cds e siti della classe > 0):
  - alfa di Cronbach: item = pazienti (conteggi di eventi nel bin / siti), soggetti = bin;
  - split-half: pazienti divisi a meta' (--n_splits volte), correlazione tra le densita' delle due meta',
    corretta con Spearman-Brown.
Il tetto e' sulla scala della densita' (R2 di un predittore della densita'), non della devianza spiegata
riportata dal pilota: per confrontarli serve l'R2 del pilota sulla densita' (eventi/siti), calcolabile dalle
predizioni out-of-fold.

Uso (env del pilota):
  python noise_ceiling_mut_classes.py --events resources/COSMIC/mutations/coding/f_breast_carcinoma.filtered.bed \
     --gtf resources/reference_data/gencode/gencode.v45.annotation.gtf.gz \
     --backbone resources/backbone/hg38_binned_20000_0.mask.filtered.sorted.bed --dataset $DS \
     --hypermutator_max 3000 --transcript_mode canonical --exclude_drivers --out_dir ceiling_mut_classes_20kb
"""
import argparse, os, sys, gzip, re
import numpy as np, pandas as pd
from scipy import sparse

p = argparse.ArgumentParser()
p.add_argument("--events", required=True); p.add_argument("--gtf", required=True)
p.add_argument("--backbone", required=True); p.add_argument("--dataset", required=True)
p.add_argument("--hypermutator_max", default="3000"); p.add_argument("--transcript_mode", choices=["canonical", "all"], default="canonical")
p.add_argument("--min_cds", type=float, default=500); p.add_argument("--n_splits", type=int, default=200)
p.add_argument("--exclude_drivers", action="store_true", help="calcola anche il tetto senza i bin driver (OG/TSG a distanza 0)")
p.add_argument("--events_classified", default=None,
               help="(consigliato) events_classified.tsv di gene_level_drivers_classes.py: classi dalla RIANNOTAZIONE sui siti, "
                    "come la pipeline con event_annotation = reference (stessa soglia degli ipermutati). Senza: classi dalle righe COSMIC")
p.add_argument("--seed", type=int, default=1); p.add_argument("--out_dir", default="ceiling_mut_classes")
a = p.parse_args(); os.makedirs(a.out_dir, exist_ok=True); rng = np.random.default_rng(a.seed)
strip = lambda s: pd.Series(s, dtype=str).str.replace("^chr", "", regex=True).to_numpy()
LOG = []
def log(m): print(m, file=sys.stderr); LOG.append(m)

# ---------------- classificazione (stesse regole di compute_coding_mutations_counts.R) ----------------
TRUNC = {"stop_gained", "splice_acceptor_variant", "splice_donor_variant"}; MIS = {"missense_variant"}
OTHER = {"start_lost", "stop_lost", "stop_retained_variant", "coding_sequence_variant", "protein_altering_variant",
         "incomplete_terminal_codon_variant", "start_retained_variant"}
SYN = {"synonymous_variant"}; INDEL = {"frameshift_variant", "inframe_insertion", "inframe_deletion"}
CLASS_OF_SEV = {2: "syn", 3: "other_coding", 4: "mis", 5: "trunc"}
def severity(desc):
    t = set(str(desc).split(","))
    return 5 if t & TRUNC else 4 if t & MIS else 3 if t & OTHER else 2 if t & SYN else 1
def coding_any(desc): return bool(set(str(desc).split(",")) & (TRUNC | MIS | OTHER | SYN | INDEL))

if a.events_classified is not None:
    # NOTE (28/10/2026): classi dalla riannotazione (file di gene_level_drivers_classes.py, stessa classificazione della pipeline);
    # il file contiene tutti gli eventi dei campioni tenuti dopo la soglia degli ipermutati, con la classe per le SNV.
    E = pd.read_csv(a.events_classified, sep="\t", dtype={"COSMIC_SAMPLE_ID": str, "mut_id": str, "chr": str})
    log(f"eventi classificati letti da {a.events_classified}: {len(E)} (classi dalla riannotazione sui siti)")
    samples = np.array(sorted(E["COSMIC_SAMPLE_ID"].unique())); s_idx = {s: i for i, s in enumerate(samples)}
    E = E[E["cls"].isin(["syn", "mis", "trunc"])].copy()
    log(f"eventi per classe: {E['cls'].value_counts().to_dict()} su {len(samples)} pazienti")
if a.events_classified is None:
    cols = ["COSMIC_SAMPLE_ID", "GENOMIC_MUTATION_ID", "CHROMOSOME", "GENOME_START", "GENOME_STOP", "TRANSCRIPT_ACCESSION",
            "MUTATION_DESCRIPTION", "GENOMIC_WT_ALLELE", "GENOMIC_MUT_ALLELE"]
    R = pd.read_csv(a.events, sep="\t", dtype=str, quoting=3, usecols=cols, keep_default_na=False)
    R["mut_id"] = np.where(R["GENOMIC_MUTATION_ID"] == "", R["CHROMOSOME"] + ":" + R["GENOME_START"] + ":" + R["GENOME_STOP"], R["GENOMIC_MUTATION_ID"])
    R["pos"] = pd.to_numeric(R["GENOME_START"], errors="coerce"); R["stop"] = pd.to_numeric(R["GENOME_STOP"], errors="coerce")
    R = R[R["pos"].notna() & (R["COSMIC_SAMPLE_ID"] != "")]
    udesc = pd.Series(R["MUTATION_DESCRIPTION"].unique())
    sev_map = dict(zip(udesc, udesc.map(severity))); cod_map = dict(zip(udesc, udesc.map(coding_any)))
    R["sev"] = R["MUTATION_DESCRIPTION"].map(sev_map); R["coding_any"] = R["MUTATION_DESCRIPTION"].map(cod_map)
    B = set("ACGT")
    R["snv"] = (R["pos"] == R["stop"]) & R["GENOMIC_WT_ALLELE"].str.upper().isin(B) & R["GENOMIC_MUT_ALLELE"].str.upper().isin(B)
    if a.transcript_mode == "canonical":
        canon = set()
        with gzip.open(a.gtf, "rt") as f:
            for line in f:
                if '\ttranscript\t' in line and 'tag "Ensembl_canonical"' in line:
                    m = re.search(r'transcript_id "([^".]+)', line)
                    if m: canon.add(m.group(1))
        log(f"trascritti Ensembl_canonical: {len(canon)}")
        R["in_mode"] = R["TRANSCRIPT_ACCESSION"].str.split(".").str[0].isin(canon)
    else:
        R["in_mode"] = True
    R["sev_mode"] = np.where(R["in_mode"], R["sev"], 0)
    E = R.groupby(["COSMIC_SAMPLE_ID", "mut_id"], sort=False).agg(chr=("CHROMOSOME", "first"), pos=("pos", "first"),
            snv=("snv", "first"), coding_any=("coding_any", "any"), sev=("sev_mode", "max")).reset_index()
    del R
    per_sample = E.groupby("COSMIC_SAMPLE_ID")["coding_any"].sum()
    if str(a.hypermutator_max).lower() not in ("none", ""):
        keep_s = per_sample[per_sample <= float(a.hypermutator_max)].index
        log(f"campioni: {len(per_sample)}, esclusi come ipermutati: {len(per_sample) - len(keep_s)} (soglia {a.hypermutator_max} eventi codificanti)")
        E = E[E["COSMIC_SAMPLE_ID"].isin(keep_s)]
    samples = np.array(sorted(E["COSMIC_SAMPLE_ID"].unique())); s_idx = {s: i for i, s in enumerate(samples)}
    E["cls"] = [CLASS_OF_SEV.get(s) if (snv and s >= 2) else None for s, snv in zip(E["sev"], E["snv"])]
    E = E[E["cls"].isin(["syn", "mis", "trunc"])]
    log(f"eventi per classe: {E['cls'].value_counts().to_dict()} su {len(samples)} pazienti")

# ---------------- bin ----------------
bb = pd.read_csv(a.backbone, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chr", "start", "end", "binID"], dtype={"chr": str, "binID": str})
bb["chr"] = strip(bb["chr"]); bb["start"] += 1
E["chr"] = strip(E["chr"]); E["bin"] = -1
bb = bb.reset_index(drop=True)
for ch, b in bb.groupby("chr"):
    m = (E["chr"] == ch).to_numpy()
    if not m.any(): continue
    b = b.sort_values("start"); i = np.searchsorted(b["start"].to_numpy(), E.loc[m, "pos"].to_numpy(), side="right") - 1
    ok = (i >= 0) & (E.loc[m, "pos"].to_numpy() <= b["end"].to_numpy()[np.clip(i, 0, None)])
    E.loc[m, "bin"] = np.where(ok, b.index.to_numpy()[np.clip(i, 0, None)], -1)
E = E[E["bin"] >= 0]

D = pd.read_csv(a.dataset, sep="\t", dtype={"binID": str}).set_index("binID")
need = ["cds_bp"] + [f"sites_{c}" for c in ("syn", "mis", "trunc")] + [f"mut_n_{c}" for c in ("syn", "mis", "trunc")]
miss = [c for c in need if c not in D.columns]
if miss: sys.exit(f"colonne mancanti nel dataset: {miss}")

# ---------------- controllo di coerenza con la pipeline ----------------
bins = bb["binID"].to_numpy()
coh = []
for c in ("syn", "mis", "trunc"):
    here = pd.Series(np.bincount(E.loc[E["cls"] == c, "bin"], minlength=len(bb)), index=bins)
    pipe = D[f"mut_n_{c}"].reindex(bins).fillna(0).astype(int)
    common = D.index.intersection(bins)
    diff = int((here.loc[common] != pipe.loc[common]).sum())
    coh.append(dict(cls=c, bins=len(common), bins_diff=diff, total_here=int(here.loc[common].sum()), total_pipeline=int(pipe.loc[common].sum())))
    log(f"COERENZA {c}: bin diversi {diff} su {len(common)}; totale qui {int(here.loc[common].sum())}, pipeline {int(pipe.loc[common].sum())}")
pd.DataFrame(coh).to_csv(f"{a.out_dir}/coherence_with_pipeline.tsv", sep="\t", index=False)
if any(r["bins_diff"] for r in coh):
    log("ATTENZIONE: la classificazione NON coincide con quella della pipeline: il tetto qui sotto non va usato")

# ---------------- tetto ----------------
def alpha(M):
    """Cronbach: item = colonne (pazienti), soggetti = righe (bin). M sparse (CSR), valori = eventi / siti."""
    k = M.shape[1]; n = M.shape[0]
    T = np.asarray(M.sum(axis=1)).ravel(); vt = T.var(ddof=1)
    s1 = np.asarray(M.sum(axis=0)).ravel(); s2 = np.asarray(M.multiply(M).sum(axis=0)).ravel()
    vi = ((s2 - s1 ** 2 / n) / (n - 1)).sum()
    return k / (k - 1) * (1 - vi / vt) if vt > 0 and k > 1 else np.nan

def split_half(M, n):
    k = M.shape[1]; out = []
    for _ in range(n):
        g = np.zeros(k, bool); g[rng.permutation(k)[: k // 2]] = True
        f1 = np.asarray(M[:, g].sum(axis=1)).ravel(); f2 = np.asarray(M[:, ~g].sum(axis=1)).ravel()
        if f1.std() == 0 or f2.std() == 0: continue
        r = np.corrcoef(f1, f2)[0, 1]; out.append(2 * r / (1 + r))
    return np.array(out)

driver = np.zeros(len(D), bool)
if a.exclude_drivers:
    for c in ("dist.to.closest.OG", "dist.to.closest.TSG"):
        if c in D.columns: driver |= (D[c].fillna(np.inf) == 0).to_numpy()
driver_bins = set(D.index[driver])

rows = []
for c in ("syn", "mis", "trunc"):
    ec = E[E["cls"] == c]
    X = sparse.csr_matrix((np.ones(len(ec)), (ec["bin"].to_numpy(), ec["COSMIC_SAMPLE_ID"].map(s_idx).to_numpy())), shape=(len(bb), len(samples)))
    X.sum_duplicates()
    sites = D[f"sites_{c}"].reindex(bins).to_numpy(float); cds = D["cds_bp"].reindex(bins).to_numpy(float)
    in_ds = pd.Series(bins).isin(D.index).to_numpy()
    base = in_ds & (cds >= a.min_cds) & (sites > 0)                              # stessi bin del pilota
    variants = [("tutti i bin", base)]
    if a.exclude_drivers: variants.append(("senza driver", base & ~pd.Series(bins).isin(driver_bins).to_numpy()))
    for label, keep in variants:
        M = sparse.diags(1.0 / sites[keep]) @ X[keep]                           # densita' per sito, per paziente
        sh = split_half(M, a.n_splits)
        r = dict(cls=c, bins=label, n_bins=int(keep.sum()), n_patients=len(samples), events=int(X[keep].sum()),
                 events_per_bin=float(X[keep].sum() / max(keep.sum(), 1)), alpha=alpha(M),
                 split_half_SB=float(np.mean(sh)) if len(sh) else np.nan,
                 split_half_lo=float(np.percentile(sh, 2.5)) if len(sh) else np.nan,
                 split_half_hi=float(np.percentile(sh, 97.5)) if len(sh) else np.nan)
        rows.append(r)
        log(f"TETTO {c:5s} {label:12s}: bin {r['n_bins']}, eventi {r['events']} ({r['events_per_bin']:.2f}/bin), "
            f"alfa {r['alpha']:.3f}, split-half {r['split_half_SB']:.3f} [{r['split_half_lo']:.3f}, {r['split_half_hi']:.3f}]")
pd.DataFrame(rows).to_csv(f"{a.out_dir}/noise_ceiling_mut_classes.tsv", sep="\t", index=False)
open(f"{a.out_dir}/log.txt", "w").write("\n".join(LOG) + "\n")
print(f"Output in {a.out_dir}/", file=sys.stderr)