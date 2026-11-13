#!/usr/bin/env python
"""
BURNER - driver per gene con il NUOVO target delle mutazioni (eventi per classe, offset dei siti della classe)
=============================================================================================================
Tutto dentro il pilota, nessuna regola della pipeline. gene_level_drivers.py resta invariato (target vecchio).

Per ogni gene g e classe c (sinonime, missense, troncanti) e per le non sinonime (missense + troncanti):
  OSSERVATO  O[g,c] = eventi (campione x mutazione, sostituzioni singole) della classe c NEL TRASCRITTO CANONICO del
             gene g. Classificazione con le stesse regole di compute_coding_mutations_counts.R (ipermutati esclusi
             oltre --hypermutator_max eventi codificanti). Un evento in due geni sovrapposti conta in ciascuno, con
             la classe che ha nel trascritto canonico di quel gene.
  NOTE 28/10/2026: la classe degli eventi viene dalla RIANNOTAZIONE sui siti (come la pipeline con event_annotation =
             reference), non dalle righe COSMIC; gene_scores_<classe>.tsv contiene anche i siti della classe nel gene.
  SITI       S[g,c,b] = sostituzioni possibili della classe c nel trascritto canonico di g che cadono nel bin b,
             contate dalla sequenza di riferimento con le stesse regole di compute_coding_sites.R (ORF = CDS +
             stop_codon; primo codone e stop = altra codificante; siti di splicing essenziali +-1/+-2 = troncanti;
             esclusi trascritti con cds_start_NF / cds_end_NF / seleno).
  ATTESO     E[g,c] = somma sui bin di tasso[b,c] x S[g,c,b], con tasso[b,c] = predizione OOF del pilota per mut_c /
             sites_c del bin. Se parte dei siti del gene cade in bin senza predizione (CDS < 500 bp), l'atteso e'
             estrapolato dalla parte coperta; sotto --min_cov il gene resta nell'universo ma senza punteggio (in fondo).
  PUNTEGGIO  residuo = -log10 P(X >= O), X ~ binomiale negativa (media E, dispersione theta stimata con il metodo dei
             momenti sui geni con atteso, escluso l'1% piu' estremo: nessuna lista di driver usata). Coda calcolata
             in modo esatto (funzione beta incompleta; somma in scala log quando va in underflow).
CONTROLLI DI COERENZA con la pipeline (stampati e salvati): (1) eventi per bin e classe = mut_n_* del dataset;
(2) siti per bin e classe (unione dei trascritti, classe piu' grave) = sites_* del dataset.
Valutazione: AP (pari merito come un unico gradino), AUROC, IC con bootstrap sui geni, stesso universo per osservato e
residuo; driver = liste IntOGen (OG + TSG). Le sinonime fanno da controllo: non dovrebbero ritrovare i driver.

Uso (env del pilota, ~30 GB di memoria per il genoma):
  python gene_level_drivers_classes.py --oof analyses/pilot_classes/20000bp/genome_nopos/oof_predictions.tsv \
    --dataset $DS --backbone resources/backbone/hg38_binned_20000_0.mask.filtered.sorted.bed \
    --events resources/COSMIC/mutations/coding/f_breast_carcinoma.filtered.bed \
    --gtf resources/reference_data/gencode/gencode.v45.annotation.gtf.gz --genome resources/reference_data/genome/hg38.fa \
    --og resources/IntOGen/OG/hg38_v45_BRCA.sorted.bed --tsg resources/IntOGen/TSG/hg38_v45_BRCA.sorted.bed \
    --out_dir analyses/pilot_classes/20000bp/drivers_gene_level
"""
import argparse, os, sys, gzip, re, itertools
import numpy as np, pandas as pd
from scipy import stats, special
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

p = argparse.ArgumentParser()
for k in ("oof", "dataset", "backbone", "events", "gtf", "genome", "og", "tsg"): p.add_argument(f"--{k}", required=True)
p.add_argument("--hypermutator_max", default="3000"); p.add_argument("--min_cov", type=float, default=0.3)
p.add_argument("--n_boot", type=int, default=1000); p.add_argument("--seed", type=int, default=1)
p.add_argument("--top_n", type=int, default=30); p.add_argument("--out_dir", default="drivers_gene_level_classes")
p.add_argument("--wsites_bin", default=None, help="con --exposure wsites: bin_weighted_sites.tsv di context_weighted_sites.py")
p.add_argument("--wsites_gene", default=None, help="gene_bin_weighted_sites.tsv di context_weighted_sites.py (per l'atteso 'wsites' e per --exposure wsites)")
p.add_argument("--exposure", choices=["sites", "cds", "wsites"], default="sites",
               help="deve coincidere con --class_exposure del pilota: sites = atteso del gene = tasso per sito del bin x siti "
                    "della classe del gene; cds = tasso per bp di CDS del bin (pred / cds_bp) x bp dell'ORF del gene nel bin")
a = p.parse_args(); os.makedirs(a.out_dir, exist_ok=True)
rng = np.random.default_rng(a.seed); TIE = np.random.default_rng(a.seed + 1)
LOG = []
def log(m): print(m, file=sys.stderr); LOG.append(str(m))
strip = lambda x: re.sub("^chr", "", str(x))
CLASSES = ("syn", "mis", "trunc"); SEV2CLS = {2: "syn", 3: "other_coding", 4: "mis", 5: "trunc"}

# =================== GTF: trascritti canonici e loro gene ===================
tx_gene, tx_tags, tx_strand, tx_chr, rows = {}, {}, {}, {}, []
with gzip.open(a.gtf, "rt") as f:
    for line in f:
        if line.startswith("#"): continue
        t = line.rstrip("\n").split("\t")
        if t[2] not in ("transcript", "exon", "CDS", "stop_codon"): continue
        m = re.search(r'transcript_id "([^".]+)', t[8])
        if not m: continue
        tx = m.group(1)
        if t[2] == "transcript":
            g = re.search(r'gene_name "([^"]+)"', t[8]); tx_gene[tx] = g.group(1) if g else tx
            tx_tags[tx] = set(re.findall(r'tag "([^"]+)"', t[8])); tx_strand[tx] = t[6]; tx_chr[tx] = strip(t[0])
        else:
            rows.append((tx, t[2], int(t[3]), int(t[4])))
canon = {tx for tx, tg in tx_tags.items() if "Ensembl_canonical" in tg}
G = pd.DataFrame(rows, columns=["tx", "type", "start", "end"]); del rows
cds_tx = set(G.loc[G["type"] == "CDS", "tx"])
use = {tx for tx in canon if tx in cds_tx}
n0 = len(use); use = {tx for tx in use if not ({"cds_start_NF", "cds_end_NF", "seleno"} & tx_tags[tx])}
log(f"trascritti canonici con CDS: {n0}; esclusi (CDS incompleta o selenoproteina): {n0 - len(use)}")
G = G[G["tx"].isin(use)]

# =================== genoma ===================
def read_fasta(path, chroms):
    seqs, cur, buf = {}, None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if cur is not None and cur in chroms: seqs[cur] = "".join(buf).upper()
                cur = strip(line[1:].split()[0]); buf = []
            elif cur in chroms: buf.append(line.strip())
    if cur is not None and cur in chroms: seqs[cur] = "".join(buf).upper()
    return seqs
GENOME = read_fasta(a.genome, {tx_chr[t] for t in use})
log(f"cromosomi letti dal genoma: {len(GENOME)}; esposizione dell'atteso del modello: {a.exposure}")

# =================== siti per trascritto (stesse regole di compute_coding_sites.R) ===================
B = "ACGT"; COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
CODE = {"".join(c): aa for c, aa in zip(itertools.product("TCAG", repeat=3), "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
CODONS = sorted(CODE); CI = {c: i for i, c in enumerate(CODONS)}
LUT = np.zeros((64, 3, 4, 3), np.int8)                    # [codone, posizione, base alternativa, contesto] -> gravita' 2..5
for c in CODONS:
    for pp in range(3):
        for bi, b in enumerate(B):
            if c[pp] == b: continue
            alt = c[:pp] + b + c[pp + 1:]
            for ctx in range(3):                          # 0 interno, 1 primo codone, 2 codone di stop
                if ctx > 0 or CODE[c] == "*": s = 3
                elif CODE[alt] == "*": s = 5
                elif CODE[alt] == CODE[c]: s = 2
                else: s = 4
                LUT[CI[c], pp, bi, ctx] = s
def tx_sites(tx, d):
    """Ritorna (pos, base_alt_genomica, gravita') delle sostituzioni possibili del trascritto."""
    ch = tx_chr[tx]; seq = GENOME.get(ch)
    if seq is None: return None
    orf = d[d["type"].isin(["CDS", "stop_codon"])].sort_values("start")
    pos = np.unique(np.concatenate([np.arange(s, e + 1) for s, e in zip(orf["start"], orf["end"])]))
    if len(pos) % 3 or len(pos) < 6: return None
    minus = tx_strand[tx] == "-"
    tpos = pos[::-1] if minus else pos
    bases = [seq[q - 1] for q in tpos]
    if minus: bases = [COMP.get(x, "N") for x in bases]
    if any(x not in COMP for x in bases): return None
    L = len(bases); cod = ["".join(bases[i:i + 3]) for i in range(0, L, 3)]
    if any(c not in CI for c in cod): return None
    ncod = L // 3; ci = np.repeat(np.arange(ncod), 3); pin = np.tile(np.arange(3), ncod)
    ctx = np.zeros(ncod, int); ctx[0] = 1; ctx[-1] = 2
    codi = np.array([CI[c] for c in cod])
    P = np.repeat(tpos, 4); Bt = np.tile(np.arange(4), L)              # base alternativa sul filamento del trascritto
    ref = np.repeat(np.array([B.index(x) for x in bases]), 4); keep = Bt != ref
    sev = LUT[codi[ci].repeat(4), pin.repeat(4), Bt, ctx[ci].repeat(4)]
    Bg = (3 - Bt) if minus else Bt                                       # complemento: A<->T (0<->3), C<->G (1<->2)
    out_p, out_b, out_s = [P[keep]], [Bg[keep]], [sev[keep]]
    ex = d[d["type"] == "exon"].sort_values("start")
    for (s1, e1), (s2, e2) in zip(zip(ex["start"], ex["end"]), list(zip(ex["start"], ex["end"]))[1:]):
        i1, i2 = e1 + 1, s2 - 1
        if i2 - i1 + 1 < 4: continue
        for q in (i1, i1 + 1, i2 - 1, i2):
            r = seq[q - 1]
            for bi, b in enumerate(B):
                if b != r: out_p.append(np.array([q])); out_b.append(np.array([bi])); out_s.append(np.array([5], np.int8))
    return np.concatenate(out_p), np.concatenate(out_b), np.concatenate(out_s).astype(np.int8)

# bin
bb = pd.read_csv(a.backbone, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chr", "start", "end", "binID"], dtype={"chr": str, "binID": str})
bb["chr"] = bb["chr"].map(strip); bb["start"] += 1; bb = bb.reset_index(drop=True)
BB = {ch: (b["start"].to_numpy(), b["end"].to_numpy(), b.index.to_numpy()) for ch, b in bb.sort_values(["chr", "start"]).groupby("chr")}
def to_bin(ch, pos):
    if ch not in BB: return np.full(len(pos), -1)
    s, e, ix = BB[ch]; i = np.searchsorted(s, pos, side="right") - 1
    ok = (i >= 0) & (pos <= e[np.clip(i, 0, None)]); return np.where(ok, ix[np.clip(i, 0, None)], -1)

gene_sites = []                                           # (gene, bin, classe, n siti)
CDS_KB = {}                                               # gene -> kb dell'ORF canonico (CDS + stop), per l'atteso "solo lunghezza"
gene_cds = []                                             # (gene, bin, bp dell'ORF)
union = {}                                                # cromosoma -> chiave (pos*4+base) -> gravita' massima
GK = {}                                                   # gene -> [(cromosoma, chiavi, gravita')] dei suoi trascritti canonici
n_bad = 0
for tx, d in G.groupby("tx"):
    r = tx_sites(tx, d)
    if r is None: n_bad += 1; continue
    P, Bg, S = r; ch = tx_chr[tx]
    orf_ = d[d["type"].isin(["CDS", "stop_codon"])]
    orf_pos = np.unique(np.concatenate([np.arange(x, y + 1) for x, y in zip(orf_["start"], orf_["end"])]))
    CDS_KB[tx_gene[tx]] = max(CDS_KB.get(tx_gene[tx], 0), len(orf_pos) / 1000)
    ob = to_bin(ch, orf_pos); ob = ob[ob >= 0]                       # bp dell'ORF del gene per bin (esposizione "cds")
    if len(ob):
        bc_ = np.bincount(ob); gene_cds += [(tx_gene[tx], int(b_), int(bc_[b_])) for b_ in np.nonzero(bc_)[0]]
    key = P.astype(np.int64) * 4 + Bg
    k_u, inv = np.unique(key, return_inverse=True); s_u = np.zeros(len(k_u), np.int8); np.maximum.at(s_u, inv, S)   # dentro il gene
    union.setdefault(ch, []).append((k_u, s_u))
    GK.setdefault(tx_gene[tx], []).append((ch, k_u, s_u))        # siti del gene (per la classe degli eventi per gene)
    bins = to_bin(ch, k_u // 4)
    for c, sv in (("syn", 2), ("mis", 4), ("trunc", 5)):
        m = (s_u == sv) & (bins >= 0)
        if m.any():
            bc = np.bincount(bins[m]); nz = np.nonzero(bc)[0]
            gene_sites += [(tx_gene[tx], int(b_), c, int(bc[b_])) for b_ in nz]
log(f"trascritti esclusi (ORF non multipla di 3, basi non ACGT, cromosoma assente): {n_bad}")
_ex = G[G["type"] == "exon"].assign(bp=lambda d_: d_["end"] - d_["start"] + 1)
GENE_STRUCT = pd.DataFrame({"n_exons": _ex.groupby("tx").size(), "exon_kb": _ex.groupby("tx")["bp"].sum() / 1000})
GENE_STRUCT["gene"] = GENE_STRUCT.index.map(tx_gene); GENE_STRUCT = GENE_STRUCT.groupby("gene").max()
GS = pd.DataFrame(gene_sites, columns=["gene", "bin", "cls", "sites"]); del gene_sites
GC = pd.DataFrame(gene_cds, columns=["gene", "bin", "bp"]); del gene_cds
# siti per bin come nella pipeline: unione dei trascritti, classe piu' grave, una volta per sostituzione
bin_sites = np.zeros((len(bb), 3), np.int64); UNION = {}
for ch, parts in union.items():
    k = np.concatenate([x[0] for x in parts]); s = np.concatenate([x[1] for x in parts])
    o = np.lexsort((-s, k)); k, s = k[o], s[o]; first = np.r_[True, k[1:] != k[:-1]]
    k, s = k[first], s[first]; bins = to_bin(ch, k // 4); UNION[ch] = (k, s)
    for j, sv in enumerate((2, 4, 5)):
        m = (s == sv) & (bins >= 0); bin_sites[:, j] += np.bincount(bins[m], minlength=len(bb))
del union

# =================== eventi (stesse regole di compute_coding_mutations_counts.R) ===================
TRUNC = {"stop_gained", "splice_acceptor_variant", "splice_donor_variant"}; MIS = {"missense_variant"}
OTHER = {"start_lost", "stop_lost", "stop_retained_variant", "coding_sequence_variant", "protein_altering_variant",
         "incomplete_terminal_codon_variant", "start_retained_variant"}
SYN = {"synonymous_variant"}; INDEL = {"frameshift_variant", "inframe_insertion", "inframe_deletion"}
def severity(desc):
    t = set(str(desc).split(","))
    return 5 if t & TRUNC else 4 if t & MIS else 3 if t & OTHER else 2 if t & SYN else 1
cols = ["COSMIC_SAMPLE_ID", "GENOMIC_MUTATION_ID", "CHROMOSOME", "GENOME_START", "GENOME_STOP", "TRANSCRIPT_ACCESSION",
        "MUTATION_DESCRIPTION", "GENOMIC_WT_ALLELE", "GENOMIC_MUT_ALLELE"]
R = pd.read_csv(a.events, sep="\t", dtype=str, quoting=3, usecols=cols, keep_default_na=False)
R["mut_id"] = np.where(R["GENOMIC_MUTATION_ID"] == "", R["CHROMOSOME"] + ":" + R["GENOME_START"] + ":" + R["GENOME_STOP"], R["GENOMIC_MUTATION_ID"])
R["pos"] = pd.to_numeric(R["GENOME_START"], errors="coerce"); R["stop"] = pd.to_numeric(R["GENOME_STOP"], errors="coerce")
R = R[R["pos"].notna() & (R["COSMIC_SAMPLE_ID"] != "")]
ud = pd.Series(R["MUTATION_DESCRIPTION"].unique())
R["sev"] = R["MUTATION_DESCRIPTION"].map(dict(zip(ud, ud.map(severity))))
R["coding_any"] = R["MUTATION_DESCRIPTION"].map(dict(zip(ud, ud.map(lambda d: bool(set(str(d).split(",")) & (TRUNC | MIS | OTHER | SYN | INDEL))))))
R["snv"] = (R["pos"] == R["stop"]) & R["GENOMIC_WT_ALLELE"].str.upper().isin(set(B)) & R["GENOMIC_MUT_ALLELE"].str.upper().isin(set(B))
R["tx"] = R["TRANSCRIPT_ACCESSION"].str.split(".").str[0]; R["canon"] = R["tx"].isin(canon)
per_sample = R.groupby(["COSMIC_SAMPLE_ID", "mut_id"])["coding_any"].any().groupby(level=0).sum()   # eventi codificanti per campione
if str(a.hypermutator_max).lower() not in ("none", ""):
    keep_s = set(per_sample[per_sample <= float(a.hypermutator_max)].index)
    log(f"campioni: {len(per_sample)}, esclusi come ipermutati: {len(per_sample) - len(keep_s)}")
    R = R[R["COSMIC_SAMPLE_ID"].isin(keep_s)]
# FIX 28/10/2026 - RIANNOTAZIONE (come compute_coding_mutations_counts.R con event_annotation = reference): la classe di una
# sostituzione singola e' quella del SITO (posizione + base mutata) nei trascritti canonici, dalla sequenza di riferimento,
# non quella delle righe COSMIC (COSMIC per alcuni geni non riporta mai il canonico Ensembl: AKT1, CASP8, RUNX1...).
# Base di riferimento COSMIC diversa dal genoma: evento escluso dalle classi (numero nel log).
ev = R.groupby(["COSMIC_SAMPLE_ID", "mut_id"], sort=False).agg(chr=("CHROMOSOME", "first"), pos=("pos", "first"), snv=("snv", "first"),
                                                            wt=("GENOMIC_WT_ALLELE", "first"), mt=("GENOMIC_MUT_ALLELE", "first")).reset_index()
ev["chr"] = ev["chr"].map(strip); ev["wt"] = ev["wt"].str.upper(); ev["mt"] = ev["mt"].str.upper()
ev_all = ev[["COSMIC_SAMPLE_ID", "mut_id", "chr", "pos", "snv"]].copy()          # tutti gli eventi dei campioni tenuti (per il tetto)
ev = ev[ev["snv"]].copy(); ev["pos"] = ev["pos"].astype(np.int64)
ev["key"] = ev["pos"] * 4 + ev["mt"].map({b: i for i, b in enumerate(B)}).astype(np.int64)
ev["ref_ok"] = [GENOME[c][p - 1] == w if c in GENOME and p <= len(GENOME[c]) else False for c, p, w in zip(ev["chr"], ev["pos"], ev["wt"])]
log(f"riannotazione: SNV con base di riferimento COSMIC diversa dal genoma (esclusi dalle classi): {int((~ev['ref_ok']).sum())}")
def lookup(keys, k, sv):                                       # gravita' del sito per ogni chiave (0 = non e' un sito)
    i = np.searchsorted(k, keys); i = np.clip(i, 0, max(len(k) - 1, 0))
    return np.where((len(k) > 0) & (k[i] == keys), sv[i], 0) if len(k) else np.zeros(len(keys), int)
# (a) classe per BIN (come la pipeline): sito nell'unione dei trascritti canonici, classe piu' grave
ev["sev"] = 0
for ch, idx in ev.groupby("chr").groups.items():
    if ch in UNION: ev.loc[idx, "sev"] = lookup(ev.loc[idx, "key"].to_numpy(), *UNION[ch])
ev.loc[~ev["ref_ok"], "sev"] = 0
ev["bin"] = -1
for ch, idx in ev.groupby("chr").groups.items(): ev.loc[idx, "bin"] = to_bin(ch, ev.loc[idx, "pos"].to_numpy(np.int64))
# eventi classificati, per noise_ceiling_mut_classes.py --events_classified (stessa classificazione, niente codice duplicato)
ev_all = ev_all.merge(ev[["COSMIC_SAMPLE_ID", "mut_id", "sev"]].assign(cls=lambda d_: d_["sev"].map(SEV2CLS)).drop(columns="sev"),
                      on=["COSMIC_SAMPLE_ID", "mut_id"], how="left")
ev_all.to_csv(f"{a.out_dir}/events_classified.tsv", sep="\t", index=False)
evb = ev[ev["bin"] >= 0]
bin_events = np.zeros((len(bb), 3), np.int64)
for j, sv in enumerate((2, 4, 5)): bin_events[:, j] = np.bincount(evb.loc[evb["sev"] == sv, "bin"], minlength=len(bb))
# (b) classe per GENE: sito nel trascritto canonico DEL gene (un evento in due geni sovrapposti conta in ciascuno)
evg = ev[ev["ref_ok"]]; byc = {ch: d_ for ch, d_ in evg.groupby("chr")}; eg_rows = []
for gene, parts in GK.items():
    for ch, k_u, s_u in parts:
        if ch not in byc: continue
        d_ = byc[ch]; sv = lookup(d_["key"].to_numpy(), k_u, s_u); m = sv > 0
        if m.any(): eg_rows.append(pd.DataFrame({"COSMIC_SAMPLE_ID": d_["COSMIC_SAMPLE_ID"].to_numpy()[m], "mut_id": d_["mut_id"].to_numpy()[m],
                                                 "gene": gene, "sev": sv[m]}))
eg = pd.concat(eg_rows).groupby(["COSMIC_SAMPLE_ID", "mut_id", "gene"])["sev"].max().reset_index() if eg_rows else \
     pd.DataFrame(columns=["COSMIC_SAMPLE_ID", "mut_id", "gene", "sev"])
O = pd.DataFrame({c: eg[eg["sev"] == sv].groupby("gene").size() for c, sv in (("syn", 2), ("mis", 4), ("trunc", 5))}).fillna(0)
# eventi NON codificanti del gene (introni vicini agli esoni, UTR), per l'informazione locale "+nc": restano definiti dalle
# righe COSMIC sul canonico (non ci sono siti per loro); per i geni senza il canonico in COSMIC sono 0 (limite noto).
Rc = R[R["canon"]]; Rg = Rc[Rc["snv"] & Rc["tx"].isin(use)].copy(); Rg["gene"] = Rg["tx"].map(tx_gene)
eg_nc = Rg.groupby(["COSMIC_SAMPLE_ID", "mut_id", "gene"])["sev"].max().reset_index()
O_NC = eg_nc[eg_nc["sev"] == 1].groupby("gene").size()
del GENOME
del R, Rc, Rg

# =================== controlli di coerenza con la pipeline ===================
D = pd.read_csv(a.dataset, sep="\t", dtype={"binID": str}).set_index("binID")
bins_id = bb["binID"].to_numpy(); common = np.isin(bins_id, D.index)
coh = []
for j, c in enumerate(CLASSES):
    for what, mine, col in (("eventi", bin_events[:, j], f"mut_n_{c}"), ("siti", bin_sites[:, j], f"sites_{c}")):
        pipe = D[col].reindex(bins_id).fillna(0).to_numpy(np.int64)
        nd = int((mine[common] != pipe[common]).sum())
        coh.append(dict(what=what, cls=c, bins_diff=nd, total_here=int(mine[common].sum()), total_pipeline=int(pipe[common].sum())))
        log(f"COERENZA {what:6s} {c:5s}: bin diversi {nd} su {int(common.sum())}; totale qui {int(mine[common].sum())}, pipeline {int(pipe[common].sum())}")
pd.DataFrame(coh).to_csv(f"{a.out_dir}/coherence_with_pipeline.tsv", sep="\t", index=False)
if any(r["bins_diff"] for r in coh): log("ATTENZIONE: differenze con la pipeline: i risultati per gene NON vanno usati prima di capirle")

# =================== atteso per gene ===================
WB = WG = None
if a.exposure == "wsites" and not (a.wsites_bin and a.wsites_gene): sys.exit("--exposure wsites richiede --wsites_bin e --wsites_gene")
if a.wsites_gene:
    WG = pd.read_csv(a.wsites_gene, sep="\t", dtype={"binID": str})
    WG["bin"] = pd.Series(np.arange(len(bb)), index=bb["binID"].to_numpy()).reindex(WG["binID"]).to_numpy()
    # controllo: i siti NON pesati per gene e classe devono coincidere con quelli contati qui
    chk = WG.groupby(["gene", "cls"])["sites"].sum(); mine = GS.groupby(["gene", "cls"])["sites"].sum()
    j = chk.index.intersection(mine.index)
    log(f"COERENZA siti per gene (file pesato vs conteggio qui): coppie gene-classe diverse {int((chk[j] != mine[j]).sum())} su {len(j)}")
if a.wsites_bin: WB = pd.read_csv(a.wsites_bin, sep="\t", dtype={"binID": str})
oof = pd.read_csv(a.oof, sep="\t", dtype={"binID": str})
bpos = pd.Series(np.arange(len(bb)), index=bins_id)
EXP = {}
for j, c in enumerate(CLASSES):
    o = oof[oof["omic"] == f"mut_{c}"]
    if o.empty: sys.exit(f"nelle predizioni manca l'omica mut_{c}")
    rate = np.full(len(bb), np.nan); b_ = bpos.reindex(o["binID"]).to_numpy()
    if a.exposure == "wsites":
        expo_bin = WB.set_index("binID")[f"wsites_{c}"]
    else:
        expo_bin = D[f"sites_{c}"] if a.exposure == "sites" else D["cds_bp"]       # stessa esposizione dell'offset del pilota
    rate[b_] = o["pred"].to_numpy(float) / expo_bin.reindex(o["binID"]).to_numpy(float)
    gs_c = GS[GS["cls"] == c]
    if a.exposure == "wsites":
        g = WG[WG["cls"] == c][["gene", "bin", "wsites"]].rename(columns={"wsites": "w"}).copy()
    else:
        g = (gs_c.rename(columns={"sites": "w"}) if a.exposure == "sites" else
             GC[GC["gene"].isin(set(gs_c["gene"]))].rename(columns={"bp": "w"})).copy()
    g["rate"] = rate[g["bin"].to_numpy()]
    tot = g.groupby("gene")["w"].sum()
    cov = g[g["rate"].notna()].groupby("gene")["w"].sum().reindex(tot.index).fillna(0) / tot
    e = (g["rate"] * g["w"]).groupby(g["gene"]).sum(min_count=1).reindex(tot.index) / cov
    EXP[c] = pd.DataFrame({"E": e.where(cov >= a.min_cov), "sites": gs_c.groupby("gene")["sites"].sum().reindex(tot.index), "cov": cov})
genes = sorted(set().union(*[set(EXP[c].index) for c in CLASSES]))

# =================== statistica: binomiale negativa ===================
def theta_mom(Ov, Ev):
    """Dispersione con i momenti: Var(O) = E + E^2/theta; escluso l'1% con lo z di Poisson piu' alto (niente liste)."""
    z = (Ov - Ev) / np.sqrt(np.clip(Ev, 1e-12, None)); k = z <= np.nanquantile(z, 0.99)
    inv = (np.sum((Ov[k] - Ev[k]) ** 2) - np.sum(Ov[k])) / np.sum(Ev[k] ** 2)
    return np.inf if inv <= 0 else 1 / inv
def nb_upper_score(Ov, Ev, theta):
    """-log10 P(X >= O), X ~ NB(media E, dispersione theta) (Poisson se theta = inf). Esatta: P(X >= O) = I_{1-p}(O, theta)
    con p = theta / (theta + E); dove va in underflow, somma della pmf in scala log su un intervallo ampio."""
    out = np.zeros(len(Ov))
    for i, (o_, e_) in enumerate(zip(Ov, Ev)):
        if o_ <= 0: out[i] = 0.0; continue
        if np.isinf(theta): lp = stats.poisson.logsf(o_ - 1, e_)
        else:
            pp = theta / (theta + e_); v = special.betainc(o_, theta, 1 - pp)
            lp = np.log(v) if v > 1e-280 else np.nan
            if not np.isfinite(lp):
                k = np.arange(o_, o_ + 20000)
                lp = special.logsumexp(stats.nbinom.logpmf(k, theta, pp))
        if not np.isfinite(lp):                                   # Poisson con underflow: somma in scala log
            k = np.arange(o_, o_ + 20000); lp = special.logsumexp(stats.poisson.logpmf(k, e_))
        out[i] = -min(lp, 0.0) / np.log(10)
    return out

SC = {}
def score_block(name, Ov, Ev):
    ok = np.isfinite(Ev); th = theta_mom(Ov[ok], Ev[ok])
    res = np.full(len(Ov), np.nan); res[ok] = nb_upper_score(Ov[ok], Ev[ok], th)
    log(f"{name}: geni con atteso {int(ok.sum())} su {len(Ov)}; dispersione theta = {th:.3g}")
    return res, th
for c in CLASSES:
    Ov = O[c].reindex(genes).fillna(0).to_numpy(float); Ev = EXP[c]["E"].reindex(genes).to_numpy(float)
    r_, th = score_block(c, Ov, Ev); SC[c] = pd.DataFrame({"O": Ov, "E": Ev, "residual": r_}, index=genes); SC[c].attrs["theta"] = th
Ov = SC["mis"]["O"] + SC["trunc"]["O"]; Ev = SC["mis"]["E"] + SC["trunc"]["E"]
r_, th = score_block("nonsyn", Ov.to_numpy(float), Ev.to_numpy(float))
SC["nonsyn"] = pd.DataFrame({"O": Ov, "E": Ev, "residual": r_}, index=genes)

# =================== valutazione ===================
def read_list(path, names):
    t = pd.read_csv(path, sep="\t", header=None, dtype=str, comment="#"); best, fr = None, 0
    for col in t.columns[3:]:
        f_ = t[col].isin(names).mean()
        if f_ > fr: best, fr = col, f_
    return set(t[best][t[best].isin(names)]) if best is not None else set()
gene_names = set(tx_gene.values()); drivers = read_list(a.og, gene_names) | read_list(a.tsg, gene_names)
log(f"driver IntOGen (OG + TSG): {len(drivers)}; nell'universo dei geni valutati: {len(drivers & set(genes))}")
def order(s): return np.lexsort((TIE.random(len(s)), -s))
def pr_steps(s, l):
    o = np.argsort(-s, kind="mergesort"); ss, lo = s[o], l[o]
    last = np.r_[np.where(np.diff(ss) != 0)[0], len(ss) - 1]; tp = np.cumsum(lo)[last]
    return tp / max(l.sum(), 1), tp / (last + 1)
def ap(s, l):
    if not l.any(): return np.nan
    rec, prec = pr_steps(s, l); return float(np.sum(np.diff(np.r_[0, rec]) * prec))
def auroc(s, l):
    r = pd.Series(s).rank().to_numpy(); n1 = l.sum(); n0 = len(l) - n1
    return (r[l].sum() - n1 * (n1 + 1) / 2) / (n1 * n0) if n1 and n0 else np.nan
def bottom(s):                                                # senza punteggio = in fondo (conta come non trovato)
    s = np.asarray(s, float); return np.where(np.isfinite(s), s, np.nanmin(s[np.isfinite(s)]) - 1 if np.isfinite(s).any() else 0)
rows, curves = [], {}
for name, df in SC.items():
    l = df.index.isin(drivers)
    for sc in ("observed", "residual"):
        s = bottom(df["O"] if sc == "observed" else df["residual"])
        bs = []
        pos_, neg_ = np.where(l)[0], np.where(~l)[0]
        for _ in range(a.n_boot):
            i = np.r_[rng.choice(pos_, len(pos_)), rng.choice(neg_, len(neg_))]; bs.append((auroc(s[i], l[i]), ap(s[i], l[i])))
        bs = np.array(bs); top = order(s)
        rows.append(dict(target=name, score=sc, n_genes=len(s), n_drivers=int(l.sum()), prevalence=l.mean(), AUROC=auroc(s, l),
                         AUROC_lo=np.nanpercentile(bs[:, 0], 2.5), AUROC_hi=np.nanpercentile(bs[:, 0], 97.5), AP=ap(s, l),
                         AP_lo=np.nanpercentile(bs[:, 1], 2.5), AP_hi=np.nanpercentile(bs[:, 1], 97.5),
                         drivers_top20=int(l[top[:20]].sum()), drivers_top100=int(l[top[:100]].sum())))
        curves[(name, sc)] = pr_steps(s, l)
    if name in EXP: df = df.assign(sites=EXP[name]["sites"].reindex(df.index))          # siti della classe nel gene (NOTE 28/10/2026)
    out = df.assign(is_driver=l).sort_values("residual", ascending=False)
    out.to_csv(f"{a.out_dir}/gene_scores_{name}.tsv", sep="\t")
    nt = min(a.top_n, len(df))
    tt = pd.DataFrame({"rank": range(1, nt + 1),
                       "top_observed": [g + ("*" if g in drivers else "") for g in df.index[order(bottom(df["O"]))][:nt]],
                       "top_residual": [g + ("*" if g in drivers else "") for g in df.index[order(bottom(df["residual"]))][:nt]]})
    tt.to_csv(f"{a.out_dir}/top_genes_{name}.tsv", sep="\t", index=False)
M = pd.DataFrame(rows); M.to_csv(f"{a.out_dir}/gene_level_metrics.tsv", sep="\t", index=False)
pd.set_option("display.width", 220); log(M.round(3).to_string(index=False))

# =================== CONFRONTO TRA ATTESI: quanto fa l'esposizione e quanto le covariate ===================
# Stessa statistica (coda della binomiale negativa, theta stimato per ciascun atteso) e stesso universo (i geni con
# l'atteso del modello); cambia SOLO l'atteso E:
#   costante  E uguale per tutti           -> classifica = conteggio osservato
#   lunghezza tasso medio per kb x kb ORF   -> solo la lunghezza
#   siti      tasso medio per sito x siti   -> esposizione completa (lunghezza + composizione dei codoni) = "solo offset"
#   modello   tasso del bin x siti          -> esposizione + covariate
# Misure: AP, AUROC, driver nei primi 20/100, theta nei NON driver (piu' alto = meno variabilita' non spiegata),
# devianza di Poisson nei non driver e devianza spiegata rispetto all'atteso "siti", Spearman con la lunghezza nei non
# driver, sovrapposizione dei primi 30 con la classifica del modello.
def pois_dev(Ov, Ev):
    Ev = np.clip(Ev, 1e-12, None); t = np.where(Ov > 0, Ov * np.log(np.where(Ov > 0, Ov, 1) / Ev), 0.0)
    return 2 * np.sum(t - (Ov - Ev))
cmp_rows, cmp_tops = [], {}
VAR = {}                                                           # punteggi per gene di ogni atteso: VAR[classe][atteso] (NaN = senza punteggio)
for name in ["syn", "mis", "trunc", "nonsyn"]:
    df = SC[name]; ok = np.isfinite(df["E"].to_numpy(float))
    sub = df[ok].copy(); l = sub.index.isin(drivers); nd = ~l
    Ov = sub["O"].to_numpy(float)
    if name == "nonsyn":
        sites_g = EXP["mis"]["sites"].reindex(sub.index).fillna(0) + EXP["trunc"]["sites"].reindex(sub.index).fillna(0)
    else:
        sites_g = EXP[name]["sites"].reindex(sub.index).fillna(0)
    kb = pd.Series(CDS_KB).reindex(sub.index).to_numpy(float)
    Es = {"constant": np.full(len(sub), Ov.sum() / len(sub)),
          "length": Ov.sum() / np.nansum(kb) * kb,
          "sites": Ov.sum() / sites_g.sum() * sites_g.to_numpy(float),
          "model": sub["E"].to_numpy(float)}
    if WG is not None:                                            # solo contesto: tasso medio per sito pesato x siti pesati del gene
        cl_ = ["mis", "trunc"] if name == "nonsyn" else [name]
        wg_ = WG[WG["cls"].isin(cl_)].groupby("gene")["wsites"].sum().reindex(sub.index).fillna(0).to_numpy(float)
        Es = {**{k: v for k, v in Es.items() if k != "model"}, "wsites": Ov.sum() / wg_.sum() * wg_, "model": Es["model"]}
    sc_model = None; tops = {}
    for en, Ev in Es.items():
        Ev = np.where(np.isfinite(Ev) & (Ev > 0), Ev, np.nan); good = np.isfinite(Ev)
        th = theta_mom(Ov[good & nd], Ev[good & nd])
        scr = np.full(len(sub), np.nan); scr[good] = nb_upper_score(Ov[good], Ev[good], th)
        VAR.setdefault(name, {})[en] = pd.Series(scr.copy(), index=sub.index); scr = bottom(scr)
        top = order(scr); tops[en] = list(sub.index[top[:30]])
        rho = pd.Series(scr[nd]).corr(pd.Series(kb[nd]), method="spearman")
        cmp_rows.append(dict(target=name, expectation=en, n_genes=len(sub), n_drivers=int(l.sum()), AP=ap(scr, l), AUROC=auroc(scr, l),
                             drivers_top20=int(l[top[:20]].sum()), drivers_top100=int(l[top[:100]].sum()),
                             theta_nondrivers=th, poisson_dev_nondrivers=pois_dev(Ov[good & nd], Ev[good & nd]),
                             spearman_length_nondrivers=rho))
    base = [r for r in cmp_rows if r["target"] == name and r["expectation"] == "sites"][0]["poisson_dev_nondrivers"]
    for r in cmp_rows:
        if r["target"] == name:
            r["dev_explained_vs_sites"] = 1 - r["poisson_dev_nondrivers"] / base
            r["top30_overlap_with_model"] = len(set(tops[r["expectation"]]) & set(tops["model"]))
    cmp_tops[name] = pd.DataFrame({en: [g + ("*" if g in drivers else "") for g in t] for en, t in tops.items()},
                                  index=range(1, len(tops["model"]) + 1))
    cmp_tops[name].to_csv(f"{a.out_dir}/expectations_top30_{name}.tsv", sep="\t", index_label="rank")
# ---- informazione LOCALE (specifica delle mutazioni), come fanno gli strumenti dedicati ----
# Per ogni gene, un fattore moltiplicativo m_g stimato da eventi NEUTRALI dello stesso gene, con contrazione verso 1
# (gamma-Poisson): m_g = (O_neutri + theta) / (E_neutri + theta), dove E_neutri e' l'atteso dei neutri dal modello di
# riferimento e theta la loro dispersione nei non driver. Pochi eventi -> m vicino a 1; molti -> m segue i dati.
#   +syn : neutri = sinonime del gene (come dNdScv); atteso = atteso delle sinonime (stessa base: model o length)
#   +nc  : neutri = eventi in introni vicini agli esoni e UTR del gene (come la "regional non-coding rate" di
#          Lawrence et al.); atteso = regressione di Poisson su n. di esoni, kb di esoni e kb di ORF (non c'e' un
#          modello regionale per i non codificanti)
# Le sinonime non sono valutate con +syn (sarebbe circolare).
from sklearn.linear_model import PoissonRegressor
def shrink(On, En, nd):
    ok = np.isfinite(En) & (En > 0); th = theta_mom(On[ok & nd], En[ok & nd])
    th = 1e6 if not np.isfinite(th) else th
    return np.where(ok, (On + th) / (En + th), 1.0), th
genes_all = SC["syn"].index
On_syn = SC["syn"]["O"].to_numpy(float)
kb_all = pd.Series(CDS_KB).reindex(genes_all).to_numpy(float)
st = GENE_STRUCT.reindex(genes_all)
Xnc = np.c_[np.log1p(st["n_exons"].fillna(1)), np.log1p(st["exon_kb"].fillna(0)), np.log1p(np.nan_to_num(kb_all))]
On_nc = O_NC.reindex(genes_all).fillna(0).to_numpy(float)
nd_all = ~genes_all.isin(drivers)
pr = PoissonRegressor(alpha=0, max_iter=1000).fit(Xnc[nd_all], On_nc[nd_all])
En_nc = pr.predict(Xnc)
m_nc, th_nc = shrink(On_nc, En_nc, nd_all)
log(f"\nINFORMAZIONE LOCALE: eventi non codificanti nei geni {int(On_nc.sum())}; theta non codificanti {th_nc:.3g}")
base_syn = {"model": SC["syn"]["E"].to_numpy(float),
            "length": On_syn.sum() / np.nansum(kb_all) * kb_all}
if WG is not None:
    ws_syn = WG[WG["cls"] == "syn"].groupby("gene")["wsites"].sum().reindex(genes_all).fillna(0).to_numpy(float)
    base_syn["wsites"] = On_syn.sum() / ws_syn.sum() * ws_syn
m_syn = {}
for bn, Eb in base_syn.items():
    m_syn[bn], th_s = shrink(On_syn, Eb, nd_all); log(f"theta sinonime (base {bn}): {th_s:.3g}")
loc_rows = []
for name in ["mis", "trunc", "nonsyn"]:
    df = SC[name]; ok = np.isfinite(df["E"].to_numpy(float)); sub = df[ok]
    idx = genes_all.get_indexer(sub.index); l = sub.index.isin(drivers); nd = ~l
    Ov = sub["O"].to_numpy(float); kb = kb_all[idx]
    bases = {"model": sub["E"].to_numpy(float), "length": Ov.sum() / np.nansum(kb) * kb}
    if WG is not None:
        cl_ = ["mis", "trunc"] if name == "nonsyn" else [name]
        wg_ = WG[WG["cls"].isin(cl_)].groupby("gene")["wsites"].sum().reindex(sub.index).fillna(0).to_numpy(float)
        bases["wsites"] = Ov.sum() / wg_.sum() * wg_
    for bn, Eb in bases.items():
        for loc, mult in (("", np.ones(len(sub))), ("+syn", m_syn[bn][idx]), ("+nc", m_nc[idx])):
            Ev = Eb * mult; Ev = Ev * Ov.sum() / np.nansum(Ev) if bn in ("length", "wsites") else Ev    # come nel confronto tra attesi
            good = np.isfinite(Ev) & (Ev > 0); th = theta_mom(Ov[good & nd], Ev[good & nd])
            scr = np.full(len(sub), np.nan); scr[good] = nb_upper_score(Ov[good], Ev[good], th)
            if loc: VAR.setdefault(name, {})[bn + loc] = pd.Series(scr.copy(), index=sub.index)
            scr = bottom(scr); top = order(scr)
            loc_rows.append(dict(target=name, expectation=bn + loc, AP=ap(scr, l), AUROC=auroc(scr, l),
                                 drivers_top20=int(l[top[:20]].sum()), drivers_top100=int(l[top[:100]].sum()),
                                 theta_nondrivers=th, poisson_dev_nondrivers=pois_dev(Ov[good & nd], Ev[good & nd])))
LOC = pd.DataFrame(loc_rows)
for name in LOC.target.unique():
    b = LOC[(LOC.target == name) & (LOC.expectation == "length")]["poisson_dev_nondrivers"].iloc[0]
    LOC.loc[LOC.target == name, "dev_explained_vs_length"] = 1 - LOC.loc[LOC.target == name, "poisson_dev_nondrivers"] / b
LOC.to_csv(f"{a.out_dir}/local_information_comparison.tsv", sep="\t", index=False)
log("\nCONFRONTO CON INFORMAZIONE LOCALE (stesso universo: geni con l'atteso del modello)")
log(LOC.round(3).to_string(index=False))

CMP = pd.DataFrame(cmp_rows); CMP.to_csv(f"{a.out_dir}/expectations_comparison.tsv", sep="\t", index=False)
# NOTE (28/10/2026): punteggi per gene di TUTTI gli attesi delle due tabelle (constant, length, sites, [wsites], model, e le
# varianti +syn / +nc), per poterli confrontare con dNdScv (compare_with_dndscv.py --ours_variants). Stessi geni dei
# gene_scores_<classe>.tsv; NaN = gene senza punteggio per quell'atteso (le tabelle usano solo i geni con l'atteso del
# modello, quindi anche le varianti senza covariate sono NaN fuori da quell'insieme: stesso svantaggio del modello).
for name, d in VAR.items():
    V = pd.DataFrame(d).reindex(SC[name].index); V.insert(0, "O", SC[name]["O"]); V["is_driver"] = V.index.isin(drivers)
    V.to_csv(f"{a.out_dir}/gene_scores_variants_{name}.tsv", sep="\t")
log("\nCONFRONTO TRA ATTESI (stessa statistica, stesso universo: geni con l'atteso del modello)")
log(CMP.round(3).to_string(index=False))

# figura: AP per bersaglio (ognuno con la sua scala) e curve precision-recall (etichette in inglese)
NAME = {"syn": "Synonymous", "mis": "Missense", "trunc": "Truncating", "nonsyn": "Non-synonymous (mis + trunc)"}
fig, axes = plt.subplots(2, 4, figsize=(16, 7.5))
for j, name in enumerate(["syn", "mis", "trunc", "nonsyn"]):
    r = M[M.target == name].set_index("score"); ax = axes[0, j]
    for k, (sc, col) in enumerate([("observed", "#7f8c8d"), ("residual", "#c0392b")]):
        v = r.loc[sc, "AP"]; ax.bar(k, v, 0.6, color=col, yerr=[[v - r.loc[sc, "AP_lo"]], [r.loc[sc, "AP_hi"] - v]], capsize=3)
        ax.text(k, r.loc[sc, "AP_hi"], f"{v:.3f}", ha="center", va="bottom", fontsize=8)
    prev = r["prevalence"].iloc[0]; ax.axhline(prev, color="k", ls="--", lw=0.8)
    ax.set_ylim(0, max(r["AP_hi"].max(), prev) * 1.35); ax.set_xticks([0, 1]); ax.set_xticklabels(["observed", "residual"])
    ax.set_title(f"{NAME[name]}\n{int(r['n_drivers'].iloc[0])} driver genes", fontsize=9)
    if j == 0: ax.set_ylabel("average precision (95% CI)")
    ax = axes[1, j]
    for sc, col in [("observed", "#7f8c8d"), ("residual", "#c0392b")]:
        rec, prec = curves[(name, sc)]; ax.step(rec, prec, where="post", color=col, lw=1.3, label=f"{sc} (AP = {r.loc[sc, 'AP']:.3f})")
    ax.axhline(prev, color="k", ls="--", lw=0.8, label=f"random ({prev:.3f})"); ax.set_ylim(0, 1); ax.legend(fontsize=7)
    ax.set_xlabel("recall (driver genes found)")
    if j == 0: ax.set_ylabel("precision")
fig.suptitle("Known drivers at gene level, new mutation target: observed events vs residual (negative binomial tail)\n"
             "Synonymous mutations are a negative control (no selection expected)", fontsize=10)
fig.tight_layout()
for ext in ("png", "pdf"): fig.savefig(f"{a.out_dir}/fig_gene_level_classes.{ext}", dpi=200, bbox_inches="tight")
open(f"{a.out_dir}/log.txt", "w").write("\n".join(LOG) + "\n")
print(f"Output in {a.out_dir}/", file=sys.stderr)