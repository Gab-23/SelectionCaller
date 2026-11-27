#!/usr/bin/env python
"""
BURNER (27/11/2026) - gene-structure map of coding mutations (step 1 of Gab's idea: pool the sparse selection signal of
all genes on ONE shared map of gene structure). This script only BUILDS the map and DESCRIBES it: no model, no
background, no reduced representation, no pipeline rule.

MAP (the same 32 cells for every gene; built from the Ensembl_canonical transcript of each gene; no external data)
  Exon categories are defined on the TRANSCRIPT exons (mRNA, UTR exons included), as the NMD rule is:
    first     first exon                                               -> 10 bins of relative position
    internal  exons 2 .. n-1, outside the pen55 window                 -> 10 bins (all internal exons together)
    pen55     last 55 nt of the mRNA before the LAST exon-exon junction -> 1 bin
              (NMD boundary; when the penultimate exon is shorter than 55 nt the window extends into earlier exons)
    last      last exon                                                -> 10 bins
    splice    essential splice sites (+-1/+-2 intronic bases of every junction, as compute_coding_sites.R)
                                                                       -> 1 bin, truncating sites only
  Single-exon gene: every coding site is "last". Two-exon gene: last 55 nt of the first exon -> pen55, rest -> first.
  Relative bin: inside each gene and category, the ORF nucleotides (CDS + stop codon) are ordered 5'->3' and their rank
  is cut into 10 equal parts. A category that is only UTR in a gene has no coding sites in that gene.
  Classes: syn / mis / trunc (first codon and stop codon are "other coding", not modelled, as in the pipeline).

SITES and EVENTS per gene: same rules as gene_level_drivers_classes.py
  - sites = possible single-base substitutions of the gene's canonical transcript, classified on the reference sequence
    with the rules of compute_coding_sites.R (transcripts with cds_start_NF / cds_end_NF / seleno excluded);
  - event = (sample, mutation), single-base substitutions only; class = class of the SITE (position + mutated base) in the
    gene's canonical transcript; COSMIC reference base must match the genome; samples with more than
    --hypermutator_max coding events (substitutions and indels) are excluded;
  - an event in two overlapping genes counts in each gene (with its class in that gene).
  In the MAP both sites and events are restricted to positions inside the bins of --backbone (the 20 kb backbone after
  the blacklist filter), so the map uses the same genome territory as the 20 kb dataset.

OUTPUTS (--out_dir)
  gene_map_long.tsv         gene, cell, cls, sites, events (one row per gene x cell x class with sites > 0)
  cell_summary.tsv          per cell and class: genes with sites, sites, events, genes with >= 1 event, events per
                            million sites, and ratio_vs_syn = (events / sites of the class) / (events / sites of syn) in the
                            SAME cell, with a 95% Wald interval on the log scale (splice cell: vs syn of all exonic cells
                            pooled). DESCRIPTIVE dN/dS-like ratio: not corrected for trinucleotide context, no background.
  transcript_structure.tsv  per gene: number of exons, ORF nucleotides per category, flags (pen55 window extends beyond the
                            penultimate exon; last exon without coding sites)
  check_vs_gene_scores.tsv  regression check against gene_scores_<cls>.tsv of gene_level_drivers_classes.py (if
                            --gene_scores_dir is given): sites summed over the map cells must equal column 'sites';
                            events WITHOUT the backbone restriction must equal column 'O' (that script does not restrict
                            the observed events to the bins, while its sites are restricted).
  log.txt

Usage (repository root on the cluster, pilot env; ~30 GB of memory for the genome):
  python gene_structure_map.py \
    --gtf resources/reference_data/gencode/gencode.v45.annotation.gtf.gz \
    --genome resources/reference_data/genome/hg38.fa \
    --backbone resources/backbone/hg38_binned_20000_0.mask.filtered.sorted.bed \
    --events resources/COSMIC/mutations/coding/f_breast_carcinoma.filtered.bed \
    --gene_scores_dir analyses/v4/BRCA/20000bp/mut_neutral_training/A/drivers_gene_level \
    --out_dir analyses/v4/BRCA/20000bp/gene_structure_map
"""
import argparse, os, sys, gzip, re, itertools
import numpy as np, pandas as pd

p = argparse.ArgumentParser()
for k in ("gtf", "genome", "backbone", "events"):
    p.add_argument(f"--{k}", required=True)
p.add_argument("--gene_scores_dir", default=None, help="folder with gene_scores_<cls>.tsv for the regression check")
p.add_argument("--hypermutator_max", default="3000")
p.add_argument("--n_bins", type=int, default=10, help="relative-position bins in first / internal / last")
p.add_argument("--pen_nt", type=int, default=55, help="mRNA nucleotides before the last exon-exon junction (NMD window)")
p.add_argument("--out_dir", default="gene_structure_map")
a = p.parse_args()
os.makedirs(a.out_dir, exist_ok=True)

LOG = []
def log(m):
    print(m, file=sys.stderr)
    LOG.append(str(m))
strip = lambda x: re.sub("^chr", "", str(x))

B = "ACGT"
BI = {b: i for i, b in enumerate(B)}
CLASS_NAMES = ["syn", "mis", "trunc"]
SEV_TO_CLASS = np.full(8, -1, np.int64)                 # severity (2 syn, 3 other, 4 mis, 5 trunc) -> class index
SEV_TO_CLASS[2], SEV_TO_CLASS[4], SEV_TO_CLASS[5] = 0, 1, 2

# ---------------- map cells ----------------
NB = a.n_bins
FIRST0, INTERNAL0, PEN, LAST0 = 0, NB, 2 * NB, 2 * NB + 1
SPLICE = 3 * NB + 1
N_CELLS = 3 * NB + 2
CELL_NAMES = ([f"first_{i + 1:02d}" for i in range(NB)] + [f"internal_{i + 1:02d}" for i in range(NB)] + ["pen55"] +
              [f"last_{i + 1:02d}" for i in range(NB)] + ["splice"])
CAT_FIRST, CAT_INTERNAL, CAT_PEN, CAT_LAST = 0, 1, 2, 3
CAT_NAMES = ["first", "internal", "pen55", "last"]

# ---------------- GTF: canonical transcripts and their gene (as gene_level_drivers_classes.py) ----------------
tx_gene, tx_tags, tx_strand, tx_chr, rows = {}, {}, {}, {}, []
with gzip.open(a.gtf, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        t = line.rstrip("\n").split("\t")
        if t[2] not in ("transcript", "exon", "CDS", "stop_codon"):
            continue
        m = re.search(r'transcript_id "([^".]+)', t[8])
        if not m:
            continue
        tx = m.group(1)
        if t[2] == "transcript":
            g = re.search(r'gene_name "([^"]+)"', t[8])
            tx_gene[tx] = g.group(1) if g else tx
            tx_tags[tx] = set(re.findall(r'tag "([^"]+)"', t[8]))
            tx_strand[tx] = t[6]
            tx_chr[tx] = strip(t[0])
        else:
            rows.append((tx, t[2], int(t[3]), int(t[4])))
canon = {tx for tx, tg in tx_tags.items() if "Ensembl_canonical" in tg}
G = pd.DataFrame(rows, columns=["tx", "type", "start", "end"])
del rows
cds_tx = set(G.loc[G["type"] == "CDS", "tx"])
use = {tx for tx in canon if tx in cds_tx}
n0 = len(use)
use = {tx for tx in use if not ({"cds_start_NF", "cds_end_NF", "seleno"} & tx_tags[tx])}
log(f"canonical transcripts with CDS: {n0}; excluded (incomplete CDS or selenoprotein): {n0 - len(use)}")
G = G[G["tx"].isin(use)]

# ---------------- genome (codes 0..3, N = 4) ----------------
CODE = np.full(256, 4, np.uint8)
for b, i in BI.items():
    CODE[ord(b)] = i
    CODE[ord(b.lower())] = i
GEN = {}
need = {tx_chr[t] for t in use}
cur, buf = None, []
with open(a.genome) as f:
    for line in f:
        if line.startswith(">"):
            if cur in need:
                GEN[cur] = CODE[np.frombuffer("".join(buf).encode(), np.uint8)]
            cur = strip(line[1:].split()[0])
            buf = []
        elif cur in need:
            buf.append(line.strip())
if cur in need:
    GEN[cur] = CODE[np.frombuffer("".join(buf).encode(), np.uint8)]
log(f"chromosomes read from the genome: {len(GEN)}")

# ---------------- genetic code and severity table (as compute_coding_sites.R) ----------------
GCODE = {"".join(c): aa for c, aa in zip(itertools.product("TCAG", repeat=3),
                                         "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
LUT = np.zeros((4, 4, 4, 3, 4, 3), np.int8)              # [b1, b2, b3, position in codon, alt, codon context] -> severity
for c1, c2, c3 in itertools.product(range(4), repeat=3):
    cod = B[c1] + B[c2] + B[c3]
    for pp in range(3):
        for al in range(4):
            if cod[pp] == B[al]:
                continue
            alt = cod[:pp] + B[al] + cod[pp + 1:]
            for cx in range(3):                          # 0 internal codon, 1 first codon, 2 last (stop) codon
                if cx > 0 or GCODE[cod] == "*":
                    LUT[c1, c2, c3, pp, al, cx] = 3      # other coding
                elif GCODE[alt] == "*":
                    LUT[c1, c2, c3, pp, al, cx] = 5      # truncating (nonsense)
                elif GCODE[alt] == GCODE[cod]:
                    LUT[c1, c2, c3, pp, al, cx] = 2      # synonymous
                else:
                    LUT[c1, c2, c3, pp, al, cx] = 4      # missense


def orf_categories(tpos, exons_tx_order, minus, pen_nt):
    """Exon category and mRNA facts of every ORF nucleotide.
    tpos: ORF genomic positions in transcript order (5'->3'); exons_tx_order: (start, end) of the transcript exons in
    transcript order; minus: transcript on the - strand. Returns (category per nucleotide, info dict) or (None, reason)."""
    starts = np.array([s for s, e in exons_tx_order])
    ends = np.array([e for s, e in exons_tx_order])
    n_exons = len(starts)
    exon_len = ends - starts + 1
    # genomic position of every mRNA nucleotide, 5'->3', and the exon (1..n) it belongs to
    pieces = [np.arange(e, s - 1, -1) if minus else np.arange(s, e + 1) for s, e in zip(starts, ends)]
    mrna_gpos = np.concatenate(pieces)
    mrna_exon = np.repeat(np.arange(1, n_exons + 1), exon_len)
    order = np.argsort(mrna_gpos)
    j = np.clip(np.searchsorted(mrna_gpos[order], tpos), 0, len(mrna_gpos) - 1)
    if not (mrna_gpos[order][j] == tpos).all():
        return None, "ORF nucleotides outside the transcript exons"
    mrna_index = order[j] + 1                            # 1-based position on the mRNA
    exon_of_nt = mrna_exon[order[j]]
    last_junction = int(exon_len[:-1].sum())             # mRNA index of the last nucleotide before the last junction
    cat = np.full(len(tpos), CAT_INTERNAL, np.int64)
    cat[exon_of_nt == 1] = CAT_FIRST
    in_window = (exon_of_nt < n_exons) & (mrna_index > last_junction - pen_nt)
    cat[in_window] = CAT_PEN                             # pen55 has priority over first / internal
    cat[exon_of_nt == n_exons] = CAT_LAST                # single-exon genes end up here too
    info = {"n_exons": n_exons,
            "window_beyond_penultimate": bool(n_exons >= 2 and exon_len[-2] < pen_nt),
            "last_exon_without_orf": bool(n_exons >= 2 and not (exon_of_nt == n_exons).any())}
    return cat, info


def relative_cells(cat, n_bins):
    """Map cell of every ORF nucleotide: category + rank-based relative bin inside the category (transcript order)."""
    cell = np.full(len(cat), -1, np.int64)
    for c, base in ((CAT_FIRST, FIRST0), (CAT_INTERNAL, INTERNAL0), (CAT_LAST, LAST0)):
        idx = np.nonzero(cat == c)[0]                    # already in 5'->3' order
        if len(idx) == 0:
            continue
        rank = np.arange(len(idx))
        cell[idx] = base + (rank * n_bins) // len(idx)
    cell[cat == CAT_PEN] = PEN
    return cell


def tx_sites_cells(tx, d):
    """Possible substitutions of one transcript: (pos, alt base on + strand, severity, map cell), or (None, reason)."""
    seq = GEN.get(tx_chr[tx])
    if seq is None:
        return None, "chromosome missing from the genome"
    orf = d[d["type"].isin(["CDS", "stop_codon"])]
    pos = np.unique(np.concatenate([np.arange(s, e + 1) for s, e in zip(orf["start"], orf["end"])]))
    if len(pos) % 3 or len(pos) < 6:
        return None, "ORF not a multiple of 3 or shorter than 6 bp"
    minus = tx_strand[tx] == "-"
    tpos = pos[::-1] if minus else pos
    b = seq[tpos - 1].astype(np.int16)
    if (b > 3).any():
        return None, "non-ACGT base in the ORF"
    if minus:
        b = 3 - b                                        # complement (A0<->T3, C1<->G2)
    # severity of every substitution (same code as context_weighted_sites.py)
    L = len(b)
    ncod = L // 3
    cb = b.reshape(ncod, 3)
    ci = np.repeat(np.arange(ncod), 3)
    pin = np.tile(np.arange(3), ncod)
    cx = np.zeros(ncod, int)
    cx[0] = 1
    cx[-1] = 2
    P = np.repeat(tpos, 4)
    At = np.tile(np.arange(4), L)                        # alternative base on the transcript strand
    keep = At != np.repeat(b, 4)
    sev = LUT[cb[ci, 0].repeat(4), cb[ci, 1].repeat(4), cb[ci, 2].repeat(4), pin.repeat(4), At, cx[ci].repeat(4)]
    Ag = (3 - At) if minus else At                       # alternative base on the + strand
    # map cell of every ORF nucleotide
    ex_tx = d[d["type"] == "exon"].sort_values("start", ascending=not minus)
    cat, info = orf_categories(tpos, list(zip(ex_tx["start"], ex_tx["end"])), minus, a.pen_nt)
    if cat is None:
        return None, info
    cell_nt = relative_cells(cat, NB)
    info["orf_nt_per_category"] = np.bincount(cat, minlength=4)
    out_p, out_a, out_s, out_c = [P[keep]], [Ag[keep]], [sev[keep]], [np.repeat(cell_nt, 4)[keep]]
    # essential splice sites (+-1/+-2 intronic) of every junction, introns >= 4 bp: truncating, cell "splice"
    ex = d[d["type"] == "exon"].sort_values("start")
    es, ee = ex["start"].to_numpy(), ex["end"].to_numpy()
    for k in range(len(es) - 1):
        i1, i2 = ee[k] + 1, es[k + 1] - 1
        if i2 - i1 + 1 < 4:
            continue
        sp = np.array([i1, i1 + 1, i2 - 1, i2])
        ref = seq[sp - 1]
        Pp = np.repeat(sp, 4)
        Aa = np.tile(np.arange(4), 4)
        ok = Aa != np.repeat(ref, 4)
        out_p.append(Pp[ok])
        out_a.append(Aa[ok])
        out_s.append(np.full(ok.sum(), 5, np.int8))
        out_c.append(np.full(ok.sum(), SPLICE, np.int64))
    res = (np.concatenate(out_p).astype(np.int64), np.concatenate(out_a).astype(np.int64),
           np.concatenate(out_s).astype(np.int8), np.concatenate(out_c).astype(np.int64))
    return res, info


# ---------------- backbone bins (as gene_level_drivers_classes.py) ----------------
bb = pd.read_csv(a.backbone, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chr", "start", "end", "binID"],
                 dtype={"chr": str, "binID": str})
bb["chr"] = bb["chr"].map(strip)
bb["start"] += 1
bb = bb.reset_index(drop=True)
BB = {ch: (x["start"].to_numpy(), x["end"].to_numpy()) for ch, x in bb.sort_values(["chr", "start"]).groupby("chr")}
def in_backbone(ch, pos):
    if ch not in BB:
        return np.zeros(len(pos), bool)
    s, e = BB[ch]
    i = np.searchsorted(s, pos, side="right") - 1
    return (i >= 0) & (pos <= e[np.clip(i, 0, None)])

# ---------------- sites per gene and cell ----------------
genes = sorted({tx_gene[t] for t in use})
GI = {g: i for i, g in enumerate(genes)}
sites = np.zeros((len(genes), N_CELLS, 3), np.int64)      # sites inside the backbone, per gene x cell x class
GK = {}                                                  # gene -> [(chromosome, sorted keys, severity, cell)]
struct_rows, bad = [], {}
for tx, d in G.groupby("tx"):
    res, info = tx_sites_cells(tx, d)
    if res is None:
        bad[info] = bad.get(info, 0) + 1
        continue
    P, A_, S, C = res
    ch = tx_chr[tx]
    gi = GI[tx_gene[tx]]
    # one entry per substitution (position + base) in the transcript, most severe class (as the pipeline)
    key = P * 4 + A_
    o = np.lexsort((-S, key))
    key, S, C = key[o], S[o], C[o]
    first = np.r_[True, key[1:] != key[:-1]]
    key, S, C = key[first], S[first], C[first]
    GK.setdefault(tx_gene[tx], []).append((ch, key, S, C))
    cls = SEV_TO_CLASS[S]
    m = (cls >= 0) & in_backbone(ch, key // 4)
    np.add.at(sites, (gi, C[m], cls[m]), 1)
    struct_rows.append({"gene": tx_gene[tx], "transcript": tx, "n_exons": info["n_exons"],
                        **{f"orf_nt_{n}": int(v) for n, v in zip(CAT_NAMES, info["orf_nt_per_category"])},
                        "window_beyond_penultimate": info["window_beyond_penultimate"],
                        "last_exon_without_orf": info["last_exon_without_orf"]})
for reason, n in bad.items():
    log(f"transcripts excluded: {n} ({reason})")
ST = pd.DataFrame(struct_rows)
ST.to_csv(f"{a.out_dir}/transcript_structure.tsv", sep="\t", index=False)
log(f"transcripts mapped: {len(ST)} ({ST['gene'].nunique()} genes); single-exon: {int((ST['n_exons'] == 1).sum())}, "
    f"two-exon: {int((ST['n_exons'] == 2).sum())}, three or more: {int((ST['n_exons'] >= 3).sum())}")
log(f"pen55 window extends beyond the penultimate exon (penultimate exon < {a.pen_nt} nt): "
    f"{int(ST['window_beyond_penultimate'].sum())} transcripts; last exon without ORF (stop codon in an earlier exon): "
    f"{int(ST['last_exon_without_orf'].sum())} transcripts")
tot_nt = ST[[f"orf_nt_{n}" for n in CAT_NAMES]].sum()
log("share of ORF nucleotides per category: " + ", ".join(f"{n} {tot_nt[f'orf_nt_{n}'] / tot_nt.sum():.3f}" for n in CAT_NAMES))

# ---------------- events (same rules as gene_level_drivers_classes.py) ----------------
TRUNC = {"stop_gained", "splice_acceptor_variant", "splice_donor_variant"}
MIS = {"missense_variant"}
OTHER = {"start_lost", "stop_lost", "stop_retained_variant", "coding_sequence_variant", "protein_altering_variant",
         "incomplete_terminal_codon_variant", "start_retained_variant"}
SYN = {"synonymous_variant"}
INDEL = {"frameshift_variant", "inframe_insertion", "inframe_deletion"}
CODING_ANY = TRUNC | MIS | OTHER | SYN | INDEL
cols = ["COSMIC_SAMPLE_ID", "GENOMIC_MUTATION_ID", "CHROMOSOME", "GENOME_START", "GENOME_STOP", "TRANSCRIPT_ACCESSION",
        "MUTATION_DESCRIPTION", "GENOMIC_WT_ALLELE", "GENOMIC_MUT_ALLELE"]
R = pd.read_csv(a.events, sep="\t", dtype=str, quoting=3, usecols=cols, keep_default_na=False)
R["mut_id"] = np.where(R["GENOMIC_MUTATION_ID"] == "",
                       R["CHROMOSOME"] + ":" + R["GENOME_START"] + ":" + R["GENOME_STOP"], R["GENOMIC_MUTATION_ID"])
R["pos"] = pd.to_numeric(R["GENOME_START"], errors="coerce")
R["stop"] = pd.to_numeric(R["GENOME_STOP"], errors="coerce")
R = R[R["pos"].notna() & (R["COSMIC_SAMPLE_ID"] != "")]
ud = pd.Series(R["MUTATION_DESCRIPTION"].unique())
coding_map = dict(zip(ud, ud.map(lambda x: bool(set(str(x).split(",")) & CODING_ANY))))
R["coding_any"] = R["MUTATION_DESCRIPTION"].map(coding_map)
R["snv"] = ((R["pos"] == R["stop"]) & R["GENOMIC_WT_ALLELE"].str.upper().isin(set(B)) &
            R["GENOMIC_MUT_ALLELE"].str.upper().isin(set(B)))
per_sample = R.groupby(["COSMIC_SAMPLE_ID", "mut_id"])["coding_any"].any().groupby(level=0).sum()
if str(a.hypermutator_max).lower() not in ("none", ""):
    keep_s = set(per_sample[per_sample <= float(a.hypermutator_max)].index)
    log(f"samples: {len(per_sample)}, excluded as hypermutators: {len(per_sample) - len(keep_s)} "
        f"(> {a.hypermutator_max} coding events)")
    R = R[R["COSMIC_SAMPLE_ID"].isin(keep_s)]
ev = R.groupby(["COSMIC_SAMPLE_ID", "mut_id"], sort=False).agg(
    chr=("CHROMOSOME", "first"), pos=("pos", "first"), snv=("snv", "first"),
    wt=("GENOMIC_WT_ALLELE", "first"), mt=("GENOMIC_MUT_ALLELE", "first")).reset_index()
del R
ev["chr"] = ev["chr"].map(strip)
ev = ev[ev["snv"]].copy()
ev["pos"] = ev["pos"].astype(np.int64)
ev["wt"] = ev["wt"].str.upper()
ev["mt"] = ev["mt"].str.upper()
ev["key"] = ev["pos"] * 4 + ev["mt"].map(BI).astype(np.int64)
ev["ref_ok"] = [c in GEN and p_ <= len(GEN[c]) and GEN[c][p_ - 1] == BI[w]
                for c, p_, w in zip(ev["chr"], ev["pos"], ev["wt"])]
log(f"single-base substitution events: {len(ev)}; COSMIC reference base different from the genome (excluded): "
    f"{int((~ev['ref_ok']).sum())}")
ev = ev[ev["ref_ok"]]
ev["in_backbone"] = False
for ch, idx in ev.groupby("chr").groups.items():
    ev.loc[idx, "in_backbone"] = in_backbone(ch, ev.loc[idx, "pos"].to_numpy(np.int64))

# class and cell of every event in every gene whose canonical transcript has that site
by_chr = {ch: x for ch, x in ev.groupby("chr")}
eg_parts = []
for gene, parts in GK.items():
    for ch, k, s, c in parts:
        if ch not in by_chr:
            continue
        x = by_chr[ch]
        keys = x["key"].to_numpy()
        i = np.clip(np.searchsorted(k, keys), 0, len(k) - 1)
        hit = k[i] == keys
        if not hit.any():
            continue
        eg_parts.append(pd.DataFrame({"sample": x["COSMIC_SAMPLE_ID"].to_numpy()[hit], "mut_id": x["mut_id"].to_numpy()[hit],
                                      "gene": gene, "sev": s[i[hit]], "cell": c[i[hit]],
                                      "in_backbone": x["in_backbone"].to_numpy()[hit]}))
if not eg_parts:
    sys.exit("ERROR: no event falls on a site of a canonical transcript: check --events, --gtf and --genome (same assembly?)")
EG = pd.concat(eg_parts, ignore_index=True)
# one event per gene: most severe class among the gene's canonical transcripts (as gene_level_drivers_classes.py)
EG = EG.sort_values("sev", ascending=False).drop_duplicates(["sample", "mut_id", "gene"])
EG["cls"] = SEV_TO_CLASS[EG["sev"].to_numpy()]
EG = EG[EG["cls"] >= 0]
events = np.zeros((len(genes), N_CELLS, 3), np.int64)
inb = EG[EG["in_backbone"]]
np.add.at(events, (inb["gene"].map(GI).to_numpy(), inb["cell"].to_numpy(), inb["cls"].to_numpy()), 1)
n_out = (~EG["in_backbone"]).groupby(EG["cls"]).sum()
log("events at positions outside the backbone bins (not in the map): " +
    ", ".join(f"{CLASS_NAMES[c]} {int(n_out.get(c, 0))}" for c in range(3)))
# events can exceed sites (the same site mutated in several samples), but an event in a cell without sites means the
# event and site mappings disagree
orphan = int(((events > 0) & (sites == 0)).sum())
if orphan:
    sys.exit(f"ERROR: {orphan} gene x cell x class entries have events but no sites: event and site mapping disagree")

# ---------------- outputs: long table ----------------
gi_, ci_, ki_ = np.nonzero(sites)
LONG = pd.DataFrame({"gene": np.array(genes)[gi_], "cell": np.array(CELL_NAMES)[ci_], "cls": np.array(CLASS_NAMES)[ki_],
                     "sites": sites[gi_, ci_, ki_], "events": events[gi_, ci_, ki_]})
LONG.to_csv(f"{a.out_dir}/gene_map_long.tsv", sep="\t", index=False)
log(f"map: {len(genes)} genes x {N_CELLS} cells; gene x cell x class entries with sites: {len(LONG)}")
for k, c in enumerate(CLASS_NAMES):
    has = sites[:, :, k] > 0
    zero_share = float((events[:, :, k][has] == 0).mean()) if has.any() else float("nan")
    log(f"  {c:5s}: sites {int(sites[:, :, k].sum())}, events {int(events[:, :, k].sum())}, "
        f"share of gene x cell entries (with sites) without events {zero_share:.3f}")

# ---------------- outputs: per-cell summary ----------------
cell_sites = sites.sum(axis=0)                           # cell x class
cell_events = events.sum(axis=0)
genes_with_sites = (sites > 0).sum(axis=0)
genes_with_events = (events > 0).sum(axis=0)
syn_pooled_rate = cell_events[:, 0].sum() / max(cell_sites[:, 0].sum(), 1)
syn_pooled_events = cell_events[:, 0].sum()
summ = []
for ci in range(N_CELLS):
    for k, c in enumerate(CLASS_NAMES):
        if cell_sites[ci, k] == 0:
            continue
        rate = cell_events[ci, k] / cell_sites[ci, k]
        row = {"cell": CELL_NAMES[ci], "cls": c, "genes_with_sites": int(genes_with_sites[ci, k]),
               "sites": int(cell_sites[ci, k]), "events": int(cell_events[ci, k]),
               "genes_with_events": int(genes_with_events[ci, k]), "events_per_Msite": 1e6 * rate,
               "ratio_vs_syn": np.nan, "ratio_lo": np.nan, "ratio_hi": np.nan, "syn_reference": ""}
        if c != "syn":
            # reference: synonymous of the same cell; the splice cell has no synonymous sites -> pooled synonymous
            if cell_sites[ci, 0] > 0:
                syn_rate, syn_ev, row["syn_reference"] = cell_events[ci, 0] / cell_sites[ci, 0], cell_events[ci, 0], "same cell"
            else:
                syn_rate, syn_ev, row["syn_reference"] = syn_pooled_rate, syn_pooled_events, "all cells pooled"
            if cell_events[ci, k] > 0 and syn_ev > 0:
                ratio = rate / syn_rate
                se = np.sqrt(1 / cell_events[ci, k] + 1 / syn_ev)    # Wald, log scale, Poisson counts
                row.update(ratio_vs_syn=ratio, ratio_lo=ratio * np.exp(-1.96 * se), ratio_hi=ratio * np.exp(1.96 * se))
        summ.append(row)
SUMM = pd.DataFrame(summ)
SUMM.to_csv(f"{a.out_dir}/cell_summary.tsv", sep="\t", index=False)
log("per-cell summary (events, ratio vs synonymous):\n" +
    SUMM.pivot(index="cell", columns="cls", values="events").reindex(CELL_NAMES).fillna(0).astype(int)
        .join(SUMM[SUMM["cls"] != "syn"].pivot(index="cell", columns="cls", values="ratio_vs_syn").add_prefix("ratio_"))
        .round(3).to_string())

# ---------------- regression check against gene_level_drivers_classes.py ----------------
if a.gene_scores_dir:
    O_all = EG.groupby(["gene", "cls"]).size()            # events WITHOUT the backbone restriction (as column 'O')
    checks = []
    for k, c in enumerate(CLASS_NAMES):
        path = f"{a.gene_scores_dir}/gene_scores_{c}.tsv"
        if not os.path.exists(path):
            log(f"check: {path} not found, class {c} not checked")
            continue
        GSC = pd.read_csv(path, sep="\t", index_col=0)
        GSC = GSC[GSC["sites"].fillna(0) > 0]              # genes with sites of this class (the file holds all classes' genes)
        mine_sites = pd.Series(sites[:, :, k].sum(axis=1), index=genes)
        mine_sites = mine_sites[mine_sites > 0]
        mine_O = O_all.xs(k, level="cls") if k in O_all.index.get_level_values("cls") else pd.Series(dtype=float)
        common = GSC.index.intersection(mine_sites.index)
        d_sites = int((GSC.loc[common, "sites"].to_numpy() != mine_sites[common].to_numpy()).sum())
        d_O = int((GSC.loc[common, "O"].to_numpy() != mine_O.reindex(common).fillna(0).to_numpy()).sum())
        checks.append({"cls": c, "genes_gene_scores": len(GSC), "genes_map": len(mine_sites), "genes_common": len(common),
                       "genes_sites_diff": d_sites, "sites_total_gene_scores": int(GSC.loc[common, "sites"].sum()),
                       "sites_total_map": int(mine_sites[common].sum()), "genes_O_diff": d_O,
                       "O_total_gene_scores": int(GSC.loc[common, "O"].sum()),
                       "O_total_map_unrestricted": int(mine_O.reindex(common).fillna(0).sum())})
        log(f"CHECK {c:5s}: genes {len(common)} in common (gene_scores {len(GSC)}, map {len(mine_sites)}); "
            f"genes with different sites {d_sites}, with different events {d_O}")
    CH = pd.DataFrame(checks)
    CH.to_csv(f"{a.out_dir}/check_vs_gene_scores.tsv", sep="\t", index=False)
    if len(CH) and (CH["genes_sites_diff"].sum() or CH["genes_O_diff"].sum() or (CH["genes_common"] != CH["genes_map"]).any()
                    or (CH["genes_common"] != CH["genes_gene_scores"]).any()):
        log("WARNING: the map does not reproduce gene_scores: do NOT use it before understanding the differences")

open(f"{a.out_dir}/log.txt", "w").write("\n".join(LOG) + "\n")
print(f"Output in {a.out_dir}/", file=sys.stderr)
