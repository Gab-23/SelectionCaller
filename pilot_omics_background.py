"""
BURNER - pilot del background agnostico per omica, con modelli leggeri.

Per ogni omica (CNA ampl, CNA del, ATAC logFC woCNAs, mutazioni codificanti) e per ogni
cromosoma (modelli cromosoma-specifici):
  * stesse feature agnostiche per tutte le omiche (nessuna colonna di alterazione, niente
    distanza da OG/TSG che restano per la validazione);
  * XGBoost con iperparametri fissi (niente BayesSearch: e' un pilota);
  * CV a blocchi intercalati con cuscinetto (INTERLEAVED_KFOLD_IDXS);
  * predizioni e SHAP out-of-fold (pred_contribs di XGBoost, esatti per gli alberi);
  * feature "esca" (decoy): copie di feature reali spostate circolarmente lungo il cromosoma,
    stessa distribuzione e autocorrelazione ma nessun allineamento col target -> soglia per
    dire che una feature e' davvero "usata".
Poi: curva R2 vs cuscinetto per le CNA; nucleo condiviso vs moduli specifici di omica.

NUOVO (20/10/2026):
  --mode chrom   un modello per cromosoma (come prima); unita' di replica = cromosoma
  --mode genome  un modello per omica su tutti i cromosomi (fold intercalate dentro ogni
                 cromosoma, cuscinetto mai tra cromosomi); unita' di replica = fold.
                 Un modello genome-wide non puo' memorizzare il profilo di un singolo
                 cromosoma tramite le coordinate.
  baseline lineare: ridge (imputazione mediana + standardizzazione, alpha fisso) sulle
                 STESSE fold, per misurare quanto aggiunge un modello non lineare.
  --loss glm     (21/10/2026) perdita per famiglia: CNA binomiale (frazione di pazienti), mutazioni Poisson
                 con offset log(kb CDS) sul conteggio dei campioni, ATAC gaussiana. Metrica: devianza
                 spiegata rispetto al modello nullo (solo media, o solo offset per le mutazioni); il residuo
                 salvato come 'resid' diventa il residuo di Pearson (osservato - atteso)/sqrt(varianza attesa).
                 Baseline lineare = XGBoost gblinear con la stessa perdita (GLM penalizzato).
  --orient_eigen gc   orienta gli autovettori Hi-C per cromosoma (segno tale che la
                 correlazione con GC_content sia positiva). Provvisorio: da togliere
                 (--orient_eigen none) quando l'orientamento sara' fatto in pipeline.

Uso (env py_ml_env):
  python pilot_omics_background.py --dataset <feature_matrix> --out_dir pilot_out [--chroms 1 8 17]
"""
import argparse, os, sys, time
import numpy as np, pandas as pd
import xgboost as xgb
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

sys.path.insert(0, os.getcwd())
try:
    from workflow.scripts.helpers.INTERLEAVED_KFOLD_IDXS import INTERLEAVED_KFOLD_IDXS
except ImportError:
    sys.exit("Serve workflow/scripts/helpers/INTERLEAVED_KFOLD_IDXS.py (lancia dalla root della repo).")

p = argparse.ArgumentParser()
p.add_argument("--dataset", required=True)
p.add_argument("--out_dir", default="pilot_out")
p.add_argument("--omics", nargs="+", default=["cna_ampl", "cna_del", "atac", "mut_syn", "mut_mis", "mut_trunc"],
               help="mut_syn / mut_mis / mut_trunc: eventi per classe (mut_n_*) con offset log(siti della classe, sites_*), "
                    "solo con --loss glm; mut = target vecchio (mut_col, offset log kb CDS), tenuto per confronto")
p.add_argument("--chroms", nargs="*", default=None, help="default: tutti")
p.add_argument("--bin_size", type=int, default=20000)
p.add_argument("--n_splits", type=int, default=5)
p.add_argument("--block_bins", type=int, default=100, help="dimensione dei blocchi intercalati (bin)")
p.add_argument("--max_buffer_bins", type=int, default=100, help="tetto al cuscinetto stimato da l (bin)")
p.add_argument("--buffer_sweep_mb", nargs="*", type=float, default=[0, 1, 2, 5, 10], help="curva R2 vs cuscinetto (solo CNA)")
p.add_argument("--n_decoys", type=int, default=10)
p.add_argument("--shadow_mode", choices=["random", "per_feature"], default="random",
               help="random = esche da feature scelte a caso (come prima); per_feature = per OGNI feature --n_shadows copie "
                    "di se' stessa spostate circolarmente nel cromosoma (alla Boruta): un gruppo e' 'usato' (used_shadow) se la "
                    "sua importanza supera quella delle copie spostate dei SUOI membri (stessa dimensione del gruppo)")
p.add_argument("--n_shadows", type=int, default=10, help="copie spostate per feature con --shadow_mode per_feature (con 3 copie, nei test, ~10%% di gruppi irrilevanti passavano per caso; con 10 nessuno)")
p.add_argument("--floor_importance", action="store_true",
               help="calcola gli SHAP anche nei modelli del pavimento e applica gli stessi criteri di uso: un gruppo usato "
                    "anche con il target spostato deve la sua importanza alla struttura spaziale, non all'occorrenza")
p.add_argument("--used_quantile", type=float, default=0.95, help="feature usata se importanza > questo quantile dei decoy")
p.add_argument("--used_frac_chroms", type=float, default=0.5, help="...in almeno questa frazione delle repliche (cromosomi o fold) con R2>0")
p.add_argument("--group_corr", type=float, default=0.8, help="feature con |rho| >= soglia raggruppate")
p.add_argument("--mut_col", default="auto", help="auto: mut_n_samples se presente (robusto agli ipermutati), altrimenti length_coding_agnostic")
p.add_argument("--exonic_col", default="auto", help="colonna bp codificanti (auto: cds_bp se presente); 'none' = proxy num.overlapping.genes_protein_coding > 0")
p.add_argument("--mut_target", choices=["density", "count"], default="density",
               help="density: log1p(eventi per kb di CDS); count: log1p(eventi). density richiede la colonna dei bp codificanti")
p.add_argument("--min_exonic_bp", type=float, default=500, help="bin delle mutazioni: almeno questi bp codificanti")
p.add_argument("--exclude_features", nargs="*", default=[], help="feature da togliere (es. dist_to_centromere dist_to_telomere per un controllo senza coordinate)")
p.add_argument("--mode", choices=["chrom", "genome"], default="chrom")
p.add_argument("--loss", choices=["squared", "glm"], default="squared",
               help="squared: errore quadratico (come prima). glm: perdita della famiglia giusta per ogni omica -"
                    " CNA binomiale sulla frazione di pazienti, mutazioni Poisson sul conteggio con offset log(kb CDS),"
                    " ATAC resta gaussiana (logFC continuo)")
p.add_argument("--orient_eigen", choices=["gc", "none"], default="gc")
p.add_argument("--ridge_alpha", type=float, default=10.0, help="baseline lineare: alpha della ridge su feature standardizzate")
p.add_argument("--n_floor", type=int, default=0,
               help="PAVIMENTO: numero di ripetizioni con il target spostato circolarmente dentro ogni cromosoma "
                    "(stessa CV, stesse feature): conserva l'autocorrelazione e rompe il legame con le feature. "
                    "0 = non calcolato")
p.add_argument("--floor_scope", choices=["chrom", "genome"], default="chrom",
               help="chrom: shift dentro ogni cromosoma (conserva la media di ogni cromosoma: il pavimento include quanto "
                    "si spiega solo dall'identita' del cromosoma); genome: shift dell'intero vettore concatenato (rompe "
                    "anche le differenze tra cromosomi). In mode chrom i due coincidono")
p.add_argument("--exposure_file", default=None, help="con --class_exposure wsites: bin_weighted_sites.tsv di context_weighted_sites.py")
p.add_argument("--no_predict_low_cds", action="store_true",
               help="NOTE 28/10/2026 - di default i modelli delle classi di mutazioni (mut_syn/mis/trunc) si ALLENANO sui bin con "
                    ">= --min_exonic_bp di CDS, ma PREDICONO anche i bin con meno CDS e almeno un sito della classe (out-of-fold: "
                    "ogni bin aggiunto e' predetto dal modello della fold del suo blocco, stessa regola di INTERLEAVED_KFOLD_IDXS). "
                    "Le metriche restano sui soli bin di training. Con questa opzione: comportamento precedente (nessun bin aggiunto)")
p.add_argument("--class_exposure", choices=["sites", "cds", "wsites"], default="sites",
               help="offset dei target per classe (mut_syn / mut_mis / mut_trunc): sites = log(siti della classe, sites_*); "
                    "cds = log(bp di CDS del bin, --exonic_col), uguale per tutte le classi (controllo: esposizione per lunghezza); "
                    "wsites = log(siti pesati per il contesto trinucleotidico, wsites_* da --exposure_file)")
p.add_argument("--n_jobs", type=int, default=4)
p.add_argument("--seed", type=int, default=1)
a = p.parse_args()
os.makedirs(a.out_dir, exist_ok=True)
rng = np.random.default_rng(a.seed)
t0 = time.time()

# ---- dati e feature agnostiche ---------------------------------------------------------------
D = pd.read_csv(a.dataset, sep="\t")
if a.class_exposure == "wsites":                                       # siti pesati per contesto, calcolati fuori dalla pipeline
    if not a.exposure_file: sys.exit("--class_exposure wsites richiede --exposure_file")
    _X = pd.read_csv(a.exposure_file, sep="\t", dtype={"binID": str})
    D = D.merge(_X, on="binID", how="left")
    print(f"Esposizione pesata per contesto: {_X.shape[0]} bin da {a.exposure_file}; bin del dataset senza valore: {int(D['wsites_syn'].isna().sum())}")
D["chr"] = D["binID"].str.split("_").str[0]
D["idx"] = D["binID"].str.split("_").str[1].astype(int)
D = D.sort_values(["chr", "idx"], key=lambda s: s.map(lambda v: int(v) if str(v).isdigit() else v) if s.name == "chr" else s).reset_index(drop=True)
ALT = [c for c in D.columns if c in ("ampl.frequency", "del.frequency", "logFC", "FDR")
       or c.startswith("length_coding") or c.startswith("length_noncoding") or c.startswith("mut_n_")]
if a.mut_col == "auto":
    a.mut_col = "mut_n_samples" if "mut_n_samples" in D.columns else "length_coding_agnostic"
if a.exonic_col == "auto":
    a.exonic_col = "cds_bp" if "cds_bp" in D.columns else None
elif a.exonic_col == "none":
    a.exonic_col = None
if a.mut_target == "density" and a.exonic_col is None:
    print("ATTENZIONE: --mut_target density richiede i bp codificanti: uso count"); a.mut_target = "count"
if a.loss == "glm" and a.exonic_col:
    print(f"Mutazioni: target = {a.mut_col} (conteggio, Poisson con offset log(kb CDS); --mut_target ignorato), bp codificanti = {a.exonic_col}")
else:
    print(f"Mutazioni: target = {a.mut_col} ({a.mut_target}), bp codificanti = {a.exonic_col}")
VALIDATION_ONLY = [c for c in ("dist.to.closest.OG", "dist.to.closest.TSG") if c in D.columns]
# la colonna dei bp esonici serve solo a filtrare i bin delle mutazioni: mai come feature
EXONIC = [a.exonic_col] if a.exonic_col else []
# siti per classe (sites_*): ESPOSIZIONE dei target per classe, usata come offset: mai come feature, altrimenti il
# modello potrebbe usarla al posto dell'offset
EXPOSURE = [c for c in D.columns if c.startswith(("sites_", "wsites_"))]
if a.exonic_col and a.exonic_col not in D.columns: sys.exit(f"--exonic_col {a.exonic_col} non e' una colonna del dataset")
FEATS = [c for c in D.columns if c not in ALT + VALIDATION_ONLY + EXONIC + EXPOSURE + list(a.exclude_features) + ["binID", "chr", "idx"]]
for c in [c for c in FEATS if c.startswith("dist")]:
    D[c] = D[c].fillna(2e6)                                    # NA su cromosomi interi = nessun gene della lista: cap
FEATS = [c for c in FEATS if D[c].isna().mean() < 0.5]
print(f"[{time.time()-t0:.0f}s] {len(D)} bin, {len(FEATS)} feature agnostiche; escluse: {len(ALT)} di alterazione + {VALIDATION_ONLY} + {EXONIC} + {EXPOSURE} + {list(a.exclude_features)}")

# orientamento degli autovettori Hi-C per cromosoma (segno arbitrario nell'output di Hi-C)
if a.orient_eigen == "gc":
    eig = [c for c in FEATS if c.startswith("eigen_")]
    n_flip = 0
    for ch, idx in D.groupby("chr").groups.items():
        for c in eig:
            r = D.loc[idx, c].corr(D.loc[idx, "GC_content"], method="spearman")
            if pd.notna(r) and r < 0:
                D.loc[idx, c] = -D.loc[idx, c]; n_flip += 1
    print(f"Orientamento autovettori (corr con GC > 0): {n_flip} segni invertiti su {len(eig) * D['chr'].nunique()} (autovettore x cromosoma)")

TARGETS = {"cna_ampl": ("ampl.frequency", None), "cna_del": ("del.frequency", None),
           "atac": ("logFC", None), "mut": (a.mut_col, "log1p")}
# eventi per classe (campione x mutazione, sostituzioni singole) e loro esposizione (sostituzioni possibili)
MUT_CLASSES = {"mut_syn": ("mut_n_syn", "sites_syn"), "mut_mis": ("mut_n_mis", "sites_mis"), "mut_trunc": ("mut_n_trunc", "sites_trunc")}
for om_, (col_, _) in MUT_CLASSES.items(): TARGETS[om_] = (col_, None)
if any(o in MUT_CLASSES for o in a.omics):
    print("Mutazioni per classe: esposizione (offset) = " + {"sites": "siti della classe (sites_*)", "cds": f"bp di CDS del bin ({a.exonic_col})",
                                                             "wsites": "siti della classe pesati per il contesto (wsites_*)"}[a.class_exposure])
if any(o in MUT_CLASSES for o in a.omics) and a.class_exposure == "cds" and not a.exonic_col:
    sys.exit("--class_exposure cds richiede la colonna dei bp di CDS (--exonic_col)")
if any(o in MUT_CLASSES for o in a.omics) and a.loss != "glm":
    sys.exit("i target per classe (mut_syn / mut_mis / mut_trunc) richiedono --loss glm (Poisson con offset dei siti)")
if a.exonic_col is None and "mut" in a.omics:
    print("ATTENZIONE: nessuna colonna di bp esonici -> per le mutazioni uso i bin con geni codificanti (proxy, include introni)")
chroms = a.chroms or sorted(D["chr"].unique(), key=lambda v: int(v) if v.isdigit() else 99)

def corr_len(y, max_lag=3000):
    y = y - y.mean(); v = np.dot(y, y)
    if v == 0: return 1
    for lag in range(1, min(max_lag, len(y) - 1)):
        if np.dot(y[:-lag], y[lag:]) / v < 0.5: return lag
    return max_lag

def deviance(y, mu, family):
    """Devianza (a meno di costanti moltiplicative, che si semplificano nel rapporto)."""
    if family == "gaussian":
        return np.sum((y - mu) ** 2)
    if family == "poisson":
        mu = np.clip(mu, 1e-12, None)
        return 2 * np.sum(np.where(y > 0, y * np.log(np.where(y > 0, y, 1) / mu), 0.0) - (y - mu))
    if family == "binomial":                                   # y = frazione di pazienti; N costante -> si semplifica
        mu = np.clip(mu, 1e-9, 1 - 1e-9)
        t1 = np.where(y > 0, y * np.log(np.where(y > 0, y, 1) / mu), 0.0)
        t0 = np.where(y < 1, (1 - y) * np.log(np.where(y < 1, 1 - y, 1) / (1 - mu)), 0.0)
        return 2 * np.sum(t1 + t0)
    raise ValueError(family)

XGB_OBJ = {"gaussian": "reg:squarederror", "poisson": "count:poisson", "binomial": "binary:logistic"}

def extra_fold_assignment(extra_binIDs, n_splits, block_bins):
    """Fold di ogni bin AGGIUNTO (non usato in training), con la STESSA regola deterministica di INTERLEAVED_KFOLD_IDXS:
    fold = ((numero del bin - 1) // block_bins) % n_splits, dal numero nel binID "chr_n". Cosi' il bin aggiunto e' predetto dal
    modello della fold del SUO blocco, che non ha visto quel blocco.
    FIX 28/10/2026: prima si usava la fold del bin di training piu' vicino; vicino ai confini tra blocchi poteva essere la
    fold di un blocco adiacente (stimato: <= 1% dell'atteso totale). LIMITE NOTO (non corretto): il cuscinetto di ogni
    modello e' calcolato attorno ai bin di TEST di training, non attorno ai bin aggiunti; un bin di training a distanza
    < cuscinetto da un bin aggiunto (ma lontano dai bin di test) puo' restare nel training del modello che lo predice."""
    n = pd.Series(extra_binIDs, dtype=str).str.rsplit("_", n=1).str[1].astype(np.int64).to_numpy()
    return ((n - 1) // block_bins) % n_splits

def fit_oof(X, y, binIDs, buffer_bins, want_shap=True, want_ridge=True, family="gaussian", offset=None, init_rate=None, extra=None):
    """OOF con fold intercalate + cuscinetto. family: gaussian / poisson (con offset = log esposizione) / binomial.
    init_rate (solo poisson): se True, il margine di partenza e' offset + log(tasso medio del training), invece del
    solo offset. Con base_margin XGBoost ignora base_score e parte da un tasso di 1 per unita' di esposizione: va bene
    per i kb di CDS (tassi ~1-4 per kb), ma per i SITI il tasso e' ~1e-3, e con max_delta_step (0.7) e eta (0.05) il
    modello impiegherebbe ~200 alberi solo per raggiungere il tasso medio (verificato su dati simulati: devianza
    spiegata negativa). La costante non cambia il modello statistico, solo il punto di partenza dell'ottimizzazione.
    Ritorna pred (media attesa sulla scala del target), shap (scala del margine), R2 xgb, R2 lineare, pred lineare,
    fold di ogni bin, e un dict con la devianza spiegata (xgb e lineare) rispetto al nullo stimato sul training."""
    folds = INTERLEAVED_KFOLD_IDXS(binIDs, a.n_splits, a.block_bins, buffer_bins, MIN_TRAIN_FRACTION=0.05)
    off = np.zeros(len(y)) if offset is None else np.asarray(offset, float)
    if extra is not None:                                      # bin aggiunti: (X_extra, binIDs_extra, offset_extra), solo poisson
        ex_fold = extra_fold_assignment(extra[1], a.n_splits, a.block_bins); pred_extra = np.full(len(extra[1]), np.nan)
    pred = np.full(len(y), np.nan); pred_lin = np.full(len(y), np.nan); null = np.full(len(y), np.nan)
    shap = np.zeros(X.shape) if want_shap else None
    fold_of = np.full(len(y), -1)
    params = dict(objective=XGB_OBJ[family], max_depth=5, eta=0.05, subsample=0.8, colsample_bytree=0.8,
                  min_child_weight=10, tree_method="hist", nthread=a.n_jobs, seed=a.seed)
    for k, (tr, te) in folds.items():
        fold_of[te] = int(k)
        if family == "gaussian":
            m = xgb.XGBRegressor(n_estimators=300, max_depth=5, learning_rate=0.05, subsample=0.8, colsample_bytree=0.8,
                                 min_child_weight=10, tree_method="hist", n_jobs=a.n_jobs, random_state=a.seed)
            m.fit(X[tr], y[tr]); pred[te] = m.predict(X[te]); booster = m.get_booster()
            dte = xgb.DMatrix(X[te]); null[te] = y[tr].mean()
        else:
            bm = off
            if family == "poisson" and init_rate:
                bm = off + np.log(max(y[tr].sum(), 1e-12) / np.exp(off[tr]).sum())   # tasso medio del training (no test)
            dtr = xgb.DMatrix(X[tr], label=y[tr], base_margin=bm[tr] if family == "poisson" else None)
            dte = xgb.DMatrix(X[te], base_margin=bm[te] if family == "poisson" else None)
            booster = xgb.train(params, dtr, num_boost_round=300)
            pred[te] = booster.predict(dte)                    # media attesa (conteggio atteso / frazione attesa)
            if extra is not None and family == "poisson" and (ex_fold == int(k)).any():
                mk = ex_fold == int(k); ebm = extra[2][mk] + (bm[tr][0] - off[tr][0])     # stessa costante del training
                pred_extra[mk] = booster.predict(xgb.DMatrix(extra[0][mk], base_margin=ebm))
            if family == "poisson":                            # nullo: tasso costante per unita' di esposizione
                null[te] = np.exp(off[te]) * y[tr].sum() / np.exp(off[tr]).sum()
            else:
                null[te] = y[tr].mean()
        if want_shap:
            shap[te] = booster.predict(dte, pred_contribs=True)[:, :-1]   # ultima colonna = bias (scala del margine)
        if want_ridge:
            imp_sc = make_pipeline(SimpleImputer(strategy="median"), StandardScaler()).fit(X[tr])
            Ztr, Zte = imp_sc.transform(X[tr]), imp_sc.transform(X[te])
            if family == "gaussian":
                lin = Ridge(alpha=a.ridge_alpha).fit(Ztr, y[tr]); pred_lin[te] = lin.predict(Zte)
            else:                                              # GLM penalizzato (L2) con la stessa perdita
                ltr = xgb.DMatrix(Ztr, label=y[tr], base_margin=bm[tr] if family == "poisson" else None)
                lte = xgb.DMatrix(Zte, base_margin=bm[te] if family == "poisson" else None)
                lb = xgb.train(dict(booster="gblinear", objective=XGB_OBJ[family], reg_lambda=a.ridge_alpha / len(tr),
                                    eta=0.5, nthread=a.n_jobs, seed=a.seed), ltr, num_boost_round=200)
                pred_lin[te] = lb.predict(lte)
    r2 = lambda p: 1 - np.sum((y - p) ** 2) / np.sum((y - y.mean()) ** 2)
    d0 = deviance(y, null, family)
    dev = dict(dev_expl=1 - deviance(y, pred, family) / d0,
               dev_expl_lin=(1 - deviance(y, pred_lin, family) / d0) if want_ridge else np.nan)
    if extra is not None: dev["pred_extra"] = pred_extra
    return pred, shap, r2(pred), (r2(pred_lin) if want_ridge else np.nan), pred_lin, fold_of, dev

def pearson_resid(y, mu, family):
    if family == "poisson":  return (y - mu) / np.sqrt(np.clip(mu, 1e-12, None))
    if family == "binomial": return (y - mu) / np.sqrt(np.clip(mu * (1 - mu), 1e-12, None))
    return y - mu

def r2_by_chrom(y, p, chroms):
    out = {}
    for c in np.unique(chroms):
        m = chroms == c
        out[c] = 1 - np.sum((y[m] - p[m]) ** 2) / np.sum((y[m] - y[m].mean()) ** 2)
    return out

def make_decoys(X, chroms):
    """Esche: feature reali spostate circolarmente DENTRO ogni cromosoma (offset 10-90%)."""
    src = rng.choice(X.shape[1], size=a.n_decoys, replace=True)
    dec = np.empty((X.shape[0], a.n_decoys))
    for i, j in enumerate(src):
        for c in np.unique(chroms):
            m = np.where(chroms == c)[0]; n = len(m)
            dec[m, i] = np.roll(X[m, j], int(rng.integers(int(.1 * n), max(int(.1 * n) + 1, int(.9 * n)))))
    return dec

def make_shadows(X, names, chroms):
    """Per ogni feature, a.n_shadows copie di se' stessa spostate circolarmente DENTRO ogni cromosoma (offset 10-90%)."""
    out, nm = [], []
    for j in range(a.n_shadows):
        S_ = np.empty_like(X)
        for c in np.unique(chroms):
            m = np.where(chroms == c)[0]; n = len(m)
            for f in range(X.shape[1]):
                S_[m, f] = np.roll(X[m, f], int(rng.integers(int(.1 * n), max(int(.1 * n) + 1, int(.9 * n)))))
        out.append(S_); nm += [f"SHADOW{j}__{f}" for f in names]
    return np.column_stack(out), nm

def circular_shift_index(chroms):
    if a.floor_scope == "genome" and len(np.unique(chroms)) > 1:
        n = len(chroms); return np.roll(np.arange(n), int(rng.integers(int(.1 * n), int(.9 * n))))
    return _shift_within(chroms)

def _shift_within(chroms):
    """Indice che sposta circolarmente i valori DENTRO ogni cromosoma (offset casuale 10-90% della lunghezza).
    I bin sono gia' ordinati per posizione: lo shift conserva l'autocorrelazione del target."""
    idx = np.arange(len(chroms))
    for c in np.unique(chroms):
        m = np.where(chroms == c)[0]; n = len(m)
        k = int(rng.integers(int(.1 * n), max(int(.1 * n) + 1, int(.9 * n))))
        idx[m] = np.roll(m, k)
    return idx

floor_rows, floor_oof, imps_floor = [], [], []
perf, imps, oof, sweep, perf_chr = [], [], [], [], []
units = [[c] for c in chroms] if a.mode == "chrom" else [chroms]          # chrom: un modello per cromosoma
for om in a.omics:
    col, tr_ = TARGETS[om]
    if col not in D.columns: print(f"salto {om}: colonna {col} assente"); continue
    for unit in units:
        d = D[D["chr"].isin(unit)]
        d_all = d
        if om == "mut" or om in MUT_CLASSES:                      # stessi bin per tutti i target delle mutazioni
            d = d[d[a.exonic_col] >= a.min_exonic_bp] if a.exonic_col else d[d["num.overlapping.genes_protein_coding"] > 0]
        if om in MUT_CLASSES:
            sc = MUT_CLASSES[om][1]
            if sc not in d.columns: print(f"salto {om}: colonna {sc} assente"); continue
            n0 = len(d); d = d[d[sc] > 0]                              # stessi bin in tutte le esposizioni (confronto equo)
            if a.class_exposure == "wsites": d = d[d["w" + sc] > 0]
            if len(d) < n0: print(f"  {om} {unit[0] if a.mode == 'chrom' else 'genome'}: {n0 - len(d)} bin senza siti di classe esclusi")
        d = d[d[col].notna()]
        if len(d) < 50 * a.n_splits: continue
        y = d[col].to_numpy(float)
        family, offset = "gaussian", None
        if a.loss == "glm" and om.startswith("cna"):
            family = "binomial"                                  # frazione di pazienti con l'evento nel bin
            if y.min() < 0 or y.max() > 1: sys.exit(f"{col}: la perdita binomiale richiede una frazione in [0, 1]")
        elif a.loss == "glm" and om in MUT_CLASSES:
            family = "poisson"                                   # eventi della classe, esposizione = siti della classe
            expo_col = {"sites": MUT_CLASSES[om][1], "cds": a.exonic_col, "wsites": "w" + MUT_CLASSES[om][1]}[a.class_exposure]
            offset = np.log(d[expo_col].to_numpy(float))
        elif a.loss == "glm" and om == "mut":
            family = "poisson"                                   # conteggio di campioni mutati, esposizione = kb di CDS
            if a.exonic_col: offset = np.log(d[a.exonic_col].to_numpy(float) / 1000.0)
        else:
            if om == "mut" and a.mut_target == "density":
                y = y / (d[a.exonic_col].to_numpy(float) / 1000.0)              # eventi per kb di CDS
            y = np.log1p(y) if tr_ == "log1p" else y
        Xr = d[FEATS]; keep = [c for c in FEATS if Xr[c].nunique() > 1]
        X = Xr[keep].to_numpy(float); ch_arr = d["chr"].to_numpy()
        if a.shadow_mode == "per_feature":
            SH, sh_names = make_shadows(X, keep, ch_arr); Xd = np.column_stack([X, SH]); names = keep + sh_names
        else:
            Xd = np.column_stack([X, make_decoys(X, ch_arr)]); names = keep + [f"DECOY_{i}" for i in range(a.n_decoys)]
        ell = int(np.median([corr_len(y[ch_arr == c]) for c in np.unique(ch_arr)])); buf = min(ell, a.max_buffer_bins)
        label = unit[0] if a.mode == "chrom" else "genome"
        extra = None
        if om in MUT_CLASSES and family == "poisson" and not a.no_predict_low_cds and a.exonic_col:
            sc_ = MUT_CLASSES[om][1]
            de = d_all[(d_all[a.exonic_col] < a.min_exonic_bp) & (d_all[sc_] > 0) & d_all[col].notna()]
            if a.class_exposure == "wsites": de = de[de["w" + sc_] > 0]
            if len(de):
                Xe = de[FEATS][keep].to_numpy(float)
                Xe = np.column_stack([Xe, np.full((len(de), Xd.shape[1] - X.shape[1]), np.nan)])   # esche: mancanti nei bin aggiunti
                extra = (Xe, de["binID"].to_numpy(), np.log(de[expo_col].to_numpy(float)))
        try:
            pred, shap, r2, r2_lin, pred_lin, fold_of, dev = fit_oof(Xd, y, d["binID"].to_numpy(), buf, family=family, offset=offset,
                                                                     init_rate=om in MUT_CLASSES, extra=extra)
        except ValueError as e:
            print(f"  {om} {label}: {e}"); continue
        perf.append(dict(omic=om, unit=label, family=family, n_bins=len(y), ell_bins=ell, buffer_bins=buf, R2=r2, R2_ridge=r2_lin,
                         dev_expl=dev["dev_expl"], dev_expl_lin=dev["dev_expl_lin"]))
        for c, v in r2_by_chrom(y, pred, ch_arr).items():
            perf_chr.append(dict(omic=om, chr=c, R2=v, R2_ridge=r2_by_chrom(y, pred_lin, ch_arr)[c]))
        # importanza: per cromosoma (mode chrom) oppure per fold (mode genome) -> unita' di replica
        rep_ids = np.array([label] * len(y)) if a.mode == "chrom" else fold_of.astype(str)
        for rid in np.unique(rep_ids):
            mabs = np.abs(shap[rep_ids == rid]).mean(axis=0)
            imps.append(pd.DataFrame(dict(omic=om, chr=rid, feature=names, imp=mabs / mabs.sum())))
        oof.append(pd.DataFrame(dict(omic=om, binID=d["binID"].to_numpy(), obs=y, pred=pred,
                                     resid=pearson_resid(y, pred, family), resid_raw=y - pred, pred_ridge=pred_lin, family=family,
                                     in_training=True)))
        if extra is not None:                                  # bin con poca CDS: predetti, non usati per allenare ne' per le metriche
            pe = dev["pred_extra"]; ye = de[col].to_numpy(float); ok_ = np.isfinite(pe)
            oof.append(pd.DataFrame(dict(omic=om, binID=de["binID"].to_numpy()[ok_], obs=ye[ok_], pred=pe[ok_],
                                         resid=pearson_resid(ye[ok_], pe[ok_], family), resid_raw=ye[ok_] - pe[ok_],
                                         pred_ridge=np.nan, family=family, in_training=False)))
            print(f"  {om} {label}: bin con < {a.min_exonic_bp:g} bp di CDS predetti (non di training): {int(ok_.sum())}; "
                  f"eventi osservati {ye[ok_].sum():.0f}, attesi {pe[ok_].sum():.1f} (rapporto {ye[ok_].sum() / max(pe[ok_].sum(), 1e-9):.2f})")
        for rep in range(a.n_floor):                            # PAVIMENTO: target spostato, stessa CV e stesse feature
            sh = circular_shift_index(ch_arr)
            ys = y[sh]; offs = None if offset is None else offset[sh]          # conteggio ed esposizione si spostano insieme
            try:
                ps, shs, r2s, _, _, fold_s, devs = fit_oof(Xd, ys, d["binID"].to_numpy(), buf, want_shap=a.floor_importance,
                                                           want_ridge=False, family=family, offset=offs, init_rate=om in MUT_CLASSES)
            except ValueError as e:
                print(f"  pavimento {om} {label} rep {rep}: {e}"); continue
            floor_rows.append(dict(omic=om, unit=label, rep=rep, family=family, scope=a.floor_scope, R2=r2s, dev_expl=devs["dev_expl"]))
            if a.floor_importance:                              # stesse unita' di replica del modello vero, per ripetizione
                rep_s = np.array([f"{label}_r{rep}"] * len(ys)) if a.mode == "chrom" else np.char.add(f"r{rep}_f", fold_s.astype(str))
                for rid in np.unique(rep_s):
                    mabs = np.abs(shs[rep_s == rid]).mean(axis=0)
                    imps_floor.append(pd.DataFrame(dict(omic=om, chr=rid, feature=names, imp=mabs / mabs.sum())))
            floor_oof.append(pd.DataFrame(dict(omic=om, rep=rep, binID=d["binID"].to_numpy(), src_binID=d["binID"].to_numpy()[sh],
                                               obs=ys, pred=ps, family=family,
                                               expo_kb=np.exp(offs) if offs is not None else np.nan)))
        if a.n_floor:
            fr = pd.DataFrame([r for r in floor_rows if r["omic"] == om and r["unit"] == label])
            if len(fr): print(f"  pavimento {om} {label}: R2 medio {fr['R2'].mean():.3f}, devEspl medio {fr['dev_expl'].mean():.3f} ({len(fr)} rip.)")
        if om.startswith("cna"):                                # curva R2 vs cuscinetto (senza decoy/SHAP/ridge)
            for mb in a.buffer_sweep_mb:
                bb = int(round(mb * 1e6 / a.bin_size))
                try:
                    _, _, r2b, _, _, _, devb = fit_oof(X, y, d["binID"].to_numpy(), bb, want_shap=False, want_ridge=False,
                                                       family=family, offset=offset)
                    r2b = devb["dev_expl"] if family != "gaussian" else r2b      # stessa metrica del run principale
                except ValueError: r2b = np.nan
                sweep.append(dict(omic=om, chr=label, buffer_mb=mb, R2=r2b))
        print(f"[{time.time()-t0:.0f}s] {om} {label} [{family}]: n={len(y)} l={ell} buffer={buf} R2_xgb={r2:.3f} R2_lin={r2_lin:.3f}"
              f" devEspl_xgb={dev['dev_expl']:.3f} devEspl_lin={dev['dev_expl_lin']:.3f}")

pd.DataFrame(perf_chr).to_csv(f"{a.out_dir}/performance_per_chr.tsv", sep="\t", index=False)
if floor_rows:
    FL = pd.DataFrame(floor_rows); FL.to_csv(f"{a.out_dir}/floor.tsv", sep="\t", index=False)
    pd.concat(floor_oof).to_csv(f"{a.out_dir}/floor_oof_predictions.tsv.gz", sep="\t", index=False)
    print("\nPavimento (target spostato circolarmente), media sulle ripetizioni:\n",
          FL.groupby("omic")[["R2", "dev_expl"]].mean().round(3))
perf = pd.DataFrame(perf); imp = pd.concat(imps); oof = pd.concat(oof)
perf.to_csv(f"{a.out_dir}/performance.tsv", sep="\t", index=False)
imp.to_csv(f"{a.out_dir}/importance_per_chr.tsv", sep="\t", index=False)
oof.to_csv(f"{a.out_dir}/oof_predictions.tsv", sep="\t", index=False)
print("\nR2 per omica (mode chrom: mediana tra cromosomi; mode genome: modello unico):\n",
      perf.groupby("omic")[["R2", "R2_ridge", "dev_expl", "dev_expl_lin"]].median().round(3))
if a.mode == "genome":
    pc = pd.DataFrame(perf_chr)
    print("\nGenome-wide, R2 per cromosoma (mediana):\n", pc.groupby("omic")[["R2", "R2_ridge"]].median().round(3))

# ---- curva R2 vs cuscinetto ----------------------------------------------------------------
if sweep:
    sw = pd.DataFrame(sweep); sw.to_csv(f"{a.out_dir}/cna_buffer_sweep.tsv", sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(6, 4))
    for om, g in sw.groupby("omic"):
        m = g.groupby("buffer_mb")["R2"].median(); ax.plot(m.index, m.values, marker="o", label=om)
    ax.axhline(0, color="grey", lw=.8); ax.set_xlabel("cuscinetto (Mb)"); ax.set_ylabel("R2 out-of-fold (mediana cromosomi)")
    ax.set_title("CNA: R2 al crescere del cuscinetto"); ax.legend(); fig.tight_layout(); fig.savefig(f"{a.out_dir}/cna_buffer_sweep.png", dpi=150)
    print("\nCurva R2/cuscinetto (mediana):\n", sw.groupby(["omic", "buffer_mb"])["R2"].median().unstack().round(3))

# ---- gruppi di feature correlate (l'importanza si divide tra feature collineari) ------------------
rho = D[FEATS].rank().corr().abs().fillna(0).to_numpy(copy=True); np.fill_diagonal(rho, 1)
grp = fcluster(linkage(squareform(1 - rho, checks=False), "average"), t=1 - a.group_corr, criterion="distance")
groups = pd.DataFrame(dict(feature=FEATS, group=grp))
lab = groups.groupby("group")["feature"].apply(lambda s: s.iloc[0] + (f" (+{len(s)-1})" if len(s) > 1 else ""))
groups["group_label"] = groups["group"].map(lab)
groups.to_csv(f"{a.out_dir}/feature_groups.tsv", sep="\t", index=False)

# ---- nucleo vs moduli ------------------------------------------------------------------------------
if a.mode == "chrom":
    good = perf[perf["R2"] > 0][["omic", "unit"]].rename(columns={"unit": "chr"})   # solo cromosomi dove il modello impara
    I = imp.merge(good, on=["omic", "chr"])
else:
    I = imp[imp["omic"].isin(perf.loc[perf["R2"] > 0, "omic"])]                 # repliche = fold
IS_DECOY = lambda f: f.str.startswith("DECOY_") | f.str.startswith("SHADOW")
def usage(I):
    """Criterio di uso per gruppo. 'used' = come prima: somma dell'importanza del gruppo > quantile --used_quantile delle
    esche, in almeno --used_frac_chroms delle repliche. Con --shadow_mode per_feature anche 'used_shadow': importanza del
    gruppo > massimo, sulle copie j, della somma delle copie spostate j dei SUOI membri (confronto a pari dimensione)."""
    thr = I[IS_DECOY(I["feature"])].groupby(["omic", "chr"])["imp"].quantile(a.used_quantile).rename("thr")
    G_ = I[~IS_DECOY(I["feature"])].merge(groups, on="feature")
    G_ = G_.groupby(["omic", "chr", "group_label"], as_index=False)["imp"].sum().merge(thr.reset_index(), on=["omic", "chr"])
    G_["used"] = G_["imp"] > G_["thr"]
    agg = dict(frac_chroms_used=("used", "mean"), mean_imp=("imp", "mean"))
    if a.shadow_mode == "per_feature":
        Sh = I[I["feature"].str.startswith("SHADOW")].copy()
        Sh["copy"] = Sh["feature"].str.extract(r"^SHADOW(\d+)__")[0]; Sh["feature"] = Sh["feature"].str.replace(r"^SHADOW\d+__", "", regex=True)
        Sh = Sh.merge(groups, on="feature").groupby(["omic", "chr", "group_label", "copy"])["imp"].sum()
        thr_g = Sh.groupby(level=["omic", "chr", "group_label"]).max().rename("thr_shadow").reset_index()
        G_ = G_.merge(thr_g, on=["omic", "chr", "group_label"], how="left")
        G_["used_shadow"] = G_["imp"] > G_["thr_shadow"]; agg["frac_used_shadow"] = ("used_shadow", "mean")
    S_ = G_.groupby(["omic", "group_label"]).agg(**agg).reset_index()
    S_["used"] = S_["frac_chroms_used"] >= a.used_frac_chroms
    if a.shadow_mode == "per_feature": S_["used_shadow"] = S_["frac_used_shadow"] >= a.used_frac_chroms
    return G_, S_
G, S = usage(I)
if a.floor_importance and imps_floor:
    _, SF = usage(pd.concat(imps_floor))
    keepc = [c for c in ("frac_chroms_used", "mean_imp", "used", "frac_used_shadow", "used_shadow") if c in SF.columns]
    S = S.merge(SF[["omic", "group_label"] + keepc].rename(columns={c: c + "_floor" for c in keepc}), on=["omic", "group_label"], how="left")
    for c in ("used", "used_shadow"):
        if c in S.columns:                                     # robusto = usato col target vero e NON col target spostato
            S[c + "_robust"] = S[c] & ~S[c + "_floor"].fillna(False).astype(bool)
    S["imp_ratio_real_floor"] = S["mean_imp"] / S["mean_imp_floor"]
    pd.concat(imps_floor).to_csv(f"{a.out_dir}/importance_floor_per_rep.tsv", sep="\t", index=False)
    print("\nUso delle feature anche nel pavimento (target spostato) - gruppi usati nel modello vero:")
    cols_ = ["omic", "group_label", "used", "used_floor"] + (["used_shadow", "used_shadow_floor"] if "used_shadow" in S.columns else []) + ["imp_ratio_real_floor"]
    print(S[S["used"]][cols_].sort_values(["omic", "imp_ratio_real_floor"]).round(2).to_string(index=False))
W = S.pivot(index="group_label", columns="omic", values="used")
missing = [o for o in a.omics if o not in W.columns]
if missing: print(f"ATTENZIONE: omiche senza repliche con R2>0, contate come 'non usano nulla': {missing}")
W = W.reindex(columns=[o for o in a.omics if o in perf["omic"].unique()]).fillna(False).astype(bool)   # nucleo = TUTTE le omiche richieste
oms = list(W.columns)
W["n_omics"] = W[oms].sum(axis=1)
W["category"] = np.where(W["n_omics"] == len(oms), "core (tutte le omiche)",
                 np.where(W["n_omics"] == 1, "modulo " + W[oms].idxmax(axis=1).astype(str),
                 np.where(W["n_omics"] == 0, "non usata", "condivisa parziale")))
W.to_csv(f"{a.out_dir}/core_vs_modules.tsv", sep="\t")
S.to_csv(f"{a.out_dir}/usage_per_omic.tsv", sep="\t", index=False)
print("\nCategorie di gruppi di feature:\n", W["category"].value_counts())
for cat in sorted(W["category"].unique()):
    if cat != "non usata": print(f"  {cat}: {', '.join(W.index[W['category'] == cat][:15])}")

# somiglianza tra omiche dei profili di importanza (Spearman, mediana sui cromosomi)
P = G.pivot_table(index=["chr", "group_label"], columns="omic", values="imp")
sim = pd.concat([g.droplevel(0).corr(method="spearman").stack().rename(ch) for ch, g in P.groupby(level=0)], axis=1).median(axis=1).unstack()
sim.to_csv(f"{a.out_dir}/importance_similarity.tsv", sep="\t"); print("\nSomiglianza dei profili di importanza (Spearman mediano):\n", sim.round(2))

# heatmap: gruppi usati da almeno un'omica x omiche
M = S.pivot(index="group_label", columns="omic", values="mean_imp").reindex(columns=oms).fillna(0).loc[W.index[W["n_omics"] > 0]]
U = W.loc[M.index, oms]
order = W.loc[M.index].sort_values(["n_omics"] + oms, ascending=False).index; M, U = M.loc[order, oms], U.loc[order, oms]
fig, ax = plt.subplots(figsize=(1.2 * len(oms) + 4, 0.28 * len(M) + 1.5))
im = ax.imshow(M.to_numpy(), aspect="auto", cmap="viridis")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        if U.iat[i, j]: ax.text(j, i, "*", ha="center", va="center", color="white", fontsize=9)
ax.set_xticks(range(len(oms))); ax.set_xticklabels(oms); ax.set_yticks(range(len(M))); ax.set_yticklabels(M.index, fontsize=7)
ax.set_title("Importanza relativa (|SHAP| OOF) per gruppo di feature\n* = sopra i decoy nella maggioranza dei cromosomi")
fig.colorbar(im, ax=ax, fraction=.03); fig.tight_layout(); fig.savefig(f"{a.out_dir}/core_vs_modules_heatmap.png", dpi=150)
print(f"\nFatto in {time.time()-t0:.0f}s. Output in {a.out_dir}/")