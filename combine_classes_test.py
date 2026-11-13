#!/usr/bin/env python
"""
BURNER - Combinare missense e troncanti invece di sommarle: e' qui il divario con dNdScv?
=========================================================================================
Legge i file per gene gia' prodotti da gene_level_drivers_classes.py (nessun rilancio) e calcola, per ogni gene, un
punteggio NON SINONIMO combinando le due classi in tre modi; valuta l'AP sui driver nello stesso universo del confronto
con dNdScv (geni presenti nel file di dNdScv), accanto a: somma delle classi (il nostro "nonsyn" di oggi) e dNdScv.
Metodi (p = coda superiore della binomiale negativa per classe, cioe' 10^-punteggio dei file):
  fisher   -2 (ln p_mis + ln p_trunc) ~ chi2 con 4 gradi di liberta'. Premia evidenza in entrambe le classi.
  minp     min(p_mis, p_trunc) * 2 (Bonferroni). Premia evidenza forte in UNA classe (es. oncosoppressori).
  lrt      rapporto di verosimiglianza alla dNdScv con il NOSTRO atteso: per classe c, O_c ~ NB(E_c * omega_c, theta_c);
           sotto H0 omega_c = 1; sotto H1 omega_c libero, stimato per massima verosimiglianza (per un'osservazione NB
           con theta fisso il massimo e' a media = O, cioe' omega = O/E). Solo selezione POSITIVA: se O_c <= E_c il
           contributo della classe e' 0. Statistica = somma sulle due classi di 2[l(O; O) - l(O; E)]; p da chi2 con 2
           gradi di liberta' (approssimazione: ignora che omega >= 1 e' un confine, quindi e' prudente).
           theta_c: metodo dei momenti come nello script dei driver (1% dei geni piu' estremi escluso).
Attesi: "model" (atteso del nostro background) per tutti i metodi; "model+syn" solo per fisher e minp, perche' i file
salvano il punteggio di model+syn ma non il suo atteso corretto (per lrt con +syn va salvato l'atteso nello script dei driver).
Uso: python combine_classes_test.py --ours_dir <...>/drivers_gene_level --dndscv <sel_cv_...tsv> --label "BRCA 20kb" --out_dir DIR
"""
import argparse, os
import numpy as np, pandas as pd
from scipy import stats, special
from sklearn.metrics import average_precision_score, roc_auc_score
p = argparse.ArgumentParser()
p.add_argument("--ours_dir", required=True); p.add_argument("--dndscv", required=True)
p.add_argument("--label", default=""); p.add_argument("--out_dir", default="combine_classes")
a = p.parse_args(); os.makedirs(a.out_dir, exist_ok=True)

def theta_mom(O, E, trim=0.01):
    """Come theta_mom di gene_level_drivers_classes.py: momenti, escluso l'1% di geni con z di Poisson piu' alto."""
    O = np.asarray(O, float); E = np.asarray(E, float); ok = np.isfinite(E) & (E > 0); O, E = O[ok], E[ok]
    z = (O - E) / np.sqrt(E); keep = z <= np.quantile(z, 1 - trim); O, E = O[keep], E[keep]
    den = np.sum((O - E) ** 2) - np.sum(O)
    return np.inf if den <= 0 else np.sum(E ** 2) / den

def nb_logpmf(o, mu, theta):
    if not np.isfinite(theta): return stats.poisson.logpmf(o, mu)
    return (special.gammaln(o + theta) - special.gammaln(theta) - special.gammaln(o + 1)
            + theta * np.log(theta / (theta + mu)) + o * np.log(mu / (theta + mu)))

G = {c: pd.read_csv(f"{a.ours_dir}/gene_scores_{c}.tsv", sep="\t", index_col=0) for c in ("mis", "trunc", "nonsyn")}
VAR = {c: pd.read_csv(f"{a.ours_dir}/gene_scores_variants_{c}.tsv", sep="\t", index_col=0) for c in ("mis", "trunc", "nonsyn")}
D = pd.read_csv(a.dndscv, sep="\t").drop_duplicates("gene_name").set_index("gene_name")
U = G["nonsyn"].index.intersection(D.index)                               # universo del confronto con dNdScv
drv = G["nonsyn"].loc[U, "is_driver"].astype(bool).to_numpy()

def pval(score): return np.power(10.0, -score.reindex(U).astype(float))      # punteggio -> p (NaN se senza atteso)
S = {}
for base, src in (("model", lambda c: G[c]["residual"]), ("model+syn", lambda c: VAR[c]["model+syn"])):
    pm, pt = pval(src("mis")), pval(src("trunc"))
    both = pm.notna() & pt.notna()
    fis = pd.Series(np.nan, index=U); fis[both] = stats.chi2.sf(-2 * (np.log(pm[both].clip(lower=1e-300)) + np.log(pt[both].clip(lower=1e-300))), 4)
    mnp = (np.minimum(pm, pt) * 2).clip(upper=1)
    S[f"fisher ({base})"] = -np.log10(fis.clip(lower=1e-300)); S[f"minp ({base})"] = -np.log10(mnp.clip(lower=1e-300))
    S[f"sum of classes ({base})"] = (G["nonsyn"]["residual"] if base == "model" else VAR["nonsyn"]["model+syn"]).reindex(U)
# LRT con l'atteso del modello
lrt = pd.Series(0.0, index=U); has = pd.Series(True, index=U)
for c in ("mis", "trunc"):
    O = G[c]["O"].reindex(U).astype(float); E = G[c]["E"].reindex(U).astype(float)
    th = theta_mom(G[c]["O"], G[c]["E"]); print(f"theta {c} (momenti, come lo script dei driver): {th:.3g}")
    ok = E.notna() & (E > 0); has &= ok
    o, e = O[ok].to_numpy(), E[ok].to_numpy(); pos = o > e
    stat = np.zeros(len(o)); stat[pos] = 2 * (nb_logpmf(o[pos], o[pos], th) - nb_logpmf(o[pos], e[pos], th))
    lrt[ok] += np.maximum(stat, 0)
lrt[~has] = np.nan
S["lrt (model)"] = -np.log10(pd.Series(stats.chi2.sf(lrt, 2), index=U).clip(lower=1e-300)).where(has)
S["dNdScv"] = -np.log10(D.loc[U, "pallsubs_cv"].astype(float).clip(lower=1e-300))

rows = []
for name, s in S.items():
    x = s.reindex(U).astype(float).fillna(-1).to_numpy()                  # senza punteggio: in fondo
    rows.append(dict(method=name, AP=average_precision_score(drv, x), AUROC=roc_auc_score(drv, x),
                     drivers_top20=int(drv[np.argsort(-x)[:20]].sum()), drivers_top100=int(drv[np.argsort(-x)[:100]].sum())))
T = pd.DataFrame(rows)
print(f"\n{a.label}: universo {len(U)} geni, driver {drv.sum()}")
print(T.round(3).to_string(index=False))
T.to_csv(f"{a.out_dir}/combine_classes_{a.label.replace(' ', '_')}.tsv", sep="\t", index=False)
pd.DataFrame(S).assign(is_driver=drv).to_csv(f"{a.out_dir}/gene_scores_combined_{a.label.replace(' ', '_')}.tsv", sep="\t")