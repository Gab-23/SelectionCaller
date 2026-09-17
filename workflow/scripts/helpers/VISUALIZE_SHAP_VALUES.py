# workflow/scripts/helpers/VISUALIZE_SHAP_VALUES.py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def VISUALIZE_SHAP_VALUES(SHAP_VALUES, BIN_SIZE_BP, OUTPUT_PATH, ROLLING_WINDOW_BINS=50):
    # SHAP_VALUES: DataFrame indicizzato per binID ("{chrom}_{bin_num}", stessa convenzione
    # di BLOCK_DATASET.py), colonne = nomi feature. Una pagina per feature in un unico PDF
    # multi-pagina - coerente con come il resto della pipeline salva i plot (un solo output
    # per regola, mai una directory con conteggio di file variabile e ignoto a parse-time).
    #
    # A 20kb di risoluzione un cromosoma puo' avere decine di migliaia di bin - una linea
    # per-bin grezza sarebbe illeggibile, quindi ogni pannello sovrappone il segnale grezzo
    # (sottile, alpha basso) a una media mobile (ROLLING_WINDOW_BINS bin) - lo smoothing e'
    # solo per la visualizzazione, non tocca mai i dati salvati in SHAP_VALUES.

    bin_nums = np.array([int(b.split("_")[1]) for b in SHAP_VALUES.index])
    order = np.argsort(bin_nums)  # il dataset e' gia' chr-sorted (assunzione di BLOCK_DATASET),
                                   # ma riordiniamo qui per sicurezza: non e' garantito che
                                   # out_of_fold_shap_values mantenga l'ordine originale dopo
                                   # il pd.concat dei fold
    positions_mb = (bin_nums[order] * BIN_SIZE_BP) / 1e6

    with PdfPages(OUTPUT_PATH) as pdf:
        for feature in SHAP_VALUES.columns:
            values = SHAP_VALUES[feature].values[order]
            smoothed = pd.Series(values).rolling(
                ROLLING_WINDOW_BINS, center=True, min_periods=1
            ).mean()

            fig, ax = plt.subplots(figsize=(14, 3))
            ax.plot(positions_mb, values, color="lightgray", linewidth=0.5, alpha=0.6, label="raw")
            ax.plot(positions_mb, smoothed, color="black", linewidth=1.2,
                     label=f"rolling mean ({ROLLING_WINDOW_BINS} bins)")
            ax.axhline(0, color="red", linestyle="--", linewidth=0.8)
            ax.set_xlabel("Genomic position (Mb)")
            ax.set_ylabel("SHAP value")
            ax.set_title(feature)
            ax.legend(loc="upper right", fontsize="small")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)