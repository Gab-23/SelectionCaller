import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from workflow.scripts.helpers.SMOOTH_SIGNAL import SMOOTH_SIGNAL


def VISUALIZE_SHAP_VALUES(SHAP_VALUES, RAW_VALUES, BIN_SIZE_BP, OUTPUT_PATH, ROLLING_WINDOW_BINS=50):
    # SHAP_VALUES: DataFrame indicizzato per binID ("{chrom}_{bin_num}"), colonne = feature.
    # RAW_VALUES: DataFrame indicizzato per binID, con (almeno) le stesse colonne feature di
    # SHAP_VALUES, contenente i valori GREZZI (non SHAP) di quelle feature - tipicamente
    # DATASET.set_index("binID") nel chiamante. Puo' avere piu' colonne/righe di SHAP_VALUES
    # (es. feature non di consenso), non e' un problema: si accede per nome/indice, non per
    # posizione.
    #
    # Una pagina per feature, due pannelli: (1) traccia genomica di |SHAP| (dove conta la
    # feature), (2) scatter SHAP-con-segno vs valore grezzo (se la relazione e' quella
    # biologicamente attesa - es. monotona - non solo dove conta ma se conta nel verso giusto).

    bin_nums = np.array([int(b.split("_")[1]) for b in SHAP_VALUES.index])
    order = np.argsort(bin_nums)
    positions_mb = (bin_nums[order] * BIN_SIZE_BP) / 1e6
    ordered_binIDs = SHAP_VALUES.index[order]

    with PdfPages(OUTPUT_PATH) as pdf:
        for feature in SHAP_VALUES.columns:
            signed_values = SHAP_VALUES[feature].values[order]
            abs_values = np.abs(signed_values)
            smoothed = SMOOTH_SIGNAL(abs_values, ROLLING_WINDOW_BINS)

            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 4), gridspec_kw={"width_ratios": [2.5, 1]})

            # panel 1: genomic landscape - ALWAYS |SHAP| (locality, not direction)
            ax1.plot(positions_mb, abs_values, color="lightgray", linewidth=0.5, alpha=0.6, label="raw")
            ax1.plot(positions_mb, smoothed, color="black", linewidth=1.2,
                      label=f"rolling mean ({ROLLING_WINDOW_BINS} bins)")
            ax1.axhline(0, color="red", linestyle="--", linewidth=0.8)
            ax1.set_xlabel("Genomic position (Mb)")
            ax1.set_ylabel("|SHAP value|")
            ax1.set_title(feature)
            ax1.legend(loc="upper right", fontsize="small")

            # panel 2: scatter - SIGNED SHAP (direction of the effect)
            if feature in RAW_VALUES.columns:
                raw_feature_values = RAW_VALUES.loc[ordered_binIDs, feature].values
                valid = ~(pd.isna(raw_feature_values) | pd.isna(signed_values))
                ax2.scatter(raw_feature_values[valid], signed_values[valid], s=3, alpha=0.3, color="steelblue")
                ax2.axhline(0, color="red", linestyle="--", linewidth=0.8)
                ax2.set_xlabel(f"{feature} (raw value)")
                ax2.set_ylabel("SHAP value (signed)")
                ax2.set_title("SHAP vs raw value")
            else:
                ax2.text(0.5, 0.5, "raw value not available\n(feature missing from RAW_VALUES)",
                          ha="center", va="center", fontsize=9)
                ax2.set_xticks([]); ax2.set_yticks([])

            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)