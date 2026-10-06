import numpy as np
import pandas as pd
import math
import os
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse

from workflow.scripts.helpers.SMOOTH_SIGNAL import *
from workflow.scripts.helpers.STRINGENCY_TO_PROMINENCE import *

parser = argparse.ArgumentParser()

parser.add_argument('--input_path_shap_values', dest = 'input_path_shap_values', required = True, type = str, help = '')
parser.add_argument('--stringency', dest = 'stringency', required = True, type = float, help = '')
parser.add_argument('--apply_smoothing', dest = 'apply_smoothing', required = False, type = bool, default = True, help = '')
parser.add_argument('--smoothing_size', dest = 'smoothing_size', required = True, type = int, help = '')
parser.add_argument('--peak_buffer_bins', dest = 'peak_buffer_bins', required = True, type = int, help = '')
parser.add_argument('--importance_aggregation', dest = 'importance_aggregation', required = True, type = str, help = '')
parser.add_argument('--bin_size_bp', dest = 'bin_size_bp', required = True, type = int, help = '')
parser.add_argument('--min_importance', dest = 'min_importance', required = False, type = float, default = None, help = '')
parser.add_argument('--output_path_feature_locality', dest = 'output_path_feature_locality', required = True, type = str, help = '')
parser.add_argument('--output_path_landscape', dest = 'output_path_landscape', required = True, type = str, help = '')
parser.add_argument('--output_path_summary', dest = 'output_path_summary', required = True, type = str, help = '')

args = parser.parse_args()

shap_values = pd.read_csv(args.input_path_shap_values, sep="\t", header=0, index_col=0)
distance = max(1, args.smoothing_size)

bin_nums = np.array([int(b.split("_")[1]) for b in shap_values.index])
order = np.argsort(bin_nums)
positions_mb = (bin_nums[order] * args.bin_size_bp) / 1e6

for path in [args.output_path_feature_locality, args.output_path_landscape, args.output_path_summary]:
    os.makedirs(os.path.dirname(path), exist_ok=True)

feature_locality_dict = {}
signals_for_plot = {}

for feature in shap_values.columns:
    x_abs_raw = np.abs(shap_values[feature].values[order])
    x_abs = x_abs_raw
    if args.apply_smoothing:
        x_abs = SMOOTH_SIGNAL(values = x_abs, smoothing_size = args.smoothing_size)

    aggregated_x_abs = pd.Series(x_abs).aggregate(func = args.importance_aggregation)

    K, prominence = STRINGENCY_TO_PROMINENCE(x_abs, args.stringency)
    peaks, properties = find_peaks(x = x_abs, prominence = prominence, distance = distance)
    n_peaks = len(peaks)

    if n_peaks > 0:
        peak_mask = np.zeros(len(x_abs), dtype=bool)
        for p in peaks:
            lo = max(0, p - args.peak_buffer_bins)
            hi = min(len(x_abs), p + args.peak_buffer_bins + 1)
            peak_mask[lo:hi] = True
        peak_importance = x_abs[peak_mask].mean()
        background_importance = x_abs[~peak_mask].mean()
    else:
        peak_importance = 0.0
        background_importance = x_abs.mean()

    feature_locality_dict[feature] = {
        "importance": aggregated_x_abs,
        "n_peaks": n_peaks,
        "peak_importance": peak_importance,
        "background_importance": background_importance,
    }
    signals_for_plot[feature] = {"x_abs_raw": x_abs_raw, "x_abs": x_abs, "peaks": peaks}

feature_locality = pd.DataFrame.from_dict(feature_locality_dict, orient="index")

min_importance = args.min_importance
if min_importance is None:
    min_importance = feature_locality["importance"].quantile(0.25)
    print(f"--min_importance not specified, using the first quartile computed from this run's data: {min_importance}")

below_cutoff = feature_locality["importance"] < min_importance
valid = ~below_cutoff

feature_locality["locality_score"] = np.nan
for feature in feature_locality.index[valid]:
    peak_importance = feature_locality.loc[feature, "peak_importance"]
    background_importance = feature_locality.loc[feature, "background_importance"]
    if background_importance == 0:
        feature_locality.loc[feature, "locality_score"] = float("nan")
        continue
    ratio = peak_importance / background_importance
    if ratio - 1 <= 0:
        feature_locality.loc[feature, "locality_score"] = float("nan")
    else:
        feature_locality.loc[feature, "locality_score"] = math.log2(ratio - 1)

feature_locality.to_csv(args.output_path_feature_locality, sep="\t", header=True, index=True)

with PdfPages(args.output_path_landscape) as pdf:
    for feature in shap_values.columns:
        x_abs_raw = signals_for_plot[feature]["x_abs_raw"]
        x_abs = signals_for_plot[feature]["x_abs"]
        peaks = signals_for_plot[feature]["peaks"]
        n_peaks = len(peaks)

        row = feature_locality.loc[feature]
        loc_str = f"{row['locality_score']:.3f}" if pd.notna(row["locality_score"]) else "NaN"
        title_stats = f"{feature} | importance={row['importance']:.2e}, n_peaks={n_peaks}, locality_score={loc_str}"

        fig, ax = plt.subplots(figsize=(14, 3))
        ax.plot(positions_mb, x_abs_raw, color="lightgray", linewidth=0.5, alpha=0.6, label="raw")
        if args.apply_smoothing:
            ax.plot(positions_mb, x_abs, color="black", linewidth=1.2, label=f"rolling mean ({args.smoothing_size} bins)")
        if n_peaks > 0:
            ax.scatter(positions_mb[peaks], x_abs[peaks], color="red", s=40, zorder=5, label=f"{n_peaks} detected peak(s)")
            for p in peaks:
                lo = max(0, p - args.peak_buffer_bins)
                hi = min(len(x_abs), p + args.peak_buffer_bins + 1)
                ax.axvspan(positions_mb[lo], positions_mb[hi-1], color="red", alpha=0.08)
        ax.axhline(0, color="red", linestyle="--", linewidth=0.5, alpha=0.5)
        ax.set_xlabel("Genomic position (Mb)")
        ax.set_ylabel("|SHAP value|")
        ax.set_title(title_stats, fontsize=9)
        ax.legend(loc="upper right", fontsize="small")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

df_valid = feature_locality.dropna(subset=["locality_score"]).copy()
df_sorted = df_valid.sort_values("locality_score", ascending=True)

fig, axes = plt.subplots(1, 2, figsize=(18, max(6, 0.28 * len(df_valid))))
ax = axes[0]
sc = ax.scatter(df_valid["background_importance"], df_valid["peak_importance"], c=df_valid["locality_score"],
                 cmap="coolwarm", s=80, edgecolor="black", linewidth=0.5)
lims = [min(df_valid["background_importance"].min(), df_valid["peak_importance"].min()),
        max(df_valid["background_importance"].max(), df_valid["peak_importance"].max())]
ax.plot(lims, lims, color="gray", linestyle=":", linewidth=1, label="peak = background")
for feat, r in df_valid.iterrows():
    ax.annotate(feat, (r["background_importance"], r["peak_importance"]), fontsize=6, alpha=0.8, xytext=(3, 3), textcoords="offset points")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("background importance (log scale)")
ax.set_ylabel("peak importance (log scale)")
ax.set_title(f"peak vs background (color = locality_score)\n{below_cutoff.sum()} excluded below min_importance={min_importance:.2e}, buffer={args.peak_buffer_bins} bins")
ax.legend(fontsize=7, loc="upper left")
plt.colorbar(sc, ax=ax, label="locality_score")

ax2 = axes[1]
vmax = max(abs(df_sorted["locality_score"].min()), abs(df_sorted["locality_score"].max()))
norm = plt.Normalize(vmin=-vmax, vmax=vmax)
ax2.barh(df_sorted.index, df_sorted["locality_score"], color=plt.cm.coolwarm(norm(df_sorted["locality_score"])))
ax2.axvline(0, color="black", linewidth=0.8)
ax2.set_xlabel("locality_score = log2(ratio - 1)")
ax2.set_title("Feature ranking by locality_score")
ax2.tick_params(axis="y", labelsize=7)

fig.tight_layout()
fig.savefig(args.output_path_summary, dpi=130)
plt.close(fig)