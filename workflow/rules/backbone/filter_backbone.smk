# Bin esclusi se la FRAZIONE della loro lunghezza coperta dalla blacklist (intervalli fusi) e' >= blacklist_max_fraction.
# Prima: esclusione per qualunque sovrapposizione (bedtools intersect -v), che con bin grandi faceva perdere interi
# megabasi per piccole regioni in blacklist (a 1 Mb: 82% dei bin tenuti contro il 93% a 20 kb; persi driver come AKT1,
# FOXA1, KMT2C). Le frazioni sono bimodali (quasi tutte < 10% o > 90%): la soglia 0.5 ("bin per la maggior parte pulito")
# sta nella valle tra i due picchi; scelta senza guardare la lista dei driver.
# Parametro solo da config (non entra nei nomi dei file): dopo un cambio usare --forcerun filter_backbone.
BLACKLIST_MAX_FRACTION = config.get("blacklist_max_fraction", ["0.5"])
if len(BLACKLIST_MAX_FRACTION) != 1 or not (0 < float(BLACKLIST_MAX_FRACTION[0]) <= 1):
    raise ValueError("blacklist_max_fraction is a config-only parameter: give exactly one value in (0, 1]")

rule filter_backbone:
    output:
        backbone_full="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.bed",
        backbone_mask="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.mask.filtered.bed"
    input:
        backbone_full="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.bed",
        backbone_mask="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.mask.bed",
        blacklist="resources/blacklist/{assembly_code}.bed",
        chrom_sizes="resources/reference_data/chrom_sizes/{assembly_code}.chrom.sizes"
    params:
        max_frac=str(BLACKLIST_MAX_FRACTION[0])
    log:
        "logs/filter_backbone/{assembly_code}_{window_mode}_{window_size}_{step}.log"
    conda:
        "../../envs/filter_backbone.yml"
    shell:
        """
        set -euo pipefail
        merged=$(mktemp)
        trap 'rm -f "$merged"' EXIT
        # intervalli della blacklist fusi (stesso ordine dei cromosomi del file, gia' ordinato come richiesto da -sorted)
        bedtools merge -i {input.blacklist} > "$merged" 2>>{log}
        for io in "{input.backbone_full}:{output.backbone_full}" "{input.backbone_mask}:{output.backbone_mask}"; do
            in="${{io%%:*}}"; out="${{io##*:}}"
            # bedtools coverage aggiunge 4 colonne: n. intervalli, basi coperte, lunghezza del bin, frazione coperta
            bedtools coverage -a "$in" -b "$merged" -g {input.chrom_sizes} -sorted 2>>{log} \
              | awk -v t={params.max_frac} 'BEGIN{{OFS="\\t"}} $NF < t {{NF=NF-4; print}}' > "$out"
            echo "$(basename "$in"): bin $(wc -l < "$in"), tenuti $(wc -l < "$out") (frazione in blacklist < {params.max_frac})" >> {log}
        done
        """