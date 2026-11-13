rule compute_coverage_atacseq:
    # Conteggio delle letture per bin.
    # FIX (21/10/2026, OOM): senza -sorted, bedtools coverage carica in memoria tutte le letture
    # del BAM (misurato: ~1 GB per milione di letture), quindi un BAM ATAC da decine di milioni
    # di letture richiede decine di GB, moltiplicati per i job in parallelo -> OOM a 128 GB.
    # Con -sorted (chromsweep) la memoria resta di pochi MB.
    # -sorted richiede che backbone e BAM seguano lo stesso ordine dei cromosomi: il backbone e'
    # ordinato in modo lessicografico (bedtools sort), il BAM nell'ordine del suo header. Quindi:
    # 1) si ricava il file genome dall'header del BAM, 2) si riordina il backbone su quell'ordine,
    # 3) si conta con -sorted, 4) si riordina l'output in modo lessicografico, cosi' il file finale
    # e' IDENTICO (contenuto e ordine delle righe) a quello prodotto dalla versione precedente.
    output:
        "resources/atacseq/files/coverage/{SRA_id}/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}.bed"
    input:
        backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.mask.filtered.sorted.bed",
        bam="resources/atacseq/files/bam/{SRA_id}/{SRA_id}_{assembly_code}_{mapq_thr}.fixed.bam"
    log:
        "logs/compute_coverage_atacseq/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{SRA_id}.log"
    conda:
        "../../envs/bedtools_samtools.yml"
    shell:
        """
        set -euo pipefail
        TMP=$(mktemp -d)
        trap 'rm -rf "$TMP"' EXIT

        samtools view -H {input.bam} | \
        awk -F'\\t' '$1=="@SQ"{{n="";l=""; for(i=2;i<=NF;i++){{if($i~/^SN:/)n=substr($i,4); if($i~/^LN:/)l=substr($i,4)}}; print n"\\t"l}}' \
        > "$TMP/genome" 2>{log}

        bedtools sort -i {input.backbone} -g "$TMP/genome" > "$TMP/backbone.bed" 2>>{log}

        bedtools coverage -sorted -g "$TMP/genome" -counts -a "$TMP/backbone.bed" -b {input.bam} 2>>{log} | \
        bedtools sort -i - > {output} 2>>{log}
        """