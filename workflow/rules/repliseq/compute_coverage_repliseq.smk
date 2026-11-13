rule compute_coverage_repliseq:
    # Conteggio delle letture per bin.
    # FIX (27/10/2026, OOM): stesso problema e stessa soluzione di compute_coverage_atacseq (21/10/2026). Senza -sorted,
    # bedtools coverage carica in memoria tutte le letture del BAM; con piu' campioni in parallelo il job della pipeline
    # e' arrivato a 157 GB ed e' stato ucciso (OUT_OF_MEMORY) proprio su questa regola. Con -sorted (chromsweep) la
    # memoria resta di pochi MB.
    # -sorted richiede che backbone e BAM seguano lo stesso ordine dei cromosomi: il backbone e' ordinato in modo
    # lessicografico (bedtools sort), il BAM nell'ordine del suo header. Quindi: 1) file genome dall'header del BAM,
    # 2) backbone riordinato su quell'ordine, 3) conteggio con -sorted, 4) output riordinato in modo lessicografico,
    # cosi' il file finale e' IDENTICO (contenuto e ordine delle righe) a quello della versione precedente.
    output:
        "resources/repliseq/files/coverage/{SRA_id}/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}.bed"
    input:
        backbone="resources/backbone/{assembly_code}_{window_mode}_{window_size}_{step}.filtered.sorted.bed",
        bam="resources/repliseq/files/bam/{SRA_id}/{assembly_code}_{mapq_thr}.noblack.sorted.bam"
    log:
        "logs/compute_coverage_repliseq/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{SRA_id}.log"
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