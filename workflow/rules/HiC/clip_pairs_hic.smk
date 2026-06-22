rule clip_pairs_hic:
    output:
        pairs="resources/hic/files/pairs/{SRA_id}/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}.filtered.dedup.sorted.clipped.pairs"
    input:
        pairs="resources/hic/files/pairs/{SRA_id}/{assembly_code}_{mapq_thr}.filtered.dedup.sorted.pairs",
        actual_chrom_sizes="resources/reference_data/chrom_sizes/{assembly_code}_{window_mode}_{window_size}_{step}.actual.chrom.sizes" 
    log:
        "logs/clip_pairs_hic/{assembly_code}_{window_mode}_{window_size}_{step}_{mapq_thr}_{SRA_id}.log"
    shell:
        """
        awk '
            BEGIN {{
                while ((getline line < "{input.actual_chrom_sizes}") > 0) {{
                    split(line, a, "\t")
                    size[a[1]] = a[2]
                }}
            }}
            /^#/ {{ print; next }}
            ($2 in size) && ($4 in size) && ($3 <= size[$2]) && ($5 <= size[$4]) {{ print }}
        ' {input.pairs} > {output.pairs} 2>{log}
        """