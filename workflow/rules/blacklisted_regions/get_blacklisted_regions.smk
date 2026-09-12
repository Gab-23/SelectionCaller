rule get_blacklisted_regions:
    output:
        "resources/blacklist/{assembly_code}.bed"
    input:
        "resources/reference_data/chrom_sizes/{assembly_code}.chrom.sizes"
    log:
        "logs/get_blacklisted_regions/{assembly_code}.log"
    conda:
        "../../envs/filter_backbone.yml"
    shell:
        # NOTE (fix): homologated with get_reference_genome's download-robustness fix - the
        # original piped `wget -qO- ... | gunzip | ...` directly, which both hides wget's own
        # errors (-q) and silently produces a truncated/corrupt intermediate if the download
        # is incomplete, with no retry. Now downloads to a real temp file first (with retry/
        # timeout/continue), then decompresses as a separate, clearly-attributable step.
        """
        TMP_GZ=$(mktemp --suffix=.bed.gz)
        wget --tries=5 --timeout=60 --continue -O "$TMP_GZ" \
        https://github.com/Boyle-Lab/Blacklist/raw/refs/heads/master/lists/{wildcards.assembly_code}-blacklist.v2.bed.gz \
        2>{log}
        gunzip -c "$TMP_GZ" | bedtools sort -g {input} > {output} 2>>{log}
        rm -f "$TMP_GZ"
        """
