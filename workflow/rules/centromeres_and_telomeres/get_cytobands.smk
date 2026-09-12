rule get_cytobands:
    output:
        "resources/reference_data/cytobands/{assembly_code}_cytoBand.txt"
    log:
        "logs/get_cytobands/{assembly_code}.log"
    shell:
        # keep only chr1..chr22 (drop chrX/chrY/chrM and alt/random/patch contigs, whose
        # names contain "_"). Output is plain text (.txt, no .gz): the original command
        # already decompressed via gunzip before writing, the ".gz" in the filename was
        # just misleading, not an actual encoding - dropped rather than re-compressed.
        # NOTE (fix): homologated with get_reference_genome's download-robustness fix - the
        # original piped `wget -qO- ... | gunzip | ...` directly, which both hides wget's own
        # errors (-q) and silently produces a truncated/corrupt intermediate if the download
        # is incomplete, with no retry. Now downloads to a real temp file first (with retry/
        # timeout/continue), then decompresses as a separate, clearly-attributable step.
        """
        TMP_GZ=$(mktemp --suffix=.txt.gz)
        wget --tries=5 --timeout=60 --continue -O "$TMP_GZ" \
        https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.assembly_code}/database/cytoBand.txt.gz \
        2>{log}
        gunzip -c "$TMP_GZ" | grep -vE '_|X|M|Y' > {output}
        rm -f "$TMP_GZ"
        """
