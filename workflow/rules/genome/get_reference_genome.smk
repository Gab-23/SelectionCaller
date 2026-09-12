rule get_reference_genome:
    output:
        "resources/reference_data/genome/{assembly_code}.fa"
    log:
        "logs/get_reference_genome/{assembly_code}.log"
    shell:
        # NOTE (fix): the original piped `wget -qO- ... | gunzip` directly - confirmed on a
        # real cluster run this failed with "gzip: stdin: unexpected end of file" (a truncated
        # download, likely a network hiccup partway through this ~900MB-1GB file), but the
        # ACTUAL cause was invisible: -q suppresses every wget diagnostic, so only gunzip's
        # downstream symptom ever reached the log. Fixed by downloading to a real file first
        # (with --tries/--timeout for automatic retry on transient failures, and --continue so
        # a retry resumes instead of restarting from scratch), decompressing as a separate
        # step, and dropping -q so wget's own errors are now actually visible in the log.
        # NOTE (fix, part 2): also using hgdownload.cse.ucsc.edu, an old/inconsistent UCSC
        # hostname - get_chrom_sizes.smk and get_cytobands.smk both already use
        # hgdownload.soe.ucsc.edu (per your own prior fix to get_chrom_sizes for this exact
        # reason). Switched to match - very plausibly the actual root cause of the truncated
        # download above, not just an unlucky transient network blip.
        """
        wget --tries=5 --timeout=60 --continue -O {output}.gz \
        https://hgdownload.soe.ucsc.edu/goldenPath/{wildcards.assembly_code}/bigZips/latest/{wildcards.assembly_code}.fa.gz \
        2>> {log}
        gunzip -f {output}.gz 2>> {log}
        """
