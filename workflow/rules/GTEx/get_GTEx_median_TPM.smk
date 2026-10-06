# GTEx release -> download location. Links live here, not in the config (the config only
# holds local paths and parameters: GTEx_releases selects the release).
GTEx_BASE_URL = "https://storage.googleapis.com/adult-gtex/bulk-gex"
GTEx_median_TPM_files = {
    "v8": "v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
    "v11": "v11/rna-seq/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz",
}

def GTEx_median_TPM_url(release):
    if release not in GTEx_median_TPM_files:
        raise ValueError(
            f"ERROR, GTEx release '{release}' has no download location "
            f"(known: {', '.join(sorted(GTEx_median_TPM_files))}; add it to GTEx_median_TPM_files "
            f"in rules/GTEx/get_GTEx_median_TPM.smk)"
        )
    return f"{GTEx_BASE_URL}/{GTEx_median_TPM_files[release]}"


rule get_GTEx_median_TPM:
    # Downloads the GTEx tissue-level median gene TPM matrix (GCT v1.2) for one GTEx release.
    # Gene-level quantification only: no genomic coordinates are used downstream, so the file
    # is the same whatever assembly_code the backbone is built on.
    output:
        "resources/GTEx/{GTEx_release}/gene_median_tpm.gct.gz"
    params:
        url=lambda wildcards: GTEx_median_TPM_url(wildcards.GTEx_release)
    wildcard_constraints:
        GTEx_release="v[0-9]+"
    log:
        "logs/get_GTEx_median_TPM/{GTEx_release}.log"
    shell:
        # Download to a .part file and move it into place only after two integrity checks:
        # (1) valid gzip stream (catches truncated downloads and HTML error pages);
        # (2) first line of the decompressed file is the GCT version tag "#1.2" (catches the
        #     wrong file behind a valid URL). Every failure branch removes the .part file
        #     explicitly, because Snakemake's strict mode (set -euo pipefail) would otherwise
        #     exit before any cleanup. "|| true" keeps head's SIGPIPE out of pipefail.
        # No --continue: the file is ~10 MB and resuming would risk appending to a stale
        # .part from a different URL.
        """
        echo "URL: {params.url}" > {log}
        if ! wget --tries=5 --timeout=60 -O {output}.part {params.url} 2>> {log}; then
            echo "ERROR: download failed" >> {log}
            rm -f {output}.part
            exit 1
        fi
        if ! gzip -t {output}.part 2>> {log}; then
            echo "ERROR: {output}.part is not a valid gzip file" >> {log}
            rm -f {output}.part
            exit 1
        fi
        header=$(gzip -cd {output}.part 2>/dev/null | head -n 1 | tr -d '\\r' || true)
        if [ "$header" != "#1.2" ]; then
            echo "ERROR: {output}.part is not a GCT v1.2 file (first line: '$header')" >> {log}
            rm -f {output}.part
            exit 1
        fi
        mv {output}.part {output}
        """