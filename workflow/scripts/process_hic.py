import sys
import argparse
import warnings
import numpy
import pandas
warnings.filterwarnings("ignore")


def orient_eigenvectors(df_eigen, phasing_track):
    """Orientamento degli autovettori (fix 20/10/2026).

    Il segno di un autovettore restituito da numpy.linalg.eigh e' arbitrario: senza un
    orientamento, lo stesso compartimento puo' risultare positivo su un cromosoma e negativo
    su un altro (osservato: segno invertito in circa meta' delle coppie autovettore x
    cromosoma). Per ogni autovettore si sceglie il segno che rende POSITIVA la correlazione di
    Spearman con una traccia di riferimento (qui il contenuto in GC del bin), come fa il
    "phasing" di cooltools. Per eigen_1 la convenzione fa coincidere i valori positivi con il
    compartimento ricco in GC (A); per gli autovettori successivi e' solo una convenzione,
    ma e' la stessa su tutti i cromosomi.

    df_eigen: DataFrame con le colonne eigen_*; phasing_track: Series allineata (stesso indice).
    Restituisce (df orientato, dict colonna -> {"rho": rho prima dell'orientamento, "flipped": bool}).
    """
    out = df_eigen.copy()
    report = {}
    track_rank = phasing_track.rank()
    for col in [c for c in out.columns if c.startswith("eigen_")]:
        ok = out[col].notna() & phasing_track.notna()
        if ok.sum() < 3:
            report[col] = {"rho": numpy.nan, "flipped": False}
            continue
        rho = numpy.corrcoef(out.loc[ok, col].rank(), track_rank[ok])[0, 1]   # Spearman = Pearson sui ranghi
        flipped = bool(numpy.isfinite(rho) and rho < 0)
        if flipped:
            out[col] = -out[col]
        report[col] = {"rho": float(rho) if numpy.isfinite(rho) else numpy.nan, "flipped": flipped}
    return out, report


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('--input_path_matrix', dest='input_path_matrix', required=True, type=str)
    parser.add_argument('--input_path_bins',   dest='input_path_bins',   required=True, type=str)
    parser.add_argument('--input_path_gc',     dest='input_path_gc',     required=False, type=str, default=None,
                        help='tabella binID / GC_content (data/GC/...) usata per orientare gli autovettori')
    parser.add_argument('--phasing',           dest='phasing',           required=False, type=str, default='gc',
                        choices=['gc', 'none'], help="gc: orienta ogni autovettore con corr(GC) > 0; none: segno grezzo di eigh")
    parser.add_argument('--output_path_eigen', dest='output_path_eigen', required=True, type=str)
    #parser.add_argument('--output_path_insulation', dest='output_path_insulation', required=True, type=str, help='')

    args = parser.parse_args()

    CONTACT_MATRIX = args.input_path_matrix
    BINS_FILE      = args.input_path_bins
    OUTFILE_EIGEN  = args.output_path_eigen
    #OUTFILE_INSULATION = args.output_path_insulation

    if args.phasing == 'gc' and args.input_path_gc is None:
        sys.exit("ERROR: --phasing gc richiede --input_path_gc")

    # Load bins (chr, start, end, binID) — cooler-indexed
    bins = pandas.read_csv(BINS_FILE, sep="\t", header=None, names=["chr", "start", "end", "binID"])

    contact = numpy.loadtxt(CONTACT_MATRIX, delimiter="\t") # load contact matrix
    print("NaNs/infs before cleaning:", numpy.isnan(contact).sum(), numpy.isinf(contact).sum())
    contact = numpy.nan_to_num(contact, nan=0.0, posinf=0.0, neginf=0.0) # set NaNs to 0
    contact = (contact + contact.T) / 2 # force matrix to be symmetric
    print("NaNs/infs after cleaning:", numpy.isnan(contact).sum(), numpy.isinf(contact).sum())

    # Sanity check
    assert contact.shape[0] == len(bins), \
        f"Matrix rows ({contact.shape[0]}) != bins ({len(bins)})"

    row_sums = contact.sum(axis=1) # compute rowSums
    bad_bins = row_sums == 0 # get all-0 columns
    valid = ~bad_bins # get valid bins
    print(f"Bad bins: {bad_bins.sum()} / {len(bad_bins)}")
    contact_valid = contact[numpy.ix_(valid, valid)] # subset matrix
    eigenvalues, eigenvectors = numpy.linalg.eigh(contact_valid)  # ascending, so reverse
    n_valid = contact_valid.shape[0]
    n_eigen = min(10, n_valid)  # fix: with fewer than 10 valid bins (possible for a small toy
                                # chromosome, or if an unusually large fraction of bins are
                                # dropped as "bad"), eigh() can only return as many
                                # eigenvectors as there are valid bins - requesting a fixed 10
                                # unconditionally caused a shape-mismatch crash when fewer were
                                # available. Real genome-scale Hi-C data has far more than 10
                                # valid bins per chromosome, so this only ever engages as a
                                # defensive fallback, not a change to normal behavior.
    top10 = eigenvectors[:, ::-1][:, :n_eigen] # get top eigenvectors (up to 10)
    n_bins = contact.shape[0] # get number of bins
    full_eigvec = numpy.full((n_bins, n_eigen), numpy.nan) # initialize with NaNs
    full_eigvec[valid, :] = top10 # fill in eigenvectors

    # Assemble output: binID + eigenvectors
    eigen_cols = [f"eigen_{i+1}" for i in range(n_eigen)]
    df_eigen = pandas.DataFrame(full_eigvec, columns=eigen_cols)

    # Orientamento (fix 20/10/2026): vedi orient_eigenvectors
    if args.phasing == 'gc':
        gc = pandas.read_csv(args.input_path_gc, sep="\t", header=0)
        if not {"binID", "GC_content"}.issubset(gc.columns):
            sys.exit(f"ERROR: {args.input_path_gc} deve avere le colonne binID e GC_content")
        gc_track = bins[["binID"]].merge(gc[["binID", "GC_content"]], on="binID", how="left")["GC_content"]
        print(f"Phasing con GC: {gc_track.notna().sum()} / {len(bins)} bin con GC disponibile", file=sys.stderr)
        df_eigen, report = orient_eigenvectors(df_eigen, gc_track.reset_index(drop=True))
        for col, r in report.items():
            print(f"{col}: rho(GC) prima = {r['rho']:.3f}, segno invertito = {r['flipped']}", file=sys.stderr)

    df_out   = pandas.concat([bins[["binID"]].reset_index(drop=True), df_eigen], axis=1)
    df_out.to_csv(OUTFILE_EIGEN, sep="\t", index=False, na_rep="NA")


if __name__ == "__main__":
    main()