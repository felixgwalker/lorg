from .sequence_analyzer import extract_window, compute_repeat_score, compute_complexity_score
from .feature_annotator import gene_density_score, regulatory_proximity_score


def score_window(sequences, genes, chrom, win_start, win_end):
    seq = extract_window(sequences, chrom, win_start, win_end)

    gd_score = gene_density_score(genes, chrom, win_start, win_end)
    reg_score = regulatory_proximity_score(sequences, chrom, win_start, win_end)
    rep_score = compute_repeat_score(seq)
    cplx_score = compute_complexity_score(seq)
    total = gd_score + reg_score + rep_score + cplx_score

    return {
        "chrom": chrom,
        "start": win_start,
        "end": win_end,
        "gene_density_score": gd_score,
        "regulatory_score": reg_score,
        "repeat_score": rep_score,
        "complexity_score": cplx_score,
        "total_score": total,
    }


def scan_locus(sequences, genes, locus, window_size=1000, step=200, flank=50000):
    chrom = locus["chrom"]
    center = (locus["start"] + locus["end"]) // 2
    scan_start = max(0, center - flank)
    chrom_len = len(sequences.get(chrom, ""))
    scan_end = min(chrom_len if chrom_len else center + flank, center + flank)

    results = []
    pos = scan_start
    while pos + window_size <= scan_end:
        row = score_window(sequences, genes, chrom, pos, pos + window_size)
        row["locus_name"] = locus["name"]
        results.append(row)
        pos += step

    return results
