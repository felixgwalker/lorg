import os


def parse_bed(bed_path):
    loci = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                parts = line.split()
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3] if len(parts) > 3 else f"{chrom}:{start}-{end}"
            loci.append({"chrom": chrom, "start": start, "end": end, "name": name})
    return loci


def demo_loci():
    return [
        {"chrom": "chr1", "start": 50000, "end": 51000, "name": "locus_A"},
        {"chrom": "chr1", "start": 100000, "end": 101000, "name": "locus_B"},
        {"chrom": "chr2", "start": 75000, "end": 76000, "name": "locus_C"},
    ]


def read_locus_context(fasta_path, chrom, start, end, flank_bp=5000):
    """Read target sequence plus flanking context from a FASTA file.

    Returns a dict with keys:
        chrom, target_start, target_end,
        context_start, context_end,
        sequence (full context including flanks),
        target_sequence (just the requested interval).
    """
    sequences = {}
    current_chrom = None
    parts = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if current_chrom is not None:
                    sequences[current_chrom] = "".join(parts)
                current_chrom = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.upper())
    if current_chrom is not None:
        sequences[current_chrom] = "".join(parts)

    chrom_seq = sequences.get(chrom, "")
    chrom_len = len(chrom_seq)

    context_start = max(0, start - flank_bp)
    context_end = min(chrom_len, end + flank_bp)

    full_context = chrom_seq[context_start:context_end]
    target_seq = chrom_seq[max(0, start):min(chrom_len, end)]

    return {
        "chrom": chrom,
        "target_start": start,
        "target_end": end,
        "context_start": context_start,
        "context_end": context_end,
        "sequence": full_context,
        "target_sequence": target_seq,
    }


def make_demo_locus(chrom="chr1", start=50000, end=51000, flank_bp=5000):
    """Return a locus context dict built from synthetic demo data.

    Uses the same synthetic genome as make_demo_sequences() so that
    downstream functions receive realistic (non-empty) sequence.
    """
    from .sequence_analyzer import make_demo_sequences  # local import to avoid circular
    sequences = make_demo_sequences()

    chrom_seq = sequences.get(chrom, "")
    chrom_len = len(chrom_seq)

    context_start = max(0, start - flank_bp)
    context_end = min(chrom_len, end + flank_bp)

    full_context = chrom_seq[context_start:context_end]
    target_seq = chrom_seq[max(0, start):min(chrom_len, end)]

    return {
        "chrom": chrom,
        "target_start": start,
        "target_end": end,
        "context_start": context_start,
        "context_end": context_end,
        "sequence": full_context,
        "target_sequence": target_seq,
    }
