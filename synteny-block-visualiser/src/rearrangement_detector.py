"""Classify chained synteny blocks and summarize rearrangements."""


def _chrom_genome_prefix(chrom):
    """Return the genome prefix (e.g. 'g1' from 'g1_chr1') for chromosome matching."""
    parts = chrom.split("_")
    if len(parts) >= 2:
        return parts[0]
    return chrom


def _chrom_base(chrom):
    """Return the base chromosome name (e.g. 'chr1' from 'g1_chr1')."""
    parts = chrom.split("_", 1)
    if len(parts) == 2:
        return parts[1]
    return chrom


def classify_blocks(chains, genome1, genome2, min_block_length=1000, k=15):
    """Classify each chain as collinear, inversion, or translocation.

    Rules:
      - If chrom1 base name != chrom2 base name (inter-chromosomal): translocation
      - Else if strand == "-": inversion
      - Else: collinear

    Only chains whose genomic span >= min_block_length are kept.
    """
    blocks = []
    for chain in chains:
        chrom1 = chain["chrom1"]
        chrom2 = chain["chrom2"]
        strand = chain.get("strand", "+")
        g1_start = chain["g1_start"]
        g1_end = chain["g1_end"] + k
        g2_start = chain["g2_start"]
        g2_end = chain["g2_end"] + k

        g1_span = abs(g1_end - g1_start)
        g2_span = abs(g2_end - g2_start)
        block_len = max(g1_span, g2_span)

        if block_len < min_block_length:
            continue

        # Determine block type
        base1 = _chrom_base(chrom1)
        base2 = _chrom_base(chrom2)
        if base1 != base2:
            block_type = "translocation"
        elif strand == "-":
            block_type = "inversion"
        else:
            block_type = "collinear"

        # Ensure g2 coordinates are in ascending order for collinear/translocation display
        if g2_start > g2_end:
            g2_start, g2_end = g2_end, g2_start

        identity = _compute_identity(chain["seeds"], k)

        blocks.append({
            "g1_chrom": chrom1,
            "g1_start": g1_start,
            "g1_end": g1_end,
            "g2_chrom": chrom2,
            "g2_start": g2_start,
            "g2_end": g2_end,
            "type": block_type,
            "strand": strand,
            "n_seeds": chain["n_seeds"],
            "identity": round(identity, 4),
            "seeds": chain["seeds"],
        })
    return blocks


def _compute_identity(seeds, k=15):
    """Estimate sequence identity from seed count (all seed k-mers are exact matches)."""
    if not seeds:
        return 0.0
    # All k-mers in seeds are exact matches; identity is approximated as 1.0
    # since seeds only contain identical k-mers by construction.
    return 1.0


def summarize_rearrangements(blocks):
    """Return a list of dicts describing inversions and translocations."""
    inversions = [b for b in blocks if b["type"] == "inversion"]
    translocations = [b for b in blocks if b["type"] == "translocation"]
    rows = []
    for b in inversions + translocations:
        rows.append({
            "type": b["type"],
            "g1_chrom": b["g1_chrom"],
            "g1_start": b["g1_start"],
            "g1_end": b["g1_end"],
            "g2_chrom": b["g2_chrom"],
            "g2_start": b["g2_start"],
            "g2_end": b["g2_end"],
            "n_seeds": b["n_seeds"],
            "identity": b["identity"],
        })
    return rows
