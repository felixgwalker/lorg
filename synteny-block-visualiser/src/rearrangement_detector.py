def _compute_identity(seeds, k=15):
    if not seeds:
        return 0.0
    matches = 0
    total = 0
    for pos1, pos2, kmer in seeds:
        matches += k
        total += k
    return 1.0 if total == 0 else matches / total


def classify_blocks(chains, genome1, genome2, min_block_length=1000, k=15):
    blocks = []
    for chain in chains:
        chrom1 = chain["chrom1"]
        chrom2 = chain["chrom2"]
        g1_start = chain["g1_start"]
        g1_end = chain["g1_end"] + k
        g2_start = chain["g2_start"]
        g2_end = chain["g2_end"] + k

        g1_span = g1_end - g1_start
        g2_span = g2_end - g2_start
        block_len = max(g1_span, g2_span)

        if block_len < min_block_length:
            continue

        g2_increasing = (g2_end - g2_start) > 0

        if chrom1 == chrom2:
            block_type = "collinear" if g2_increasing else "inversion"
        else:
            g1_prefix = chrom1.split("_")[0] if "_" in chrom1 else chrom1
            g2_prefix = chrom2.split("_")[0] if "_" in chrom2 else chrom2
            if g1_prefix == g2_prefix:
                block_type = "collinear" if g2_increasing else "inversion"
            else:
                block_type = "translocation"

        identity = _compute_identity(chain["seeds"], k)

        blocks.append({
            "g1_chrom": chrom1,
            "g1_start": g1_start,
            "g1_end": g1_end,
            "g2_chrom": chrom2,
            "g2_start": g2_start,
            "g2_end": g2_end,
            "type": block_type,
            "n_seeds": chain["n_seeds"],
            "identity": round(identity, 4),
            "seeds": chain["seeds"],
        })
    return blocks


def summarize_rearrangements(blocks):
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
