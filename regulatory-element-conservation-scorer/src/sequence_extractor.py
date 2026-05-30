import random


def load_fasta(fasta_path):
    sequences = {}
    current = None
    parts = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if current is not None:
                    sequences[current] = "".join(parts)
                current = line[1:].split()[0]
                parts = []
            else:
                parts.append(line.upper())
    if current is not None:
        sequences[current] = "".join(parts)
    return sequences


def extract_element_sequence(sequences, element, offset=0):
    chrom = element["chrom"]
    seq = sequences.get(chrom, "")
    if not seq:
        return ""
    start = max(0, element["start"] + offset)
    end = min(len(seq), element["end"] + offset)
    return seq[start:end]


def extract_sequences(regions, genome_fasta):
    """Extract sequences for a list of BED regions from a genome FASTA.

    Parameters
    ----------
    regions : list of RegionRecord dicts (chrom, start, end, name, …)
    genome_fasta : str — path to a FASTA file, or None for demo/synthetic mode

    Returns
    -------
    list of (name, seq) tuples, one per region.
    Coordinates are clipped to chromosome end gracefully.
    When genome_fasta is None (or the file cannot be opened), synthetic
    sequences are generated from the region length so the rest of the
    pipeline can run without a real genome.
    """
    if genome_fasta is None:
        return _synthetic_sequences(regions)

    try:
        genome = load_fasta(genome_fasta)
    except (OSError, IOError):
        return _synthetic_sequences(regions)

    result = []
    for reg in regions:
        name = reg.get("name") or reg.get("id") or f"{reg['chrom']}:{reg['start']}-{reg['end']}"
        chrom = reg["chrom"]
        chrom_seq = genome.get(chrom, "")
        if not chrom_seq:
            result.append((name, ""))
            continue
        start = max(0, reg["start"])
        end = min(len(chrom_seq), reg["end"])
        result.append((name, chrom_seq[start:end]))
    return result


def _synthetic_sequences(regions):
    """Generate deterministic synthetic sequences when no genome is available."""
    random.seed(42)
    bases = list("ACGT")
    result = []
    for reg in regions:
        name = reg.get("name") or reg.get("id") or f"{reg['chrom']}:{reg['start']}-{reg['end']}"
        length = max(1, reg["end"] - reg["start"])
        seq = "".join(random.choices(bases, k=length))
        result.append((name, seq))
    return result


def demo_species_sequences(elements, species_list):
    random.seed(123)
    bases = list("ACGT")

    def rand_seq(n):
        return "".join(random.choices(bases, k=n))

    conservation_levels = {
        "enhancer_01": [1.0, 0.95, 0.85, 0.30],
        "promoter_01": [1.0, 0.98, 0.96, 0.92],
        "enhancer_02": [1.0, 0.90, 0.75, 0.60],
        "silencer_01": [1.0, 0.80, 0.55, 0.25],
        "enhancer_03": [1.0, 0.92, 0.88, 0.35],
        "promoter_02": [1.0, 0.97, 0.93, 0.88],
        "insulator_01": [1.0, 0.85, 0.70, 0.50],
        "enhancer_04": [1.0, 0.78, 0.45, 0.20],
        "promoter_03": [1.0, 0.99, 0.95, 0.91],
        "enhancer_05": [1.0, 0.88, 0.65, 0.40],
        "silencer_02": [1.0, 0.72, 0.48, 0.22],
        "enhancer_06": [1.0, 0.94, 0.80, 0.55],
    }

    species_seqs = {sp: {} for sp in species_list}

    for elem in elements:
        eid = elem["id"]
        length = elem["end"] - elem["start"]
        ref_seq = list(rand_seq(length))
        levels = conservation_levels.get(eid, [1.0, 0.8, 0.6, 0.4])

        for sp_idx, sp in enumerate(species_list):
            cons = levels[sp_idx] if sp_idx < len(levels) else 0.5
            mutated = list(ref_seq)
            for i in range(length):
                if random.random() > cons:
                    mutated[i] = random.choice([b for b in bases if b != ref_seq[i]])
            species_seqs[sp][eid] = "".join(mutated)

    return species_seqs
