import numpy as np


def parse_fasta(path):
    """Parse a FASTA file and return {seq_name: sequence_str}."""
    sequences = {}
    current_id = None
    current_seq = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
    if current_id is not None:
        sequences[current_id] = "".join(current_seq)
    return sequences


def make_demo_sequences(length=10000, divergence=0.05, seed=42):
    """Return two synthetic sequences (proxy, target) of *length* bp with ~divergence SNP rate.

    The proxy sequence is random nucleotide sequence; the target is derived from
    it by introducing SNPs and small indels to reach approximately *divergence*
    overall sequence divergence (~5% by default).

    Returns
    -------
    tuple[str, str]
        (proxy_seq, target_seq) — plain strings, no FASTA headers.
    """
    rng = np.random.default_rng(seed)
    bases = list("ACGT")

    proxy_arr = rng.choice(bases, size=length)
    target_arr = list(proxy_arr)

    # Introduce SNPs to reach ~4% divergence.
    n_snps = int(length * 0.04)
    snp_positions = rng.choice(length, size=n_snps, replace=False)
    for pos in snp_positions:
        orig = target_arr[pos]
        alts = [b for b in bases if b != orig]
        target_arr[pos] = alts[rng.integers(0, len(alts))]

    # Introduce small indels (~1% additional divergence).
    n_indels = int(length * 0.005)
    indel_positions = sorted(
        rng.choice(max(length - 20, 1), size=n_indels, replace=False).tolist()
    )
    offset = 0
    for pos in indel_positions:
        adj = pos + offset
        adj = min(adj, max(len(target_arr) - 1, 0))
        if rng.random() < 0.5:
            indel_len = int(rng.integers(1, 10))
            insert = list(rng.choice(bases, size=indel_len))
            target_arr = target_arr[:adj] + insert + target_arr[adj:]
            offset += indel_len
        else:
            del_len = int(rng.integers(1, 8))
            del_len = min(del_len, len(target_arr) - adj)
            if del_len > 0:
                target_arr = target_arr[:adj] + target_arr[adj + del_len:]
                offset -= del_len

    return "".join(proxy_arr), "".join(target_arr)


def generate_demo_genomes(n_chroms=3, chrom_len=10000, n_snvs=50, n_indels=10, rng=None):
    if rng is None:
        rng = np.random.default_rng(42)

    bases = list("ACGT")
    proxy = {}
    target = {}

    for chrom_idx in range(n_chroms):
        chrom_name = f"chr{chrom_idx+1}"
        proxy_seq = list(rng.choice(bases, size=chrom_len))
        target_seq = list(proxy_seq)

        snv_positions = rng.choice(chrom_len, size=n_snvs, replace=False)
        for pos in snv_positions:
            orig = target_seq[pos]
            alts = [b for b in bases if b != orig]
            target_seq[pos] = alts[rng.integers(0, len(alts))]

        indel_positions = sorted(rng.choice(chrom_len - 20, size=n_indels, replace=False))
        offset = 0
        for pos in indel_positions:
            adj_pos = pos + offset
            adj_pos = min(adj_pos, len(target_seq) - 1)
            if rng.random() < 0.5:
                indel_len = int(rng.integers(1, 15))
                insert = list(rng.choice(bases, size=indel_len))
                target_seq = target_seq[:adj_pos] + insert + target_seq[adj_pos:]
                offset += indel_len
            else:
                del_len = int(rng.integers(1, 10))
                del_len = min(del_len, len(target_seq) - adj_pos)
                target_seq = target_seq[:adj_pos] + target_seq[adj_pos + del_len:]
                offset -= del_len

        proxy[chrom_name] = "".join(proxy_seq)
        target[chrom_name] = "".join(target_seq)

    return proxy, target
