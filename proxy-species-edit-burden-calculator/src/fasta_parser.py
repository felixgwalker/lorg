import numpy as np


def parse_fasta(path):
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
