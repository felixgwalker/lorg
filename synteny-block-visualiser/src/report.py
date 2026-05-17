import os
import csv


def write_synteny_blocks(blocks, output_dir):
    path = os.path.join(output_dir, "synteny_blocks.csv")
    fieldnames = [
        "g1_chrom", "g1_start", "g1_end",
        "g2_chrom", "g2_start", "g2_end",
        "type", "n_seeds", "identity",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for b in blocks:
            writer.writerow({f: b[f] for f in fieldnames})
    return path


def write_rearrangements(rearrangements, output_dir):
    path = os.path.join(output_dir, "rearrangements.csv")
    fieldnames = [
        "type", "g1_chrom", "g1_start", "g1_end",
        "g2_chrom", "g2_start", "g2_end",
        "n_seeds", "identity",
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for r in rearrangements:
            writer.writerow({f: r[f] for f in fieldnames})
    return path
