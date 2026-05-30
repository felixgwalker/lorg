import os
import csv


def write_per_element_tsv(scored_elements, species_list, output_dir):
    """Write a per-element conservation report as a TSV file.

    Columns: element_id, chrom, start, end, classification,
             combined_score, mean_sequence_identity, motif_retention,
             plus one identity_<species> column for every species.

    Returns the output path.
    """
    path = os.path.join(output_dir, "per_element_conservation.tsv")
    base_fields = [
        "element_id", "chrom", "start", "end",
        "classification", "combined_score",
        "mean_sequence_identity", "motif_retention",
    ]
    sp_fields = [f"identity_{sp}" for sp in species_list]
    fieldnames = base_fields + sp_fields

    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for elem in scored_elements:
            row = {}
            for f in base_fields:
                row[f] = elem.get(f, "")
            for sp in species_list:
                row[f"identity_{sp}"] = elem["per_species_identity"].get(sp, "")
            writer.writerow(row)
    return path


def write_conservation_scores(scored_elements, species_list, output_dir):
    path = os.path.join(output_dir, "conservation_scores.csv")
    base_fields = ["element_id", "mean_sequence_identity", "motif_retention",
                   "combined_score", "classification"]
    sp_fields = [f"identity_{sp}" for sp in species_list]
    fieldnames = base_fields + sp_fields

    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for elem in scored_elements:
            row = {f: elem[f] for f in base_fields}
            for sp in species_list:
                row[f"identity_{sp}"] = elem["per_species_identity"].get(sp, "")
            writer.writerow(row)
    return path


def write_motif_matrix(scored_elements, species_list, output_dir):
    path = os.path.join(output_dir, "motif_matrix.csv")
    if not scored_elements:
        return path

    all_motifs = list(scored_elements[0]["per_species_motifs"].get(species_list[0], {}).keys())
    fieldnames = ["element_id", "species"] + all_motifs

    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for elem in scored_elements:
            for sp in species_list:
                motifs = elem["per_species_motifs"].get(sp, {})
                row = {"element_id": elem["element_id"], "species": sp}
                row.update({m: int(motifs.get(m, False)) for m in all_motifs})
                writer.writerow(row)
    return path


def write_poorly_conserved(scored_elements, output_dir):
    path = os.path.join(output_dir, "poorly_conserved.csv")
    diverged = [e for e in scored_elements if e["classification"] == "diverged"]
    fieldnames = ["element_id", "combined_score", "mean_sequence_identity", "motif_retention"]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for elem in diverged:
            writer.writerow({f: elem[f] for f in fieldnames})
    return path
