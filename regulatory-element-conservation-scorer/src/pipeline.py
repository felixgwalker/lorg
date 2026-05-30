import os
import sys

from .bed_parser import parse_bed, make_demo_bed
from .sequence_extractor import load_fasta, extract_element_sequence, extract_sequences, demo_species_sequences
from .conservation_scorer import score_element_conservation, score_conservation
from .report import write_conservation_scores, write_motif_matrix, write_poorly_conserved, write_per_element_tsv
from .plot import plot_conservation_heatmap


def run_pipeline(args):
    os.makedirs(args.output_dir, exist_ok=True)

    if args.demo:
        elements = make_demo_bed()
        species_list = ["human", "mouse", "zebrafish", "fruitfly"]
        ref_species = "human"
        species_seqs = demo_species_sequences(elements, species_list)
    else:
        if not args.bed or not args.fastas:
            print("ERROR: --bed and --fastas are required unless --demo is used.", file=sys.stderr)
            return {"status": "error", "message": "missing inputs"}

        elements = parse_bed(args.bed)
        fasta_paths = args.fastas
        species_list = [os.path.splitext(os.path.basename(p))[0] for p in fasta_paths]
        ref_species = species_list[0]

        genome_seqs = [load_fasta(p) for p in fasta_paths]
        species_seqs = {}
        for sp, seqs in zip(species_list, genome_seqs):
            species_seqs[sp] = {}
            for elem in elements:
                seq = extract_element_sequence(seqs, elem)
                species_seqs[sp][elem["id"]] = seq

    scored_elements = []
    for elem in elements:
        result = score_element_conservation(elem["id"], species_seqs, species_list, ref_species)
        result["chrom"] = elem["chrom"]
        result["start"] = elem["start"]
        result["end"] = elem["end"]

        # Attach ConservationResult metrics computed across all species sequences
        all_seqs_for_elem = [
            species_seqs.get(sp, {}).get(elem["id"], "")
            for sp in species_list
        ]
        cr = score_conservation(all_seqs_for_elem)
        result["entropy_conservation_score"] = cr.entropy_conservation_score
        result["conserved_positions_fraction"] = cr.conserved_positions_fraction
        result["mean_pairwise_identity_pct"] = cr.mean_pairwise_identity

        scored_elements.append(result)

    scores_path = write_conservation_scores(scored_elements, species_list, args.output_dir)
    motif_path = write_motif_matrix(scored_elements, species_list, args.output_dir)
    poor_path = write_poorly_conserved(scored_elements, args.output_dir)
    tsv_path = write_per_element_tsv(scored_elements, species_list, args.output_dir)

    plot_paths = []
    if not getattr(args, "no_plot", False):
        plot_paths = plot_conservation_heatmap(scored_elements, species_list, args.output_dir)

    conserved = sum(1 for e in scored_elements if e["classification"] == "conserved")
    partial = sum(1 for e in scored_elements if e["classification"] == "partially_conserved")
    diverged = sum(1 for e in scored_elements if e["classification"] == "diverged")

    print(f"Scored {len(scored_elements)} regulatory elements across {len(species_list)} species.")
    print(f"  Conserved: {conserved}  Partially conserved: {partial}  Diverged: {diverged}")
    print(f"Outputs written to: {args.output_dir}")

    return {
        "status": "ok",
        "n_elements": len(scored_elements),
        "n_species": len(species_list),
        "conserved": conserved,
        "partially_conserved": partial,
        "diverged": diverged,
        "outputs": [scores_path, motif_path, poor_path, tsv_path] + plot_paths,
    }


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Regulatory Element Conservation Scorer")
    parser.add_argument("--bed", help="BED file of regulatory elements")
    parser.add_argument("--fastas", nargs="+", help="FASTA files, one per species (first = reference)")
    parser.add_argument("--output-dir", default="output", dest="output_dir")
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--no-plot", action="store_true", dest="no_plot")
    args = parser.parse_args()
    result = run_pipeline(args)
    return 0 if result.get("status") == "ok" else 1
