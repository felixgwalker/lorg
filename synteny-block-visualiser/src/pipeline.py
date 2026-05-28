import os
import sys

try:
    from .fasta_reader import parse_fasta, make_demo_genomes
    from .kmer_indexer import build_kmer_index, find_unique_seeds
    from .chain_builder import chain_seeds
    from .rearrangement_detector import classify_blocks, summarize_rearrangements
    from .report import write_synteny_blocks, write_rearrangements
    from .dot_plotter import plot_dotplot
    from .ribbon_plotter import plot_ribbon
except ImportError:
    from fasta_reader import parse_fasta, make_demo_genomes
    from kmer_indexer import build_kmer_index, find_unique_seeds
    from chain_builder import chain_seeds
    from rearrangement_detector import classify_blocks, summarize_rearrangements
    from report import write_synteny_blocks, write_rearrangements
    from dot_plotter import plot_dotplot
    from ribbon_plotter import plot_ribbon


def run_pipeline(args):
    os.makedirs(args.output_dir, exist_ok=True)

    k = getattr(args, "kmer_size", 15)
    min_chain = getattr(args, "min_chain_score", 3)
    min_block_len = getattr(args, "min_block_length", 1000)

    if args.demo:
        genome1, genome2 = make_demo_genomes()
    else:
        if not args.fasta1 or not args.fasta2:
            print("ERROR: --fasta1 and --fasta2 are required unless --demo is used.", file=sys.stderr)
            return {"status": "error", "message": "missing inputs"}
        genome1 = parse_fasta(args.fasta1)
        genome2 = parse_fasta(args.fasta2)

    print(f"Genome 1: {len(genome1)} sequences, total {sum(len(v) for v in genome1.values())} bp")
    print(f"Genome 2: {len(genome2)} sequences, total {sum(len(v) for v in genome2.values())} bp")
    print(f"Building k-mer index (k={k})...")

    index1 = build_kmer_index(genome1, k=k)
    print(f"  {len(index1)} unique k-mers indexed from genome 1.")

    print("Finding unique seeds...")
    seeds = find_unique_seeds(genome2, index1, k=k)
    print(f"  {len(seeds)} unique seeds found.")

    print("Chaining seeds...")
    chains = chain_seeds(seeds, min_chain_score=min_chain)
    print(f"  {len(chains)} chains assembled.")

    print("Classifying blocks...")
    blocks = classify_blocks(chains, genome1, genome2, min_block_length=min_block_len, k=k)
    rearrangements = summarize_rearrangements(blocks)

    n_col = sum(1 for b in blocks if b["type"] == "collinear")
    n_inv = sum(1 for b in blocks if b["type"] == "inversion")
    n_tra = sum(1 for b in blocks if b["type"] == "translocation")
    print(f"  {len(blocks)} blocks: {n_col} collinear, {n_inv} inversions, {n_tra} translocations")

    blocks_path = write_synteny_blocks(blocks, args.output_dir)
    rear_path = write_rearrangements(rearrangements, args.output_dir)

    plot_paths = []
    if not getattr(args, "no_plot", False):
        print("Generating dot plot...")
        plot_paths += plot_dotplot(blocks, genome1, genome2, args.output_dir)
        print("Generating ribbon diagram...")
        plot_paths += plot_ribbon(blocks, genome1, genome2, args.output_dir)

    print(f"Outputs written to: {args.output_dir}")

    return {
        "status": "ok",
        "n_blocks": len(blocks),
        "collinear": n_col,
        "inversions": n_inv,
        "translocations": n_tra,
        "outputs": [blocks_path, rear_path] + plot_paths,
    }


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Synteny Block Visualiser")
    parser.add_argument("--fasta1", help="First genome FASTA")
    parser.add_argument("--fasta2", help="Second genome FASTA")
    parser.add_argument("--output-dir", default="output", dest="output_dir")
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--no-plot", action="store_true", dest="no_plot")
    parser.add_argument("--kmer-size", type=int, default=15, dest="kmer_size")
    parser.add_argument("--min-chain-score", type=int, default=3, dest="min_chain_score")
    parser.add_argument("--min-block-length", type=int, default=1000, dest="min_block_length")
    args = parser.parse_args()
    result = run_pipeline(args)
    return 0 if result.get("status") == "ok" else 1
