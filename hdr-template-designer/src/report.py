"""Write HDR template outputs: FASTA, QC CSV, arm variants CSV."""

from pathlib import Path
import pandas as pd


def write_template_fasta(template: dict, output_dir: Path, name: str = "hdr_template") -> Path:
    out_path = output_dir / "hdr_template.fa"
    full_seq = template["full_template"]
    with open(out_path, "w") as fh:
        fh.write(f">{name} left_arm={template['left_arm_len']}bp edit={template['edit_seq']} "
                 f"right_arm={template['right_arm_len']}bp\n")
        for i in range(0, len(full_seq), 60):
            fh.write(full_seq[i:i + 60] + "\n")
    return out_path


def write_qc_report(qc_checks: list[dict], output_dir: Path) -> Path:
    out_path = output_dir / "qc_report.csv"
    df = pd.DataFrame(qc_checks)
    df.to_csv(out_path, index=False)
    return out_path


def write_arm_variants(variants: list[dict], output_dir: Path) -> Path:
    out_path = output_dir / "arm_variants.csv"
    df = pd.DataFrame(variants)
    df.to_csv(out_path, index=False)
    return out_path
