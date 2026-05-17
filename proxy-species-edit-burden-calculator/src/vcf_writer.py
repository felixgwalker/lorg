import os
from datetime import date


def write_vcf(variants, output_dir, proxy_name="PROXY", target_name="TARGET"):
    path = os.path.join(output_dir, "variants.vcf")
    today = date.today().strftime("%Y%m%d")

    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##fileDate={today}\n")
        fh.write(f"##source=proxy-species-edit-burden-calculator\n")
        fh.write(f"##reference={proxy_name}\n")
        fh.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant class">\n')
        fh.write('##INFO=<ID=IMPACT,Number=1,Type=String,Description="Impact category">\n')
        fh.write('##INFO=<ID=WEIGHT,Number=1,Type=Integer,Description="Edit burden weight">\n')
        fh.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for i, v in enumerate(variants):
            chrom = v.get("chrom", "chrUn")
            pos = v.get("pos", 0) + 1
            ref = v.get("ref", "N") or "N"
            alt = v.get("alt", ".") or "."
            vc = v.get("variant_class", "SNV")
            impact = v.get("impact_category", "intergenic")
            weight = v.get("weight", 1)
            info = f"SVTYPE={vc};IMPACT={impact};WEIGHT={weight}"
            fh.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\n")

    return path
