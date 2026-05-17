IMPACT_SCORES = {
    "exonic_nonsynonymous": 10,
    "splice_site": 8,
    "exonic_synonymous": 3,
    "intronic": 2,
    "intergenic": 1,
}


def _overlaps_feature(pos, features, margin=0):
    for feat in features:
        if (feat["start"] - margin) <= pos <= (feat["end"] + margin):
            return feat
    return None


def _parse_gff3(gff_path):
    features = []
    with open(gff_path, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            if parts[2] in ("exon", "CDS"):
                features.append({
                    "chrom": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                    "feature": parts[2],
                    "strand": parts[6],
                })
    return features


def annotate_impact(variants, gff_path=None):
    features_by_chrom = {}
    if gff_path:
        try:
            all_features = _parse_gff3(gff_path)
            for feat in all_features:
                chrom = feat["chrom"]
                if chrom not in features_by_chrom:
                    features_by_chrom[chrom] = []
                features_by_chrom[chrom].append(feat)
        except Exception:
            features_by_chrom = {}

    annotated = []
    for v in variants:
        pos = v.get("pos", 0)
        chrom = v.get("chrom", "")
        chrom_features = features_by_chrom.get(chrom, [])

        if chrom_features:
            splice_feat = _overlaps_feature(pos, chrom_features, margin=5)
            exon_feat = _overlaps_feature(pos, chrom_features, margin=0)
            if splice_feat and not exon_feat:
                impact_category = "splice_site"
            elif exon_feat:
                if v.get("type") == "SNV":
                    impact_category = "exonic_nonsynonymous"
                else:
                    impact_category = "exonic_nonsynonymous"
            else:
                impact_category = "intergenic"
        else:
            if v.get("type") == "SNV":
                impact_category = "intergenic"
            elif v.get("variant_class", "").startswith("SV"):
                impact_category = "exonic_nonsynonymous"
            elif v.get("variant_class", "").startswith("LARGE"):
                impact_category = "intronic"
            else:
                impact_category = "intergenic"

        v["impact_category"] = impact_category
        v["impact_score"] = IMPACT_SCORES.get(impact_category, 1)
        annotated.append(v)

    return annotated
