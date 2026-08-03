#!/bin/bash
HLA_TSV="pvacseq_hla_final.tsv"
VCF_DIR="/workdir_bctl/vep_snv_with_GT_and_GX"
OUT_TSV="vcf_sample_hla_mapping.tsv"
UNMATCHED="unmatched_vcf_samples.txt"

echo -e "vcf_sample\thla" > "$OUT_TSV"
> "$UNMATCHED"

for vcf in "$VCF_DIR"/*.filtered.vep.gt.expr.vcf; do
    vcf_sample=$(basename "$vcf" .filtered.vep.gt.expr.vcf)
    num=$(echo "$vcf_sample" | sed -E 's/^(BL|S)//; s/D$//')
    num_padded=$(printf "%03d" "$num" 2>/dev/null)

    hla=$(awk -F'\t' -v n="$num" -v np="$num_padded" '
        $1 == "WB"n || $1 == "Wb"n || $1 == "Plasmanormal"n || $1 == "Plasmanormal"np {print $2; exit}
    ' "$HLA_TSV")

    if [ -z "$hla" ]; then
        echo "$vcf_sample (numeric ID: $num)" >> "$UNMATCHED"
    else
        echo -e "${vcf_sample}\t${hla}" >> "$OUT_TSV"
    fi
done

echo "Done: $OUT_TSV"
echo "Unmatched samples (need manual review):"
cat "$UNMATCHED"
