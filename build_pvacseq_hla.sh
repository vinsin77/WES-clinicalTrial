#!/bin/bash
# Build final pVACseq HLA allele string per sample from hlahd_alleles.tsv

IN_TSV="hlahd_alleles.tsv"
OUT_TSV="pvacseq_hla_final.tsv"

echo -e "sample\thla" > "$OUT_TSV"

tail -n +2 "$IN_TSV" | while IFS=$'\t' read -r sample a b c drb1 dq dp; do
    # Add HLA- prefix to class I alleles (comma-separated within each field)
    a_pref=$(echo "$a" | sed 's/^/HLA-/; s/,/,HLA-/g')
    b_pref=$(echo "$b" | sed 's/^/HLA-/; s/,/,HLA-/g')
    c_pref=$(echo "$c" | sed 's/^/HLA-/; s/,/,HLA-/g')

    # Class II fields (drb1, dq, dp) already have no prefix - keep as-is
    fields=("$a_pref" "$b_pref" "$c_pref" "$drb1" "$dq" "$dp")

    # Join all non-empty fields with commas
    final=""
    for f in "${fields[@]}"; do
        [ -z "$f" ] && continue
        if [ -z "$final" ]; then
            final="$f"
        else
            final="${final},${f}"
        fi
    done

    echo -e "${sample}\t${final}" >> "$OUT_TSV"
done

echo "Done: $OUT_TSV"
