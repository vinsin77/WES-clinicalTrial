#!/bin/bash
# Parse HLA-HD final result files -> pVACseq-formatted allele table
# Takes the top candidate pair per locus, reduces to 2-field resolution,
# pairs DQ/DP alpha+beta chains, and flags highly ambiguous calls for review.

HLAHD_OUT="/globalscratch/workdir_ceci/pvacseq_project/hlahd_output"
OUT_TSV="hlahd_alleles.tsv"
FLAGGED="ambiguous_calls.txt"
AMBIG_THRESHOLD=10   # flag if more than this many candidate pairs at a locus

echo -e "sample\tA\tB\tC\tDRB1\tDQ_pairs\tDP_pairs" > "$OUT_TSV"
> "$FLAGGED"

# extract column 2 or 3 of a given locus line, reduced to 2-field resolution
get_allele() {
    local file=$1 locus=$2 col=$3
    local val
    val=$(awk -F'\t' -v loc="$locus" '$1==loc {print; exit}' "$file" | cut -f"$col")
    if [ -z "$val" ] || [ "$val" == "-" ] || [[ "$val" == *"Not typed"* ]]; then
        echo ""
    else
        echo "$val" | sed 's/^HLA-//' | awk -F: '{print $1":"$2}'
    fi
}

# count candidate pairs at a locus (ambiguity level)
count_pairs() {
    local file=$1 locus=$2
    local nf
    nf=$(awk -F'\t' -v loc="$locus" '$1==loc {print NF; exit}' "$file")
    [ -z "$nf" ] && { echo 0; return; }
    echo $(( (nf - 1) / 2 ))
}

for f in "$HLAHD_OUT"/*/result/*_final.result.txt; do
    sample=$(basename "$f" _final.result.txt)

    a1=$(get_allele "$f" A 2); a2=$(get_allele "$f" A 3)
    b1=$(get_allele "$f" B 2); b2=$(get_allele "$f" B 3)
    c1=$(get_allele "$f" C 2); c2=$(get_allele "$f" C 3)
    drb1_1=$(get_allele "$f" DRB1 2); drb1_2=$(get_allele "$f" DRB1 3)

    dqa1_1=$(get_allele "$f" DQA1 2); dqa1_2=$(get_allele "$f" DQA1 3)
    dqb1_1=$(get_allele "$f" DQB1 2); dqb1_2=$(get_allele "$f" DQB1 3)
    dpa1_1=$(get_allele "$f" DPA1 2); dpa1_2=$(get_allele "$f" DPA1 3)
    dpb1_1=$(get_allele "$f" DPB1 2); dpb1_2=$(get_allele "$f" DPB1 3)

    # Build A/B/C/DRB1 comma-separated strings (skip empty second allele = homozygous)
    a_str=$(printf "%s" "$a1"); [ -n "$a2" ] && [ "$a2" != "$a1" ] && a_str="${a_str},${a2}"
    b_str=$(printf "%s" "$b1"); [ -n "$b2" ] && [ "$b2" != "$b1" ] && b_str="${b_str},${b2}"
    c_str=$(printf "%s" "$c1"); [ -n "$c2" ] && [ "$c2" != "$c1" ] && c_str="${c_str},${c2}"
    drb1_str=$(printf "%s" "$drb1_1"); [ -n "$drb1_2" ] && [ "$drb1_2" != "$drb1_1" ] && drb1_str="${drb1_str},${drb1_2}"

    # Build paired DQ / DP strings (haplotype1 alpha-beta, haplotype2 alpha-beta)
    dq_str=""
    [ -n "$dqa1_1" ] && [ -n "$dqb1_1" ] && dq_str="${dqa1_1}-${dqb1_1}"
    if [ -n "$dqa1_2" ] && [ -n "$dqb1_2" ]; then
        pair2="${dqa1_2}-${dqb1_2}"
        [ "$pair2" != "$dq_str" ] && dq_str="${dq_str},${pair2}"
    fi

    dp_str=""
    [ -n "$dpa1_1" ] && [ -n "$dpb1_1" ] && dp_str="${dpa1_1}-${dpb1_1}"
    if [ -n "$dpa1_2" ] && [ -n "$dpb1_2" ]; then
        pair2="${dpa1_2}-${dpb1_2}"
        [ "$pair2" != "$dp_str" ] && dp_str="${dp_str},${pair2}"
    fi

    echo -e "${sample}\t${a_str}\t${b_str}\t${c_str}\t${drb1_str}\t${dq_str}\t${dp_str}" >> "$OUT_TSV"

    # Flag ambiguous loci for manual review
    for locus in A B C DRB1 DQA1 DQB1 DPA1 DPB1; do
        n=$(count_pairs "$f" "$locus")
        if [ "$n" -gt "$AMBIG_THRESHOLD" ]; then
            echo -e "${sample}\t${locus}\t${n} candidate pairs" >> "$FLAGGED"
        fi
    done
done

echo "Done. Output: $OUT_TSV"
echo "Flagged ambiguous calls (review manually): $FLAGGED"
wc -l "$FLAGGED"
