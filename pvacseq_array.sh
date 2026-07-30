#!/bin/bash
#SBATCH --job-name=pvacseq_array
#SBATCH --array=1-113%10          # 113 samples, 10 running concurrently — adjust as needed
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/pvacseq_%A_%a.out
#SBATCH --error=logs/pvacseq_%A_%a.err

# activate micromamba environment
export MAMBA_ROOT_PREFIX=~/micromamba
eval "$(~/micromamba/micromamba shell hook -s bash)"
micromamba activate pvactools_env

export CUDA_VISIBLE_DEVICES=""

TSV=/path/to/all_samples_hla.tsv
VCF_DIR=/path/to/vep_snv_with_GT_and_GX
OUT_DIR=/path/to/pvacseq_output
IEDB_DIR=~/iedb

# +1 skips the header row; SLURM_ARRAY_TASK_ID selects the data row
LINE=$(tail -n +2 "$TSV" | sed -n "${SLURM_ARRAY_TASK_ID}p")
sample=$(echo "$LINE" | cut -f1)
hla=$(echo "$LINE" | cut -f2)

pvacseq run \
    ${VCF_DIR}/${sample}.filtered.vep.gt.expr.vcf \
    TUMOR \
    ${hla} \
    MHCflurry MHCflurryEL NetMHCpan NetMHCpanEL NetMHCIIpan NetMHCIIpanEL \
    ${OUT_DIR}/${sample} \
    -e1 8,9,10,11 \
    -e2 12,13,14,15,16,17,18 \
    --pass-only \
    -t 4 \
    --iedb-install-directory ${IEDB_DIR}
