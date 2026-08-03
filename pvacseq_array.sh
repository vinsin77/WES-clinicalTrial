#!/bin/bash
#SBATCH --job-name=pvacseq_array
#SBATCH --array=1-N%10            # replace N with actual row count from vcf_sample_hla_mapping.tsv (minus header)
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/workdir_ceci/pvacseq_log/pvacseq_%A_%a.out
#SBATCH --error=/workdir_ceci/pvacseq_log/pvacseq_%A_%a.err

# activate micromamba environment
export MAMBA_ROOT_PREFIX=~/micromamba
eval "$(~/micromamba/micromamba shell hook -s bash)"
micromamba activate pvactools_env

export CUDA_VISIBLE_DEVICES=""

TSV=/globalscratch/work_dir/vcf_sample_hla_mapping.tsv
VCF_DIR=/workdir_ceci/vep_snv_with_GT_and_GX
OUT_DIR=/globalscratch/work_dir/pvacseq_project/pvacseq_output
IEDB_DIR=/globalscratch/work_dir/pvacseq_project/iedb

# +1 skips the header row; SLURM_ARRAY_TASK_ID selects the data row
LINE=$(tail -n +2 "$TSV" | sed -n "${SLURM_ARRAY_TASK_ID}p")
sample=$(echo "$LINE" | cut -f1)
hla=$(echo "$LINE" | cut -f2)

echo "Processing sample: $sample"
echo "HLA alleles: $hla"

pvacseq run \
    "${VCF_DIR}/${sample}.filtered.vep.gt.expr.vcf" \
    TUMOR \
    "${hla}" \
    MHCflurry MHCflurryEL NetMHCpan NetMHCpanEL NetMHCIIpan NetMHCIIpanEL \
    "${OUT_DIR}/${sample}" \
    -e1 8,9,10,11 \
    -e2 12,13,14,15,16,17,18 \
    --pass-only \
    -t 4 \
    --iedb-install-directory "${IEDB_DIR}"

echo "Done: $sample"
