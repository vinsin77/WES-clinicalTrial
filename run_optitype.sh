#!/bin/bash
#SBATCH --job-name=optitype
#SBATCH --array=1-103
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=batch
#SBATCH --output=/my_dir/logs/optitype_%A_%a.out
#SBATCH --error=/my_dir/logs/optitype_%A_%a.err

# Load modules
module load BWA/0.7.18-GCCcore-13.3.0
module load SAMtools/1.21-GCC-13.3.0
module load GLPK/5.0-GCCcore-13.3.0
module load Python/3.12.3-GCCcore-13.3.0

# Set PATH
export PATH=/my_dir/.local/bin:$PATH

# Get sample name from array index
#SAMPLE=$(ls /my_dir/pvacseq_hla_bams/hla_reads_*_1.fastq | sed -n "${SLURM_ARRAY_TASK_ID}p" | basename | sed 's/hla_reads_//' | sed 's/_1.fastq//')
FILES=(/my_dir/pvacseq_hla_bams/hla_reads_*_1.fastq)

FILE=${FILES[$((SLURM_ARRAY_TASK_ID-1))]}
SAMPLE=$(basename "$FILE" _1.fastq)
SAMPLE=${SAMPLE#hla_reads_}
# Create output directory
#mkdir -p /my_dir/optitype_output/${SAMPLE}
# Create output directory
mkdir -p /my_dir/optitype_output/${SAMPLE}

# Run OptiType
optitype run \
  -i /my_dir/pvacseq_hla_bams/hla_reads_${SAMPLE}_1.fastq \
  -i /my_dir/pvacseq_hla_bams/hla_reads_${SAMPLE}_2.fastq \
  --dna \
  -o /my_dir/optitype_output/${SAMPLE} \
  --prefix ${SAMPLE}
