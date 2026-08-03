#!/bin/bash
#SBATCH --job-name=hlahd
#SBATCH --array=1-N%K              # replace N with sample count; %K throttles concurrent tasks (optional, e.g. %10)
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=batch
#SBATCH --output=/workdir_bctl/logs_hlahd/hlahd_%A_%a.out
#SBATCH --error=/workdir_bctl/logs_hlahd/hlahd_%A_%a.err

# Load modules
module load SAMtools/1.21-GCC-13.3.0
module load GCC/13.3.0
# Activate micromamba env for bowtie2
source ~/.bashrc
micromamba activate pvactools_env

# Add HLA-HD to PATH
export PATH=/workdir_bctl/hlahd.1.7.1/bin:$PATH

# Get sample name from array index
FILES=(/globalscratch/workdir_ceci/pvacseq_project/pvacseq_hla_bams/hla_reads_*_1.fastq)
FILE=${FILES[$((SLURM_ARRAY_TASK_ID-1))]}
SAMPLE=$(basename "$FILE" _1.fastq)
SAMPLE=${SAMPLE#hla_reads_}

echo "Processing sample: $SAMPLE"

# Run HLA-HD
hlahd.sh -t 8 -m 100 \
  /globalscratch/workdir_ceci/pvacseq_project/pvacseq_hla_bams/hla_reads_${SAMPLE}_1.fastq \
  /globalscratch/workdir_ceci/pvacseq_project/pvacseq_hla_bams/hla_reads_${SAMPLE}_2.fastq \
  /workdir_bctl/hlahd.1.7.1/HLA_gene.split.txt \
  /workdir_bctl/hlahd.1.7.1/dictionary \
  ${SAMPLE} \
  /globalscratch/workdir_ceci/pvacseq_project/hlahd_output

echo "Done: $SAMPLE"
