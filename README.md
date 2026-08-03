# WES-clinicalTrial
Bladder Cancer WES data analysis- Somatic mutations

# pVACseq pipeline (Neoantigen Prediction Analysis):
- The documentation is: https://pvactools.readthedocs.io/en/latest/pvacseq.html

1-Input file preparation (use: pVACseq_analysis_bash_commands):
- 
- Using WES VCF files
- Annotating VCFs with VEP
- Adding coverage data to your VCF
- Adding expression data to your VCF
    •using TPM.rds file converting from hgnc_symbol to ensembl_gene_id (use: pvacseq_converting_genes.R, pvacseq_converting_genes_for_all_samples.R)
    •adding expression column to the VCFs (use: pVACseq_analysis_bash_commands)

2-Preparing fastq files from germline BAM files (use: commands_for_bam_to_fastq_files):
 -
 - First, creating HLA region-only BAMs from the original BAM files
 - Second, creating collated.BAM files from HLA region-only BAMs
 - Then, creating fastq1 and fastq2 files from collated.BAM files
  
3- Using HLAHD (HLA typing algorithm, MHC Class I and Class II); https://pubmed.ncbi.nlm.nih.gov/28419628/ 
-
 -  Request access to HLA-HD (do this on your own laptop, not the server)
   •	Open a browser and go to: https://www.genome.med.kyoto-u.ac.jp/HLA-HD /
   •	Fill out the request form: name, email, institution, and intended use (put "academic research, neoantigen prediction"). 
   •	You'll receive an email (sometimes within minutes, sometimes up to a day or two) with a download link and possibly a password.
   •	Click the link and download hlahd.1.7.1.tar.gz to your laptop (e.g., into your Downloads folder). Note the exact filename — the version numbers           change.
 - Get the file onto the university server
   •	scp C:/Users/xxx/Downloads/hlahd.1.7.1.tar.gz xxxx@so-bctl-gpu01.hpda.ulb.ac.be:/work_dir/pvacseq/
   •	then move to CECI cluster: scp /work_dir/pvacseq/hlahd.1.7.1.tar.gz xxx@lemaitre4:/work_dir/
   •    extract the archive: tar -xzvf hlahd.1.7.1.tar.gz
 - Running HLAHD (use: HLAHD_commands)
 - Running SLURM (use: run_hlahd.sh)
   •    sbatch run_hlahd.sh
   •	squeue -u <user>
   •	sacct -j <job_id> --format=JobID,JobName,State,ExitCode | grep -v COMPLETED


