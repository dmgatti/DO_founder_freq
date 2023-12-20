#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10
#SBATCH --mem 16G
#SBATCH --time 0-6:00
#SBATCH --array 1-20

################################################################################
# Using the 48 seqWell DO low-coverage sequencing samples to call variants
# in the DO. Combine the sample BAMs into a single BAM file. Get thet set 
# diffference of the Sanger SNPs and the called variants. 
# 
# Daniel Gatti
# dan.gatti@jax.org
# 2023-11-16
################################################################################

set -e -u -o pipefail

##### VARIABLES #####

CHR=${SLURM_ARRAY_TASK_ID}

if [ ${CHR} -eq 20 ]
then
  CHR=X
fi

# Base directory for project.
BASE_DIR=/compsci/gedi/DO_founder_freq

# BAM directory.
BAM_DIR=${BASE_DIR}/data/bams

# List of BAM files.
BAM_FILES=`ls ${BAM_DIR}/*.bam`

# Output directory.
OUT_DIR=${BASE_DIR}/results/do_variants

# Combined BAM file.
COMBINED_BAM=${OUT_DIR}/seqwell_do_combined.bam

# Singularity container directory.
CONTAINER_DIR=${BASE_DIR}/containers

# Samtools container.
SAMTOOLS=${CONTAINER_DIR}/quay.io-biocontainers-samtools-1.14--hb421002_0.img

# Bcftools container.
BCFTOOLS=${CONTAINER_DIR}/quay.io-biocontainers-bcftools-1.10.2--h4f4756c_3.img

# Mouse GRCm39 reference FASTA file.
REF_FASTA=/compsci/gedi/reference/Mus_musculus.GRCm39.dna.primary_assembly.fa

# Number of threads to use.
N_THREADS=10

# Output file for this chromosome.
VCF_FILE=${OUT_DIR}/do_variants_chr${CHR}.vcf.gz


##### MAIN #####

echo CHR ${CHR}

module load singularity

# Create the output directory.
mkdir -p ${OUT_DIR}

# Create file with paths to BAM files.
ls -1 ${BAM_DIR}/*sorted.bam > bamlist.txt

# Merge the BAM files.
#singularity  exec -B /compsci ${SAMTOOLS} samtools merge \
#                                         -b bamlist.txt \
#                                         -o ${COMBINED_BAM} \
#                                         --threads 8
# Sort the BAM file.
#singularity exec -B /compsci ${SAMTOOLS} samtools sort \
#                                      -O bam \
#                                      -o 
#                                      ${COMBINED_BAM}

# Index the BAM file.
#singularity exec -B /compsci ${SAMTOOLS} samtools index \
#                                         -@ 8 \
#                                         ${COMBINED_BAM}

# Sort and index the sample BAMs.
#for BAM in ${BAM_FILES}
#do

  # Sort the BAM file.
#  singularity exec -B /compsci ${SAMTOOLS} samtools sort \
#                                           -O bam \
#                                           -o ${BAM/skipped_downsampled/sorted} \
#                                           ${BAM}
  # Index the BAM file.
#  singularity exec -B /compsci ${SAMTOOLS} samtools index \
#                                           -@ 8 \
#                                           ${BAM/skipped_downsampled/sorted}

#done

# List of sorted BAM files.
SORTED_BAM_FILES=`ls ${BAM_DIR}/*sorted.bam`

# Call variants in the individual BAM files.
singularity exec -B /compsci ${BCFTOOLS} bcftools mpileup \
                                         --fasta-ref ${REF_FASTA} \
                                         --output-type u \
                                         --max-depth 1000 \
                                         --threads ${N_THREADS} \
                                         --ignore-RG \
                                         --regions ${CHR}:1-200000000 \
                                         ${SORTED_BAM_FILES} |
singularity exec -B /compsci ${BCFTOOLS} bcftools call \
                                         --multiallelic-caller \
                                         --threads ${N_THREADS} \
                                         --output-type z \
                                         --output ${VCF_FILE} \
                                         --variants-only

# Index the VCF.
singularity exec -B /compsci ${SAMTOOLS} tabix --preset vcf ${VCF_FILE}

# Call variants in the combined BAM file.
VCF_FILE=${OUT_DIR}/do_variants_chr${CHR}_combined.vcf.gz

singularity exec -B /compsci ${BCFTOOLS} bcftools mpileup \
                                         --fasta-ref ${REF_FASTA} \
                                         --output-type u \
                                         --max-depth 1000 \
                                         --threads ${N_THREADS} \
                                         --ignore-RG \
                                         --regions ${CHR}:1-200000000 \
                                         ${COMBINED_BAM} |
singularity exec -B /compsci ${BCFTOOLS} bcftools call \
                                         --multiallelic-caller \
                                         --threads ${N_THREADS} \
                                         --output-type z \
                                         --output ${VCF_FILE} \
                                         --variants-only

