#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 10
#SBATCH --mem 128G
#SBATCH --time 0-4:00
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

# Base directory for project.
BASE_DIR=/compsci/gedi/DO_founder_freq

# BAM directory.
BAM_DIR=${BASE_DIR}/data/bams

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

##### MAIN #####

module load singularity

# Create the output directory.
mkdir -p ${OUT_DIR}

# Create file with paths to BAM files.
ls -1 ${BAM_DIR}/*.bam > bamlist.txt

# Merge the BAM files.
#singularity  exec -B /compsci ${SAMTOOLS} samtools merge \
#                                         -b bamlist.txt \
#                                         -o ${COMBINED_BAM} \
#                                         --threads 8
# Sort the BAM file.
#singularity exec ${SAMTOOLS} samtools sort \
#                                      -O bam \
#                                      -o 
#                                      ${COMBINED_BAM}

# Index the BAM file.
#singularity exec -B /compsci ${SAMTOOLS} samtools index \
#                                         -@ 8 \
#                                         ${COMBINED_BAM}

# Call variants in the combined BAM file.
singularity exec -B /compsci ${BCFTOOLS} bcftools mpileup \
                                         -Ou \
                                         -f ${REF_FASTA} \
                                         --max-depth 1000 \
                                         --threads 8 \
                                         ${COMBINED_BAM} | \
singularity exec -B /compsci ${BCFTOOLS} bcftools call \
                                         -mv \
                                         --output-type z \
                                         --output ${OUT_DIR}/seqwell_do_calls.vcf.gz \
                                         --threads 8

