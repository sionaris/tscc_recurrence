#!/bin/bash
#SBATCH --time=23:59:59
#SBATCH --partition=batch
#SBATCH --nodes=1                # Use 5 nodes
#SBATCH --ntasks=1               # 1 task per node
#SBATCH --cpus-per-task=12       # 12 CPU per task
#SBATCH --mem=100G               # 100 GB memory for the node
#SBATCH --job-name=tscc_rnaseq   # Job name
#SBATCH --output=logs/%x_%A_%a.out # Standard output
#SBATCH --error=logs/%x_%A_%a.err  # Standard error
#SBATCH --array=0-79             # Array jobs to process batches (80 batches to cover 80 samples)

cd path/to/project

# Load required modules
module load gcc/13.2.0-iqpfkya openjdk/11.0.20 samtools/1.19.2 apptainer

# Define the number of samples per batch
SAMPLES_PER_BATCH=1  # 1 sample per batch to fit the memory constraints

# Calculate start and end index for this array job
BATCH_START=$(( SLURM_ARRAY_TASK_ID * SAMPLES_PER_BATCH + 2 ))  # +2 to skip the header
BATCH_END=$(( (SLURM_ARRAY_TASK_ID + 1) * SAMPLES_PER_BATCH + 1 ))  # +1 to include the correct range

# Create a temporary CSV file for the current batch
TEMP_CSV=$(mktemp batch_${SLURM_ARRAY_TASK_ID}_XXXX.csv)

# Add the header to the temporary CSV file
echo "sample,fastq_1,strandedness" > $TEMP_CSV

# Extract the batch samples from samplesheet.csv and append to the temporary CSV file
sed -n "${BATCH_START},${BATCH_END}p" samplesheet.csv >> $TEMP_CSV

# Set apptainer variables
APPTAINER_CACHEDIR=path/to/apptainer_cache
NXF_APPTAINER_CACHEDIR=path/to/apptainer_cache
NXF_APPTAINER_HOME_MOUNT=true

# Run Nextflow pipeline on the current batch using the temporary CSV file
path/to/nextflow run \
    nf-core/rnaseq \
    --input $TEMP_CSV \
    --outdir NF-processed_samples/batch_${SLURM_ARRAY_TASK_ID} \
    --gtf Homo_sapiens.GRCh38.112.gtf \
    --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --igenomes_ignore \
    --genome 'GRCh38' \
    -profile apptainer

# Clean up the temporary CSV file
rm $TEMP_CSV
