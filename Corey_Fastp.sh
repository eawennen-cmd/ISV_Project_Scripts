#!/bin/bash
#SBATCH --partition=acpu
#SBATCH --job-name=fastpLoop
#SBATCH --output=%x.%j.out
#SBATCH --time=5:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=17
#SBATCH --mail-type=ALL
#SBATCH --mail-user=c837142097@colostate.edu

module purge
source activate base
conda activate fastp

mkdir -vp trimmed/htmls

for FILE in *R1_001.fastq.gz; do
    echo "Processing R1: ${FILE}"
    
    # Cleaned up the redundant sed line to keep the script neat
    TWO="${FILE/_R1/_R2}"
    echo "Processing R2: ${TWO}"
    
    # FIX: Moved closing parenthesis inside the curly braces
    TRIM1="trimmed/trimmed.$(basename "${FILE}")"
    TRIM2="trimmed/trimmed.$(basename "${TWO}")"
    
    # Wrapped $cmd variables in quotes to handle any unexpected characters safely
    cmd="fastp -i ${FILE} -I ${TWO} -o ${TRIM1} -O ${TRIM2} -h trimmed/htmls/$(basename "${TWO}").html -j trimmed/htmls/$(basename "${TWO}").json -w $((SLURM_NTASKS-1))"
    
    echo "Running: ${cmd}"
    eval "${cmd}"
done