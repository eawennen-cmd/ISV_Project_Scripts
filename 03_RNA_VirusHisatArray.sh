#!/bin/bash
#SBATCH --partition=acpu
#SBATCH --job-name=HisatVirus
#SBATCH --output=/scratch/alpine/c837142097@colostate.edu/Virus_Alignments/slurmlogs/outlog/%x.%A_%a.out
#SBATCH --error=/scratch/alpine/c837142097@colostate.edu/Virus_Alignments/slurmlogs/errlog/%x.%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=c837142097@colostate.edu

# =============================================================================
# 03_RNA_VirusHisatArray.sh
#
# Purpose:
#   SLURM array job. Takes the "Neither" read pairs from script 02 (reads that
#   aligned to neither Cxt nor RVFV) and aligns them, independently, to each of
#   three insect-specific virus indexes:
#       dsRNA, ssRNAminus, ssRNAplus
#   Every virus class gets the SAME input, so a pair can show up in more than
#   one class (e.g. a conserved region). Classes are not filtered against each
#   other.
#
# Input:
#   ${TRIM_DIR}/CulexTrimmed.txt   same sample list used by script 02, so
#                                  array task ID -> sample is identical
#   ${NEITHER_DIR}/<sample>unaligned_R1.fastq.gz / _R2.fastq.gz
#
# Output (under ${RUN_DIR} = ${BASE_DIR}/Data/VirusAlignments):
#   <class>/<sample><class>.bam (+ .bai)   coordinate-sorted; pairs with at
#                                          least one mapped mate (both mates
#                                          unmapped pairs are dropped)
#   Summaries/<class>/<sample><class>_summary.txt    HISAT2 alignment summary
#   Summaries/<class>/<sample><class>_idxstats.txt   mapped records per virus
#                                          sequence -> shows WHICH virus hit
#
# Usage:
#   sbatch --array=0-20 03_RNA_VirusHisatArray.sh
#   N = (lines in CulexTrimmed.txt) - 1. Test with --array=0 first.
#
# Notes:
#   - RUN_DIR is hard-coded (not derived from the script name). Under sbatch
#     the script name resolves to "slurm_script", which would mix these outputs
#     into the folder script 02 already uses.
#   - Index prefixes include ".fa" (files are Virus_Index_dsRNA.fa.1.ht2 etc.),
#     so the -x prefix is Virus_Index_dsRNA.fa
#   - --no-spliced-alignment is used (viral genomes have no introns).
#   - idxstats counts mapped READ RECORDS (each mate counted), not pairs.
#   - Requires samtools >= 1.10 in hisat2_env (view -G).
# =============================================================================

echo "[$0] $SLURM_JOB_NAME $@"
module purge

# --- Load conda ---
SCRATCH=/scratch/alpine/c837142097@colostate.edu
source ${SCRATCH}/miniconda3/etc/profile.d/conda.sh || { echo "conda.sh not found"; exit 1; }
conda activate hisat2 || { echo "activate failed"; exit 1; }
date

# --- Paths ---
BASE_DIR="/scratch/alpine/c837142097@colostate.edu/Virus_Alignments"
TRIM_DIR="${BASE_DIR}/Cxt_FQ/trimmed"
NEITHER_DIR="${BASE_DIR}/Cxt_FQ/slurm_script/3_NeitherAligned"
RUN_DIR="${BASE_DIR}"
SUM_DIR="${RUN_DIR}/Summaries_03"
TMP_DIR="${RUN_DIR}/sorttmp"

VIRUS_IDX_DIR="/scratch/alpine/c837142097@colostate.edu/Virus_Alignments/Index"

# Virus class -> HISAT2 index prefix
declare -A IDX
IDX[dsRNA]="${VIRUS_IDX_DIR}/dsRNA_Indexes/Virus_Index_dsRNA.fa"
IDX[ssRNAminus]="${VIRUS_IDX_DIR}/ssRNAminus_Indexes/Virus_Index_ssRNAminus.fa"
IDX[ssRNAplus]="${VIRUS_IDX_DIR}/ssRNAplus_Indexes/Virus_Index_ssRNAplus.fa"
CLASSES=(dsRNA ssRNAminus ssRNAplus)

THREADS=${SLURM_CPUS_PER_TASK:-8}
SORT_THREADS=4

# --- Sanity checks ---
filename="${TRIM_DIR}/CulexTrimmed.txt"

[ -d "${NEITHER_DIR}" ] || { echo "ERROR: NEITHER_DIR does not exist: ${NEITHER_DIR}"; exit 1; }
[ -f "${filename}" ]    || { echo "ERROR: Sample list not found: ${filename}"; exit 1; }

for CLASS in "${CLASSES[@]}"; do
    P="${IDX[$CLASS]}"
    if [ ! -f "${P}.1.ht2" ] && [ ! -f "${P}.1.ht2l" ]; then
        echo "ERROR: ${CLASS} index not found or not readable at prefix: ${P}"
        echo "       Check that c837142097 has granted you read access, or copy the index folder to your scratch."
        exit 1
    fi
done

command -v samtools >/dev/null || { echo "ERROR: samtools not in hisat2_env"; exit 1; }
samtools view --help 2>&1 | grep -q -- '-G' || { echo "ERROR: samtools too old (need >= 1.10 for view -G)"; exit 1; }

mkdir -p "${TMP_DIR}"
for CLASS in "${CLASSES[@]}"; do
    mkdir -vp "${RUN_DIR}/${CLASS}" "${SUM_DIR}/${CLASS}"
done

linenum=0
while read -r line
do
    if [ "$SLURM_ARRAY_TASK_ID" -eq "$linenum" ]
    then
      TWO=${line/_R1_/_R2_}
      OUT=${TWO%R2_001.fastq.gz}   # sample prefix (keeps trailing "_"), same as script 02
      echo "${OUT}"

      R1="${NEITHER_DIR}/${OUT}unaligned_R1.fastq.gz"
      R2="${NEITHER_DIR}/${OUT}unaligned_R2.fastq.gz"
      [ -f "${R1}" ] && [ -f "${R2}" ] \
        || { echo "ERROR: Neither FASTQs missing for ${OUT} in ${NEITHER_DIR}"; exit 1; }

      for CLASS in "${CLASSES[@]}"; do
          echo "--- ${CLASS} ---"
          BAM="${RUN_DIR}/${CLASS}/${OUT}${CLASS}.bam"

          # HISAT2 writes SAM to stdout; keep pairs with >=1 mapped mate, sort to BAM
          set -o pipefail
          hisat2 --phred33 --rna-strandness RF --no-spliced-alignment -p ${THREADS} \
            -x "${IDX[$CLASS]}" \
            -1 "${R1}" -2 "${R2}" \
            --summary-file "${SUM_DIR}/${CLASS}/${OUT}${CLASS}_summary.txt" \
            | samtools view -u -G 12 - \
            | samtools sort -@ ${SORT_THREADS} -T "${TMP_DIR}/${OUT}${CLASS}" -o "${BAM}" - \
            || { echo "ERROR: ${CLASS} alignment failed for ${OUT}"; exit 1; }

          samtools index "${BAM}" || { echo "ERROR: index failed for ${BAM}"; exit 1; }
          samtools idxstats "${BAM}" > "${SUM_DIR}/${CLASS}/${OUT}${CLASS}_idxstats.txt"

          MAPPED=$(awk '{s+=$3} END{print s+0}' "${SUM_DIR}/${CLASS}/${OUT}${CLASS}_idxstats.txt")
          echo "${OUT} ${CLASS}: ${MAPPED} mapped read records"
      done

      date
    fi
    linenum=$((linenum + 1))
done < "$filename"
