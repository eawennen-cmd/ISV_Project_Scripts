# ISV_Project_Scripts

SLURM scripts for finding insect-specific viruses (ISVs) in *Culex* (Cxt) RNA-seq data on the CURC Alpine cluster. Reads are trimmed, filtered against the host genome and Rift Valley fever virus (RVFV), and the leftover reads are aligned to three virus-class indexes.

## Pipeline overview

```
raw paired-end FASTQ
   |  01b_runFastp_RNA_HTML.sh
   v
trimmed FASTQ + fastp QC reports
   |  02_RNA_Cxthisatarray_unaligned.sh
   |    Stage 1: trimmed reads -> Cxt genome
   |    Stage 2: Cxt-unaligned pairs -> RVFV
   v
3_NeitherAligned (pairs that aligned to neither Cxt nor RVFV)
   |  03_RNA_VirusHisatArray.sh
   v
BAMs + summaries for dsRNA, ssRNA-minus, ssRNA-plus virus indexes
```

## Scripts

| Script | Type | Purpose |
|---|---|---|
| `01b_runFastp_RNA_HTML.sh` | single job | Quality-trims every R1/R2 pair with fastp and writes HTML/JSON QC reports |
| `02_RNA_Cxthisatarray_unaligned.sh` | array job | Two-stage HISAT2: host genome, then RVFV. Keeps host-aligned, RVFV-aligned, and "neither" datasets |
| `03_RNA_VirusHisatArray.sh` | array job | Aligns the "neither" pairs to dsRNA, ssRNAminus, and ssRNAplus virus indexes |
| `old_scripts/` | archive | Earlier versions, kept for reference. Not part of the current pipeline |

## Directory layout

Base directory: `/scratch/alpine/<user>/Virus_Alignments`

```
Virus_Alignments/
├── Cxt_FQ/
│   ├── *_R1_001.fastq.gz, *_R2_001.fastq.gz   raw reads (input to 01b)
│   ├── trimmed/                               fastp output + CulexTrimmed.txt
│   │   └── htmls/                             fastp .html and .json reports
│   └── slurm_script/                          script 02 output (see note below)
│       ├── 1_CxtAligned/
│       ├── 1_CxtUnaligned/
│       ├── 2_RVFVAligned/
│       ├── 3_NeitherAligned/
│       └── Summaries/
├── Index/
│   ├── CxtK1_genome_Index/DNAGenome.*.ht2
│   ├── RVFV_Indexes/RVFV_index.*.ht2
│   ├── dsRNA_Indexes/Virus_Index_dsRNA.fa.*.ht2
│   ├── ssRNAminus_Indexes/Virus_Index_ssRNAminus.fa.*.ht2
│   └── ssRNAplus_Indexes/Virus_Index_ssRNAplus.fa.*.ht2
├── dsRNA/, ssRNAminus/, ssRNAplus/            script 03 BAMs (+ .bai)
├── Summaries_03/<class>/                      script 03 summaries
└── slurmlogs/outlog, slurmlogs/errlog         SLURM logs
```

**Note on `slurm_script/`:** under `sbatch`, `basename "${BASH_SOURCE[0]}"` resolves to `slurm_script`, so script 02 writes to `Cxt_FQ/slurm_script/` rather than a folder named after the script. Script 03 reads from that path, so keep them consistent.

## Environments

- Script 01b: conda env `cellSquito` (fastp). Sets `LD_LIBRARY_PATH` to `${CONDA_PREFIX}/lib` to fix the `libisal.so.2` error.
- Scripts 02 and 03: conda env `hisat2` (HISAT2 and samtools >= 1.10, needed for `samtools view -G`).
- Each script sources `${SCRATCH}/miniconda3/etc/profile.d/conda.sh` before activating, because conda is not initialized in non-interactive SLURM jobs.

## Step 0: create log folders (once)

SLURM opens its log files before the script runs, so these must exist before the first submission:

```bash
mkdir -p /scratch/alpine/$USER/Virus_Alignments/slurmlogs/outlog
mkdir -p /scratch/alpine/$USER/Virus_Alignments/slurmlogs/errlog
mkdir -p /scratch/alpine/$USER/slurmlogs/outlog /scratch/alpine/$USER/slurmlogs/errlog
```

The last line is for script 01b, which logs to `/scratch/alpine/<user>/slurmlogs/`.

## Step 1: trim reads (01b)

**Input:** `Cxt_FQ/*_R1_001.fastq.gz` and the matching `*_R2_001.fastq.gz`

**Output:** `Cxt_FQ/trimmed/trimmed.<sample>_R1_001.fastq.gz` and `_R2_001.fastq.gz`, plus `trimmed/htmls/<sample>_R2_001.fastq.gz.html` and `.json`

```bash
sbatch 01b_runFastp_RNA_HTML.sh
```

Resources: 17 tasks (16 fastp threads), 10 h limit. R1 files with no matching R2 are skipped with an error message.

## Step 2: build the sample list

Both array scripts read `Cxt_FQ/trimmed/CulexTrimmed.txt`, one R1 filename per line. The array task ID is the 0-indexed line number.

```bash
cd /scratch/alpine/$USER/Virus_Alignments/Cxt_FQ/trimmed
ls | grep -v / | grep "_R1_001.fastq.gz" | grep "trimmed" > CulexTrimmed.txt
wc -l CulexTrimmed.txt
```

If there are N lines, use `--array=0-(N-1)`. Regenerate the list whenever trimmed files are added or removed. The list is alphabetical, so task IDs can shift to different samples. Do not mix task IDs from runs that used different versions of the list.

## Step 3: host and RVFV filtering (02)

For each sample:

1. **Stage 1:** align trimmed reads to the Cxt genome (`--rna-strandness RF`, spliced alignment on).
2. Split the raw SAM by flags. Pairs with at least one mapped mate go to `1_CxtAligned/<sample>.sam`. Pairs with both mates unmapped (flag 12) go to `1_CxtUnaligned/*.fastq.gz`.
3. **Stage 2:** align the Cxt-unaligned pairs to the RVFV index (`--no-spliced-alignment`).
4. Split again: RVFV-aligned reads go to `2_RVFVAligned/<sample>.sam`, and pairs that matched neither go to `3_NeitherAligned/*.fastq.gz`.
5. The raw SAMs are deleted after splitting.

```bash
sbatch --array=0 02_RNA_Cxthisatarray_unaligned.sh      # test one sample first
sbatch --array=1-20 02_RNA_Cxthisatarray_unaligned.sh   # then the rest
```

Resources: 16 CPUs, 18 h limit. Budget scratch space for one uncompressed raw SAM per running task.

**Count check per sample:**

- Unique read names in the `1_CxtAligned` SAM + pairs in `1_CxtUnaligned` FASTQ = total input pairs.
- RVFV SAM + `3_NeitherAligned` FASTQ = Cxt-unaligned pairs.

**Why not `--un-conc-gz`:** it also captures discordant and one-mate-aligned pairs, which would put host-mapped reads in the "unaligned" set. This pipeline splits on SAM flags instead.

**Caveat:** a viral pair with one spuriously host-mapped mate stays in the Cxt SAM and never reaches the RVFV or neither sets.

### One-time RVFV index build

RVFV has three segments (RefSeq NC_014397, NC_014396, NC_014395), all in one FASTA:

```bash
mkdir -p Index/RVFV_Indexes && cd Index/RVFV_Indexes
curl -o RVFV_all.fasta \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_014397,NC_014396,NC_014395&rettype=fasta&retmode=text"
grep ">" RVFV_all.fasta        # expect 3 sequences
hisat2-build -p 4 RVFV_all.fasta RVFV_index
ls RVFV_index.*.ht2            # expect 8 files
```

## Step 4: insect-specific virus alignment (03)

Takes `3_NeitherAligned` pairs and aligns them independently to each index (dsRNA, ssRNAminus, ssRNAplus). Every class gets the same input, so a pair can appear in more than one class.

**Output per sample and class:**

- `<class>/<sample><class>.bam` (+ `.bai`), coordinate-sorted, pairs with at least one mapped mate
- `Summaries_03/<class>/<sample><class>_summary.txt` (HISAT2 summary)
- `Summaries_03/<class>/<sample><class>_idxstats.txt` (mapped records per virus sequence, which shows which virus was hit)

```bash
sbatch --array=0 03_RNA_VirusHisatArray.sh
sbatch --array=1-20 03_RNA_VirusHisatArray.sh
```

Resources: 8 CPUs, 6 h limit. `idxstats` counts mapped read records (each mate counted), not pairs.

## Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `sbatch: error: Batch script contains DOS line breaks` | Windows line endings. Run `sed -i 's/\r$//' script.sh` |
| `gzip: unexpected end of file` / `fewer reads in file specified with -1 than in -2` | Truncated intermediate FASTQ. Check with `gzip -t file.fastq.gz`, delete the affected sample's outputs, and rerun that array task |
| Job fails instantly with no log | `slurmlogs/outlog` or `errlog` is missing. See Step 0 |
| `conda: invalid choice: 'deactivate'` | conda not initialized in the batch shell. Source `conda.sh` first (already done in these scripts) |
| Output folders recreated empty | Script run from a different location or the output path changed. Check `RUN_DIR` and the `slurm_script` note above |
| Failed array tasks | `sacct -j <jobid> -X --format=JobID,State,ExitCode,Elapsed,MaxRSS \| grep -v COMPLETED`, then rerun with `--array=<ids>` |

## Notes

- Scripts 02 and 03 stop at the first error and print `ERROR: ...`. Check the `.out` and `.err` logs in `slurmlogs/`.
- Run the scripts from `Virus_Alignments/` or anywhere, since all paths are absolute. Edit `SCRATCH` and `BASE_DIR` at the top of each script if you move the project.
