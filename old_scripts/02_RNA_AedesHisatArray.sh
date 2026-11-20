#!/bin/bash
#SBATCH --partition=acpu
#SBATCH --job-name=HisatArray
#SBATCH --output=%x.%j.out
#SBATCH --time=12:00:00
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
#SBATCH --mail-user=camcorey@colostate.edu


#create CxtHisat folder before starting

#Run this first to generate a list of files.
#ls | grep -v / | grep ".fastq.gz" | grep "R1" | grep "trimmed" > AedesTrimmed.txt

#wc -l AedesTrimmed.txt

#sbatch --array=0-<N> <script.sh>
#Where N is the number of lines in your AedesTrimmed.txt file minus 1. You can get that easily using the following command



# to run this script use sbatch --array=0-9 02_RNA_AedesHisatArray.sh



echo "[$0] $SLURM_JOB_NAME $@" # log the command line

module purge
source activate base
conda activate rna-seq

date # timestamp

filename=AedesTrimmed.txt # $1

linenum=0
while read -r line
do
    if [ $SLURM_ARRAY_TASK_ID -eq $linenum ]
    then
      TWO=${line/_R1_/_R2_}
      OUT=${TWO%R2_001.fastq.gz} #Remove end. 
      echo ${TWO}
      echo $OUT
      hisat2 --phred33 --rna-strandness RF -p 16 -x /scratch/alpine/camcorey@colostate.edu/CR_RNASeq/index/Cxt_Main/CtarK1_index -1 ${line} -2 ${TWO} -S CxtHisat/$OUT.sam
    fi
    linenum=$((linenum + 1))
done < $filename
