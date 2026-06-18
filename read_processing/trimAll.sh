#!/bin/bash
#SBATCH --partition=epyc-64
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=02:00:00
#SBATCH	--error=/scratch1/mruggeri/AcerHapPanel/log/%x_%a.err
#SBATCH --output=/scratch1/mruggeri/AcerHapPanel/log/%x_%a.out
#SBATCH --mail-type=END
#SBATCH --mail-user=mruggeri@usc.edu
#SBATCH --array=1-36

# enter your job environment parameters here
module load gcc/11.3.0 trimgalore/0.6.6
module load py-cutadapt/3.5
module load fastqc

config=/scratch1/mruggeri/AcerHapPanel/config.txt
export NAME=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# enter your job specific code below this line
trim_galore --nextseq 20 --adapter AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --clip_R1 2 --clip_R2 2 --retain_unpaired --paired $NAME'_R1.fastq' $NAME'_R2.fastq' 

