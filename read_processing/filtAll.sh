#!/bin/bash
#SBATCH --partition=epyc-64
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=03:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=mruggeri@usc.edu
#SBATCH	--error=/scratch1/mruggeri/AcerHapPanel/log/%x_%a.err
#SBATCH --output=/scratch1/mruggeri/AcerHapPanel/log/%x_%a.out
#SBATCH --array=1-36


# enter your job environment parameters here
config=/scratch1/mruggeri/AcerHapPanel/config.txt
export NAME=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# enter your job specific code below this line 
cat $NAME'_R1_val_1.fq' $NAME'_R1_unpaired_1.fq' | fastq_quality_filter -q 20 -p 90 > ../clean/$NAME'_R1.clean'
cat $NAME'_R2_val_2.fq' $NAME'_R2_unpaired_2.fq' | fastq_quality_filter -q 20 -p 90 > ../clean/$NAME'_R2.clean'

# now repair reads using perl script, note needs to be run in directory with clean reads so
cd /scratch1/mruggeri/AcerHapPanel/clean/
rePair.pl $NAME'_R1.clean' $NAME'_R2.clean'
