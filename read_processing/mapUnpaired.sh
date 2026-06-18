#!/bin/bash
#SBATCH --partition=epyc-64   
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2G
#SBATCH --time=01:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=mruggeri@usc.edu
#SBATCH	--error=/scratch1/mruggeri/AcerHapPanel/log/map/%x_%a.err
#SBATCH --output=/scratch1/mruggeri/AcerHapPanel/log/map/%x_%a.out
#SBATCH --array=1-36

# enter your job environment parameters here
export REF=/project/ckenkel_26/RefSeqs/AcerGenome/2023update
module load gcc/11.3.0 intel/19.0.4 bwa/0.7.17 samtools/1.17
config=/scratch1/mruggeri/AcerHapPanel/config.txt
export NAME=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# enter your job specific code below this line 
bwa mem -t 32 $REF/AcerSfitBwa ../clean/Unp_$NAME'_R1.clean_'$NAME'_R2.clean' > Unp_$NAME'.bwa.sam'

# sym scaffolds start with 'k127' so can use that to separate mapped reads
# first for host
grep -v 'k127_' Unp_$NAME'.bwa.sam' > Unp_$NAME'.host.sam'
samtools view -b -S Unp_$NAME'.host.sam' > Unp_$NAME'.host.bam'
samtools sort Unp_$NAME'.host.bam' -@ 32 -O BAM -o Unp_$NAME'.host.sorted.bam'

# now for sym
grep -v 'contig\|LG' Unp_$NAME'.bwa.sam' > Unp_$NAME'.sym.sam'
samtools view -b -S Unp_$NAME'.sym.sam' > Unp_$NAME'.sym.bam'
samtools sort Unp_$NAME'.sym.bam' -@ 32 -O BAM -o Unp_$NAME'.sym.sorted.bam'
