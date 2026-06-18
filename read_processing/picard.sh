#!/bin/bash
#SBATCH --partition=epyc-64   
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=mruggeri@usc.edu
#SBATCH	--error=/scratch1/mruggeri/AcerHapPanel/log/dup/%x_%a.err
#SBATCH --output=/scratch1/mruggeri/AcerHapPanel/log/dup/%x_%a.out
#SBATCH --array=1-36

module load gcc/11.3.0 samtools/1.17 picard/2.26.2
#config=/scratch1/mruggeri/AcerHapPanel/config.txt
#export NAME=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
export NAME=AR_113

# merge paired and single end bams
#samtools merge -o $NAME'.host.merge.bam' $NAME'.host.sorted.bam' Unp_$NAME'.host.sorted.bam'
samtools merge -o $NAME'.sym.merge.bam' $NAME'.sym.sorted.bam' Unp_$NAME'.sym.sorted.bam'

# reformat files for GATK
#java -Xmx4G -jar $PICARD AddOrReplaceReadGroups I=$NAME'.host.merge.bam' O=$NAME'.host.rg.bam' RGID=group1 RGLB=lib1 RGPU=unit1 RGPL=illumina RGSM=$NAME
java -Xmx4G -jar $PICARD AddOrReplaceReadGroups I=$NAME'.sym.merge.bam' O=$NAME'.sym.rg.bam' RGID=group1 RGLB=lib1 RGPU=unit1 RGPL=illumina RGSM=$NAME

# mark PCR duplicates
#java -Xmx6G -jar $PICARD MarkDuplicates I=$NAME'.host.rg.bam' O=../final_bams/$NAME'.host.marked_duplicates.bam' M=../log/dup/$NAME'.host.marked_dup_metrics.txt'
java -Xmx6G -jar $PICARD MarkDuplicates I=$NAME'.sym.rg.bam' O=../final_bams/$NAME'.sym.marked_duplicates.bam' M=../log/dup/$NAME'.sym.marked_dup_metrics.txt'

# index file
#samtools index ../final_bams/$NAME'.host.marked_duplicates.bam'
samtools index ../final_bams/$NAME'.sym.marked_duplicates.bam'
