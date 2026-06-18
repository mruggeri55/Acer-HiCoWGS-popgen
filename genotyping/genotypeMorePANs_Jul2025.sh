##### genotyping more PAN #####
# Jul 21 2025
SCRIPTS=/scratch1/mruggeri/AcerHapPanel/scripts/genotype

# first call genotypes for PAN samples at loci in final vcf
# make a list of loci to genotype
VCF=/project/ckenkel_26/Acer_WGS/vcfs/snps.5bpind.DPqual.AC3.noCU.HWEfilt.bi.maf5_rm.outliers_LD.prune.vcf.gz
module load bcftools
bcftools query -f '%CHROM\t%POS\n' $VCF > snps.5bpind.DPqual.AC3.noCU.HWEfilt.bi.maf5_rm.outliers_LD.prune.pos

# run the following as job 

# #!/bin/bash
# #SBATCH --partition=epyc-64   
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=32
# #SBATCH --mem-per-cpu=2G
# #SBATCH --time=02:00:00
# #SBATCH --mail-type=END
# #SBATCH --mail-user=mruggeri@usc.edu
# #SBATCH --error=/scratch1/mruggeri/AcerHapPanel/log/geno/%x_%a.err
# #SBATCH --output=/scratch1/mruggeri/AcerHapPanel/log/geno/%x_%a.out
# #SBATCH --array=1-32
# 
# this script is for genotyping additional Panama samples from Vollmer paper for demo history analysis
# genotyping just for loci already in main dataset (see bcftools -R option)
# 
# module load bcftools/1.19
# config=/scratch1/mruggeri/AcerHapPanel/config_morePAN.txt
# export NAME=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
# 
# IN=/project2/ckenkel_26/Acer_WGS/bams/Vollmer_PAN_n.32_post.GATK_bam
# # bamOUT=/scratch1/mruggeri/AcerHapPanel/final_sanger/morePAN
# OUT=/scratch1/mruggeri/AcerHapPanel/vcfSanger/morePAN
# REF=/project/ckenkel_26/RefSeqs/AcerGenome/Sanger_Inst_GCA_964034985.1/GCA_964034985.1_Acer_S.fitti.fa
# POS=/scratch1/mruggeri/AcerHapPanel/vcfSanger/morePAN/snps.5bpind.DPqual.AC3.noCU.HWEfilt.bi.maf5_rm.outliers_LD.prune.pos
# 
# # note including -R to only genotype loci in main vcf
# bcftools mpileup -f $REF -q 40 -Q 20 -R $POS -Ou -A $IN/$NAME'_GATKd_sort.bam' \
# -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
# | bcftools call -m -Oz -f GQ,GP -o $OUT/$NAME'.vcf.gz' 

# module load htslib/1.19.1
# tabix -p vcf $OUT/$NAME'.vcf.gz'
sbatch -J morePAN $SCRIPTS/SingleCall_morePAN.sh

# 31 Jul 2025
# merge new PAN samples with core VCF
IN=/scratch1/mruggeri/AcerHapPanel/vcfSanger/morePAN
VCF=/project/ckenkel_26/Acer_WGS/vcfs/snps.5bpind.DPqual.AC3.noCU.HWEfilt.bi.maf5_rm.outliers_LD.prune.vcf.gz

cd $IN
echo -e `ls *.vcf.gz` $VCF | tr ' ' '\n' > mergeMorePan.list

# run the following as a job -- mergeMorePan.sh
#!/bin/bash
#SBATCH --partition=epyc-64   
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2G
#SBATCH --time=06:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=mruggeri@usc.edu
#SBATCH	--error=/scratch1/mruggeri/AcerHapPanel/log/geno/%x_%a.err
#SBATCH --output=/scratch1/mruggeri/AcerHapPanel/log/geno/%x_%a.out
module load bcftools/1.19
IN=/scratch1/mruggeri/AcerHapPanel/vcfSanger/morePAN
cd $IN
bcftools merge -m both -i IDV:sum,DP:sum,DP4:sum,AD:sum,MQ:avg -l mergeMorePan.list --force-samples -o core_morePAN_merge.vcf.gz


###### STARTING OVER #########
# Sep 17, 2025
# previously genotyped PAN at only loci in final vcf but need to re-genotype PAN at all loci
# HOWEVER going to first use the above merged file to check for clones so that I only include unique PAN

####### so relatedness using temp vcf of only PAN
module load plink2
plink2 --vcf core_morePAN_merge.vcf.gz --allow-extra-chr --make-king-table
awk '$6>0.3 {print}' *.kin0
awk '$6>0.3 {print $1}' *.kin0 > temp1 # get list of PAN inds that have a clone in dataset
awk '$6>0.3 {print $2}' *.kin0 > temp2
cat temp1 temp2 | sort -u > temp # list of PAN inds who have clones
awk -F'[.]' '{print $1}' mergeMorePan.list | head -n -1 > temp.list.all # list of all samples
grep -x -v -f temp temp.list.all > unique.PANs
grep -x -v -f ../../pop_gen/lists/PAN.list unique.PANs
shuf -n 6 unique.PANs > unique.PANs.n6 # randomly choose 6 extra PAN
grep -f unique.PANs.n6 *.kin0 # just triple checking no clones
cat unique.PANs.n6 ../../pop_gen/lists/PAN.list > n10.final.PAN.list

####### re-genotype extra 6 panama at all loci
# need to make config file /scratch1/mruggeri/AcerHapPanel/config_uniqPAN6.txt
cat -n unique.PANs.n6 > /scratch1/mruggeri/AcerHapPanel/config_uniqPAN6.txt
# add header ArrayTaskID	sample
# then can run the following -- SingleCall_uniqPAN6.sh
SCRIPTS=/scratch1/mruggeri/AcerHapPanel/scripts

# #!/bin/bash
# #SBATCH --partition=epyc-64   
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=32
# #SBATCH --mem-per-cpu=2G
# #SBATCH --time=02:00:00
# #SBATCH --mail-type=END
# #SBATCH --mail-user=mruggeri@usc.edu
# #SBATCH	--error=/scratch1/mruggeri/AcerHapPanel/log/geno/%x_%a.err
# #SBATCH --output=/scratch1/mruggeri/AcerHapPanel/log/geno/%x_%a.out
# #SBATCH --array=1-6
# 
# # this script is for genotyping additional Panama samples from Vollmer paper for demo history analysis
# # genotyping just for loci already in main dataset (see bcftools -R option)
# module load bcftools/1.19
# config=/scratch1/mruggeri/AcerHapPanel/config_uniqPAN6.txt
# export NAME=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
# 
# IN=/project2/ckenkel_26/Acer_WGS/bams/Vollmer_PAN_n.32_post.GATK_bam
# #bamOUT=/scratch1/mruggeri/AcerHapPanel/final_sanger/morePAN
# OUT=/scratch1/mruggeri/AcerHapPanel/vcfSanger/sample
# REF=/project2/ckenkel_26/RefSeqs/AcerGenome/Sanger_Inst_GCA_964034985.1/GCA_964034985.1_Acer_S.fitti.fa
# 
# cd $IN
# bcftools mpileup -f $REF -q 40 -Q 20 -Ou -A $NAME'_GATKd_sort.bam' \
# -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
# | bcftools call -m -Oz -f GQ,GP -o $OUT/$NAME'.vcf.gz' 
# 
# module load htslib/1.19.1
# tabix -p vcf $OUT/$NAME'.vcf.gz'

sbatch -J SingleCall_PAN6_allloci $SCRIPTS/genotype/SingleCall_uniqPAN6.sh

####### then merge by chromosome
# need to exclude SRR24007594 (PAN_7594) and CU21E_1049 from config file
# first need to make new list of samples to merge excluding the above and including the new PANs
grep -v "PAN_7594" vcf_list | grep -v "CU21E_1049" > vcf_list_clean
cat vcf_list_clean ../morePAN/lists/unique.PANs.n6 > vcf_morePAN.list 
# double check everything has full file name

# also note need to run script for major contigs (chromosomes) and then minor contigs (CAX) separately
# contents of mergeByCHR_Sep2025.sh
# #!/bin/bash
# #SBATCH --partition=epyc-64   
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=32
# #SBATCH --mem-per-cpu=2G
# #SBATCH --time=06:00:00
# #SBATCH --mail-type=END
# #SBATCH --mail-user=mruggeri@usc.edu
# #SBATCH	--error=/scratch1/mruggeri/AcerHapPanel/log/geno/%x_%a.err
# #SBATCH --output=/scratch1/mruggeri/AcerHapPanel/log/geno/%x_%a.out
# #SBATCH --array=1-15
# module load gcc/11.3.0 openblas/0.3.21 bcftools/1.14 htslib/1.17
# config=/scratch1/mruggeri/AcerHapPanel/vcfSanger/lists/chrom_noCAX.txt
# CHR=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
# 
# LIST=/scratch1/mruggeri/AcerHapPanel/vcfSanger/lists/vcf_morePAN.list 
# IN=/scratch1/mruggeri/AcerHapPanel/vcfSanger/sample
# OUT=/scratch1/mruggeri/AcerHapPanel/vcfSanger/byCHR_morePAN
# 
# cd $IN
# # merge sample vcfs into single vcf
# bcftools merge -r $CHR -l $LIST -m both -i IDV:sum,DP:sum,DP4:sum,AD:sum,MQ:avg -o $OUT/$CHR'_merge.vcf.gz'

sbatch -J mergeByChr $SCRIPTS/genotype/mergeByCHR_Sep2025.sh

# run CAX as separate job -- mergeByCAX_Sep2025.sh
# #!/bin/bash
# #SBATCH --partition=epyc-64   
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=32
# #SBATCH --mem-per-cpu=2G
# #SBATCH --time=01:00:00
# #SBATCH --mail-type=END
# #SBATCH --mail-user=mruggeri@usc.edu
# #SBATCH	--error=/scratch1/mruggeri/AcerHapPanel/log/geno/%x_%a.err
# #SBATCH --output=/scratch1/mruggeri/AcerHapPanel/log/geno/%x_%a.out
# 
# # merge minor contigs (CAX) across samples
# # vcf_morePAN.list includes more PAN and excludes SRR24007594 (PAN_7594) and CU21E_1049 -- change if needed
# # enter your job environment parameters here
# module load bcftools/1.19
# CHR=CAX
# LIST=/scratch1/mruggeri/AcerHapPanel/vcfSanger/lists/vcf_morePAN.list 
# IN=/scratch1/mruggeri/AcerHapPanel/vcfSanger/sample
# OUT=/scratch1/mruggeri/AcerHapPanel/vcfSanger/byCHR_morePAN
# 
# # enter your job specific code below this line 
# cd $IN
# # merge sample vcfs into single vcf
# bcftools merge -r `cat $CHR` -l $LIST -m both -i IDV:sum,DP:sum,DP4:sum,AD:sum,MQ:avg -o $OUT/$CHR'_merge.vcf.gz' 
sbatch -J mergeCAX_morePAN ../../scripts/genotype/mergeByCAX_Sep2025.sh

######### concatenate chromosomes together
# first make list of files to concatenate
ls *.vcf.gz > ../lists/cat_list
# double check same order as chromosome list -- if not reorder

# run the following as job -- concatCHR_Sep2025.sh
#!/bin/bash
#SBATCH --partition=epyc-64   
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2G
#SBATCH --time=06:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=mruggeri@usc.edu
#SBATCH	--error=/scratch1/mruggeri/AcerHapPanel/log/geno/%x_%a.err
#SBATCH --output=/scratch1/mruggeri/AcerHapPanel/log/geno/%x_%a.out
# 
# # enter your job environment parameters here
# IN=/scratch1/mruggeri/AcerHapPanel/vcfSanger/byCHR_morePAN
# OUT=/scratch1/mruggeri/AcerHapPanel/vcfSanger/final_merge/morePAN
# LIST=/scratch1/mruggeri/AcerHapPanel/vcfSanger/lists/cat_list
# 
# module load bcftools/1.19
# # enter your job specific code below this line
# cd $IN 
# bcftools concat -f $LIST -o $OUT/all_loci_merged_morePAN.vcf.gz
sbatch -J concatMorePAN ../../scripts/genotype/concatCHR_Sep2025.sh

# now pull out snps (excluding those within 5bp of indels)
# extract snps and remove those within 5bp on indels
salloc --time 3:00:00 --mem=16GB
module load bcftools/1.19
bcftools filter --SnpGap 5 all_loci_merged_morePAN.vcf.gz | bcftools view -v snps -o snps.5bpind.merged.morePAN.vcf.gz
# this is the second time I got broken pipe error but looks like it output same exact thing
bcftools query -f '%CHROM\t%POS\n' snps.5bpind.merged.morePAN.vcf.gz | wc -l #5,255,871 snps

# will run as job with different file name to check -- ok this finished
module load bcftools/1.19
bcftools filter --SnpGap 5 all_loci_merged_morePAN.vcf.gz | bcftools view -v snps -o snps.5bpind.merged.morePAN.v2.vcf.gz
# 5,255,871 snps
# 5,485,192 snps before indel filtering (see bcftools stats)

bcftools view -i 'F_MISSING<0.1' snps.5bpind.merged.morePAN.v2.vcf.gz | bcftools query -f '%CHROM\t%POS\n' | wc -l
# 4,505,949 snps after filtering for missingness
bcftools view -i 'MEAN(FMT/DP)>=10 && MEAN(FMT/DP)<=80 && QUAL>20 && GQ>20 && MQ>40 && FS<10 && F_MISSING<0.1' snps.5bpind.merged.morePAN.v2.vcf.gz | bcftools query -f '%CHROM\t%POS\n' | wc -l
# 4,291,465 snps after missingness and depth filtering
bcftools view -i 'MEAN(FMT/DP)>=10 && MEAN(FMT/DP)<=80 && QUAL>20 && GQ>20 && MQ>40 && FS<10 && F_MISSING<0.1' --max-alleles 2 snps.5bpind.merged.morePAN.v2.vcf.gz | bcftools query -f '%CHROM\t%POS\n' | wc -l
# 4,231,760 bi-allelic snps
bcftools view -i 'MEAN(FMT/DP)>=10 && MEAN(FMT/DP)<=80 && QUAL>20 && GQ>20 && MQ>40 && FS<10 && F_MISSING<0.1' --max-alleles 2 snps.5bpind.merged.morePAN.v2.vcf.gz | bcftools filter -e 'AC=0 || AC==AN' | bcftools query -f '%CHROM\t%POS\n' | wc -l
# 4,227,460 bi-alellics excluding monomorphic

# kinship again just in case but these PANs should be unique
bcftools view --max-alleles 2 \
	-i 'MEAN(FMT/DP)>=10 && MEAN(FMT/DP)<=80 && QUAL>20 && GQ>20 && MQ>40 && FS<10 && F_MISSING<0.1' -Ou snps.5bpind.merged.morePAN.v2.vcf.gz |\
	bcftools filter \
	-e 'AC=0 || AC==AN' \
	-o temp.snps.5bpind.merged.morePAN.qual.miss.depth.nomono.vcf.gz
VCF=temp.snps.5bpind.merged.morePAN.qual.miss.depth.nomono.vcf.gz
module load plink2
plink2 --vcf $VCF --allow-extra-chr --make-king-table
awk '$6>0.3 {print}' *.kin0

# get general stats
bcftools stats temp.snps.5bpind.merged.morePAN.qual.miss.depth.nomono.vcf.gz > bcf_stats_snps_post_filter.txt
bcftools stats all_loci_merged_morePAN.vcf.gz > bcf_stats_all_loci.txt

# get missingness per sample
paste \
<(bcftools query -f '[%SAMPLE\t]\n' $VCF | head -1 | tr '\t' '\n') \
<(bcftools query -f '[%GT\t]\n' $VCF | awk -v OFS="\t" '{for (i=1;i<=NF;i++) if ($i == "./.") sum[i]+=1 } END {for (i in sum) print i, sum[i] / NR }' | sort -k1,1n | cut -f 2) \
>> miss_by_sample_post_filter.txt
# great do not see any bias by missingness between PAN and others and all realy low (but note already filtered loci with >10% missingness)
VCF=snps.5bpind.merged.morePAN.v2.vcf.gz
paste \
<(bcftools query -f '[%SAMPLE\t]\n' $VCF | head -1 | tr '\t' '\n') \
<(bcftools query -f '[%GT\t]\n' $VCF | awk -v OFS="\t" '{for (i=1;i<=NF;i++) if ($i == "./.") sum[i]+=1 } END {for (i in sum) print i, sum[i] / NR }' | sort -k1,1n | cut -f 2) \
>> miss_by_sample.txt

# depth per sample
VCF=temp.snps.5bpind.merged.morePAN.qual.miss.depth.nomono.vcf.gz
vcftools --gzvcf $VCF --depth --out mean_depth_per_sample_post_filter.txt
cat mean_depth_per_sample_post_filter.txt.idepth | grep "SRR" | awk '{sum +=$3} END {print sum/NR}' 
# 47.1187 for PAN only
cat mean_depth_per_sample_post_filter.txt.idepth | grep -v "SRR" | awk '{sum +=$3} END {print sum/NR}'
# 30.1171 for newly sequenced samples
# sftp depth and missingness files to plot in R

