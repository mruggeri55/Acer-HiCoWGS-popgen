##### calculating diversity stats for haplotype panel ######
  
#####################################
# LD decay -- for populations
####################################
# 7 Oct 2025 -- running for K of 5 with more PAN included
# will subsample each pop to 8 individuals so sample size does not effect LD results

# first create a soft filtered vcf with snps only -- filtered by depth, missingness, quality, bi-allelics
# note going with mean max depth 50 for now after looking at depth distribution
salloc
module load bcftools/1.19
VCF=/scratch1/mruggeri/AcerHapPanel/vcfSanger/final_merge/morePAN/snps.5bpind.merged.morePAN.v2.vcf.gz
bcftools view --max-alleles 2 \
	-i 'MEAN(FMT/DP)>=10 && MEAN(FMT/DP)<=50 && QUAL>20 && GQ>20 && MQ>40 && FS<10 && F_MISSING<0.1' -Ou $VCF |\
	bcftools filter \
	-e 'AC=0 || AC==AN' \
	-Oz --write-index \
	-o snps.5bpind.minmeanDP10.maxmeanDP50.qual.miss10.nomono.vcf.gz
# make list of samples
lists=/scratch1/mruggeri/AcerHapPanel/pop_gen/lists
bcftools view -h $VCF | grep "CHROM" | sed 's/\t/\n/g' | tail -n +10 > morePAN.list
# also need to change sample names because plink won't allow underscores or dashes
sed 's/_//g' $lists/morePAN.list | sed 's/-//g' > $lists/morePAN_4plink.txt
VCF=/scratch1/mruggeri/AcerHapPanel/vcfSanger/final_merge/morePAN/snps.5bpind.minmeanDP10.maxmeanDP50.qual.miss10.nomono.vcf.gz
bcftools reheader --samples $lists/morePAN_4plink.txt -o snps.5bpind.minmeanDP10.maxmeanDP50.qual.miss10.nomono.4plink.vcf.gz $VCF
bcftools view -h snps.5bpind.minmeanDP10.maxmeanDP50.qual.miss10.nomono.4plink.vcf.gz | grep "CHROM"
bcftools view -h $VCF | grep "CHROM"
bcftools query -f '%CHROM\t%POS\n' $VCF | wc -l #3,956,510 snps
# make new lists for each population
cd $lists
# removing admixed individual and individual grouping with MEX/BZ from FL
grep -v CRF_Acer-099 FLkeys.list | grep -v CRF_Acer-102 > FLkeys.pop.list 
# adding in FL sample grouping with MEX/BZ
sed '$a\CRF_Acer-102' MEX_BZ.list > MEX_BZ.pop.list
# full list of Panama samples
grep "SRR" morePAN.list > PAN.pop.list
# also make pop list by copying others
cp AR_CU.list AR_CU.pop.list
cp DR_JAM.list DR_JAM.pop.list

#!/bin/bash
#SBATCH --partition=epyc-64  
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --error=/scratch1/mruggeri/AcerHapPanel/pop_gen/LD/%x_%j.err
#SBATCH --output=/scratch1/mruggeri/AcerHapPanel/pop_gen/LD/%x_%j.out
#SBATCH --mail-type=END
#SBATCH --mail-user=mruggeri@usc.edu
#SBATCH --array=1-5
# enter your job specific code below this line 
config=/scratch1/mruggeri/AcerHapPanel/pop_gen/lists/config_region.txt
export POP=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
module load plink2
VCF=/scratch1/mruggeri/AcerHapPanel/vcfSanger/final_merge/morePAN/snps.5bpind.minmeanDP10.maxmeanDP50.qual.miss10.nomono.4plink.vcf.gz
OUT=/scratch1/mruggeri/AcerHapPanel/pop_gen/LD/by_pop
lists=/scratch1/mruggeri/AcerHapPanel/pop_gen/lists
# reformat sample names because cannot have underscores or dashes in plink (already reformatted in vcf above)
cat $lists/$POP'.pop.list' | sed 's/_//g' | sed 's/-//g' | awk '{print $0"\t"$0}' | shuf -n 8 > $lists/$POP'.pop.list8.4plink'
plink --vcf $VCF --keep $lists/$POP'.pop.list8.4plink' \
  --allow-extra-chr \
  --ld-window 999999 \
  --ld-window-kb 100 \
  --ld-window-r2 0 \
  --out $OUT/$POP'_ld_decay' \
  --r2 \
  --thin 0.01 
# to calculate distance between positions
# col 1 of ld summary is dist and col 2 is r2
# will also sort in increasing order    
cat $OUT/$POP'_ld_decay.ld' | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' | \ 
sort -k1,1n > $OUT/$POP'_ld.summary'
# then in R get mean r2 for 100bp intervals
module load r
R
library(dplyr)
library(stringr)
library(tools)
files <- list.files(pattern = "\\.summary$")
for (file in files){
	dfr <- read.delim(file,sep="",check.names=F,stringsAsFactors=F)
	colnames(dfr) <- c("dist","rsq")
	# get mean r2 from 100bp intervals
	dfr$distc <- cut(dfr$dist,breaks=seq(from=min(dfr$dist)-1,to=max(dfr$dist)+1,by=100))
	dfr1 <- dfr %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
	dfr1 <- dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))
	write.csv(dfr1,paste(file_path_sans_ext(file),"r2_decay_100bp_sum.csv",sep="_"))
	}
q()
# now concatenate all together and add in filename
for f in `ls *_ld_r2_decay_100bp_sum.csv`; do
	awk -v filename="$(basename "$f")" 'FNR > 1 {print $0 "," filename}' "$f" >> ld_r2_decay_100bp_sum_by_pop.csv ;
	done
# added header using nano then sftped to computer to visualize in R

###### PIXY ######## -- the following is adapted from Carly's code
# Oct 22, 2025
# first create conda env to install pixy
mamba create --prefix /project2/ckenkel_26/mruggeri/software/conda_env/pixy
conda activate /project2/ckenkel_26/mruggeri/software/conda_env/pixy
conda install -c conda-forge pixy
conda install -c bioconda htslib

salloc --time=4:00:00 --cpus-per-task=4 --mem=32G
# need to filter VCF for SOR, DP, miss, etc
# following best practices here https://pixy.readthedocs.io/en/latest/guide/pixy_guide.html
VCFs=/project2/ckenkel_26/Acer_WGS/vcfs/Sep2025_morePAN
module load bcftools
bcftools view -i 'MEAN(FMT/DP)>=10 && MEAN(FMT/DP)<=50 && QUAL>20 && GQ>20 && MQ>40 && SOR<4 && F_MISSING<0.1' --write-index \
	$VCFs/snps.5bpind.merged.morePAN.v2_SOR.add.vcf.gz \
	-o $VCFs/snps.5bpind.merged.morePAN.minDP10.maxDP50.qual.SOR.miss10.vcf.gz
# also need to include invariant sites
# extract those and then can combine both files back together -- need to run the following as a job (ffilt_invariant.sh)
VCFs2=/scratch1/mruggeri/AcerHapPanel/vcfSanger/final_merge/morePAN
bcftools view -i 'ALT="." && MEAN(FMT/DP)>=10 && MEAN(FMT/DP)<=50 && QUAL>20 && MQ>40 && F_MISSING<0.1' \
 $VCFs2/all_loci_merged_morePAN.vcf.gz --write-index -Oz -o $VCFs2/invariant_loci_merged.morePAN.minDP10.maxDP50.qual.miss10.vcf.gz
# concatenate the two files together
bcftools concat --allow-overlaps \
$VCFs/snps.5bpind.merged.morePAN.minDP10.maxDP50.qual.SOR.miss10.vcf.gz $VCFs2/invariant_loci_merged.morePAN.minDP10.maxDP50.qual.miss10.vcf.gz \
-Oz --write-index -o $VCFs2/snps.invariant.morePAN.minDP10.maxDP50.qual.SOR.miss10.nosort.vcf.gz
#-Ou | bcftools sort --write-index -o $VCFs2/snps.invariant.morePAN.minDP10.maxDP50.qual.SOR.miss10.vcf.gz
# also make list of loci -- removing excess het 
bcftools query -f '%CHROM\t%POS\n' $VCFs2/snps.invariant.morePAN.minDP10.maxDP50.qual.SOR.miss10.nosort.vcf.gz > all_sites.txt 
# also excluding snps with excess heterozygosity (most likely sequencing errors)
rmSNPs=/project2/ckenkel_26/Acer_WGS/misc_from_PD/snps.excluded.by.RUTH.filt_morePANn46.tsv 
grep -v -f $rmSNPs all_sites.txt > sites_rmHWE.txt # note this takes a while because all_sites large file ~1.5 hrs

# need sample tab popID file for pixy
lists=/scratch1/mruggeri/AcerHapPanel/pop_gen/lists
cd $lists
awk '{filename = FILENAME; sub(/\.pop.list+$/,"",filename); print $1 "\t" filename, $2}' AR_CU.pop.list DR_JAM.pop.list FLkeys.pop.list MEX_BZ.pop.list PAN.pop.list \
> sample2pop.txt

# now run pixy 
conda activate /project2/ckenkel_26/mruggeri/software/conda_env/pixy # activate environment before running job
pixy --help # note version 2 update is not online yet so see here for stats options
pixy --version #v2.0.0

# run the following as job -- pixy.sh
OUT=/scratch1/mruggeri/AcerHapPanel/pop_gen/pixy
VCFs2=/scratch1/mruggeri/AcerHapPanel/vcfSanger/final_merge/morePAN
lists=/scratch1/mruggeri/AcerHapPanel/pop_gen/lists
sites=/scratch1/mruggeri/AcerHapPanel/pop_gen/pixy/sites_rmHWE.txt

pixy --stats pi watterson_theta tajima_d \
 --vcf $VCFs2/snps.invariant.morePAN.minDP10.maxDP50.qual.SOR.miss10.nosort.vcf.gz \
 --populations $lists/sample2pop.txt \
 --include_multiallelic_snps \
 --sites_file $sites \
 --window_size 100000 --n_cores 4 \
 --output_prefix 100kb --output_folder $OUT

# took ~3hrs

#### subsetting to 8 individuals per pop to see if patterns hold
lists=/scratch1/mruggeri/AcerHapPanel/pop_gen/lists
cd $lists
for file in AR_CU.pop.list DR_JAM.pop.list FLkeys.pop.list MEX_BZ.pop.list PAN.pop.list; do \
    prefix=$(basename "$file" .pop.list); \
    awk -v prefix="$prefix" '{print $1 "\t" prefix}' "$file" | shuf -n 8 >> sample2pop_sub8.txt; done

# activate conda environment
conda activate /project2/ckenkel_26/mruggeri/software/conda_env/pixy
# then rerun pixy using this sample file
# adding prefix 100kb_sub8

#### rerunning pixy as one metapopulation for genomic scans
awk '{print $1}' $lists/sample2pop.txt | awk -v word="global" '{print $0 "\t" word}' > $lists/sample2singlepop.txt
nano ../../scripts/pop_gen/pixy.sh
# change to 10kb windows
# output prefix 10kb_global
conda activate /project2/ckenkel_26/mruggeri/software/conda_env/pixy
sbatch -J pixy_global ../../scripts/pop_gen/pixy.sh
# also redoing with 1kb windows
sbatch -J pixy_global_1kb ../../scripts/pop_gen/pixy.sh


##### COUNT PRIVATE ALLELES #######
salloc
VCF=/project2/ckenkel_26/Acer_WGS/vcfs/Sep2025_morePAN/snps.5bpind.merged.morePAN.minDP10.maxDP50.qual.SOR.miss10.vcf.gz
VCF2=/project2/ckenkel_26/Acer_WGS/vcfs/Sep2025_morePAN/full.panel_5bpind.morePAN.HWE.DPqual.AC3.bi.SOR.vcf.gz
LIST=/scratch1/mruggeri/AcerHapPanel/pop_gen/lists

module load bcftools/1.19

bcftools view -S $LIST/AR_CU.pop.list -x $VCF | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AN\t%AC\t%AD\n' > AR_CU_private
# 193,634
# many of these only have one individual with alternate allele
bcftools query -f '%CHROM\t%POS\n' $VCF | wc -l
# 3,822,703 variants so ~5% of snps private to AR_CU
# honestly kind of a lot

bcftools view -S $LIST/AR_CU.pop.list -x $VCF2 | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AN\t%AC\n' > AR_CU_private2
# 15,534 
# note this is restricting allele count to 3 individuals
# also bi-allelics only
bcftools query -f '%CHROM\t%POS\n' $VCF2 | wc -l
# 2,527,725 variants so ~0.6% private, very different scale

# what if we post-filter one with multi-allelic to count only those where AC >2 or 3
awk '$6 >= 2 {count++} END {print count}' AR_CU_private
# 64,034
awk '$6 >= 3 {count++} END {print count}' AR_CU_private
# 25,430
# so most are unique to one individual

##### run above for each pop
salloc
VCF=/project2/ckenkel_26/Acer_WGS/vcfs/Sep2025_morePAN/snps.5bpind.merged.morePAN.minDP10.maxDP50.qual.SOR.miss10.vcf.gz
VCF2=/project2/ckenkel_26/Acer_WGS/vcfs/Sep2025_morePAN/full.panel_5bpind.morePAN.HWE.DPqual.AC3.bi.SOR.vcf.gz
LIST=/scratch1/mruggeri/AcerHapPanel/pop_gen/lists
module load bcftools/1.19

for POP in `cat $LIST/pop.list`; do 
bcftools view -S $LIST/$POP'.pop.list' -x $VCF | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AN\t%AC\t%AD\n' > $POP'_private.txt';
done
# note AD not re-calculated for pop subset, but alternative allele depth should be accurate because only present in this subpop

wc -l *.txt
#   193634 AR_CU_private.txt
#   250587 DR_JAM_private.txt
#   217095 FLkeys_private.txt
#   277012 MEX_BZ_private.txt
#   285390 PAN_private.txt

# now filter by alt allele depth of at least 10x
awk '{ split($7, a, ","); if (a[2] >= 10) print }' AR_CU_private.txt | wc -l
# 130,174

for file in `ls *.txt`; do 
awk '{ split($7, a, ","); if (a[2] >= 10) print }' $file | wc -l;
done
# 130174
# 168317
# 145482
# 193483
# 219806

# made excel file with private alleles for these different count thresholds
# going to make new column for AF since populations have different number of samples
for file in `ls *.txt`; do 
awk -v OFS='\t' '{ $(NF+1) = $6 / $5; print }' $file > $file'.af'
done
# so now can filter based on AF column
for file in `ls *.af`; do 
awk '$8 >= 0.1 {count++} END {print count}' $file;
done
# did this across different thresholds
# you can find summary here "/Users/maria/Desktop/Kenkel_lab/Acerv_CBASS_GWAS/haplotype_panel/pop_gen/2025/morePAN/private_alleles/private_allele_sum.xlsx"

####### try setting GPs<0.99 to missing
salloc
VCF=/project2/ckenkel_26/Acer_WGS/vcfs/Sep2025_morePAN/snps.5bpind.merged.morePAN.minDP10.maxDP50.qual.SOR.miss10.vcf.gz
LIST=/scratch1/mruggeri/AcerHapPanel/pop_gen/lists
module load bcftools/1.19

bcftools +setGT $VCF -- -t q -n ./. -e 'FORMAT/GP < 0.99' > /project2/ckenkel_26/Acer_WGS/vcfs/Sep2025_morePAN/snps.5bpind.merged.morePAN.minDP10.maxDP50.qual.SOR.miss10.GPfilt.vcf.gz
# Filled 263307090 alleles

VCF_GPfilt=/project2/ckenkel_26/Acer_WGS/vcfs/Sep2025_morePAN/snps.5bpind.merged.morePAN.minDP10.maxDP50.qual.SOR.miss10.GPfilt.vcf.gz
for POP in `cat $LIST/pop.list`; do 
bcftools view -S $LIST/$POP'.pop.list' -x $VCF_GPfilt | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AN\t%AC\t%AD\n' > $POP'_private_GPfilt.txt';
done
# did not change anything
# check if genotypes got set to missing
bcftools view -S $LIST/AR_CU.pop.list -x $VCF_GPfilt | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AN\t%AC\t%AD\t[%GT]\n' > temp
bcftools query -r "OZ035966.1:2198" $VCF -f '%CHROM\t%POS\t%REF\t%ALT\t%AN\t%AC\t%AD\t[%GP:%GT:%GQ:%AD]\n'
# nope, 0.99 GP can still be not well supported i.e. 1 alt / 3 reads still gets GP >0.99 

# try setting to missing if less than 10x
bcftools +setGT $VCF -- -t q -n ./. -i 'FMT/DP<10' > /project2/ckenkel_26/Acer_WGS/vcfs/Sep2025_morePAN/snps.5bpind.merged.morePAN.minDP10.maxDP50.qual.SOR.miss10.GTdpfilt.vcf.gz
# Filled 10177026 alleles
VCF_GTdpfilt=/project2/ckenkel_26/Acer_WGS/vcfs/Sep2025_morePAN/snps.5bpind.merged.morePAN.minDP10.maxDP50.qual.SOR.miss10.GTdpfilt.vcf.gz
bcftools view -S $LIST/AR_CU.pop.list -x $VCF_GTdpfilt | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AN\t%AC\t%AD\t[%GT]\n' > temp
# 176,766 priv all in AR_CU so honestly didnt reduce by much ~$20K

#### 4 Feb 2026
# filter based on alternate allele depth>10x and add in allele frequency
for file in `ls *_private.txt`; do 
awk '{ split($7, a, ","); if (a[2] >= 10) print }' $file | \
awk -v OFS='\t' '{ $(NF+1) = $6 / $5; print }' > $file'.AAD10x.af'
done
# so now can filter based on AF column
for file in `ls *.AAD10x.af`; do 
awk '$8 >= 0.1 {count++} END {print count}' $file;
done
# did this across different thresholds





