# redoing EEMS after adding more PAN
salloc --time 2:00:00 --mem=16GB
# new vcf after removing PAN clone -- "neutral" so will have removed HWE outliers
# also should have already filtered for physical linkage, multi-allelics, and maf
VCF=/project2/ckenkel_26/Acer_WGS/vcfs/Sep2025_morePAN/Acer_main.NEU.PANEL_n46.vcf.gz

# make diffs file
module load plink2
plink2 --vcf $VCF --allow-extra-chr --make-bed --out snpsForPopGen_plink
# also need to fix chromosome codes
awk '{$1="0";print $0}' snpsForPopGen_plink.bim > snpsForPopGen_plink.bim.tmp
mv snpsForPopGen_plink.bim.tmp snpsForPopGen_plink.bim
# now run bed2diffs to make diff file
bed2diffs_v1 --bfile snpsForPopGen_plink

# get info for params
# how many loci?
module load bcftools
bcftools query -Hf '%CHROM\t%POS\n' $VCF | wc -l 
#228059 loci
# copy and modify coord file to include extra panama samples
cp ../sample_to_coord.txt .
cat sample_to_coord.txt | awk '{print $3,$4}' > snpsForPopGen_plink.coord
# copy over outer file -- this is the same, just to draw outer polygon
cp ../snpsForPopGen_plink.outer .
# make output directory
mkdir out
# make param file and add/edit the following
nano param.txt
# datapath = snpsForPopGen_plink
# mcmcpath = out
# nIndiv = 46
# nSites = 228059
# nDemes = 500
# diploid = true
# numMCMCIter = 5000000
# numBurnIter = 1000000
# numThinIter = 9999

# run EEMS.sh
# runeems_snps --params param.txt --seed 123
sbatch -J EEMS EEMS.sh
# took a 1.2-2.2- hours

# note out is d500
# running d400 and d300 and naming out dir deme size in param file

# plot results in R
module load geos gdal proj udunits r
R

library(rEEMSplots)
setwd('/scratch1/mruggeri/AcerHapPanel/pop_gen/EEMS/morePAN/')

mcmcpath = "out/"
plotpath = "out/d500"
projection_none <- "+proj=longlat +datum=WGS84"
projection_mercator <- "+proj=merc +datum=WGS84"

eems.plots(mcmcpath, plotpath, longlat = F, projection.in = projection_none, projection.out = projection_mercator, out.png=FALSE, add.grid=F, add.demes=T, add.outline=T, add.map=T, add.abline = T, add.r.squared = T)

q()

