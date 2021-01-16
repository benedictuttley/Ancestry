#!/bin/bash

# Pipeline to address the following question:
# "Using the test dataset provided, identify a subset of female samples who have similar ancestry to the GBR reference genotypes in the 1000-genomes project."
# Calculate PCs using the all dataset
# Calculate means and SD for individuals that are GB

PLINK=/home/benedict/Bioinformatics/ancestry-study-coursework/plink
echo "path to plink executable: ${PLINK}"

### MODULE 1 ###
echo "<<< MODULE 1 >>>"

# STEP 1

# Enforce common strand by examining the .bim files and correcting
# To correct this we will flip the genotype to the reverse complement

# Two Datasets:
    # gbr-1000g* - reference file for GBR ancestry
    # Contains 91 individuals
    # 'Within Family ID (IID)' matches corresponding FID uniformly
    # 25,473,793 markers

# met583-test:
    # Unknown genotype dataset with origins in the USA
    # Contains 3563 individuals
    # 'Within Family ID (IID)' is 1 uniformly
    # 575,867 markers

# Data state on receipt:
    # Marker names have been updated to align with reference genome hg19
    # indels updated to have I/D coding
    # non-standard variation codes have been removed (CT/GA  -> T/A)
    # Duplicate SNPs have been removed

# Limit the data to common markers - minimum allele frequency of 5%
echo "Step 1"
$PLINK --bfile ../data/ftp.1000genomes.ebi.ac.uk/all-1000g-phase3-chrall-mac5-v2 --maf 0.05 --make-bed --out ../temp/1000g-gbr-maf5 #TODO: remove GBR references
# 5,495,334 markers remaining ... 6,405,953

# Step 2
# Identify and remove ambiguous SNPs
# These are those that have a W or S genotype
# AT, TA, CG, GC
awk '($5=="T"&&$6=="A")||($5=="A"&&$6=="T")||($5=="C"&&$6=="G")||($5=="G"&&$6=="C")' \
 ../temp/1000g-gbr-maf5.bim | \
 cut -f 2 > ../temp/1000g-gbr-maf5-ambig.exclude

$PLINK --bfile ../temp/1000g-gbr-maf5 \
--exclude ../temp/1000g-gbr-maf5-ambig.exclude \
--make-bed \
--out ../temp/1000g-gbr-maf5-noWS
# 4,659,586 markers remaining ... 5,427,892

# STEP 3
# Limit data to ancestry informative markers - using the *.aims SNP list - Filtration step
# Shows a divergence in allele frequency between populations
# all-1000g-phase3-chrall-mac5-v2.aims contains 6,079,790 aims

echo "Step 3"
$PLINK --bfile ../temp/1000g-gbr-maf5-noWS \
--extract ../data/ftp.1000genomes.ebi.ac.uk/all-1000g-phase3-chrall-mac5-v2.aims \
--make-bed \
--out ../temp/1000g-gbr-maf5-noWS-aims
# 2,539,713 markers remaining ... 3,373,982

### MODULE 2 ###
echo "<<< MODULE 2 >>>"
# Step 1
# Limit data to overlapping markers
# Merge in the test and reference data sets
# Limit to overlapping markers
# Create an *.extract file of common SNPs
echo "Step 1"

comm -1 -2 \
<(cut -f 2 ../temp/1000g-gbr-maf5-noWS-aims.bim | sort ) \
<(cut -f 2 ../data/met583-test/met583-test-qc-v9.bim  | sort ) \
>  ../temp/overlap.extract
    
$PLINK --bfile ../temp/1000g-gbr-maf5-noWS-aims \
--extract ../temp/overlap.extract \
--make-bed \
--out ../temp/1000g-gbr-maf5-noWS-aims-overlap
# 300,658 markers remaining ... 331,914 

# Step 2
# Remove regions of LRLD
# These regions should be excluded when performing PCA on genotype data
# Price et al. (2008) Long-Range LD Can Confound Genome Scans in Admixed Populations. Am. J. Hum. Genet. 86, 127-147
echo "Step 2"
 $PLINK --bfile ../temp/1000g-gbr-maf5-noWS-aims-overlap \
    --make-set ../data/high-ld-b37.ranges \
    --write-set \
    --out ../temp/1000g-gbr-maf5-noWS-aims-overlap-noLD

$PLINK --bfile ../temp/1000g-gbr-maf5-noWS-aims-overlap \
    --exclude ../temp/1000g-gbr-maf5-noWS-aims-overlap-noLD.set \
    --make-bed \
    --out ../temp/1000g-gbr-maf5-noWS-aims-overlap-noLD
# 322,925 markers remaining 

# Step 3
# Identify LD independent SNP set
    
# For each pair of SNPs, drop one if r2 > 0.2
# Window size: 50 snps
# Step size: 10 snps
echo "Step 3"
$PLINK --bfile ../temp/1000g-gbr-maf5-noWS-aims-overlap-noLD \
    --indep-pairwise 50 10 0.2 \
    --out ../temp/1000g-gbr-maf5-noWS-aims-overlap-noLD

### MODULE 3 ###
echo "<<< MODULE 3 >>>"
# Step 1
echo "Step 1"
# Limit reference dataset to ld independent markers
$PLINK --bfile ../data/ftp.1000genomes.ebi.ac.uk/all-1000g-phase3-chrall-mac5-v2 \
--extract ../temp/1000g-gbr-maf5-noWS-aims-overlap-noLD.prune.in \
--make-bed \
--out ../temp/all-1000g-phase3-chrall-mac5-v2-prune
    
    
# Limit test dataset to ld independent markers
$PLINK --bfile ../data/met583-test/met583-test-qc-v9 \
--extract ../temp/1000g-gbr-maf5-noWS-aims-overlap-noLD.prune.in \
--make-bed \
--out ../temp/met583-test-qc-v9-prune


# Step 2
# Merge test and reference datasets
echo "Step 2"
$PLINK --bfile ../temp/all-1000g-phase3-chrall-mac5-v2-prune \
--bmerge ../temp/met583-test-qc-v9-prune \
--make-bed \
--out ../temp/combined
# 66,922 markers remaining

# Step 3
    # Validate merge success

# 3+ alleles present:
    # This can occur when data to be merged is not harmonised to a common strand
    # PLINK has function to filp strand issues of a defined set of SNPs
    # Only need to fix the strand for one dataset
    # Then re-attempt the merge
echo "Step 3"
$PLINK --bfile ../temp/all-1000g-phase3-chrall-mac5-v2-prune \
--flip ../temp/combined-merge.missnp \
--make-bed \
--out ../temp/all-1000g-phase3-chrall-mac5-v2-prune-flip

$PLINK --bfile ../temp/all-1000g-phase3-chrall-mac5-v2-prune-flip \
--bmerge ../temp/met583-test-qc-v9-prune \
--make-bed \
--out ../temp/combined
# 89,822 markers remaining
### MODULE 4 ###
echo "<<< MODULE 4 >>>"

# Step 1
# Perform PCA analysis on the combined dataset
# This is a standard method used to cluster individuals by ancestry
# It is used to detect population structure


echo "Step 1"
$PLINK --bfile ../temp/combined \
--pca 10 header tabs \
--out ../results/combined_alt

### REST OF TASK TO BE PERFORMED IN R OR PYTHON ###

# Tips on creating the PCA plot
# PC1 vs PC2
# PC2 vs PC3
# PC1 vs PC3
# Color datapoints based on population and super-population
# Population variable needs to be mpped to the fid iid in the *.eigenvec file
# Need to add a test population value for the test data in the *eigenvec file


# MODULE 6
# Identify subset of female participants in the met583-test data
# Tips:
# Use the X chromosome (23 in plink) to infer sex
# Output: met683-test-female.keep
# 1==Male, 2==Female

$PLINK --bfile ../data/met583-test-qc-v9 --check-sex
awk 'BEGIN{printf("%10s%10s\n", "FID", "IID")};{ if ($4 == 2) {printf("%10s%10s\n", $1, $2)}}' plink.sexcheck > met583-test-female.keep

# Housekeeping
echo "Removing intermediate files..."
rm ../temp/1000g-gbr-maf5.*
rm ../temp/combined.*
echo "Done"