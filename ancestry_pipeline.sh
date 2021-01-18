#!/bin/bash

PLINK=/home/benedict/Bioinformatics/ancestry-study-coursework/plink
echo "path to plink executable: ${PLINK}"

########## MODULE 1 START ##########
echo "<<< MODULE 1 >>>"

# STEP 1: Limit the data to common markers - minimum allele frequency of 5%
echo "[M1,S1] Remove markers with minimum allele frequency < 5%"
$PLINK --bfile ../data/ftp.1000genomes.ebi.ac.uk/all-1000g-phase3-chrall-mac5-v2 --maf 0.05 --make-bed --out ../temp/1000g-all-maf5

# STEP 2: Identify and remove ambiguous SNPs - Those that have a W or S genotypes (AT, TA, CG, GC)
echo "[M1,S2] Remove ambiguous W and S genotype SNPs"
awk '($5=="T"&&$6=="A")||($5=="A"&&$6=="T")||($5=="C"&&$6=="G")||($5=="G"&&$6=="C")' \
    ../temp/1000g-all-maf5.bim | \
    cut -f 2 > ../temp/1000g-all-maf5-ambig.exclude

$PLINK --bfile ../temp/1000g-all-maf5 \
--exclude ../temp/1000g-all-maf5-ambig.exclude \
--make-bed \
--out ../temp/1000g-all-maf5-noWS

# STEP 3: Limit data to ancestry informative markers - Filtration step
echo "[M1,S3] Limit data to ancestry informative markers"
$PLINK --bfile ../temp/1000g-all-maf5-noWS \
--extract ../data/ftp.1000genomes.ebi.ac.uk/all-1000g-phase3-chrall-mac5-v2.aims \
--make-bed \
--out ../temp/1000g-all-maf5-noWS-aims
########## MODULE 1 END ##########


########## MODULE 2 START ##########
echo "<<< MODULE 2 >>>"

# STEP 1: Identify SNPs common to both reference and test datasets 
echo "[M2,S1] Identify SNPs common to both reference and test datasets "
comm -1 -2 \
    <(cut -f 2 ../temp/1000g-all-maf5-noWS-aims.bim | sort ) \
    <(cut -f 2 ../data/met583-test/met583-test-qc-v9.bim  | sort ) >  ../temp/overlap.extract
    
$PLINK --bfile ../temp/1000g-all-maf5-noWS-aims \
    --extract ../temp/overlap.extract \
    --make-bed \
    --out ../temp/1000g-all-maf5-noWS-aims-overlap

# STEP 2: Remove regions of LRLD - Price et al. (2008) Long-Range LD Can Confound Genome Scans in Admixed Populations. Am. J. Hum. Genet. 86, 127-147
echo "[M2,S2] Remove SNPs in regions of long-range linkage disequilibrium"
 $PLINK --bfile ../temp/1000g-all-maf5-noWS-aims-overlap \
    --make-set ../data/high-ld-b37.ranges \
    --write-set \
    --out ../temp/1000g-all-maf5-noWS-aims-overlap-noLD

$PLINK --bfile ../temp/1000g-all-maf5-noWS-aims-overlap \
    --exclude ../temp/1000g-all-maf5-noWS-aims-overlap-noLD.set \
    --make-bed \
    --out ../temp/1000g-all-maf5-noWS-aims-overlap-noLD


# STEP 3: Identify LD independent SNP set
echo "[M2,S3] Identify LD independent SNP set"
$PLINK --bfile ../temp/1000g-all-maf5-noWS-aims-overlap-noLD \
    --indep-pairwise 50 10 0.2 \
    --out ../temp/1000g-all-maf5-noWS-aims-overlap-noLD
########## MODULE 2 END ##########


########## MODULE 3 START ##########
echo "<<< MODULE 3 >>>"

# STEP 1: Limit reference dataset to ld independent markers
echo "[M3,S1] Limit referenece dataset to ld independent markers"
$PLINK --bfile ../data/ftp.1000genomes.ebi.ac.uk/all-1000g-phase3-chrall-mac5-v2 \
--extract ../temp/1000g-all-maf5-noWS-aims-overlap-noLD.prune.in \
--make-bed \
--out ../temp/all-1000g-phase3-chrall-mac5-v2-prune
       
# STEP 2: Limit test dataset to ld independent markers
echo "[M3,S2] Limit test dataset to ld independent markers"
$PLINK --bfile ../data/met583-test/met583-test-qc-v9 \
--extract ../temp/1000g-all-maf5-noWS-aims-overlap-noLD.prune.in \
--make-bed \
--out ../temp/met583-test-qc-v9-prune

# STEP 3: Merge test and reference datasets
echo "[M3,S3] Merge test and reference datasets"
$PLINK --bfile ../temp/all-1000g-phase3-chrall-mac5-v2-prune \
--flip ../temp/combined-merge.missnp \
--make-bed \
--out ../temp/all-1000g-phase3-chrall-mac5-v2-prune-flip

$PLINK --bfile ../temp/all-1000g-phase3-chrall-mac5-v2-prune-flip \
--bmerge ../temp/met583-test-qc-v9-prune \
--make-bed \
--out ../temp/combined
########## MODULE 3 END ##########


########## MODULE 4 START ##########
echo "<<< MODULE 4 >>>"

# STEP 1: Perform PCA analysis on the combined dataset
echo "[M4,S1] Perform PCA computing top 10 principle components"
$PLINK --bfile ../temp/combined \
--pca 10 header tabs \
--out ../results/combined_alt
########## MODULE 4 END ##########


########## MODULE 5 START ##########
echo "<<< MODULE 5 >>>"

echo "Modules 1-4 complete. Please run R script 'vis.R' for modules 5 and 6"
########## MODULE 5 END ############


########## MODULE 6 START ##########
$PLINK --bfile ../data/met583-test-qc-v9 --check-sex
awk 'BEGIN{printf("%10s%10s\n", "FID", "IID")};{ if ($4 == 2) {printf("%10s%10s\n", $1, $2)}}' plink.sexcheck > met583-test-female.keep
########## MODULE 6 END ##########


########## HOUSEKEEPING ##########
echo "Removing intermediate files..."
rm ../temp/1000g-gbr-maf5.*
rm ../temp/combined.*
echo "Done"
########## HOUSEKEEPING ##########