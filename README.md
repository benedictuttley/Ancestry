# Ancestry Inference
Pipeline to address the following question:
"Using the test dataset provided, identify a subset of female samples who have similar ancestry to the GBR reference genotypes in the 1000-genomes project."

## Input
### 1000 Genomes Reference Dataset
* all-1000g* - reference file for GBR ancestry
* Contains  individuals
* 'Within Family ID (IID)' matches corresponding FID uniformly
* 25,473,793 markers

### Test Dataset
* Unknown genotype dataset with origins in the USA
* Contains 3563 individuals
* 'Within Family ID (IID)' is 1 uniformly
* 575,867 markers

## Output

### Plots
* scree.png
* PC1_vs_PC2.png
* PC2_vs_PC3.png
* PC1_vs_PC3.png
### keep
* GBR-like.keep
* GBR-like-females.keep

## Preliminary Data Harmonisation
* Marker names have been updated to align with reference genome hg19
* indels updated to have I/D coding
* non-standard variation codes have been removed (CT/GA  -> T/A)
* Duplicate SNPs have been removed

## Running the Pipeline
* Six modules
* Type `bash ancestry_pipeline.sh` to run modules 1-4
* Run vis.R R script to run modules 5-6