
#!/bin/bash

# to run interactively: ./run.sh
# to submit: qsub run.sh

# Modify according to need
#$ -S /bin/bash
#$ -M lailhh@gis.a-star.edu.sg
#$ -V
#$ -cwd

# Set number of CPU here
#$ -pe OpenMP 4

#$ -l mem_free=10G,h_rt=12:00:00


### Changes accordingly
# master file
mf=$1

# working directory
wdir=$2

# number of mutations to be input to PyClone
num=$3

# whether to select mutations from copy neutral regions only
copy_neutral=$4 # "yes" or "no"

# random seed for selecting mutations
random_seed=$5

### script directory
script_dir=/mnt/projects/lailhh/workspace/pipelines/PyClone/PLANET_Pyclone_pipeline/scripts/S01_prepareInputs

### activate conda environment
conda activate my_R-3.6

### since phyloWGS assume all sequenza results are in the same folder
# --> create soft links for sequenza and 
# get a list of copy-neutral segments for each DNA lib
Rscript $script_dir/S01a_sequenza_CopyNeuRegs_perSam.r $mf $wdir


### for each patient,
# intersect copy-neutral segments from their samples 
# --> those interesected regions will be used as bed file to filter mutations in mutect
Rscript $script_dir/S01b_sequenza_CopyNeuRegs_perPat.r $mf $wdir


### for each sample, get bed file of coding regions using its MAF file (funcotator results)
Rscript $script_dir/S02a_mutect_CodingRegs_perSam.r $mf $wdir


### for each patient, combine coding regions from their samples 
Rscript $script_dir/S02b_mutect_CodingRegs_perPat.r $mf $wdir


### remove mutations with invalid chromosomes (i.e. those containing "_" in chromosome names) in mpileup VCF,
# sort the VCF file, then
# compress the sorted file, then
# index the compressed file
Rscript $script_dir/S03a_compress_index_mpileup_VCF.r $mf $wdir


### from mpileup VCF and coding bed file (from mutect MAF)
# use vcftools to separate mutations into 2 groups:
# - mutations in coding regions 
# - mutations in non-oding regions 
Rscript $script_dir/S03b_mpileup_CodingRegs.r $mf $wdir


### for each DNA lib, 
# using copy-neutral bed files (from sequenza results) and 
# mpileup of coding mutations and 
# mpileup of non-coding mutations
# --> output 2 files: 
# mutations in copy-neutral + coding regions
# mutations in copy-neutral + non-coding regions
Rscript $script_dir/S03c_mpileup_CodingRegs_CopyNeutral.r $mf $wdir


### for each patients, randomly select mutations from mpileup VCF
Rscript $script_dir/S03d_randomly_select_mutations.r $mf $wdir $num $copy_neutral $random_seed
