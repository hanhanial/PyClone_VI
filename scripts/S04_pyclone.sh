#!/bin/bash
#$ -V
#$ -pe OpenMP 4
#$ -l mem_free=4g,h_rt=120:00:00
#$ -cwd

# change accordingly
wdir=/mnt/projects/lailhh/workspace/pipelines/PyClone/testing/d20210429/
vcfList=$wdir/S03_list_of_vcfs.csv
script_dir=/mnt/projects/lailhh/workspace/pipelines/PyClone/PyClone_VI/scripts/helper_scripts
sequenzadir=/mnt/projects/zhaiww1/planet/data_storage/DNA_data/hg38/sequenza_results/cnvkit
mutectdir=/mnt/projects/zhaiww1/planet/data_storage/DNA_data/hg38/mutect2_results


# create output directories
odir=$wdir/S04_pyclone
mkdir -p $odir

cd $odir

/mnt/projects/lailhh/workspace/softwares/nextflow/nextflow run /mnt/projects/lailhh/workspace/pipelines/PyClone/PyClone_VI/scripts/S04_pyclone.nf -resume -config $script_dir/nextflow.config --wdir $odir --script_dir $script_dir --sequenzadir $sequenzadir --mutectdir $mutectdir --vcfList $vcfList
