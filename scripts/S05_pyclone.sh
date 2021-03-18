#!/bin/bash
#$ -V
#$ -pe OpenMP 4
#$ -l mem_free=16g,h_rt=24:00:00
#$ -cwd

# change accordingly
wdir=/mnt/projects/lailhh/workspace/pipelines/PyClone/testing/d20210318/
vcfList=$wdir/S04_list_of_vcfs.csv
script_dir=/mnt/projects/lailhh/workspace/pipelines/PyClone/PyClone_VI/scripts/helper_scripts
sequenzadir=/mnt/projects/zhaiww1/planet/data_storage/DNA_data/hg38/sequenza_results/cnvkit
mutectdir=/mnt/projects/zhaiww1/planet/data_storage/DNA_data/hg38/mutect2_results


# create output directories
odir=$wdir/S05_pyclone
mkdir -p $odir

cd $odir

/mnt/projects/lailhh/workspace/softwares/nextflow/nextflow run /mnt/projects/lailhh/workspace/pipelines/PyClone/PyClone_VI/scripts/S05_pyclone.nf \ 
  -resume -C $script_dir/nextflow.config \
  --wdir $odir \
  --script_dir $script_dir --sequenzadir $sequenzadir \
  --mutectdir $mutectdir --vcfList $vcfList
