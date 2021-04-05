// working directory
params.wdir = "/mnt/projects/lailhh/workspace/pipelines/PyClone/testing/d20210318/S05_pyclone"

// directory to helper scripts
params.script_dir = "/mnt/projects/lailhh/workspace/pipelines/PyClone/PyClone_VI/scripts/helper_scripts"

// number of random muts to be selected
params.num_mut = 3000

// seed for randomly selecting muts
params.random_seed = 1234

// sequenza
params.sequenzadir = "/mnt/projects/zhaiww1/planet/data_storage/DNA_data/hg38/sequenza_results/cnvkit"

// mutect
params.mutectdir = "/mnt/projects/zhaiww1/planet/data_storage/DNA_data/hg38/mutect2_results"

// list of VCF files
params.vcfList = "/mnt/projects/lailhh/workspace/pipelines/PyClone/testing/d20210318/S03_list_of_vcfs.csv"

// version of citup
params.citup_type = "iter"

// Rscript
params.R_path = "/mnt/projects/chuakp/khipin_lungtcr/softwares/anaconda3/envs/r34/bin/Rscript"

Channel
  .fromPath(params.vcfList)
  .splitCsv(header:true)
  .map{ row -> tuple(row.Unified_ID, row.DNA_lib, row.coding_VCF, row.noncoding_VCF) }
  .set { VCFs_ch }


process pileup_pos {
    maxForks 25
    publishDir 'pileup_results', mode: 'copy'

    input:
    set Unified_ID, DNA_lib, coding_VCF, noncoding_VCF from VCFs_ch

    output:
    set Unified_ID, file("${Unified_ID}_coding_mutant_reads.tsv"), file("${Unified_ID}_coding_read_depth.tsv"), file("${Unified_ID}_noncoding_mutant_reads.tsv"), file("${Unified_ID}_noncoding_read_depth.tsv") into pat_muts
    set Unified_ID, file("${Unified_ID}.sectors") into pat_sectors

    script:
    """
    python ${params.script_dir}/var_total_reads_from_vcf.py \
      -v $coding_VCF -p $Unified_ID -m ${params.mutectdir} -t coding
    python ${params.script_dir}/var_total_reads_from_vcf.py \
      -v $noncoding_VCF -p $Unified_ID -m ${params.mutectdir} -t noncoding
    mv ${Unified_ID}_noncoding.sectors ${Unified_ID}.sectors
    """
} 


process prepare_seg_purity {
    maxForks 25
    publishDir 'seg_purity', mode: 'copy'

    input:
    set Unified_ID, file(sectors) from pat_sectors
    
    output:
    set Unified_ID, file("${Unified_ID}.seg") into pat_segs
    set Unified_ID, file("${Unified_ID}.purity") into pat_purity

    """
    python ${params.script_dir}/seg_purity_file_from_sectors.py \
        -s $sectors -p $Unified_ID -d ${params.sequenzadir}
    """
}


pat_muts_segs = pat_muts.join(pat_segs)

// Assign major and minor copy number information to mutatons to prepare for pyclone

process assign_CN {
    maxForks 25
    publishDir 'assign_CN', mode: 'copy'

    input:
    set Unified_ID, file(coding_var_reads_file), file(coding_depth_file), file(noncoding_var_reads_file), file(noncoding_depth_file), file(segs_file) from pat_muts_segs

    output:
    set Unified_ID, file("${Unified_ID}_pyclone_input/*.tsv") into assign_CN_ch

    shell:
    """
    ${params.script_dir}/assignCN_to_maf_multiSec.sh $Unified_ID $segs_file $coding_var_reads_file $coding_depth_file $noncoding_var_reads_file $noncoding_depth_file ${params.script_dir}
    
    ${params.script_dir}/remove_muts_with_no_CN.sh $Unified_ID
    """
}



process pyclone_input {
  maxForks 25
  publishDir 'pyclone_input', mode: 'copy'
  
  input:
  set Unified_ID, file(pat_sams) from assign_CN_ch
  
  output:
  set Unified_ID, file("${Unified_ID}/*.tsv") into pyclone_input_ch

  script:
  """
  ${params.R_path} ${params.script_dir}/randomly_select_mutations.r ${Unified_ID} ${params.num_mut} ${params.random_seed} ${params.wdir}  
  """
}



process run_pyclone {
  maxForks 25
  conda '/home/lailhh/miniconda3/envs/pyclone-vi'
  clusterOptions='-V -pe OpenMP 1 -l mem_free=20g,h_rt=96:00:00'
  publishDir 'pyclone_results', mode: 'copy'

  input:
  set Unified_ID, file(pat_input) from pyclone_input_ch

  output:
  set Unified_ID, file("${Unified_ID}/*") into pyclone_results_ch

  script:
  """
  ${params.script_dir}/pyclone_pipeline.sh ${Unified_ID} ${Unified_ID}.tsv
  """  
}



/*
process run_citup {
    maxForks 25
    conda '/home/lailhh/miniconda3/envs/citup'
    clusterOptions='-V -pe OpenMP 24 -l mem_free=20g,h_rt=48:00:00'
    errorStrategy 'ignore'
    publishDir "citup_results/${pat_id}", mode: 'copy'

    input:
    set Unified_ID, file(pat_folder) from pyclone_results_ch

    output:
    set Unified_ID, file("${Unified_ID}_*.*") into citup_results_ch

    script:
    if ( params.citup_type == 'iter')
        """
        ${params.R_path} ${params.script_dir}/preparing_citup_inputs.r ${Unified_ID}
        ${params.script_dir}/run_citup.sh ${Unified_ID}

        # Get tree
        python ${params.script_dir}/get_optimal_tree.py ${Unified_ID}_results.h5 ${Unified_ID} iter
        """
    else
        """
        ${params.R_path} ${params.script_dir}/preparing_citup_inputs.R ${Unified_ID}
        ${params.script_dir}/run_citup_qip.sh ${Unified_ID}

        # Get tree
        python ${params.script_dir}/get_optimal_tree.py ${Unified_ID}_results.h5 ${Unified_ID} QIP
        """
}



process done {
    publishDir 'pipeline_complete_indicator', mode: 'copy'

    input:
    set Unified_ID, file(pat_files) from citup_results_ch 

    output:
    file("${Unified_ID}.done")

    script:
    """
    touch ${Unified_ID}.done
    """
}


*/

