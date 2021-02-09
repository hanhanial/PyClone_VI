#!/usr/bin/env nextflow

params.script_dir = "/mnt/projects/lailhh/workspace/pipelines/PyClone/PLANET_Pyclone_pipeline/scripts/S02_PyClone/scripts"

// Run parameter, change accordingly
params.merged_maf = "$baseDir/test.maf"
params.purity = "$baseDir/test.purity"
params.segments = "$baseDir/test.seg"
params.citup_type = "iter"
params.genome_dir = "/mnt/projects/lailhh/workspace/pipelines/PyClone/PLANET_Pyclone_pipeline/scripts/S02_PyClone/reference/human_g1k_v37.fasta"

// [hechuan] params
// mpileup
params.mpileupdir = "/mnt/projects/lailhh/workspace/Liver_TCR/PyClone/PyClone_package/mpileup_results"

// sequenza
params.sequenzadir = "/mnt/projects/zhaiww1/planet/data_storage/DNA_data/hg38/sequenza_results/cnvkit"

// mutect
params.mutectdir = "/mnt/projects/zhaiww1/planet/data_storage/DNA_data/hg38/mutect2_results"

// path to your Rscript
params.R_path = "/mnt/projects/chuakp/khipin_lungtcr/softwares/anaconda3/envs/r34/bin/Rscript"


println """\
	
	########################################################################

	Use a post-processed merged MAF file for all samples together
	with merged segments output (Any program as long as with major
	and minor copy number) and purity/ploidy estimates.

	Merged MAF file: ${params.merged_maf}
	Merged purity file: ${params.purity}
	Merged segment files: ${params.segments}

	Important Note (These will fix most errors):
	Script assumed the sample to have 4 characters prefix, i.e. A511_xxx
	The headers required below are case-sensitive.

	MAF file requires extra column of VAF, alt_freq and ref_freq

	Segments files header should be:
	Tumor_Sample_Barcode, chromosome, start.pos, end.pos, CNt, A, B

	Purity file requires Tumor_Sample_Barcode, purity and ploidy 
	(Order is important).

	For get_pileup_pos step, make sure the bams directory for each
	sample can be read correctly. Modify accordingly. Here the script
	assumes the bams can be found in /path/to/mutectdir/samplename/*.bam

	Note that for indels, the sectors that do not have the calls will
	have read depth set to the average read depth of other sectors.

	Specify NF parameter citup_type to "qip" if you want to use the QIP
	version of Citup. In this mode, the bash script will assign a different
	"cluster id" to each of the cluster fed into citup.

	Note that the QIP version requires you to set the cluster assignment to
	be the maxnode so for the cluster we get from PyClone, we can't just set
	dummy assignment as per what I am doing here. E.g. if max node is 10, 
	each cluster needs to be assigned from 0 to 9.

	Iter version (Default) will run clustering prior to phylogeny inference. So 
	subclones with similar CCF *can* be clustered together.
	Normally only an issue if you have a lot of subclones.

	By default the binomial model of PyClone is used as I assume WES/WGS data
	which has low coverage. Change "density" parameter in pyclone_pipeline.sh
	script if you want to use beta binomial model (deep sequencing).

	Usual workflow:
	1. Copy qsub_pyclone.sh, nextflow.config, run_pyclone.nf to the working directory.
	2. Change qsub_pyclone.sh to point to the correct segments, purity, and MAF file.
	3. Change mutectdir in run_pyclone.nf script to point to correct BAM locations.
	4. (Optional) Add "errorStrategy='ignore'" if you want the pipeline to run in
	   the case of failures on any patient. Otherwise the pipeline will stop
	   whenever there's a failure.
	5. Do "qsub -V qsub_pyclone.sh" to run the pipeline. Check the stdout file
	   to monitor progress and error.

	########################################################################

"""

// [Hechuan] start from mpileup output

mpileup_vcfs=Channel.fromPath(params.mpileupdir+'/*.vcf')
mpileup_vcfs
    .map {vcf-> tuple(vcf.baseName.substring(0,6),vcf)}
    .set {pat_vcfs}

process pileup_pos {
    maxForks 25
    publishDir 'pileup_results', mode: 'copy'

    input: 
    set pat_id, file(vcf) from pat_vcfs

    output:
	set val(pat_id), file("${pat_id}_mutant_reads.tsv"), file("${pat_id}_read_depth.tsv") into pat_muts
	set val(pat_id), file("${pat_id}.sectors") into pat_sectors

	script:
	"""
    python ${params.script_dir}/var_total_reads_from_vcf.py \
        -v $vcf -p $pat_id -m ${params.mutectdir}
    """
}

process prepare_seg_purity {
    maxForks 25
    publishDir 'seg_purity', mode: 'copy'

    input:
    set pat_id, file(sectors) from pat_sectors

    output:
    set val(pat_id), file("${pat_id}.seg") into pat_segs
    set val(pat_id), file("${pat_id}.purity") into pat_purity

    """
    python ${params.script_dir}/seg_purity_file_from_sectors.py \
        -s $sectors -p $pat_id -d ${params.sequenzadir}
    """
} 

pat_muts_segs = pat_muts.join(pat_segs)

// Assign major and minor copy number information to mutatons to prepare for pyclone

process assign_CN {
    maxForks 25
    publishDir 'pyclone_input', mode: 'copy'

    input:
    set pat_id, file(var_reads_file), file(depth_file), file(segs_file) from pat_muts_segs

    output:
    set val(pat_id), file("${pat_id}_pyclone_input/*.tsv") into pyclone_input_tmp
	
    script:
    """
    python ${params.script_dir}/assignCN_to_maf_multiSec.py -m ${var_reads_file} -s ${segs_file} -p ${pat_id} -d ${depth_file}
    ${params.script_dir}/remove_muts_with_no_CN.sh $pat_id
    """
}






pyclone_input = pyclone_input_tmp.join(pat_purity)

process run_pyclone_pipeline {
    maxForks 25
    conda '/home/lailhh/.conda/envs/PyClone'
    clusterOptions='-V -pe OpenMP 1 -l mem_free=20g,h_rt=96:00:00'
    publishDir 'pyclone_results', mode: 'copy'
	
    input:
    set pat_id, file(sample_input), file(purity) from pyclone_input

    output:
    set val(pat_id), file("output/${pat_id}/*") into pyclone_results_ch

    script:
    """
    ${params.script_dir}/pyclone_pipeline.sh ${pat_id} ${purity}
    """
}



// Run citup on the clusters' CCF. The citup script set a maximum of 10 
// clones as it will otherwise be extremely slow

process run_citup {
    maxForks 25
    conda '/home/lailhh/.conda/envs/PyClone'
    clusterOptions='-V -pe OpenMP 24 -l mem_free=20g,h_rt=48:00:00'
    errorStrategy 'ignore'
    publishDir "citup_results/${pat_id}", mode: 'copy'

    input:
    set pat_id, file(pat_folder) from pyclone_results_ch

    output:
    set val(pat_id), file("${pat_id}_*.*") into citup_results_ch

    script:
    if ( params.citup_type == 'iter')
	"""
	${params.R_path} ${params.script_dir}/calculate_ITH_clonal_subclonal.R ${pat_id}
	${params.script_dir}/run_citup.sh ${pat_id}
		
	# Get tree
	python ${params.script_dir}/get_optimal_tree.py ${pat_id}_results.h5 ${pat_id} iter
	"""
    else
	"""
	${params.R_path} ${params.script_dir}/calculate_ITH_clonal_subclonal.R ${pat_id}
	${params.script_dir}/run_citup_qip.sh ${pat_id}
		
	# Get tree
	python ${params.script_dir}/get_optimal_tree.py ${pat_id}_results.h5 ${pat_id} QIP
	"""
}



process draw_mapscape {
    maxForks 25
    clusterOptions='-V -pe OpenMP 1 -l mem_free=2g,h_rt=1:00:00'
    errorStrategy 'ignore'
    publishDir "mapscape_results", mode: 'copy'

    input:
    set pat_id, file(citup_tree) from citup_results_ch

    output:
    set val(pat_id), file("*_mapscape.html") into mapscape_ch

    script:
    if ( params.citup_type == 'iter')
	"""
	${params.R_path} ${params.script_dir}/plot_tree.R ${pat_id} iter
	"""
    else
	"""
	${params.R_path} ${params.script_dir}/plot_tree.R ${pat_id} qip
	"""
}


// Touch when done

process done {
    publishDir 'pipeline_complete_indicator'

    input:
    set pat_id, file(finished_files) from mapscape_ch

    output:
    file("${pat_id}.done")

    script:
    """
    echo -e "${pat_id}\t${finished_files}" > ${pat_id}.done
    """
}
