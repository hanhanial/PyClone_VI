"""
Use this script to assign copy number information to the mutations generated
using the prepare_treeomics_pileup.py and treeomics workflow script for
multi-sectors data. Takes in the read depth and variant counts information
for all mutations called across different sectors of the same patient and
assign copy number information to the variants using the segments
information.

Segments file needs the following columns (Following Sequenza's naming
convention):
Tumor_Sample_Barcode chromosome start.pos end.pos CNt A B

Mutation TSV file (Or any tab delimited mutation file) uses the following
columns:
Tumor_Sample_Barcode Start_position Chromosome Hugo_Symbol

followed by samples from the 5th columns onwards.

If --singleSample is specified, will run the single sample version of the
script to prepare input for pyclone using single sample. In this case, -m is
followed by the mutation MAF file.

The script checks the segments where the mutations belong to and assign the
CN of that segments to the mutations. Note that X and Y chromosomes are
ignored. Indels are not considered as well. There are some edge cases where
the segments have zero major copy number which will lead to pyclone's error.
These segments together with the mutations are discarded. Note that pyclone
only looks at the intersect of the mutations between all sectors, which is
why the pileup is needed for mutations not called in any sector.

(The MAF in the script name has nothing to do with MAF file, just legacy
script naming...)


"""

import sys
import os
import glob
import pandas as pd
import getopt
import argparse
import re

def segments_to_dict(segments):
    """
    This convert the segments to a dictionary with the chromosomes as the
    keys and the segments as its value.
    """
    seg_dict = {}
    for chrom in segments.chromosome.unique():
        seg_dict[chrom] = segments[segments['chromosome'] == chrom]
    return(seg_dict)

def search_overlap_singleSample(mutation, seg_dict):
    """
    Search for the segment overlapping the mutation. If an indel or if cannot
    find overlapping segments, return empty dataframe. Indel less than 20bp
    is allowed.
    """
    chrom_to_search = mutation['Chromosome']
    mut_start = mutation['Start_position']
    mut_end = mutation['End_position']
    if (mut_end - mut_start > 25):
        print("{}:{}:{} is a long indel (25bp)".format(mutation['Hugo_Symbol'], mutation['Chromosome'], mutation['Start_position']))
        return(pd.DataFrame)
    # reduce to segments with start position earlier than the mutation position.
    true_start = seg_dict[chrom_to_search][seg_dict[chrom_to_search]['start.pos'] < mut_start]
    # Find the segment by checking if the end position of the segment is higher
    # than the mutation position
    overlap_segment = true_start[true_start['end.pos'] > mut_end]
    if (overlap_segment.empty):
        # print("No overlap!")
        # print(chrom_to_search, mut_start, mut_end)
        return(pd.DataFrame)
    return(overlap_segment)

def search_overlap(mutation, seg_dict):
    """
    Search for the segment overlapping the mutation. This function does not
    handle indels, as it assumes indels are already removed in the input.
    E.g. this is used mostly for the multi-sectoring version of main function
    below since that removes indels prior to calculating overlap.
    """
    chrom_to_search = str(mutation['Chromosome'])
    try:
        mut_start = mutation['Position']
    except KeyError:
        mut_start = mutation['Start_position']
    except:
        print("Something wrong with the mutation start position. Check if the column is proper named as Position or Start_position")

    # reduce to segments with start position earlier than the mutation position.
    true_start = seg_dict[chrom_to_search][seg_dict[chrom_to_search]['start.pos'] < mut_start]
    # Find the segment by checking if the end position of the segment is higher
    # than the mutation position
    overlap_segment = true_start[true_start['end.pos'] > mut_start]
    if (overlap_segment.empty):
        # print("No overlap!")
        # print(chrom_to_search, mut_start, mut_end)
        return(pd.DataFrame)
    return(overlap_segment)

def main(var_file, rd_file, segment_file, patient_names, vaf_threshold=0.05, filterSegments = False):
    """
    Main function to assign CN using the functions above. Iterates through
    sample and output combined MAF file with CNV information.
    """
    patient = patient_names

    patient_varcount = pd.read_csv(var_file, low_memory=False, delimiter="\t")
    patient_readdepth = pd.read_csv(rd_file, low_memory=False, delimiter="\t")

    # Sanity check to see if the columns are identical
    unmatch = patient_varcount.loc[patient_varcount.loc[:, 'Chromosome'] != patient_readdepth.loc[:, 'Chromosome']]
    if (unmatch.empty != True):
        print("Something wrong with sample, columns order do not match!")

    tumor_sample = patient_varcount.columns[4:]
    info_col = patient_varcount.columns[:4]

    # Can use this to remove indels
    # patient_readdepth = patient_readdepth.loc[patient_readdepth['Change'].str.contains('-')!=True]
    # patient_varcount = patient_varcount.loc[patient_varcount['Change'].str.contains('-')!=True]

    # Make sure there's no zero read depth position for any sector, as Pyclone
    # will assume that the mutatation has identical VAF at that sector
    tmp = (patient_readdepth.loc[:, tumor_sample] == 0)
    patient_readdepth = patient_readdepth.loc[tmp.any(axis=1)==False]
    patient_varcount = patient_varcount.loc[tmp.any(axis=1)==False]

    # Transform RD to ref count which is just the difference between RD and varcount
    patient_readdepth.iloc[:, 4:] = patient_readdepth.iloc[:, 4:] - patient_varcount.iloc[:, 4:]
    # Get VAF and filter out those with < 0.05 VAF called in any sector.
    patient_VAF = patient_varcount.iloc[:, 4:] / patient_readdepth.iloc[:, 4:]
    patient_VAF = (patient_VAF < vaf_threshold)

    # Remove the mutations where the condition is true for ALL segments, i.e. it has to be below
    # 0.05 for all sectors. If it's above 0.05 in any sector, keep the mutations. This will keep most
    # of the private mutations.
    filter_VAF_index = (patient_VAF.all(axis=1) == False)
    num_filtered = filter_VAF_index.loc[filter_VAF_index == False, ]
    print("Patient {} has {} mutations with average VAF < {} removed".format(patient, num_filtered.shape[0], vaf_threshold))
    # Filter out the variants
    patient_readdepth = patient_readdepth.loc[filter_VAF_index, ]
    patient_varcount = patient_varcount.loc[filter_VAF_index, ]

    all_segments = pd.read_csv(segment_file, low_memory=False, delimiter='\t')

    if not os.path.exists("{}_mutations_withCN".format(patient)):
        os.makedirs("{}_mutations_withCN".format(patient))
    if not os.path.exists("{}_pyclone_input".format(patient)):
        os.makedirs("{}_pyclone_input".format(patient))

    for sample in tumor_sample:
        # The treeomics input has this weird problem of not accepting dash
        # in the name, so the output from my script in preparing treeomics
        # input has underscore instead. Change it back here.
        samplename = re.sub(r'_', r'-', sample)
        print(samplename)
        col_to_get = list(info_col)
        col_to_get.extend([sample])
        var_pat = patient_varcount.loc[:, col_to_get]
        var_pat.rename(columns={sample:"var_counts"}, inplace=True)
        ref_pat = patient_readdepth.loc[:, col_to_get]
        ref_pat.rename(columns={sample:"ref_counts"}, inplace=True)
        merge_sample_mut = var_pat.merge(ref_pat, how="left")
        merge_sample_mut.loc[:, 'normal_cn'] = 2
        merge_sample_mut.loc[:, 'mutation_id'] = merge_sample_mut.loc[:, 'Gene'].map(str) + "_" + merge_sample_mut.loc[:, 'Chromosome'].map(str) + ":" + merge_sample_mut.loc[:, 'Position'].map(str)
        sample_segments = all_segments[all_segments['Tumor_Sample_Barcode'] == samplename]

        seg_dict = segments_to_dict(sample_segments)

        overlap_seg = pd.DataFrame()
        filtered_seg = pd.DataFrame()
        for _, mut_row in merge_sample_mut.iterrows():
            # Skip X and Y chromosome
            if (mut_row['Chromosome'] == "X" or mut_row['Chromosome'] == "Y"):
                continue

            # Search for the segment
            buf = search_overlap(mut_row, seg_dict)
            # Skip if no overlapping segments
            if (buf.empty):
                continue
            # Filter segments with unreliable calls. This is according to Canopy's guideline. However, I set CNt to 8 instead of 6 since
            # LUAD tends to have higher ploidy than the other cancer types.
            elif filterSegments:
                print("--filterSegments specified. Will filter segments of low quality.")
                if (buf.iloc[0]['numMarker'] < 100) or (buf.iloc[0]['end.pos'] - buf.iloc[0]['start.pos'] < 5000000) or (buf.iloc[0]['CNt'] >= 8):
                    if (filtered_seg.empty):
                        filtered_seg = buf.iloc[0].to_frame()
                    else:
                        filtered_seg = pd.concat([filtered_seg, buf.iloc[0]], axis=1)
            else:
                # Get copy number for mutations
                assigned_row = mut_row.copy(deep=True)
                assigned_row['CNt'] = buf.iloc[0]['CNt']
                assigned_row['major_cn'] = buf.iloc[0]['A']
                assigned_row['minor_cn'] = buf.iloc[0]['B']
                # Initialize dataframe for merging.
                if (overlap_seg.empty):
                    overlap_seg = assigned_row.to_frame()
                else:
                    overlap_seg = pd.concat([overlap_seg, assigned_row], axis=1)

        overlap_seg = overlap_seg.transpose()
        overlap_seg.to_csv("./{}_mutations_withCN/{}_SNV_withCN.maf".format(patient, samplename),sep="\t", index=False)

        filtered_seg = filtered_seg.transpose()
        print("Sample {} has {} segments with marker<100 or smaller than 5 Mb or >= 8 copy number (Canopy guideline)".format(sample, filtered_seg.shape[0]))
        filtered_seg.to_csv("./{}_mutations_withCN/{}_filtered_seg.maf".format(patient, samplename),sep="\t", index=False)

        towrite = overlap_seg.loc[:, ['mutation_id', 'ref_counts', 'var_counts', 'normal_cn', 'minor_cn', 'major_cn']]
        # Remove those with major CN = 0. Most likely false positive. Note this will however remove the mutations across all
        # sectors when Pyclone run analysis
        weird_mut = towrite.loc[towrite.loc[:, 'major_cn'] == 0]
        print("{} mutations for sample {} are located in regions with major_cn 0!".format(weird_mut.shape[0], samplename))
        towrite = towrite.loc[towrite.loc[:, 'major_cn'] != 0]
        towrite['ref_counts'] = towrite['ref_counts'].map(int)
        towrite['var_counts'] = towrite['var_counts'].map(int)

        towrite.to_csv("./{}_pyclone_input/{}.tsv".format(patient, samplename), sep='\t', index=False)

def main_SS(maf_file, segment_file, vaf_threshold = 1.05, filterSegments = False):
    """
    Main function to assign CN using the functions above. Iterates through
    sample and output combined MAF file with CNV information.
    """
    all_mutations = pd.read_csv(maf_file, low_memory=False, delimiter='\t')
    all_segments = pd.read_csv(segment_file, low_memory=False, delimiter='\t')

    if not os.path.exists("./sample_mutations_withCN"):
        os.makedirs("./sample_mutations_withCN")
    if not os.path.exists("./pyclone_input"):
        os.makedirs("./pyclone_input")

    for i, sample in enumerate(all_mutations.Tumor_Sample_Barcode.unique()):
        print("Processing sample {}: {}".format(i+1, sample))

        # Subset the mutations and segments to those belonging to the patient
        sample_mutations = all_mutations[all_mutations['Tumor_Sample_Barcode'] == sample]
        sample_segments = all_segments[all_segments['Tumor_Sample_Barcode'] == sample]

        patient_VAF = sample_mutations.loc[:, 'VAF']
        filter_VAF_index = (patient_VAF > vaf_threshold)

        # Remove the mutations where the condition is true for ALL segments, i.e. it has to be below
        # 0.05 for all sectors. If it's above 0.05 in any sector, keep the mutations. This will keep most
        # of the private mutations.
        num_filtered = filter_VAF_index.loc[filter_VAF_index == False, ]
        print("Patient {} has {} mutations with average VAF < {} removed".format(sample, num_filtered.shape[0], vaf_threshold))
        # Filter out the variants
        sample_mutations = sample_mutations.loc[filter_VAF_index, ]
        # Get the segments dictionary for the patient.
        seg_dict = segments_to_dict(sample_segments)

        overlap_seg = pd.DataFrame()
        filtered_seg = pd.DataFrame()
        for _, mut_row in sample_mutations.iterrows():
            # Skip X and Y chromosome
            if (mut_row['Chromosome'] == "X" or mut_row['Chromosome'] == "Y"):
                continue

            # Search for the segment
            buf = search_overlap_singleSample(mut_row, seg_dict)
            # Skip if no overlapping segments
            if (buf.empty):
                continue
            elif filterSegments:
                print("--filterSegments specified. Will filter segments of low quality.")
                if (buf.iloc[0]['numMarker'] < 100) or (buf.iloc[0]['end.pos'] - buf.iloc[0]['start.pos'] < 5000000) or (buf.iloc[0]['CNt'] >= 8):
                    if (filtered_seg.empty):
                        filtered_seg = buf.iloc[0].to_frame()
                    else:
                        filtered_seg = pd.concat([filtered_seg, buf.iloc[0]], axis=1)
            else:
                # Get copy number for mutations
                assigned_row = mut_row.copy(deep=True)
                assigned_row['CNt'] = buf.iloc[0]['CNt']
                assigned_row['Major_CN'] = buf.iloc[0]['A']
                assigned_row['Minor_CN'] = buf.iloc[0]['B']
                assigned_row['adjustedCN'] = buf.iloc[0]['adjustedCN']
                # Initialize dataframe for merging.
                if (overlap_seg.empty):
                    overlap_seg = assigned_row.to_frame()
                else:
                    overlap_seg = pd.concat([overlap_seg, assigned_row], axis=1)

        overlap_seg = overlap_seg.transpose()
        overlap_seg.to_csv("./sample_mutations_withCN/{}_SNV_withCN.maf".format(sample),sep="\t", index=False)

        filtered_seg = filtered_seg.transpose()
        print("Sample {} has {} segments with marker<100 or smaller than 5 Mb or >= 8 copy number (Canopy guideline)".format(sample, filtered_seg.shape[0]))
        filtered_seg.to_csv("./sample_mutations_withCN/{}_filtered_seg.maf".format(sample),sep="\t", index=False)

        pyclone_input = overlap_seg.loc[:, ['Hugo_Symbol', 'Chromosome',
        'Start_position', 'ref_count', 'alt_count', 'VAF', 'Major_CN',
        'Minor_CN']]
        pyclone_input['mutation_id'] = pyclone_input['Hugo_Symbol'].map(str) + "_" + pyclone_input['Chromosome'].map(str) + ":" + pyclone_input['Start_position'].map(str)
        pyclone_input['normal_cn'] = 2
        towrite = pyclone_input.loc[:, ['mutation_id', 'ref_count', 'alt_count', 'normal_cn', 'Minor_CN', 'Major_CN']]
        towrite.columns = ['mutation_id', 'ref_counts', 'var_counts', 'normal_cn', 'minor_cn', 'major_cn']
        towrite['ref_counts'] = towrite['ref_counts'].map(int)
        towrite['var_counts'] = towrite['var_counts'].map(int)
        towrite.to_csv("./pyclone_input/{}_mutations.tsv".format(sample), sep='\t', index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-m", "--mutation_tsv",
                        help="mutant reads tsv file")
    parser.add_argument("-d", "--depth_tsv",
                        help="reads depth tsv file")
    parser.add_argument("-s", "--segments_file",
                        help="CNV file (e.g. merged sequenza's segments file) for all sectors")
    parser.add_argument("-p", "--patient_prefix",
                        help="Prefix for sample, e.g. patient name A511.")
    parser.add_argument("-f", "--vaf_threshold", type=float, default=0.05,
                        help="VAF threshold to filter variants (Default = 0.05)")
    parser.add_argument("-u", "--singleSample", action="store_true",
                        help="if specified, run the single sample version to prepare pyclone input for individual sample")
    parser.add_argument("--filterSegments", action="store_true",
                        help="if specified, filter the segments using Canopy suggested criteria.")

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    if (args.singleSample):
        print("Running single sample version of the script.")
        main_SS(
            args.mutation_tsv,
            args.segments_file,
            vaf_threshold=args.vaf_threshold,
            filterSegments = args.filterSegments)
    else:
        main(
            args.mutation_tsv,
            args.depth_tsv,
            args.segments_file,
            args.patient_prefix,
            vaf_threshold=args.vaf_threshold,
            filterSegments = args.filterSegments)
