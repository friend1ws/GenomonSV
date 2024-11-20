#! /usr/bin/env python

from .run import *
import argparse

def create_parser():

    ####################
    # top level parser
    parser = argparse.ArgumentParser(prog = "GenomonSV", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--version", action = "version", version = "GenomonSV-0.8.3")

    subparsers = parser.add_subparsers()

    ####################


    ####################
    # GenomonSV parse
    parse_parser = subparsers.add_parser("parse", 
                                         help = "Parse and cluster supporting read pairs for candidate structural variations",
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
     
    parse_parser.add_argument("bam_file", metavar = "input.bam", type = str,
                              help = "Path to input bam file")

    parse_parser.add_argument("output_prefix", metavar = "output_prefix", type = str,
                              help = "Prefix of the pathes of output files")

    parse_parser.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

    parse_junction_group = parse_parser.add_argument_group("parse_junction_condition",
                                                           "Parameters used for parsing breakpoint containing read pairs from bam files")

    parse_junction_group.add_argument("--junction_abnormal_insert_size", type = int, default = 1000,
                                      help = "Size of abnormal insert size. Used for checking the consistency of paired read of breakpoint containing reads (default: %(default)s)")

    parse_junction_group.add_argument("--junction_min_major_clipping_size", type = int, default = 15,
                                      help = "Minimum number of clipped bases for junction read (default: %(default)s)")

    parse_junction_group.add_argument("--junction_max_minor_clipping_size", type = int, default = 15,
                                      help = "If the clipped bases numbers of both sides are greater than this value, then the read is filtered out (default: %(default)s)")

    parse_junction_group.add_argument("--junction_min_mapping_qual", type = int, default = 0,
                                      help = "The minimum acceptable mapping qualitiy of the primary junction read (default: %(default)s)")
   
    parse_junction_group.add_argument("--reference_genome", type = str, default = "",
                                      help = "The path to the reference genomoe sequence")

    parse_junction_group.add_argument("--sort_option", metavar = "-S 1G", type = str, default = "-S 1G", 
                                      help = "Options for Linux sort command (default: '-S 1G')")

    cluster_junction_group = parse_parser.add_argument_group("cluster_junction_condition",
                                                             "Parameters used for clustering breakpoint containing read pairs")

    cluster_junction_group.add_argument("--junction_check_margin_size", type = int, default = 1000,
                                        help = "This value should be at least 1000. But very large value will lead to slow computation (default: %(default)s)")

    cluster_junction_group.add_argument("--junction_check_maximum_unique_pairs", type = int, default = 10000,
                                        help = "The number of unique junction pairs at each local region. \
                                                Set mainly for avoiding local regions with extremely high depths and huge amount of junction reads (default: %(default)s)")


    parse_improper_group = parse_parser.add_argument_group("parse_improper_condition",
                                                           "Parameters used for parsing improperly aligned read pairs from bam files")

    parse_improper_group.add_argument("--improper_abnormal_insert_size", type = int, default = 2000,
                                      help = "The size of abnormal insert size used for identifying improper read pairs (default: %(default)s)")

    parse_improper_group.add_argument("--improper_min_mapping_qual", type = int, default = 30,
                                      help = "The minimum acceptable mapping qualities of reads within improper read pairs (default: %(default)s)")

    parse_improper_group.add_argument("--improper_max_clipping_size", type = int, default = 5,
                                      help = "The maximum acceptable clipping size of reads within improper read pairs(default: %(default)s) ")


    cluster_improper_group = parse_parser.add_argument_group("cluster_improper_condition",
                                                             "Parameters used for parsing improperly aligned read pairs")

    cluster_improper_group.add_argument("--junction_dist_margin", type = int, default = 500,
                                        help = "Possible junction position margin from the junction direction (default: %(default)s)")

    cluster_improper_group.add_argument("--junction_opposite_dist_margin", type = int, default = 30,
                                        help = "Possible junction position in the opposite side (default: %(default)s)")

    cluster_improper_group.add_argument("--improper_check_margin_size", type = int, default = 1500,
                                        help = "This should be sufficiently greater than insert size, but the computational time will increase when too large (default: %(default)s)")

    cluster_improper_group.add_argument("--improper_check_maximum_unique_pairs", type = int, default = 10000,
                                        help = "The number of unique improper pairs at each local region. \
                                                Set mainly for avoiding local regions with extremely high depths and huge amount of junction reads (default: %(default)s)")

    parse_parser.set_defaults(func = genomonSV_parse)
    ####################   


    ####################
    # GenomonSV filt
    filt_parser = subparsers.add_parser("filt", 
                                         help = "Filter candidate structural variations")

    filt_parser.add_argument("bam_file", metavar = "input.bam", type = str,
                             help = "Path to input bam file")

    filt_parser.add_argument("output_prefix", metavar = "output_prefix", type = str,
                             help = "Prefix of the pathes of output files")

    filt_parser.add_argument("reference_genome", metavar = "reference.fa", type = str,
                             help = "The path to the reference genomoe sequence")

    # filt_parser.add_argument("annotation_dir", metavar = "annotation_dir", type = str,
    #                          help = "the path to database directory")

    filt_parser.add_argument("--matched_control_bam", metavar = "matched_control.bam", default = "", type = str,
                             help = "Path to matched control bam file")

    filt_parser.add_argument("--non_matched_control_junction", metavar = "merged.junction.control.bedpe.gz", default = "", type = str,
                             help = "Path to the non-matched control panel data generated by merge stage")

    filt_parser.add_argument("--matched_control_label", default = "", type = str,
                             help = "Matched control sample label (junction read pairs of this label is ignored in the non-matched control filtering step)")

    filt_parser.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg19",
                             help = "The genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")

    filt_parser.add_argument("--grc", default = False, action = 'store_true',
                             help = "Use Genome Reference Consortium nomenclature rather than UCSC (default: %(default)s)")

    filt_parser.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

    filt_parser.add_argument('--thread_num', default = 1, type=int,
                             help = "The number of threads")

    filt_parser.add_argument("--blat", default = False, action = 'store_true', help = "use blat in the realignment function")
    filt_parser.add_argument("--edlib", default = False, action = 'store_true', help = "use edlib in the realignment function")

    filt_parser.add_argument("--sort_option", metavar = "-S 1G", type = str, default = "-S 1G", 
                             help = "Options for Linux sort command (default: '-S 1G')")


    filter_condition_group = filt_parser.add_argument_group("filter_condition",
                                                             "Parameters used in various filtering steps in GenomonSV filt command")

    filter_condition_group.add_argument("--min_junc_num", type = int, default = 3,
                                       help = "Minimum required number of supporting junction read pairs (default: %(default)s)")

    filter_condition_group.add_argument("--min_sv_size", type = int, default = 10,
                                       help = "Minimum structural variation size (default: %(default)s)")

    filter_condition_group.add_argument("--min_inversion_size", type = int, default = 500,
                                        help = "Minimum inversion size (default: %(default)s)")

    filter_condition_group.add_argument("--control_panel_num_thres", type = int, default = 2,
                                        help = "Acceptable number of supporting junction read pairs in any of non-matched control samples (default: %(default)s)")

    filter_condition_group.add_argument("--control_panel_check_margin", type = int, default = 40,
                                        help = "Positions of junction slip very large (sometimes over 20bp). The author thinks this value should be at least 30bp, but may be more (default: %(default)s)")

    filter_condition_group.add_argument("--min_support_num", type = int, default = 3,
                                        help = "Mininum required number of supporting read pairs (both junction read pairs and improper read pairs) (default: %(default)s)")

    filter_condition_group.add_argument("--min_mapping_qual", type = int, default = 40,
                                        help = "Threshold of mapping quality (if the median of supporting reads is below this value, then filtered out) (default: %(default)s)")

    filter_condition_group.add_argument("--min_overhang_size", type = int, default = 50,
                                        help = "Minimum region size arround each break-point which have to be covered by at least one aligned short read (default: %(default)s)")

    filter_condition_group.add_argument("--close_check_margin", type = int, default = 25,
                                        help = "Check size for removing duplicated candidates (produced by alignment artifacts). The author believes 25 should be enough (default: %(default)s)")

    filter_condition_group.add_argument("--close_check_thres", type = int, default = 2,
                                        help = "If the #{breakpoint-containing read pairs} of the putative duplicated candidate is below this value, remove it (default: %(default)s)")


    realignment_condition_group = filt_parser.add_argument_group("realignment_condition",
                                                                 "parameters used in realignment filtering steps in GenomonSV filt command")

    realignment_condition_group.add_argument("--max_read_num", type = int, default = 5000,
                                             help = "when the total read number exceeds max_read_num, then max_read_num reads are randomly selected (default: %(default)s)")

    realignment_condition_group.add_argument("--search_length", type = int, default = 1000,
                                             help = "Read pairs within the specified length from the breakpoint will be checked (default: %(default)s)")

    realignment_condition_group.add_argument("--search_margin", type = int, default = 5,
                                             help = "Read pairs in the opposite side to the breakpoint within the specified length will be checked (default: %(default)s)")

    realignment_condition_group.add_argument("--split_refernece_thres", type = int, default = 1000,
                                             help = "Threshould for splitting the validation sequences (with variant) (default: %(default)s)")

    realignment_condition_group.add_argument("--validate_sequence_length", type = int, default = 1000,
                                            help = "Length of sequences for validation (each from the breakpoint) (default: %(default)s)")

    realignment_condition_group.add_argument("--short_tandem_reapeat_thres", type = int, default = 500,
                                             help = "Threshold for identifying short tandem duplication length (since STD has slightly different approach for realignment validation), \
                                                     this value should be significantly larger than the read length (default: %(default)s)")

    realignment_condition_group.add_argument("--blat_option", type = str, default = "-stepSize=5 -repMatch=2253",
                                             help = "Option used in blat")

    realignment_condition_group.add_argument("--min_tumor_variant_read_pair", type = int, default = 3,
                                             help = "Minimum required number of read pairs in tumor sample (default: %(default)s)")

    realignment_condition_group.add_argument("--min_tumor_allele_freq", type = float, default = 0.01,
                                             help = "Minimum required allele frequency of the tumor sample (default: %(default)s)")

    realignment_condition_group.add_argument("--max_control_variant_read_pair", type = int, default = 1,
                                             help = "Maximum allowed number of read pairs in matched control sample (default: %(default)s)")

    realignment_condition_group.add_argument("--max_control_allele_freq", type = float, default = 0.10,
                                             help = "Maximum allowed allele frequency of normal sample (default: %(default)s)")

    realignment_condition_group.add_argument("--max_fisher_pvalue", type = float, default = 0.10,
                                             help = "Maximum allowed fisher's exact test p-value (default: %(default)s)")



    filt_parser.set_defaults(func = genomonSV_filt)
    ####################


    ####################
    # GenomonSV merge 
    merge_parser = subparsers.add_parser("merge",
                                         help = "Merge clustered junction files for filtering")

    merge_parser.add_argument("control_info_file", metavar = "control_info.txt", type = str,
                              help = "Tab-delimited file on non-matched control. 1st column: sample labels, 2nd column: output prefix generated at the parse stage")

    merge_parser.add_argument("merge_output_file", type = str,
                              help = "The path to the merged junction file")

    merge_parser.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")

    merge_parser.add_argument("--merge_check_margin_size", type = int, default = 100,
                             help = "This value should be at least 50 (more than the length of possible inserted sequence between break points)")

    merge_parser.add_argument("--sort_option", metavar = "-S 1G", type = str, default = "-S 1G", 
                             help = "Options for Linux sort command (default: '-S 1G')")

    merge_parser.set_defaults(func = genomonSV_merge)
    ####################

    return parser

