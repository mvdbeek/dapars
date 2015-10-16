import argparse
import numpy as np
import os
import sys
import datetime

import scipy as sp

from bisect import bisect

from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector

import math


def time_now():  # return time
    curr_time = datetime.datetime.now()
    return curr_time.strftime("%c")


def Convert_wig_into_bp_coverage(extracted_coverage, extracted_3UTR_region, strand_info):
    bp_coverage = np.zeros(extracted_3UTR_region[-1] - extracted_3UTR_region[0])
    relative_start = extracted_3UTR_region[0]
    for i in range(len(extracted_coverage)):
        start = extracted_3UTR_region[i] - relative_start
        curr_region_end = extracted_3UTR_region[i + 1] - relative_start
        bp_coverage[start:curr_region_end] = extracted_coverage[i]
    if strand_info == '-':
        bp_coverage = bp_coverage[::-1]

    return bp_coverage


def create_output_dir(output_directory):
    d = os.path.dirname(output_directory)
    if not os.path.exists(d):
        os.makedirs(d)
    temp_dir = os.path.join(d, '/tmp/')
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    return temp_dir

def parse_args(argv):
    """
    Returns floating point values except for input files.
    My initial approach will not filter anything. (FDR. fold_change, PDUI, Num_least ...)
    :param argv:
    :return:
    """
    parser = argparse.ArgumentParser(description='Determines the usage of proximal polyA usage')
    parser.add_argument("-c", "--control_alignments", nargs="+", required=True,
                        help="Alignment files in BAM format from control condition")
    parser.add_argument("-t", "--treatment_alignments", nargs="+", required=True,
                        help="Alignment files in BAM format from treatment condition")
    parser.add_argument("-u", "--utr_bed_file", nargs=1, required=True, type=file,
                        help="Bed file describing longest 3UTR positions")
    parser.add_argument("-o", "--output_file", nargs=1, required=True, type=argparse.FileType('w'),
                        help="file containing output")
    return parser.parse_args(argv)

def process_alignments(All_Wig_files, utr_dict):
    """
    Takes in an iterable and produces a 3\'UTR coverage dict and a numpy array of total sample depth.
    :param All_Wig_files:
    :param utr_dict:
    :return:
    """
    """
    :param All_Wig_files:
    :param utr_dict:
    :return:
    """
    def make_coverage_dict(wig_file):
        """Get coverage_dict for wig_file"""
        coverage_dict = {}
        depth = 0
        for line in wig_file:
            if '#' not in line and line[0:3] == 'chr':
                fields = line.strip('\n').split('\t')
                chrom_name = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                count = int(fields[-1])
                depth += int(fields[-1]) * (end - start)
                if chrom_name not in coverage_dict:
                    coverage_dict[chrom_name] = [[0], [0]]
                if start > coverage_dict[chrom_name][0][-1]:
                    coverage_dict[chrom_name][0].append(start)
                    coverage_dict[chrom_name][1].append(0)
                coverage_dict[chrom_name][0].append(end)
                coverage_dict[chrom_name][1].append(count)
        coverage_dict[chrom_name][1].append(0)
        return coverage_dict, depth

    def get_utr_coverage(coverage_dict, utr_dict):
        for curr_3UTR_event_id in utr_dict:
            curr_3UTR_structure = utr_dict[curr_3UTR_event_id]
            curr_chr = curr_3UTR_structure[0]
            if curr_chr in coverage_dict:
                curr_chr_coverage = coverage_dict[curr_chr]
                start = curr_3UTR_structure[1]
                end = curr_3UTR_structure[2]
                left_region_index = bisect(curr_chr_coverage[0], start)
                right_region_index = bisect(curr_chr_coverage[0], end)

                extracted_coverage = curr_chr_coverage[1][left_region_index:right_region_index + 1]
                extracted_3UTR_region = curr_chr_coverage[0][left_region_index:right_region_index]
                extracted_3UTR_region.insert(0, start)
                extracted_3UTR_region.append(end)
                if curr_3UTR_event_id not in coverage_dict_all_samples:
                    coverage_dict_all_samples[curr_3UTR_event_id] = []
                coverage_dict_all_samples[curr_3UTR_event_id].append(
                    [extracted_coverage, extracted_3UTR_region])
        return coverage_dict_all_samples

    ##Load coverage for all samples
    depth_all_samples = []
    coverage_dict_all_samples = {}
    for wig_file in All_Wig_files:
        coverage_dict, depth = make_coverage_dict(wig_file)
        depth_all_samples.append(depth)
        coverage_dict_all_samples = get_utr_coverage(coverage_dict, utr_dict, coverage_dict_all_samples)
    return coverage_dict_all_samples, np.array(depth_all_samples)

class UtrFinder():
    '''
    This seems to be the main caller.
    '''
    def __init__(self, args):
        self.control_alignments = args.control_alignments
        self.treatment_alignments = args.treatment_alignments
        self.output_directory = args.output_directory
        self.utr = args.utr
        self.result_file = args.result_file
        #self.control_surpassing_coverage = args.control_surpassing_coverage
        #self.treatment_surpassing_coverage = args.treatment_surpassing_coverage
        #self.coverage_cutoff = args.coverage_cutoff
        #self.fdr_cutoff = args.fdr_cutoff
        #self.fold_change_cutoff = args.fold_change_cutoff
        #self.pdui_cutoff = args.pdui_cutoff
        #self.num_controls = len(control_alignments)
        self.all_alignments = self.control_alignments + self.treatment_alignments
        self.utr_dict = self.get_utr_dict(self)

    def get_utr_dict(self, shift=0.2):
        """Iterate over bed file. Shift 3' by param shift (default 20% of UTR length)
        If the shifted UTR is more than 50 nucleotides apart create a dictionary entry.
        dictionary structure is:
        {"gene":["chr1", 101, 221, (chr1, 100, 200)  ]}
        """
        utr_dict= {}
        for line in self.utr:
            fields = line.strip('\n').split('\t')
            chr = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            gene = fields[3]
            curr_strand = fields[5]
            UTR_pos = (chr, start, end)
            end_shift = int(round(abs(start - end) * shift))  # Shift 3' end by 20%
            if curr_strand == '+':
                end = end - end_shift
            else:
                start = start + end_shift
            start = start
            end = int(end)
            if start + 50 < end:
                utr_dict[gene] = [chr, start, end, fields[-1], UTR_pos]

def remaining_crap():
    pass
    #temp_dir = create_output_dir(output_directory)
    #Output_all_prediction_file = os.path.join(output_directory, Output_result_file, '_result_temp.txt')


    ##Prepare output directory

    #Output_result = open(Output_all_prediction_file, 'w')

    #num_samples = len(All_Sample_files)

    ##Debug
    print >> sys.stderr, "[%s] Loading coverage ..." % time_now()
    All_samples_Target_3UTR_coverages, All_samples_sequencing_depths, UTR_events_dict = process_alignments(All_Wig_files)

    All_sample_coverage_weights = All_samples_sequencing_depths / np.mean(All_samples_sequencing_depths)
    print >> sys.stderr, "[%s] Loading coverage finished ..." % time_now()
    ##Write the first line
    first_line = ['Gene', 'fit_value', 'Predicted_Proximal_APA', 'Loci']
    for i in range(num_group_1):
        curr_long_exp = 'A_%s_long_exp' % str(i + 1)
        curr_short_exp = 'A_%s_short_exp' % str(i + 1)
        curr_ratio = 'A_%s_PDUI' % str(i + 1)
        first_line.extend([curr_long_exp, curr_short_exp, curr_ratio])
    for i in range(num_samples - num_group_1):
        curr_long_exp = 'B_%s_long_exp' % str(i + 1)
        curr_short_exp = 'B_%s_short_exp' % str(i + 1)
        curr_ratio = 'B_%s_PDUI' % str(i + 1)
        first_line.extend([curr_long_exp, curr_short_exp, curr_ratio])
    first_line.append('PDUI_Group_diff')

    Output_result.writelines('\t'.join(first_line) + '\n')

    for curr_3UTR_id in UTR_events_dict:
        curr_3UTR_structure = UTR_events_dict[curr_3UTR_id]
        start = curr_3UTR_structure[1]
        end = curr_3UTR_structure[2]
        curr_strand = curr_3UTR_structure[-2]
        UTR_pos = curr_3UTR_structure[-1]
        if curr_3UTR_id in All_samples_Target_3UTR_coverages:
            curr_3UTR_coverage_wig = All_samples_Target_3UTR_coverages[curr_3UTR_id]
            curr_3UTR_all_samples_bp_coverage = []
            for curr_sample_curr_3UTR_coverage_wig in curr_3UTR_coverage_wig:
                curr_3UTR_curr_sample_bp_coverage = Convert_wig_into_bp_coverage(curr_sample_curr_3UTR_coverage_wig[0],
                                                                                 curr_sample_curr_3UTR_coverage_wig[1],
                                                                                 curr_strand)
                curr_3UTR_all_samples_bp_coverage.append(curr_3UTR_curr_sample_bp_coverage)

            select_mean_squared_error, selcted_break_point, UTR_abundances = De_Novo_3UTR_Coverage_estimation_Genome_for_TCGA_multiple_samples(
                curr_3UTR_all_samples_bp_coverage, start, end, curr_strand, All_sample_coverage_weights)

            if str(select_mean_squared_error) != "Na":
                Long_3UTR_exp_all = np.array(UTR_abundances[0])
                Short_3UTR_exp_all = np.array(UTR_abundances[1])
                num_non_zero = sum((Long_3UTR_exp_all + Short_3UTR_exp_all) > 0)
                if num_non_zero == num_samples:
                    All_Long_inclusion_ratios = []
                    line_write = [curr_3UTR_id, "%.1f" % select_mean_squared_error, str(selcted_break_point), UTR_pos]
                    for i in range(num_samples):
                        curr_sample_ratio = float(UTR_abundances[0][i]) / (
                        float(UTR_abundances[0][i]) + float(UTR_abundances[1][i]))  ##long 3'UTR percentage
                        All_Long_inclusion_ratios.append(curr_sample_ratio)
                        line_write.append("%.2f" % UTR_abundances[0][i])
                        line_write.append("%.2f" % UTR_abundances[1][i])
                        line_write.append("%.2f" % curr_sample_ratio)

                    Group1_IR = All_Long_inclusion_ratios[:num_group_1]
                    Group2_IR = All_Long_inclusion_ratios[num_group_1:]
                    inclusion_ratio_Group_diff = np.mean(np.array(Group1_IR)) - np.mean(np.array(Group2_IR))

                    line_write.append("%.2f" % inclusion_ratio_Group_diff)

                    Output_result.writelines('\t'.join(line_write) + '\n')

    Output_result.close()

    print >> sys.stderr, "[%s] Filtering the result ..." % time_now()

    Output_Motif_filtered_result_file = output_directory + Output_result_file + '_All_Prediction_Results.txt'
    # UTR_APA_Result_filtering(Output_all_prediction_file,Genome_seq_fasta,Output_Motif_filtered_result_file)

    DaPars_Filtering(Output_all_prediction_file, num_samples, num_group_1, Output_Motif_filtered_result_file)

    try:
        os.remove(Output_all_prediction_file)
    except OSError:
        pass

    try:
        os.rmdir(temp_dir)
    except OSError:
        pass

    print >> sys.stderr, "[%s] Finished!" % time_now()


# def DaPars_Filtering_debug():
def DaPars_Filtering(input_file, num_samples, num_group1, output_file):
    # cfg_file = 'CFIm25_Configure.txt'
    # control_alignment,treatment_alignment,output_directory,utr,Output_result_file,control_surpassing_coverage,treatment_surpassing_coverage,coverage_cutoff,fdr_cutoff,fold_change_cutoff,pdui_cutoff = parse_cfgfile(cfg_file)

    # input_file = 'CFIm25_KD_vs_Control_3UTR_All_prediction.txt'
    # input_file = 'Wagner_3UTR_New_Nov_5_2012_All_prediction.txt'
    # output_file = 'filtered.txt'
    # num_samples = 2
    # num_group1 = 1



    #     if fdr_cutoff != '':
    #         global FDR_cutoff
    #         FDR_cutoff = fdr_cutoff
    #         print FDR_cutoff
    #     if pdui_cutoff != '':
    #         global PDUI_cutoff
    #         PDUI_cutoff = pdui_cutoff
    #         print PDUI_cutoff





    output_write = open(output_file, 'w')
    num_line = 0

    result_dict = {}
    All_P_values = []
    Selected_events_id = []
    All_mean_abundance = []

    for line in open(input_file, 'r'):
        if num_line > 0:
            fields = line.strip('\n').split('\t')
            group1_coverages = np.zeros(2)
            group2_coverages = np.zeros(2)
            num_group1_pass = 0
            group1_PDUIs = 0
            for i in range(num_group1):
                curr_long = fields[4 + i * 3]
                curr_short = fields[5 + i * 3]
                if curr_long != 'NA':
                    curr_long = float(curr_long)
                    curr_short = float(curr_short)
                    if curr_long + curr_short >= Coverage_cutoff:
                        group1_PDUIs = group1_PDUIs + float(fields[6 + i * 3])
                        num_group1_pass += 1
                        group1_coverages[0] = group1_coverages[0] + curr_long
                        group1_coverages[1] = group1_coverages[1] + curr_short
                    else:
                        fields[4 + i * 3] = 'NA'
                        fields[5 + i * 3] = 'NA'
                        fields[6 + i * 3] = 'NA'

            num_group2_pass = 0
            group2_PDUIs = 0
            for i in range(num_samples - num_group1):
                curr_long = fields[4 + (i + num_group1) * 3]
                curr_short = fields[5 + (i + num_group1) * 3]
                if curr_long != 'NA':
                    curr_long = float(curr_long)
                    curr_short = float(curr_short)
                    if curr_long + curr_short >= Coverage_cutoff:
                        group2_PDUIs = group2_PDUIs + float(fields[6 + (i + num_group1) * 3])
                        num_group2_pass += 1
                        group2_coverages[0] = group2_coverages[0] + curr_long
                        group2_coverages[1] = group2_coverages[1] + curr_short
                    else:
                        fields[4 + (i + num_group1) * 3] = 'NA'
                        fields[5 + (i + num_group1) * 3] = 'NA'
                        fields[6 + (i + num_group1) * 3] = 'NA'

            if num_group1_pass >= Num_least_in_group1 and num_group2_pass >= Num_least_in_group2:
                Final_group_diff = str(group1_PDUIs / num_group1_pass - group2_PDUIs / num_group2_pass)

                All_mean_abundance.append([group1_PDUIs / num_group1_pass, group2_PDUIs / num_group2_pass])

                fields[-1] = str(Final_group_diff)
                ratio_val, P_val = sp.stats.fisher_exact(
                    [group1_coverages / num_group1_pass, group2_coverages / num_group2_pass])

                All_P_values.append(P_val)
                Selected_events_id.append(fields[0])
                # print P_val
                # print ratio_val
            else:
                fields[-1] = 'NA'

            result_dict[fields[0]] = fields

        else:
            first_line = line.strip('\n').split('\t')

        num_line += 1


    ##Filtering
    stats = importr('stats')
    All_p_adjust = stats.p_adjust(FloatVector(All_P_values), method='BH')
    first_line.insert(-1, 'Group_A_Mean_PDUI')
    first_line.insert(-1, 'Group_B_Mean_PDUI')
    first_line.extend(['P_val', 'adjusted.P_val', 'Pass_Filter'])
    output_write.writelines('\t'.join(first_line) + '\n')
    for curr_event_id in result_dict:
        mean_PDUI_group1 = 'NA'
        mean_PDUI_group2 = 'NA'
        curr_P_val = 'NA'
        curr_FDR_val = 'NA'
        Pass_filter = 'N'
        curr_fields = result_dict[curr_event_id]
        if curr_event_id in Selected_events_id:
            sel_ind = Selected_events_id.index(curr_event_id)
            curr_P_val = str(All_P_values[sel_ind])
            curr_FDR_val = str(All_p_adjust[sel_ind])

            mean_PDUI_group1 = All_mean_abundance[sel_ind][0]
            mean_PDUI_group2 = All_mean_abundance[sel_ind][1]

            if float(curr_FDR_val) <= FDR_cutoff and abs(float(curr_fields[-1])) >= PDUI_cutoff and abs(
                    math.log((mean_PDUI_group1 + 1e-5) / (mean_PDUI_group2 + 1e-5), 2)) >= Fold_change_cutoff:
                Pass_filter = 'Y'

        curr_fields.insert(-1, str(mean_PDUI_group1))
        curr_fields.insert(-1, str(mean_PDUI_group2))
        curr_fields.append(curr_P_val)
        curr_fields.append(curr_FDR_val)
        curr_fields.append(Pass_filter)

        output_write.writelines('\t'.join(curr_fields) + '\n')

    output_write.close()


def get_version():
    return "0.9.1"


def De_Novo_3UTR_Coverage_estimation_Genome_for_TCGA_multiple_samples(All_Samples_curr_3UTR_coverages, UTR_start,
                                                                      UTR_end, curr_strand, weight_for_second_coverage):
    '''For UTR-APA new
       Load one chromosome by chromosome
       Just for TCGA data analysis. So no peak evenness checking
       Jan-17-2013
       2-28-2013
    '''
    coverage_threshold = 20
    search_point_start = 200
    search_point_end = int(abs((UTR_end - UTR_start)) * 0.1)

    num_samples = len(All_Samples_curr_3UTR_coverages)
    ##read coverage
    Region_Coverages = []
    Region_mean_Coverages = []
    Region_first_100_coverage_all_samples = []
    for i in range(num_samples):
        curr_Region_Coverage_raw = All_Samples_curr_3UTR_coverages[i]  ##strand is reversed in load
        curr_Region_Coverage = curr_Region_Coverage_raw / weight_for_second_coverage[i]
        Region_mean_Coverages.append(np.mean(curr_Region_Coverage_raw))
        Region_Coverages.append(curr_Region_Coverage)
        curr_first_100_coverage = np.mean(curr_Region_Coverage_raw[0:99])
        Region_first_100_coverage_all_samples.append(curr_first_100_coverage)
    if sum(np.array(
            Region_first_100_coverage_all_samples) >= coverage_threshold) >= num_samples and UTR_end - UTR_start >= 150:
        if curr_strand == "+":
            search_region = range(UTR_start + search_point_start, UTR_end - search_point_end + 1)
        else:
            search_region = range(UTR_end - search_point_start, UTR_start + search_point_end - 1, -1)

        search_start = search_point_start
        search_end = UTR_end - UTR_start - search_point_end
        Mean_squared_error_list = []
        Estimated_3UTR_abundance_list = []
        for curr_point in range(search_start, search_end + 1):
            curr_search_point = curr_point
            All_samples_result = [[], [], []]
            for curr_sample_region_coverage in Region_Coverages:
                Mean_Squared_error, Long_UTR_abun, Short_UTR_abun = Estimation_abundance(curr_sample_region_coverage,
                                                                                         curr_search_point)
                All_samples_result[0].append(Mean_Squared_error)
                All_samples_result[1].append(Long_UTR_abun)
                All_samples_result[2].append(Short_UTR_abun)

            Mean_Squared_error = np.mean(np.array(All_samples_result[0]))
            Mean_squared_error_list.append(Mean_Squared_error)
            Estimated_3UTR_abundance_list.append([All_samples_result[1], All_samples_result[2]])

        if len(Mean_squared_error_list) > 0:
            min_ele_index = Mean_squared_error_list.index(min(Mean_squared_error_list))

            select_mean_squared_error = Mean_squared_error_list[min_ele_index]
            UTR_abundances = Estimated_3UTR_abundance_list[min_ele_index]
            selcted_break_point = search_region[min_ele_index]

        else:
            select_mean_squared_error = 'Na'
            UTR_abundances = 'Na'
            selcted_break_point = 'Na'

    else:
        select_mean_squared_error = 'Na'
        UTR_abundances = 'Na'
        selcted_break_point = 'Na'

    return select_mean_squared_error, selcted_break_point, UTR_abundances


def Estimation_abundance(Region_Coverage, break_point):
    Long_UTR_abun = np.mean(Region_Coverage[break_point:])
    Short_UTR_abun = np.mean(Region_Coverage[0:break_point] - Long_UTR_abun)
    if Short_UTR_abun < 0:
        Short_UTR_abun = 0
    Coverage_diff = Region_Coverage[0:break_point] - Long_UTR_abun - Short_UTR_abun
    Coverage_diff = np.append(Coverage_diff, Region_Coverage[break_point:] - Long_UTR_abun)
    Mean_Squared_error = np.mean(Coverage_diff ** 2)

    return Mean_Squared_error, Long_UTR_abun, Short_UTR_abun

if __name__ == '__main__':

    args = parse_args(sys.argv)
    find_utr = UtrFinder(args)

    ##debug
    # DaPars_Filtering_debug()
