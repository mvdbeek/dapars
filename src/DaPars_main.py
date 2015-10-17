import argparse
import numpy as np
import os
import sys
import datetime
from bisect import bisect
from collections import OrderedDict


def time_now():  # return time
    curr_time = datetime.datetime.now()
    return curr_time.strftime("%c")


def wig_to_bp_coverage(extracted_coverage, extracted_3UTR_region, strand_info):
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


def parse_args():
    """
    Returns floating point values except for input files.
    My initial approach will not filter anything. (FDR. fold_change, PDUI, Num_least ...)
    :param argv:
    :return:
    """
    parser = argparse.ArgumentParser(description='Determines the usage of proximal polyA usage')
    parser.add_argument("-c", "--control_alignments", nargs="+", required=True, type=file,
                        help="Alignment files in BAM format from control condition")
    parser.add_argument("-t", "--treatment_alignments", nargs="+", required=True, type=file,
                        help="Alignment files in BAM format from treatment condition")
    parser.add_argument("-u", "--utr_bed_file", required=True, type=file,
                        help="Bed file describing longest 3UTR positions")
    parser.add_argument("-o", "--output_file", required=True, type=argparse.FileType('w'),
                        help="file containing output")
    return parser.parse_args()


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
        """Get coverage_dict for wig_file".
        covarage dict structure is:
        {"chr1": [ [position], [count] ]}
        Depth is total covered based (WEIRD!)"""
        coverage_dict = {}
        depth = 0
        for line in wig_file:
            if '#' not in line and line[0:3] == 'chr':
                fields = line.strip('\n').split('\t')
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                count = int(fields[-1])
                depth += int(fields[-1]) * (end - start)
                if chrom not in coverage_dict:
                    coverage_dict[chrom] = [[0], [0]]
                if start > coverage_dict[chrom][0][-1]:  # No strand info in wig files: buhuuuu!
                    coverage_dict[chrom][0].append(start)
                    coverage_dict[chrom][1].append(0)
                coverage_dict[chrom][0].append(end)
                coverage_dict[chrom][1].append(count)
        coverage_dict[chrom][1].append(0)
        return coverage_dict, depth

    def get_utr_coverage(coverage_dict, utr_dict, coverage_dict_all_samples):
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
        self.utr = args.utr_bed_file
        self.result_file = args.output_file
        self.all_alignments = self.control_alignments + self.treatment_alignments
        self.utr_dict = self.get_utr_dict(0.2)
        self.coverage_dict_all_samples, self.coverage_weights = get_coverage_dicts(self)
        result_d = self.calculate_apa_ratios()

    def get_utr_dict(self, shift):
        """Iterate over bed file. Shift 3' by param shift (default 20% of UTR length)
        If the shifted UTR is more than 50 nucleotides apart create a dictionary entry.
        dictionary structure is:
        {"gene":["chr1", 101, 221, (chr1, 100, 200)  ]}
        """
        utr_dict = OrderedDict()
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
        return utr_dict

    def calculate_apa_ratios(self):
        utr_dict = self.utr_dict
        coverage_dict_all_samples = self.coverage_dict_all_samples
        coverage_weights = self.coverage_weights
        num_samples = coverage_weights.size
        result_d = OrderedDict()
        for utr, structure in utr_dict.iteritems():
            start = structure[1]
            end = structure[2]
            strand = structure[-2]
            UTR_pos = structure[-1]
            if utr in coverage_dict_all_samples:
                coverage_wig = coverage_dict_all_samples[utr]
                utr_bp_coverage = []
                for wig in coverage_wig:
                    bp_coverage = wig_to_bp_coverage(wig[0], wig[1], strand)
                    utr_bp_coverage.append(bp_coverage)

                mse, breakpoint, abundances = estimate_coverage_extended_utr(
                    utr_bp_coverage, start, end, strand, coverage_weights)

                if not str(mse) == "Na":
                    long_exp_all = np.array(abundances[0])
                    short_exp_all = np.array(abundances[1])
                    num_non_zero = sum((long_exp_all + short_exp_all) > 0)  # TODO: This introduces bias
                    if num_non_zero == num_samples:
                        percentages_long = []
                        result_d[utr] = [mse, breakpoint, UTR_pos]
                        for i in range(num_samples):
                            ratio = float(abundances[0][i]) / (
                                float(abundances[0][i]) + float(abundances[1][i]))  # long 3'UTR percentage
                            percentages_long.append(ratio)
                            result_d[utr].extend((abundances[0][i], abundances[1][i], ratio))

                        Group1_IR = percentages_long[:len(self.control_alignments)]
                        Group2_IR = percentages_long[len(self.control_alignments):]
                        inclusion_ratio_Group_diff = np.mean(np.array(Group1_IR)) - np.mean(np.array(Group2_IR))

                        result_d[utr].extend([inclusion_ratio_Group_diff])
        return result_d


def get_coverage_dicts(self):
    """
    """
    coverage_dict_all_samples, depth_all_samples = process_alignments(self.all_alignments, self.utr_dict)
    coverage_weights = depth_all_samples / np.mean(depth_all_samples)
    return coverage_dict_all_samples, coverage_weights


def write_header():
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

    self.result_file.writelines('\t'.join(first_line) + '\n')


def estimate_coverage_extended_utr(All_Samples_curr_3UTR_coverages, UTR_start,
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
                Mean_Squared_error, Long_UTR_abun, Short_UTR_abun = estimate_abundance(curr_sample_region_coverage,
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


def estimate_abundance(Region_Coverage, break_point):
    Long_UTR_abun = np.mean(Region_Coverage[break_point:])
    Short_UTR_abun = np.mean(Region_Coverage[0:break_point] - Long_UTR_abun)
    if Short_UTR_abun < 0:
        Short_UTR_abun = 0
    Coverage_diff = Region_Coverage[0:break_point] - Long_UTR_abun - Short_UTR_abun
    Coverage_diff = np.append(Coverage_diff, Region_Coverage[break_point:] - Long_UTR_abun)
    mse = np.mean(Coverage_diff ** 2)

    return mse, Long_UTR_abun, Short_UTR_abun


if __name__ == '__main__':
    args = parse_args()
    find_utr = UtrFinder(args)

    ##debug
    # DaPars_Filtering_debug()
