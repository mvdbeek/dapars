import argparse
import numpy as np
import os
import sys
import datetime
from bisect import bisect
from collections import OrderedDict
import pysam


def time_now():  # return time
    curr_time = datetime.datetime.now()
    return curr_time.strftime("%c")


def get_bp_coverage(alignment, utr, is_reverse=False):
    """
    Use pysam pileup method on indexed bam alignment files to retrieve strand-specific coverage in numpy array format.
    :param alignment:
    :param utr:
    :param is_reverse:
    :return:
    """
    """
    :param alignment:
    :param utr:
    :param is_reverse:
    :return:
    """
    pos_n = {pu.pos:sum([ int(pur.alignment.is_reverse == is_reverse) for pur in pu.pileups ]) for pu in alignment.pileup(*utr)}
    pus = np.array([pos_n.get(i, 0) for i in xrange(utr[1], utr[2])])
    if not is_reverse:
        return pus
    else:
        return pus[::-1]

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
    parser.add_argument("-t", "--treatment_alignments", nargs="+", required=True,
                        help="Alignment files in BAM format from treatment condition")
    parser.add_argument("-u", "--utr_bed_file", required=True,
                        help="Bed file describing longest 3UTR positions")
    parser.add_argument("-o", "--output_file", required=True, type=argparse.FileType('w'),
                        help="file containing output")
    return parser.parse_args()


class UtrFinder():
    '''
    This seems to be the main caller.
    '''

    def __init__(self, args):
        self.control_alignments = [pysam.AligmentFile(file) for file in args.control_alignments]
        self.treatment_alignments = [pysam.AligmentFile(file) for file in  args.treatment_alignments]
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
            if structure[-2] == "+"
                is_reverse = False
            else:
                is_reverse = True
            UTR_pos = structure[-1]
            if utr in coverage_dict_all_samples:
                coverage_wig = coverage_dict_all_samples[utr]
                utr_bp_coverage = []
                for bam in self.all_alignments:
                    bp_coverage = get_bp_coverage(bam, utr, is_reverse)
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
                                   UTR_end, is_reverse, weight_for_second_coverage):
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
        if not is_reverse:
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
