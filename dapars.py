import argparse
import os
import numpy as np
import pysam
from collections import OrderedDict, namedtuple
import filter_utr


def get_bp_coverage(alignment, chr, start, end, is_reverse=False):
    """
    Use pysam pileup method on indexed bam alignment files to retrieve strand-specific coverage in numpy array format.
    """
    pos_n = {pu.pos:sum([ int(pur.alignment.is_reverse == is_reverse) for pur in pu.pileups ]) for pu in alignment.pileup(chr, start, end)}
    pus = np.array([pos_n.get(i, 0) for i in xrange(start, end)])
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
    parser.add_argument("-c", "--control_alignments", nargs="+", required=True,
                        help="Alignment files in BAM format from control condition")
    parser.add_argument("-t", "--treatment_alignments", nargs="+", required=True,
                        help="Alignment files in BAM format from treatment condition")
    parser.add_argument("-u", "--utr_bed_file", required=True, type=file,
                        help="Bed file describing longest 3UTR positions")
    parser.add_argument("-o", "--output_file", required=True, type=argparse.FileType('w'),
                        help="file containing output")
    return parser.parse_args()


class UtrFinder():
    """
    This seems to be the main caller.
    """

    def __init__(self, args):
        self.control_alignments = [pysam.AlignmentFile(file) for file in args.control_alignments]
        self.treatment_alignments = [pysam.AlignmentFile(file) for file in args.treatment_alignments]
        self.available_chromosomes = self.get_available_chromosome()
        self.utr = args.utr_bed_file
        self.gtf_fields = filter_utr.get_gtf_fields()
        self.result_file = args.output_file
        self.all_alignments = self.control_alignments + self.treatment_alignments
        self.num_samples = len(self.all_alignments)
        self.utr_dict = self.get_utr_dict(0.2)
        self.utr_coverages = self.get_utr_coverage()
        self.coverage_weights = self.get_coverage_weights()
        result_d = self.calculate_apa_ratios()
        test = result_d

    def get_available_chromosome(self):
        return self.control_alignments[0].references

    def get_utr_dict(self, shift):
        utr_dict = OrderedDict()
        for line in self.utr:
            if not line.startswith("#"):
                filter_utr.get_utr_dict( line, self.gtf_fields, utr_dict )
                gene, utr_d = utr_dict.popitem()
                utr_d = utr_d[0]
                end_shift = int(round(abs(utr_d["start"] - utr_d["end"]) * shift))
                if utr_d["strand"] == "+":
                    utr_d["new_end"] = utr_d["end"] - end_shift
                    utr_d["new_start"] = utr_d["start"]
                else:
                    utr_d["new_end"] = utr_d["end"]
                    utr_d["new_start"] = utr_d["start"] + end_shift
                if utr_d["new_start"] + 50 < utr_d["new_end"]:
                    utr_dict[gene] = utr_d
        return utr_dict

    def get_utr_dict_from_bed(self, shift):
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

    def get_utr_coverage(self):
        """
        Returns a dict:
        { UTR : [coverage_aligment1, ...]}
        """
        utr_coverages = {}
        for utr, utr_d in self.utr_dict.iteritems():
            if utr_d["chr"] in self.available_chromosomes:
                if utr_d["strand"] == "+":
                    is_reverse = False
                else:
                    is_reverse = True
                utr_coverage = []
                for bam in self.all_alignments:
                    bp_coverage = get_bp_coverage(bam, utr_d["chr"], utr_d["new_start"], utr_d["new_end"], is_reverse)
                    utr_coverage.append(bp_coverage)
                utr_coverages[utr] = utr_coverage
        return utr_coverages

    def get_coverage_weights(self):
        """
        Return weights for normalizing coverage.
        """
        coverage_per_alignment = []
        for utr in self.utr_coverages.itervalues():
            utr_coverage = []
            for pileup in utr:
                utr_coverage.append(np.sum(pileup))
            coverage_per_alignment.append(utr_coverage)
        coverages = np.array([ sum(x) for x in zip(*coverage_per_alignment) ])
        coverage_weights = coverages / np.mean(coverages)  # TODO: proabably median is better suited?
        return coverage_weights

    def get_result_tuple(self):
        static_desc = ["chr", "start", "end", "strand", "gene", "breakpoint", "control_mean_percent", "treatment_mean_percent" ]
        samples_desc = []
        for statistic in ["coverage_long", "coverage_short", "percent_long"]:
            for i, sample in enumerate(self.control_alignments):
                samples_desc.append("control_{i}_{statistic}".format(i=i+1, statistic = statistic))
            for i, sample in enumerate(self.treatment_alignments):
                samples_desc.append("treatment_{i}_{statistic}".format(i=i+1, statistic = statistic))
        return namedtuple("result", static_desc + samples_desc, rename=True)



    def calculate_apa_ratios(self):
        result_d = OrderedDict()
        for utr, utr_d in self.utr_dict.iteritems():
            result_tuple = self.get_result_tuple()
            if utr_d["strand"] == "+":
                is_reverse = False
            else:
                is_reverse = True
            utr_coverage = self.utr_coverages[utr]
            mse, breakpoint, abundances = estimate_coverage_extended_utr(utr_coverage,
                                                                         utr_d["new_start"],
                                                                         utr_d["new_end"],
                                                                         is_reverse,
                                                                         self.coverage_weights)
            if not str(mse) == "Na":
                long_coverage_vector = abundances[0]
                short_coverage_vector = abundances[1]
                num_non_zero = sum((np.array(long_coverage_vector) + np.array(short_coverage_vector)) > 0)  # TODO: This introduces bias
                if num_non_zero == self.num_samples:
                    percentage_long = []
                    for i in range(self.num_samples):
                        ratio = float(long_coverage_vector[i]) / (long_coverage_vector[i] + short_coverage_vector[i])  # long 3'UTR percentage
                        percentage_long.append(ratio)
                    control_mean_percent = np.mean(np.array(percentage_long[:len(self.control_alignments)]))
                    treatment_mean_percent = np.mean(np.array(percentage_long[len(self.control_alignments):]))
                    res = result_tuple(utr_d["chr"], utr_d["start"], utr_d["end"], utr_d["strand"],
                                       utr, breakpoint, control_mean_percent, treatment_mean_percent,
                                       *zip(percentage_long, long_coverage_vector, short_coverage_vector))
                    result_d[utr] = res
        return result_d


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


def estimate_coverage_extended_utr(utr_coverage, UTR_start,
                                   UTR_end, is_reverse, coverage_weigths):
    """
    We are searching for a breakpoint in coverage?!
    utr_coverage is a list with items corresponding to numpy arrays of coverage for a sample.
    """
    coverage_threshold = 15
    search_point_start = 200
    search_point_end = int(abs((UTR_end - UTR_start)) * 0.1)  # TODO: This is 10% of total UTR end. Why?
    num_samples = len(utr_coverage)
    ##read coverage filtering
    normalized_utr_coverage = [coverage/ coverage_weigths[i] for i, coverage in enumerate( utr_coverage )]
    start_coverage = [np.mean(coverage[0:99]) for coverage in utr_coverage]  # filters threshold on mean coverage over first 100 nt
    is_above_threshold = sum(np.array(start_coverage) >= coverage_threshold) >= num_samples  # This filters on the raw threshold. Why?
    is_above_length = UTR_end - UTR_start >= 150
    if (is_above_threshold) and (is_above_length):
        if not is_reverse:
            search_region = range(UTR_start + search_point_start, UTR_end - search_point_end + 1)
        else:
            search_region = range(UTR_end - search_point_start, UTR_start + search_point_end - 1, -1)
        search_start = search_point_start
        search_end = UTR_end - UTR_start - search_point_end
        mse_list = []
        Estimated_3UTR_abundance_list = []
        for curr_point in range(search_start, search_end + 1):
            curr_search_point = curr_point
            All_samples_result = [[], [], []]
            for curr_sample_region_coverage in normalized_utr_coverage:
                Mean_Squared_error, Long_UTR_abun, Short_UTR_abun = estimate_abundance(curr_sample_region_coverage,
                                                                                       curr_search_point)
                All_samples_result[0].append(Mean_Squared_error)
                All_samples_result[1].append(Long_UTR_abun)
                All_samples_result[2].append(Short_UTR_abun)
            Mean_Squared_error = np.mean(np.array(All_samples_result[0]))
            mse_list.append(Mean_Squared_error)
            Estimated_3UTR_abundance_list.append([All_samples_result[1], All_samples_result[2]])
        if len(mse_list) > 0:
            min_ele_index = mse_list.index(min(mse_list))
            select_mean_squared_error = mse_list[min_ele_index]
            UTR_abundances = Estimated_3UTR_abundance_list[min_ele_index]
            selected_break_point = search_region[min_ele_index]
        else:
            select_mean_squared_error = 'Na'
            UTR_abundances = 'Na'
            selected_break_point = 'Na'
    else:
        select_mean_squared_error = 'Na'
        UTR_abundances = 'Na'
        selected_break_point = 'Na'

    return select_mean_squared_error, selected_break_point, UTR_abundances


def estimate_abundance(Region_Coverage, breakpoint):
    """
    get abundance of long utr vs short utr with breakpoint specifying the position of long and short utr.
    """
    long_utr_vector = Region_Coverage[breakpoint:]
    short_utr_vector = Region_Coverage[0:breakpoint]
    mean_long_utr = np.mean(long_utr_vector)
    mean_short_utr = np.mean(Region_Coverage[0:breakpoint])
    square_mean_centered_short_utr_vector = (short_utr_vector - mean_short_utr) **2
    square_mean_centered_long_utr_vector = (long_utr_vector - mean_long_utr) ** 2
    mse = np.mean(np.append(square_mean_centered_short_utr_vector, square_mean_centered_long_utr_vector))

    return mse, mean_long_utr, mean_short_utr


if __name__ == '__main__':
    args = parse_args()
    find_utr = UtrFinder(args)
