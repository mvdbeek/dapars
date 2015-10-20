import argparse
import os
import csv
import numpy as np
from collections import OrderedDict, namedtuple
import filter_utr
import subprocess


def get_bp_coverage(alignment, chr, start, end, is_reverse=False):
    """
    Use pysam pileup method on indexed bam alignment files to retrieve strand-specific coverage in numpy array format.
    This is way too slow, unfortunately.
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
        self.control_alignments = [file for file in args.control_alignments]
        self.treatment_alignments = [file for file in args.treatment_alignments]
        #self.available_chromosomes = self.get_available_chromosome()
        self.utr = args.utr_bed_file
        self.gtf_fields = filter_utr.get_gtf_fields()
        self.result_file = args.output_file
        self.all_alignments = self.control_alignments + self.treatment_alignments
        self.alignment_names = { file: os.path.basename(file) for file in self.all_alignments }
        self.num_samples = len(self.all_alignments)
        self.utr_dict = self.get_utr_dict(0.2)
        self.dump_utr_dict_to_bedfile()
        print "done getting utr_dict"
        self.coverage_files = self.run_bedtools_coverage()
        print "produced coverage files"
        self.utr_coverages = self.read_coverage_result()
        print "got coverage_dict"
        self.coverage_weights = self.get_coverage_weights()
        self.result_tuple = self.get_result_tuple()
        self.result_d = self.calculate_apa_ratios()
        self.write_results()


    def dump_utr_dict_to_bedfile(self):
        w = csv.writer(open("tmp_bedfile.bed", "w"), delimiter="\t")
        for gene, utr in self.utr_dict.iteritems():
            w.writerow([utr["chr"], utr["new_start"]-1, utr["new_end"], gene, ".", utr["strand"]])

    def run_bedtools_coverage(self):
        """
        Use bedtools coverage to generate pileup data for all alignment files for the regions specified in utr_dict.
        """
        coverage_files = []
        for alignment_file in self.all_alignments:
            cmd = "sort -k1,1 -k2,2n tmp_bedfile.bed | "
            cmd = cmd + "/usr/bin/bedtools coverage -d -s -abam {alignment_file} -b stdin |" \
                        " cut -f 4,7,8 > coverage_file_{alignment_name}".format(
                alignment_file = alignment_file, alignment_name= self.alignment_names[alignment_file] )
            print cmd
            subprocess.call([cmd], shell=True)
            coverage_files.append("gene_position_coverage_{alignment_name}".format(
                alignment_name = self.alignment_names[alignment_file]))
        return coverage_files

    def read_coverage_result(self):
        """
        Read coverages back in and store as dictionary of numpy arrays
        """
        coverage_dict = { gene: { file: np.zeros(utr_d["new_end"]+1-utr_d["new_start"]) for file in self.all_alignments } for gene, utr_d in self.utr_dict.iteritems() }
        for alignment_name in self.alignment_names.itervalues():
            with open("coverage_file_{alignment_name}".format(alignment_name = alignment_name)) as coverage_file:
                for line in coverage_file:
                    gene, position, coverage= line.strip().split("\t")
                    coverage_dict[gene][alignment_name][int(position)-1] = coverage
        for utr_d in self.utr_dict.itervalues():
            if utr_d["strand"] == "-":
                for alignment_name in self.alignment_names.values():
                    coverage_dict[gene][alignment_name] = coverage_dict[gene][alignment_name][::-1]
        return coverage_dict

#    def get_available_chromosome(self):
#        return self.control_alignments[0].references

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
        utr_coverage is still confusing.
        """
        coverage_per_alignment = []
        for utr in self.utr_coverages.itervalues():  # TODO: be smarter about this.
            utr_coverage = []
            for vector in utr.itervalues():
                utr_coverage.append(np.sum(vector))
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
                    stats = zip(long_coverage_vector, short_coverage_vector, percentage_long)
                    stats = [item for sublist in stats for item in sublist]
                    res = self.result_tuple(utr_d["chr"], utr_d["start"], utr_d["end"], utr_d["strand"],
                                       utr, breakpoint, control_mean_percent, treatment_mean_percent,
                                       *stats)
                    result_d[utr] = res
        return result_d

    def write_results(self):
        w = csv.writer(self.result_file, delimiter='\t')
        w.writerow((self.result_tuple._fields))    # field header
        w.writerows( self.result_d.values())


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
    normalized_utr_coverage = [coverage/ coverage_weigths[i] for i, coverage in enumerate( utr_coverage.values() )]
    start_coverage = [np.mean(coverage[0:99]) for coverage in utr_coverage.values()]  # filters threshold on mean coverage over first 100 nt
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

