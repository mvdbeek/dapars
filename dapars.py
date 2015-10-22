import argparse
import os
import csv
import numpy as np
from collections import OrderedDict, namedtuple
import filter_utr
import subprocess
import time
from multiprocessing import Pool

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
        self.start = time.time()
        self.result_d = self.calculate_apa_ratios()
        print "numpy took {0}".format(time.time() - self.start)
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
            cmd = cmd + "bedtools coverage -d -s -abam {alignment_file} -b stdin |" \
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
                samples_desc.append("control_{i}_{statistic}".format(i=i, statistic = statistic))
            for i, sample in enumerate(self.treatment_alignments):
                samples_desc.append("treatment_{i}_{statistic}".format(i=i, statistic = statistic))
        return namedtuple("result", static_desc + samples_desc)



    def calculate_apa_ratios(self):
        result_d = OrderedDict()
        import time
        start_time = time.time()
        arg_d = {"result_tuple": self.result_tuple,
                 "coverage_weights":self.coverage_weights,
                 "num_samples":self.num_samples,
                 "num_control":len(self.control_alignments),
                 "num_treatment":len(self.treatment_alignments),
                 "result_d":result_d}
        #[ calculate_all_utr(self.utr_coverages[utr], utr, utr_d, **arg_d) for utr, utr_d in self.utr_dict.iteritems() ]
        pool = Pool(2)
        tasks = [ (self.utr_coverages[utr], utr, utr_d, self.result_tuple._fields, self.coverage_weights, self.num_samples,
                    len(self.control_alignments), len(self.treatment_alignments), result_d) for utr, utr_d in self.utr_dict.iteritems() ]
        result = [ pool.apply_async(calculate_all_utr, t) for t in tasks]
        [res.get() for res in result]
        for utr, utr_d in result_d:
            result_d[utr] = self.result_tuple(*utr_d)
        stop_time = time.time() - start_time
        print "took {0} seconds".format(stop_time)
        return result_d

    def write_results(self):
        w = csv.writer(self.result_file, delimiter='\t')
        w.writerow((self.result_tuple._fields))    # field header
        w.writerows( self.result_d.values())

def calculate_all_utr(utr_coverage, utr, utr_d, result_tuple_fields, coverage_weights, num_samples, num_control, num_treatment, result_d):
    res = dict(zip(result_tuple_fields, result_tuple_fields))
    if utr_d["strand"] == "+":
        is_reverse = False
    else:
        is_reverse = True
    mse, breakpoint, abundances = estimate_coverage_extended_utr(utr_coverage,
                                                                 utr_d["new_start"],
                                                                 utr_d["new_end"],
                                                                 is_reverse,
                                                                 coverage_weights)
    if not str(mse) == "Na":
        long_coverage_vector = abundances[0]
        short_coverage_vector = abundances[1]
        num_non_zero = sum((np.array(long_coverage_vector) + np.array(short_coverage_vector)) > 0)  # TODO: This introduces bias
        if num_non_zero == num_samples:
            percentage_long = []
            for i in range(num_samples):
                ratio = float(long_coverage_vector[i]) / (long_coverage_vector[i] + short_coverage_vector[i])  # long 3'UTR percentage
                percentage_long.append(ratio)
            for i in range(num_control):
                res["control_{i}_coverage_long".format(i=i)] = float(long_coverage_vector[i])
                res["control_{i}_coverage_short".format(i=i)] = float(short_coverage_vector[i])
                res["control_{i}_percent_long".format(i=i)] = percentage_long[i]
            for k in range(num_treatment):
                i = k + num_control
                res["treatment_{i}_coverage_long".format(i=k)] = float(long_coverage_vector[i])
                res["treatment_{i}_coverage_short".format(i=k)] = float(short_coverage_vector[i])
                res["treatment_{i}_percent_long".format(i=k)] = percentage_long[i]
            control_mean_percent = np.mean(np.array(percentage_long[:num_control]))
            treatment_mean_percent = np.mean(np.array(percentage_long[num_control:]))
            res["chr"] = utr_d["chr"]
            res["start"] = utr_d["start"]
            res["end"] = utr_d["end"]
            res["strand"] = utr_d["strand"]
            res["breakpoint"] = breakpoint
            res["control_mean_percent"] = control_mean_percent
            res["treatment_mean_percent"] = treatment_mean_percent
            result_d[utr] = res


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
        normalized_utr_coverage = np.array(normalized_utr_coverage)
        breakpoints = range(search_start, search_end + 1)
        mse_list = [ estimate_mse(normalized_utr_coverage,bp, num_samples) for bp in breakpoints ]
        if len(mse_list) > 0:
            min_ele_index = mse_list.index(min(mse_list))
            breakpoint = breakpoints[min_ele_index]
            UTR_abundances = estimate_abundance(normalized_utr_coverage, breakpoint, num_samples)
            select_mean_squared_error = mse_list[min_ele_index]
            selected_break_point = breakpoint
        else:
            select_mean_squared_error = 'Na'
            UTR_abundances = 'Na'
            selected_break_point = 'Na'
    else:
        select_mean_squared_error = 'Na'
        UTR_abundances = 'Na'
        selected_break_point = 'Na'

    return select_mean_squared_error, selected_break_point, UTR_abundances


def estimate_mse(cov, bp, num_samples):
    """
    get abundance of long utr vs short utr with breakpoint specifying the position of long and short utr.
    """
    long_utr_vector = cov[:num_samples, bp:]
    short_utr_vector = cov[:num_samples, 0:bp]
    mean_long_utr = np.mean(long_utr_vector, 1)
    mean_short_utr = np.mean(short_utr_vector, 1)
    square_mean_centered_short_utr_vector = (short_utr_vector[:num_samples] - mean_short_utr[:, np.newaxis] )**2
    square_mean_centered_long_utr_vector = (long_utr_vector[:num_samples] - mean_long_utr[:, np.newaxis])**2
    mse = np.mean(np.append(square_mean_centered_short_utr_vector[:num_samples], square_mean_centered_long_utr_vector[:num_samples]))
    return mse

def estimate_abundance(cov, bp, num_samples):
    long_utr_vector = cov[:num_samples, bp:]
    short_utr_vector = cov[:num_samples, 0:bp]
    mean_long_utr = np.mean(long_utr_vector, 1)
    mean_short_utr = np.mean(short_utr_vector, 1)
    return mean_long_utr, mean_short_utr


if __name__ == '__main__':
    args = parse_args()
    find_utr = UtrFinder(args)

