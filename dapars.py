import argparse
import os
import csv
import numpy as np
from scipy import stats
from collections import OrderedDict, namedtuple
import filter_utr
import subprocess
from multiprocessing import Pool
import warnings
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from tabulate import tabulate

def directory_path(str):
    if os.path.exists(str):
        return str
    else:
        os.mkdir(str)
        return str

def parse_args():
    """
    Returns floating point values except for input files.
    My initial approach will not filter anything. (FDR. fold_change, PDUI, Num_least ...)
    """
    parser = argparse.ArgumentParser(prog='DaPars', description='Determines the usage of proximal polyA usage')
    parser.add_argument("-c", "--control_alignments", nargs="+", required=True,
                        help="Alignment files in BAM format from control condition")
    parser.add_argument("-t", "--treatment_alignments", nargs="+", required=True,
                        help="Alignment files in BAM format from treatment condition")
    parser.add_argument("-u", "--utr_bed_file", required=True, type=file,
                        help="Bed file describing longest 3UTR positions")
    parser.add_argument("-o", "--output_file", required=True, type=argparse.FileType('w'),
                        help="file containing output")
    parser.add_argument("-cpu", required=False, type=int, default=1,
                        help="Number of CPU cores to use.")
    parser.add_argument("-s", "--search_start", required=False, type=int, default=50,
                        help="Start search for breakpoint n nucleotides downstream of UTR start")
    parser.add_argument("-ct", "--coverage_threshold", required=False, type=float, default=20,
                        help="minimum coverage in each aligment to be considered for determining breakpoints")
    parser.add_argument("-b", "--breakpoint_bed", required=False, type=argparse.FileType('w'),
                        help="Write bedfile with coordinates of breakpoint positions to supplied path.")
    parser.add_argument("-v", "--version", action='version', version='%(prog)s 0.2.0')
    parser.add_argument("-p", "--plot_path", default=None, required=False, type=directory_path,
                        help="If plot_path is specified will write a coverage plot for every UTR in that directory.")
    parser.add_argument("-html", "--html_file", default=None, required=False, type=argparse.FileType('w'),
                        help="Write an html file to the specified location. Only to be used within a galaxy wrapper")
    return parser.parse_args()


class UtrFinder():
    """
    This seems to be the main caller.
    """

    def __init__(self, args):
        self.control_alignments = [file for file in args.control_alignments]
        self.treatment_alignments = [file for file in args.treatment_alignments]
        self.n_cpus = args.cpu
        self.search_start = args.search_start
        self.coverage_threshold = args.coverage_threshold
        self.plot_path = args.plot_path
        self.html_file = args.html_file
        self.utr = args.utr_bed_file
        self.gtf_fields = filter_utr.get_gtf_fields()
        self.result_file = args.output_file
        self.all_alignments = self.control_alignments + self.treatment_alignments
        self.alignment_names = OrderedDict(( file, os.path.basename(file)) for file in self.all_alignments )
        self.num_samples = len(self.all_alignments)
        self.utr_dict = self.get_utr_dict(0.2)
        self.dump_utr_dict_to_bedfile()
        print "Established dictionary of 3\'UTRs"
        self.coverage_files = self.run_bedtools_coverage()
        self.utr_coverages = self.read_coverage_result()
        print "Established dictionary of 3\'UTR coverages"
        self.coverage_weights = self.get_coverage_weights()
        self.result_tuple = self.get_result_tuple()
        self.result_d = self.calculate_apa_ratios()
        self.write_results()
        if args.breakpoint_bed:
            self.bed_output = args.breakpoint_bed
            self.write_bed()
        if self.plot_path:
            self.write_html()

    def dump_utr_dict_to_bedfile(self):
        w = csv.writer(open("tmp_bedfile.bed", "w"), delimiter="\t")
        for gene, utr in self.utr_dict.iteritems():
            w.writerow([utr["chr"], utr["new_start"]-1, utr["new_end"], gene, ".", utr["strand"]])

    def run_bedtools_coverage(self):
        """
        Use bedtools coverage to generate pileup data for all alignment files for the regions specified in utr_dict.
        """
        cmds = []
        for alignment_file in self.all_alignments:
            cmd = "sort -k1,1 -k2,2n tmp_bedfile.bed | "
            cmd = cmd + "bedtools coverage -d -s -abam {alignment_file} -b tmp_bedfile.bed |" \
                        " cut -f 4,7,8 > coverage_file_{alignment_name}".format(
                alignment_file = alignment_file, alignment_name= self.alignment_names[alignment_file] )
            cmds.append(cmd)
        subprocesses = [subprocess.Popen([cmd], shell=True) for cmd in cmds]
        [p.wait() for p in subprocesses]
        coverage_files = ["gene_position_coverage_{alignment_name}".format(
                alignment_name = self.alignment_names[alignment_file]) for alignment_file in self.all_alignments ]
        return coverage_files

    def read_coverage_result(self):
        """
        Read coverages back in and store as dictionary of numpy arrays
        """
        coverage_dict = { gene: [np.zeros(utr_d["new_end"]+1-utr_d["new_start"]) for name in self.alignment_names.values() ] for gene, utr_d in self.utr_dict.items() }
        for i, alignment_name in enumerate(self.alignment_names.values()):
            with open("coverage_file_{alignment_name}".format(alignment_name = alignment_name)) as coverage_file:
                for line in coverage_file:
                    gene, position, coverage= line.strip().split("\t")
                    coverage_dict[gene][i][int(position)-1] = coverage
        for utr_d in self.utr_dict.values():
            if utr_d["strand"] == "-":
                for i, alignment_name in enumerate(self.alignment_names.values()):
                    coverage_dict[gene][i] = coverage_dict[gene][i][::-1]
        return coverage_dict

    def get_utr_dict(self, shift):
        """
        The utr end is extended by UTR length * shift, to discover novel distal polyA sites.
        Set to 0 to disable.
        """
        utr_dict = OrderedDict()
        for line in self.utr:
            if not line.startswith("#"):
                filter_utr.get_feature_dict( line=line, gtf_fields=self.gtf_fields, utr_dict=utr_dict, feature="UTR" )
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

    def get_coverage_weights(self):
        """
        Return weights for normalizing coverage.
        utr_coverage is still confusing.
        """
        coverage_per_alignment = []
        for utr in self.utr_coverages.itervalues():  # TODO: be smarter about this.
            utr_coverage = []
            for vector in utr:
                utr_coverage.append(np.sum(vector))
            coverage_per_alignment.append(utr_coverage)
        coverages = np.array([ sum(x) for x in zip(*coverage_per_alignment) ])
        coverage_weights = coverages / np.mean(coverages)  # TODO: proabably median is better suited? Or even no normalization!
        return coverage_weights

    def get_result_tuple(self):
        static_desc = ["chr", "start", "end", "strand", "gene", "t_stat", "p_value", "breakpoint",
                       "breakpoint_type", "control_mean_percent", "treatment_mean_percent" ]
        samples_desc = []
        for statistic in ["coverage_long", "coverage_short", "percent_long"]:
            for i, sample in enumerate(self.control_alignments):
                samples_desc.append("control_{i}_{statistic}".format(i=i, statistic = statistic))
            for i, sample in enumerate(self.treatment_alignments):
                samples_desc.append("treatment_{i}_{statistic}".format(i=i, statistic = statistic))
        return namedtuple("result", static_desc + samples_desc)

    def calculate_apa_ratios(self):
        result_d = OrderedDict()
        arg_d = {"result_tuple": self.result_tuple,
                 "coverage_weights":self.coverage_weights,
                 "num_samples":self.num_samples,
                 "num_control":len(self.control_alignments),
                 "num_treatment":len(self.treatment_alignments),
                 "result_d":result_d}
        pool = Pool(self.n_cpus)
        tasks = [ (self.utr_coverages[utr], self.plot_path, utr, utr_d, self.coverage_weights, len(self.control_alignments),
                   len(self.treatment_alignments), self.search_start, self.coverage_threshold) \
                  for utr, utr_d in self.utr_dict.iteritems() ]
        processed_tasks = [ pool.apply_async(calculate_all_utr, t) for t in tasks]
        result_list = [res.get() for res in processed_tasks]
        #result_list = [calculate_all_utr(*t) for t in tasks]  # uncomment for easier debugging
        for res_control, res_treatment in result_list:
            if not res_control:
                continue
            for i, result in enumerate(res_control):
                if isinstance(result, dict):
                    t = self.result_tuple(**result)
                    result_d[result["gene"]+"_bp_control_{i}".format(i=i)] = t
            for i, result in enumerate(res_treatment):
                if isinstance(result, dict):
                    t = self.result_tuple(**result)
                    result_d[result["gene"]+"_bp_treatment_{i}".format(i=i)] = t
        return result_d

    def write_results(self):
        w = csv.writer(self.result_file, delimiter='\t')
        header = list(self.result_tuple._fields)
        header[0] = "#chr"
        w.writerow(header)    # field header
        w.writerows( self.result_d.values())

    def write_html(self):
        output_lines = [(gene_str_to_link(result.gene), result.breakpoint, result.breakpoint_type, result.p_value ) for result in self.result_d.itervalues()]
        if self.html_file:
            self.html_file.write(tabulate(output_lines, headers=["gene", "breakpoint", "breakpoint_type", "p_value"], tablefmt="html"))
        else:
            with open(os.path.join(self.plot_path, "index.html"), "w") as html_file:
                html_file.write(tabulate(output_lines, headers=["gene", "breakpoint", "breakpoint_type", "p_value"], tablefmt="html"))

    def write_bed(self):
        w = csv.writer(self.bed_output, delimiter='\t')
        bed = [(result.chr, result.breakpoint, int(result.breakpoint)+1, result.gene+"_"+result.breakpoint_type, 0, result.strand) for result in self.result_d.itervalues()]
        w.writerows(bed)


def calculate_all_utr(utr_coverage, plot_path, utr, utr_d, coverage_weights, num_control, num_treatment, search_start, coverage_threshold):
    if utr_d["strand"] == "+":
        is_reverse = False
    else:
        is_reverse = True
    control_breakpoints, control_abundances, treatment_breakpoints, treatment_abundances  = \
        optimize_breakpoint(plot_path, utr, utr_coverage, utr_d["new_start"], utr_d["new_end"], coverage_weights, search_start, coverage_threshold, num_control)
    res_control = breakpoints_to_result(utr, utr_d, control_breakpoints, "control_breakpoint", control_abundances, is_reverse,
                             num_control, num_treatment)
    res_treatment = breakpoints_to_result(utr, utr_d, treatment_breakpoints, "treatment_breakpoint", treatment_abundances, is_reverse,
                             num_control, num_treatment)
    return res_control, res_treatment


def breakpoints_to_result(utr, utr_d, breakpoints, breakpoint_type,
                         abundances, is_reverse, num_control, num_treatment):
    """
    Takes in a result dictionary res and fills the necessary fields
    """
    if not breakpoints:
        return False
    result = []
    for breakpoint, abundance in zip(breakpoints, abundances):
        res = {}
        long_coverage_vector = abundance[0]
        short_coverage_vector = abundance[1]
        percentage_long = long_coverage_vector/(long_coverage_vector+short_coverage_vector)
        for i in range(num_control):
            res["control_{i}_coverage_long".format(i=i)] = float(long_coverage_vector[i])
            res["control_{i}_coverage_short".format(i=i)] = float(short_coverage_vector[i])
            res["control_{i}_percent_long".format(i=i)] = percentage_long[i]
        for k in range(num_treatment):
            i = k + num_control
            res["treatment_{i}_coverage_long".format(i=k)] = float(long_coverage_vector[i])
            res["treatment_{i}_coverage_short".format(i=k)] = float(short_coverage_vector[i])
            res["treatment_{i}_percent_long".format(i=k)] = percentage_long[i]
        res["t_stat"], res["p_value"] = stat_test(percentage_long[:num_control], percentage_long[num_control:])
        control_mean_percent = np.mean(np.array(percentage_long[:num_control]))
        treatment_mean_percent = np.mean(np.array(percentage_long[num_control:]))
        res["chr"] = utr_d["chr"]
        res["start"] = utr_d["start"]
        res["end"] = utr_d["end"]
        res["strand"] = utr_d["strand"]
        if is_reverse:
            breakpoint = utr_d["new_end"] - breakpoint
        else:
            breakpoint = utr_d["new_start"] + breakpoint
        res["breakpoint"] = breakpoint
        res["breakpoint_type"] = breakpoint_type
        res["control_mean_percent"] = control_mean_percent
        res["treatment_mean_percent"] = treatment_mean_percent
        res["gene"] = utr
        result.append(res)
    return result


def optimize_breakpoint(plot_path, utr, utr_coverage, UTR_start, UTR_end, coverage_weigths, search_start, coverage_threshold, num_control):
    """
    We are searching for a point within the UTR that minimizes the mean squared error, if the coverage vector was divided
    at that point. utr_coverage is a list with items corresponding to numpy arrays of coverage for a sample.
    """
    num_samples = len(utr_coverage)
    normalized_utr_coverage = np.array(utr_coverage)/np.expand_dims(coverage_weigths, axis=1)
    start_coverage = [np.mean(coverage[0:99]) for coverage in utr_coverage]  # filters threshold on mean coverage over first 100 nt
    is_above_threshold = sum(np.array(start_coverage) >= coverage_threshold) >= num_samples  # This filters on the raw threshold. Why?
    is_above_length = UTR_end - UTR_start >= 150
    if (is_above_threshold) and (is_above_length):
        search_end = UTR_end - UTR_start
        breakpoints = range(search_start, search_end + 1)
        mse_list = [ estimate_mse(normalized_utr_coverage, bp, num_samples, num_control) for bp in breakpoints ]
        mse_list = [mse_list[0] for i in xrange(search_start)] + mse_list
        if plot_path:
            plot_coverage_breakpoint(plot_path, utr, mse_list, normalized_utr_coverage, num_control)
        if len(mse_list) > 0:
            return mse_to_breakpoint(mse_list, normalized_utr_coverage, num_samples)
    return False, False, False, False


def plot_coverage_breakpoint(plot_path, utr, mse_list, normalized_utr_coverage, num_control):
    """

    """
    fig = plt.figure(figsize=(8, 8))
    gs = gridspec.GridSpec(2, 1)
    ax1 = plt.subplot(gs[0, :])
    ax2 = plt.subplot(gs[1, :])
    ax1.set_title("mean-squared error plot")
    ax1.set_ylabel("mean-squared error")
    ax1.set_xlabel("nt after UTR start")
    ax2.set_title("coverage plot")
    ax2.set_xlabel("nt after UTR start")
    ax2.set_ylabel("normalized nucleotide coverage")
    mse_control = [ condition[0] for condition in mse_list]
    mse_treatment = [ condition[1] for condition in mse_list]
    minima_control = get_minima(np.array(mse_control))
    minima_treatment = get_minima(np.array(mse_treatment))
    control = normalized_utr_coverage[:num_control]
    treatment = normalized_utr_coverage[num_control:]
    ax1.plot(mse_control, "b-")
    ax1.plot(mse_treatment, "r-")
    [ax2.plot(cov, "b-") for cov in control]
    [ax2.plot(cov, "r-") for cov in treatment]
    [ax2.axvline(val, color="b", alpha=0.25) for val in minima_control]
    ax2.axvline(mse_control.index(min(mse_control)), color="b", alpha=1)
    [ax2.axvline(val, color="r", alpha=0.25) for val in minima_treatment]
    ax2.axvline(mse_treatment.index(min(mse_treatment)), color="r", alpha=1)
    fig.add_subplot(ax1)
    fig.add_subplot(ax2)
    gs.tight_layout(fig)
    fig.savefig(os.path.join(plot_path, "{utr}.svg".format(utr=utr)))


def mse_to_breakpoint(mse_list, normalized_utr_coverage, num_samples):
    """
    Take in mse_list with control and treatment mse and return breakpoint and utr abundance for all local minima
    in mse_list
    """
    mse_control = np.array([mse[0] for mse in mse_list])
    mse_treatment = np.array([mse[1] for mse in mse_list])
    control_breakpoints = list(get_minima(mse_control))
    treatment_breakpoints = list(get_minima(mse_treatment))
    control_abundances = [estimate_abundance(normalized_utr_coverage, bp, num_samples) for bp in control_breakpoints]
    treatment_abundances = [estimate_abundance(normalized_utr_coverage, bp, num_samples) for bp in treatment_breakpoints]
    return control_breakpoints, control_abundances, treatment_breakpoints, treatment_abundances

def get_minima(a):
    """
    get minima for numpy array a
    """
    return np.where(np.r_[True, a[1:] < a[:-1]] & np.r_[a[:-1] < a[1:], True])[0]+1

def estimate_mse(cov, bp, num_samples, num_control):
    """
    get abundance of long utr vs short utr with breakpoint specifying the position of long and short utr.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        long_utr_vector = cov[:num_samples, bp:]
        short_utr_vector = cov[:num_samples, 0:bp]
        mean_long_utr = np.mean(long_utr_vector, 1)
        mean_short_utr = np.mean(short_utr_vector, 1)
        square_mean_centered_short_utr_vector = (short_utr_vector[:num_samples] - mean_short_utr[:, np.newaxis] )**2
        square_mean_centered_long_utr_vector = (long_utr_vector[:num_samples] - mean_long_utr[:, np.newaxis])**2
        mse_control = np.mean(np.append(square_mean_centered_long_utr_vector[:num_control], square_mean_centered_short_utr_vector[:num_control]))
        mse_treatment = np.mean(np.append(square_mean_centered_long_utr_vector[num_control:], square_mean_centered_short_utr_vector[num_control:]))
        return mse_control, mse_treatment


def estimate_abundance(cov, bp, num_samples):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        long_utr_vector = cov[:num_samples, bp:]
        short_utr_vector = cov[:num_samples, 0:bp]
        mean_long_utr = np.mean(long_utr_vector, 1)
        mean_short_utr = np.mean(short_utr_vector, 1)
        return mean_long_utr, mean_short_utr

def stat_test(a,b):
    return stats.ttest_ind(a,b)

def gene_str_to_link(str):
    return "<a href=\"{str}.svg\" type=\"image/svg+xml\" target=\"_blank\">{str}</a>".format(str=str)

if __name__ == '__main__':
    args = parse_args()
    find_utr = UtrFinder(args)

