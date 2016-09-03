#!/usr/bin/env python
import sys, os, argparse,subprocess, time
from Bio import SeqIO
from signal import signal, SIGPIPE, SIG_DFL
def extract_from_bam(bam_results):
    out_dict = {}
    p = subprocess.Popen(["genomeCoverageBed", "-ibam", bam_results], stdout=subprocess.PIPE)
    outfile, error = p.communicate()
    count = 0
    if p.returncode !=0:
        sys.stderror.write(outfile)
    else:
        fh = outfile.split('\n')[:-1]
        for line in fh:
            cols = line.split()
            count = count +1
            try:
                d = out_dict[cols[0]]
            except KeyError:
                d = {}
                out_dict[cols[0]] = d
            if int(cols[1]) == 0:
                d["percentage_covered"] = 100 - float(cols[4]) * 100.0
            else:
                d["cov_mean"] = d.get("cov_mean", 0) + (int(cols[1]) * float(cols[4]))
            sys.stdout.write("\rNumber of record processed for coverage: %i" % count)
            sys.stdout.flush()
    return out_dict


def create_cov_file(length_dictionary,bedtools_dict,outname):
    print
    """
    ---------------------------------------------------------
                    Building Coverage File
    ---------------------------------------------------------"""
    contigs = ["contig"]
    coverage = ["average_coverage"]
    length = ["length"]
    del bedtools_dict['genome']
    count = 0
    for contig in bedtools_dict.viewkeys():
        contig_dict = bedtools_dict[contig]
        get_length = length_dictionary[contig]
        contig_dict.setdefault("length",str(get_length))
        average_coverage = contig_dict.get('cov_mean')
        if average_coverage is None:
            average_coverage = '0'
        average_coverage=str(average_coverage)
        contigs.append(contig)
        coverage.append(average_coverage)
        length.append(contig_dict.get("length"))
        count = count +1
        sys.stdout.write("\rFinal Contigs Processed: %i" % count)
        sys.stdout.flush()
    y = zip(contigs,coverage,length)
    with open(outname, 'w') as file:
        file.writelines('\t'.join(i) + '\n' for i in y)

def find_contig_length(fastafile):
    len_dict = {}
    count = 0
    for record in SeqIO.parse(fastafile, "fasta"):
        len_dict.setdefault(str(record.id),str(len(record.seq)))
        count = count + 1
        sys.stdout.write("\rNumber of sequences Processed for length: %i" % count)
        sys.stdout.flush()
    return len_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='contig-coverage-bam.py', usage='%(prog)s -b Bam File -r Fasta File -o Output File')
    parser.add_argument("-f", dest="fastafile", help="Fasta File")
    parser.add_argument("-b", dest="bamfiles", help="BAM file")
    parser.add_argument("-o",dest="outfile",help="outfile name (e.g file.coverage)")
    args = parser.parse_args()
    signal(SIGPIPE,SIG_DFL)
    print """
    ---------------------------------------------------------
            Finding Length information for each Contig
    ---------------------------------------------------------
    """
    x =find_contig_length(args.fastafile)
    print """
    ---------------------------------------------------------
     Extracting Coverage Information from provided BAM file
    ---------------------------------------------------------
    """
    y = extract_from_bam(args.bamfiles)
    create_cov_file(x,y,args.outfile)