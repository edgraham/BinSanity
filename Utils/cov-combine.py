#! /usr/bin/env python

import argparse, csv, glob, time
import numpy as np



def get_contigs(c):
    all_coverage = glob.glob("*%s" % (c))
    all_contigs =[]
    for f in all_coverage:
        full = list(csv.reader(open(str(f), 'rb'), delimiter='\t'))
        del full[0]
        cov = []
        contig = []
        for lists in full:
            contig.append(lists[0])
            cov.append(lists[1])
        for name in contig:
            if name not in all_contigs:
                all_contigs.append(name)
    return all_contigs
    

def get_coverage(a,b,c):
    coverage_log_all = {}
    all_ = list(a)
    all_coverage = glob.glob("*%s" % (c))
    for f in all_coverage:
        z = list(csv.reader(open(str(f), 'rb'), delimiter='\t'))
        del z[0]
        cov = []
        coverage_log = []
        contig = []
        for lists in z:
            contig.append(lists[0])
            cov.append(lists[1])
        for n in cov:
            log = np.log10((float(n))+1)
            coverage_log.append(str(log))
    
        data = zip(contig,coverage_log)
        data = dict(data)
        for name in all_:
            if name not in contig:
                if name in coverage_log_all:
                    coverage_log_all[name].append('0')
                elif name not in coverage_log_all:
                    coverage_log_all.setdefault(name, [])
                    coverage_log_all[name].append('0')
        for key in data.viewkeys():
            if key in all_: 
                if key in coverage_log_all:
                    coverage_log_all[key].append(data[key])
                elif key not in coverage_log_all:
                    coverage_log_all.setdefault(key, [])
                    coverage_log_all[key].append(data[key])

                
    out = open(str(b), 'w')
    for key, value in coverage_log_all.iteritems():
        out.write('%s\t%s\n' % (key,'\t'.join(value)))
    out.close()
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='cov-combine.py', usage='%(prog)s -o Output File')
    parser.add_argument("-o", dest="inputoutput", help="Specify the output file")
    parser.add_argument("-c", dest="inputCoverage", default="coverage", help="Specify suffix linking coverage profiles")
    args = parser.parse_args()
    if args.inputoutput is None:
        parser.error('-o output fasta file needed')
    if args.inputCoverage is None:
        parser.error('-c suffic linking coverage profiles needed')
    else:
        time = time.time()
        print"""
        ---------------------------------------------------------------------
        """

        get_coverage(get_contigs(args.inputCoverage),args.inputoutput,args.inputCoverage)
        print ("""
        --------------------------------------------------------------------
                 Finished combined coverage profiles in %s seconds
        ____________________________________________________________________""" % (time.time() - time))