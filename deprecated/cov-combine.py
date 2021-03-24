#! /usr/bin/env python

import argparse, csv, glob, time, pandas
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
            log = ((float(n)))
            coverage_log.append(str(log))
    
        data = list(zip(contig,coverage_log))
        data = dict(data)
        for name in all_:
            if name not in contig:
                if name in coverage_log_all:
                    coverage_log_all[name].append('0')
                elif name not in coverage_log_all:
                    coverage_log_all.setdefault(name, [])
                    coverage_log_all[name].append('0')
        for key in data.keys():
            if key in all_: 
                if key in coverage_log_all:
                    coverage_log_all[key].append(data[key])
                elif key not in coverage_log_all:
                    coverage_log_all.setdefault(key, [])
                    coverage_log_all[key].append(data[key])

                
    out = open(str(b), 'w')
    for key, value in coverage_log_all.items():
        out.write('%s\t%s\n' % (key,'\t'.join(value)))
    out.close()
    
def get_log(file_,outfile):
    cov_file = list(csv.reader(open(str(file_), 'rb'), delimiter='\t'))
    for list_ in cov_file:
       for num_ in range(1,len(list_)):
            list_[num_] = np.log10(((float(list_[num_])))+1)
    pd = pandas.DataFrame(cov_file)
    pd.to_csv(outfile,sep="\t",index=False,header=False)
    
def get_multiple(file_,outfile,num):
    cov_file = list(csv.reader(open(str(file_), 'rb'), delimiter='\t'))
    for list_ in cov_file:
        for num_ in range(1,len(list_)):
            list_[num_] = float(list_[num_])*int(num)
    pd = pandas.DataFrame(cov_file)
    pd.to_csv(outfile,sep="\t",index=False,header=False)  

def get_squareroot(file_,outfile):
    cov_file = list(csv.reader(open(str(file_), 'rb'), delimiter='\t'))
    for list_ in cov_file:
       for num_ in range(1,len(list_)):
            list_[num_] = np.sqrt(((float(list_[num_])))+1)
    pd = pandas.DataFrame(cov_file)
    pd.to_csv(outfile,sep="\t",index=False,header=False)
              
if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='cov-combine.py', usage='%(prog)s -o Output File')
    parser.add_argument("-o", dest="inputoutput", help="Specify the output file")
    parser.add_argument("-c", dest="inputCoverage", default="coverage", help="Specify suffix linking coverage profiles")
    parser.add_argument("-t",dest="transform", default = "log", help ="""Indicate what type of data transformation you want in the final file (default is log):\n
    log --> Log transform \n
    None --> Raw Coverage Values \n
    X5 --> Multiplication by 5 \n
    X10 --> Multiplication by 10 \n
    SQR --> Square root \n
    
    We recommend using a log transformation for initial testing. Other transformations can be useful in cases where there is an extremely low range distribution of coverages and when coverage values are low
    """)
    args = parser.parse_args()
    if args.inputoutput is None:
        parser.error('-o output fasta file needed')
    if args.inputCoverage is None:
        parser.error('-c suffic linking coverage profiles needed')
    else:
        start_time = time.time()
        print("""
        ---------------------------------------------------------------------
        Getting Coverage.....
        """)

        get_coverage(get_contigs(args.inputCoverage),args.inputoutput,args.inputCoverage)
        
        x = str(args.inputoutput)
        if args.transform == "log":
            print("Transforming your combined covergae profile")
            get_log(x, x+".lognorm")
        elif args.transform =="X5":
            print("Transforming your combined covergae profile")            
            get_multiple(x,x+".x5", int(5))
        elif args.transform =="x10":
            print("Transforming your combined covergae profile")            
            get_multiple(x,x+".x5", int(10))
        elif args.transform =="SQR":
            print("Transforming your combined covergae profile")            
            get_squareroot(x, x+".sqrt")
            
            
      
        
        print(("""
        --------------------------------------------------------------------
                 Finished combined coverage profiles in %s seconds
        ____________________________________________________________________""" % (time.time() - start_time)))