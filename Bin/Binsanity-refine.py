#! /usr/bin/env python

import collections
import os, sys, time, argparse
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from sklearn.cluster import AffinityPropagation
from Bio.SeqUtils import GC
import csv
##########################Get GC Count#######################################
def GC_Fasta_file(b,name):
    """Makes a tab-delimeted output indicating the GC count associated with each contig"""
    GC_dict = {}
    input_file = open(os.path.join(b,name), 'r') 
    output_file = open('GC_count','w')
    output_file.write("%s\t%s\n" % ('Gene','GC'))
    for cur_record in SeqIO.parse(input_file, "fasta") :
        gene_name = cur_record.id 
        GC_percent = float(GC(cur_record.seq))
        GC_dict.setdefault(gene_name,[])
        GC_dict[gene_name].append(GC_percent) 
        output_line = '%s\t%i\n' % \
        (gene_name, GC_percent) 
        output_file.write(output_line)

        
    output_file.close() 
    input_file.close()
###########################Get tetramer frequencies#####################################

def kmer_list(dna, k):
    """Makes list of k-mers based on input of k and the dna sequence"""
    result = []
    dna = dna.upper()
    dna_edit = Seq(str(dna))
    reverse_complement = (dna_edit).reverse_complement()
    for x in range(len(dna)+1-k):
        result.append(dna[x:x+k])
    for x in range(len(reverse_complement)+1-k):
        result.append(reverse_complement[x:x+k])
    result= [s for s in result if not s.strip('AGTC')]
    return result


def kmer_counts(kmer,b,name):
    """Builds dictonary with each contig and its respective tetramer counts"""
    kmer_dict = {}
    for record in SeqIO.parse(os.path.join(b,name), "fasta"):
       id_=record.id
       seq = record.seq
       length = len(seq)
       tetra = kmer_list(seq,kmer)           
       c = collections.Counter(tetra)
       c = dict(c)
       val = list(c.values())
       val_edit = []
       for freq in val:
           freq= (float(freq)/float(length))
           val_edit.append(freq)
    
       keys = []
       for key in c:
   	    keys.append(str(key))
       c_edit = dict(zip(keys,val_edit))
       kmer_dict.setdefault(str(id_),[])
       kmer_dict[str(id_)].append(c_edit)
    return kmer_dict
   	                
def output(a):
    """Builds tab-delimeted file based on dictionary of tetramers"""
    topleft = 'contigs' 
    headers = sorted(set(key
                        for row in a.values()
                        for key in row[0]))
    writer = csv.writer(open('tetramer-frequencies','wb'), delimiter='\t')
    writer.writerow([topleft] + headers)
    for key in a:
        row = [key]
        for header in headers:
            row.append(a[key][0].get(header, 0))
        writer.writerow(row)
#########################Combine tetramer frequencies and GC counts in tab delimeted format and normalize################################
def get_contigs(c):
    """Pulls GC file and Tetramer file to create a merged tab-delimeted profile"""
    file_names = []
    contig_names = []
    for filename in os.listdir('.'):
        if str(filename).startswith('GC_count'):
            file_names.append(str(filename))
        elif str(filename).startswith('tetramer-frequencies'):
            file_names.append(str(filename))
        elif str(filename).endswith(c):
            file_names.append(str(filename))
    for filename in file_names:
        if filename.startswith('GC_count') or filename.startswith('tetramer-frequencies'):
            full = list(csv.reader(open(str(filename),'rb'), delimiter='\t'))
            del full[0]
            contig = []
            for lists in full:
                contig.append(lists[0])
            for name in contig: 
                if name not in contig_names:
                    contig_names.append(name)
    contig_dict = {}
    for filename in file_names:
        if filename.startswith('GC_count') or filename.startswith('tetramer-frequencies'):
            full = list(csv.reader(open(str(filename),'rb'), delimiter='\t'))
            del full[0]
            for value in list(full[1]):
                if list(full[1]).index(value) > 0:
                    for lists in full:
                        if lists[0] in contig_names:
                            if lists[0] in contig_dict:
                                contig_dict[lists[0]].append(str(np.log10(float(lists[list(full[1]).index(value)])+1)))
                            elif lists[0] not in full:
                                contig_dict.setdefault(lists[0],[])
                                contig_dict[lists[0]].append(str(np.log10(float(lists[list(full[1]).index(value)])+1)))
    full = list(csv.reader(open(str(c),'rb'), delimiter='\t'))
    for value in list(full[1]):
        if list(full[1]).index(value) > 0:
            for lists in full:
                if lists[0] in contig_names:
                    if lists[0] in contig_dict:
                        contig_dict[lists[0]].append(str(float(lists[list(full[1]).index(value)])))
                    elif lists[0] not in full:
                        contig_dict.setdefault(lists[0],[])
                        contig_dict[lists[0]].append(str(float(lists[list(full[1]).index(value)])))                
        out = open("tetra-GC-out", "w")
        for key,value in contig_dict.iteritems():
            out.write('%s\t%s\n' % (key,'\t'.join(value)))
        out.close()

##########################Run AP for refinement##############################    
def create_fasta_dict(fastafile):
        """makes dictionary using fasta files to be binned"""
	fasta_dict = {}
	for record in SeqIO.parse(fastafile, "fasta"):
		fasta_dict[record.id] = record.seq
	return fasta_dict


def get_cov_data(h):
    """Makes a dictionary of coverage values for affinity propagation"""
    all_cov_data = {}
    #for line in open(str(sys.argv[1]), "r"):
    for line in open(str(h), "r"):
   	line = line.rstrip()
	cov_data = line.split()
	all_cov_data[cov_data[0]] = cov_data
    return all_cov_data


def cov_array(a,b,name,size):
    """Computes a coverage array based on the contigs in files and the contig names associated with the coverage file"""
    names = []
    c = 0
    cov_array = []
    unused = []
    for record in SeqIO.parse(os.path.join(b,name), "fasta"):
        if len(record.seq) >= int(size):
  	    if record.id in a.keys():
 	         data = a[record.id]
	         names.append(data[0])	
	         data.remove(data[0])
	         line = " ".join(data)
	         if c == 1:
                    temparray = np.fromstring(line, dtype=float, sep=' ')
                    cov_array = np.vstack((cov_array, temparray))
                 if c == 0:
                    cov_array = np.fromstring(line, dtype=float, sep=' ')
                    c += 1
        else:
            unused.append(str(record.id))   
    print cov_array.shape
    return cov_array, names, unused			     

def affinity_propagation(array,names,unused,d,e,f,g,b,name):	
            """Uses affinity propagation to make putative bins"""		     
            apclust = AffinityPropagation(damping=float(d), max_iter=int(e), convergence_iter=int(f), copy=True, preference=float(g), affinity='euclidean', verbose=False).fit_predict(array)
            print
            """
            -------------------------------------------------------
                            --Creating Bins--
            -------------------------------------------------------"""
            outfile_data = {}
            i = 0
            while i < len(names):
                if apclust[i] in outfile_data.keys():
			outfile_data[apclust[i]].append(names[i])
		if apclust[i] not in outfile_data.keys():
			outfile_data[apclust[i]] = [names[i]]
		i += 1
            out_name = name.split(".fna")[0]
            unbinned_file = open(str(out_name)+".unbinned-refine.fna", "w" )
            with open(os.path.join(b,name),"r") as input2_file: 
	  	fasta_dict = create_fasta_dict(input2_file)                
	  	for id_ in list(unused):
	  	    unbinned_file.write(">"+str(id_)+"\n"+str(fasta_dict[id_])+"\n") 
	  	count = 0                               
                for k in outfile_data:
	  	    if len(outfile_data[k]) >= 5:
	  	        output_file = open(str(out_name)+".refined_bin_%s.fna" % (k), "w" )
	  	        for x in outfile_data[k]:
	  	            output_file.write(">"+str(x)+"\n"+str(fasta_dict[x])+"\n")
	  	        output_file.close()
	  	        count = count + 1
	  	    elif len(outfile_data[k]) < 5:
	  	        if any((len(fasta_dict[x])>100000) for x in outfile_data[k]):
	  	            output_file = open(str(out_name)+".refined_bin_%s.fna" % (k), "w" )
	  	            for x in outfile_data[k]:
	  	                output_file.write(">"+str(x)+"\n"+str(fasta_dict[x])+"\n")
	  	            output_file.close()
	  	            count = count +1
	  	        else:
	  	            for x in outfile_data[k]:
	  	                unbinned_file.write(">"+str(x)+"\n"+str(fasta_dict[x])+"\n")
	  	    print "Cluster "+str(k)+": "+str(len(outfile_data[k]))
	  	print ("""
	  	Total Number of Bins: %i""" % count)
            unbinned_file.close()       
            
class Logger(object):
    def __init__(self):
        self.terminal = sys.stdout
        self.log = open("cluster-log.txt", "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  

    def flush(self):
        pass  

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='Binsanity-refine.py', usage='%(prog)s  -f [Path To Contig File] -l [name of file] {optional [-x Contig Size Cut Off] [-plot Y] [-p Preference] [-m Max Iterations] [-v Convergence Iterations] [-d Damping factor]}',description="""Script designed to use Affinity Propagation to split
    metagenomic data into bins using contig coverage values. It takes as input a coverage file and files containing the contigs to be binned, then outputs clusters of contigs in putative bins.""")
    parser.add_argument("-f", dest="inputContigFiles", help="Specify directory containing your contigs")
    parser.add_argument("-p", type=float, dest="preference", default=-100, help="Specify a preference (default is -500) Note: decreasing the preference leads to more lumping, increasing will lead to more splitting")
    parser.add_argument("-m", type=int, dest="maxiter", default=4000, help="Specify a max number of iterations (default is 2000)")
    parser.add_argument("-v", type=int, dest="conviter",default=400, help="Specify the convergence iteration number (default is 200), e.g Number of iterations with no change in the number of estimated clusters that stops the convergence.")
    parser.add_argument("-d",default=0.95, type=float, dest="damp", help="Specify a damping factor between 0.5 and 1, default is 0.9")
    parser.add_argument("-l",dest="name",help="Specify file you want to refine")
    parser.add_argument("-x",dest="ContigSize", type=int, default=1000,help="Specify the contig size cut-off (Default 1000 bp)")
    parser.add_argument("-c",dest="inputCoverage",help="Identify coverage file")
    parser.add_Argument("-t", dest="inputKmer", type=int, default=4, help="Specify k-mer desired for analysis, default is 4")

    args = parser.parse_args()
    if (args.inputContigFiles is None) and not args.name:
        print "Please indicate directory in which file to be refined is contained"
    elif (args.name is None) and not args.inputContigFiles:
        print " Please indicate what fasta, fa, or fna file you want to refine"
    elif args.inputCoverage is None:
        print "Please indicate suffix linking coverage profile"
    else:
        start_time = time.time()
        sys.stdout = Logger()
        print """
        -------------------------------------------------------
                     Calculating GC content
        -------------------------------------------------------"""
        GC_time = time.time()
        GC_Fasta_file(args.inputContigFiles, args.name)
        print "GC content calculated in %s seconds" % (time.time() - GC_time)
        
        print """
        -------------------------------------------------------
                     Calculating tetramer frequencies
        -------------------------------------------------------"""
        kmer_time = time.time()
        output(kmer_counts(args.inputKmer,args.inputContigFiles,args.name))
        print "tertamer frequency calculated in %s seconds" % (time.time()- kmer_time)
        
        print """
        ------------------------------------------------------
          Creating Combined Composition and Coverage Profile
        -------------------------------------------------------"""
        combine_time = time.time()
        print "Combined profile created in %s seconds" % (time.time()- combine_time)
       
        get_contigs(args.inputCoverage)
        print """
        -------------------------------------------------------
                Reclustering on Composition and Coverage
        -------------------------------------------------------"""
        print "Preference: " + str(args.preference)
        print "Maximum Iterations: " + str(args.maxiter)
        print "Convergence Iterations: " + str(args.conviter)
        print "Contig Cut-Off: " + str(args.ContigSize)
        print "Damping Factor: " + str(args.damp) + "\n"
        
        
        val1, val2, val3 = cov_array((get_cov_data('tetra-GC-out')), args.inputContigFiles, args.name,args.ContigSize)

        print """
        
        -------------------------------------------------------
                        Re-clustering contigs
        -------------------------------------------------------""" 
        affinity_propagation(val1, val2,val3, args.damp, args.maxiter, args.conviter, args.preference,args.inputContigFiles,args.name)
  
        print("""
        --------------------------------------------------------
              --- Ammended Bins Computed in %s seconds ---
        --------------------------------------------------------""" % (time.time() - start_time))