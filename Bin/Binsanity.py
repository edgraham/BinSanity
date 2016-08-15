#! /usr/bin/env python
import os, sys, time, argparse
import numpy as np
from Bio import SeqIO
from sklearn.cluster import AffinityPropagation

def create_fasta_dict(fastafile):
        """makes dictionary using fasta files to be binned"""
	fasta_dict = {}
	for record in SeqIO.parse(fastafile, "fasta"):
		fasta_dict[record.id] = record.seq
	return fasta_dict

def get_cov_data(h):
    """Makes a dictionary of coverage values for affinity propagation"""
    all_cov_data = {}
    for line in open(str(h), "r"):
   	line = line.rstrip()
	cov_data = line.split()
	all_cov_data[cov_data[0]] = cov_data
    return all_cov_data

def cov_array(a,b,link,size):
    """Computes a coverage array based on the contigs in files and the contig names associated with the coverage file"""
    count = 0
    for filename in os.listdir(b):
        file_names = []
        if str(filename).endswith(link):
            file_names.append(filename)
   	    names = []
   	    c = 0
   	    cov_array = []
   	    unused = []
   	    for record in SeqIO.parse(os.path.join(b, filename), "fasta"):
   	        count = count + 1
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
            return cov_array, names, file_names, unused			     

def affinity_propagation(array,names,files,d,e,f,g,unused,b):	
            """Uses affinity propagation to make putative bins"""		     
            apclust = AffinityPropagation(damping=float(d), max_iter=int(e), convergence_iter=int(f), copy=True, preference=int(g), affinity='euclidean', verbose=False).fit_predict(array)
            print
            """
            -------------------------------------------------------
                            --Creating Bins--
            -------------------------------------------------------
            
            """
            outfile_data = {}
            i = 0
            while i < len(names):
                if apclust[i] in outfile_data.keys():
			outfile_data[apclust[i]].append(names[i])
		if apclust[i] not in outfile_data.keys():
			outfile_data[apclust[i]] = [names[i]]
		i += 1
            unbinned_file = open("afprop.unbinned.fna", "w" )
            for filename in files:
                out_name = filename.split(".")[0]
            with open(os.path.join(b,filename),"r") as input2_file: 
	  	fasta_dict = create_fasta_dict(input2_file)                
	  	for id_ in list(unused):
	  	    unbinned_file.write(">"+str(id_)+"\n"+str(fasta_dict[id_])+"\n")
	  	count = 0                
                for k in outfile_data:
	  	    if len(outfile_data[k]) >= 5:
	  	        output_file = open(str(out_name)+".bin_%s.fna" % (k), "w" )
	  	        for x in outfile_data[k]:
	  	            output_file.write(">"+str(x)+"\n"+str(fasta_dict[x])+"\n")
	  	        output_file.close()
	  	        count = count + 1
	  	    elif len(outfile_data[k]) < 5:
	  	        if any((len(fasta_dict[x])>50000) for x in outfile_data[k]):
	  	            output_file = open(str(out_name)+".bin_%s.fna" % (k), "w" )
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


###################################################################

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='Binsanity.py', usage='%(prog)s -c [Coverage File] -f [Path To Contig File] -l [Suffix Linking Contig files] {optional [-x Contig Size Cut Off] [-plot Y] [-p Preference] [-m Max Iterations] [-v Convergence Iterations] [-d Damping factor]}',description="""Script designed to use Affinity Propagation to split
    metagenomic data into bins using contig coverage values. It takes as input a coverage file and files containing the contigs to be binned, then outputs clusters of contigs in putative bins.""")
    parser.add_argument("-c", dest="inputCovFile", help="Specify a Coverage File")
    parser.add_argument("-f", dest="inputContigFiles", help="Specify directory containing your contigs")
    parser.add_argument("-p", type=float, dest="preference", default=-3, help="Specify a preference (default is -3) Note: decreasing the preference leads to more lumping, increasing will lead to more splitting. If your range of coverages are low you will want to decrease the preference, if you have 10 or less replicates increasing the preference could benefit you.")
    parser.add_argument("-m", type=int, dest="maxiter", default=4000, help="Specify a max number of iterations (default is 2000)")
    parser.add_argument("-v", type=int, dest="conviter",default=400, help="Specify the convergence iteration number (default is 200), e.g Number of iterations with no change in the number of estimated clusters that stops the convergence.")
    parser.add_argument("-d",default=0.95, type=float, dest="damp", help="Specify a damping factor between 0.5 and 1, default is 0.9")
    parser.add_argument("-l",dest="link", default="fa",help="Specify the suffix linking your fasta contig files, default is fa")
    parser.add_argument("-x",dest="ContigSize", type=int, default=1000,help="Specify the contig size cut-off (Default 1000 bp)")

    args = parser.parse_args()
    if args.inputCovFile is None:
        print "Please indicate -c coverage file"
    elif args.inputContigFiles is None:
        print "Please indicate -f directory containing your contigs"
    elif args.inputCovFile is None:
        if (args.inputContigFiles is None):
            parser.print_help()
    elif args.inputContigFiles and not args.link:
        parser.error('-l Suffix Linking Contig Files Needed')


    else:
        start_time = time.time()
        sys.stdout = Logger()
        print """
        -------------------------------------------------------
                          Running Binsanity
                    
                    ---Computing Coverage Array ---
        -------------------------------------------------------
        """
        print "Preference: " + str(args.preference)
        print "Maximum Iterations: " + str(args.maxiter)
        print "Convergence Iterations: " + str(args.conviter)
        print "Contig Cut-Off: " + str(args.ContigSize)
        print "Damping Factor: " + str(args.damp) + "\n"
        
        
        val1, val2, val3, val4 = cov_array((get_cov_data(args.inputCovFile)), args.inputContigFiles, args.link,args.ContigSize)

        print """
        
        
        -------------------------------------------------------
                      ---Clustering Contigs---
        -------------------------------------------------------
        
        """ 
        affinity_propagation(val1, val2, val3, args.damp, args.maxiter, args.conviter, args.preference,val4,args.inputContigFiles)
  
        print("""
        
        
        --------------------------------------------------------
              --- Putative Bins Computed in %s seconds ---
        --------------------------------------------------------""" % (time.time() - start_time))