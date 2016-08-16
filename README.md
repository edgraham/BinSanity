# BinSanity v0.1.1

Program implements Affinity Propagation to cluster contigs into putative genomes. BinSanity uses contig coverage as an input, while BinSanity-refine incorporates tetranucleotide frequencies, GC content, and an optional input of coverage. All relevant scripts to produce inputs for BinSanity are provided here.

##  BinSanity ##
###Dependencies###
Versions used at time of last update to script are provided in parenthesis.

* [Numpy](http://www.numpy.org/) (v1.11.1)
* [SciKit](http://scikit-learn.org/stable/install.html) (v0.17.1)
* [Biopython](http://biopython.org/wiki/Download) (v1.66)
* [BedTools] (http://bedtools.readthedocs.io/en/latest/content/installation.html)(v2.17.0)

###Input Files###
* Fasta file
* Combined Coverage Profile:

$contig  coverage-1 coverage-2 coverage-3
con-1   121        89         95
con-2   14         29         21
.......

* 


##Issues##

If an issue arises in the process of utilizing BinSanity please provide please create an issue and we will address is as soon as possible. To expedite a response please provide any associated error messages. 

##Citation##
In Prep:
Elaina Graham, John Heidelberg, & Benjamin Tully. 2016. BinSanity: Unsupervised Clustering of Environmental Microbial Assemblies Using Coverage and Affinity Propagation. 
