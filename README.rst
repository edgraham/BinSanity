# BinSanity v0.1.3#

Program implements Affinity Propagation to cluster contigs into putative genomes. BinSanity uses contig coverage as an input, while BinSanity-refine incorporates tetranucleotide frequencies, GC content, and an optional input of coverage. All relevant scripts to produce inputs for BinSanity are provided here.

## BinSanity ##
###Dependencies###
>Versions used at time of last update to script are provided in parenthesis.

* [Numpy](http://www.numpy.org/) (v1.11.1)
* [SciKit](http://scikit-learn.org/stable/install.html) (v0.17.1)
* [Biopython](http://biopython.org/wiki/Download) (v1.66)
* [BedTools](http://bedtools.readthedocs.io/en/latest/content/installation.html) (v2.17.0)
* [Pandas] (http://pandas.pydata.org/) (v0.13.1)

>Programs used to prepare input files for BinSanity and associated utility scripts include:

* [Bowtie2] (https://sourceforge.net/projects/bowtie-bio/) (v2.2.5)
* [Samtools] (http://www.htslib.org/) (v1.2)

>Program used by us to putatively analyze bin completion and redundancy for refinement:
* [CheckM] (https://github.com/Ecogenomics/CheckM)
#Installation#
```
pip install BinSanity
```
#Workflow#

