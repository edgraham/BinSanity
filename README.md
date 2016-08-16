# BinSanity v0.1.1

Program implements Affinity Propagation to cluster contigs into putative genomes. BinSanity uses contig coverage as an input, while BinSanity-refine incorporates tetranucleotide frequencies, GC content, and an optional input of coverage. All relevant scripts to produce inputs for BinSanity are provided here.

## BinSanity ##
###Dependencies###
>Versions used at time of last update to script are provided in parenthesis.

* [Numpy](http://www.numpy.org/) (v1.11.1)
* [SciKit](http://scikit-learn.org/stable/install.html) (v0.17.1)
* [Biopython](http://biopython.org/wiki/Download) (v1.66)
* [BedTools](http://bedtools.readthedocs.io/en/latest/content/installation.html) (v2.17.0)

>Programs used to prepare input files for BinSanity and associated utility scripts include:

* [Bowtie2] (https://sourceforge.net/projects/bowtie-bio/) (v2.2.5)
* [Samtools] (http://www.htslib.org/) (v1.2)

>Program used by us to putatively analyze bin completion and redundancy for refinement:
* [CheckM] (https://github.com/Ecogenomics/CheckM)

###Input Files###
* Fasta file
* Combined Coverage Profile
```
contig  coverage-1 coverage-2 coverage-3
con-1   121        89         95
con-2   14         29         21
.......
```

###Script Usage's###
To generate input files for BinSanity the scripts `contig-coverage-bam.py` and `cov-combine.py` are provided:
* `contig-coverage-bam` generates a `.coverage` file that produces a tab delimited file containing average contig coverage from a `.BAM` file. In our tests we used Bowtie2 to produce a `.SAM` file.  To maintain consistency we used the `.coverage`suffix for all files output from this script.
```
contig-coverage-bam -f [fasta-file] -b [Bam-file] -o [out.coverage] 
```
The standard output:
```
 ---------------------------------------------------------
            Finding Length information for each Contig
 ---------------------------------------------------------
 Number of sequences processed for length: 2343
 
  ---------------------------------------------------------
    Extracting Coverage Information from provided BAM file
  ---------------------------------------------------------
  
  Number of records processed for coverage: 2343
  
  ---------------------------------------------------------
                    Building Coverage File
  ---------------------------------------------------------
  Final Contigs processed: 2343

```
An example of the output file is shown below:
```
$less sample-1.coverage

contig cov length
con-1  250 3049
con-2  215 1203
con-3  123 4032
con-4  110 5021
....
```

* `cov-combine.py` generates a combined coverage profile from the individual coverage files produces by `contig-coverage-bam`
```
cov-combine -c [suffix-linking-coverage-files] -o [output-file]
```
Standard output will read:

```
--------------------------------------------------------------------
        Finished combined coverage profiles in 650 seconds
____________________________________________________________________
```
* `Binsanity` clusters contigs based on the input of a combined coverage profile. Preference, damping factor, contig cut-off, convergence iterations, and maximum iterations can be optionally adjusted using `-p`, `-d`, `-x`, `-v`, and `-m` respectively.
```
Binsanity -f [directory-with-fasta-file] -l [fna,fa,fasta] -c [combined-cov]

-------------------------------------------------------
               Running Binsanity
          ---Computing Coverage Array ---
-------------------------------------------------------
Preference: -3
Maximum iterations: 4000
Convergence Iterations: 400
Conitg Cut-Off: 1000
Damping Factor: 0.95
.......
```
* `Binsanity-refine` is used to refine bins with high redundancy or low completion by reclustering only those contigs from bins indicated and incorporation of tetranucleotide frequenices and GC-content. K-mer, preference, damping factor, contig cut-off, convergence iterations, and maximum iterations can be optionally adjusted using `-t`, `-p`, `-d`, `-x`, `-v`, and `-m` respectively.

```
Binsanity-refine -f [directory-with-fasta-file] -l [fasta-file-identifier] -c [combined_cov_file]
```

##Example Problem##
>The Infant Gut Metagenome collected and curated by [Sharon et al 2013](https://t.co/6h8LmNpxpk) was clustered by us to test BinSanity. To confirm you have BinSanity working we have provided a folder `Example` containing the fasta file (`INFANT-GUT-ASSEMBLY.fa`) containing contigs for the Infant Gut Metagenome provided by [Eren et al. 2015](https://doi.org/10.7717/peerj.1319). All files associated with our BinSanity run are also provided, which includes the combined coverage profile (produced using Bowtie2 v2.2.5 on defaults, `contig-coverage-bam.py`, and `cov-combine.py`.

To run the test use the following command while in the `Example` directory:

```
Binsanity -f . -l .fa -p -10
```
The output should be as follows:
```
        -------------------------------------------------------
                         Running Bin-sanity
                    
                    ---Computing Coverage Array ---
        -------------------------------------------------------
        
Preference: -10.0
Maximum Iterations: 4000
Convergence Iterations: 400
Contig Cut-Off: 1000
Damping Factor: 0.95


        -------------------------------------------------------
                         Running Bin-sanity
                    
                    ---Computing Coverage Array ---
        -------------------------------------------------------
        
Preference: -10.0
Maximum Iterations: 4000
Convergence Iterations: 400
Contig Cut-Off: 1000
Damping Factor: 0.95

(4189, 11)
        
        -------------------------------------------------------
                    ---Clustering Contigs---
        -------------------------------------------------------
        
        

Cluster 0: 5
Cluster 1: 14
Cluster 2: 75
Cluster 3: 105
Cluster 4: 54
Cluster 5: 20
Cluster 6: 34
Cluster 7: 43
Cluster 8: 27
Cluster 9: 105
Cluster 10: 35
Cluster 11: 10
Cluster 12: 39
Cluster 13: 30
Cluster 14: 727
Cluster 15: 256
Cluster 16: 574
Cluster 17: 7
Cluster 18: 620
Cluster 19: 508
Cluster 20: 350
Cluster 21: 551

	  	Total Number of Bins: 22

        
        
        --------------------------------------------------------
              --- Putative Bins Computed in 233.362998962 seconds ---
        --------------------------------------------------------
```

##Issues##

If an issue arises in the process of utilizing BinSanity please provide please create an issue and we will address is as soon as possible. To expedite a response please provide any associated error messages. 

##Citation##
In Prep:
Elaina Graham, John Heidelberg, & Benjamin Tully. 2016. BinSanity: Unsupervised Clustering of Environmental Microbial Assemblies Using Coverage and Affinity Propagation. 


