
#BinSanity v0.2.0#

Program implements Affinity Propagation to cluster contigs into putative genomes. BinSanity uses contig coverage as an input, while BinSanity-refine incorporates tetranucleotide frequencies, GC content, and an optional input of coverage. All relevant scripts to produce inputs for BinSanity are provided here.

##BinSanity ##
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
* Install the latest stable version via pip
```
pip install BinSanity
```

##Script Usage's##

First you need to generate input files for Binsanity (e.g the coverage profile).
To generate input files for BinSanity the script `Binsanity-profile` is provided:
* `Binsanity-profile` generates a `.cov` file containing both average contig coverage from a `.BAM` file calculated via multiBamCov in Bedtools. In our tests we used Bowtie2 to produce a `.SAM` file.  To maintain consistency we used the `.cov`suffix for all files output from this script.
*There are multiple transformation options identified by the flag `--transform`. We recommend the Scaled option.
* Scale --> Scaled by multiplying by 100 and log transformed
* log --> Log transform
* None --> Raw Coverage Values
* X5 --> Multiplication by 5 
* X10 --> Multiplication by 10
* SQR --> Square root<br />
<p> Other transformations can be useful in cases where there is an extremely low range distribution of coverages and when coverage values are low

```
Binsanity-profile -o Infant-gut-assembly --contigs contigs_of_interest.txt -i igm.fa -s /path/to/bam/files --transform Scale
```

The standard output:
```
 ---------------------------------------------------------
                       Formating BED File
 ---------------------------------------------------------

  ---------------------------------------------------------
                    Generating Read Counts
  ---------------------------------------------------------

  
  ---------------------------------------------------------
                   Combining all Readcount Files
  ---------------------------------------------------------
   ---------------------------------------------------------
                   Transforming Coverage Profile
  ---------------------------------------------------------
```
* This script will output two files. The raw `Infant-gut-assembly.cov` file and the transformed `Infant-gut-assembly.cov.x100.lognorm` coverage profile.

```
$less Infant-gut-assembly.cov.x100.lognorm

con-1  1.2 0.4
con-2  1.0 0.4
con-3  1.3 4.2
con-4  1.1 5.1
....
```
###Running BinSanity###
* Three scripts are available to run Binsanity: `Binsanity`,`Binsanity-refine`,`Binsanity-wf`
* `Binsanity` runs the coverage based clustering and takes as input a transformed coverage profile.
* `Binsanity-refine` is used to refine bins using tetranucleotide frequencies, GC content, and Coverage.
* `Binsanity-wf` is an automated workflow that combines the Binsanity and Binsanity-refine step. ```

####Running the Binsanity Workflow####
The help menu for `Binsanity-wf` is shown below:
```
usage: Binsanity-wf [-h] [-c INPUTCOVFILE] [-f INPUTCONTIGFILES]
                    [-p PREFERENCE] [-m MAXITER] [-v CONVITER] [-d DAMP]
                    [-l FASTAFILE] [-x CONTIGSIZE] [-o OUTPUTDIR]
                    [--threads THREADS] [--kmer KMER]
                    [--refine-preference INPUTREFINEDPREF] [--version]

Script designed to use Affinity Propagation to split
    metagenomic data into bins using contig coverage values.
    It takes as input a coverage file and files containing
    the contigs to be binned, then outputs clusters of contigs in putative bins.

optional arguments:
  -h, --help            show this help message and exit
  -c INPUTCOVFILE       Specify a Coverage File

  -f INPUTCONTIGFILES   Specify directory containing your contigs

  -p PREFERENCE         Specify a preference (default is -3)
                            Note: decreasing the preference leads to more lumping,
                            increasing will lead to more splitting. If your range
                            of coverages are low you will want to decrease the
                            preference, if you have 10 or less replicates increasing
                            the preference could benefit you.

  -m MAXITER            Specify a max number of iterations (default is 4000)

  -v CONVITER           Specify the convergence iteration number (default is 400)
                            e.g Number of iterations with no change in the number
                            of estimated clusters that stops the convergence.

  -d DAMP               Specify a damping factor between 0.5 and 1, default is 0.95

  -l FASTAFILE          Specify the fasta file containing contigs you want to cluster

  -x CONTIGSIZE         Specify the contig size cut-off (Default 1000 bp)

  -o OUTPUTDIR          Give a name to the directory BinSanity results will be output in [Default is 'BINSANITY-RESULTS']

  --threads THREADS     Indicate how many threads you want dedicated to the subprocess CheckM

  --kmer KMER           Indicate a number for the kmer calculation, the default is 4

  --refine-preference INPUTREFINEDPREF
                        Specify a preference for refinement. The Default is -25

  --version             show program's version number and exit
```
Using the Infant Gut Metagenome (available in the Examples file) the command would be as follows:
```
Binsanity-wf -f /path/to/fasta -l ifm.fa -c Infant_gut_assembly.cov.x100.lognorm 
```

The default preference for the initial binning and refinement step are -3 and -25 respectively. We feel these are the best settings in most cases, but modifications can be made relative to the sample type and expected level of diversity.
<p>
This workflow will do the following:
<p>
* run`Binsanity` solely with coverage. 
* run CheckM to determine completion and redundancy (Values used to make the distinction between completion and redundant are given below, we provide the uncoupled scripts so that the user can optionally use their own methods to discern completion and redundnacy)
* run `Binsanity-refine` to recluster redundant bins and refine bins with low completion.
#####Setting Completion and Redundnacy Estimates for refinement#####

For the purposes of our analysis we used CheckM as a means of generally indicating high and low redundancy bins to use the refinement script on. To speed up this process a script was written `checkm_analysis` to parse the output of checkM qa and separate Binsanity produced bins into categories of high redundancy, low completion, high completion, and strain redundacy.<p>

Currently the thresholds written into the script place bins into categories using the following parameters:<p>
* High completion: greater than 95% complete with < 10% redundancy, greater than 80% with <5% redundancy,  or > 50% with < 2% redundacy
* Low completion: less than 50% complete with > 90% strain heterogeneity
* High Redundancy: 80% complete with >10% redundacy, or 50% complete > 5% redundacy
<p>
The program is written in to the script `Binsanity-wf`, but can also be called as a stand alone script available in the Utils. It is run as: <p>
`checkm_analysis -checkM [checkm_qa tab delimited output]`
<p>
It should be noted that selection of the high and low redundancy values are an arbitrary cut off and the values of generally accepted redundancy, completion, and strain heterogeneity are up for debate so it is recommended that if you use the script that you decide what the best cut off values are for your purposes.<p>
CheckM is also only one means of evaluating bins. This script is provided as a means to make refinement using BinSanity slightly simpler by quickly moving bins produced during a first pass of BinSanity into smaller categories for further analysis (Note this isn't really necessar if you have a small enough data, but for example in cases where we have produced 100's of bins using BinSanity it becomes increasingly more time consuming to manually separate the high and low redundancy bins.)

##Example Problem##
>The Infant Gut Metagenome collected and curated by [Sharon et al. (2013)](http://dx.doi.org/10.1101/gr.142315.112) was clustered by us to test BinSanity. To confirm you have BinSanity working we have provided a folder `Example` containing the fasta file (`INFANT-GUT-ASSEMBLY.fa`) containing contigs for the Infant Gut Metagenome provided by [Eren et al. (2015)](https://doi.org/10.7717/peerj.1319). All files associated with our BinSanity run are also provided, which includes the combined coverage profile (produced using Bowtie2 v2.2.5 on defaults, `contig-coverage-bam.py`, and `cov-combine.py`.

To run the test use the following command using the igm.fa and Infant_gut_assembly.cov.x100.lognorm

```
Binsanity -f . -l igm.fa -p -10
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

If an issue arises in the process of utilizing BinSanity please create an issue and we will address is as soon as possible. To expedite a response please provide any associated error messages. 
As this project is actively being improved any comments or suggestions are welcome.

* [Binsanity Forum](https://groups.google.com/forum/#!forum/binsanity)

##Citation##
Graham, E., Heidelberg, J. & Tully, B. BinSanity: Unsupervised Clustering of Environmental Microbial Assemblies Using Coverage and Affinity Propagation. bioRxiv, doi:10.1101/069567 (2016).


