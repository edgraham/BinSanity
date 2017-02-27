#BinSanity v0.2.5.5#
BinSanity contains a suite a scripts designed to cluster contigs generated from metagenomic assembly into putative genomes.
##Scripts:##
> More detailed descriptions of each script is given under the **Usage** section
* `Binsanity`
  * BinSanity implements Affinity Propagation to cluster contigs into putative genomes using contig coverage as an input * 
`Binsanity-refine`
  * BinSanity-refine incorporates tetranucleotide frequencies, GC%, and optionally incorporates the coverage profile * `Binsanity-wf`
  * Binsanity-wf runs Binsanity and Binsanity-refine sequentially to optimize cluster results * `Binsanity-profile`
  * Binsanity-profile uses [featureCounts](http://bioinf.wehi.edu.au/featureCounts/) to produce the coverage profiles requires in 
Binsanity, Binsanity-refine, and Binsanity-wf * `Binsanity-lc`
  * Binsanity-lc is an experimental script for large metagenomic assemblies where Binsanity and Binsanity-refine become to memory 
intensive. It uses K-means to subset contigs before implementing Binsanity **(*BetaVersion*)** * `checkm_analysis`
  * checkm_analysis uses [CheckM](http://ecogenomics.github.io/CheckM/) to evaluate completion, redundancy, and strain heterogeneity and 
subsets clusters to aid downstream refinement efforts. * `transform-coverage-profile`
  * If Binsanity-profile is not used to generate the coverage profile transform-coverage-profile is provided to transform a raw coverage 
matrix via one of the provided methods.
##Dependencies##
> Versions used at time of last update to script are provided in parenthesis.
* [Numpy](http://www.numpy.org/) (v1.11.1) * [SciKit](http://scikit-learn.org/stable/install.html) (v0.17.1) * 
[Biopython](http://biopython.org/wiki/Download) (v1.66) * [BedTools](http://bedtools.readthedocs.io/en/latest/content/installation.html) 
(v2.17.0) * [Pandas] (http://pandas.pydata.org/) (v0.13.1) * [FeatureCounts](http://bioinf.wehi.edu.au/featureCounts/)(v1.5.0-p2)
>Programs used to prepare input files for BinSanity and associated utility scripts include:
* [Bowtie2] (https://sourceforge.net/projects/bowtie-bio/) (v2.2.5) * [Samtools] (http://www.htslib.org/) (v1.2)
>Program used by us to putatively analyze bin completion and redundancy for refinement:
* [CheckM] (https://github.com/Ecogenomics/CheckM)
<p>

#Installation#
* Usage of BinSanity requires installation of numpy ``` $ pip install numpy ``` * If you want to use the BinSanity workflow 
`Binsanity-wf` you will also need to download [HMMER](http://hmmer.org/) and [pplacer](http://matsen.fhcrc.org/pplacer/) * HMMER can be 
downloaded like this: ``` $ wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2.tar.gz $ tar -zxvf hmmer-3.1b2.tar.gz $ cd 
hmmer-3.1b2 $ ./configure && make && sudo make install $ cd easel && make check && sudo make install ```
> to download pplacer follow the instructions [here](http://matsen.github.io/pplacer/compiling.html) or do the following:
* Go to the [pplacer webpage](http://matsen.fredhutch.org/pplacer/) and click on `latest release`. This will take you to the github page 
and contain source code. Download the file `pplacer-linux-v1.1.alpha19.zip`. In that file are pre-compiled version of `pplacer`, `guppy`, 
and `rppr`. The location of these need to be exported to your path. ``` $ export PATH=/path/to/directory/with/pplacer:$PATH ```
> If you want to use `Binsanity-profile` you will need to download [featureCounts](http://bioinf.wehi.edu.au/featureCounts/)
* To do this download the latest [subread](https://sourceforge.net/projects/subread/files/) package (subread-1.x.x-source.tar.gz) ``` $ 
tar zxvf subread-1.*.*-source.tar.gz $ cd subread-1.*.*-source/src $ make -f Makefile.Linux ; cd ../ ``` * This will produce a directory 
called `bin` in the `subread-1.x.x-source` file containing the executables for `featureCounts`. These should be copied into your path. 
``` $ sudo cp -r bin ``` 
**Install the latest stable version of the BinSanity suite of scripts via pip** 

``` 
$ sudo pip install BinSanity 
```
## FAQ ##
<p> **Why does Binsanity use more memory than other programs like CONCOCT or MetaBat?** <p> Binsanity implements Affinity Propagation 
(AP), an algorithm that identifies exemplars among data points and forms clusters of data points around these exemplars. This is done by 
considering every data point as a potential exemplar and exchanges messages until a good set of clusters emerges. AP is a [deterministic 
algorithm](https://en.wikipedia.org/wiki/Deterministic_algorithm) in which time and memory requirements scale linearly with the number of 
similarities. <p> BinSanity's accuracy is due in part to the biphasic approach in which BinSanity separates coverage and composition 
during clustering, but also relies heavily on the implementation of AP. Unfortunately our attempts to use other clustering algorithms 
that were less computatonally intensive ultimately sacrificed accuracy. We are currently working on methods to solve the current memory 
limitations. <p> **How much memory will BinSanity Require?** <p> On a Dell PowerEdge R920 with 1TB of available RAM and Intel Xeon 2.3GHz 
processors, it took 191 minutes and ~ 54 GB RAM to run 27,643 contigs. Due to the linear increase of memory we have chosen to cap contigs 
at 100,000 by choosing appropiate size cut-offs for use of this machine. Please contact us with any questions regarding this or 
suggestions on the best way to implement BinSanity using whatever computer cluster you have access to. <p> **If Binsanity-lc uses less 
memory then why not implement this all the time?** <p> Binsanity-lc reduces memory complexity by subseting the contigs into groups based 
on roughly clustering contigs using k-means. In some cases Binsanity-lc will produce identical results to Binsanity, but in some cases 
Binsanity-lc will lead to incorrect contig assignment. Unlike Affinity Propagation, K-means clustering requires human input of 
information criteria that dictate the ultimate number of clusters (N). You could estimate this number by using single copy genes to 
estiamte how many genomes you may have in an assembly (such as [here](http://merenlab.org/2015/12/07/predicting-number-of-genomes/)) and 
use this as a guide to initialize clustering. Methods that require a priori identification of cluster numbers (N) in some cases can 
mis-cluster contigs because they can end up forcing a contig to fit in one of N number of bins when a fit may not exists. So in essence 
the most memory efficient route isn't always the best one. The computational intensity of Affinity Propagation may make the method more 
difficult to implement, but ultimately maintains a consistent level of accuracy.
##Script Usage##
First you need to generate input files for Binsanity (e.g the coverage profile). To generate input files for BinSanity the script 
`Binsanity-profile` is provided: * `Binsanity-profile` generates a `.cov` file containing both average contig coverage from a `.BAM` file 
calculated via featureCounts. In our tests we used Bowtie2 to produce a `.SAM` file, and the converted to a `bam` file using SamTools.  
To maintain consistency we used the `.cov` suffix for all files output from this script. The script also takes as input the contig ids to 
be used in generating a profile. This is to reduce the size of the file by excluding contigs that wont be clustered (for example if you 
wanted to cluster only contigs greater than 1000 bp you would have a list of contig ids greater than 1000bp). * There are multiple 
transformation options identified by the flag `--transform`. We recommend the scaled option.
  * scale --> Scaled by multiplying by 100 and log transformed
  * log --> Log transform
  * None --> Raw Coverage Values
  * X5 --> Multiplication by 5
  * X10 --> Multiplication by 10
  * X100 --> multiplication by 100
  * SQR --> Square root<br /> * The `--ids` flag should indicate a file containing the contig ids. This can easily be generated via bash, 
or through our utilites script `get-ids`. This has parameters `-f` to identify directory containing assembly file, `-l` to identify name 
of fasta file,`-o` to identify name for output file containing contig ids, `-x` containign the size cutoff for inclusion in the final 
coverage profile. ``` contig-1 contig-2 contig-3 ..... ``` <p> Other transformations can be useful in cases where there is an extremely 
low range distribution of coverages and when coverage values are low ``` $ Binsanity-profile -i assembly.fa -s directory/to/BAM/files 
--ids ids.txt ``` * Note that it will read in ALL BAM files in the directory indicated and try to incorporate those into the coverage 
profile * This script will output two files. The raw `.cov` file and the transformed `.cov.x100.lognorm` coverage profile. ``` $less 
assmebly.cov.x100.lognorm contig-1 1.2 0.4 contig-2 1.0 0.4 contig-3 1.3 4.2 contig-4 1.1 5.1 .... ```
###Running BinSanity###
* Four scripts are available to run Binsanity: `Binsanity`,`Binsanity-refine`,`Binsanity-wf`, `Binsanity-lc`
  * `Binsanity` runs the coverage based clustering and takes as input a transformed coverage profile.
  * `Binsanity-refine` is used to refine bins using tetranucleotide frequencies, GC content, and Coverage.
  * `Binsanity-wf` is an automated workflow that combines the Binsanity and Binsanity-refine step.
  * `Binsanity-lc` combines both K-means and Affinity Propagation to reduce memory intensity (**BETA VERSION**)
####Running the Binsanity Workflow####
The help menu for `Binsanity-wf` is shown below (Notation is shared for `Binsanity`, `Binsanity-refine`, and `Binsanity-lc`): ``` usage: 
Binsanity-wf [-h] [-c INPUTCOVFILE] [-f INPUTCONTIGFILES]
                    [-p PREFERENCE] [-m MAXITER] [-v CONVITER] [-d DAMP]
                    [-l FASTAFILE] [-x CONTIGSIZE] [-o OUTPUTDIR]
                    [--threads THREADS] [--kmer KMER]
                    [--refine-preference INPUTREFINEDPREF] [--version]
    Binsanity clusters contigs based on coverage and refines these 'bins' using tetramer frequencies and GC content.
    ----------------------------------------------------------------------------------------------------------------
    Binsanity-wf is a workflow script that runs Binsanity and Binsanity-refine sequenctially.
    The following is including in the workflow:
    STEP 1. Run Binsanity
    STEP 2: Run CheckM to estimate completeness for Refinement
    STEP 3: Run Binsanity-refine
    STEP 4: Create Final BinSanity Clusters
  
  Arguments:
  -h, --help show this help message and exit
  -c INPUTCOVFILE Specify a Coverage File
  -f INPUTCONTIGFILES Specify directory containing your contigs
  -p PREFERENCE Specify a preference (default is -3)
                            Note: decreasing the preference leads to more lumping,
                            increasing will lead to more splitting. If your range
                            of coverages are low you will want to decrease the
                            preference, if you have 10 or less replicates increasing
                            the preference could benefit you.
  -m MAXITER Specify a max number of iterations (default is 4000)
  -v CONVITER Specify the convergence iteration number (default is 400)
                            e.g Number of iterations with no change in the number
                            of estimated clusters that stops the convergence.
  -d DAMP Specify a damping factor between 0.5 and 1, default is 0.95
  -l FASTAFILE Specify the fasta file containing contigs you want to cluster
  -x CONTIGSIZE Specify the contig size cut-off (Default 1000 bp)
  -o OUTPUTDIR Give a name to the directory BinSanity results will be output in [Default is 'BINSANITY-RESULTS']
  --threads THREADS Indicate how many threads you want dedicated to the subprocess CheckM
  --kmer KMER Indicate a number for the kmer calculation, the default is 4
  --refine-preference INPUTREFINEDPREF
                        Specify a preference for refinement. The Default is -25
  --version show program's version number and exit ``` Using the Infant Gut Metagenome (available in the Examples file) the command would 
be as follows: ``` $ Binsanity-wf -f /path/to/fasta -l igm.fa -c Infant_gut_assembly.cov.x100.lognorm ``` The default preference for the 
initial binning and refinement step are -3 and -25 respectively. In most cases these work well, but modifications can be made relative to 
the sample type and expected level of diversity. <p> This workflow will do the following: <p> * run`Binsanity` solely with coverage. * 
run CheckM to determine completion and redundancy (Values used to make the distinction between completion and redundant are given below, 
we provide the uncoupled scripts so that the user can optionally use their own methods to discern completion and redundancy) * run 
`Binsanity-refine` to recluster redundant bins and refine bins with low completion.
#### Setting Completion and Redundancy Estimates for refinement ####
For the purposes of our analysis we used CheckM as a means of generally indicating high and low redundancy bins to use the refinement 
script on. To speed up this process a script was written `checkm_analysis` to parse the output of checkM qa and separate Binsanity 
produced bins into categories of high redundancy, low completion, high completion, and strain redundacy.<p> Currently the thresholds 
written into the script place bins into categories using the following parameters:<p> * High completion: greater than 95% complete with < 
10% redundancy, greater than 80% with <5% redundancy, or > 50% with < 2% redundacy * Low completion: less than 50% complete with 
<5%redundancy * Strain redundancy: greater than 90% complete with > 90% strain heterogeneity * High Redundancy: 80% complete with >10% 
redundacy, or 50% complete > 5% redundacy <p> The program is written in to the script `Binsanity-wf`, but can also be called as a stand 
alone script available in the Utils. It is run as: <p> `checkm_analysis -checkM [checkm_qa tab delimited output]` <p> It should be noted 
that selection of the high and low redundancy values are an arbitrary cut off and the values of generally accepted redundancy, 
completion, and strain heterogeneity are up for debate so it is recommended that if you use the script that you decide what the best cut 
off values are for your purposes. For example if you are looking for members of the [candidate phyla radiation 
(CPR)](doi:10.1038/nature14486) the completion estimates may place the bin below accepted thresholds within this context of 'high 
completion' and 'low completion'. <p> CheckM is also only one means of evaluating bins. This script is provided as a means to make 
refinement using BinSanity slightly simpler by quickly moving bins produced during a first pass of BinSanity into smaller categories for 
further analysis (Note this isn't really necessary if you have a small enough data, but for example in cases where we have produced 100's 
of bins using BinSanity it becomes increasingly more time consuming to manually separate the high and low redundancy bins.)
##Example Problem##
>The Infant Gut Metagenome collected and curated by [Sharon et al. (2013)](http://dx.doi.org/10.1101/gr.142315.112) was clustered by us 
>to test BinSanity. To confirm you have BinSanity working we have provided a folder `Example` containing the fasta file 
>(`INFANT-GUT-ASSEMBLY.fa`) containing contigs for the Infant Gut Metagenome provided by [Eren et al. 
>(2015)](https://doi.org/10.7717/peerj.1319). All files associated with our BinSanity run are also provided, which includes the combined 
>coverage profile (produced using Bowtie2 v2.2.5 on defaults, `contig-coverage-bam.py`, and `cov-combine.py`.
To run the test use the following command using the igm.fa and Infant_gut_assembly.cov.x100.lognorm ``` $ Binsanity -f . -l igm.fa -p -10 
-c Infant_gut_assembly.cov.x100.lognorm ``` The output should be as follows: ```
        -------------------------------------------------------
                         Running Bin-sanity
                    
                    ---Computing Coverage Array ---
        -------------------------------------------------------
        
Preference: -10.0 Maximum Iterations: 4000 Convergence Iterations: 400 Contig Cut-Off: 1000 Damping Factor: 0.95 (4189, 11)
        
        -------------------------------------------------------
                    ---Clustering Contigs---
        ------------------------------------------------------- Cluster 0: 5 Cluster 1: 14 Cluster 2: 75 Cluster 3: 105 Cluster 4: 54 
Cluster 5: 20 Cluster 6: 34 Cluster 7: 43 Cluster 8: 27 Cluster 9: 105 Cluster 10: 35 Cluster 11: 10 Cluster 12: 39 Cluster 13: 30 
Cluster 14: 727 Cluster 15: 256 Cluster 16: 574 Cluster 17: 7 Cluster 18: 620 Cluster 19: 508 Cluster 20: 350 Cluster 21: 551
  Total Number of Bins: 22
        --------------------------------------------------------
              --- Putative Bins Computed in 233.362998962 seconds ---
        -------------------------------------------------------- ```
##Issues##
If an issue arises in the process of utilizing BinSanity please create an issue and we will address is as soon as possible. To expedite a 
response please provide any associated error messages. As this project is actively being improved any comments or suggestions are 
welcome. * [Binsanity Forum](https://groups.google.com/forum/#!forum/binsanity)
##Citation##
Graham, E., Heidelberg, J. & Tully, B. BinSanity: Unsupervised Clustering of Environmental Microbial Assemblies Using Coverage and 
Affinity Propagation. bioRxiv, doi:10.1101/069567 (2016).
