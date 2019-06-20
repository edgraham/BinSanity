FROM continuumio/miniconda2:latest
MAINTAINER Shengwei Hou, housw2010@gmail.com

# update
RUN apt-get -qq update && \
    apt-get install -y --no-install-recommends gcc make build-essential libtool zlib1g-dev ncbi-blast+ hmmer prodigal mummer && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# install pplacer
RUN conda install --yes -c bioconda pplacer subread bedtools bowtie2 samtools && \
    conda clean -ya

# install python packages
RUN pip install numpy scipy pandas scikit-learn matplotlib biopython pysam dendropy BinSanity

# install checkM
ENV checkM_DIR /checkM
RUN mkdir -p $checkM_DIR
WORKDIR $checkM_DIR
RUN apt-get install -y --no-install-recommends git && \
    git clone https://github.com/Ecogenomics/CheckM.git && \
    cd CheckM && python setup.py install && \
    apt-get remove -y --auto-remove git && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# configure checkM database
ENV checkM_db /db/CheckM
RUN mkdir -p $checkM_db
RUN echo -e "\n/db/CheckM\n" | checkm data setRoot
WORKDIR $checkM_db
RUN apt-get install -y --no-install-recommends wget && \
    wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz && \
    tar zxfv checkm_data_2015_01_16.tar.gz && rm checkm_data_2015_01_16.tar.gz && \
    apt-get remove -y --auto-remove wget && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# update PATH
ENV PATH $PATH:/opt/conda/bin

# Entry
WORKDIR /mnt
CMD [ "/bin/bash" ]

