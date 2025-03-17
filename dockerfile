FROM ubuntu:22.04

################### System requirements ###################
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    parallel \
#    curl \
#    zlib1g-dev \
#    openjdk-11-jre \
    build-essential

RUN mkdir -p ~/miniconda3 && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh && \
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3 && \
    rm ~/miniconda3/miniconda.sh

# Set the runtime environment path for conda
ENV PATH="~/miniconda3/bin:$PATH"

# Initialize conda 
SHELL ["/bin/bash", "--login", "-c"]
RUN conda init bash

# Create a conda environment
RUN conda create -n biotools37 python=3.7 -y
RUN conda create -n biotools38 python=3.8 -y
#    echo "conda activate biotools" >> ~/.bashrc

################## Tools & Requirements ###################
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

RUN conda install -n biotools38 -y \
    fastp \
    spades
#    refseq_masher \
#    ncbi-datasets-cli \

RUN conda install -n biotools37 -y \
    quast \
    perl-bio-tools-run-alignment-tcoffee \
    perl-bioperl=1.7.2 \
    snippy=4.6.0 \
    gubbins

###################### Docker setup #######################
WORKDIR /data

# Put the shell script into docker container and give permissions
COPY pipeline.sh /data/pipeline.sh
RUN chmod +x /data/pipeline.sh

CMD ["conda", "run", "-n", "biotools", "/data/pipeline.sh"]