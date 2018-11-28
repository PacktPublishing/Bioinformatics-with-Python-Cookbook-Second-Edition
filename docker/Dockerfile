FROM continuumio/anaconda3:5.2.0
MAINTAINER Tiago Antao <tiagoantao@gmail.com>

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update
RUN apt-get upgrade -y

RUN apt-get install -y git wget build-essential unzip
#RUN apt-get install -y samtools mafft muscle raxml tabix

#R
#RUN apt-get install -y r-bioc-biobase

#RUN apt-get install -y graphviz libgraphviz-dev pkg-config  # phylo/biopython
#RUN apt-get install -y swig  # simupop
#RUN apt-get install -y libx11-dev
#RUN apt-get install -y libgsl0ldbl
#RUN apt-get install -y libgsl0-dev libopenblas-dev liblapacke-dev

RUN apt-get clean


RUN conda config --add channels bioconda
RUN conda install --yes biopython=1.70
RUN conda install --yes statsmodels pysam plink gffutils genepop trimal
RUN conda install --yes  simuPOP

RUN root/miniconda/bin/pip install rpy2
#RUN root/miniconda/bin/pip install pygraphviz eigensoft
RUN root/miniconda/bin/pip install seaborn
RUN root/miniconda/bin/pip install pyvcf
RUN root/miniconda/bin/pip install dendropy
RUN root/miniconda/bin/pip install pexpect
RUN root/miniconda/bin/pip install reportlab
RUN root/miniconda/bin/pip install networkx

RUN pip install pygenomics


EXPOSE 9875

RUN git clone https://github.com/PacktPublishing/Bioinformatics-with-Python-Cookbook-Second-Edition.git
WORKDIR /PacktPublising/notebooks

RUN echo setterm -foreground magenta >> /etc/bash.bashrc
CMD jupyter-notebook --ip=0.0.0.0 --no-browser --port=9875
