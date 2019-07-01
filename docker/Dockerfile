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
RUN conda config --add channels r
RUN conda config --add channels conda-forge
RUN conda install --yes biopython=1.70
RUN conda install --yes statsmodels pysam plink gffutils genepop trimal
RUN conda install --yes simuPOP
RUN conda install --yes pip
RUN conda install --yes rpy2

RUN conda install --yes pygraphviz eigensoft
RUN conda install --yes seaborn pexpect pyvcf dendropy networkx reportlab

RUN pip install pygenomics


EXPOSE 9875

RUN git clone https://github.com/PacktPublishing/Bioinformatics-with-Python-Cookbook-Second-Edition.git
WORKDIR /PacktPublising/notebooks

RUN echo setterm -foreground magenta >> /etc/bash.bashrc
CMD jupyter-notebook --ip=0.0.0.0 --no-browser --allow-root --port=9875
