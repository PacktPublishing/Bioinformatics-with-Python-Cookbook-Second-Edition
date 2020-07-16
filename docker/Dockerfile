FROM continuumio/anaconda3:5.2.0
MAINTAINER Tiago Antao <tiagoantao@gmail.com>
#ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get upgrade -y && apt-get install -y git wget build-essential unzip graphviz libgraphviz-dev pkg-config swig libx11-dev libgsl0-dev libopenblas-dev liblapacke-dev
#RUN apt-get install -y samtools mafft muscle raxml tabix

# R
#RUN apt-get install -y r-bioc-biobase

RUN conda update -n base -c defaults conda --force && conda update -y --all && conda install -y -c conda-forge statsmodels pip rpy2 r-ggplot2 r-gridextra pygraphviz seaborn pexpect networkx reportlab tzlocal simupop biopython
RUN conda install -y -c bioconda gffutils pyvcf dendropy genepop plink trimal eigensoft pysam
RUN pip install pygenomics ipython
EXPOSE 9875

RUN git clone https://github.com/PacktPublishing/Bioinformatics-with-Python-Cookbook-Second-Edition.git
WORKDIR /Bioinformatics-with-Python-Cookbook-Second-Edition

RUN echo setterm -foreground magenta >> /etc/bash.bashrc
CMD jupyter-notebook --ip=0.0.0.0 --no-browser --allow-root --port=9875
