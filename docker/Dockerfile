FROM continuumio/anaconda3:5.2.0
MAINTAINER Tiago Antao <tiagoantao@gmail.com>
#ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get upgrade -y && apt-get install -y git wget build-essential unzip graphviz libgraphviz-dev pkg-config swig libx11-dev libgsl0-dev libopenblas-dev liblapacke-dev
#RUN apt-get install -y samtools mafft muscle raxml tabix

# R
#RUN apt-get install -y r-bioc-biobase

#RUN conda update -n base -c defaults conda --force && conda update -y --all
RUN conda upgrade -n base conda
RUN conda config --add channels r
RUN conda config --add channels bioconda
RUN conda config --add channels r
RUN conda install -y pip ipython python=3.6
RUN conda install -y -c r r-ggplot2 r-gridextra rpy2
RUN conda install -y -c conda-forge statsmodels
RUN conda install -y -c conda-forge seaborn pexpect networkx reportlab tzlocal simupop biopython
RUN apt-get install -y plink1.9
RUN conda install -y -c bioconda gffutils pyvcf dendropy genepop trimal eigensoft pysam
RUN conda install -y pygraphviz pandas
RUN pip install pygenomics
EXPOSE 9875

RUN git clone https://github.com/PacktPublishing/Bioinformatics-with-Python-Cookbook-Second-Edition.git
WORKDIR /Bioinformatics-with-Python-Cookbook-Second-Edition

RUN echo setterm -foreground magenta >> /etc/bash.bashrc
CMD jupyter-notebook --ip=0.0.0.0 --no-browser --allow-root --port=9875
