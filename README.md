# CauseVarTE
An integrated pipeline for an end-to-end analysis of mobile element insertions in whole-genome resequencing data


 The goal of the pipeline is to minimize manual intervention in identifying the non-reference mobile element insertions in the human genome and their potential role in diseases. The pipeline identifies insertions using previously published tools (e.g. MELT (Gardener et al. 2017 Genome Research) and TranSurVeyor (Rajaby and Sung 2018 Nucleic Acid Research). The insertion breakpoints identified by both tools are further validated and annotated using AnnotSV (Geoffroy et al. 2018 Bioinformatics). The candidate insertions segregating with diseases and those insertions in the regions of interest are identified for further evaluation.


Usage:

Python3 Scriptname -conf configfile

Before running the script, update the configfile with the path of the tools and other details. 

Please email jainyt@genetics.utah.edu for questions or help 


