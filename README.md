# Agoglia_HumanChimpanzee2020
Code used for data analysis and figure making for the human-chimpanzee cell fusion project/manuscript. For questions regarding code or analysis, please open a new issue. For questions regarding the project as a whole, access to cell lines, etc., please contact the corresponding authors:

Hunter B Fraser hbfraser@stanford.edu
Sergiu P Pasca  spasca@stanford.edu

Raw and processed data for this project is freely available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144825


Mapping
-------

This folder contains all of the code used to map reads, remove duplicates, correct for mapping bias, and quantify allele specific expression in data from human-chimpanzee hybrid cells. This pipeline was also used for parental cell samples for consistency. 

  - Snakefile: the mapping pipeline used for this analysis
  - sbatch_snakemake.sbatch: how to run this pipeline on Sherlock, potentially useful for the Fraser lab
  - find_intersecting_SNPs.py, filter_remapped_reads.py: from the Hornet pipeline (https://github.com/TheFraserLab/Hornet)
  - GetGeneASEbyReads.py: from the ASEr package (https://github.com/TheFraserLab/ASEr)
  - rmdup_for_ase.py: a wrapper for Samtools' rmdup command that randomizes which duplicate read is discarded (from Ryo Kita)
  
  Also required: STAR (https://github.com/alexdobin/STAR), Samtools (http://www.htslib.org/)
  
  
Figures
-------

This folder contains all of the code used to generate the plots presented in Agoglia et al. 

  - karyotype.R: used to generate Extended Data Figure 2a,b,d
  - karyo.sbactch: how to run this script on Sherlock, potentially useful for the Fraser lab
  - Figure1.R: 
  - Figure2.R: 
  - Figure2_Monocle.R: 
  - Figure3.R: 
  - Figure3_parental.R: 
  - Figure4.R: 
  - scKaryo.R: 
 
