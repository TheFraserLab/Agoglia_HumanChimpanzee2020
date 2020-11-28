# Agoglia_HumanChimpanzee2020
Code used for data analysis and figure making for the human-chimpanzee cell fusion project/manuscript. For questions regarding code or analysis, please open a new issue. For questions regarding the project as a whole, access to cell lines, etc., please contact the corresponding authors:

Hunter B Fraser (hbfraser@stanford.edu) and Sergiu P Pasca (spasca@stanford.edu)

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

  - karyotype.R: Used to generate Extended Data Figure 2a,b,d
  - karyo.sbactch: How to run this script on Sherlock, potentially useful for the Fraser lab
  - Figure1.R: Analysis and plots for Figure 1e-k and Extended Data Figures 1f, 2c,e, 3a,c,d: hybrid iPS cells
  - Figure2.R: Analysis and plots for Figures 2c-f, 4b and Extended Data Figures 4a-c, 5, 8e: hybrid cortical spheroids
  - Figure2_scKaryo.R: Analysis and plots for Extended Data Figure 4h-j: single cell karyotype analysis
  - Figure2_Monocle.R: Analysis and plots for Figure 2g-i and Extended Data Figure 3h: monocle analysis
  - Figure3.R: Analysis and plots for Figure 3d-e and Extended Data Figures 4f-g, 7a-e, g: allele specific expression in CS
  - Figure3_parental.R: Analysis and plots for Figure 3c,f-l and Extended Data Figure 6d-e, 7f,h, 8b-d: hCS and cCS timecourse data
  - Figure4.R: Analysis and plots for Figure 4a,b,e,f,i,k,l and Extended Data Figures 8e-g, 9a-i,k, 10a-b, d-e: ASE genes and SSTR2
 
