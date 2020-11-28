## Analysis and plots for Figure 3d-e and Extended Data Figures 4f-g, 7a-e, g: allele specific expression in CS
## Everything here was performed in Rstudio on a Macbook Pro with 16G RAM, unless otherwise noted
## Clear cache in Rstudio
rm(list = ls(all.names = TRUE))

## Load required libraries
library(ggplot2)
library(gplots)
library(gridExtra)
library(Rtsne)
library(DESeq2)
library(dplyr)
library(reshape)
library(scales)
library(VennDiagram)
library(Seurat)
library(graphics)
library(SDMTools)
# library(plyr) # will need to load and unload this later, as it interferes with certain dplyr features

peek <- function(data){
  rowmax = min(5, nrow(data))
  colmax = min(5, ncol(data))
  print(data[1:rowmax, 1:colmax])
}

## Directories where the data is / where to deposit files and plots
data_dir = "/Users/rachelagoglia/Desktop/Revisions/ProcessedData"
analysis_dir = "/Users/rachelagoglia/Desktop/Revisions/AnalysisFiles"
plot_dir = "/Users/rachelagoglia/Desktop/Revisions/Plots"

## Load data 
# Gene meta data
geneInfo = read.table(paste(data_dir, "geneInfo.GRCh38.txt", sep = "/"))
names(geneInfo) = c("gene", "length", "chrom")
geneInfo$genesizeKb = as.vector(geneInfo$length/1000)
dim(geneInfo)
head(geneInfo)

# RNAseq data
hyCS = read.table(paste(data_dir, "hyCS.GRCh38.txt", sep = "/"), header = TRUE)
rownames(hyCS) = geneInfo$gene
dim(hyCS)
hyCS_tot = hyCS[, seq(1,75,by=3)]
hyCS_meta = read.table(paste(data_dir, "hyCS.meta.txt", sep = "/"), header = TRUE)
hyCS_meta$Line = as.character(hyCS_meta$Line)
hyCS_meta$Induction = as.character(hyCS_meta$Induction)

# ID genes with mapping bias so that they can be removed from subsequent allele-specific analyses
# Chimp genome gene meta data
geneInfo.pt = read.table(paste(data_dir, "geneInfo.PanTro5.txt", sep = "/"))
names(geneInfo.pt) = c("gene", "length", "chrom")
geneInfo.pt$genesizeKb = as.vector(geneInfo.pt$length/1000)
dim(geneInfo.pt)
head(geneInfo.pt)

hyCS_pt = read.table(paste(data_dir, "hyCS.PanTro5.txt", sep = "/"), header = TRUE)
rownames(hyCS_pt) = geneInfo.pt$gene
dim(hyCS_pt)
hyCS_pt_tot = hyCS_pt[, seq(1,75,by=3)]

# DESeq2 (to ID genes with biased ASE)
hyCS_c = hyCS[,seq(3, length(hyCS), 3)]
hyCS_h = hyCS[,seq(2, length(hyCS), 3)]

hyCS_pt_h = hyCS_pt[,seq(2, length(hyCS_pt), 3)]
hyCS_pt_c = hyCS_pt[,seq(3, length(hyCS_pt), 3)]


deseqASE <- function(counts_matrix, cond_species, cond_line, cond_induction = 0){
  require(DESeq2)
  if(cond_induction == 0){
    cond <- ~ cond_line + cond_species
    cond_reduced <- ~ cond_line 
    coldata <- data.frame(cond_line, cond_species)
  } else {
    cond <- ~ cond_line + cond_induction + cond_species
    cond_reduced <- ~ cond_line + cond_induction
    coldata <- data.frame(cond_line,cond_induction, cond_species) 
  }
  dds <- DESeqDataSetFromMatrix(counts_matrix, coldata, cond)
  dds <- DESeq(dds,test="LRT",betaPrior=FALSE,full=cond,reduced=cond_reduced)
  res <- results(dds)
  res
}

d50_GRCh38_res = deseqASE(cbind(hyCS_h[,1:7], hyCS_c[,1:7]),
                          c(rep("Human", 7), rep("Chimp", 7)),
                          rep(hyCS_meta$Line[1:7], 2),
                          rep(hyCS_meta$Induction[1:7], 2))
d50_PanTro5_res = deseqASE(cbind(hyCS_pt_h[,1:7], hyCS_pt_c[,1:7]),
                           c(rep("Human", 7), rep("Chimp", 7)),
                           rep(hyCS_meta$Line[1:7], 2),
                           rep(hyCS_meta$Induction[1:7], 2))

d100_GRCh38_res = deseqASE(cbind(hyCS_h[,8:16], hyCS_c[,8:16]),
                           c(rep("Human", 9), rep("Chimp", 9)),
                           rep(hyCS_meta$Line[8:16], 2),
                           rep(hyCS_meta$Induction[8:16], 2))  
d100_PanTro5_res =  deseqASE(cbind(hyCS_pt_h[,8:16], hyCS_pt_c[,8:16]),
                             c(rep("Human", 9), rep("Chimp", 9)),
                             rep(hyCS_meta$Line[8:16], 2),
                             rep(hyCS_meta$Induction[8:16], 2)) 

d150_PanTro5_res = deseqASE(cbind(hyCS_pt_h[,17:23], hyCS_pt_c[,17:23]),
                            c(rep("Human", 7), rep("Chimp", 7)),
                            rep(hyCS_meta$Line[17:23], 2),
                            rep(hyCS_meta$Induction[17:23], 2))
d150_GRCh38_res = deseqASE(cbind(hyCS_h[,17:23], hyCS_c[17:23]),
                           c(rep("Human", 7), rep("Chimp", 7)),
                           rep(hyCS_meta$Line[17:23], 2),
                           rep(hyCS_meta$Induction[17:23], 2)) 

# All one induction so leave out that variable
d200_PanTro5_res = deseqASE(cbind(hyCS_pt_h[,24:25], hyCS_pt_c[,24:25]),
                            c(rep("Human", 2), rep("Chimp", 2)),
                            rep(hyCS_meta$Line[24:25], 2)) 
d200_GRCh38_res = deseqASE(cbind(hyCS_h[,24:25], hyCS_c[24:25]),
                           c(rep("Human", 2), rep("Chimp", 2)),
                           rep(hyCS_meta$Line[24:25], 2)) 

allCS_PanTro5_res = deseqASE(cbind(hyCS_pt_h, hyCS_pt_c),
                             c(rep("Human", 25), rep("Chimp", 25)),
                             rep(hyCS_meta$Line, 2),
                             rep(hyCS_meta$Induction, 2))   
allCS_GRCh38_res = deseqASE(cbind(hyCS_h, hyCS_c),
                            c(rep("Human", 25), rep("Chimp", 25)),
                            rep(hyCS_meta$Line, 2),
                            rep(hyCS_meta$Induction, 2)) 


# Call mapping bias from full results because this has the most power
hc_genes = rownames(hyCS)[which(rownames(hyCS) %in% rownames(hyCS_pt))]
h_points = allCS_GRCh38_res[hc_genes, "log2FoldChange"]
c_points = allCS_PanTro5_res[hc_genes, "log2FoldChange"]
all_biased_genes = hc_genes[which(abs(allCS_GRCh38_res[hc_genes, "log2FoldChange"]-allCS_PanTro5_res[hc_genes, "log2FoldChange"])>1)]
all_unbiased_genes = hc_genes[which(abs(allCS_GRCh38_res[hc_genes, "log2FoldChange"]-allCS_PanTro5_res[hc_genes, "log2FoldChange"])<=1)]
length(all_biased_genes)    # 291
length(all_unbiased_genes)  # 14014

write.table(file = paste(analysis_dir, "all_biased_genes.rev.txt", sep = "/"), x = as.data.frame(all_biased_genes), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(file = paste(analysis_dir, "all_unbiased_genes.rev.txt", sep = "/"), x = as.data.frame(all_unbiased_genes), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# Table of biased genes
write.table(data.frame(cbind(all_biased_genes, h_points[which(hc_genes %in% all_biased_genes)], c_points[which(hc_genes %in% all_biased_genes)])), 
            file=paste(analysis_dir, "CS.biased.rev.txt", sep = "/"), 
            sep="\t", quote = FALSE,na = "", row.names = FALSE, col.names = FALSE)


# Human only results (unfiltered)
hyCS.DESeq2.ASE = data.frame(cbind(d50_GRCh38_res[which(!(rownames(d50_GRCh38_res) %in% all_biased_genes)), c("log2FoldChange", "padj")], 
                                   d100_GRCh38_res[which(!(rownames(d100_GRCh38_res) %in% all_biased_genes)), c("log2FoldChange", "padj")], 
                                   d150_GRCh38_res[which(!(rownames(d150_GRCh38_res) %in% all_biased_genes)), c("log2FoldChange", "padj")], 
                                   d200_GRCh38_res[which(!(rownames(d200_GRCh38_res) %in% all_biased_genes)), c("log2FoldChange", "padj")]))

colnames(hyCS.DESeq2.ASE) = c("D50_HumanLFC", "D50_HumanPadj", 
                              "D100_HumanLFC", "D100_HumanPadj", 
                              "D150_HumanLFC", "D150_HumanPadj", 
                              "D200_HumanLFC", "D200_HumanPadj")
dim(hyCS.DESeq2.ASE)

# Merged results (filtered for genes annotated in both genomes)
hyCS.DESeq2.ASE.merge = data.frame(cbind(d50_GRCh38_res[hc_genes[which(!(hc_genes %in% all_biased_genes))], c("log2FoldChange", "padj")], 
                                         d50_PanTro5_res[hc_genes[which(!(hc_genes %in% all_biased_genes))], c("log2FoldChange", "padj")],
                                         d100_GRCh38_res[hc_genes[which(!(hc_genes %in% all_biased_genes))], c("log2FoldChange", "padj")], 
                                         d100_PanTro5_res[hc_genes[which(!(hc_genes %in% all_biased_genes))], c("log2FoldChange", "padj")],
                                         d150_GRCh38_res[hc_genes[which(!(hc_genes %in% all_biased_genes))], c("log2FoldChange", "padj")], 
                                         d150_PanTro5_res[hc_genes[which(!(hc_genes %in% all_biased_genes))], c("log2FoldChange", "padj")],
                                         d200_GRCh38_res[hc_genes[which(!(hc_genes %in% all_biased_genes))], c("log2FoldChange", "padj")], 
                                         d200_PanTro5_res[hc_genes[which(!(hc_genes %in% all_biased_genes))], c("log2FoldChange", "padj")]))
colnames(hyCS.DESeq2.ASE.merge) = c("D50_HumanLFC", "D50_HumanPadj", "D50_ChimpLFC", "D50_ChimpPadj",
                                    "D100_HumanLFC", "D100_HumanPadj", "D100_ChimpLFC", "D100_ChimpPadj",
                                    "D150_HumanLFC", "D150_HumanPadj", "D150_ChimpLFC", "D150_ChimpPadj",
                                    "D200_HumanLFC", "D200_HumanPadj", "D200_ChimpLFC", "D200_ChimpPadj")
dim(hyCS.DESeq2.ASE.merge)

# Unbiased merged results (only genes in both annotations, only strictly unbiased genes)
hyCS.DESeq2.highConfidenceASE = data.frame(cbind(d50_GRCh38_res[all_unbiased_genes, c("log2FoldChange", "padj")], 
                                                 d50_PanTro5_res[all_unbiased_genes, c("log2FoldChange", "padj")],
                                                 d100_GRCh38_res[all_unbiased_genes, c("log2FoldChange", "padj")], 
                                                 d100_PanTro5_res[all_unbiased_genes, c("log2FoldChange", "padj")],
                                                 d150_GRCh38_res[all_unbiased_genes, c("log2FoldChange", "padj")], 
                                                 d150_PanTro5_res[all_unbiased_genes, c("log2FoldChange", "padj")],
                                                 d200_GRCh38_res[all_unbiased_genes, c("log2FoldChange", "padj")], 
                                                 d200_PanTro5_res[all_unbiased_genes, c("log2FoldChange", "padj")]))
colnames(hyCS.DESeq2.highConfidenceASE) = c("D50_HumanLFC", "D50_HumanPadj", "D50_ChimpLFC", "D50_ChimpPadj",
                                            "D100_HumanLFC", "D100_HumanPadj", "D100_ChimpLFC", "D100_ChimpPadj",
                                            "D150_HumanLFC", "D150_HumanPadj", "D150_ChimpLFC", "D150_ChimpPadj",
                                            "D200_HumanLFC", "D200_HumanPadj", "D200_ChimpLFC", "D200_ChimpPadj")
dim(hyCS.DESeq2.highConfidenceASE)

write.table(file = paste(analysis_dir, "hyCS.DESEq2.ASE.rev.txt", sep = "/"), x = hyCS.DESeq2.ASE, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
write.table(file = paste(analysis_dir, "hyCS.DESEq2.ASE.merge.rev.txt", sep = "/"), x = hyCS.DESeq2.ASE.merge, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
write.table(file = paste(analysis_dir, "hyCS.DESEq2.ASE.highConf.rev.txt", sep = "/"), x = hyCS.DESeq2.highConfidenceASE, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")


# Genes and chromosomes to remove from subsequent analyses
chromKeep = which(geneInfo$chrom!="chr20" & geneInfo$chrom!="chr18" & geneInfo$chrom!="chrX" & geneInfo$chrom!="chrY")
geneKeep = rownames(hyCS)[chromKeep][which(!(rownames(hyCS)[chromKeep] %in% all_biased_genes))]
geneKeepIdx = which(geneInfo$gene %in% geneKeep) 

write.table(file = paste(analysis_dir, "geneKeep.rev.txt", sep = "/"), x = cbind(geneKeep, geneKeepIdx), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
geneKeepInfo = read.table(file = paste(analysis_dir, "geneKeep.rev.txt", sep = "/"), header = TRUE, fill = TRUE)
geneKeep = as.vector(geneKeepInfo$geneKeep)
geneKeepIdx = as.vector(geneKeepInfo$geneKeepIdx)

ase_hc = read.table(paste(analysis_dir, "hyCS.DESEq2.ASE.highConf.rev.txt", sep = "/"), header = TRUE, row.names = 1, sep = "\t", fill = TRUE)
dim(ase_hc)

write.table(file = paste(analysis_dir, "hyCS.DESEq2.ASE.highConf.noAneuploidy.noSex.rev.txt", sep = "/"), 
            x = ase_hc[which(rownames(ase_hc) %in% geneKeep),], col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")


# Calculate TPM and CPM
hyCS_tot_keep = hyCS_tot[geneKeep,]
hyCS_tot_cpm = sweep(as.matrix(hyCS_tot_keep), 2, as.double(colSums(hyCS_tot_keep)/1000000), `/`)
hyCS_tot_fpkm = sweep(as.matrix(hyCS_tot_cpm), 1, as.double(geneInfo$genesizeKb[geneKeepIdx]), `/`)
hyCS_tot_rpk = sweep(as.matrix(hyCS_tot_keep), 1, as.double(geneInfo$genesizeKb[geneKeepIdx]), `/`)
hyCS_tot_tpm = sweep(as.matrix(hyCS_tot_rpk), 2, as.double(colSums(hyCS_tot_rpk)/1000000), `/`)
dim(hyCS_tot_tpm)
peek(hyCS_tot_tpm)

# Only the human and chimp counts
hyCS_hc = hyCS[, c(seq(2, length(hyCS), 3), seq(3, length(hyCS), 3))]
dim(hyCS_hc)
hyCS_hc_keep = hyCS_hc[geneKeep,]
hyCS_hc_cpm = sweep(as.matrix(hyCS_hc_keep), 2, as.double(colSums(hyCS_hc)/1000000), `/`)
hyCS_hc_rpk = sweep(as.matrix(hyCS_hc_keep), 1, as.double(geneInfo$genesizeKb[geneKeepIdx]), `/`)
hyCS_hc_tpm = sweep(as.matrix(hyCS_hc_rpk), 2, as.double(colSums(hyCS_hc_rpk)/1000000), `/`)
dim(hyCS_hc_tpm)
peek(hyCS_hc_tpm)

# Save allelic TPM (for plots in figure 4)
write.table(file = paste(analysis_dir, "hyCS_allelic_TPM.rev.txt", sep = "/"),
            hyCS_hc_tpm, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")

# Save ase per sample (for plots in figure 4)
hyCS_h = hyCS_h[geneKeep,]
hyCS_c = hyCS_c[geneKeep,]
ase = log2(hyCS_h/hyCS_c)
dim(ase)
write.table(file = paste(analysis_dir, "hyCS_ase_perSample.rev.txt", sep = "/"),
            ase, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")

# Replicate heat map
hylabels = paste(hyCS_meta$Line, hyCS_meta$Age, sep = "_")

# Extended Data Figure 4f
pdf(paste(plot_dir, "HybridTimecourseHeatmap.rev.pdf", sep = "/")) 
fontsize = 0.75
hyCS_cor = cor(hyCS_tot_cpm, use="complete.obs")
dim(hyCS_cor)
colnames(hyCS_cor) = hylabels
rownames(hyCS_cor) = hylabels
col<- colorRampPalette(c("turquoise4", "white", "darkorange2"))(200)
par(mai=c(1,0.5,0.5,0.25))
heatmap.2(x = hyCS_cor, col = col, srtCol = 55, symm = FALSE, dendrogram="none", trace="none", density.info="none", keysize = 0.9, breaks = seq(0.5, 1, length.out=201), cexCol = fontsize, cexRow = fontsize)
dev.off()


# PCA of hybrid CS with iPS cells
# With hybrids split
hyCS_c = hyCS[, seq(3,length(hyCS),by=3)]
dim(hyCS_c)
head(hyCS_c)

hyCS_h = hyCS[, seq(2,length(hyCS),by=3)]
dim(hyCS_h)
head(hyCS_h)

hyCS_hy = cbind(hyCS_c, hyCS_h)
dim(hyCS_hy)

hyCS_hy_cpm = sweep(as.matrix(hyCS_hy), 2, as.vector(colSums(hyCS_hy)/1000000), `/`)
dim(hyCS_hy_cpm)

hyCS_hy_noAnSex = hyCS_hy_cpm[geneKeep,]
dim(hyCS_hy_noAnSex)
head(hyCS_hy_noAnSex)

# iPS split species
iPS = read.table(paste(data_dir, "iPSC.GRCh38.txt", sep = "/"), header = TRUE)
rownames(iPS) = geneInfo$gene
dim(iPS)
iPS_meta = read.table(paste(data_dir, "iPSC.meta.txt", sep = "/"), header = TRUE)
iPS_meta$Species = as.character(iPS_meta$Species)
iPS_meta$Line = as.character(iPS_meta$Line)
iPS_tot = iPS[, seq(1,66,by=3)]
iPS_tot_keep = iPS_tot[geneKeep,]
iPS_tot_cpm = sweep(as.matrix(iPS_tot_keep), 2, as.double(colSums(iPS_tot_keep)/1000000), `/`)


iPS_hum_h = iPS[, seq(2, 65, by = 3)]
iPS_hum_h = iPS_hum_h[,seq(1,6)]
dim(iPS_hum_h)
head(iPS_hum_h)

iPS_chi_c = iPS[, seq(3, 66, by = 3)]
iPS_chi_c = iPS_chi_c[,seq(7,12)]
dim(iPS_chi_c)
head(iPS_chi_c)

iPS_c = iPS[, seq(3, 66, by = 3)]
iPS_c = iPS_c[, seq(13,22)]
dim(iPS_c)
head(iPS_c)

iPS_h = iPS[, seq(2, 65, by = 3)]
iPS_h = iPS_h[, seq(13,22)]
dim(iPS_h)
head(iPS_h)

iPS_hy_par = cbind(iPS_c, iPS_h, iPS_hum_h, iPS_chi_c)
head(iPS_hy_par)

iPS_hy_par = sweep(as.matrix(iPS_hy_par), 2, as.vector(colSums(iPS_hy_par)/1000000), `/`)
head(iPS_hy_par)

iPS_hy_par_noAnSex = iPS_hy_par[geneKeepIPS,]
dim(iPS_hy_par_noAnSex)
head(iPS_hy_par_noAnSex)

iPS_hy_par_noAnSex = iPS_hy_par[geneKeep,]

# Merge the two
hyCS_hy_par_noAnSex = cbind(hyCS_hy_noAnSex, iPS_hy_par_noAnSex)

# Expression threshold
hyCS_hy_par_comp = hyCS_hy_par_noAnSex[complete.cases(hyCS_hy_par_noAnSex),]

hyCS_hy_par_1 = hyCS_hy_par_comp[which(rowSums(hyCS_hy_par_comp)>1),]
dim(hyCS_hy_par_1)
head(hyCS_hy_par_1)

hyCS.hy_par.pca = prcomp(t(hyCS_hy_par_1), scale=TRUE, center=TRUE)
summary(hyCS.hy_par.pca)
hyCS.hy_par.pca_data = as.data.frame(hyCS.hy_par.pca$x)
hyCS_hy_par_species = c(rep("Hy-Chimp", 25), rep("Hy-Human", 25), rep("Hy-Chimp", 10), rep("Hy-Human", 10), rep("Human", 6), rep("Chimp", 6))


# Extended Data Figure 4g
pdf(paste(plot_dir, "hyCS_PCA_hyPar.rev.pdf", sep = "/"), height = 2, 2.5, useDingbats = FALSE) 
fontsize = 8
p = ggplot(hyCS.hy_par.pca_data, aes(PC1, PC2)) + 
  geom_point(aes(col=hyCS_hy_par_species, shape=c(rep("CS", 50), rep("iPS", 32))), size = 1)  + 
  ggtitle("PCA") +
  xlab("PC1 (32% of variance)") +
  ylab("PC2 (9% of variance)") +
  theme_bw() +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values = c("blue4","firebrick4", "deepskyblue", "indianred1")) +
  #labs(color = "Species") +
  theme(legend.position = "none")
p 
dev.off()

# PCA's of only organoid samples
dim(hyCS_tot_cpm) 

# Hybrids with parents
# Human
hCS_200 = read.table(paste(data_dir, "hCS_200.GRCh38.txt", sep = "/"), header = TRUE)
hCS_200_meta = read.table(paste(data_dir, "hCS_200.meta.txt", sep = "/"), header = TRUE)

hCS_35 = read.table(paste(data_dir, "hCS_35.GRCh38.txt", sep = "/"), header = TRUE)
hCS_35_meta = read.table(paste(data_dir, "hCS_35.meta.txt", sep = "/"), header = TRUE)

# Merge the two timecourses into one
hCS = as.data.frame(cbind(hCS_35, hCS_200))
rownames(hCS) = geneInfo$gene

hCS_meta = as.data.frame(rbind(hCS_35_meta, hCS_200_meta))
hCS_meta$Species = as.character(hCS_meta$Species)
hCS_meta$Line = as.character(hCS_meta$Line)

# Chimp
cCS_200 = read.table(paste(data_dir, "cCS_200.GRCh38.txt", sep = "/"), header = TRUE)
cCS_200_meta = read.table(paste(data_dir, "cCS_200.meta.txt", sep = "/"), header = TRUE)

cCS_35 = read.table(paste(data_dir, "cCS_35.GRCh38.txt", sep = "/"), header = TRUE)
cCS_35_meta = read.table(paste(data_dir, "cCS_35.meta.txt", sep = "/"), header = TRUE)

# Merge the two timecourses into one
cCS = as.data.frame(cbind(cCS_35, cCS_200))
rownames(cCS) = geneInfo$gene

cCS_meta = as.data.frame(rbind(cCS_35_meta, cCS_200_meta))
cCS_meta$Species = as.character(cCS_meta$Species)
cCS_meta$Line = as.character(cCS_meta$Line)


hCS_tot = hCS[, seq(1,123,by=3)]
hCS_tot = hCS_tot[geneKeep,]
hCS_tot_cpm = sweep(as.matrix(hCS_tot), 2, as.double(colSums(hCS_tot)/1000000), `/`)
cCS_tot = cCS[, seq(1,117,by=3)]
cCS_tot = cCS_tot[geneKeep,]
cCS_tot_cpm = sweep(as.matrix(cCS_tot), 2, as.double(colSums(cCS_tot)/1000000), `/`)
hcCS_tot_cpm = as.data.frame(cbind(hCS_tot_cpm, cCS_tot_cpm))

hyCS_par_cpm = cbind(hyCS_tot_cpm, hcCS_tot_cpm)
dim(hyCS_par_cpm)

time_hy_par = c(hyCS_meta$Age, hCS_meta$Age, cCS_meta$Age)


hyCS_par_cpm_comp = hyCS_par_cpm[complete.cases(hyCS_par_cpm),]
hyCS_par_cpm_1 = hyCS_par_cpm_comp[which(rowSums(hyCS_par_cpm_comp)>1),]
dim(hyCS_par_cpm_1)

hyCS_par.pca = prcomp(t(hyCS_par_cpm_1), scale=TRUE, center=TRUE)
summary(hyCS_par.pca)
hyCS_par.pca_data = as.data.frame(hyCS_par.pca$x)

hyCS_par_species = c(rep("Hybrid", 25), rep("Human", 41),rep("Chimp", 39))

# Extended Data Figure 7a,b
pdf(paste(plot_dir, "HybridParentTimecoursePCA.rev.pdf", sep = "/"), height = 1.65, width = 2.5, useDingbats = FALSE)
fontsize = 8
pointsize = 0.6

p1 = ggplot(hyCS_par.pca_data, aes(PC1, PC2)) + 
  geom_point(aes(col=hyCS_par_species, shape = hyCS_par_species), size = pointsize)  + 
  theme(title = element_blank()) +
  xlab("PC1 (14% of variance)") +
  ylab("PC2 (12% of variance)") +
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values = c("dodgerblue2","salmon", "mediumpurple2")) +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize)) + 
  labs(color = "Species") +
  theme(legend.position = "none")
p1
p1 = ggplot(hyCS_par.pca_data, aes(PC2, PC3)) + 
  geom_point(aes(col=hyCS_par_species, shape = hyCS_par_species), size = pointsize)  + 
  theme(title = element_blank()) +
  xlab("PC2 (12% of variance)") +
  ylab("PC3 (8% of variance)") +
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values = c("dodgerblue2","salmon", "mediumpurple2")) +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize)) + 
  labs(color = "Species") +
  theme(legend.position = "none")
p1

p2 = ggplot(hyCS_par.pca_data, aes(PC1, PC2)) + 
  geom_point(aes(col=as.factor(time_hy_par), shape = hyCS_par_species), size = pointsize)  + 
  theme(title = element_blank()) +
  xlab("PC1 (14% of variance)") +
  ylab("PC2 (12% of variance)") +
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values = c("brown2", "chocolate1", "darkgoldenrod1", "forestgreen", "dodgerblue2", "darkslateblue", "darkorchid3")) +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize)) + 
  labs(color = "Age", shape = "Species") +
  theme(legend.position = "none") #+
#theme(legend.position = "bottom", legend.box = "horizontal")
p2
dev.off()


# Heat maps of marker gene expression
# Heatmap of marker gene expression (forebrain and mesenchyme)
markers = c("FOXG1", "EMX1", "EOMES", "TBR1", "DCN", "LUM","COL1A1", "COL1A2")
hyCS.DESeq2.ASE = read.table(file = paste(analysis_dir, "hyCS.DESEq2.ASE.rev.txt", sep = "/"), header = TRUE, row.names = 1)

df = t(hyCS.DESeq2.ASE[markers, c(1,3,5,7)])
dim(df)

row.names(df) = c(50,100,150,200)
df
dfm = melt(cbind(df))
dfm
colnames(dfm) = c("Day", "Gene", "ASE")
dfm$Day <- factor(dfm$Day, levels = c(50,100,150,200))
dfm$Gene <- factor(dfm$Gene, levels = c("FOXG1", "EMX1", "EOMES",  "TBR1", "DCN", "LUM", "COL1A2", "COL1A1"))


# Figure 3d
pdf(paste(plot_dir, "HybridMarkersHeatmap.rev.pdf", sep = "/"), height = 1.35, width = 1.5) 
fontsize = 8
p= ggplot(dfm, aes(factor(Day), Gene)) + 
  geom_tile(aes(fill = ASE), color = "black", height = 1, width=1) +
  scale_fill_gradientn(colours = c("dodgerblue4", "white", "indianred3"), limits = c(-10,10)) +
  theme_bw()+
  theme(axis.title.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(axis.ticks.x=element_blank()) +
  xlab("Timepoint (Days)")+
  labs(fill = "Log2(Human/Chimp)") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text.x=element_text(size=fontsize, angle = 45, hjust = 1), axis.text.y = element_text(size = fontsize), axis.title=element_text(size=fontsize)) +
  theme(title = element_blank())
p

dev.off()



markerMapHy <- function(marker_genes, comp = "comp"){
  require(plyr)
  fontsize = 8.0
  mark = data.frame(t(hyCS_tot_tpm[marker_genes,]))
  
  mark$Timepoint = hyCS_meta$Age
  mark_summary = group_by(mark, Timepoint)
  mark_summary = summarise_all(mark_summary, mean)
  mark_summary
  
  df = as.data.frame(mark_summary)[,2:length(mark_summary)]
  row.names(df) = mark_summary$Timepoint
  
  dfm = melt(cbind(df, ind = rownames(df)))
  dfm <- ddply(dfm, .(variable), transform, rescale = rescale(log10(value)))
  dfm$ind <- factor(dfm$ind, levels = c(50,100,150,200))
  
  
  low = min(dfm$rescale) - 0.1
  high = max(dfm$rescale) + 0.1
  
  p1 =ggplot(dfm, aes(factor(ind), variable, width = 3, height = 1)) + 
    geom_tile(aes(fill = rescale), color = "black", height = 1, width=1) +
    scale_fill_gradient(low = "white", high = "darkgreen", limits = c(low, high)) +
    theme_bw()+
    ylab("") +
    theme(axis.ticks.x=element_blank()) +
    xlab("Timepoint (Days)")+
    theme(axis.title.y=element_blank(),
          axis.ticks.y=element_blank()) +
    labs(fill = "TPM") +
    theme(legend.position = "top", legend.box = "horizontal") +
    theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize)) + 
    theme(legend.text=element_text(size=fontsize))+
    theme(plot.title = element_text(hjust = 0.5, size=fontsize)) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank())
  grid.arrange(p1,p1,ncol=1)
}

# Extended Data Figure 7d
pdf(paste(plot_dir, "HybridTimecourseMarkerGenes_legend.rev.pdf", sep = "/"), height = 1.75, width = 3)

# Astrogenesis
markerMapHy(c("CLU", "HEPACAM", "GFAP", "AGT"))

# Cell division
markerMapHy(c("MKI67", "TOP2A", "NEK2", "ASPM"))

# Neurogenesis
markerMapHy(c("TUBB3", "MLLT11", "MAP2", "RBFOX3"))

dev.off()

detach(package:plyr)


## CIBERSORT
write.table(file = paste(analysis_dir, "hyCS_ciber_mixture_011620.txt", sep = "/"), x = hyCS_tot_fpkm, quote=FALSE, col.names = TRUE, row.names = TRUE, sep = '\t')

# CIBERSORT results
hyb_cib = read.table(paste(analysis_dir, "hyCS_ciber_result_011620.txt", sep = "/"), header=TRUE, row.names = 1)
hyb_cib[1,]
colnames(hyb_cib) = c("mes", "epi","neuron", "astro", "pro", "p.value", "pearson", "rmse")
hyb_cib$Timepoint = hyCS_meta$Age

hy_time_ord = paste("D", hyCS_meta$Age, sep = "")
hy_line_ord = hyCS_meta$Line
hylabels = paste(hy_time_ord, hy_line_ord, sep = "_")

# Group
hyb_cib_summary = group_by(hyb_cib, Timepoint)

# Calculate mean and standard error for each species/timepoint (normalized as a percentage of the neuralectoderm)
hyb_cib_summary = summarise(hyb_cib_summary, mean_pro = mean(pro/(pro+neuron+astro)), sd_pro = sd(pro/(pro+neuron+astro)), se_pro = sd(pro/(pro+neuron+astro))/sqrt(n()), mean_n = mean(neuron/(pro+neuron+astro)), sd_n = sd(neuron/(pro+neuron+astro)), se_n = sd(neuron/(pro+neuron+astro))/sqrt(n()), mean_a = mean(astro/(pro+neuron+astro)), sd_a = sd(astro/(pro+neuron+astro)), se_a = sd(astro/(pro+neuron+astro))/sqrt(n()))
hyb_cib_summary

# NOT USED
pdf(paste(plot_dir, "HybridTimecourseProportions.pdf", sep = "/"), height = 2, 1.6, useDingbats = FALSE)
fontsize = 8
pointsize = 0.3
p1 <- ggplot() +
  geom_bar(data = hyb_cib_summary, aes(x=factor(Timepoint), y=mean_pro), fill = "grey29", alpha = 0.8, stat = "identity", position = "dodge", width = 0.7, colour = "black") +
  geom_errorbar(data = hyb_cib_summary, aes(x=factor(Timepoint), ymin = mean_pro - se_pro, ymax = mean_pro + se_pro), width=0.2, position = position_dodge(0.7)) +
  geom_point(data=hyb_cib, aes(x=factor(Timepoint), y=pro/(pro+astro+neuron)), colour = "black", position = position_dodge(0.7), size = pointsize) +
  #geom_smooth(data = hyb_cib_summary[which(hyb_cib_summary$species == "human"),],aes(x=as.numeric(factor(Timepoint)), y=mean_pro), colour = "indianred3", size = 0.5, group = 1, linetype = "longdash", se=FALSE) +
  #geom_smooth(data = hyb_cib_summary[which(hyb_cib_summary$species == "chimp"),],aes(x=as.numeric(factor(Timepoint)), y=mean_pro), colour = "dodgerblue4", size = 0.5, group = 1, linetype = "longdash", se=FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 5)) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(size=fontsize))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab("Percentage") +
  ggtitle("Cycling progenitors") +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))+ theme(legend.title = element_blank()) 

p2 <- ggplot() +
  geom_bar(data = hyb_cib_summary, aes(x=factor(Timepoint), y=mean_n), fill = "steelblue4", alpha = 0.8, stat = "identity", position = "dodge", width = 0.7, colour = "black") +
  geom_errorbar(data = hyb_cib_summary, aes(x=factor(Timepoint), ymin = mean_n - se_n, ymax = mean_n + se_n), width=0.2, position = position_dodge(0.7)) +
  geom_point(data=hyb_cib, aes(x=factor(Timepoint), y=neuron/(pro+astro+neuron)), colour = "black", position = position_dodge(0.7), size = pointsize) +
  #geom_smooth(data = hyb_cib_summary[which(hyb_cib_summary$species == "human"),],aes(x=as.numeric(factor(Timepoint)), y=mean_n), colour = "indianred3", size = 0.5, group = 1, linetype = "longdash", se=FALSE) +
  #geom_smooth(data = hyb_cib_summary[which(hyb_cib_summary$species == "chimp"),],aes(x=as.numeric(factor(Timepoint)), y=mean_n), colour = "dodgerblue4", size = 0.5, group = 1, linetype = "longdash", se=FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 5)) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(size=fontsize))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab("Percentage") +
  ggtitle("Neurons") +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))+ theme(legend.title = element_blank()) 

p3 <- ggplot() +
  geom_bar(data = hyb_cib_summary, aes(x=factor(Timepoint), y=mean_a), fill = "red4", alpha = 0.8, stat = "identity", position = "dodge", width = 0.7, colour = "black") +
  geom_errorbar(data = hyb_cib_summary, aes(x=factor(Timepoint), ymin = mean_a - se_a, ymax = mean_a + se_a), width=0.2, position = position_dodge(0.7)) +
  geom_point(data=hyb_cib, aes(x=factor(Timepoint), y=astro/(pro+astro+neuron)), colour = "black", position = position_dodge(0.7), size = pointsize) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 5)) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(size=fontsize))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab("Percentage") +
  ggtitle("Astroglia") +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))+ theme(legend.title = element_blank()) 
grid.arrange(p1,p2,p3,ncol=1)
dev.off()


# Raw proportions
# Extended Data Figure 7c
pdf(paste(plot_dir, "HybridTimecourseProportions.perSample.pdf", sep = "/"), height = 4.5, width = 6)
fontsize = 8

t = as.data.frame(t(as.matrix(hyb_cib[, 1:5])))
t
rownames(t) = c(4,5,2,1,3)

tm <- melt(cbind(t, ind = rownames(t)), id.vars = c('ind'))
tm$Timepoint = rep(hyCS_meta$Age, each=5)
tm$Timepoint = factor(tm$Timepoint, levels = c(50,100, 150, 200), ordered = TRUE)
tm = tm[order(tm$Timepoint),]
tm$variable <- factor(tm$variable, levels = unique(tm$variable))
ggplot(tm,aes(x = variable, y = value, fill = ind)) + 
  geom_bar(position = "fill",stat = "identity", alpha = 0.85) +
  scale_y_continuous(labels = percent_format()) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c("red4", "steelblue4", "gray29", "gold", "burlywood"),
                    name = "Cell Type",
                    labels=c("Astroglia", "Neuron", "Cycling progenitors", "Mesenchyme", "Epithelial")) +
  theme(plot.title = element_text(size=fontsize))+
  xlab("Sample") + 
  ylab("Percentage") +
  ggtitle("Hybrid samples: all cell types") +
  theme(legend.position = "bottom", legend.box = "horizontal")+
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize),axis.text.x = element_text(angle = 65, hjust = 1)) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_x_discrete(labels = hylabels)
dev.off()

# Human
h_cib = read.table(paste(analysis_dir, "hCS_ciber_results_011620.txt", sep = "/"), header=TRUE, row.names = 1)

# Chimp
c_cib = read.table(paste(analysis_dir, "cCS_ciber_results_011620.txt", sep = "/"), header=TRUE, row.names = 1)

# Merge species
hc_cib = rbind(h_cib, c_cib)
colnames(hc_cib) = c("mes", "epi","neuron", "astro", "pro", "p.value", "pearson", "rmse")
hc_cib

# Add metadata
hc_cib$Species = c(hCS_meta$Species, cCS_meta$Species)
hc_cib$Timepoint = c(hCS_meta$Age, cCS_meta$Age)
hc_cib_summary <- hc_cib
hc_cib_summary = group_by(hc_cib_summary, Timepoint, Species)

# Calculate mean and standard error for each species/timepoint (normalized as a percentage of the neuralectoderm)
# detach(package:plyr) # unload it if it was loaded
hc_cib_summary = summarise(hc_cib_summary, 
                           mean_pro = mean(pro/(pro+neuron+astro)), 
                           sd_pro = sd(pro/(pro+neuron+astro)), 
                           se_pro = sd(pro/(pro+neuron+astro))/sqrt(n()), 
                           mean_n = mean(neuron/(pro+neuron+astro)), 
                           sd_n = sd(neuron/(pro+neuron+astro)), 
                           se_n = sd(neuron/(pro+neuron+astro))/sqrt(n()), 
                           mean_a = mean(astro/(pro+neuron+astro)), 
                           sd_a = sd(astro/(pro+neuron+astro)), 
                           se_a = sd(astro/(pro+neuron+astro))/sqrt(n()))
hc_cib_summary

# Figure 3e
pdf(paste(plot_dir, "ParentalTimecourseProportions.pdf", sep = "/"), height = 3.5, width = 2.5, useDingbats = FALSE)
pointsize = 0.4
p1 <- ggplot()+
  geom_bar(data = hc_cib_summary, aes(x=factor(Timepoint), y=mean_pro, fill=factor(Species)), alpha = 0.8, stat = "identity", position = "dodge", width = 0.5, colour = "black") +
  geom_errorbar(data = hc_cib_summary, aes(x=factor(Timepoint), fill=factor(Species), ymin = mean_pro - se_pro, ymax = mean_pro + se_pro), width=0.2, position = position_dodge(0.5)) +
  geom_smooth(data = hc_cib_summary[which(hc_cib_summary$Species == "Human"),],aes(x=as.numeric(factor(Timepoint)), y=mean_pro), colour = "indianred3", size = 0.5, group = 1, linetype = "longdash", se=FALSE) +
  stat_smooth(data = hc_cib_summary[which(hc_cib_summary$Species == "Human"),],aes(x=as.numeric(factor(Timepoint)), y=mean_pro),
              alpha = 0.3,fill = 'salmon', se = FALSE, geom = "area", method = "loess")+
  geom_smooth(data = hc_cib_summary[which(hc_cib_summary$Species == "Chimp"),],aes(x=as.numeric(factor(Timepoint)), y=mean_pro), colour = "dodgerblue4", size = 0.5, group = 1, linetype = "longdash", se=FALSE) +
  stat_smooth(data = hc_cib_summary[which(hc_cib_summary$Species == "Chimp"),],aes(x=as.numeric(factor(Timepoint)), y=mean_pro),
              alpha = 0.3,fill = 'dodgerblue2', se = FALSE, geom = "area", method = "loess")+
  geom_point(data=hc_cib, aes(x=factor(Timepoint), y=pro/(pro+astro+neuron), fill=factor(Species)), colour = "grey15", position = position_dodge(0.5), size = pointsize) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 5), limits = c(0, 0.99)) + 
  theme_bw() + 
  xlab("") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0, size=8))+
  # theme(axis.title.x=element_blank(),
  #        axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) +
  ylab("Percentage") +
  scale_fill_manual(values=c("dodgerblue3", "salmon"))+
  ggtitle("Cycling progenitors") +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8))+ theme(legend.title = element_blank()) +
  theme(legend.position = "none")

p2 <- ggplot() +
  geom_bar(data = hc_cib_summary, aes(x=factor(Timepoint), y=mean_n, fill=factor(Species)), alpha = 0.8,stat = "identity", position = "dodge", width = 0.5, colour = "black") +
  geom_errorbar(data = hc_cib_summary, aes(x=factor(Timepoint), fill=factor(Species), ymin = mean_n - se_n, ymax = mean_n + se_n), width=0.2, position = position_dodge(0.5)) +
  geom_smooth(data = hc_cib_summary[which(hc_cib_summary$Species == "Human"),],aes(x=as.numeric(factor(Timepoint)), y=mean_n), colour = "indianred3", size = 0.5, group = 1,linetype = "longdash", se=FALSE) +
  stat_smooth(data = hc_cib_summary[which(hc_cib_summary$Species == "Human"),],aes(x=as.numeric(factor(Timepoint)), y=mean_n),
              alpha = 0.3,fill = 'salmon', se = FALSE, geom = "area", method = "loess")+
  geom_smooth(data = hc_cib_summary[which(hc_cib_summary$Species == "Chimp"),],aes(x=as.numeric(factor(Timepoint)), y=mean_n), colour = "dodgerblue4", size = 0.5, group = 1, linetype = "longdash", se=FALSE) +
  stat_smooth(data = hc_cib_summary[which(hc_cib_summary$Species == "Chimp"),],aes(x=as.numeric(factor(Timepoint)), y=mean_n),
              alpha = 0.3,fill = 'dodgerblue2', se = FALSE, geom = "area", method = "loess")+
  geom_point(data=hc_cib, aes(x=factor(Timepoint), y=neuron/(pro+astro+neuron), fill=factor(Species)), colour = "grey15", position = position_dodge(0.5), size = pointsize) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 5), limits = c(0, 0.99)) + 
  theme_bw() +
  xlab("") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0, size=8))+
  # theme(axis.title.x=element_blank(),
  #        axis.text.x=element_blank(),
  #        axis.ticks.x=element_blank()) +
  ylab("Percentage") +
  scale_fill_manual(values=c("dodgerblue3", "salmon"))+
  ggtitle("Neurons") +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8))+ theme(legend.title = element_blank()) +
  theme(legend.position = "none")

p3 <- ggplot() +
  geom_bar(data = hc_cib_summary, aes(x=factor(Timepoint), y=mean_a, fill=factor(Species)), alpha = 0.8, stat = "identity", position = "dodge", width = 0.5, colour = "black") +
  geom_errorbar(data = hc_cib_summary, aes(x=factor(Timepoint), fill=factor(Species), ymin = mean_a - se_a, ymax = mean_a + se_a), width=0.2, position = position_dodge(0.5)) +
  geom_smooth(data = hc_cib_summary[which(hc_cib_summary$Species == "Human"),],aes(x=as.numeric(factor(Timepoint)), y=mean_a), colour = "indianred3", size = 0.5, group = 1, se=FALSE, linetype = "longdash") +
  stat_smooth(data = hc_cib_summary[which(hc_cib_summary$Species == "Human"),],aes(x=as.numeric(factor(Timepoint)), y=mean_a),
              alpha = 0.3,fill = 'salmon', se = FALSE, geom = "area", method = "loess")+
  geom_smooth(data = hc_cib_summary[which(hc_cib_summary$Species == "Chimp"),],aes(x=as.numeric(factor(Timepoint)), y=mean_a), colour = "dodgerblue4", size = 0.5, group = 1, se=FALSE, linetype = "longdash") +
  stat_smooth(data = hc_cib_summary[which(hc_cib_summary$Species == "Chimp"),],aes(x=as.numeric(factor(Timepoint)), y=mean_a),
              alpha = 0.3,fill = 'dodgerblue2', se = FALSE, geom = "area", method = "loess")+
  geom_point(data=hc_cib, aes(x=factor(Timepoint), y=astro/(pro+astro+neuron), fill=factor(Species)), colour = "grey15", position = position_dodge(0.5), size = pointsize) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 5), limits = c(0, 0.99)) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(plot.title = element_text(hjust = 0, size=8))+
  xlab("Timepoint (Days)") + 
  ylab("Percentage") +
  scale_fill_manual(values=c("dodgerblue3", "salmon"))+
  ggtitle("Astroglia") +
  theme(axis.text=element_text(size=8), axis.title=element_text(size=8)) + theme(legend.title = element_blank()) +
  theme(legend.position = "none")

grid.arrange(p1,p2,p3,ncol=1)

dev.off()

# Plot all cell type proportions for every sample individually (to show the variability)
# Extended Data Figure 7e,g
pdf(paste(plot_dir, "ParentalTimecourseProportions.perSample.pdf", sep = "/"), height = 4.5, width = 6)
fontsize = 8

# Human
hlabels = paste(hCS_meta$Line, "_D", hCS_meta$Age, sep = "")
t = as.data.frame(t(as.matrix(h_cib[, 1:5])))
rownames(t) = c(4,5,2,1,3)

tm <- melt(cbind(t, ind = rownames(t)), id.vars = c('ind'))
tm$Timepoint = rep(hCS_meta$Age, each=5)
tm$Timepoint = factor(tm$Timepoint, levels = c(15,25,35,50,100, 150, 200), ordered = TRUE)
tm = tm[order(tm$Timepoint),]
h_time_raw = factor(hCS_meta$Age, levels = c(15,25,35,50,100, 150, 200), ordered = TRUE)
hlabels_ord = hlabels[order(h_time_raw)]
tm$variable <- factor(tm$variable, levels = unique(tm$variable))
ph <- ggplot(tm,aes(x = variable, y = value, fill = ind)) + 
  geom_bar(position = "fill",stat = "identity", alpha = 0.85) +
  scale_y_continuous(labels = percent_format()) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c("red4", "steelblue4", "gray29", "gold", "burlywood"),
                    name = "Cell Type",
                    labels=c("Astrocyte", "Neuron", "Radial Glia", "Mesenchyme", "Epithelial")) +
  theme(plot.title = element_text(hjust = 0.5, size=fontsize))+
  xlab("Sample") + 
  ylab("Percentage") +
  theme(legend.position = "right", legend.box = "vertical")+
  theme(axis.text=element_text(size=fontsize-2), axis.title=element_text(size=fontsize),axis.text.x = element_text(angle = 65, hjust = 1)) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_x_discrete(labels = hlabels_ord)

# Chimp
clabels = paste(cCS_meta$Line, "_D", cCS_meta$Age, sep = "")
t = as.data.frame(t(as.matrix(c_cib[, 1:5])))
rownames(t) = c(4,5,2,1,3)

tm <- melt(cbind(t, ind = rownames(t)), id.vars = c('ind'))
tm$Timepoint = rep(cCS_meta$Age, each=5)
tm$Timepoint = factor(tm$Timepoint, levels = c(15,25,35,50,100, 150, 200), ordered = TRUE)
tm = tm[order(tm$Timepoint),]
c_time_raw = factor(cCS_meta$Age, levels = c(15,25,35,50,100, 150, 200), ordered = TRUE)
clabels_ord = clabels[order(c_time_raw)]
tm$variable <- factor(tm$variable, levels = unique(tm$variable))
pc <- ggplot(tm,aes(x = variable, y = value, fill = ind)) + 
  geom_bar(position = "fill",stat = "identity", alpha = 0.85) +
  scale_y_continuous(labels = percent_format()) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_fill_manual(values = c("red4", "steelblue4", "gray29", "gold", "burlywood"),
                    name = "Cell Type",
                    labels=c("Astrocyte", "Neuron", "Radial Glia", "Mesenchyme", "Epithelial")) +
  theme(plot.title = element_text(hjust = 0.5, size=fontsize))+
  xlab("Sample") + 
  ylab("Percentage") +
  theme(legend.position = "right", legend.box = "vertical")+
  theme(axis.text=element_text(size=fontsize-2), axis.title=element_text(size=fontsize),axis.text.x = element_text(angle = 65, hjust = 1)) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_x_discrete(labels = clabels_ord)

grid.arrange(ph,pc, nrow=2)
dev.off()




