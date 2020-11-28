## Analysis and plots for Fig. 1 and Extended Data Fig. 1, 2, 3: hybrid iPS cells 
## Everything here was performed in Rstudio on a Macbook Pro with 16G RAM, unless otherwise noted
## Clear cache in Rstudio
rm(list = ls(all.names = TRUE))

## Load required libraries
library(ggplot2)
library(gridExtra)
library(gplots)
library(DESeq2)
library(sva)

# Useful
peek <- function(data){
  rowmax = min(5, nrow(data))
  colmax = min(5, ncol(data))
  print(data[1:rowmax, 1:colmax])
}

## Directories where the data is / where to deposit files and plots
data_dir = "/Users/rachelagoglia/Desktop/Manuscript_backedUp/Revisions/ProcessedData"
analysis_dir = "/Users/rachelagoglia/Desktop/Manuscript_backedUp/Revisions/AnalysisFiles"
plot_dir = "/Users/rachelagoglia/Desktop/Manuscript_backedUp/Revisions/Plots"

## Load data 
# Gene meta data
geneInfo = read.table(paste(data_dir, "geneInfo.GRCh38.txt", sep = "/"))
names(geneInfo) = c("gene", "length", "chrom")
geneInfo$genesizeKb = as.vector(geneInfo$length/1000)
dim(geneInfo)
head(geneInfo)

# RNAseq data
iPS = read.table(paste(data_dir, "iPSC.GRCh38.txt", sep = "/"), header = TRUE)
rownames(iPS) = geneInfo$gene
dim(iPS)
iPS_meta = read.table(paste(data_dir, "iPSC.meta.txt", sep = "/"), header = TRUE)
iPS_meta$Species = as.character(iPS_meta$Species)
iPS_meta$Line = as.character(iPS_meta$Line)

# DESeq2 for ASE
# Human genome
counts_matrix = iPS[, c(seq(2, 65, by = 3), seq(3, 66, by = 3))]
counts_matrix = counts_matrix[,c(seq(13,22), seq(35, 44))]
num_samples = ncol(counts_matrix) / 2
cond_species <- factor(c(rep("Human", num_samples), rep("Chimp", num_samples))) 
cond_sample <- factor(rep(1:num_samples, 2)) 
cond_line <- factor(rep(rep(c("Hy1.30", "Hy1.25",  "Hy1.29", "Hy2.16", "Hy2.9"), each = 2), 2))

cond <- ~ cond_line + cond_species
cond_reduced <- ~ cond_line 
coldata <- data.frame(cond_line, cond_sample, cond_species)
dds <- DESeqDataSetFromMatrix(counts_matrix, coldata, cond)

dds <- DESeq(dds,test="LRT",betaPrior=FALSE,full=cond,reduced=cond_reduced)

iPS_h_res <- results(dds)
resultsNames(dds)

output_name = paste(analysis_dir, "iPS_ASE.DESeq2.GRCh38.v6.txt", sep = "/")

write.table(data.frame(iPS_h_res), 
            file=output_name, 
            sep="\t", quote = FALSE,na = "", row.names = TRUE, col.names = TRUE)


# Chimp genome gene meta data
geneInfo.pt = read.table(paste(data_dir, "geneInfo.PanTro5.txt", sep = "/"))
names(geneInfo.pt) = c("gene", "length", "chrom")
geneInfo.pt$genesizeKb = as.vector(geneInfo.pt$length/1000)
dim(geneInfo.pt)
head(geneInfo.pt)

# Chimp genome RNAseq data
iPS_pt = read.table(paste(data_dir, "iPSC.PanTro5.txt", sep = "/"), header = TRUE)
rownames(iPS_pt) = geneInfo.pt$gene
dim(iPS_pt)

counts_matrix = iPS_pt[, c(seq(2, 65, by = 3), seq(3, 66, by = 3))]
counts_matrix = counts_matrix[,c(seq(13,22), seq(35, 44))]
num_samples = ncol(counts_matrix) / 2
cond_species <- factor(c(rep("Human", num_samples), rep("Chimp", num_samples))) 
cond_sample <- factor(rep(1:num_samples, 2)) 
cond_line <- factor(rep(rep(c("Hy1.30", "Hy1.25",  "Hy1.29", "Hy2.16", "Hy2.9"), each = 2), 2))

cond <- ~ cond_line + cond_species
cond_reduced <- ~ cond_line 
coldata <- data.frame(cond_line, cond_sample, cond_species)
dds <- DESeqDataSetFromMatrix(counts_matrix, coldata, cond)

dds <- DESeq(dds,test="LRT",betaPrior=FALSE,full=cond,reduced=cond_reduced)

iPS_c_res <- results(dds)

output_name = paste(analysis_dir, "iPS_ASE.DESeq2.PanTro5.v6.txt", sep = "/")

write.table(data.frame(iPS_c_res), 
            file=output_name, 
            sep="\t", quote = FALSE,na = "", row.names = TRUE, col.names = TRUE)


# Mapping bias for ASE
hc_genes = rownames(iPS)[which(rownames(iPS) %in% rownames(iPS_pt))]
h_points = data.frame(iPS_h_res)[hc_genes, "log2FoldChange"]
c_points = data.frame(iPS_c_res)[hc_genes, "log2FoldChange"]

# Extended Data Figure 3d
pdf(paste(plot_dir, "iPS_mappingBiasASE.v6.pdf", sep = "/"), height = 2.3, 2.3, useDingbats = FALSE) 
pointsize = 0.3
fontsize = 8
p1 = ggplot() +
  geom_point(aes(x = h_points, y = c_points), size = pointsize) +      
  theme(title = element_blank()) + 
  xlab("GRCh38 ASE") +
  ylab("PanTro5 ASE") + 
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(slope=1, intercept = 0) +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))+
  geom_point(aes(h_points[which(abs(h_points-c_points)>1)], c_points[which(abs(h_points-c_points)>1)]), col = "red", size = pointsize) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p1
dev.off()


filt = intersect(which(is.finite(h_points)), which(is.finite(c_points)))
cor(h_points[filt], c_points[filt])*cor(h_points[filt], c_points[filt]) # 0.786


# Remove genes with mapping bias and aneuploid chromosomes
iPS_biased_genes = hc_genes[which(abs(h_points-c_points)>1)]
iPS_unbiased_genes = hc_genes[which(abs(h_points-c_points)<=1)]

length(iPS_biased_genes)
write.table(data.frame(cbind(iPS_biased_genes, h_points[which(hc_genes %in% iPS_biased_genes)], c_points[which(hc_genes %in% iPS_biased_genes)])), 
            file=paste(analysis_dir, "iPS.biased.rev.txt", sep = "/"), 
            sep="\t", quote = FALSE,na = "", row.names = FALSE, col.names = FALSE)



# Human only results (unfiltered)
iPS.DESeq2.ASE = data.frame(cbind(iPS_h_res[which(!(rownames(iPS_h_res) %in% iPS_biased_genes)), c("log2FoldChange", "padj")]))
colnames(iPS.DESeq2.ASE) = c("iPS_HumanLFC", "iPS_HumanPadj")   
dim(iPS.DESeq2.ASE)
write.table(data.frame(iPS.DESeq2.ASE), 
            file=paste(analysis_dir, "iPS.DESeq2.rev.txt", sep = "/"), 
            sep="\t", quote = FALSE,na = "", row.names = TRUE, col.names = TRUE)

# Merged results (filtered for genes annotated in both genomes)
iPS.DESeq2.ASE.merge = data.frame(cbind(iPS_h_res[hc_genes[which(!(hc_genes %in% iPS_biased_genes))], c("log2FoldChange", "padj")], 
                                        iPS_c_res[hc_genes[which(!(hc_genes %in% iPS_biased_genes))], c("log2FoldChange", "padj")]))
colnames(iPS.DESeq2.ASE.merge) = c("iPS_HumanLFC", "iPS_HumanPadj", "iPS_ChimpLFC", "iPS_ChimpPadj")                                         
dim(iPS.DESeq2.ASE.merge)
write.table(data.frame(iPS.DESeq2.ASE.merge), 
            file=paste(analysis_dir, "iPS.DESeq2.merge.rev.txt", sep = "/"), 
            sep="\t", quote = FALSE,na = "", row.names = TRUE, col.names = TRUE)

# Unbiased merged results (only genes in both annotations, only strictly unbiased genes)
highConf_iPS = data.frame(cbind(iPS_h_res[iPS_unbiased_genes, c("log2FoldChange", "padj")], iPS_c_res[iPS_unbiased_genes, c("log2FoldChange", "padj")]))
colnames(highConf_iPS) = c("HumanLFC", "HumanPadj", "ChimpLFC", "ChimpPadj")
dim(highConf_iPS)
write.table(data.frame(highConf_iPS), 
            file=paste(analysis_dir, "iPS.DESeq2.merge.highConfASE.rev.txt", sep = "/"), 
            sep="\t", quote = FALSE,na = "", row.names = TRUE, col.names = TRUE)


# Aneuploid chromosomes
chromKeepIPS = which(geneInfo$chrom!="chr20" & geneInfo$chrom!="chrX" & geneInfo$chrom!="chrY")
geneKeepIPS = rownames(iPS)[chromKeepIPS][which(!(rownames(iPS)[chromKeepIPS] %in% iPS_biased_genes))]
geneKeepIdxIPS = which(geneInfo$gene %in% geneKeepIPS) 

write.table(file = paste(analysis_dir, "iPSgeneKeep.rev.txt", sep = "/"), x = cbind(geneKeepIPS, geneKeepIdxIPS), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(data.frame(highConf_iPS[which(rownames(highConf_iPS) %in% geneKeepIPS),]), 
            file=paste(analysis_dir, "iPS.DESeq2.merge.highConfASE.noAneuploid.noSex.rev.txt", sep = "/"), 
            sep="\t", quote = FALSE,na = "", row.names = TRUE, col.names = TRUE)


keepIPS = read.table(paste(analysis_dir, "iPSgeneKeep.rev.txt", sep = "/"), header = TRUE)
geneKeepIPS = as.character(keepIPS$geneKeepIPS)
geneKeepIdxIPS = keepIPS$geneKeepIdxIPS

# Calculate TPM and CPM
iPS_tot = iPS[, seq(1,66,by=3)]
iPS_tot_keep = iPS_tot[geneKeepIPS,]
iPS_tot_cpm = sweep(as.matrix(iPS_tot_keep), 2, as.double(colSums(iPS_tot_keep)/1000000), `/`)
iPS_tot_fpkm = sweep(as.matrix(iPS_tot_cpm), 1, as.double(geneInfo$genesizeKb[geneKeepIdxIPS]), `/`)
iPS_tot_rpk = sweep(as.matrix(iPS_tot_keep), 1, as.double(geneInfo$genesizeKb[geneKeepIdxIPS]), `/`)
iPS_tot_tpm = sweep(as.matrix(iPS_tot_rpk), 2, as.double(colSums(iPS_tot_rpk)/1000000), `/`)


geneKeepIPS_pt = geneInfo.pt$gene[which(geneInfo.pt$gene %in% geneKeepIPS)]
geneKeepIPSIdx_pt = which(geneInfo.pt$gene %in% geneKeepIPS_pt)

iPS_pt_tot = iPS_pt[, seq(1,66,by=3)]
iPS_pt_tot_keep = iPS_pt_tot[geneKeepIPS_pt,]
iPS_pt_tot_cpm = sweep(as.matrix(iPS_pt_tot_keep), 2, as.double(colSums(iPS_pt_tot_keep)/1000000), `/`)
iPS_pt_tot_fpkm = sweep(as.matrix(iPS_pt_tot_cpm), 1, as.double(geneInfo.pt$genesizeKb[geneKeepIPSIdx_pt]), `/`)
iPS_pt_tot_rpk = sweep(as.matrix(iPS_pt_tot_keep), 1, as.double(geneInfo.pt$genesizeKb[geneKeepIPSIdx_pt]), `/`)
iPS_pt_tot_tpm = sweep(as.matrix(iPS_pt_tot_rpk), 2, as.double(colSums(iPS_pt_tot_rpk)/1000000), `/`)

# Plot
hc_genes = which(rownames(iPS_tot_tpm) %in% rownames(iPS_pt_tot_tpm))
ch_genes = which(rownames(iPS_pt_tot_tpm) %in% rownames(iPS_tot_tpm))
h_filt = iPS_tot_tpm[hc_genes, 13]
c_filt = iPS_pt_tot_tpm[ch_genes,13]

# Extended Data Figure 3c
pdf(paste(plot_dir, "iPS_mappingBias.rev.pdf", sep = "/"), height = 2.3, 2.3, useDingbats = FALSE) 
fontsize = 8
pointsize = 0.3
p1 = ggplot()+
  geom_point(aes(x = log10(h_filt), y= log10(c_filt)), size = pointsize) +      
  theme(title = element_blank()) + 
  xlab("GRCh38 Log10(TPM)") +
  ylab("PanTro5 Log10(TPM)") + 
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(slope=1, intercept = 0) +
  theme_bw() + 
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
p1
dev.off()

filt = intersect(which(is.finite(log10(c_filt))), which(is.finite(log10(h_filt))))
cor(log10(c_filt)[filt], log10(h_filt)[filt])*cor(log10(c_filt)[filt], log10(h_filt)[filt]) # 0.854

# Bar plots for pluripotency markers
fontsize = 8
barwidth = 0.6
pointsize = 0.5

iPS_bar <-function(gene, data = iPS_tot_tpm){
  require(dplyr)
  temp = data.frame(t(as.matrix(rbind(data[which(row.names(data) == gene),], iPS_meta$Species, iPS_meta$Line))))
  names(temp) = c("exp", "Species", "Sample")
  #temp$Sample = rownames(temp)
  temp$Species = factor(temp$Species, levels = c("Human", "Chimp", "Hybrid"))
  temp$Sample = factor(temp$Sample, levels=c("H1", "H2", "H3", "C1", "C2", "C3", "Hy1-25", "Hy1-29", "Hy1-30", "Hy2-9", "Hy2-16"), ordered = TRUE)
  temp = temp[order(temp$Species, temp$Sample),]
  temp2 = group_by(temp, Sample, Species)
  #temp2 = group_by(temp, Species)
  temp2 = summarise(temp2, 
                    mean_exp = mean(as.numeric(as.character(exp))), 
                    se_exp = sd(as.numeric(as.character(exp)))/sqrt(n()))
  p <- ggplot() + 
    #geom_bar(data = temp2, aes(x=Species, y=mean_exp, fill=Species), colour="black", stat = "identity",width=barwidth, position = position_dodge(width=barwidth), alpha = 0.8) +
    #geom_errorbar(data = temp2, aes(x=Species, ymin = mean_exp - se_exp, ymax = mean_exp + se_exp), width=0.2, position = position_dodge(0.2)) +
    #geom_point(data=temp, aes(x=Species, y=as.numeric(as.character(exp)), fill=Species), colour = "grey15", position = position_dodge(0.7), size = pointsize) +
    
    geom_bar(data = temp2, aes(x=Sample, y=mean_exp, fill=Species), colour="black", stat = "identity",width=barwidth, position = position_dodge(width=barwidth), alpha = 0.8) +
    geom_point(data=temp, aes(x=Sample, y=as.numeric(as.character(exp)), fill=Species), colour = "black", position = position_dodge(0.7), size = pointsize) +
    #geom_errorbar(data = temp2, aes(x=Sample, ymin = mean_exp - se_exp, ymax = mean_exp + se_exp), width=0.2, position = position_dodge(0.2)) +
    ggtitle(gene) + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values=c("salmon", "dodgerblue2", "mediumpurple2"))+
    theme(plot.title = element_text(size=fontsize))+
    xlab("Species") + 
    ylab("TPM") +
    theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize), axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank()) +
    scale_x_discrete(labels = temp2$Sample) +
    theme(legend.position = "none")
  p  
}

# Figure 1 e,f,g,h
pdf(paste(plot_dir, "iPS_Pluripotency_long.rev.pdf", sep = "/"), height = 1.5, width = 5.5, useDingbats = FALSE)
p1 = iPS_bar("NANOG")
p2 = iPS_bar("SOX2")
p3 = iPS_bar("MYCL")
p4 = iPS_bar("KLF4")
grid.arrange(p2,p1,p3,p4, ncol = 4)
dev.off()


# Hybrid allelic expression
iPS_h = iPS[, seq(2, 66, by = 3)]
iPS_h_keep = iPS_h[geneKeepIPS,]
iPS_h_rpk = sweep(as.matrix(iPS_h_keep), 1, as.double(geneInfo$genesizeKb[geneKeepIdxIPS]), `/`)
iPS_h_tpm = sweep(as.matrix(iPS_h_rpk), 2, as.double(colSums(iPS_h_rpk)/1000000), `/`)
iPS_h_tpm = iPS_h_tpm[,seq(13,22)]
dim(iPS_h_tpm)
head(iPS_h_tpm)

iPS_c = iPS[, seq(3, 66, by = 3)]
iPS_c_keep = iPS_c[geneKeepIPS,]
iPS_c_rpk = sweep(as.matrix(iPS_c_keep), 1, as.double(geneInfo$genesizeKb[geneKeepIdxIPS]), `/`)
iPS_c_tpm = sweep(as.matrix(iPS_c_rpk), 2, as.double(colSums(iPS_c_rpk)/1000000), `/`)
iPS_c_tpm = iPS_c_tpm[,seq(13,22)]
dim(iPS_c_tpm)
head(iPS_c_tpm)

iPS_hc_species = c(rep("Human", 10), rep("Chimp", 10))

iPS_hc_tpm = cbind(iPS_h_tpm, iPS_c_tpm)
dim(iPS_hc_tpm)

iPS_hc_sample = rep(c( "Hy1-30", "Hy1-25", "Hy1-29", "Hy2-16", "Hy2-9"), each = 2)

fontsize = 8
barwidth = 0.5
pointsize = 0.5

iPS_hy_bar <-function(gene, data = iPS_hc_tpm){
  temp = data.frame(t(as.matrix(rbind(data[which(row.names(data) == gene),], iPS_hc_species, iPS_hc_sample))))
  names(temp) = c("exp", "Species", "Sample")
  temp$Species = factor(temp$Species, levels = c("Human", "Chimp"))
  temp$Sample = factor(temp$Sample, levels=c("Hy1-25", "Hy1-29", "Hy1-30", "Hy2-9", "Hy2-16"), ordered = TRUE)
  temp = temp[order(temp$Species, temp$Sample),]
  temp2 = group_by(temp, Sample, Species)
  temp2 = summarise(temp2, mean_exp = mean(as.numeric(as.character(exp))), se_exp = sd(as.numeric(as.character(exp)))/sqrt(n()))
  p <- ggplot() + 
    geom_bar(data = temp2, aes(x=Sample, y=mean_exp, fill=Species), colour="black", stat = "identity",width=barwidth, position = position_dodge(width=barwidth), alpha = 0.8) +
    #geom_errorbar(data = temp2, aes(x=Sample, ymin = mean_exp - se_exp, ymax = mean_exp + se_exp, fill = Species), width=0.2, position = position_dodge(barwidth)) +
    geom_point(data=temp, aes(x=Sample, y=as.numeric(as.character(exp)), fill=Species), colour = "black", position = position_dodge(barwidth), size = pointsize) +
    ggtitle(paste("Allelic expression:", gene)) + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values=c("indianred1", "dodgerblue2"))+
    theme(plot.title = element_text(size=fontsize))+
    xlab("Species") + 
    ylab("TPM") +
    theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize), axis.text.x = element_text(angle = 65, hjust = 1)) +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank()) +
    scale_x_discrete(labels = unique(temp$Sample)) +
    theme(legend.position = "none")
  p  
}

# Extended Data Figure 2 c,e
pdf(paste(plot_dir, "iPS_chrX_mito.final.pdf", sep = "/"), height = 2.5, width = 2.5, useDingbats = FALSE)
p1 = iPS_bar("RNR1")
p2 = iPS_hy_bar("RNR1")
grid.arrange(p1,p2, ncol=1)

p1 = iPS_bar("RNR2")
p2 = iPS_hy_bar("RNR2")
grid.arrange(p1,p2, ncol=1)


# Check XIST expression - need the X chromosome for this
geneKeepIPSX = union(geneKeepIPS, geneInfo$gene[which(geneInfo$chrom=="chrX")])
geneKeepIdxIPSX = which(geneInfo$gene %in% geneKeepIPSX)

iPS_tot_keepX = iPS_tot[geneKeepIPSX,]
iPS_tot_rpkX = sweep(as.matrix(iPS_tot_keepX), 1, as.double(geneInfo$genesizeKb[geneKeepIdxIPSX]), `/`)
iPS_tot_tpmX = sweep(as.matrix(iPS_tot_rpkX), 2, as.double(colSums(iPS_tot_rpkX)/1000000), `/`)
dim(iPS_tot_tpmX)

iPS_h_keepX = iPS_h[geneKeepIPSX,]
iPS_h_rpkX = sweep(as.matrix(iPS_h_keepX), 1, as.double(geneInfo$genesizeKb[geneKeepIdxIPSX]), `/`)
iPS_h_tpmX = sweep(as.matrix(iPS_h_rpkX), 2, as.double(colSums(iPS_h_rpkX)/1000000), `/`)
iPS_h_tpmX = iPS_h_tpmX[,seq(13,22)]
dim(iPS_h_tpmX)
head(iPS_h_tpmX)

iPS_c_keepX = iPS_c[geneKeepIPSX,]
iPS_c_rpkX = sweep(as.matrix(iPS_c_keepX), 1, as.double(geneInfo$genesizeKb[geneKeepIdxIPSX]), `/`)
iPS_c_tpmX = sweep(as.matrix(iPS_c_rpkX), 2, as.double(colSums(iPS_c_rpkX)/1000000), `/`)
iPS_c_tpmX = iPS_c_tpmX[,seq(13,22)]
dim(iPS_c_tpmX)
head(iPS_c_tpmX)


iPS_hc_tpmX = cbind(iPS_h_tpmX, iPS_c_tpmX)
dim(iPS_hc_tpmX)

p1 = iPS_bar("XIST", data = iPS_tot_tpmX)
p2 = iPS_hy_bar("XIST", data = iPS_hc_tpmX)
grid.arrange(p1,p2, ncol=1)
dev.off()


## Replication and dimensionality reduction
# Correlations between samples
# Extended Data Figure 3a
pdf(paste(plot_dir, "iPS_Heatmap.rev.pdf", sep = "/")) 
fontsize = 0.75
iPS_cor = cor(iPS_tot_tpm, use="complete.obs")
dim(iPS_cor)
colnames(iPS_cor) = iPS_meta$Line
rownames(iPS_cor) = iPS_meta$Line
col<- colorRampPalette(c("white", "orange", "red"))(200)
par(mai=c(1,0.5,0.5,0.25))
heatmap.2(x = iPS_cor, col = col, symm = FALSE, dendrogram="row", trace="none", density.info="none", keysize = 0.9, srtCol = 65, cexCol = fontsize, cexRow = fontsize)
dev.off()

# PCA
# All genes, removing sex chromosomes and chromosome 20
dim(iPS_tot_cpm)


# Expression threshold
iPS_tot_cpm_comp = iPS_tot_cpm[complete.cases(iPS_tot_cpm),]

iPS_tot_cpm_1 = iPS_tot_cpm_comp[which(rowSums(iPS_tot_cpm_comp)>1),]
dim(iPS_tot_cpm_1)

iPS.pca = prcomp(t(iPS_tot_cpm_1), scale=TRUE, center=TRUE)
summary(iPS.pca)
iPS.pca_data = as.data.frame(iPS.pca$x)

# Figure 1i
pdf(paste(plot_dir, "iPS_PCA.rev.pdf", sep = "/"), height = 2, 2, useDingbats = FALSE) 
fontsize = 8
p1 = ggplot(iPS.pca_data, aes(PC1, PC2)) + 
  geom_point(aes(col=iPS_meta$Species), size = 1)  + 
  ggtitle("PCA") +
  xlab("PC1 (23% of variance)") +
  ylab("PC2 (17% of variance)") +
  theme_bw() +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values = c("dodgerblue2","salmon", "mediumpurple2")) +
  #labs(color = "Species")
  theme(legend.position = "none")
p1 
dev.off()


# PCA with split hybrids and parents
# Parental allelic expression
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

# Expression threshold
iPS_hy_par_comp = iPS_hy_par_noAnSex[complete.cases(iPS_hy_par_noAnSex),]

iPS_hy_par_1 = iPS_hy_par_comp[which(rowSums(iPS_hy_par_comp)>1),]
dim(iPS_hy_par_1)
head(iPS_hy_par_1)

iPS.hy_par.pca = prcomp(t(iPS_hy_par_1), scale=TRUE, center=TRUE)
summary(iPS.hy_par.pca)
iPS.hy_par.pca_data = as.data.frame(iPS.hy_par.pca$x)
iPS_hy_par_species = c(rep("Hy-Chimp", 10), rep("Hy-Human", 10), rep("Human", 6), rep("Chimp", 6))

# Figure 1j
pdf(paste(plot_dir, "iPS_PCA_hyPar.rev.pdf", sep = "/"), height = 2, 2, useDingbats = FALSE) 
fontsize = 8
p = ggplot(iPS.hy_par.pca_data, aes(PC1, PC2)) + 
  geom_point(aes(col=iPS_hy_par_species), size = 1)  + 
  ggtitle("PCA") +
  xlab("PC1 (24% of variance)") +
  ylab("PC2 (13% of variance)") +
  theme_bw() +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values = c("blue4","firebrick4", "deepskyblue", "indianred1")) +
  #labs(color = "Species") +
  theme(legend.position = "none")
p 
dev.off()


# Cis trans plots
# Differential expression in parent iPS data
counts_matrix = iPS_tot[,1:12]
cond_species = c(rep("Human", 6), rep("Chimp", 6))
cond <- ~ cond_species
coldata <- data.frame(cond_species)
dds <- DESeqDataSetFromMatrix(counts_matrix, coldata, cond)
dds <- DESeq(dds)
ips_par_res = results(dds)

dim(ips_par_res)


iPS_h_res = read.table(paste(analysis_dir, "iPS_ASE.DESeq2.GRCh38.v6.txt", sep = "/"), header = TRUE, row.names = 1, fill=TRUE, sep = "\t")
dim(iPS_h_res)

cis = iPS_h_res$log2FoldChange[geneKeepIdxIPS]
trans = ips_par_res$log2FoldChange[geneKeepIdxIPS]

# Figure 1k
pdf(paste(plot_dir, "iPS_cis_trans.rev.pdf", sep = "/"), height = 2, 2, useDingbats = FALSE) 
fontsize = 8
ggplot() +  
  ggtitle("Cis/Trans Test") +
  ylab("Log2(Human / Chimp)") +
  xlab("Log2(Hybrid ASE)") + 
  geom_point(aes(y = trans, x = cis), col = ifelse(rownames(ips_par_res)[geneKeepIdxIPS] == "NRP2", "red", "black"), size = 0.1) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlim(-10,10) +
  ylim(-10,10) +
  theme_bw() +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_abline(slope=1, intercept = 0)
dev.off()

filt = intersect(which(is.finite(trans)), which(is.finite(cis)))
cor(cis[filt], trans[filt]) * cor(cis[filt], trans[filt]) # 0.540

# Percent cis/trans
percentCis = abs(cis)/(abs(cis) + abs(trans-cis))
length(percentCis)
mean(percentCis, na.rm = TRUE) # 55%


# Pluritest
pluri = read.csv(paste(analysis_dir, "PluriTest_Results.csv", sep = "/"), header = TRUE)
head(pluri)

pluri$Line2 = c(1,2,3,1,2,3,1,2,3,1,2,3,1,1,2,2,3,3,4,4,5,5)

# Extended Data Figure 1f
pdf(paste(plot_dir, "pluriTest.rev.pdf", sep = "/"), height = 2.5, 3.5, useDingbats = FALSE) 
fontsize = 8
pointsize = 1
p1 = ggplot(pluri, aes(novelClassic, pluriClassic)) + 
  geom_point(aes(col=Species, shape = Lab), size = pointsize)  + 
  ggtitle("PluriTest") +
  xlab("Novelty") +
  ylab("Pluripotency") +
  geom_hline(yintercept = 20, linetype = "dashed", col = "grey60") +
  geom_vline(xintercept = 1.6, linetype = "dashed", col = "grey60") +
  ylim(-40, 40) +
  xlim(1.3,2.1) +
  theme_bw() +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values = c("dodgerblue2","salmon", "mediumpurple2")) +
  labs(color = "Species", shape = "Line") 
  #theme(legend.position = "none")
p1 
 dev.off()


# PCA with samples from Gilad lab **NOT USED
gilad = read.table(paste(analysis_dir, "HumanChimp_iPSC_Gilad.txt", sep = "/"), header = TRUE, fill=TRUE, sep = "\t")
gilad = gilad[,2:26]
rownames(gilad) = geneInfo$gene

gilad_filt = gilad[which(rownames(gilad) %in% rownames(iPS_tot_cpm)),]
gilad_cpm = sweep(as.matrix(gilad_filt), 2, as.double(colSums(gilad_filt)/1000000), `/`)

gilad_merge_cpm = cbind(gilad_cpm, iPS_tot_cpm)
dim(gilad_merge_cpm)

gilad_merge_cpm_comp = gilad_merge_cpm[complete.cases(gilad_merge_cpm),]

gilad_merge_cpm_1 = gilad_merge_cpm_comp[which(rowSums(gilad_merge_cpm_comp)>1),]
dim(gilad_merge_cpm_1)

gilad.pca = prcomp(t(gilad_merge_cpm_1), scale=TRUE, center=TRUE)
summary(gilad.pca)
gilad.pca_data = as.data.frame(gilad.pca$x)

gilad_merge_species = c(rep("Human", 14), rep("Chimp", 11), iPS_meta$Species)
gilad_merge_study = c(rep("Marchetto", 4), rep("GallegoRomero", 10), rep("Marchetto", 4), rep("GallegoRomero", 7), rep("Agoglia", 22))

pdf(paste(plot_dir, "gilad_PCA.rev.pdf", sep = "/"), height = 2, 2, useDingbats = FALSE) 
fontsize = 8
pointsize = 0.7
p1 = ggplot(gilad.pca_data, aes(PC1, PC2)) + 
  geom_point(aes(col=gilad_merge_species, shape = gilad_merge_study), size = pointsize)  + 
  ggtitle("PCA") +
  xlab("PC1 (26% of variance)") +
  ylab("PC2 (15% of variance)") +
  theme_bw() +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values = c("dodgerblue2","salmon", "mediumpurple2")) +
  #labs(color = "Species", shape = "Study") 
  theme(legend.position = "none")
p1 
dev.off()



