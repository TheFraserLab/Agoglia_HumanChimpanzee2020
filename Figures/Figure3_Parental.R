## Analysis and plots for Figure 3c,f-l and Extended Data Figure 6d-e, 7f,h, 8b-d: hCS and cCS timecourse data
## Everything here was performed in Rstudio on a Macbook Pro with 16G RAM, unless otherwise noted
## Clear cache in Rstudio
rm(list = ls(all.names = TRUE))

## Load required libraries
library(ggplot2)
library(gridExtra)
library(gplots)
library(dplyr)
library(reshape)
library(scales)
library(DESeq2)
#library(plyr) # will need to load and unload this later, as it interferes with certain dplyr features

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

# Timecourse RNA seq data (counts per sample, performed in two rounds)
# Human
hCS_200 = read.table(paste(data_dir, "hCS_200.GRCh38.txt", sep = "/"), header = TRUE)
hCS_200_meta = read.table(paste(data_dir, "hCS_200.meta.txt", sep = "/"), header = TRUE)
dim(hCS_200)

hCS_35 = read.table(paste(data_dir, "hCS_35.GRCh38.txt", sep = "/"), header = TRUE)
hCS_35_meta = read.table(paste(data_dir, "hCS_35.meta.txt", sep = "/"), header = TRUE)
dim(hCS_35)

# Merge the two timecourses into one
hCS = as.data.frame(cbind(hCS_35, hCS_200))
rownames(hCS) = geneInfo$gene
dim(hCS)
peek(hCS)
hCS_tot = hCS[, seq(1,length(hCS),by=3)]

hCS_meta = as.data.frame(rbind(hCS_35_meta, hCS_200_meta))
dim(hCS_meta)
head(hCS_meta)
hCS_meta$Species = as.character(hCS_meta$Species)
hCS_meta$Line = as.character(hCS_meta$Line)

# Chimp
cCS_200 = read.table(paste(data_dir, "cCS_200.GRCh38.txt", sep = "/"), header = TRUE)
cCS_200_meta = read.table(paste(data_dir, "cCS_200.meta.txt", sep = "/"), header = TRUE)
dim(cCS_200)

cCS_35 = read.table(paste(data_dir, "cCS_35.GRCh38.txt", sep = "/"), header = TRUE)
cCS_35_meta = read.table(paste(data_dir, "cCS_35.meta.txt", sep = "/"), header = TRUE)
dim(cCS_35)

# Merge the two timecourses into one
cCS = as.data.frame(cbind(cCS_35, cCS_200))
rownames(cCS) = geneInfo$gene
dim(cCS)
peek(cCS)
cCS_tot = cCS[, seq(1,length(cCS),by=3)]

cCS_meta = as.data.frame(rbind(cCS_35_meta, cCS_200_meta))
dim(cCS_meta)
head(cCS_meta)
cCS_meta$Species = as.character(cCS_meta$Species)
cCS_meta$Line = as.character(cCS_meta$Line)

# DESeq2 for cell fate plots and comparison with hyCS data later
deseqParent <- function(counts_matrix, cond_species, cond_induction){
  require(DESeq2)
  cond <- ~ cond_induction + cond_species
  cond_reduced <- ~ cond_induction 
  coldata <- data.frame(cond_induction, cond_species)
  dds <- DESeqDataSetFromMatrix(counts_matrix, coldata, cond)
  dds <- DESeq(dds,test="LRT",betaPrior=FALSE,full=cond,reduced=cond_reduced)
  res <- results(dds)
  res
}

d15_par_res = deseqParent(cbind(hCS_tot[,1:6], cCS_tot[,1:6]),
                       c(rep("Human", 6), rep("Chimp", 6)),
                       c(as.character(hCS_meta$Induction[1:6]), as.character(cCS_meta$Induction[1:6])))

d25_par_res = deseqParent(cbind(hCS_tot[,7:12], cCS_tot[,7:12]),
                       c(rep("Human", 6), rep("Chimp", 6)),
                       c(as.character(hCS_meta$Induction[7:12]), as.character(cCS_meta$Induction[7:12])))

d35_par_res = deseqParent(cbind(hCS_tot[,13:18], cCS_tot[,13:18]),
                       c(rep("Human", 6), rep("Chimp", 6)),
                       c(as.character(hCS_meta$Induction[13:18]), as.character(cCS_meta$Induction[13:18])))

d50_par_res = deseqParent(cbind(hCS_tot[,19:24], cCS_tot[,19:24]),
                       c(rep("Human", 6), rep("Chimp", 6)),
                       c(as.character(hCS_meta$Induction[19:24]), as.character(cCS_meta$Induction[19:24])))

d100_par_res = deseqParent(cbind(hCS_tot[,25:30], cCS_tot[,25:29]),
                        c(rep("Human", 6), rep("Chimp", 5)),
                        c(as.character(hCS_meta$Induction[25:30]), as.character(cCS_meta$Induction[25:29])))

d150_par_res = deseqParent(cbind(hCS_tot[,31:36], cCS_tot[,30:34]),
                        c(rep("Human", 6), rep("Chimp", 5)),
                        c(as.character(hCS_meta$Induction[31:36]), as.character(cCS_meta$Induction[30:34])))

d200_par_res = deseqParent(cbind(hCS_tot[,37:41], cCS_tot[,35:39]),
                        c(rep("Human", 5), rep("Chimp", 5)),
                        c(as.character(hCS_meta$Induction[37:41]), as.character(cCS_meta$Induction[35:39])))

merged_par_res = deseqParent(cbind(hCS_tot, cCS_tot),
                          c(rep("Human", 41), rep("Chimp", 39)),
                          c(as.character(hCS_meta$Induction), as.character(cCS_meta$Induction)))

output_name = paste(analysis_dir, "D150_par_DESeq2.rev.txt", sep = "/")

write.table(data.frame(d150_par_res), 
            file=output_name, 
            sep="\t", quote = FALSE,na = "", row.names = TRUE, col.names = TRUE)

output_name = paste(analysis_dir, "Merged_par_DESeq2.rev.txt", sep = "/")
write.table(data.frame(merged_par_res), 
            file=output_name, 
            sep="\t", quote = FALSE,na = "", row.names = TRUE, col.names = TRUE)

# Supplementary TableS2
output_name = paste(analysis_dir, "All_par_DESeq2.rev.txt", sep = "/")
all_par_res = cbind(d15_par_res$log2FoldChange, d15_par_res$padj, 
                    d25_par_res$log2FoldChange, d25_par_res$padj,
                    d35_par_res$log2FoldChange, d35_par_res$padj, 
                    d50_par_res$log2FoldChange, d50_par_res$padj, 
                    d100_par_res$log2FoldChange, d100_par_res$padj, 
                    d150_par_res$log2FoldChange, d150_par_res$padj, 
                    d200_par_res$log2FoldChange, d200_par_res$padj)
rownames(all_par_res) = rownames(d15_par_res)

dim(all_par_res)
colnames(all_par_res) = c("Day15_LFC", "Day15_padj", 
                          "Day25_LFC", "Day25_padj", 
                          "Day35_LFC", "Day35_padj", 
                          "Day50_LFC", "Day50_padj", 
                          "Day100_LFC", "Day100_padj", 
                          "Day150_LFC", "Day150_padj", 
                          "Day200_LFC", "Day200_padj")
peek(all_par_res)
write.table(data.frame(all_par_res), 
            file=output_name, 
            sep="\t", quote = FALSE,na = "", row.names = FALSE, col.names = TRUE)


# Use the list of biased genes from the hyCS data and remove the same ones here
biased_genes = read.table( paste(analysis_dir, "all_biased_genes.txt", sep = "/"), header = TRUE)
dim(biased_genes)

par_unbiased = which(!(rownames(hCS_tot) %in% biased_genes$all_biased_genes))
length(par_unbiased)

# Calculate TPM and CPM
hCS_tot_keep = hCS_tot[par_unbiased,]
hCS_tot_cpm = sweep(as.matrix(hCS_tot_keep), 2, as.double(colSums(hCS_tot_keep)/1000000), `/`)
hCS_tot_fpkm = sweep(as.matrix(hCS_tot_cpm), 1, as.double(geneInfo$genesizeKb[par_unbiased]), `/`)
hCS_tot_rpk = sweep(as.matrix(hCS_tot_keep), 1, as.double(geneInfo$genesizeKb[par_unbiased]), `/`)
hCS_tot_tpm = sweep(as.matrix(hCS_tot_rpk), 2, as.double(colSums(hCS_tot_rpk)/1000000), `/`)

# Calculate TPM and CPM
cCS_tot_keep = cCS_tot[par_unbiased,]
cCS_tot_cpm = sweep(as.matrix(cCS_tot_keep), 2, as.double(colSums(cCS_tot_keep)/1000000), `/`)
cCS_tot_fpkm = sweep(as.matrix(cCS_tot_cpm), 1, as.double(geneInfo$genesizeKb[par_unbiased]), `/`)
cCS_tot_rpk = sweep(as.matrix(cCS_tot_keep), 1, as.double(geneInfo$genesizeKb[par_unbiased]), `/`)
cCS_tot_tpm = sweep(as.matrix(cCS_tot_rpk), 2, as.double(colSums(cCS_tot_rpk)/1000000), `/`)


# Write files for CIBERSORT ("mixture" file; see the Rscript for figure 3 for generation of the "pheno_classes" and "reference" files for this analysis)
write.table(file = paste(analysis_dir, "hCS.ciber_mixture_011320.txt", sep = "/"), x=hCS_tot_fpkm, sep = "\t", quote=FALSE)
write.table(file = paste(analysis_dir, "cCS.ciber_mixture_011320.txt", sep = "/"), x=cCS_tot_fpkm, sep = "\t", quote=FALSE)


# Heatmap
hlabels = paste0(paste0(hCS_meta$Line, rep("_D", 41), sep = ""), hCS_meta$Age, sep = "")
clabels = paste0(paste0(cCS_meta$Line, rep("_D", 39), sep = ""), cCS_meta$Age, sep = "")

# Separate for each species
# Extended Data Figure 6d,e
pdf(paste(plot_dir, "ParentalTimecourseHeatmap.rev.pdf", sep = "/")) 
fontsize = 0.75

hCS_cor = cor(hCS_tot_fpkm, use="complete.obs")
dim(hCS_cor)
colnames(hCS_cor) = hlabels
rownames(hCS_cor) = hlabels
col<- colorRampPalette(c("turquoise4", "white", "darkorange2"))(200)
par(mai=c(1,0.5,0.5,0.25))
heatmap.2(x = hCS_cor, col = col, srtCol = 55, symm = FALSE, dendrogram="row", trace="none", density.info="none", keysize = 0.9, breaks = seq(0.5, 1, length.out=201), cexCol = fontsize, cexRow = fontsize)

cCS_cor = cor(cCS_tot_fpkm, use="complete.obs")
dim(cCS_cor)
colnames(cCS_cor) = clabels
rownames(cCS_cor) = clabels
col<- colorRampPalette(c("turquoise4", "white", "darkorange2"))(200)
par(mai=c(1,0.5,0.5,0.25))
heatmap.2(x = cCS_cor, col = col, srtCol = 55, symm = FALSE, dendrogram="row", trace="none", density.info="none", keysize = 1.0, breaks = seq(0.5, 1, length.out=201), cexCol = fontsize, cexRow = fontsize)
dev.off()

# median correlation values
median(hCS_cor*hCS_cor) # 0.87
median(cCS_cor*cCS_cor) # 0.80

# medians without outliers
rowMedians(hCS_cor)
rowMedians(cCS_cor)

hCS_keep = c(seq(1,26), seq(28,41))
median(hCS_cor[hCS_keep,hCS_keep]*hCS_cor[hCS_keep,hCS_keep]) # 0.88

cCS_keep = c(seq(1,34), 36, 39)
median(hCS_cor[cCS_keep,cCS_keep]*hCS_cor[cCS_keep,cCS_keep]) # 0.89


## Dimentionality reduction and clustering
# PCA of all samples in the time course
hcCS_tot_cpm = as.data.frame(cbind(hCS_tot_cpm, cCS_tot_cpm))
dim(hcCS_tot_cpm)

# Save for bar plots in figure 4
write.table(file = paste(analysis_dir, "hcCS_tot_cpm.rev.txt", sep = "/"), x = hcCS_tot_cpm, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")

hcCS_species = c(hCS_meta$Species, cCS_meta$Species)
hcCS_time = c(hCS_meta$Age, cCS_meta$Age)

# Exclude any genes that are missing
hcCS_tot_cpm_comp = hcCS_tot_cpm[complete.cases(hcCS_tot_cpm),]

# Require at least one cpm per gene
hcCS_tot_cpm_1 = hcCS_tot_cpm_comp[which(rowSums(hcCS_tot_cpm_comp)>1),]
dim(hcCS_tot_cpm_1)

hcCS.pca = prcomp(t(hcCS_tot_cpm_1), scale=TRUE, center=TRUE)
summary(hcCS.pca)
hcCS.pca_data = as.data.frame(hcCS.pca$x)


# ID top genes contributing to PC1, plot these across samples
hcCS_species = c(hCS_meta$Species, cCS_meta$Species)
hcCS_time = c(hCS_meta$Age, cCS_meta$Age)

hc_bar <- function(gene){
  temp = data.frame(as.matrix(t(rbind(hcCS_tot_cpm[which(rownames(hcCS_tot_cpm)==gene),], hcCS_species, hcCS_time))))
  names(temp) = c("exp", "species", "timepoint")
  temp$exp = as.numeric(as.character(temp$exp))
  temp$timepoint = factor(temp$timepoint, levels = c(15,25,35,50,100, 150, 200), ordered = TRUE)
  temp = temp[order(temp$timepoint),]
  temp2 = temp %>% group_by(species, timepoint) %>% summarise(mean_exp = mean(exp), se_exp = sd(exp)/sqrt(n()))
  ymax = max(temp$exp) + max(temp2$se_exp)
  p = ggplot(temp2, aes(x=factor(timepoint), y=mean_exp, fill = species)) + 
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_errorbar(aes(ymin = mean_exp - se_exp, ymax = mean_exp + se_exp), width=0.2, position = position_dodge(0.7))+
    geom_point(data = temp, aes(x = factor(timepoint), y = exp, fill = species), size = 0.9, stat = "identity", position_dodge(width = 0.7)) +
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("dodgerblue3", "salmon")) +
    ggtitle(gene)+
    theme(plot.title = element_text(size = fontsize)) +
    theme(legend.position = "none") +
    ylim(0, ymax) +
    theme(axis.text=element_text(size=fontsize), axis.text.x=element_text(angle = 45, hjust =1), axis.title=element_blank())
  p
}

peek(hcCS.pca$rotation)
top_PC1 = names(hcCS.pca$rotation[,1][order(abs(hcCS.pca$rotation[,1]),decreasing=TRUE)][1:20])

for(i in top_PC1) {
  cat(i)
  print(hc_bar(i))
}

# CAPZA1, G3BP1** associated with astrocytoma, PPP4R2, SNRPE


top_PC2 = names(hcCS.pca$rotation[,2][order(abs(hcCS.pca$rotation[,2]),decreasing=TRUE)][1:20])
top_PC2

for(i in top_PC2) {
  cat(i)
  print(hc_bar(i))
}

# For GO enrichment
PC1_ranked = names(hcCS.pca$rotation[,2][order(hcCS.pca$rotation[,2],decreasing=TRUE)])
PC2_ranked = names(hcCS.pca$rotation[,2][order(hcCS.pca$rotation[,2],decreasing=TRUE)])

write.table(x = cbind(PC1_ranked, PC2_ranked), file = paste(analysis_dir, "PC1_PC2_ranked.txt", sep = "/"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# PC1: Top of list (genes that decrease over time) are enriched for cell cycle; 
# bottom of list (increase over time) are enriched for synaptic activity / neuron function in general /
# nervous system development




# Heat map of marker genes' expression across time in both species
markerMap <- function(marker_genes){
  require(plyr)
  fontsize = 8.0
  h_mark = data.frame(t(hCS_tot_tpm[marker_genes,]))
  
  h_mark$Timepoint = hCS_meta$Age
  h_mark_summary = group_by(h_mark, Timepoint)
  h_mark_summary = summarise_all(h_mark_summary, mean)
  
  dfh = as.data.frame(h_mark_summary)[,2:length(h_mark_summary)]
  row.names(dfh) = h_mark_summary$Timepoint
  
  # Chimp
  c_mark = data.frame(t(cCS_tot_tpm[marker_genes,]))
  
 # c_time = hc_cib$Timepoint[42:80]
  c_mark$Timepoint = cCS_meta$Age
  c_mark_summary = group_by(c_mark, Timepoint)
  c_mark_summary = summarise_all(c_mark_summary, mean)
  
  dfc = as.data.frame(c_mark_summary)[,2:length(c_mark_summary)]
  row.names(dfc) = c_mark_summary$Timepoint
  
  dfhm = melt(cbind(dfh, ind = rownames(dfh)))
  dfhm <- ddply(dfhm, .(variable), transform, rescale = rescale(log10(value)))
  dfhm$ind <- factor(dfhm$ind, levels = c(15,25,35,50,100,150,200))
    
  dfcm = melt(cbind(dfc, ind = rownames(dfc)))
  dfcm <- ddply(dfcm, .(variable), transform, rescale = rescale(log10(value)))
  dfcm$ind <- factor(dfcm$ind, levels = c(15,25,35,50,100,150,200))
    
  low = min(min(dfhm$rescale), min(dfcm$rescale)) - 0.1
  high = max(max(dfhm$rescale),max(dfcm$rescale)) + 0.1
    
  p1 =ggplot(dfhm, aes(factor(ind), variable, width = 3, height = 1)) + 
      geom_tile(aes(fill = rescale), color = "black", height = 1, width=1) +
      scale_fill_gradient(low = "white", high = "darkgreen", limits = c(low, high)) +
      theme_bw()+
      ylab("") +
      theme(axis.ticks.x=element_blank()) +
      xlab("Timepoint (Days)")+
      theme(axis.title.y=element_blank(),
            axis.ticks.y=element_blank()) +
      #labs(fill = "TPM") +
      theme(legend.position = "none", legend.box = "horizontal") +
      theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize)) + 
      #theme(legend.text=element_text(size=fontsize))+
      theme(plot.title = element_text(hjust = 0.5, size=fontsize)) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank())
    
  p2 = ggplot(dfcm, aes(factor(ind), variable), width = 3, height = 1) + 
      geom_tile(aes(fill = rescale), color = "black", height = 1, width=1) +
      scale_fill_gradient(low = "white", high = "darkgreen", limits = c(low, high)) +
      theme_bw()+
      ylab("") +
      theme(axis.ticks.x=element_blank()) +
      xlab("Timepoint (Days)")+
      theme(axis.title.y=element_blank(),
            axis.ticks.y=element_blank()) +
      labs(fill = "Row-scaled TPM") +
      theme(legend.position = "top", legend.box = "horizontal") +
      #theme(legend.position = "none") +
      theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize)) + 
      theme(legend.text=element_text(size=fontsize))+
      theme(plot.title = element_text(hjust = 0.5, size=fontsize))+
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank())
    
  grid.arrange(p1,p2, ncol =1)
}

# Extended Data Figure 7f,h
pdf(paste(plot_dir, "ParentalTimecourseMarkerGenes.rev.pdf", sep = "/"), height = 4, width = 5)

# Astrogenesis
markerMap(c("AQP4", "HEPACAM", "GFAP", "AGT"))

# Cell division
markerMap(c("MKI67", "TOP2A", "NEK2", "ASPM"))

# Neurogenesis
markerMap(c("TUBB3", "MLLT11", "MAP2", "RBFOX3"))
dev.off()

# Differences
all_par_res = read.table(paste(analysis_dir, "All_par_DESeq2.rev.txt", sep = "/"), header = TRUE, fill = TRUE, sep = "\t")
rownames(all_par_res) = geneInfo$gene
markers = c("FOXG1", "EMX1", "EOMES", "TBR1", "DCN", "LUM", "COL1A1", "COL1A2")
df = t(all_par_res[markers, c(1,3,5,7,9,11,13)])
dim(df)

row.names(df) = c(15,25,35,50,100,150,200)
df
dfm = melt(cbind(df))
dfm
colnames(dfm) = c("Day", "Gene", "LFC")
dfm$Day <- factor(dfm$Day, levels = c(15,25,35,50,100,150,200))
#dfm$Gene <- factor(dfm$Gene, levels = c("COL1A1","COL1A2", "COL3A1","BGN", "LUM","DCN","PAX6","SIX3", "EMX1", "FOXG1"))
dfm$Gene <- factor(dfm$Gene, levels = c("FOXG1", "EMX1", "EOMES", "TBR1","DCN", "LUM", "COL1A2", "COL1A1"))

# Figure 3c
pdf(paste(plot_dir, "ParentMarkersHeatmap.rev.pdf", sep = "/"), height = 1.35, width = 1.5) 
fontsize = 8
p= ggplot(dfm, aes(factor(Day), Gene)) + 
  geom_tile(aes(fill = LFC), color = "black", height = 1, width=1) +
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


# NOT USED
pdf(paste(plot_dir, "ParentalTimecourseMarkerGenes_TimeDiffsAstroNeuron.pdf", sep = "/"), height = 1, width = 2)
markers = c("AQP4", "HEPACAM", "GFAP", "AGT", "TUBB3", "MLLT11", "MAP2", "RBFOX3")
df = t(all_par_res[markers, c(1,3,5,7,9,11,13)])
dim(df)

row.names(df) = c(15,25,35,50,100,150,200)
df
dfm = melt(cbind(df))
dfm
colnames(dfm) = c("Day", "Gene", "LFC")
dfm$Day <- factor(dfm$Day, levels = c(15,25,35,50,100,150,200))
dfm$Gene <- factor(dfm$Gene, levels = c("AQP4", "HEPACAM", "GFAP", "AGT", "TUBB3", "MLLT11", "MAP2", "RBFOX3"))
fontsize = 8
p= ggplot(dfm, aes(factor(Day), Gene)) + 
  geom_tile(aes(fill = LFC), color = "black", height = 1, width=1) +
  scale_fill_gradientn(colours = c("dodgerblue4", "white", "indianred3"), limits = c(-4,4)) +
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


## NOT USED
# Enrichment tests for differential expression at day 25
msigdb = read.table("/Users/rachelagoglia/Desktop/Manuscript/AnalysisFiles/MSigDB_BP_MF_KEGG.txt", header = TRUE, fill = TRUE)
msigdb = msigdb[-c(1),]
msigdb = as.data.frame(msigdb)
dim(msigdb)
peek(msigdb)

hcCS_de = read.table(paste(analysis_dir, "all_par_DESeq2.txt", sep = "/"), header = TRUE, row.names = 1, fill = TRUE, sep = "\t")
head(hcCS_de)

hcCS_de_d25 = hcCS_de$Day25_LFC

datad25 = hcCS_de_d25[which(!is.na(hcCS_de_d25))]
rowsd25 = rownames(hcCS_de)[which(!is.na(hcCS_de_d25))]

d25_cats = c()
d25_pvals = c()
for(i in seq(1, length(msigdb))){
  genes = as.vector(msigdb[,i])[which(msigdb[,i] != "")]
  test = datad25[which(rowsd25 %in% genes)]
  background = datad25[which(!(rowsd25 %in% genes))]
  if(length(test) >= 20){
    x = wilcox.test(test,background)
    d25_cats = c(d25_cats, colnames(msigdb)[i])
    d25_pvals = c(d25_pvals, x$p.value)
  }
}
length(d25_cats)
length(d25_pvals) 

min(d25_pvals)

pvalBar<-function(cats, pvals, cutoff = 0.01, ncats = 10, color = "grey50"){
  df = as.data.frame(cbind(cats, pvals))
  df2 = df[order(df$pvals),]
  topcats = as.data.frame(df2[1:ncats,])
  topcats = topcats[which(as.numeric(as.character(topcats$pvals)) < cutoff),]
  
  neglog = -log10(as.numeric(as.character(topcats$pvals)))
  df3 = as.data.frame(cbind(as.character(topcats$cats), neglog))
  df3$V1 = factor(df3$V1, levels = df3$V1, ordered = TRUE)
  print(df3)
  
  p = ggplot(data=df3, aes(x=as.factor(V1),y=as.numeric(as.character(neglog)))) +
    geom_bar(position="dodge",stat="identity", color = "black", fill = color, width = 0.7) + 
    coord_flip() +
    theme(title = element_blank())+
    xlab("") + 
    ylab("-Log10(p-value)") +
    ylim(0,4) +
    theme(axis.text.y=element_text(size=fontsize), axis.title.y=element_text(size=fontsize)) +
    theme_bw() + 
    theme(legend.position = "none") +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank()) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  p
}

pvalBar(d25_cats, d25_pvals, color = "grey60")






# WGCNA
library("WGCNA")
disableWGCNAThreads()
options(stringsAsFactors = FALSE)

hcCS_tot_fpkm = cbind(hCS_tot_fpkm, cCS_tot_fpkm)
dim(hcCS_tot_fpkm)

net_hcCS = as.data.frame(t(hcCS_tot_fpkm))
net_gsg = goodSamplesGenes(net_hcCS, verbose = 3)
net_gsg$allOK

# Remove the bad genes / samples
net_hcCS_filt = net_hcCS[net_gsg$goodSamples, net_gsg$goodGenes]
dim(net_hcCS_filt)

# Cluster to ID outliers
sampleTree = hclust(dist(net_hcCS_filt), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
 

# Remove outlier sample
# Plot a line to show the cut
abline(h = 5800, col = "red")
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 5800, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
net_hcCS_filt = net_hcCS_filt[keepSamples, ]

nGenes = ncol(net_hcCS_filt)
nSamples = nrow(net_hcCS_filt)


# Choose soft thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(net_hcCS_filt, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.80,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


write.table(x = net_hcCS_filt, file = paste(analysis_dir,  "net_hcCS_filt.txt", sep = "/"), row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# Gene network construction // module identification
net = blockwiseModules(net_hcCS_filt, power = 8,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespechcCSDendro = FALSE,
                       saveTOMs = TRUE,
                       maxBlockSize = 50000,
                       saveTOMFileBase = "hcCS_TOM",
                       verbose = 3)

table(net$colors)
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)

plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


dev.off()

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "hcCS_moduleInfo.RData")

net_hcCS_filt = read.table(file = paste(analysis_dir,  "net_hcCS_filt.test.txt", sep = "/"), row.names = 1, header = TRUE,sep = "\t")

hcCS_moduleInfo = load(file = paste(analysis_dir, "hcCS_moduleInfo.RData",sep = "/"))
head(MEs)
head(moduleLabels)
head(moduleColors)

length(moduleLabels)
length(moduleColors)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(net_hcCS_filt, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = c()

hcCS_meta = rbind(hCS_meta, cCS_meta)
hcCS_meta_filt = hcCS_meta[which(rownames(net_hcCS) %in% rownames(net_hcCS_filt)),]

# Table of all genes and what module they are in
write.table(x = cbind(names(net_hcCS_filt), moduleColors), file = paste(analysis_dir, "hcCS_wgcna_geneModules.txt", sep = "/"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# ID modules that have a species difference across time points
species_pvals = c()
for(i in names(MEs0)){
  eigengenes = MEs0[,i]
  human_eig = eigengenes[1:41]
  chimp_eig = eigengenes[42:73]
  species_pvals = c(species_pvals, t.test(human_eig, chimp_eig)$p.value)
}

length(species_pvals)
min(species_pvals)
cat(paste0(paste(names(MEs0), species_pvals, sep = "\t"), sep = "\n"))

modulePlot <-function(module){
  data = data.frame(hcCS_meta_filt$Age, hcCS_meta_filt$Species, MEs0[,module])
  names(data) = c("Age", "Species", "Module")
  data$Module = as.numeric(as.character(data$Module))
  data$Age = factor(data$Age, levels = c(15,25,35,50,100, 150, 200), ordered = TRUE)
  data = data[order(data$Age),]
  data2 = data %>% group_by(Species, Age) %>% summarise(mean_eig = mean(Module), se_eig = sd(Module)/sqrt(n()), median_eig = median(Module))
  p = ggplot() + 
    geom_boxplot(data = data, aes(x=factor(Age), y=Module, fill = Species), position = "dodge", width = 0.7, color = "black", size = 0.2, outlier.shape=NA) +
    geom_smooth(data = as.data.frame(data2[which(data2$Species == "Human"),]), aes(x=as.numeric(factor(Age)), y=median_eig), colour = "indianred3", size = 0.5, group = 1,linetype = "dashed", se=FALSE) +
    geom_smooth(data = as.data.frame(data2[which(data2$Species == "Chimp"),]), aes(x=as.numeric(factor(Age)), y=median_eig), colour = "dodgerblue4", size = 0.5, group = 1, linetype = "dashed", se=FALSE) +
    geom_point(data = data, aes(x = factor(Age), y = Module, fill = Species), size = pointsize, stat = "identity", position_dodge(width = 0.7)) +
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("dodgerblue3", "salmon")) +
    ggtitle(module)+
    theme(plot.title = element_text(size = fontsize)) +
    theme(legend.position = "none") +
    theme(axis.text=element_text(size=fontsize), axis.text.x=element_text(angle = 45, hjust =1), axis.title=element_blank()) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank()) 
  p
}

# Ones with a significant difference across all samples (t-test)
require("dplyr")
species_mods = names(MEs0)[which(species_pvals<0.05)]
for(i in species_mods) {
  print(modulePlot(i))
}


# Correlations with time
moduleTraitCor = cor(MEs, hcCS_meta_filt, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)


labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(hcCS_meta_filt),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

timeMods = rownames(as.data.frame(moduleTraitPvalue))[which(as.data.frame(moduleTraitPvalue)$Age<0.05)]
for(i in timeMods) {
  print(modulePlot(i))
}

modulePlot("MEbisque4")
modulePlot("MEblack")
modulePlot("MEblue")
modulePlot("MEbrown")
modulePlot("MEdarkgrey")
modulePlot("MEmagenta")
modulePlot("MEmidnightblue")
modulePlot("MEred")
modulePlot("MEroyalblue")
modulePlot("MEsalmon")
modulePlot("MEturquoise")
modulePlot("MEyellow")
modulePlot("MEskyblue")
modulePlot("MEcyan")
modulePlot("MEbrown4")

# Plot the main 3 modules
# Extended Data Figure 8b
pointsize = 0.2
pdf(paste(plot_dir, "ParentalWGCNA_Modules_processes.pdf", sep = "/"), height = 4, width = 3, useDingbats = FALSE)
p1 = modulePlot("MEblue")
p2 = modulePlot("MEbrown")
p3 = modulePlot("MEred")
grid.arrange(p1,p2,p3, ncol=1)
dev.off()

# Figure 3f
pointsize = 0.2
pdf(paste(plot_dir, "ParentalWGCNA_Modules_extra.pdf", sep = "/"), height = 3, width = 3, useDingbats = FALSE)
p1 = modulePlot("MEmagenta")
p2 = modulePlot("MEsalmon")
grid.arrange(p1,p2, ncol=1)
dev.off()


# Now look in hyCS
geneKeepInfo = read.table(file = paste(analysis_dir, "geneKeep.rev.txt", sep = "/"), header = TRUE, fill = TRUE)
geneKeep = as.vector(geneKeepInfo$geneKeep)
geneKeepIdx = as.vector(geneKeepInfo$geneKeepIdx)



hyCS = read.table(paste(data_dir, "hyCS.GRCh38.txt", sep = "/"), header = TRUE)
rownames(hyCS) = geneInfo$gene
dim(hyCS)
hyCS_tot = hyCS[, seq(1,75,by=3)]
hyCS_meta = read.table(paste(data_dir, "hyCS.meta.txt", sep = "/"), header = TRUE)
hyCS_meta$Line = as.character(hyCS_meta$Line)
hyCS_meta$Induction = as.character(hyCS_meta$Induction)

hyCS_hc = hyCS[, c(seq(2, length(hyCS), 3), seq(3, length(hyCS), 3))]
dim(hyCS_hc)
hyCS_hc_keep = hyCS_hc[geneKeep,]
hyCS_hc_cpm = sweep(as.matrix(hyCS_hc_keep), 2, as.double(colSums(hyCS_hc)/1000000), `/`)




hy_mod_box <-function(mod){
  mod_genes = colnames(net_hcCS_filt)[which(moduleColors == mod)]
  hyCS_mod = hyCS_hc_cpm[which(rownames(hyCS_hc_cpm) %in% mod_genes),]
  hyCS_mod_comp = hyCS_mod[complete.cases(hyCS_mod),]
  hyCS_mod_1 = hyCS_mod_comp[which(rowSums(hyCS_mod_comp)>1),]
  mod.pca = prcomp(t(hyCS_mod_1), scale=TRUE, center=TRUE)
  mod.pca_data = as.data.frame(mod.pca$x)

  allele = c(rep("Human", 25), rep("Chimp", 25))
  age = rep(hyCS_meta$Age, 2)
  data = data.frame(age, mod.pca_data$PC1, allele)
  names(data) = c("Age", "Eigengene", "Allele")
  data$Eigengene = as.numeric(as.character(data$Eigengene))
  
  # Re-orient eigengenes for consistency with WGCNA eigengenes
  if(cor(colMeans(hyCS_mod_1), data$Eigengene)[1] < 0 ){
    data$Eigengene = (-1)*data$Eigengene
  }
  
  data$Age = factor(data$Age, levels = c(50,100, 150, 200), ordered = TRUE)
  data = data[order(data$Age),]
  data2 = data %>% group_by(Allele, Age) %>% summarise(mean_eig = mean(Eigengene), se_eig = sd(Eigengene)/sqrt(n()), median_eig = median(Eigengene))
  ymax = max(data$Eigengene) + max(data2$se_eig)
  ymin = min(data$Eigengene) - mean(data2$se_eig)
  p = ggplot() + 
    geom_boxplot(data = data, aes(x=factor(Age), y=Eigengene, fill = Allele), position = "dodge", width = 0.7, color = "black", size = 0.2, outlier.shape = NA) +
    geom_smooth(data = as.data.frame(data2[which(data2$Allele == "Human"),]), aes(x=as.numeric(factor(Age)), y=median_eig), colour = "indianred3", size = 0.5, group = 1,linetype = "dashed", se=FALSE) +
    geom_smooth(data = as.data.frame(data2[which(data2$Allele == "Chimp"),]), aes(x=as.numeric(factor(Age)), y=median_eig), colour = "dodgerblue4", size = 0.5, group = 1, linetype = "dashed", se=FALSE) +
    geom_point(data = data, aes(x = factor(Age), y = Eigengene, fill = Allele), size = pointsize, stat = "identity", position_dodge(width = 0.7)) +
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("dodgerblue3", "salmon")) +
    ggtitle(mod)+
    theme(plot.title = element_text(size = fontsize)) +
    theme(legend.position = "none") +
    theme(axis.text=element_text(size=fontsize), axis.text.x=element_text(angle = 45, hjust =1), axis.title=element_blank()) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank()) 
  p
}

# Extended Data Figure 8d
pdf(paste(plot_dir, "hyCS_allelic_WGCNA_Modules_processes.pdf", sep = "/"), height = 4, width = 1.7, useDingbats = FALSE)
p1 = hy_mod_box("blue")
p2 = hy_mod_box("brown")
p3 = hy_mod_box("red")
grid.arrange(p1,p2,p3,ncol=1)
dev.off()

# Figure 3g
pdf(paste(plot_dir, "hyCS_allelic_WGCNA_Modules_extra.pdf", sep = "/"), height = 2.5, width = 1.7, useDingbats = FALSE)
p4 = hy_mod_box("magenta")
p5 = hy_mod_box("salmon")
grid.arrange(p4,p5,ncol=1)
dev.off()

pdf(paste(plot_dir, "hyCS_allelic_WGCNA_Modules_ALL.pdf", sep = "/"), height = 4, width = 4, useDingbats = FALSE)
for(i in unique(moduleColors)) {
  print(hy_mod_box(i))
}
dev.off()


# Modules plotted in single cell data
total_filt = read.table(file=paste(analysis_dir, "total_filt.rev.txt", sep = "/"), header = TRUE)
umap = read.table(file=paste(analysis_dir, "sc_umap_coords.rev.txt", sep = "/"), header = TRUE)

sc_mod <-function(mod){
  mod_genes = colnames(net_hcCS_filt)[which(moduleColors == mod)]
  sc_mod = total_filt[which(rownames(total_filt) %in% mod_genes),]
  sc_mod_comp = sc_mod[complete.cases(sc_mod),]
  sc_mod_1 = sc_mod_comp[which(rowSums(sc_mod_comp)>1),]
  mod.pca = prcomp(t(sc_mod_1), scale=TRUE, center=TRUE)
  mod.pca_data = as.data.frame(mod.pca$x)
  eig = mod.pca_data$PC1
  
  # Re-orient eigengenes for consistency with WGCNA eigengenes
  if(cor(colMeans(sc_mod_1), eig)[1] < 0 ){
    eig = (-1)*eig
  }
  
  coords = umap
  col_by = eig
  bplot = ggplot() +
    ggtitle(mod) +
    theme(plot.title = element_text(size=fontsize)) +
    geom_point(aes(x = coords[,1], y = coords[,2], color = col_by), size = pointsize) + 
    scale_color_gradient(low = "grey90", high = "maroon4", na.value = "grey80") +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    theme(axis.title.y =element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    theme(legend.text=element_text(size=fontsize - 3),legend.title=element_text(size=fontsize - 3)) +
    theme(legend.key.size = unit(0.1, "in")) +
    labs(color = "Exp")
  #theme(title = element_blank()) 
  bplot
}

# Extended Data Figure 8c
pdf(paste(plot_dir, "sc_WGCNA_Modules.pdf", sep = "/"), height = 2, width =2, useDingbats = FALSE)
fontsize = 8
pointsize = 0.1
sc_mod("blue")
sc_mod("brown")
sc_mod("red")
sc_mod("magenta")
sc_mod("salmon")
dev.off()

# Now with neural cells only
total_filt = read.table(file=paste(analysis_dir, "total_filt.rev.txt", sep = "/"), header = TRUE)
umap = read.table(file=paste(analysis_dir, "neural_sc_umap_coords.rev.txt", sep = "/"), header = TRUE)

neural_filt = total_filt[,which(colnames(total_filt) %in% rownames(umap))]

neural_sc_mod <-function(mod){
  mod_genes = colnames(net_hcCS_filt)[which(moduleColors == mod)]
  neural_sc_mod = neural_filt[which(rownames(neural_filt) %in% mod_genes),]
  neural_sc_mod_comp = neural_sc_mod[complete.cases(neural_sc_mod),]
  neural_sc_mod_1 = neural_sc_mod_comp[which(rowSums(neural_sc_mod_comp)>1),]
  mod.pca = prcomp(t(neural_sc_mod_1), scale=TRUE, center=TRUE)
  mod.pca_data = as.data.frame(mod.pca$x)
  eig = mod.pca_data$PC1
  
  # Re-orient eigengenes for consistency with WGCNA eigengenes
  if(cor(colMeans(neural_sc_mod_1), eig)[1] < 0 ){
    eig = (-1)*eig
  }
  
  coords = umap
  col_by = eig
  bplot = ggplot() +
    ggtitle(mod) +
    theme(plot.title = element_text(size=fontsize)) +
    geom_point(aes(x = coords[,1], y = coords[,2], color = col_by), size = pointsize) + 
    scale_color_gradient(low = "grey90", high = "maroon4", na.value = "grey80") +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    theme(axis.title.y =element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    theme(legend.text=element_text(size=fontsize - 3),legend.title=element_text(size=fontsize - 3)) +
    theme(legend.key.size = unit(0.1, "in")) +
    labs(color = "Exp")
  #theme(title = element_blank()) 
  bplot
}

# NOT USED
pdf(paste(plot_dir, "neural_sc_WGCNA_Modules.pdf", sep = "/"), height = 2, width =2, useDingbats = FALSE)
fontsize = 8
pointsize = 0.1
neural_sc_mod("blue")
neural_sc_mod("brown")
neural_sc_mod("red")
neural_sc_mod("magenta")
neural_sc_mod("salmon")
dev.off()

 # Salmon module sign test
hycs_deseq2 = read.table(paste(analysis_dir, "hyCS.DESEq2.ASE.highConf.noAneuploidy.noSex.rev.txt", sep = "/"), header = TRUE, row.names = 1, sep = "\t", fill = TRUE)
salmon_genes = colnames(net_hcCS_filt)[which(moduleColors == "salmon")]
hycs_deseq2_salmon = hycs_deseq2[which(rownames(hycs_deseq2) %in% salmon_genes ),]
dim(hycs_deseq2_salmon)

head(hycs_deseq2_salmon)
fontsize = 8

data = as.data.frame(cbind(hycs_deseq2_salmon$D50_HumanLFC, hycs_deseq2_salmon$D50_HumanLFC>0))
names(data) = c("LFC", "human")
p1 = ggplot(data = data, aes(x = as.numeric(as.character(LFC)), fill = as.factor(human))) + 
  geom_histogram(position = 'identity', color = "black", size = 0.3, binwidth = 0.1, boundary = 1) + 
  ggtitle("Day 50") + 
  xlab("Log2(Human/Chimpanzee)") + 
  ylab("Frequency")  + 
  scale_fill_manual(values = c("dodgerblue2", "salmon")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(axis.text = element_text(size = fontsize), axis.title = element_text(size = fontsize), title = element_text(size = fontsize))+
  theme(legend.position = "none")
p1

length(which(hycs_deseq2_salmon$D50_HumanLFC>0)) # 103
length(which(hycs_deseq2_salmon$D50_HumanLFC<0)) # 62
wilcox.test(hycs_deseq2_salmon$D50_HumanLFC, hycs_deseq2$D50_HumanLFC) # p-value = 6.049e-06


data = as.data.frame(cbind(hycs_deseq2_salmon$D100_HumanLFC, hycs_deseq2_salmon$D100_HumanLFC>0))
names(data) = c("LFC", "human")
p2 = ggplot(data = data, aes(x = as.numeric(as.character(LFC)), fill = as.factor(human))) + 
  geom_histogram(position = 'identity', color = "black", size = 0.3, binwidth = 0.1, boundary = 1) + 
  ggtitle("Day 100") + 
  xlab("Log2(Human/Chimpanzee)") + 
  ylab("Frequency")  + 
  scale_fill_manual(values = c("dodgerblue2", "salmon")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(axis.text = element_text(size = fontsize), axis.title = element_text(size = fontsize), title = element_text(size = fontsize))+
  theme(legend.position = "none")
p2

length(which(hycs_deseq2_salmon$D100_HumanLFC>0)) # 122
length(which(hycs_deseq2_salmon$D100_HumanLFC<0)) # 44
wilcox.test(hycs_deseq2_salmon$D100_HumanLFC, hycs_deseq2$D100_HumanLFC) #p-value = 2.566e-12


data = as.data.frame(cbind(hycs_deseq2_salmon$D150_HumanLFC, hycs_deseq2_salmon$D150_HumanLFC>0))
names(data) = c("LFC", "human")
p3 = ggplot(data = data, aes(x = as.numeric(as.character(LFC)), fill = as.factor(human))) + 
  geom_histogram(position = 'identity', color = "black", size = 0.3, binwidth = 0.1, boundary = 1) + 
  ggtitle("Day 150") + 
  xlab("Log2(Human/Chimpanzee)") + 
  ylab("Frequency")  + 
  scale_fill_manual(values = c("dodgerblue2", "salmon")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(axis.text = element_text(size = fontsize), axis.title = element_text(size = fontsize), title = element_text(size = fontsize))+
  theme(legend.position = "none")
p3

length(which(hycs_deseq2_salmon$D150_HumanLFC>0)) # 115
length(which(hycs_deseq2_salmon$D150_HumanLFC<0)) # 52
wilcox.test(hycs_deseq2_salmon$D150_HumanLFC, hycs_deseq2$D150_HumanLFC) #p-value = 9.457e-11


data = as.data.frame(cbind(hycs_deseq2_salmon$D200_HumanLFC, hycs_deseq2_salmon$D200_HumanLFC>0))
names(data) = c("LFC", "human")
p4 = ggplot(data = data, aes(x = as.numeric(as.character(LFC)), fill = as.factor(human))) + 
  geom_histogram(position = 'identity', color = "black", size = 0.3, binwidth = 0.1, boundary = 1) + 
  ggtitle("Day 200") + 
  xlab("Log2(Human/Chimpanzee)") + 
  ylab("Frequency")  + 
  scale_fill_manual(values = c("dodgerblue2", "salmon")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(axis.text = element_text(size = fontsize), axis.title = element_text(size = fontsize), title = element_text(size = fontsize)) +
  theme(legend.position = "none")
p4
length(which(hycs_deseq2_salmon$D200_HumanLFC>0)) # 122
length(which(hycs_deseq2_salmon$D200_HumanLFC<0)) # 39
wilcox.test(hycs_deseq2_salmon$D200_HumanLFC, hycs_deseq2$D200_HumanLFC) #p-value = 1.486e-14

# Figure 3i,j,k,l
pdf(paste(plot_dir, "salmon_wilcoxon.pdf", sep = "/"), height = 1.5, width =6, useDingbats = FALSE)
grid.arrange(p1,p2,p3,p4, ncol = 4)
dev.off()



for(i in salmon_genes){cat(i, "\n")}

salmon_deseq150 = as.data.frame(cbind(rownames(hycs_deseq2_salmon), hycs_deseq2_salmon$D150_HumanLFC))
names(salmon_deseq150) = c("gene", "LFC")
salmon_deseq150$gene = as.character(salmon_deseq150$gene)
salmon_deseq150$LFC = as.numeric(as.character(salmon_deseq150$LFC))
salmon_deseq150 = salmon_deseq150[order(salmon_deseq150$LFC, decreasing = TRUE),]

write.table(file = paste(analysis_dir, "hyCS.wgcnaSalmon.deseq150.rev.txt", sep = "/"), x = salmon_deseq150 , quote = FALSE, col.names = F, row.names = F, sep = "\t")



astro_overlap = c("TLR4", "PTPN11", "S100B", "EGFR", "CNTF", "SERPINE2", "PLP1", "SOX6")
hycs_deseq2_astro = hycs_deseq2[which(rownames(hycs_deseq2) %in% astro_overlap ),]
dim(hycs_deseq2_astro)

# None significant
head(hycs_deseq2_astro)
hist(hycs_deseq2_astro$D200_HumanLFC)
length(which(hycs_deseq2_astro$D200_HumanLFC>0))
length(which(hycs_deseq2_astro$D200_HumanLFC<0)) 
wilcox.test(hycs_deseq2_astro$D200_HumanLFC, hycs_deseq2$D200_HumanLFC) 


glia_overlap = c("LYN", "CLU", "TLR4", "PTPN11", "S100B", "EGFR", "CNTF", "SERPINE2", "PLP1", "SOX6", "OLIG1")
hycs_deseq2_glia = hycs_deseq2[which(rownames(hycs_deseq2) %in% glia_overlap ),]
dim(hycs_deseq2_glia)

head(hycs_deseq2_glia)
hist(hycs_deseq2_glia$D50_HumanLFC)
length(which(hycs_deseq2_glia$D50_HumanLFC>0)) # 5
length(which(hycs_deseq2_glia$D50_HumanLFC<0)) # 3
wilcox.test(hycs_deseq2_glia$D50_HumanLFC, hycs_deseq2$D50_HumanLFC) # NS

hist(hycs_deseq2_glia$D100_HumanLFC)
length(which(hycs_deseq2_glia$D100_HumanLFC>0)) # 5
length(which(hycs_deseq2_glia$D100_HumanLFC<0)) # 3
wilcox.test(hycs_deseq2_glia$D100_HumanLFC, hycs_deseq2$D100_HumanLFC) # NS

hist(hycs_deseq2_glia$D150_HumanLFC)
length(which(hycs_deseq2_glia$D150_HumanLFC>0)) # 7
length(which(hycs_deseq2_glia$D150_HumanLFC<0)) # 1
wilcox.test(hycs_deseq2_glia$D150_HumanLFC, hycs_deseq2$D150_HumanLFC) # p-value = 0.06675

hist(hycs_deseq2_glia$D200_HumanLFC)
length(which(hycs_deseq2_glia$D200_HumanLFC>0)) # 7
length(which(hycs_deseq2_glia$D200_HumanLFC<0)) # 1
wilcox.test(hycs_deseq2_glia$D200_HumanLFC, hycs_deseq2$D200_HumanLFC) # p-value = 0.02675



# Looking at the Salmon module in monocle data
total_filt = read.table(file=paste(analysis_dir, "total_filt.rev.txt", sep = "/"), header = TRUE)
human_filt = read.table(file=paste(analysis_dir, "human_filt.rev.txt", sep = "/"), header = TRUE)
chimp_filt = read.table(file=paste(analysis_dir, "chimp_filt.rev.txt", sep = "/"), header = TRUE)

human_filt_cpm = sweep(as.matrix(human_filt), 2, as.double(colSums(human_filt)/1000000), `/`)
chimp_filt_cpm = sweep(as.matrix(chimp_filt), 2, as.double(colSums(chimp_filt)/1000000), `/`)


monocle_module <- function(gene){
  mod_genes = colnames(net_hcCS_filt)[which(moduleColors == mod)]
  h_points = human_filt_cpm[which(rownames(human_filt_cpm) %in% mod_genes),which(colnames(total_filt) %in% rownames(pData(rcds)))]
  c_points = chimp_filt_cpm[which(rownames(chimp_filt_cpm) %in% mod_genes),which(colnames(total_filt) %in% rownames(pData(rcds)))]
  pseudotime = pData(rcds)$Pseudotime
  hc_filt_cpm = cbind(h_points, c_points)
  
  mono_mod_comp = hc_filt_cpm[complete.cases(hc_filt_cpm),]
  mono_mod_1 = mono_mod_comp[which(rowSums(mono_mod_comp)>1),]
  mod.pca = prcomp(t(mono_mod_1), scale=TRUE, center=TRUE)
  mod.pca_data = as.data.frame(mod.pca$x)
  eig = mod.pca_data$PC1
  
  # Re-orient eigengenes for consistency with WGCNA eigengenes
  if(cor(colMeans(mono_mod_1), eig)[1] < 0 ){
    eig = (-1)*eig
  }
  
  
  all_points = as.data.frame(cbind(rep(pseudotime, 2), eig, c(rep("Human", 348), rep("Chimp", 348))))
  names(all_points) = c("Pseudotime", "Eig", "Species")
  all_points$Pseudotime = as.numeric(as.character(all_points$Pseudotime))
  all_points$Eig = as.numeric(as.character(all_points$Eig))
  
  ggplot(data = all_points, aes(x = Pseudotime, y = Eig, col = as.factor(Species))) +
    geom_point(size = pointsize) +
    scale_color_manual(values = c("dodgerblue2", "salmon")) +
    geom_smooth(data = as.data.frame(all_points[which(all_points$Species == "Human"),]), aes(x = Pseudotime, y = Eig ), size = 0.25,linetype = "dashed", col = "indianred3", se = FALSE) +
    geom_smooth(data = as.data.frame(all_points[which(all_points$Species == "Chimp"),]), aes(x = Pseudotime, y = Eig ), size = 0.25,linetype = "dashed", col = "dodgerblue4", se = FALSE) +
    theme_bw() +
    ylim(-2, 20) +
    theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(legend.position = "none")
  
}

# Figure 3h
pdf(paste(plot_dir, "monocle.allelic.salmon.rev.pdf", sep = "/"), height = 2.5, width = 2.5, useDingbats = FALSE)
pointsize = 0.03
monocle_module("salmon")
dev.off()



