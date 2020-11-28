## Analysis and plots for Figures 2c-f, 4b and Extended Data Figures 4a-c, 5, 8e: hybrid cortical spheroids
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


# Induction success rates (how many make it to D100)
ff_suc = 1/12
disp_suc = 0/9
disp_mat_suc = 9/9

rates = c(ff_suc, disp_suc, disp_mat_suc)
ratedf = data.frame(rate = rates, lab = c("Feeder Free", "Feeder", "Feeder + Matrigel"))
ratedf$lab = factor(ratedf$lab, levels = c("Feeder Free", "Feeder", "Feeder + Matrigel"), ordered = TRUE)
ratedf = ratedf[order(ratedf$lab),]

# Extended Data Figure 4c
pdf(paste(plot_dir, "successRates.rev.pdf", sep = "/"), height = 2, width = 1.5, useDingbats = FALSE)
fontsize = 8
ggplot(data = ratedf, aes(x = as.factor(lab), y = rate)) +
  geom_bar(stat = "identity", color = "black") +
  ylab("Percent Success") +
  
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  #scale_color_manual(values = c("dodgerblue2","salmon", "mediumpurple2")) +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize)) + 
  theme(axis.title.x = element_blank(), axis.text.x=element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none")
dev.off()


## Load data 
# Gene meta data
geneInfo = read.table(paste(data_dir, "geneInfo.GRCh38.txt", sep = "/"))
names(geneInfo) = c("gene", "length", "chrom")
geneInfo$genesizeKb = as.vector(geneInfo$length/1000)
dim(geneInfo)
head(geneInfo)


# Sngle cell RNAseq data
sc = read.table(paste(data_dir, "scRNAseq.GRCh38.txt", sep = "/"), header = TRUE)
rownames(sc) = geneInfo$gene
dim(sc)
sc_meta = read.table(paste(data_dir, "scRNAseq.meta.txt", sep = "/"), header = TRUE)
sc_meta$Line = as.character(sc_meta$Line)



# Calculate TPM and CPM
geneKeepInfo = read.table(file = paste(analysis_dir, "geneKeep.rev.txt", sep = "/"), header = TRUE, fill = TRUE)
geneKeep = as.vector(geneKeepInfo$geneKeep)
geneKeepIdx = as.vector(geneKeepInfo$geneKeepIdx)
sc_keep = sc[geneKeep,]
sc_cpm = sweep(as.matrix(sc_keep), 2, as.double(colSums(sc_keep)/1000000), `/`)
sc_fpkm = sweep(as.matrix(sc_cpm), 1, as.double(geneInfo$genesizeKb[geneKeepIdx]), `/`)
sc_rpk = sweep(as.matrix(sc_keep), 1, as.double(geneInfo$genesizeKb[geneKeepIdx]), `/`)
sc_tpm = sweep(as.matrix(sc_rpk), 2, as.double(colSums(sc_rpk)/1000000), `/`)

human = sc_keep[, seq(2,length(sc),by=3)]
chimp = sc_keep[, seq(3,length(sc),by=3)]
total = sc_keep[, seq(1,length(sc),by=3)]

dim(human)
dim(chimp)
dim(total)

# Aggregated total counts 
aggregates = cbind(rowSums(human), rowSums(chimp), rowSums(total))
colnames(aggregates) = c("Human", "Chimp", "Total")

data = CreateSeuratObject(counts = total, project = "hyCS", min.cells = 3, min.features = 200)
data

mitoGenes = geneInfo$gene[which(geneInfo$chrom=="chrMT")]
mitoGenes = mitoGenes[which(mitoGenes %in% rownames(data))]
data[["percent.mt"]] <- PercentageFeatureSet(data, features = mitoGenes)

VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

data <- subset(data, subset = nCount_RNA > 100000 & percent.mt < 7)

data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 1500)
top10 <- head(VariableFeatures(data), 10)
plot1 <- VariableFeaturePlot(data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
data <- RunPCA(data, features = VariableFeatures(object = data))
print(data[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(data, dims = 1:2, reduction = "pca")
DimPlot(data, reduction = "pca")
DimHeatmap(data, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(data)


data <- FindNeighbors(data, dims = 1:13)
data <- FindClusters(data, resolution = 0.4)

data <- RunUMAP(data, dims = 1:14)
DimPlot(data, reduction = "umap")

data <- RunTSNE(data, dims = 1:14)
DimPlot(data, reduction = "tsne")

save(data, file = paste(analysis_dir, "seurat_sc.Rdata", sep = "/"))
write.table(file = paste(analysis_dir, "sc_umap_coords.rev.txt", sep = "/"), x = data@reductions$umap@cell.embeddings, col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(file = paste(analysis_dir, "sc_tsne_coords.rev.txt", sep = "/"), x = data@reductions$tsne@cell.embeddings, col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(file = paste(analysis_dir, "sc_meta_data_clusters.rev.txt", sep = "/"), x = data@meta.data, col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

#*coords = read.table(file = paste(analysis_dir, "sc_umap_coords.rev.txt", sep = "/"), header = TRUE, row.names = 1)
#*seurat_clusters = read.table(file = paste(analysis_dir, "sc_meta_data_clusters.rev.txt", sep = "/"), row.names = 1, header = TRUE)


# Plot cells colored by cell type
# Extended Data Figure 5a
pdf(paste(plot_dir, "sc_umap.rev.pdf", sep = "/"), height = 3, 3, useDingbats = FALSE) 
fontsize = 8
pointsize = 0.35
coords = data@reductions$umap@cell.embeddings
bplot = ggplot() +
  geom_point(aes(x = coords[,1], y = coords[,2], color = as.factor(data@meta.data$seurat_clusters)), size = pointsize) + 
  scale_color_manual(values = c("goldenrod1", "royalblue4", "dodgerblue1", "red4", "burlywood1", "grey60", "darkgoldenrod")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position = "none")
bplot
dev.off()


# Top cluster genes
data.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
data.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top10 <- data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top100 <- data.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)

write.table(file = paste(analysis_dir, "top100ClusterGenes.all.rev.txt", sep = "/"),
            data.frame(cbind(top100$cluster, top100$gene, top100$avg_logFC, top100$p_val_adj)),
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)


#pdf(paste(plot_dir, "sc_clusterHeatmap.pdf", sep = "/"), height = 8, 8) 
#fontsize = 1
#DoHeatmap(data, features = top10$gene, size = 3, 
#          group.colors = c("goldenrod1", "royalblue4", "dodgerblue1", "red4", "burlywood1", "grey60", "darkgoldenrod"),
#          raster = FALSE) + NoLegend()
#dev.off()

# Dot Plot for marker genes
# Extended Data Figure 5c
pdf(paste(plot_dir, "sc_dotPlot.rev.pdf", sep = "/"), height = 2.5, 4, useDingbats = FALSE)
fontsize = 8
markers = as.vector(c("COL1A2", "COL3A1", "NEUROD6", "SLA", "DLX6", "DLX5", "SLC1A3", "ADGRV1", "KRT5", "KRT19", "CDK1", "TOP2A", "COL9A1", "COL2A1"))
DotPlot(data, features = markers, dot.min = 0.1, col.min = 0, col.max = 0, dot.scale = 3) +
  ylab("Cluster") +
  #geom_point(col = rep(c("goldenrod1", "royalblue4", "dodgerblue1", "red4", "burlywood1", "grey60", "darkgoldenrod"), each = 14), size = ) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 65, hjust = 1)) +
  theme(axis.title.y = element_text(size = fontsize)) +
  theme(axis.text = element_text(size = fontsize))
dev.off()


# Filter original data frames
human_filt = human[rownames(data@assays$RNA),which(colnames(total) %in% colnames(data@assays$RNA))]
chimp_filt = chimp[rownames(data@assays$RNA),which(colnames(total) %in% colnames(data@assays$RNA))]
total_filt = total[rownames(data@assays$RNA),colnames(data@assays$RNA)]
aggregates_filt = as.data.frame(cbind(rowSums(human_filt), rowSums(chimp_filt), rowSums(total_filt)))
colnames(aggregates_filt) = c("Human", "Chimp", "Total")

dim(human_filt)
dim(chimp_filt)
dim(total_filt)
dim(aggregates_filt)


write.table(x = total_filt, file=paste(analysis_dir, "total_filt.rev.txt", sep = "/"), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t" )
write.table(x = human_filt, file=paste(analysis_dir, "human_filt.rev.txt", sep = "/"), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t" )
write.table(x = chimp_filt, file=paste(analysis_dir, "chimp_filt.rev.txt", sep = "/"), col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t" )

# Meta data
meta_filt = sc_meta[which(colnames(total) %in% colnames(total_filt)),]
dim(meta_filt)

# Line
line = meta_filt$Line
hl2_16_sc = total_filt[, which(meta_filt$Line=="HL2-16")]
hl2_9_sc = total_filt[, which(meta_filt$Line=="HL2-9")]

dim(hl2_16_sc)
dim(hl2_9_sc)

# NOT USED
pdf(paste(plot_dir, "sc_umap_Line.rev.pdf", sep = "/"), height = 1.5,1.5, useDingbats = FALSE) 
pointsize = 0.1
col = (line == "HL2-9")
colbg = (line == "HL2-16")
bplot = ggplot() +
  geom_point(aes(x = coords[,1], y = coords[,2]),size = pointsize) + 
  xlab("Dimension 1") + 
  ylab("Dimension 2")  + 
  geom_point(aes(coords[,1][colbg], coords[,2][colbg]), colour = "slateblue3",size = pointsize) + 
  geom_point(aes(coords[,1][col], coords[,2][col]), colour = "mediumseagreen",size = pointsize) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position = "none")
print(bplot)
dev.off()

# Bar plot of cells per line per cluster
line = sc_meta$Line[which(colnames(total) %in% colnames(total_filt))]
clustLine = as.data.frame(cbind(line, data$seurat_clusters))#as.vector(seurat_clusters$data.seurat_clusters)))
names(clustLine) = c("Line", "Cluster")
clustLine.summary = group_by(clustLine, Cluster, Line)
library(plyr)
clustLine.summary = count(clustLine.summary)

length(which(line=="HL2-16")) # 338
length(which(line=="HL2-9")) # 368


clustLine.summary$percent = clustLine.summary$freq / rep(c(338, 368), 7)

# Extended Data Figure 5b
pdf(paste(plot_dir, "Cluster_by_Line.rev.pdf", sep = "/"), height = 2, 2) 
fontsize = 8
p = ggplot(data = as.data.frame(clustLine.summary), aes(x = as.factor(Cluster), y = percent, fill = Line)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("Percentage") +
  xlab("Cluster") +
  scale_fill_manual(values=c("slateblue3", "mediumseagreen"))+
  ggtitle("") +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))+ theme(legend.title = element_blank()) 

p
dev.off()


fontsize = 8
pointsize =0.5
genePlot <- function(gene,exp = total_filt){
  #coords = data@reductions$umap@cell.embeddings
  coords = read.table(file = paste(analysis_dir, "sc_umap_coords.rev.txt", sep = "/"), header = TRUE, row.names = 1)
  col_by = log10(as.numeric(exp[gene,]))
  bplot = ggplot() +
    ggtitle(gene) +
    theme(plot.title = element_text(size=fontsize)) +
    geom_point(aes(x = coords[,1], y = coords[,2], color = col_by), size = pointsize) + 
    scale_color_gradient(low = "grey80", high = "maroon4", na.value = "grey80") +
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


# Extended Data Figure 5d
pdf(paste(plot_dir, "sc_genePlots.pdf", sep = "/"), height = 16, 8, useDingbats = FALSE) 
fontsize = 8
pointsize = 0.3
p1 = genePlot("COL1A1")
p2 = genePlot("COL3A1")
p3 = genePlot("POSTN")
p4 = genePlot("LUM")
p5 = genePlot("FMOD")
p6 = genePlot("COL9A1")
p7 = genePlot("KRT19")
p8 = genePlot("KRT14")
p9 = genePlot("KRT5")
p10 = genePlot("PERP")
p11 = genePlot("NEUROD6")
p12 = genePlot("STMN2")
p13 = genePlot("SLA")
p14 = genePlot("MAP2")
p15 = genePlot("DLX5")
p16 = genePlot("DLX6")
p17 = genePlot("DLX1")
p18 = genePlot("DLX2")
p19 = genePlot("SATB2")
p20 = genePlot("TBR1")
p21 = genePlot("POU3F2")
p22 = genePlot("EOMES")
p23 = genePlot("CDK1")
p24 = genePlot("TOP2A")
p25 = genePlot("MKI67")
p26 = genePlot("VIM")
p27 = genePlot("PAX6")
p28 = genePlot("FOXG1")
p29 = genePlot("SOX9")
p30 = genePlot("HOPX")
p31 = genePlot("AGT")
p32 = genePlot("FAM107A")
p33 = genePlot("PON2")
p34 = genePlot("HEPACAM")
p35 = genePlot("GFAP")

grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,
             p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,
             p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,
             p31,p32,p33,p34,p35, ncol = 4)
dev.off()


#plot = DimPlot(data, reduction = "umap")
#select.cells <- CellSelector(plot = plot)

#write.table(as.data.frame(n_data$seurat_clusters), paste(analysis_dir, "seurat.neural.cells.txt", sep = "/"), quote = FALSE, col.names = TRUE, row.names = TRUE)
select.cells = read.table(paste(analysis_dir, "seurat.neural.cells.txt", sep = "/"), header = TRUE, row.names = 1)
select.cells = as.vector(rownames(select.cells))


# Make a plot showing which cells were selected as neural
# Extended Data Figure 5h
pdf(paste(plot_dir, "sc_umapSelected.rev.pdf", sep = "/"), height = 2.5, 2.5, useDingbats = FALSE) 
fontsize = 8
pointsize = 0.3
bplot = ggplot() +
  theme(plot.title = element_text(size=fontsize)) +
  geom_point(aes(x = coords[,1], y = coords[,2], color = colnames(total_filt) %in% select.cells), size = pointsize) + 
  scale_color_manual(values = c("grey80", "red2")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position = "none")
bplot
dev.off()

# Define each cell type
MES_CELLS = names(data$orig.ident[which(data$seurat_clusters==0 | data$seurat_clusters==6)])
EPI_CELLS = names(data$orig.ident[which(data$seurat_clusters==4)])
ASTRO_CELLS = names(data$orig.ident[which(data$seurat_clusters==3)])
NEURON_CELLS = names(data$orig.ident[which(data$seurat_clusters==1 | data$seurat_clusters==2)])
PROGENITOR_CELLS = names(data$orig.ident[which(data$seurat_clusters==5)])


# Percentage of expression per cell type
epi_percent = round((rowSums(total[, EPI_CELLS]) / aggregates[,3]) * 100, 1)
mes_percent = round((rowSums(total[, MES_CELLS]) / aggregates[,3]) * 100, 1)
neuron_percent = round((rowSums(total[, NEURON_CELLS]) / aggregates[,3]) * 100, 1)
progenitor_percent = round((rowSums(total[, PROGENITOR_CELLS]) / aggregates[,3]) * 100, 1)
astro_percent = round((rowSums(total[, ASTRO_CELLS]) / aggregates[,3]) * 100, 1)

neuron_percent[is.nan(neuron_percent)] = 0
progenitor_percent[is.nan(progenitor_percent)] = 0
astro_percent[is.nan(astro_percent)] = 0
mes_percent[is.nan(mes_percent)] = 0
epi_percent[is.nan(epi_percent)] = 0

# ASE by cell type
epi_ase = log2(rowSums(human[, which(colnames(total) %in% EPI_CELLS)]) / rowSums(chimp[, which(colnames(total) %in% EPI_CELLS)]))
mes_ase = log2(rowSums(human[, which(colnames(total) %in% MES_CELLS)]) / rowSums(chimp[, which(colnames(total) %in% MES_CELLS)]))
neuron_ase = log2(rowSums(human[, which(colnames(total) %in% NEURON_CELLS)]) / rowSums(chimp[, which(colnames(total) %in% NEURON_CELLS)]))
progenitor_ase = log2(rowSums(human[, which(colnames(total) %in% PROGENITOR_CELLS)]) / rowSums(chimp[, which(colnames(total) %in% PROGENITOR_CELLS)]))
astro_ase = log2(rowSums(human[, which(colnames(total) %in% ASTRO_CELLS)]) / rowSums(chimp[, which(colnames(total) %in% ASTRO_CELLS)]))

cellTypeStats = as.data.frame(cbind(progenitor_percent,astro_percent,mes_percent,epi_percent,neuron_percent,progenitor_ase,astro_ase,mes_ase,epi_ase,neuron_ase))
names(cellTypeStats) = c("progenitor_percent","astro_percent","mes_percent","epi_percent","neuron_percent","progenitor_ase","astro_ase","mes_ase","epi_ase","neuron_ase")
write.table(x = cellTypeStats, file = paste(analysis_dir, "scCellTypeStats.rev.txt", sep = "/"), sep="\t", row.names = TRUE,
            col.names = TRUE, quote=FALSE)


# ASE histograms
binsize = 1.8
p1 <- qplot(astro_ase, 
            geom = "histogram", 
            main="ASE: astroglia", 
            xlab="Log2(ASE Ratio)",
            ylab="Frequency",
            fill=I("red4"), 
            col=I("black"), 
            binwidth=binsize) +
  theme_bw() + 
  theme(axis.text = element_text(size = fontsize), axis.title = element_text(size = fontsize), title = element_text(size = fontsize)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept = 0) 

p2 <-qplot(neuron_ase, 
           geom = "histogram", 
           main="ASE: neurons", 
           xlab="Log2(ASE Ratio)",
           ylab="Frequency",
           fill=I("dodgerblue1"), 
           col=I("black"), 
           binwidth=binsize) +
  theme_bw() + 
  theme(axis.text = element_text(size = fontsize), axis.title = element_text(size = fontsize), title = element_text(size = fontsize)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept = 0)

p3 <-qplot(progenitor_ase, 
           geom = "histogram", 
           main="ASE: cycling progenitors", 
           xlab="Log2(ASE Ratio)",
           ylab="Frequency",
           fill=I("gray60"), 
           col=I("black"), 
           binwidth=binsize) +
  theme_bw() + 
  theme(axis.text = element_text(size = fontsize), axis.title = element_text(size = fontsize), title = element_text(size = fontsize)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept = 0)

p4 <-qplot(mes_ase, 
           geom = "histogram", 
           main="ASE: mesenchyme", 
           xlab="Log2(ASE Ratio)",
           ylab="Frequency",
           fill=I("goldenrod1"), 
           col=I("black"), 
           binwidth=binsize) +
  theme_bw() + 
  theme(axis.text = element_text(size = fontsize), axis.title = element_text(size = fontsize), title = element_text(size = fontsize)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept = 0)

p5 <-qplot(epi_ase, 
           geom = "histogram", 
           main="ASE: epithelial", 
           xlab="Log2(ASE Ratio)",
           ylab="Frequency",
           fill=I("burlywood1"), 
           col=I("black"),
           binwidth=binsize) +
  theme_bw() + 
  theme(axis.text = element_text(size = fontsize), axis.title = element_text(size = fontsize), title = element_text(size = fontsize)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept = 0)

# Extended Data Figure 5i
pdf(paste(plot_dir, "scASE_hist.pdf", sep = "/"), height = 1.5, width = 6) 
grid.arrange(p1,p2,p3,p4,p5, nrow = 1)
dev.off()


# Define genes that are neuralectoderm-specific
neural_genes_filter = which(progenitor_percent+astro_percent+neuron_percent >=90)

neural_genes = rownames(total)[neural_genes_filter]
length(neural_genes) # 6222

write.table(x = as.data.frame(neural_genes), file = paste(analysis_dir, "neuralGenes.rev.txt", sep = "/"), sep="\t", row.names = FALSE,
            col.names = TRUE, quote=FALSE)


## Testing differences between matrigel versus no matrigel
# SY15 data 
sy = read.table(paste(data_dir, "SY15_CS.GRCh38.txt", sep = "/"), header = TRUE)
rownames(sy) = geneInfo$gene
dim(sy)
sy_meta = read.table(paste(data_dir, "SY15_CS.meta.txt", sep = "/"), header = TRUE)
sy_meta$Line = as.character(sy_meta$Line)
sy_meta$Species = as.character(sy_meta$Species)

# Calculate TPM and CPM
sy_tot = sy[, seq(1,length(sy),by=3)]
sy_tot_keep = sy_tot[geneKeep,]
sy_tot_cpm = sweep(as.matrix(sy_tot_keep), 2, as.double(colSums(sy_tot_keep)/1000000), `/`)

ff = sy_tot_cpm[,11]

# Bulk hyCS data
hyCS = read.table(paste(data_dir, "hyCS.GRCh38.txt", sep = "/"), header = TRUE)
rownames(hyCS) = geneInfo$gene
dim(hyCS)
hyCS_tot = hyCS[, seq(1,75,by=3)]
hyCS_tot_keep = hyCS_tot[geneKeep,]
hyCS_tot_cpm = sweep(as.matrix(hyCS_tot_keep), 2, as.double(colSums(hyCS_tot_keep)/1000000), `/`)

matrigel = hyCS_tot_cpm[,14]

ff[which(!is.finite(ff))] <-NA
matrigel[which(!is.finite(matrigel))] <-NA

tempdf = as.data.frame(cbind(ff[which(matrigel>(2.8*ff + 10))], matrigel[which(matrigel>(2.8*ff + 10))]))

index = which(rownames(sy_tot_cpm) %in% c("FMOD")) #c("ACAN", "COL2A1", "COL11A1", "COL9A1", "COL1A1", "FMOD"))

# Extended Data Figure 5e
pdf(paste(plot_dir, "matrigelCompare.rev.pdf", sep = "/"), height = 2.5, 2.5, useDingbats = FALSE) 
fontsize = 8
pointsize = 0.4
p = ggplot() + 
  geom_point(aes(x = log10(ff), y = log10(matrigel)), size = pointsize) +  
  ggtitle("Effects of matrigel") +
  xlab("SY15 D50 HL1-29 (Log10 CPM)") +
  ylab("RA8 D50 HL1-29 (Log10 CPM)") + 
  theme(legend.position="none")  + 
  geom_point(color='black') +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_abline(slope=1, intercept = 0) +
  geom_point(aes(log10(tempdf$V1), log10(tempdf$V2)), color="firebrick3", size = pointsize)+
  geom_point(aes(log10(ff[index]), log10(matrigel[index])), color="darkgreen", size = pointsize) +
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))
p
dev.off()

cor(ff, matrigel)*cor(ff, matrigel) # 0.120

# Worst offenders
which(log10(matrigel) > 3.5 & log10(ff) <2)
which(log10(matrigel) > 2 & log10(ff) <0)

# NOT USED
pdf(paste(plot_dir, "matrigelComparePlots.rev.pdf", sep = "/"), height = 4, 6, useDingbats = FALSE) 
p1 = genePlot("COL2A1")
p2 = genePlot("COL11A1")
p3 = genePlot("COL1A1")
p4 = genePlot("COL9A1")
p5 = genePlot("ACAN")
p6 = genePlot("FMOD")
grid.arrange(p1,p2,p3,p4,p5,p6, ncol = 3)
dev.off()

# PCA of SY15 only
iPS = read.table(paste(data_dir, "iPSC.GRCh38.txt", sep = "/"), header = TRUE)
rownames(iPS) = geneInfo$gene
dim(iPS)
iPS_meta = read.table(paste(data_dir, "iPSC.meta.txt", sep = "/"), header = TRUE)
iPS_meta$Species = as.character(iPS_meta$Species)
iPS_meta$Line = as.character(iPS_meta$Line)
iPS_tot = iPS[, seq(1,66,by=3)]
iPS_tot_keep = iPS_tot[geneKeep,]
iPS_tot_cpm = sweep(as.matrix(iPS_tot_keep), 2, as.double(colSums(iPS_tot_keep)/1000000), `/`)
sy_iPS_cpm = cbind(iPS_tot_cpm, sy_tot_cpm)
sy_tot_cpm_comp = sy_iPS_cpm[complete.cases(sy_iPS_cpm),]

sy_tot_cpm_1 = sy_tot_cpm_comp[which(rowSums(sy_tot_cpm_comp)>1),]
dim(sy_tot_cpm_1)

sy.pca = prcomp(t(sy_tot_cpm_1), scale=TRUE, center=TRUE)
summary(sy.pca)
sy.pca_data = as.data.frame(sy.pca$x)

# Extended Data Figure 4a
pdf(paste(plot_dir, "sy_PCA.rev.pdf", sep = "/"), height = 2, 2, useDingbats = FALSE) 
fontsize = 8
p1 = ggplot(sy.pca_data, aes(PC1, PC2)) + 
  geom_point(aes(col=c(iPS_meta$Species, sy_meta$Species), shape = c(rep("iPS", 22), rep("CS", 11))), size = 1)  + 
  ggtitle("PCA") +
  xlab("PC1 (37% of variance)") +
  ylab("PC2 (11% of variance)") +
  theme_bw() +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values = c("dodgerblue2","salmon", "mediumpurple2")) +
  #labs(color = "Species")
  theme(legend.position = "none")
p1 
dev.off()

# With hybrids split
sy_hum_h = sy[, seq(2,length(sy),by=3)]
sy_hum_h = sy_hum_h[,c(2,5,8)]
dim(sy_hum_h)
head(sy_hum_h)

sy_chi_c = sy[, seq(3,length(sy),by=3)]
sy_chi_c = sy_chi_c[,c(1,4,7,10)]
dim(sy_chi_c)
head(sy_chi_c)

sy_c = sy[, seq(3,length(sy),by=3)]
sy_c = sy_c[, c(3,6,9,11)]
dim(sy_c)
head(sy_c)

sy_h = sy[, seq(2,length(sy),by=3)]
sy_h = sy_h[, c(3,6,9,11)]
dim(sy_h)
head(sy_h)

sy_hy_par = cbind(sy_c, sy_h, sy_hum_h, sy_chi_c)
head(sy_hy_par)

sy_hy_par = sweep(as.matrix(sy_hy_par), 2, as.vector(colSums(sy_hy_par)/1000000), `/`)
head(sy_hy_par)

sy_hy_par_noAnSex = sy_hy_par[geneKeep,]
dim(sy_hy_par_noAnSex)
head(sy_hy_par_noAnSex)

# iPS split species
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
sy_hy_par_noAnSex = cbind(sy_hy_par_noAnSex, iPS_hy_par_noAnSex)

# Expression threshold
sy_hy_par_comp = sy_hy_par_noAnSex[complete.cases(sy_hy_par_noAnSex),]

sy_hy_par_1 = sy_hy_par_comp[which(rowSums(sy_hy_par_comp)>1),]
dim(sy_hy_par_1)
head(sy_hy_par_1)

sy.hy_par.pca = prcomp(t(sy_hy_par_1), scale=TRUE, center=TRUE)
summary(sy.hy_par.pca)
sy.hy_par.pca_data = as.data.frame(sy.hy_par.pca$x)
sy_hy_par_species = c(rep("Hy-Chimp", 4), rep("Hy-Human", 4), rep("Human", 3), rep("Chimp", 4), rep("Hy-Chimp", 10), rep("Hy-Human", 10), rep("Human", 6), rep("Chimp", 6))


# Extended Data Figure 4b
pdf(paste(plot_dir, "sy_PCA_hyPar.rev.pdf", sep = "/"), height = 2, 2, useDingbats = FALSE) 
fontsize = 8
p = ggplot(sy.hy_par.pca_data, aes(PC1, PC2)) + 
  geom_point(aes(col=sy_hy_par_species, shape=c(rep("CS", 15), rep("iPS", 32))), size = 1)  + 
  ggtitle("PCA") +
  xlab("PC1 (35% of variance)") +
  ylab("PC2 (11% of variance)") +
  theme_bw() +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values = c("blue4","firebrick4", "deepskyblue", "indianred1")) +
  #labs(color = "Species") +
  theme(legend.position = "none")
p 
dev.off()


## Comparison with single cell data from Sloan et al.
## Publication: https://www.cell.com/neuron/fulltext/S0896-6273(17)30683-9
steven = read.table(paste(analysis_dir, "Steven_SmartSeq2.txt", sep = "/"), header = TRUE, row.names = 1)
dim(steven)

# List of genes that are in both data sets 
steven_genes = read.table(paste(analysis_dir, "steven_genes_filter.txt", sep = "/"), header=TRUE, row.names = 1)
steven_filt = steven[which(row.names(steven) %in% row.names(steven_genes)), colSums(steven)>100000]
dim(steven_filt)

# Steven metadata
steven_meta = read.csv(paste(analysis_dir, "Metadata_minibrain.csv", sep = "/"), header = TRUE, row.names = 1)
dim(steven_meta)

steven_meta_filt = steven_meta[colSums(steven)>100000,]
dim(steven_meta_filt)
head(steven_meta_filt)

all_steven_filt = total[which(row.names(total) %in% row.names(steven_genes)),]
dim(steven_filt)
dim(all_steven_filt)

steven_merged_all = as.data.frame(cbind(steven_filt[which(rownames(steven_filt) %in% rownames(all_steven_filt)),], all_steven_filt))
dim(steven_merged_all)

# Clustering both data sets
sr_data = CreateSeuratObject(counts = steven_merged_all, project = "sr", min.cells = 3, min.features = 200)
sr_data

mitoGenes = geneInfo$gene[which(geneInfo$chrom=="chrMT")]
mitoGenes = mitoGenes[which(mitoGenes %in% rownames(sr_data))]
sr_data[["percent.mt"]] <- PercentageFeatureSet(sr_data, features = mitoGenes)

VlnPlot(sr_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

sr_data <- subset(sr_data, subset = nCount_RNA > 100000 & percent.mt < 5)
sr_data <- NormalizeData(sr_data, normalization.method = "LogNormalize", scale.factor = 10000)
sr_data <- FindVariableFeatures(sr_data, selection.method = "vst", nfeatures = 1500)

top10 <- head(VariableFeatures(sr_data), 10)
plot1 <- VariableFeaturePlot(sr_data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(sr_data)
sr_data <- ScaleData(sr_data, features = all.genes)
sr_data <- RunPCA(sr_data, features = VariableFeatures(object = sr_data))
print(sr_data[["pca"]], dims = 1:5, nfeatures = 5)
ElbowPlot(sr_data)


sr_data <- FindNeighbors(sr_data, dims = 1:15)
sr_data <- FindClusters(sr_data, resolution = 0.4)

sr_data <- RunUMAP(sr_data, dims = 1:14)
DimPlot(sr_data, reduction = "umap")

sr_data <- RunTSNE(sr_data, dims = 1:14)
DimPlot(sr_data, reduction = "tsne")

save(sr_data, file = paste(analysis_dir, "steven_seurat_sc.Rdata", sep = "/"))

sr_data.markers <- FindAllMarkers(sr_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- sr_data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

#pdf(paste(plot_dir, "steven_sc_clusterHeatmap.pdf", sep = "/"), height = 10, 8) 
#fontsize = 1
#DoHeatmap(sr_data, features = top10$gene, size = 3, 
#          group.colors = c("royalblue4", "red4", "goldenrod1", "dodgerblue2", "grey60", "burlywood1", "violet", "red1", "darkgoldenrod"),
#          raster = FALSE) + NoLegend()
#dev.off()

write.table(file = paste(analysis_dir, "steven_sc_umap_coords.rev.txt", sep = "/"), x = sr_data@reductions$umap@cell.embeddings, col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(file = paste(analysis_dir, "steven_sc_tsne_coords.rev.txt", sep = "/"), x = sr_data@reductions$tsne@cell.embeddings, col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

sr_tsne = read.table(paste(analysis_dir, "steven_sc_tsne_coords.rev.txt", sep = "/"), header = T, row.names = 1)

firstRachel = which(names(sr_data$orig.ident) %in% colnames(total))[1]
source_all = c(rep("Steven", firstRachel-1), rep("Rachel", length(sr_data$orig.ident) - (firstRachel-1)))


# Extended Data Figure 5f
pdf(paste(plot_dir, "steven_sc_umap.rev.pdf", sep = "/"), height = 2.5, 2.5, useDingbats = FALSE) 
fontsize = 8
pointsize = 0.5
sr_coords = sr_tsne
bplot = ggplot() +
  theme(plot.title = element_text(size=fontsize)) +
  geom_point(aes(x = sr_coords[,1], y = sr_coords[,2], color = source_all), size = pointsize) + 
  scale_color_manual(values = c("goldenrod3", "darkgreen")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position = "none")
bplot

bplot = ggplot() +
  theme(plot.title = element_text(size=fontsize)) +
  geom_point(aes(x = sr_coords[,1], y = sr_coords[,2], color = sr_data$seurat_clusters), size = pointsize) + 
  scale_color_manual(values = c("royalblue4", "red4", "goldenrod1", "dodgerblue1", "grey60", "burlywood1", "violet", "red3", "darkgoldenrod")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position = "none")
bplot
dev.off()


# Original data frame for gene plots
steven_total_filt = steven_merged_all[rownames(sr_data@assays$RNA),colnames(sr_data@assays$RNA)]
dim(steven_total_filt)

fontsize = 8
pointsize =0.5
genePlotSteven <- function(gene,exp = steven_total_filt){
  sr_coords = sr_data@reductions$tsne@cell.embeddings
  col_by = log10(as.numeric(exp[gene,]))
  bplot = ggplot() +
    ggtitle(gene) +
    theme(plot.title = element_text(size=fontsize)) +
    geom_point(aes(x = sr_coords[,1], y = sr_coords[,2], color = col_by), size = pointsize) + 
    scale_color_gradient(low = "grey80", high = "maroon4", na.value = "grey80") +
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

# Extended Data Figure 5g
pdf(paste(plot_dir, "stevenGenePlots.rev.pdf", sep = "/"), height = 16, 8, useDingbats = FALSE) 
p1 = genePlotSteven("COL1A1")
p2 = genePlotSteven("COL3A1")
p3 = genePlotSteven("POSTN")
p4 = genePlotSteven("LUM")
p5 = genePlotSteven("FMOD")
p6 = genePlotSteven("COL9A1")
p7 = genePlotSteven("KRT19")
p11 = genePlotSteven("NEUROD6")
p12 = genePlotSteven("STMN2")
p14 = genePlotSteven("MAP2")
p15 = genePlotSteven("DLX5")
p16 = genePlotSteven("DLX6")
p17 = genePlotSteven("DLX1")
p18 = genePlotSteven("DLX2")
p19 = genePlotSteven("SATB2")
p20 = genePlotSteven("TBR1")
p21 = genePlotSteven("POU3F2")
p22 = genePlotSteven("EOMES")
p23 = genePlotSteven("CDK1")
p24 = genePlotSteven("TOP2A")
p25 = genePlotSteven("MKI67")
p26 = genePlotSteven("VIM")
p27 = genePlotSteven("PAX6")
p28 = genePlotSteven("FOXG1")
p29 = genePlotSteven("SOX9")
p30 = genePlotSteven("HOPX")
p31 = genePlotSteven("AGT")
p32 = genePlotSteven("FAM107A")
p34 = genePlotSteven("HEPACAM")
p35 = genePlotSteven("GFAP")
p36 = genePlotSteven("CRABP1")
p37 = genePlotSteven("MAL")
p38 = genePlotSteven("S100B")

grid.arrange(p1,p2,p3,p4,p5,p6,p7,
             p11,p12,p14,p15,p16,p17,p18,p19,p20,
             p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,
             p31,p32,p34,p35,p36,p37,p38, ncol = 4)


dev.off()


## Sub-clustering neuralectoderm-derived cells
head(select.cells)

neural = total[,select.cells]
dim(neural)

n_data = CreateSeuratObject(counts = neural, project = "neural", min.cells = 3, min.features = 200)
n_data

mitoGenes = geneInfo$gene[which(geneInfo$chrom=="chrMT")]
mitoGenes = mitoGenes[which(mitoGenes %in% rownames(n_data))]
n_data[["percent.mt"]] <- PercentageFeatureSet(n_data, features = mitoGenes)

# Don't filter because these have already been filtered in the first iteration
n_data <- NormalizeData(n_data, normalization.method = "LogNormalize", scale.factor = 10000)
n_data <- FindVariableFeatures(n_data, selection.method = "vst", nfeatures = 1500)
n_data
top10 <- head(VariableFeatures(n_data), 10)
all.genes <- rownames(n_data)
n_data <- ScaleData(n_data, features = all.genes)
n_data <- RunPCA(n_data, features = VariableFeatures(object = n_data))
print(n_data[["pca"]], dims = 1:5, nfeatures = 5)
ElbowPlot(n_data)


n_data <- FindNeighbors(n_data, dims = 1:12)
n_data <- FindClusters(n_data, resolution = 0.4)

n_data <- RunUMAP(n_data, dims = 1:14)
DimPlot(n_data, reduction = "umap")

n_data <- RunTSNE(n_data, dims = 1:14)
DimPlot(n_data, reduction = "tsne")

save(n_data, file = paste(analysis_dir, "neural_seurat_sc.rev.Rdata", sep = "/"))
load(paste(analysis_dir, "neural_seurat_sc.rev.Rdata", sep = "/"))
n_data.markers <- FindAllMarkers(n_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top100 <- n_data.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_logFC)

write.table(file = paste(analysis_dir, "top100ClusterGenes.neural.rev.txt", sep = "/"),
            data.frame(cbind(top100$cluster, top100$gene, top100$avg_logFC, top100$p_val_adj)),
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)


#pdf(paste(plot_dir, "neural_sc_clusterHeatmap.pdf", sep = "/"), height = 10, 8) 
#fontsize = 1
#DoHeatmap(n_data, features = top10$gene, size = 3, 
#          group.colors = c("royalblue4", "dodgerblue2", "red4", "cadetblue1", "grey60"),
#          raster = FALSE) + NoLegend()
#dev.off()

write.table(file = paste(analysis_dir, "neural_sc_umap_coords.rev.txt", sep = "/"), x = n_data@reductions$umap@cell.embeddings, col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(file = paste(analysis_dir, "neural_sc_tsne_coords.rev.txt", sep = "/"), x = n_data@reductions$tsne@cell.embeddings, col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

# Umap 
# Figure 2c
pdf(paste(plot_dir, "neural_sc_umap.rev.pdf", sep = "/"), height = 2.5, 2.5, useDingbats = FALSE) 
fontsize = 8
pointsize = 0.5
n_coords = n_data@reductions$umap@cell.embeddings

n_coords = read.table(paste(analysis_dir, "neural_sc_umap_coords.rev.txt", sep = "/"), header = TRUE, row.names = 1)


bplot = ggplot() +
  theme(plot.title = element_text(size=fontsize)) +
  geom_point(aes(x = n_coords[,1], y = n_coords[,2], color = n_data$seurat_clusters), size = pointsize) + 
  scale_color_manual(values = c("royalblue4", "dodgerblue2", "red4", "cadetblue3", "grey60")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position = "none")
bplot
dev.off()

# Meta data
neural_filt = neural[rownames(n_data@assays$RNA),colnames(n_data@assays$RNA)]
dim(neural_filt)
meta_filt_neural = sc_meta[which(colnames(total) %in% colnames(neural_filt)),]
dim(meta_filt_neural)

# Line
line_neural = meta_filt_neural$Line
hl2_16_sc_neural = neural_filt[, which(meta_filt_neural$Line=="HL2-16")]
hl2_9_sc_neural = neural_filt[, which(meta_filt_neural$Line=="HL2-9")]

dim(hl2_16_sc_neural)
dim(hl2_9_sc_neural)

# NOT USED
pdf(paste(plot_dir, "neural_sc_umap_Line.rev.pdf", sep = "/"), height = 1.3,1.3, useDingbats = FALSE) 
pointsize = 0.35
col = (line_neural == "HL2-9")
colbg = (line_neural == "HL2-16")
bplot = ggplot() +
  geom_point(aes(x = n_coords[,1], y = n_coords[,2]),size = pointsize) + 
  xlab("Dimension 1") + 
  ylab("Dimension 2")  + 
  geom_point(aes(n_coords[,1][colbg], n_coords[,2][colbg]), colour = "slateblue3",size = pointsize) + 
  geom_point(aes(n_coords[,1][col], n_coords[,2][col]), colour = "mediumseagreen",size = pointsize) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_blank(), panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position = "none")
print(bplot)
dev.off()

# Bar plot of cells per line per cluster
clustLine_neural = as.data.frame(cbind(line_neural, as.vector(n_data$seurat_clusters)))
names(clustLine_neural) = c("Line", "Cluster")
clustLine_neural.summary = group_by(clustLine_neural, Cluster, Line)
clustLine_neural.summary = count(clustLine_neural.summary)

dim(hl2_16_sc_neural) # 160
dim(hl2_9_sc_neural) # 208

clustLine_neural.summary$percent = clustLine_neural.summary$n / c(rep(160,5), rep(208, 5))

# Figure 2d
pdf(paste(plot_dir, "neural_Cluster_by_Line.rev.pdf", sep = "/"), height = 1.5, 1.8) 
fontsize = 8
p = ggplot(data = as.data.frame(clustLine_neural.summary), aes(x = as.factor(Cluster), y = percent, fill = Line)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylab("# of cells") +
  xlab("Cluster") +
  scale_fill_manual(values=c("slateblue3", "mediumseagreen"))+
  theme(title = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))+ theme(legend.title = element_blank()) 

p
dev.off()


# Dot Plot for marker genes
# Figure 2e
pdf(paste(plot_dir, "neural_sc_dotPlot.rev.pdf", sep = "/"), height = 2, 3.65, useDingbats = FALSE)
fontsize = 8
markers = as.vector(c("NEUROD6", "FAM49A", "DLX6", "DLX5", "SLC1A3", "CLU", "SCG2", "NEFL", "CDK1", "NUSAP1"))
DotPlot(n_data, features = markers, dot.min = 0.01, col.min = 0, col.max = 0, dot.scale = 3) +
  ylab("Cluster") +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 65, hjust = 1)) +
  theme(axis.title.y = element_text(size = fontsize)) +
  theme(axis.text = element_text(size = fontsize))
dev.off()


# Original data frame for gene plots
fontsize = 8
pointsize =0.5
genePlotNeural <- function(gene,exp = neural_filt){
  n_coords = n_data@reductions$umap@cell.embeddings
  col_by = log10(as.numeric(exp[gene,]))
  bplot = ggplot() +
    ggtitle(gene) +
    theme(plot.title = element_text(size=fontsize)) +
    geom_point(aes(x = n_coords[,1], y = n_coords[,2], color = col_by), size = pointsize) + 
    scale_color_gradient(low = "grey80", high = "maroon4", na.value = "grey80") +
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


# Figure 2f
pdf(paste(plot_dir, "neuralGenePlots.rev.pdf", sep = "/"), height = 10, 8, useDingbats = FALSE) 
p11 = genePlotNeural("NEUROD6")
p12 = genePlotNeural("STMN2")
p14 = genePlotNeural("MAP2")
p15 = genePlotNeural("DLX5")
p16 = genePlotNeural("DLX6")
p17 = genePlotNeural("DLX1")
p18 = genePlotNeural("DLX2")
p19 = genePlotNeural("SATB2")
p20 = genePlotNeural("TBR1")
p21 = genePlotNeural("POU3F2")
p22 = genePlotNeural("EOMES")
p23 = genePlotNeural("CDK1")
p24 = genePlotNeural("TOP2A")
p25 = genePlotNeural("MKI67")
p26 = genePlotNeural("VIM")
p27 = genePlotNeural("PAX6")
p28 = genePlotNeural("FOXG1")
p29 = genePlotNeural("SOX9")
p30 = genePlotNeural("HOPX")
p31 = genePlotNeural("AGT")
p32 = genePlotNeural("FAM107A")
p34 = genePlotNeural("HEPACAM")
p35 = genePlotNeural("GFAP")
p36 = genePlotNeural("SCG2")
p37 = genePlotNeural("NEFL")
p38 = genePlotNeural("GFRA1")

grid.arrange(p11,p12,p14,p15,p16,p17,p18,p19,p20,
             p21,p22,p23,p24,p25,p26,p27,p28,p29,p30,
             p31,p32,p34,p35,p36,p37,p38, ncol = 4)

dev.off()

# For later figures
# Extended Data Figure 8e, Figure 4b
pdf(paste(plot_dir, "neuralGenePlotsTOP.rev.pdf", sep = "/"), height = 1, 1.5, useDingbats = FALSE) 
pointsize = 0.3
genePlotNeural("PMP2")
genePlotNeural("SSTR2")
dev.off()

# Save cluster assignments for all cells, steven+all cells and neural cells
write.table(x = as.data.frame(data$seurat_clusters), file = paste(analysis_dir, "clusterAssignmentsAll.rev.txt", sep = "/"), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
write.table(x = as.data.frame(sr_data$seurat_clusters), file = paste(analysis_dir, "clusterAssignmentsSteven.rev.txt", sep = "/"), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
write.table(x = as.data.frame(n_data$seurat_clusters), file = paste(analysis_dir, "clusterAssignmentsNeural.rev.txt", sep = "/"), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")


# Files for cytoTrace
# filter cells but not genes
sc_filt = sc[, colnames(total_filt)]
sc_filt_cpm = sweep(as.matrix(sc_filt), 2, as.double(colSums(sc_filt)/1000000), `/`)
sc_filt_fpkm = sweep(as.matrix(sc_filt_cpm), 1, as.double(geneInfo$genesizeKb), `/`)
sc_filt_rpk = sweep(as.matrix(sc_filt), 1, as.double(geneInfo$genesizeKb), `/`)
sc_filt_tpm = sweep(as.matrix(sc_filt_rpk), 2, as.double(colSums(sc_filt_rpk)/1000000), `/`)


write.table(x = as.data.frame(sc_filt_tpm), file = paste(analysis_dir, "cytoTrace_total_tpm.txt", sep = "/"), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
cytoMeta = sc_meta[which(colnames(total_filt) %in% colnames(total)),]
cytoMeta$Cluster = data$seurat_clusters
rownames(cytoMeta) = colnames(sc_filt_tpm)
write.table(x = as.data.frame(as.data.frame(cytoMeta[,9], row.names = rownames(cytoMeta))), file = paste(analysis_dir, "cytoTrace_meta.txt", sep = "/"), quote = FALSE, row.names = TRUE, col.names = FALSE, sep = "\t")


# Human and chimp separately
sc_filt_human = human[, which(colnames(total) %in% colnames(total_filt))]
sc_filt_human_cpm = sweep(as.matrix(sc_filt_human), 2, as.double(colSums(sc_filt_human)/1000000), `/`)
sc_filt_human_fpkm = sweep(as.matrix(sc_filt_human_cpm), 1, as.double(geneInfo$genesizeKb[geneKeepIdx]), `/`)
sc_filt_human_rpk = sweep(as.matrix(sc_filt_human), 1, as.double(geneInfo$genesizeKb[geneKeepIdx]), `/`)
sc_filt_human_tpm = sweep(as.matrix(sc_filt_human_rpk), 2, as.double(colSums(sc_filt_human_rpk)/1000000), `/`)


sc_filt_chimp = chimp[, which(colnames(total) %in% colnames(total_filt))]
sc_filt_chimp_cpm = sweep(as.matrix(sc_filt_chimp), 2, as.double(colSums(sc_filt_chimp)/1000000), `/`)
sc_filt_chimp_fpkm = sweep(as.matrix(sc_filt_chimp_cpm), 1, as.double(geneInfo$genesizeKb[geneKeepIdx]), `/`)
sc_filt_chimp_rpk = sweep(as.matrix(sc_filt_chimp), 1, as.double(geneInfo$genesizeKb[geneKeepIdx]), `/`)
sc_filt_chimp_tpm = sweep(as.matrix(sc_filt_chimp_rpk), 2, as.double(colSums(sc_filt_chimp_rpk)/1000000), `/`)

write.table(x = as.data.frame(sc_filt_human_tpm), file = paste(analysis_dir, "cytoTrace_human_all_tpm.txt", sep = "/"), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
write.table(x = as.data.frame(sc_filt_chimp_tpm), file = paste(analysis_dir, "cytoTrace_chimp_all_tpm.txt", sep = "/"), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")


cytoMeta_human = cbind(colnames(sc_filt_human_tpm), as.numeric(as.character(data$seurat_clusters)))
cytoMeta_chimp = cbind(colnames(sc_filt_chimp_tpm), as.numeric(as.character(data$seurat_clusters)))   

write.table(x = cytoMeta_human, file = paste(analysis_dir, "cytoTrace_meta_human_all.txt", sep = "/"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(x = cytoMeta_chimp, file = paste(analysis_dir, "cytoTrace_meta_chimp_all.txt", sep = "/"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")



# Neural cells only, human and chimp
sc_neural = sc[, rownames(neural_cells)]
sc_neural_human = human[which(colnames(total) %in% rownames(neural_cells))]
sc_neural_human_cpm = sweep(as.matrix(sc_neural_human), 2, as.double(colSums(sc_neural_human)/1000000), `/`)
sc_neural_human_fpkm = sweep(as.matrix(sc_neural_human_cpm), 1, as.double(geneInfo$genesizeKb[geneKeepIdx]), `/`)
sc_neural_human_rpk = sweep(as.matrix(sc_neural_human), 1, as.double(geneInfo$genesizeKb[geneKeepIdx]), `/`)
sc_neural_human_tpm = sweep(as.matrix(sc_neural_human_rpk), 2, as.double(colSums(sc_neural_human_rpk)/1000000), `/`)

sc_neural_chimp = chimp[which(colnames(total) %in% rownames(neural_cells))]
sc_neural_chimp_cpm = sweep(as.matrix(sc_neural_chimp), 2, as.double(colSums(sc_neural_chimp)/1000000), `/`)
sc_neural_chimp_fpkm = sweep(as.matrix(sc_neural_chimp_cpm), 1, as.double(geneInfo$genesizeKb[geneKeepIdx]), `/`)
sc_neural_chimp_rpk = sweep(as.matrix(sc_neural_chimp), 1, as.double(geneInfo$genesizeKb[geneKeepIdx]), `/`)
sc_neural_chimp_tpm = sweep(as.matrix(sc_neural_chimp_rpk), 2, as.double(colSums(sc_neural_chimp_rpk)/1000000), `/`)


write.table(x = as.data.frame(sc_neural_human_tpm), file = paste(analysis_dir, "cytoTrace_human_tpm.txt", sep = "/"), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")
write.table(x = as.data.frame(sc_neural_chimp_tpm), file = paste(analysis_dir, "cytoTrace_chimp_tpm.txt", sep = "/"), quote = FALSE, row.names = TRUE, col.names = TRUE, sep = "\t")



cytoMeta_neural_human = cbind(colnames(sc_neural_human_tpm), as.numeric(as.character(n_data$seurat_clusters)))
cytoMeta_neural_chimp = cbind(colnames(sc_neural_chimp_tpm), as.numeric(as.character(n_data$seurat_clusters)))   
  
write.table(x = cytoMeta_neural_human, file = paste(analysis_dir, "cytoTrace_meta_human.txt", sep = "/"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(x = cytoMeta_neural_chimp, file = paste(analysis_dir, "cytoTrace_meta_chimp.txt", sep = "/"), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

  
  
