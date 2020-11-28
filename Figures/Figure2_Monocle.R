## Analysis and plots for Figure 2g-i and Extended Data Figure 3h: monocle analysis
## Everything here was performed in Rstudio on a Macbook Pro with 16G RAM, unless otherwise noted
## Clear cache in Rstudio
rm(list = ls(all.names = TRUE))

## Load required libraries
library(dplyr)
library(reshape)
library(scales)
library(limma)
library(Matrix)
library(tidyverse)
library(monocle)
library(gridExtra)

peek <- function(data){
  rowmax = min(5, nrow(data))
  colmax = min(5, ncol(data))
  print(data[1:rowmax, 1:colmax])
}

## Directories where the data is / where to deposit files and plots
data_dir = "/Users/rachelagoglia/Desktop/Revisions/ProcessedData"
analysis_dir = "/Users/rachelagoglia/Desktop/Revisions/AnalysisFiles"
plot_dir = "/Users/rachelagoglia/Desktop/Revisions/Plots"

## Load data (counts matrix, meta data)
neural_cells = read.table(paste(analysis_dir, "seurat.neural.cells.txt",sep = "/"))
head(neural_cells)

# Gene meta data
geneInfo = read.table(paste(data_dir, "geneInfo.GRCh38.txt", sep = "/"), row.names = 1)
names(geneInfo) = c("length", "chrom")
geneInfo$genesizeKb = as.vector(geneInfo$length/1000)
geneInfo$gene_short_name = rownames(geneInfo)
dim(geneInfo)
head(geneInfo)

sc = read.table(paste(data_dir, "scRNAseq.GRCh38.txt", sep = "/"), header = TRUE)
rownames(sc) = rownames(geneInfo)
dim(sc)

# Genes and chromosomes to remove from subsequent analyses
geneKeepInfo = read.table(file = paste(analysis_dir, "geneKeep.txt", sep = "/"), header = TRUE, fill = TRUE)
geneKeep = as.vector(geneKeepInfo$geneKeep)
geneKeepIdx = as.vector(geneKeepInfo$geneKeepIdx)

# Calculate TPM and CPM
sc_keep = sc[geneKeep,]
dim(sc_keep)
sc_cpm = sweep(as.matrix(sc_keep), 2, as.double(colSums(sc_keep)/1000000), `/`)
sc_fpkm = sweep(as.matrix(sc_cpm), 1, as.double(geneInfo$genesizeKb[geneKeepIdx]), `/`)

sc_neural = sc_fpkm[,as.character(rownames(neural_cells))]
dim(sc_neural)

meta = data.frame(cbind(as.character(rownames(neural_cells)), neural_cells$n_data.seurat_clusters))
rownames(meta) = meta$X1
colnames(meta) = c("cell_ID", "cell_type")
head(meta)


geneInfoKeep = geneInfo[geneKeepIdx,]
dim(geneInfoKeep)

# Create CellDataSet object
pd <- new("AnnotatedDataFrame", data = meta)
fd <- new("AnnotatedDataFrame", data = geneInfoKeep)
cds <- newCellDataSet(as.matrix(sc_neural), 
                      phenoData = pd, 
                      featureData = fd, 
                      lowerDetectionLimit = 0.1,
                      expressionFamily = tobit(Lower = 0.1)) 

# Convert fpkm to mRNA transcript counts
rpc_matrix <- relative2abs(cds, method = "num_genes")

# Create new CellDataSet object with the RNA counts
rcds = newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

# Estimate size factors and dispersions
rcds <- estimateSizeFactors(rcds)
rcds <- estimateDispersions(rcds)

# Detect genes that are expressed
rcds <- detectGenes(rcds, min_expr = 0.1)
head(fData(rcds))

expressed_genes <- row.names(subset(fData(rcds),
                                    num_cells_expressed >= 10))

# Assess expression levels of the cells
pData(rcds)$Total_mRNAs <- Matrix::colSums(exprs(rcds))

rcds <- rcds[,pData(rcds)$Total_mRNAs < 1e6]

upper_bound <- 10^(mean(log10(pData(rcds)$Total_mRNAs)) +
                     2*sd(log10(pData(rcds)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(rcds)$Total_mRNAs)) -
                     2*sd(log10(pData(rcds)$Total_mRNAs)))

qplot(Total_mRNAs, data = pData(rcds), 
      color = cell_type, 
      geom = "density") +
  geom_vline(xintercept = lower_bound) +
  geom_vline(xintercept = upper_bound)

rcds <- rcds[,pData(rcds)$Total_mRNAs > lower_bound &
               pData(rcds)$Total_mRNAs < upper_bound]

# Save list of which cells to keep in analysis 
monocle_cells = pData(rcds)$cell_ID
length(monocle_cells)

# Verify that expression values follow a roughly lognormal distribution
L <- log(exprs(rcds[expressed_genes,]))
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
qplot(value, geom = "density", data = melted_dens_df) +
  stat_function(fun = dnorm, size = 0.5, color = 'red') +
  xlab("Standardized log(FPKM)") +
  ylab("Density")


# Cluster cells without marker genes
disp_table <- dispersionTable(rcds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
rcds <- setOrderingFilter(rcds, unsup_clustering_genes$gene_id)
plot_ordering_genes(rcds)

plot_pc_variance_explained(rcds, return_all = F)
rcds <- reduceDimension(rcds, max_components = 2, num_dim = 4,
                        reduction_method = 'tSNE', verbose = T)

rcds <- clusterCells(rcds, num_clusters = 2) # Distance cutoff 1.156351
plot_cell_clusters(rcds, 1, 2)

# Choose genes to use for pseudotime ordering 
fData(rcds)$use_for_ordering <- fData(rcds)$num_cells_expressed > 0.05 * ncol(rcds)

plot_pc_variance_explained(rcds, return_all = F)

rcds <- reduceDimension(rcds,
                        max_components = 2,
                        norm_method = 'log',
                        num_dim = 4,
                        reduction_method = 'tSNE',
                        verbose = T)

rcds <- clusterCells(rcds, verbose = F)
plot_cell_clusters(rcds, color_by = 'as.factor(Cluster)')
plot_cell_clusters(rcds, color_by = 'as.factor(cell_type)')

plot_rho_delta(rcds, rho_threshold = 2, delta_threshold = 4 )
rcds <- clusterCells(rcds,
                     rho_threshold = 6,
                     delta_threshold = 10,
                     skip_rho_sigma = T,
                     verbose = F)

plot_cell_clusters(rcds, color_by = 'as.factor(Cluster)')
plot_cell_clusters(rcds, color_by = 'as.factor(cell_type)')

# Takes a few minutes
clustering_DEG_genes <- differentialGeneTest(rcds[expressed_genes,],
                                             fullModelFormulaStr = '~Cluster',
                                             cores = 1)


rcds_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]

rcds <- setOrderingFilter(rcds, ordering_genes = rcds_ordering_genes)
rcds <- reduceDimension(rcds, method = 'DDRTree')
rcds <- orderCells(rcds)


# Plots
plot_cell_trajectory(rcds, color_by = "State")

# Redefine the root of the tree to be the progenitor cells
rcds <- orderCells(rcds, root_state = 3)

# Figure 2g
pdf(paste(plot_dir, "monocle.rev.pdf", sep = "/"), height = 2.5, width = 2.5, useDingbats = FALSE)
fontsize = 8
pointsize = 0.5
plot_cell_trajectory(rcds, color_by = "Pseudotime", cell_size = pointsize, 
                     show_branch_points = FALSE,
                     cell_link_size = pointsize) +
  theme(legend.position="none") +
  theme(plot.title = element_blank())+
  theme(axis.text=element_blank(), axis.title=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank()) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank()) 

plot_cell_trajectory(rcds, color_by = "Pseudotime", cell_size = pointsize, 
                     show_branch_points = FALSE,
                     cell_link_size = pointsize) +
  theme(legend.position="bottom", legend.box = "horizontal") +
  theme(plot.title = element_blank())+
  theme(axis.text=element_blank(), axis.title=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank()) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank())

plot_cell_trajectory(rcds, color_by = "cell_type", cell_size = pointsize, 
                     show_branch_points = FALSE,
                     cell_link_size = pointsize) +
  scale_color_manual(values = c("royalblue4", "dodgerblue2", "red4", "cadetblue3", "grey60")) +
  theme(legend.position="none") +
  theme(plot.title = element_blank())+
  theme(axis.text=element_blank(), axis.title=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank()) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank()) 
dev.off()


# Plot salmon genes eigengene
mod_genes =colnames(net_hcCS_filt)[which(moduleColors == "salmon")]
mono_sc_mod = total_filt[which(rownames(total_filt) %in% mod_genes), which(colnames(total_filt) %in% rcds$cell_ID)]
mono_sc_mod_comp = mono_sc_mod[complete.cases(mono_sc_mod),]
mono_sc_mod_1 = mono_sc_mod_comp[which(rowSums(mono_sc_mod_comp)>1),]
mod.pca = prcomp(t(mono_sc_mod_1), scale=TRUE, center=TRUE)
mod.pca_data = as.data.frame(mod.pca$x)
eig = mod.pca_data$PC1
eig = -1*eig

# Drop one outlier value
rcds$salmon_eig = log(eig +1.4)
rcds$salmon_eig = rank(eig)

# Figure 3h
pdf(paste(plot_dir, "monocle.salmon.legend.rev.pdf", sep = "/"), height = 2.5, width = 2.5, useDingbats = FALSE)
pointsize = 0.3
plot_cell_trajectory(rcds, color_by = "salmon_eig", cell_size = pointsize, 
                     show_branch_points = FALSE,
                     cell_link_size = pointsize) +
  scale_color_gradient(low = "grey80", high = "purple4", na.value = "grey80") +
  theme(plot.title = element_blank())+
  theme(axis.text=element_blank(), axis.title=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank()) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank()) 
dev.off()



## not used
dim(total_filt)
dim(human_filt)
dim(chimp_filt)

ase_filt = log2(human_filt/chimp_filt)

ase_salmon = ase_filt[mod_genes, which(colnames(total_filt) %in% rcds$cell_ID)]
dim(ase_salmon)
test = colMeans(as.matrix(ase_salmon), na.rm = T)
hist(test)
rcds$test = test
plot_cell_trajectory(rcds, color_by = "test", cell_size = pointsize, 
                     show_branch_points = FALSE,
                     cell_link_size = pointsize) +
  scale_color_gradient(low = "blue", high = "red", na.value = "grey80") +
  theme(plot.title = element_blank())+
  theme(axis.text=element_blank(), axis.title=element_blank()) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x = element_blank()) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank()) +
  theme(legend.position = "none")


dim(ase_filt)

ase_filt_neural = ase_filt[, which(colnames(total_filt) %in% colnames(neural_filt))]
dim(ase_filt_neural)

tests = c()
results = c()
for(i in seq(1, nrow(ase_filt))){
  x = as.numeric(ase_filt[i, which(neural_clusters$n_data.seurat_clusters == 0)])
  y = as.numeric(ase_filt[i, which(neural_clusters$n_data.seurat_clusters == 2)])
  if(length(which(is.finite(x))) > 5 & length(which(is.finite(y))) > 5){
    test = wilcox.test(x, y)
    tests = c(tests, "")
    results = c(results, test$p.value)
  }

}
length(tests)
length(which(results < 0.05))

# Make a heatmap of branch dependent genes
#BEAM_res <- BEAM(rcds, branch_point = 1, cores = 1) # Takes forever!
#BEAM_res <- BEAM_res[order(BEAM_res$qval),]
#BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

#pdf(paste(plot_dir, "monocleBEAM.pdf", sep = "/"), height = 7, width = 4)
#fontsize = 8
#plot_genes_branched_heatmap(rcds[row.names(subset(BEAM_res,
#                                                  qval < 1e-6)),],
#                            branch_point = 1,
#                            num_clusters = 4,
#                            cores = 1,
#                            use_gene_short_name = T,
#                            show_rownames = T)
#dev.off()

# Save the BEAM data so we don't have to run this again
#save(BEAM_res, file = paste(analysis_dir, "monocleBEAM_res.rev.RData", sep = "/"))

# Plot specific branching genes

# Figure 2h,i
pdf(paste(plot_dir, "monocle.branchedGenes.extra.rev.pdf", sep = "/"), height = 13, width = 2, useDingbats = FALSE)
fontsize = 8
pointsize = 0.5
branch_genes <- row.names(subset(fData(rcds),
                                 gene_short_name %in% c("TOP2A",  "MKI67", "CDK1","SLC1A3", "NEUROD6", "SLA", "SATB2", "SPARC", "AGT", "TBR1","DLX5", "DLX6", "EOMES", "FEZF2", "BCL11B", "RORB", "TLE3", "TLE1", "POU3F2", "NEUROD6", "SLA", "SATB2", "SSTR2")))

plot_genes_branched_pseudotime(rcds[branch_genes,],
                               branch_point = 1,
                               color_by = "cell_type",
                               ncol = 1,
                               cell_size = pointsize) +
  scale_color_manual(values = c("royalblue4", "dodgerblue2", "red4", "cadetblue3", "grey60")) +
  theme(plot.title = element_blank())+
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize)) +
  theme(legend.position = "none")

dev.off()
