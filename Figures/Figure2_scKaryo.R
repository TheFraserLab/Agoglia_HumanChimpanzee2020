## Analysis and plots for Extended Data Figure 4h-j
## Everything here was performed in Rstudio on a Macbook Pro with 16G RAM, unless otherwise noted

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


# Sngle cell RNAseq data
sc = read.table(paste(data_dir, "scRNAseq.GRCh38.txt", sep = "/"), header = TRUE)
rownames(sc) = geneInfo$gene
dim(sc)
sc_meta = read.table(paste(data_dir, "scRNAseq.meta.txt", sep = "/"), header = TRUE)
sc_meta$Line = as.character(sc_meta$Line)

# Allelic and total reads
human_base = sc[, seq(2,length(sc),by=3)]
chimp_base = sc[, seq(3,length(sc),by=3)]
total_base = sc[, seq(1,length(sc),by=3)]

# Remove cells that don't pass filtering
filtered = read.table(paste(analysis_dir, "total_filt.rev.txt", sep = "/"))
keepers = which(colnames(total_base) %in% colnames(filtered))

# Look at percent human reads for each cell
wholeGenome = (colSums(human_base[,keepers]) / (colSums(human_base[,keepers]) + colSums(chimp_base[,keepers])) ) *100
qplot(wholeGenome) + 
  geom_histogram()

# Extended Data Figure 4h, inset bottom right
pdf(paste(plot_dir, "scKaryo_WG.rev.pdf", sep = "/"), height = 1.07, width = 1.2, useDingbats = FALSE) 
binsize = 0.5
fontsize = 6
qplot(wholeGenome, 
      geom = "histogram", 
      main="", 
      xlab="% Human Reads (whole genome)",
      ylab="# of cells",
      fill=I("mediumpurple2"), 
      col=I("black"), 
      binwidth=binsize) +
  theme_bw() + 
  theme(axis.text = element_text(size = fontsize), axis.title = element_text(size = fontsize), title = element_text(size = fontsize)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept = 0) 
dev.off()

max(wholeGenome)
min(wholeGenome)

# Allelic bias per chromosome per cell
library(dplyr)

totals = c()
totPerChrom = data.frame(init = rep("t", 706))
cols = data.frame(init = rep("t", 706))
for(i in c(seq(1, 22), "X")){
  total = sum(rowSums(total_base[which(geneInfo$chrom == paste("chr", as.character(i), sep = "")),keepers]))
  totals = c(totals, total)
  chr_human = human_base[which(geneInfo$chrom == paste("chr", as.character(i), sep = "")),keepers]
  chr_chimp = chimp_base[which(geneInfo$chrom == paste("chr", as.character(i), sep = "")),keepers]
  totChrom = colSums(chr_human) + colSums(chr_chimp)
  chr_percent = (colSums(chr_human) / (colSums(chr_human) + colSums(chr_chimp)) ) *100
  cols[paste("chr", as.character(i), sep = "")] = chr_percent
  totPerChrom[paste("chr", as.character(i), sep = "")] = as.vector(totChrom)

}

dim(cols)
head(cols)

dim(totPerChrom)
head(totPerChrom)


cols$line = sc_meta$Line[keepers]
cols = arrange(cols, cols$line)

max(totals)
totals_norm = round((totals/max(totals))*200)

totColors = colorRampPalette(c("white", "black"))(n = 200)
totColors = totColors[totals_norm]

# Extended Data Figure 4h, bottom left
pdf(paste(plot_dir, "scKaryo_perCell_key.rev.pdf", sep = "/"), height = 7, width = 10)
hcolors = colorRampPalette(c("blue", "white", "darkred"))(n = 200)
heatmap.2(as.matrix(cols[, 2:24]), dendrogram = "none", col = hcolors, tracecol = "black", Rowv = F, Colv = F, RowSideColors = c(rep("slateblue3", 338), rep("mediumseagreen", 368)))
dev.off()

# Extended Data Figure 4h
pdf(paste(plot_dir, "scKaryo_perCell.rev.pdf", sep = "/"), height = 2, width = 3)
hcolors = colorRampPalette(c("blue", "white", "darkred"))(n = 200)
heatmap.2(as.matrix(cols[, 2:24]), dendrogram = "none", col = hcolors, tracecol = NA, Rowv = F, Colv = F, ColSideColors = totColors, RowSideColors = c(rep("slateblue3", 338), rep("mediumseagreen", 368)), key = FALSE)
dev.off()

# Extended Data Figure 4h, top right key
pdf(paste(plot_dir, "scKaryo_perCell_countKey.rev.pdf", sep = "/"), height = 7, width = 10)
hcolors = colorRampPalette(c("white", "black"))(n = 200)
heatmap.2(as.matrix(cols[, 2:24]), dendrogram = "none", col = hcolors, tracecol = NA, Rowv = F, Colv = F, ColSideColors = totColors, RowSideColors = c(rep("slateblue3", 338), rep("mediumseagreen", 368)))
dev.off()


# Extended Data Figure 4i,j
pdf(paste(plot_dir, "PerChrom.pdf", sep = "/"), height = 1.5, width = 3)
for(i in seq(1, 22)){
  chr_human = human_base[which(geneInfo$chrom == paste("chr", as.character(i), sep = "")),keepers]
  chr_chimp = chimp_base[which(geneInfo$chrom == paste("chr", as.character(i), sep = "")),keepers]
  chr_percent = (colSums(chr_human) / (colSums(chr_human) + colSums(chr_chimp)) ) *100
  chr_line = sc_meta[keepers,]$Line
  
  data = data.frame(cbind(chr_percent, chr_line))
  p = ggplot(data = data, aes(x = as.numeric(as.character(chr_percent)), fill = chr_line)) + 
    geom_histogram(alpha = 0.5, position = 'identity', color = "black", size = 0.3) + 
    ggtitle(paste("chr", as.character(i), sep = "")) + 
    xlim(0, 100) +
    xlab("% of Reads") + 
    ylab("Frequency")  + 
    geom_vline(xintercept = 50, size = 0.5) +
    geom_vline(xintercept = mean(chr_percent[which(chr_line == "HL2-16")], na.rm = TRUE), size = 0.5, linetype = "dashed", color = "darkblue") +
    geom_vline(xintercept = mean(chr_percent[which(chr_line == "HL2-9")], na.rm = TRUE), size = 0.5, linetype = "dashed", color = "darkgreen") +
    scale_fill_manual(values = c("slateblue3", "mediumseagreen")) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA)) +
    theme(axis.text = element_text(size = fontsize), axis.title = element_text(size = fontsize), title = element_text(size = fontsize))
    

  print(p)
}
dev.off()
i=20

# Correlate variance in each chromosome with read depth
# Extended Data Figure 4h, bottom middle
pdf(paste(plot_dir, "PerChromVariance.pdf", sep = "/"), height = 0.9, width = 1.5, useDingbats = FALSE)
pointsize = 0.2
fontsize = 6
cor(colVars(as.matrix(cols[,2:24]), na.rm = TRUE), totals)
qplot(y=colVars(as.matrix(cols[,2:24]), na.rm = TRUE), x = totals/1000000, size = I(pointsize))+
  xlab("Total Reads (millions)") + 
  ylab("Variance") + 
  theme_bw() +
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(axis.text = element_text(size = fontsize, colour = "black"), axis.title = element_text(size = fontsize, color = "black"), title = element_text(size = fontsize))
dev.off()


# Same, but with individual measurements
pdf(paste(plot_dir, "PerChromVariance_perCell.pdf", sep = "/"), height = 0.9, width = 1.5, useDingbats = FALSE)
pointsize = 0.2
fontsize = 6
totPerChrom_long = melt(totPerChrom)
cols_long = melt(cols)
qplot(y=cols_long$value, x = totPerChrom_long$value/1000000, size = I(pointsize))+
  geom_point(aes(color=ifelse(cols_long$variable=="chr18", "red", "black"))) +
  xlab("Total Reads (millions)") + 
  ylab("% Human Reads") + 
  theme_bw() +
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent", colour = NA), plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(axis.text = element_text(size = fontsize, colour = "black"), axis.title = element_text(size = fontsize, color = "black"), title = element_text(size = fontsize))
dev.off()
