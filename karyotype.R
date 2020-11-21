#!/usr/bin/env Rscript
library(ggplot2)
library(gridExtra)

args = commandArgs(TRUE)
aseFile = args[1]
posFile = args[2]
pdfFile = args[3]

libSize = as.numeric(scan(aseFile, what="character", n=4)[4])

data = read.table(aseFile, header=TRUE)
pos = read.table(posFile, header=FALSE)

chroms = levels(data$chrom)
total_exp = data$ref_counts+data$alt_counts+data$no_ase_counts+data$ambig_ase_counts

# Filter out lowly expressed genes
data = data[total_exp>20,]
pos = pos[total_exp>20,3]

# Calculate ASE
ase = log2(data$ref_counts/data$alt_counts)

# Plot each chromosome
pdf(pdfFile)

for(i in seq(1,length(chroms))){
  chromo = chroms[i]
  chr_ase = ase[data$chrom==chromo]
  chr_pos= pos[data$chrom==chromo]
  chr_data = as.data.frame(cbind(chr_ase, chr_pos))
  
  chr_ord = chr_data[order(chr_data$chr_pos),]
  
  mean_ase = mean(chr_ord$chr_ase[is.finite(chr_ord$chr_ase)])
  
  p1 = ggplot(chr_ord, aes(chr_pos, chr_ase)) + 
    geom_point(col = "brown3", size = 0.8)  +  
    ggtitle(paste("ASE: chromosome", chromo))+
    ylab("Log2(ASE Ratio)") +
    xlab("") +
    theme(legend.position="none") +
    theme(axis.title.x=element_blank()) +
    theme(axis.text=element_text(size=9), axis.title=element_text(size=11))  +
    geom_hline(yintercept = 0) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    geom_hline(yintercept = mean_ase, linetype="longdash") +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
      
  
  # Calculate average, median ASE in a sliding window of 20 genes
  means = c()
  medians = c()
  positions = c()
  stats = c()
  window = 20
  
  if(length(chr_ase) < 25) next
  
  for(i in seq(1, (length(chr_ord$chr_ase)-(window-1)))){
    values = chr_ord$chr_ase[i:(i+(window-1))]
    position = chr_ord$chr_pos[i]
    means = c(means, mean(values, na.rm=TRUE))
    medians = c(medians, median(values, na.rm=TRUE))
    stats = c(stats, tryCatch(wilcox.test(ase, values)$p.value, error=function(e) NA))
    positions = c(positions, position)
  }
  
  tempdf = as.data.frame(cbind(positions,means, medians))
  
  p2 = ggplot(tempdf, aes(positions, means)) + 
    geom_point(col = "orange2", size = 0.8)  +  
    ggtitle("Averages")+
    ylab("Average Log2(ASE Ratio)") +
    xlab("") +
    ylim(c(-2,2)) +
    theme(legend.position="none") +
    theme(axis.title.x=element_blank()) +
    theme(axis.text=element_text(size=9), axis.title=element_text(size=11))  +
    geom_hline(yintercept = 0) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    geom_hline(yintercept = mean_ase, linetype="longdash") +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
      
  
  p3 = ggplot(tempdf, aes(positions, medians)) + 
    geom_point(col = "seagreen4", size = 0.8)  +  
    ggtitle("Medians")+
    theme(axis.title.x=element_blank()) +
    ylab("Median Log2(ASE Ratio)") +
    ylim(c(-2,2)) +
    xlab("") +
    theme(legend.position="none") +
    theme(axis.text=element_text(size=9), axis.title=element_text(size=11))  +
    geom_hline(yintercept = 0) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    geom_hline(yintercept = mean_ase, linetype="longdash") +  
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
      
  
  bonf = -log10(0.05/length(ase))
  top = bonf
  if(max(-log10(stats), na.rm=TRUE)> bonf){
    top=max(-log10(stats), na.rm=TRUE)
  }
  
  p4 = ggplot(tempdf, aes(positions, -log10(stats))) + 
    geom_point(col = "royalblue3", size = 0.8)  +  
    ggtitle("Wilcoxon test")+
    xlab("Genomic Position") +
    ylab("-Log10(p-value)") +
    ylim(c(0, top+1)) +
    theme(legend.position="none") +
    theme(axis.text=element_text(size=9), axis.title=element_text(size=11))  +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = bonf, linetype="longdash") +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
      
  
  # Display three plots (to include averages, add p2 as well)
  grid.arrange(p1, p3, p4, nrow = 3)
  
  write(paste("Done with chrom:", chromo), stderr())
}

dev.off()
