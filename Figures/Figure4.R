## Analysis and plots for Figure 4a,b,e,f,i,k,l and Extended Data Figures 8e-g, 9a-i,k, 10a-b, d-e: ASE genes and SSTR2
## Everything here was performed in Rstudio on a Macbook Pro with 16G RAM, unless otherwise noted
## Clear cache in Rstudio
rm(list = ls(all.names = TRUE))

## Load required libraries
library(ggplot2)
library(gridExtra)
library(gplots)
library(dplyr)
library(VennDiagram)
library(ggrepel)
library(reshape)

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


# Gene keep info 
geneKeepInfo = read.table(file = paste(analysis_dir, "geneKeep.rev.txt", sep = "/"), header = TRUE, fill = TRUE)
geneKeep = as.vector(geneKeepInfo$geneKeep)
geneKeepIdx = as.vector(geneKeepInfo$geneKeepIdx)


# Cis / trans for spheroid data
hyCS_deseq = read.table(paste(analysis_dir, "hyCS.DESEq2.ASE.rev.txt", sep = "/"), header = TRUE, row.names = 1)
CS_par_deseq = read.table(paste(analysis_dir, "all_par_DESeq2.rev.txt", sep = "/"), sep = "\t", fill = TRUE, na.strings = "", header = TRUE)
rownames(CS_par_deseq) = geneInfo$gene
testgenes = intersect(intersect(rownames(hyCS_deseq), rownames(CS_par_deseq)), geneKeep)

cis = hyCS_deseq$D150_HumanLFC[which(rownames(hyCS_deseq) %in% testgenes)]
trans = CS_par_deseq$Day150_LFC[which(rownames(CS_par_deseq) %in% testgenes)]

# Extended Data Figure 9e
pdf(paste(plot_dir, "CS_cis_trans.rev.pdf", sep = "/"), height = 2.5, 2.5, useDingbats = FALSE) 
fontsize = 8
ggplot() +  
  ggtitle("Cis/Trans Test") +
  ylab("Log2(Human / Chimp)") +
  xlab("Log2(Hybrid ASE)") + 
  geom_point(aes(y = trans, x = cis), size = 0.1) +
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
cor(cis[filt], trans[filt]) * cor(cis[filt], trans[filt]) 

# percent cis/trans
percentCis = abs(cis)/(abs(cis) + abs(trans-cis))
length(percentCis)
mean(percentCis, na.rm = TRUE) # 39%


# Color points by whether they are mesenchyme specific
percents = read.table(paste(analysis_dir, "scCellTypeStats.rev.txt", sep = "/"), header = TRUE, row.names = 1)
mes_genes = rownames(percents)[which(percents$mes_percent >=90)]
neural_genes = rownames(percents)[which(percents$astro_percent + percents$neuron_percent + percents$progenitor_percent >=90)]

# NOT USED
pdf(paste(plot_dir, "CS_cis_trans.mes.rev.pdf", sep = "/"), height = 2, 2, useDingbats = FALSE) 
fontsize = 8
ggplot() +  
  ggtitle("Cis/Trans Test") +
  ylab("Log2(Human / Chimp)") +
  xlab("Log2(Hybrid ASE)") + 
  geom_point(aes(y = trans, x = cis), size = 0.04, col = "white", alpha = 0.3) +
  geom_point(aes(y = trans[which(testgenes %in% mes_genes)], x = cis[which(testgenes %in% mes_genes)]), size = 0.04, col = "goldenrod3", alpha = 0.3) +
  #geom_point(aes(y = trans[which(testgenes %in% neural_genes)], x = cis[which(testgenes %in% neural_genes)]), size = 0.04, col = "dodgerblue4", alpha = 0.3) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlim(-10,10) +
  ylim(-10,10) +
  theme_bw() +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_abline(slope=1, intercept = 0)
dev.off()

# Venn diagrams
# iPS data
geneKeepInfoIPS = read.table(file = paste(analysis_dir, "IPSgeneKeep.txt", sep = "/"), header = TRUE, fill = TRUE)
geneKeepIPS = as.vector(geneKeepInfo$geneKeep)
geneKeepIdxIPS = as.vector(geneKeepInfo$geneKeepIdx)

iPS_deseq_HC = read.table(paste(analysis_dir, "iPS.DESeq2.merge.highConfASE.txt", sep = "/"), header = TRUE, row.names = 1, fill = TRUE)
iPS_deseq_HC = iPS_deseq_HC[which(rownames(iPS_deseq_HC) %in% geneKeepIPS),]
head(iPS_deseq_HC)
dim(iPS_deseq_HC)

# hyCS data
hyCS_deseq_HC = read.table(paste(analysis_dir, "hyCS.DESEq2.ASE.highConf.rev.txt", sep = "/"), header = TRUE, row.names = 1, fill = TRUE)
hyCS_deseq_HC = hyCS_deseq_HC[which(rownames(hyCS_deseq_HC) %in% geneKeep),]
head(hyCS_deseq_HC)
dim(hyCS_deseq_HC)

# Filter parent CS data for same genes in hyCS high conf. set
hcCS_d150_deseq_keep = CS_par_deseq[which(rownames(CS_par_deseq) %in% rownames(hyCS_deseq_HC)),]
dim(hcCS_d150_deseq_keep)

hcCS_d150_deseq_sig = rownames(hcCS_d150_deseq_keep)[which(hcCS_d150_deseq_keep$Day150_padj < 0.01)]
iPS_sig = rownames(iPS_deseq_HC)[which(iPS_deseq_HC$HumanPadj<0.01 & iPS_deseq_HC$ChimpPadj<0.01)]
d50_sig = rownames(hyCS_deseq_HC)[which(hyCS_deseq_HC$D50_HumanPadj<0.01 & hyCS_deseq_HC$D50_ChimpPadj<0.01)]
d100_sig = rownames(hyCS_deseq_HC)[which(hyCS_deseq_HC$D100_HumanPadj<0.01 & hyCS_deseq_HC$D100_ChimpPadj<0.01)]
d150_sig = rownames(hyCS_deseq_HC)[which(hyCS_deseq_HC$D150_HumanPadj<0.01 & hyCS_deseq_HC$D150_ChimpPadj<0.01)]

# Extended Data Figure 9a-d
pdf(paste(plot_dir, "VennDiagrams.pdf", sep = "/"), height = 4, width = 4.5)
vennPlot = venn.diagram(list("Day 50" = 1:length(d50_sig)-1, "Day 100" = (length(d50_sig)-length(intersect(d50_sig, d100_sig))):((length(d50_sig) + length(d100_sig))-length(intersect(d50_sig, d100_sig)) -1)),fill = c("forestgreen", "dodgerblue2"), fontfamily = "Helvetica", cex=1.5, cat.fontfamily = "Helvetica", cat.cex = 1.5, filename=NULL, cat.pos = c(220, 150), cat.dist = 0.05)
grid.newpage()
grid.draw(vennPlot)

vennPlot = venn.diagram(list("Day 100" = 1:length(d100_sig)-1, "Day 150" = (length(d100_sig)-length(intersect(d100_sig, d150_sig))):((length(d100_sig) + length(d150_sig))-length(intersect(d100_sig, d150_sig)) -1)),fill = c("dodgerblue2", "darkslateblue"), fontfamily = "Helvetica", cex=1.5, cat.fontfamily = "Helvetica", cat.cex = 1.5, filename=NULL, cat.pos = c(220, 150), cat.dist = 0.05)
grid.newpage()
grid.draw(vennPlot)

vennPlot = venn.diagram(list("Day 150" = 1:length(d150_sig)-1, "iPSC" = (length(d150_sig)-length(intersect(d150_sig, iPS_sig))):((length(d150_sig) + length(iPS_sig))-length(intersect(d150_sig, iPS_sig)) -1)),fill = c("darkslateblue", "gray50"), fontfamily = "Helvetica", cex=1.5, cat.fontfamily = "Helvetica", cat.cex = 1.5, filename=NULL, cat.pos = c(330, 40), cat.dist = 0.05)
grid.newpage()
grid.draw(vennPlot)

vennPlot = venn.diagram(list("Hybrid" = 1:length(d150_sig)-1, "Parent" = (length(d150_sig)-length(intersect(d150_sig, hcCS_d150_deseq_sig))):((length(d150_sig) + length(hcCS_d150_deseq_sig))-length(intersect(d150_sig, hcCS_d150_deseq_sig)) -1)),fill = c("darkslateblue", "slateblue"), fontfamily = "Helvetica", cex=1.5, cat.fontfamily = "Helvetica", cat.cex = 1.5, filename=NULL, cat.pos = c(330, 40), cat.dist = 0.05)
grid.newpage()
grid.draw(vennPlot)
dev.off()

# SFARI
sfari = read.table(paste(analysis_dir, "SFARI.txt", sep = "/"), header = TRUE, fill = TRUE, sep = "\t")
head(sfari)

# Venn diagram
# Extended Data Figure 9f
pdf(paste(plot_dir, "VennDiagramSFARI.pdf", sep = "/"), height = 4, width = 4.5)
vennPlot = venn.diagram(list("Day 150" = 1:length(d150_sig), "SFARI" = (length(d150_sig)-length(intersect(d150_sig, sfari$gene.symbol))+1):((length(d150_sig) + length(sfari$gene.symbol))-length(intersect(d150_sig, sfari$gene.symbol)))),fill = c("darkslateblue", "peru"), fontfamily = "Helvetica", cex=1.5, cat.fontfamily = "Helvetica", cat.cex = 1.5, filename=NULL, inverted = TRUE, cat.pos = c(150, 220), cat.dist = 0.05)
grid.newpage()
grid.draw(vennPlot)
dev.off()

# Extended Data Figure 9g (second plot in the series was used)
pdf(paste(plot_dir, "SFARI_histograms.pdf", sep = "/"), height = 1.7, width = 3.4)

sfari_overlap_all = intersect(d150_sig, sfari$gene.symbol)
dataHC150sig = hyCS_deseq_HC$D150_HumanLFC[which(rownames(hyCS_deseq_HC) %in% d150_sig)]
p1 <- qplot(dataHC150sig, 
            geom = "histogram", 
            main="Background", 
            xlab="Log2(ASE Ratio)",
            ylab="Frequency",
            fill=I("white"), 
            col=I("black"), 
            binwidth=0.75) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0, linetype = "longdash")

p2 <- qplot(dataHC150sig[which(d150_sig %in% sfari_overlap_all)], 
            geom = "histogram", 
            main="All SFARI genes", 
            xlab="Log2(ASE Ratio)",
            ylab="Frequency",
            fill=I("grey80"), 
            col=I("black"), 
            binwidth=0.2) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0, linetype = "longdash")
grid.arrange(p1,p2, ncol=2)

sfari_overlap_4321 = intersect(d150_sig, sfari$gene.symbol[which(sfari$gene.score<=4)])
length(sfari_overlap_4321)
p3 <- qplot(dataHC150sig[which(d150_sig %in% sfari_overlap_4321)], 
            geom = "histogram", 
            main="SFARI genes (score <=4)", 
            xlab="Log2(ASE Ratio)",
            ylab="Frequency",
            fill=I("grey80"), 
            col=I("black"), 
            binwidth=0.2) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0, linetype = "longdash")
grid.arrange(p1,p3, ncol=2)

sfari_overlap_321 = intersect(d150_sig, sfari$gene.symbol[which(sfari$gene.score<=3)])
length(sfari_overlap_321)
p4 <- qplot(dataHC150sig[which(d150_sig %in% sfari_overlap_321)], 
            geom = "histogram", 
            main="SFARI genes (score <=3)", 
            xlab="Log2(ASE Ratio)",
            ylab="Frequency",
            fill=I("grey80"), 
            col=I("black"), 
            binwidth=0.2) +
  theme_bw() + 
  ylim(0,6) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0, linetype = "longdash")
grid.arrange(p1,p4, ncol=2)

sfari_overlap_21 = intersect(d150_sig, sfari$gene.symbol[which(sfari$gene.score<=2)])
length(sfari_overlap_21)
p5 <- qplot(dataHC150sig[which(d150_sig %in% sfari_overlap_21)], 
            geom = "histogram", 
            main="SFARI genes (score <=2)", 
            xlab="Log2(ASE Ratio)",
            ylab="Frequency",
            fill=I("grey80"), 
            col=I("black"), 
            binwidth=0.2) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0, linetype = "longdash")
grid.arrange(p1,p5, ncol=2)
dev.off()

# Box plots of allelic expression in hyCS (for top genes)
hyCS_hc_tpm = read.table(paste(analysis_dir, "hyCS_allelic_TPM.rev.txt", sep = "/"), header = TRUE, row.names = 1)
head(hyCS_hc_tpm)
dim(hyCS_hc_tpm)

hyCS_meta = read.table(paste(data_dir, "hyCS.meta.txt", sep = "/"), header = TRUE)
hyCS_meta$Line = as.character(hyCS_meta$Line)
hyCS_meta$Induction = as.character(hyCS_meta$Induction)


boxwidth = 0.7
box <-function(gene){
  timepoints = rep(hyCS_meta$Age, 2)
  species = c(rep("Human", 25), rep("Chimp", 25))
  temp = data.frame(t(as.matrix(rbind(hyCS_hc_tpm[which(row.names(hyCS_hc_tpm) == gene),], timepoints, species))))
  names(temp) = c("exp", "timepoints", "species")
  temp$timepoints = factor(temp$timepoints, levels = c(50, 100, 150, 200),ordered = TRUE)
  temp = temp[order(temp$timepoints),]
  names(temp) = c("exp", "Timepoint", "Species")
  p <- ggplot() + 
    geom_boxplot(data = temp, aes(x=factor(Timepoint), y=as.numeric(as.character(exp)), fill=Species), outlier.shape=NA, position = "dodge", width = boxwidth, colour = "black", size = 0.5) +
    ggtitle(paste("Allelic expression:", gene)) + 
    theme_bw() + 
    theme(legend.position = "none") +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    geom_smooth(data = temp, aes(x=factor(Timepoint), y=as.numeric(as.character(exp)), group=Species, fill=Species, color=Species), size = 0.4, linetype = "longdash", se=FALSE) +
    #geom_smooth(data = temp[which(temp$Species == "Human"),], aes(x = as.numeric(factor(Timepoint)), y = as.numeric(as.character(exp))), formula = y~x,colour = "indianred3", size = 0.5, linetype = "longdash", se = FALSE) +
    geom_point(data = temp, aes(x=factor(Timepoint), y=as.numeric(as.character(exp)), fill = factor(Species)), position = position_dodge(width = boxwidth), col = "black", size = pointsize) +
    scale_color_manual(values=c("dodgerblue4", "tomato3"))+
    scale_fill_manual(values=c("dodgerblue3", "salmon"))+
    theme(plot.title = element_text(size=fontsize))+
    xlab("Timepoint") + 
    ylab("Expression (TPM)") +
    theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))
  p  
}

# Extended Data Figure 9h-i
pdf(paste(plot_dir, "SFARI_topGenes.rev.pdf", sep = "/"), height = 1.6, width = 1.6, useDingbats = FALSE)
fontsize = 8
pointsize = 0.2
box("GRIN2A")
box("SCN1A")
dev.off()

# Extended Data Figure 8g
pdf(paste(plot_dir, "TOPGenes_allelicBox.rev.pdf", sep = "/"), height = 1.5, width = 1.7, useDingbats = FALSE)
fontsize = 8
pointsize = 0.2
box("SSTR2")
box("PMP2")
dev.off()


# Top genes
dim(hyCS_deseq_HC) # 12912
length(d150_sig) # 2088
length(which(d150_sig %in% neural_genes)) # 183

top_genes = d150_sig[which(d150_sig %in% neural_genes)]

dim(hcCS_d150_deseq_keep)

filt_df = as.data.frame(cbind(dataHC150sig[which(d150_sig %in% top_genes)], hcCS_d150_deseq_keep$Day150_LFC[which(rownames(hcCS_d150_deseq_keep) %in% top_genes)]))
names(filt_df) = c("Hybrid", "Parent")
rownames(filt_df) = top_genes

# Extended Data Figure 9k
pdf(paste(plot_dir, "topGenesASE.D150.rev.pdf", sep = "/"), height = 2, width = 2, useDingbats = FALSE)
fontsize = 8
pointsize = 0.5
ggplot(data = filt_df, aes(x = Hybrid, y = Parent))+
  geom_point(col = ifelse(rownames(filt_df) %in% c("PMP2", "SSTR2") , "red", "black"),size = pointsize)+
  geom_abline(slope = 1, intercept = 0, linetype ="longdash") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  theme_bw() +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))+ theme(legend.title = element_blank()) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

dev.off()  

cor(filt_df$Hybrid, filt_df$Parent)*cor(filt_df$Hybrid, filt_df$Parent) # 0.43

write.table(file = paste(analysis_dir, "d150_topGenes.rev.txt", sep = "/"), x = filt_df, col.names = T, row.names = T, quote = F, sep = "\t")


# Directionality
length(which(filt_df$Hybrid<0 & filt_df$Parent>0)) # 24
length(which(filt_df$Hybrid>0 & filt_df$Parent<0)) # 16

(24+16)/183 # 22% (so 78% agree)

# PMP2 and SSTR2
# Create box plots of ASE over time
ase = read.table(paste(analysis_dir, "hyCS_ase_perSample.rev.txt", sep = "/"), header = TRUE, row.names = 1)

ase_box <-function(gene){
  temp = data.frame(t(as.matrix(rbind(ase[which(row.names(ase) == gene),], hyCS_meta$Age))))
  names(temp) = c("exp", "Timepoint")
  temp$Timepoint = factor(temp$Timepoint, levels = c(50, 100, 150, 200), ordered = TRUE)
  temp = temp[order(temp$Timepoint),]
  p <- ggplot() + 
    geom_boxplot(data = temp, aes(x=factor(Timepoint), y=as.numeric(as.character(exp)), fill=Timepoint), outlier.shape=NA) + 
    ggtitle(paste("ASE:", gene)) + theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values=c("thistle1", "plum2", "plum3", "plum4"))+
    geom_point(data = temp, aes(x=factor(Timepoint), y=as.numeric(as.character(exp))), position = position_dodge(width = boxwidth), col = "black", size = 0.5) +
    #geom_smooth(data = temp, aes(x=factor(Timepoint), y=as.numeric(as.character(exp)), group = 1), color="black", se=FALSE, size = 0.4, linetype = "longdash") +
    theme(plot.title = element_text(size=fontsize))+
    xlab("Timepoint") + ylab("Log2(ASE)")+
    theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))+
    geom_hline(yintercept = 0, linetype= "longdash") +
    theme(legend.position = "none")
  p
}

# Figure 4a
pdf(paste(plot_dir, "hyCS_ASEboxPlot.rev.pdf", sep = "/"), height = 1.5, width = 1.5, useDingbats = FALSE)
ase_box("PMP2")
ase_box("SSTR2")
dev.off()

# Parental expression as a bar plot
hcCS_tot_cpm = read.table(file = paste(analysis_dir, "hcCS_tot_cpm.rev.txt", sep = "/"), header = TRUE, row.names = 1, sep = "\t")

hCS_200_meta = read.table(paste(data_dir, "hCS_200.meta.txt", sep = "/"), header = TRUE)
hCS_35_meta = read.table(paste(data_dir, "hCS_35.meta.txt", sep = "/"), header = TRUE)
hCS_meta = as.data.frame(rbind(hCS_35_meta, hCS_200_meta))
hCS_meta$Species = as.character(hCS_meta$Species)
hCS_meta$Line = as.character(hCS_meta$Line)

cCS_200_meta = read.table(paste(data_dir, "cCS_200.meta.txt", sep = "/"), header = TRUE)
cCS_35_meta = read.table(paste(data_dir, "cCS_35.meta.txt", sep = "/"), header = TRUE)
cCS_meta = as.data.frame(rbind(cCS_35_meta, cCS_200_meta))
cCS_meta$Species = as.character(cCS_meta$Species)
cCS_meta$Line = as.character(cCS_meta$Line)

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
    geom_point(data = temp, aes(x = factor(timepoint), y = exp, fill = species), size = pointsize, stat = "identity", position_dodge(width = 0.7)) +
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
 
# NOT USED
pdf(paste(plot_dir, "ParentalExpOverTime.rev.pdf", sep = "/"), height = 1.5, width = 3, useDingbats = FALSE)
pointsize = 0.2
hc_bar("SSTR2")
hc_bar("PMP2")
dev.off()

# As a boxplot
hc_box <- function(gene){
  temp = data.frame(as.matrix(t(rbind(hcCS_tot_cpm[which(rownames(hcCS_tot_cpm)==gene),], hcCS_species, hcCS_time))))
  names(temp) = c("exp", "species", "timepoint")
  temp$exp = as.numeric(as.character(temp$exp))
  temp$timepoint = factor(temp$timepoint, levels = c(15,25,35,50,100, 150, 200), ordered = TRUE)
  temp = temp[order(temp$timepoint),]
  #temp2 = temp %>% group_by(species, timepoint) %>% summarise(mean_exp = mean(exp), se_exp = sd(exp)/sqrt(n()))
  #ymax = max(temp$exp) + max(temp2$se_exp)
  p = ggplot() + 
    geom_boxplot(data = temp, aes(x=factor(timepoint), y=exp, fill = species), outlier.shape=NA, position = "dodge", width = boxwidth, colour = "black", size = 0.5) +
    geom_point(data = temp, aes(x = factor(timepoint), y = exp, fill = species), size = pointsize,  position = position_dodge(width = boxwidth), col = "black") +
    geom_smooth(data = temp, aes(x=factor(timepoint), y=exp, group=species, fill=species, color=species), size = 0.4, linetype = "longdash", se=FALSE) +
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("dodgerblue3", "salmon")) +
    scale_color_manual(values = c("dodgerblue3", "salmon")) +
    ggtitle(gene)+
    ylab("Expression (TPM)") +
    xlab("Timepoint") +
    theme(plot.title = element_text(size = fontsize)) +
    theme(legend.position = "none") +
   # ylim(0, ymax) +
    theme(axis.text=element_text(size=fontsize), axis.title = element_text(size = fontsize))
  p
}
boxwidth = 0.7
# Extended Data Figures 8f and 10a
pdf(paste(plot_dir, "ParentalExpOverTimeBox.rev.pdf", sep = "/"), height = 1.5, width = 3, useDingbats = FALSE)
pointsize = 0.2
hc_box("SSTR2")
hc_box("PMP2")
dev.off()


geom_boxplot(data = temp, aes(x=factor(Timepoint), y=as.numeric(as.character(exp)), fill=Species), outlier.shape=NA, position = "dodge", width = boxwidth, colour = "black", size = 0.5) +
  ggtitle(paste("Allelic expression:", gene)) + 
  theme_bw() + 
  theme(legend.position = "none") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  geom_smooth(data = temp, aes(x=factor(Timepoint), y=as.numeric(as.character(exp)), group=Species, fill=Species, color=Species), size = 0.4, linetype = "longdash", se=FALSE) +
  #geom_smooth(data = temp[which(temp$Species == "Human"),], aes(x = as.numeric(factor(Timepoint)), y = as.numeric(as.character(exp))), formula = y~x,colour = "indianred3", size = 0.5, linetype = "longdash", se = FALSE) +
  geom_point(data = temp, aes(x=factor(Timepoint), y=as.numeric(as.character(exp)), fill = factor(Species)), position = position_dodge(width = boxwidth), col = "black", size = pointsize) +
  


# Single cell plots
n_coords = read.table(paste(analysis_dir, "neural_sc_umap_coords.rev.txt", sep = "/"), header = TRUE, row.names = 1)
sc_keep = sc[geneKeep,]
dim(sc_keep)
sc_cpm = sweep(as.matrix(sc_keep), 2, as.double(colSums(sc_keep)/1000000), `/`)
sc_fpkm = sweep(as.matrix(sc_cpm), 1, as.double(geneInfo$genesizeKb[geneKeepIdx]), `/`)
dim(sc_fpkm)

sc_neural = sc_fpkm[,which(paste("Rachels", colnames(sc_fpkm), sep="") %in% rownames(neural_clusters))]
dim(sc_neural)

fontsize = 8
pointsize =0.5
genePlotNeural <- function(gene,exp = sc_neural){
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

# Figure 4b
pdf(paste(plot_dir, "neuralGenePlotsTOP.rev.pdf", sep = "/"), height = 1, 1.5, useDingbats = FALSE) 
pointsize = 0.3
genePlotNeural("SSTR2")
genePlotNeural("PMP2")
dev.off()

# Extended Data Figure 8e
pdf(paste(plot_dir, "GenePlotsTOP.rev.pdf", sep = "/"), height = 2, 2, useDingbats = FALSE) 
pointsize = 0.2
genePlot("SSTR2")
genePlot("PMP2")
dev.off()


# Layers data
# Also look at these genes in the in vivo data 
kh_hum = read.table(paste(analysis_dir, "Khaitovich_Layers_Data/DS1_sectionRPKM_human.tsv", sep = "/"), header = TRUE)
kh_chi = read.table(paste(analysis_dir, "Khaitovich_Layers_Data/DS1_sectionRPKM_chimp.tsv", sep = "/"), header = TRUE)
kh_mac = read.table(paste(analysis_dir, "Khaitovich_Layers_Data/DS1_sectionRPKM_macaque.tsv", sep = "/"), header = TRUE)

dim(kh_hum)
dim(kh_chi)
dim(kh_mac)

kh_plot <- function(gene, title = gene){
  index = which(rownames(kh_hum) == gene)
  temp = as.data.frame(t(rbind(kh_hum[index,],kh_chi[index,],kh_mac[index,])))
  colnames(temp) = c("Human", "Chimp", "Macaque")
  temp$Section = as.character(seq(1,17))
  temp2 = melt(temp, id="Section")
  temp2$Species = c(rep("Human", 17), rep("Chimp", 17), rep("Macaque", 17))
  #temp2$Section = factor(temp2$Section, levels = temp2$Section, ordered = TRUE)
  p = ggplot(temp2, aes(x=as.factor(Section), y=value, color=Species, group=variable)) + geom_point() + geom_line() +
    scale_x_discrete(limits=c(seq(1,17))) +
    geom_point(size = 0.5) +
    theme(plot.title = element_text(size=fontsize))+
    xlab("Cortical section") + 
    ylab("Expression (RPKM)") +
    scale_color_manual(values = c("dodgerblue3", "salmon", "grey50")) +
    theme_bw() +
    theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    ggtitle(paste("Khaitovich 2017:", title)) +
    geom_vline(xintercept = c(1.5,3.5, 7,9, 12.5, 16.5), linetype="longdash") +
    theme(legend.position = "none")
  p
}

# Extended Data Figure 10b
pdf(paste(plot_dir, "topGenesLayers.pdf", sep = "/"), height = 1.5, width = 2, useDingbats = FALSE)
fontsize = 8
p = kh_plot("ENSG00000180616.4", "SSTR2")
p
dev.off()

# FURA2  
fura = read.table(paste(analysis_dir, "Calcium_Jan2020/Amplitudes_all_shifted.txt", sep = "/"), header = FALSE, fill = TRUE)
head(fura)
tail(fura)

group_by(fura, V2) %>%
  summarise(
    count = n(),
    mean = mean(V1, na.rm = TRUE),
    sd = sd(V1, na.rm = TRUE)
  )

fura$V2 = factor(fura$V2, levels = c("Chimp_base", "Chimp_drug", "Human_base", "Human_drug"), ordered = TRUE)
fura = fura[order(fura$V2),]

# Figure 4i
pdf(paste(plot_dir, "fura_amp.pdf", sep = "/"), height = 1.9, width = 1.7, useDingbats = FALSE)
fontsize = 8
ggplot(fura, aes(x = V2, y = V1, fill = V2)) +
  #geom_violin(color = "black", trim = FALSE)+
 # geom_point(position = "jitter", size = 0.1 )+
  geom_dotplot( binaxis='y', stackdir='center',  dotsize=1.5, binwidth = 0.015, position = "dodge", stroke=0)+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
               geom = "crossbar", width = 0.5) +
  theme_bw() +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Condition") + 
  ylab("Amplitude") +
  scale_fill_manual(values = c("dodgerblue2", "cyan1", "salmon","hotpink1")) +
  theme(legend.position = "none") 

dev.off()

wilcox.test(fura$V1[which(fura$V2=="Human_base")], fura$V1[which(fura$V2=="Human_drug")]) # 5.082e-05
wilcox.test(fura$V1[which(fura$V2=="Human_base")], fura$V1[which(fura$V2=="Chimp_base")]) # 0.01959
wilcox.test(fura$V1[which(fura$V2=="Chimp_base")], fura$V1[which(fura$V2=="Chimp_drug")]) # 0.08708
wilcox.test(fura$V1[which(fura$V2=="Chimp_drug")], fura$V1[which(fura$V2=="Human_drug")]) # 4.143e-14

mean(fura$V1[which(fura$V2=="Human_base")]) # 0.4568498
mean(fura$V1[which(fura$V2=="Human_drug")]) # 0.5628883
mean(fura$V1[which(fura$V2=="Chimp_base")]) # 0.3884296
mean(fura$V1[which(fura$V2=="Chimp_drug")]) # 0.3372768

log2(mean(fura$V1[which(fura$V2=="Human_drug")])/mean(fura$V1[which(fura$V2=="Human_base")])) # 0.3011286
log2(mean(fura$V1[which(fura$V2=="Chimp_drug")])/mean(fura$V1[which(fura$V2=="Chimp_base")])) # -0.2037199
log2(mean(fura$V1[which(fura$V2=="Human_base")])/mean(fura$V1[which(fura$V2=="Chimp_base")])) # 0.234067
log2(mean(fura$V1[which(fura$V2=="Human_drug")])/mean(fura$V1[which(fura$V2=="Chimp_drug")])) # 0.7389155

## NOT USED
# FURA data split by line (2 human, 1 chimp) ##TODO
furaLine = read.table(paste(analysis_dir, "Calcium_Jan2020/fura_all_byLine_long.txt", sep = "/"), header = FALSE, fill = TRUE)
head(furaLine)
tail(furaLine)

group_by(furaLine, V2) %>%
  summarise(
    count = n(),
    mean = mean(V1, na.rm = TRUE),
    sd = sd(V1, na.rm = TRUE)
  )

# NOT USED
pdf(paste(plot_dir, "fura_byLine.pdf", sep = "/"), height = 2, width = 3, useDingbats = FALSE)
fontsize = 8
ggplot(furaLine, aes(x = V2, y = V1, fill = V2)) +
  #geom_violin(color = "black")+
  geom_dotplot( binaxis='y', stackdir='center',  dotsize=4, binwidth = 0.2, position = "dodge", stroke=0)+
  theme_bw() +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize), axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Condition") + 
  ylab("# of spikes") +
  scale_fill_manual(values = c("dodgerblue2", "cyan", "salmon","hotpink1","salmon","hotpink1")) +
  theme(legend.position = "none") 

dev.off()

wilcox.test(furaLine$V1[which(furaLine$V2=="C3649_base")], furaLine$V1[which(furaLine$V2=="C3649_cyn")]) # 
wilcox.test(furaLine$V1[which(furaLine$V2=="H20682_base")], furaLine$V1[which(furaLine$V2=="H20682_cyn")]) #
wilcox.test(furaLine$V1[which(furaLine$V2=="H20961_base")], furaLine$V1[which(furaLine$V2=="H20961_cyn")]) # 



# Frequency data
freq = read.table(paste(analysis_dir, "Calcium_Jan2020/GluFRAP_all_NoSpikes_Filter1.1uncaging_2ormorespikes.txt", sep = "/"), header = FALSE, fill = TRUE)
head(freq)
tail(freq)

group_by(freq, V2) %>%
  summarise(
    count = n(),
    mean = mean(V1, na.rm = TRUE),
    sd = sd(V1, na.rm = TRUE)
  )


# Figure 4k
pdf(paste(plot_dir, "freq.pdf", sep = "/"), height = 1.7, width = 1.7, useDingbats = FALSE)
fontsize = 8
ggplot(freq, aes(x = V2, y = V1, fill = V2)) +
  geom_violin(color = "black", trim = FALSE, size = 0.2, draw_quantiles =0.5)+
  geom_jitter(shape=16, size = 0.2,  width = 0.2)+
 # geom_point(position = "jitter", size = 0.1 )+
  #geom_dotplot( binaxis='y', stackdir='centerwhole',  dotsize=1, binwidth = 0.2, position = "jitter", stroke=0, method = "histodot")+
  theme_bw() +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Condition") + 
  ylab("Num. of Ca spikes") +
  scale_fill_manual(values = c("dodgerblue2", "cyan1", "salmon","hotpink1")) +
  theme(legend.position = "none") 

dev.off()

wilcox.test(freq$V1[which(freq$V2=="Chimp_base")], freq$V1[which(freq$V2=="Chimp_cyn")]) # 0.1139
wilcox.test(freq$V1[which(freq$V2=="Chimp_base")], freq$V1[which(freq$V2=="Human_base")]) # 0.2715
wilcox.test(freq$V1[which(freq$V2=="Human_base")], freq$V1[which(freq$V2=="Human_cyn")]) # 0.001338
wilcox.test(freq$V1[which(freq$V2=="Human_cyn")], freq$V1[which(freq$V2=="Chimp_cyn")]) # 0.001181


# Split by cell line
freqLine = read.table(paste(analysis_dir, "Calcium_Jan2020/GluFRAP_all_byline_NoSpikes_Filter1.1uncaging_2ormorespikes.txt", sep = "/"), header = FALSE, fill = TRUE)
head(freqLine)
tail(freqLine)

group_by(freqLine, V2) %>%
  summarise(
    count = n(),
    mean = mean(V1, na.rm = TRUE),
    sd = sd(V1, na.rm = TRUE)
  )

# Figure 4l
pdf(paste(plot_dir, "freq_active_byLine.pdf", sep = "/"), height = 1.7, width = 4, useDingbats = FALSE)
fontsize = 8
ggplot(freqLine, aes(x = V2, y = V1, fill = V2)) +
  geom_violin(color = "black", trim = FALSE, size = 0.2, draw_quantiles =0.5)+
  geom_jitter(shape=16, size = 0.2,  width = 0.2)+
  #geom_dotplot( binaxis='y', stackdir='center',  dotsize=4, binwidth = 0.2, position = "dodge", stroke=0)+
  theme_bw() +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Condition") + 
  ylab("Num of Ca spikes") +
  scale_fill_manual(values = c("dodgerblue2", "cyan","dodgerblue2", "cyan", "dodgerblue2", "cyan","salmon","hotpink1","salmon","hotpink1","salmon","hotpink1")) +
  scale_color_manual(values = c("dodgerblue2", "cyan","dodgerblue2", "cyan", "dodgerblue2", "cyan","salmon","hotpink1","salmon","hotpink1","salmon","hotpink1")) +
  theme(legend.position = "none") 

dev.off()

wilcox.test(freqLine$V1[which(freqLine$V2=="C3647-base")], freqLine$V1[which(freqLine$V2=="C3647-cyn")]) # NA
wilcox.test(freqLine$V1[which(freqLine$V2=="C3649-base")], freqLine$V1[which(freqLine$V2=="C3649-cyn")]) # 0.3154
wilcox.test(freqLine$V1[which(freqLine$V2=="C3651-base")], freqLine$V1[which(freqLine$V2=="C3651-cyn")]) # 0.1444
wilcox.test(freqLine$V1[which(freqLine$V2=="H20682_base")], freqLine$V1[which(freqLine$V2=="H20682_cyn")]) # 0.05019
wilcox.test(freqLine$V1[which(freqLine$V2=="H21796_base")], freqLine$V1[which(freqLine$V2=="H21796_cyn")]) # 0.03051
wilcox.test(freqLine$V1[which(freqLine$V2=="H20961_base")], freqLine$V1[which(freqLine$V2=="H20961_cyn")]) # 0.07824

# Combined p-values (using Fisher's method)
# Human
-2*(log(0.07824) + log(0.05019) + log(0.03051)) # statistic = 18.05923, dof = 6
pchisq(18.05923, 6, lower.tail = FALSE) # 0.006085869

# Chimpanzee
-2*(log(0.3154) + log(0.1444)) # statistic = 6.178163, dof = 4, p-value = 0.1
pchisq(6.178163, 4, lower.tail = FALSE) # 0.1862322
    
    
## SSTR2 bar plots
tubb3 = read.table(paste(analysis_dir, "SSTR2_TUBB3.txt",sep = "/"), header = FALSE)
tubb3

# Figure 4f
pdf(paste(plot_dir, "SSTR2_TUBB3.pdf", sep = "/"), height = 1.2, width = 1.35, useDingbats = FALSE)
pointsize = 0.01
fontsize = 8
boxwidth = 0.55
p = ggplot(data = tubb3, aes(x = factor(V2), y = V1)) +
  geom_boxplot( position = "dodge", width = boxwidth, aes(fill = V2), color = "black", outlier.shape = NA) +
  geom_point(data = tubb3, aes(x = factor(V2), y = V1), color = "black",size = pointsize, stat = "identity", position = position_jitter(width = 0.1)) +
  scale_fill_manual(values = c("dodgerblue2", "salmon")) +
  scale_color_manual(values = c("dodgerblue4", "indianred3")) +
  ylab("SSTR2/TUBB3")+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text = element_text(size = fontsize, color = "black"), axis.title = element_text(size = fontsize), title = element_text(size = fontsize)) +
  theme(axis.title.x = element_blank())+
  theme(axis.ticks = element_line(color = "black"))+
  theme(legend.position = "none")
p
dev.off()


tubb3 = read.table(paste(analysis_dir, "TUBB3.txt",sep = "/"), header = FALSE)

# Figure 4e
pdf(paste(plot_dir, "TUBB3.pdf", sep = "/"), height = 1.2, width = 1.35, useDingbats = FALSE)
pointsize = 0.01
fontsize = 8
boxwidth = 0.55
p = ggplot(data = tubb3, aes(x = factor(V2), y = V1)) +
  geom_boxplot( position = "dodge", width = boxwidth, aes(fill = V2), color = "black", outlier.shape = NA) +
  geom_point(data = tubb3, aes(x = factor(V2), y = V1), color = "black",size = pointsize, stat = "identity", position = position_jitter(width = 0.1)) +
  scale_fill_manual(values = c("dodgerblue2", "salmon")) +
  scale_color_manual(values = c("dodgerblue4", "indianred3")) +
  ylab("TUBB3")+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text = element_text(size = fontsize, color = "black"), axis.title = element_text(size = fontsize), title = element_text(size = fontsize)) +
  theme(axis.title.x = element_blank())+
  theme(axis.ticks = element_line(color = "black"))+
  theme(legend.position = "none")
p
dev.off()




map2 = read.table(paste(analysis_dir, "MAP2.txt",sep = "/"), header = FALSE)

# Extended Data Figure 10d
pdf(paste(plot_dir, "MAP2.pdf", sep = "/"), height = 1.2, width = 1.35, useDingbats = FALSE)
pointsize = 0.01
fontsize = 8
boxwidth = 0.55
p = ggplot(data = map2, aes(x = factor(V2), y = V1)) +
  geom_boxplot( position = "dodge", width = boxwidth, aes(fill = V2), color = "black", outlier.shape = NA) +
  geom_point(data = map2, aes(x = factor(V2), y = V1), color = "black",size = pointsize, stat = "identity", position = position_jitter(width = 0.1)) +
  scale_fill_manual(values = c("dodgerblue2", "salmon")) +
  scale_color_manual(values = c("dodgerblue4", "indianred3")) +
  ylab("MAP2")+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text = element_text(size = fontsize, color = "black"), axis.title = element_text(size = fontsize), title = element_text(size = fontsize)) +
  theme(axis.title.x = element_blank())+
  theme(axis.ticks = element_line(color = "black"))+
  theme(legend.position = "none")
p
dev.off()

map2 = read.table(paste(analysis_dir, "MAP2_SSTR2.txt",sep = "/"), header = FALSE)

# Extended Data Figure 10e
pdf(paste(plot_dir, "MAP2_SSTR2.pdf", sep = "/"), height = 1.2, width = 1.35, useDingbats = FALSE)
pointsize = 0.01
fontsize = 8
boxwidth = 0.55
p = ggplot(data = map2, aes(x = factor(V2), y = V1)) +
  geom_boxplot( position = "dodge", width = boxwidth, aes(fill = V2), color = "black", outlier.shape = NA) +
  geom_point(data = map2, aes(x = factor(V2), y = V1), color = "black",size = pointsize, stat = "identity", position = position_jitter(width = 0.1)) +
  scale_fill_manual(values = c("dodgerblue2", "salmon")) +
  scale_color_manual(values = c("dodgerblue4", "indianred3")) +
  ylab("MAP2")+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text = element_text(size = fontsize, color = "black"), axis.title = element_text(size = fontsize), title = element_text(size = fontsize)) +
  theme(axis.title.x = element_blank())+
  theme(axis.ticks = element_line(color = "black"))+
  theme(legend.position = "none")
p
dev.off()





































## NOT USED
# Plot drug effect versus mean expression at D150
hcCS_tot_cpm["SSTR2", which(hcCS_time == 150)]
x = c(3.695,8.913, 77.54, 86.27, 109.98, 70.73)

log10(mean(freqLine$V1[which(freqLine$V2=="C3647-cyn")])/mean(freqLine$V1[which(freqLine$V2=="C3647-base")]))
log10(mean(freqLine$V1[which(freqLine$V2=="C3649-cyn")])/mean(freqLine$V1[which(freqLine$V2=="C3649-base")]))
log10(mean(freqLine$V1[which(freqLine$V2=="C3651-cyn")])/mean(freqLine$V1[which(freqLine$V2=="C3651-base")]))
log10(mean(freqLine$V1[which(freqLine$V2=="H20682_cyn")])/mean(freqLine$V1[which(freqLine$V2=="H20682_base")]))
log10(mean(freqLine$V1[which(freqLine$V2=="H21796_cyn")])/mean(freqLine$V1[which(freqLine$V2=="H21796_base")]))
log10(mean(freqLine$V1[which(freqLine$V2=="H20961_cyn")])/mean(freqLine$V1[which(freqLine$V2=="H20961_base")]))
y = c(0, 0.0099, -0.1187806, 0.08889948,0.1731907,0.1440872)

spec = c("Chimp", "Chimp", "Chimp", "Human", "Human", "Human")

pdf(paste(plot_dir, "freq_active_byLine_exp.pdf", sep = "/"), height = 1.7, width = 1.7, useDingbats = FALSE)
fontsize = 8
pointsize = 1
qplot(x,y, size = I(pointsize))+
  geom_point(aes(col = spec), size = pointsize) +
  theme_bw() +
  ylab("Log10(Drug/Base)") +
  xlab("SSTR2 Expression (CPM)") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("dodgerblue2", "salmon")) +
  theme(axis.text=element_text(size=fontsize), axis.title=element_text(size=fontsize), axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

cor(x, y)*cor(x, y) # 0.77



# marker genes
hcCS_line = c(hCS_meta$Line, cCS_meta$Line)
hc_bar_line<- function(gene){
  temp = data.frame(as.matrix(t(rbind(hcCS_tot_cpm[which(rownames(hcCS_tot_cpm)==gene),], hcCS_species, hcCS_time, hcCS_line))))
  names(temp) = c("exp", "species", "timepoint", "line")
  temp$exp = as.numeric(as.character(temp$exp))
  temp$timepoint = factor(temp$timepoint, levels = c(15,25,35,50,100, 150, 200), ordered = TRUE)
  temp = temp[order(temp$timepoint),]
  temp2 = temp %>% group_by(species, timepoint) %>% summarise(mean_exp = mean(exp), se_exp = sd(exp)/sqrt(n()))
  ymax = max(temp$exp) + max(temp2$se_exp)
  p = ggplot(temp2, aes(x=factor(timepoint), y=mean_exp, fill = species)) + 
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_errorbar(aes(ymin = mean_exp - se_exp, ymax = mean_exp + se_exp), width=0.2, position = position_dodge(0.7))+
    geom_point(data = temp, aes(x = factor(timepoint), y = exp, fill = species, col = factor(temp$line)), size = 1, stat = "identity", position_dodge(width = 0.7)) +
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    scale_fill_manual(values = c("dodgerblue3", "salmon")) +
    scale_color_manual(values = c("black", "black", "black", "darkred", "red", "hotpink"))+
    ggtitle(gene)+
    theme(plot.title = element_text(size = fontsize)) +
    theme(legend.position = "none") +
    ylim(0, ymax) +
    theme(axis.text=element_text(size=fontsize), axis.text.x=element_text(angle = 45, hjust =1), axis.title=element_blank())
  p
}


# Dorsal Forebrain
hc_bar_line("FOXG1")
hc_bar_line("EMX1")
hc_bar_line("NEUROD6")

# Ventral Forebrain
hc_bar_line("DLX1")
hc_bar_line("DLX2")
hc_bar_line("NKX2-1")
hc_bar_line("NKX6-2")

# Hindbrain/Spinal Cord
hc_bar_line("HOXB4")
hc_bar_line("GATA2")
hc_bar_line("GATA3")
hc_bar_line("HOXA2")
hc_bar_line("EN1")

# Meso/Mesenchyme
hc_bar_line("TBXT")
hc_bar_line("POSTN")
hc_bar_line("TBX6")
hc_bar_line("MYOD")

# Neural Crest
hc_bar_line("FOXD3")
hc_bar_line("HNK1")
hc_bar_line("ENDO")
hc_bar_line("SOX17")
hc_bar_line("GATA4")
hc_bar_line("GATA6")

# Maturation
hc_bar_line("DCX")
hc_bar_line("STMN2")
hc_bar_line("GFAP")
hc_bar_line("SOX2")
hc_bar_line("EOMES")
hc_bar_line("PAX6")
hc_bar_line("BCL11B")


