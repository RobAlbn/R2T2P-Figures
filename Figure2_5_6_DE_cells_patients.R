library(cowplot)
library(zoo)
library(reshape2)
library(limma)
library(UpSetR)
library(readxl)
library(ggpubr)
library(GenomeInfoDb)

suppressPackageStartupMessages(library("RiboseQC"))
suppressPackageStartupMessages(library("ggfortify"))
trans <- function(x) {reso=-log(x, 10);reso[is.na(reso)]=0;return(reso)}
inv <- function(x) 10^(-x)
setwd("/group/calviello/Roberto/R2T2P_scripts_for_Figures_with_container/R2T2P_Figures_GitHub_repo")

source("data/new_riboseqc_all_HT_oct23_26.R")
#orfs<-get(load("/scratch/lorenzo.calviello/lopit_R2T2P/ORFquant/combined_forORFquant.RData_final_ORFquant_results"))
load("data/gencode.v47.annotation.gtf.gz_Rannot")

multiDE_all1<-get(load("data/04h_vs_00h_multiDE_results_summary_allgenes.RData"))
multiDE_all2<-get(load("data/08h_vs_00h_multiDE_results_summary_allgenes.RData"))
multiDE_all3<-get(load("data/14h_vs_00h_multiDE_results_summary_allgenes.RData"))


target_genes_file <- read_excel("data/Table_S1_DUX4_targets.xls", skip=1, col_names=T)
target_genes <- unique(GTF_annotation$trann$gene_id[sapply(strsplit(GTF_annotation$trann$gene_id,"[.]"),"[[",1) %in% sapply(strsplit(target_genes_file$Ensembl,"[.]"),"[[",1)])
target_genes2 <- unique(GTF_annotation$trann$gene_id[match(target_genes_file$gene.name,GTF_annotation$trann$gene_name)])
target_genes<-unique(c(target_genes,target_genes2))

load("data/SRonly_UniqueORFq_Rannot")

multiDE_all<-rbind(multiDE_all1$uniq[[1]],
                   multiDE_all2$uniq[[1]],
                   multiDE_all3$uniq[[1]])
multiDE_all$gene_biotype<-as.character(multiDE_all$gene_biotype)
multiDE_all$gene_biotype[grepl(multiDE_all$gene_id,pattern="^R1|^R2")]<-multiDE_all$biotype[grepl(multiDE_all$gene_id,pattern="^R1|^R2")]
multiDE_all$gene_biotype[multiDE_all$gene_id%in%target_genes]<-"DUX4_targets"
multiDE_all$gene_biotype[multiDE_all$gene_biotype%in%c("lncRNA","pseudogene","ncRNA")]<-"non-coding"

#IMPORTANT: here filter for not overlapping DE genes
multiDE_all<-multiDE_all[!multiDE_all$gene_overlapping_noDE,]
multiDE_all<-multiDE_all[multiDE_all$biotype!="artifact",]
multiDE_all_idux4<-multiDE_all

supp2<-multiDE_all_idux4
supp2<-supp2[,!grepl(colnames(supp2),pattern="ovl|xy|len|Intron|dist|GCpc|n_exons")]
write.table(supp2, file="/group/calviello/Roberto/R2T2P_scripts_for_Figures_with_container/R2T2P_Figures_GitHub_repo/Supplementary_Tables/Supplementary_Table2.tsv",sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE)

df<-melt(table(multiDE_all$RNA_padj<.05,multiDE_all$experiment,multiDE_all$gene_biotype,multiDE_all$RNA_log2FC>0))
df<-df[df$Var1,]
df<-df[df$Var4,]
df$type<-"RNA-seq"
df1<-df
df<-melt(table(multiDE_all$Ribo_padj<.05,multiDE_all$experiment,multiDE_all$gene_biotype,multiDE_all$Ribo_log2FC>0))
df<-df[df$Var1,]
df<-df[df$Var4,]
df$type<-"Ribo-seq"
df2<-df
df<-rbind(df1,df2)
df$Var2<-gsub(as.character(df$Var2),pattern="_vs_00h",replacement="")
df$Var3<-as.character(df$Var3)
df$Var3[df$Var3=="novel_antisense_intr"]<-"novel_antisense_intronic"

df$Var3<-factor(df$Var3,levels=c("protein_coding","non-coding","DUX4_targets","novel_intergenic","novel_antisense","novel_antisense_intronic"))
df$type<-factor(as.character(df$type),levels=c("RNA-seq","Ribo-seq"))
colssi<-c("#0072B2", "#56B4E9", "#009E73", "grey30","grey", "goldenrod1")
names(colssi)<-c("novel_antisense", "novel_antisense_intronic", "novel_intergenic","protein_coding", "non-coding", "DUX4_targets")

gg_barDE <- ggplot(df, aes(x = Var2 , y = value,label=value,fill=Var3)) +
  #stat_summary(fun.data = mean_se,size=1.5) +
  #stat_summary(fun.y = mean, geom="point",colour="magenta", size=2) +
  scale_fill_manual(values = colssi,"") +
  geom_bar(stat="identity",position=position_dodge(width = .9)) +
  scale_y_log10()+
  ylab("# upregulated genes") +
  theme_bw()+
  xlab("")+
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=16), 
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16),
        panel.grid.major = element_blank(),  # remove grid lines from panel
        panel.grid.minor = element_blank(),  # remove grid lines from panel
        strip.text=element_text(size=16),legend.position="top",legend.text=element_text(size=14),
        legend.title=element_text(size=16)) +
  facet_wrap(~type,scales="fixed",nrow = 1) + geom_text(vjust=-.1, color="black",position=position_dodge(width = .9))
pdf(file = "Figures/SupplementaryFigure3.pdf",width = 11,height = 4)
gg_barDE+ coord_cartesian(ylim=c(1,5001))
dev.off()

wang_G1<-get(load("data/G1_FSHD_vs_G1_Control_multiDE_results_summary_allgenes.RData"))
wang_G2<-get(load("data/G2_FSHD_vs_G1_Control_multiDE_results_summary_allgenes.RData"))
wang_G3<-get(load("data/G3_FSHD_vs_G1_Control_multiDE_results_summary_allgenes.RData"))
wang_G4<-get(load("data/G4_FSHD_vs_G1_Control_multiDE_results_summary_allgenes.RData"))
yao_tubes<-get(load("data/FSHD_vs_Control_multiDE_results_summary_allgenes.RData"))


multiDE_all<-rbind(multiDE_all1$uniq[[1]],
                   multiDE_all2$uniq[[1]],
                   multiDE_all3$uniq[[1]])
multiDE_all$gene_biotype<-as.character(multiDE_all$gene_biotype)
multiDE_all$gene_biotype[grepl(multiDE_all$gene_id,pattern="^R1|^R2")]<-multiDE_all$biotype[grepl(multiDE_all$gene_id,pattern="^R1|^R2")]
multiDE_all$gene_biotype[multiDE_all$gene_id%in%target_genes]<-"DUX4_targets"
multiDE_all$gene_biotype[multiDE_all$gene_biotype%in%c("lncRNA","pseudogene","ncRNA")]<-"non-coding"

#IMPORTANT: here filter for not overlapping DE genes
multiDE_all<-multiDE_all[!multiDE_all$gene_overlapping_noDE,]
multiDE_all<-multiDE_all[multiDE_all$biotype!="artifact",]
multiDE_all_idux4<-multiDE_all

df<-melt(table(multiDE_all$RNA_padj<.05,multiDE_all$experiment,multiDE_all$gene_biotype,multiDE_all$RNA_log2FC>0))
df<-df[df$Var1,]
df<-df[df$Var4,]
df$type<-"RNA-seq"
df1<-df
df<-melt(table(multiDE_all$Ribo_padj<.05,multiDE_all$experiment,multiDE_all$gene_biotype,multiDE_all$Ribo_log2FC>0))
df<-df[df$Var1,]
df<-df[df$Var4,]
df$type<-"Ribo-seq"
df2<-df
df<-rbind(df1,df2)
df$Var2<-gsub(as.character(df$Var2),pattern="_vs_00h",replacement="")
df$Var3<-as.character(df$Var3)
df$Var3[df$Var3=="novel_antisense_intr"]<-"novel_antisense_intronic"

df$Var3<-factor(df$Var3,levels=c("protein_coding","non-coding","DUX4_targets","novel_intergenic","novel_antisense","novel_antisense_intronic"))
df$type<-factor(as.character(df$type),levels=c("RNA-seq","Ribo-seq"))
colssi<-c("#0072B2", "#56B4E9", "#009E73", "grey30","grey", "goldenrod1")
names(colssi)<-c("novel_antisense", "novel_antisense_intronic", "novel_intergenic","protein_coding", "non-coding", "DUX4_targets")

gg_barDE <- ggplot(df, aes(x = Var2 , y = value,label=value,fill=Var3)) +
  #stat_summary(fun.data = mean_se,size=1.5) +
  #stat_summary(fun.y = mean, geom="point",colour="magenta", size=2) +
  scale_fill_manual(values = colssi,"") +
  geom_bar(stat="identity",position=position_dodge(width = .9)) +
  scale_y_log10()+
  ylab("# upregulated genes") +
  theme_bw()+
  xlab("")+
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=16), 
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16),
        panel.grid.major = element_blank(),  # remove grid lines from panel
        panel.grid.minor = element_blank(),  # remove grid lines from panel
        strip.text=element_text(size=16),legend.position="top",legend.text=element_text(size=14),
        legend.title=element_text(size=16)) +
  facet_wrap(~type,scales="fixed",nrow = 1) + geom_text(vjust=-.1, color="black",position=position_dodge(width = .9))
# pdf(file = "Figures/SupplementaryFigure3A.pdf",width = 11,height = 4)
gg_barDE+ coord_cartesian(ylim=c(1,5001))
# dev.off()

df2<-df[!df$Var3%in%c("protein_coding","non-coding"),]
gg_barDE2 <- ggplot(df2, aes(x = Var2 , y = value,label=value,fill=Var3)) +
  #stat_summary(fun.data = mean_se,size=1.5) +
  #stat_summary(fun.y = mean, geom="point",colour="magenta", size=2) +
  scale_fill_manual(values = colssi,"") +
  geom_bar(stat="identity",position=position_dodge(width = .9)) +
  #scale_y_log10()+
  ylab("# upregulated genes") +
  theme_bw()+
  xlab("")+
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=16), 
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16),
        panel.grid.major = element_blank(),  # remove grid lines from panel
        panel.grid.minor = element_blank(),  # remove grid lines from panel
        strip.text=element_text(size=16),legend.position="top",legend.text=element_text(size=13),
        legend.title=element_text(size=16)) +
  facet_wrap(~type,scales="fixed",nrow = 1) + geom_text(vjust=-.1, color="black",position=position_dodge(width = .9))
pdf(file = "Figures/Figure2_C.pdf",width = 8,height = 4)
gg_barDE2+ coord_cartesian(ylim=c(1,270))
dev.off()

gns<-unique(multiDE_all$gene_id[which(multiDE_all$RNA_padj<.05 &multiDE_all$RNA_log2FC>0)])

df<-multiDE_all[multiDE_all$gene_id%in%gns,]
df$experiment<-as.numeric(gsub(df$experiment,pattern="h_vs_00h",replacement=""))
df$type<-"RNA-seq"

df1<-df
df01<-df1
df01$log2FC<-0
df01$experiment<-0
df1$log2FC<-df1$RNA_log2FC

gns<-unique(multiDE_all$gene_id[which(multiDE_all$Ribo_padj<.05 &multiDE_all$Ribo_log2FC>0)])

df<-multiDE_all[multiDE_all$gene_id%in%gns,]
df$experiment<-as.numeric(gsub(df$experiment,pattern="h_vs_00h",replacement=""))
df$type<-"Ribo-seq"
df2<-df
df2$log2FC<-df2$Ribo_log2FC
df02<-df2
df02$experiment<-0
df02$log2FC<-0
df<-rbind(df01,df02,df1,df2)

df<-df[,c("gene_biotype","experiment","log2FC","type")]
df$gene_biotype<-as.character(df$gene_biotype)
df$gene_biotype[df$gene_biotype=="novel_antisense_intr"]<-"novel_antisense_intronic"

df$gene_biotype<-factor(df$gene_biotype,levels=c("protein_coding","non-coding","DUX4_targets","novel_intergenic","novel_antisense","novel_antisense_intronic"))
df$type<-factor(as.character(df$type),levels=c("RNA-seq","Ribo-seq"))

colssi<-c("#0072B2", "#56B4E9", "#009E73", "grey30","grey", "goldenrod1")
names(colssi)<-c("novel_antisense", "novel_antisense_intronic", "novel_intergenic","protein_coding", "non-coding", "DUX4_targets")


gg_linefc <- ggplot(df, aes(x = experiment, y = log2FC,color=gene_biotype)) +
  stat_summary(fun.data = mean_se,size=1.5,geom="line") +
  stat_summary(fun.data = mean_se,size=1.5,geom="point")+
  stat_summary(fun.data = mean_se,size=1.5,geom="errorbar",width=0.2)+
  scale_color_manual(values = alpha(colssi,.7),"") +
  #geom_bar(stat="identity",position=position_dodge(width = .9)) +
  #scale_y_log10()+
  ylab("log2FC") +
  theme_bw()+
  xlab("Time after DUX4 induction")+
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=16), 
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16),
        panel.grid.major = element_blank(),  # remove grid lines from panel
        panel.grid.minor = element_blank(),  # remove grid lines from panel
        strip.text=element_text(size=16),legend.position="top",legend.text=element_text(size=14),
        legend.title=element_text(size=16))+ 
  facet_grid(~type,scales="free",space="free") 

pdf(file = "Figures/SupplementaryFigure4.pdf",width = 8,height = 5)
gg_linefc + scale_x_continuous(breaks = c(0,4,8,14),labels = c("0h","4h","8h","14h"))
dev.off()


df2<-df[df$type=="RNA-seq",]
gg_linefc2 <- ggplot(df2, aes(x = experiment, y = log2FC,color=gene_biotype)) +
  stat_summary(fun.data = mean_se,size=1.5,geom="line") +
  stat_summary(fun.data = mean_se,size=1.5,geom="point")+
  stat_summary(fun.data = mean_se,size=1.5,geom="errorbar",width=0.2)+
  scale_color_manual(values = alpha(colssi,.7),"") +
  #geom_bar(stat="identity",position=position_dodge(width = .9)) +
  #scale_y_log10()+
  ylab("log2FC") +
  theme_bw()+
  xlab("Time after DUX4 induction")+
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=16), 
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16),
        panel.grid.major = element_blank(),  # remove grid lines from panel
        panel.grid.minor = element_blank(),  # remove grid lines from panel
        strip.text=element_text(size=16),legend.position="top",legend.text=element_text(size=14),
        legend.title=element_text(size=16))+ 
  facet_grid(~type,scales="free",space="free") 


pdf(file = "Figures/Figure2_D.pdf",width = 6,height = 5)
gg_linefc2 + scale_x_continuous(breaks = c(0,4,8,14),labels = c("0h","4h","8h","14h"))
dev.off()



multiDE_all<-multiDE_all_idux4
chipseq_peaks <- rtracklayer::import("data/hg38_peaks.bed", format="BED")
proms <- promoters(GTF_annotation$txs,upstream = 1000,downstream=300)

common <- intersect(seqlevels(chipseq_peaks), seqlevels(proms))
chipseq_peaks <- keepSeqlevels(chipseq_peaks, common, pruning.mode="coarse")

proms <-split(proms,GTF_annotation$trann$gene_id[match(GTF_annotation$txs$tx_name,GTF_annotation$trann$transcript_id)])

gns<-unique(multiDE_all$gene_id[which(multiDE_all$RNA_padj<.05 &multiDE_all$RNA_log2FC>0)])
biots<-multiDE_all$gene_biotype[match(gns,multiDE_all$gene_id)]
allg<-unique(GTF_annotation$trann$gene_id)
freqs_all <- c()

for (g in unique(biots)){
  nong<-allg[!allg%in%gns[biots==g]]
  peaks<-chipseq_peaks[!chipseq_peaks%over%proms[nong]]
  promok<-proms[gns[biots==g]]
  freqs_all<-rbind(freqs_all,c(sum(promok%over%peaks),length(gns[biots==g])))
}
rownames(freqs_all)<-unique(biots)

df<-melt(freqs_all[,1]/freqs_all[,2]*100)
df$biotype<-rownames(df)
df$type<-"Upregulated in RNA-seq"
df1<-df


gns<-unique(multiDE_all$gene_id[which(multiDE_all$Ribo_padj<.05 &multiDE_all$Ribo_log2FC>0)])
biots<-multiDE_all$gene_biotype[match(gns,multiDE_all$gene_id)]
allg<-unique(GTF_annotation$trann$gene_id)
freqs_all <- c()

for (g in unique(biots)){
  nong<-allg[!allg%in%gns[biots==g]]
  peaks<-chipseq_peaks[!chipseq_peaks%over%proms[nong]]
  promok<-proms[gns[biots==g]]
  freqs_all<-rbind(freqs_all,c(sum(promok%over%peaks),length(gns[biots==g])))
}
rownames(freqs_all)<-unique(biots)

df<-melt(freqs_all[,1]/freqs_all[,2]*100)
df$biotype<-rownames(df)
df$type<-"Upregulated in Ribo-seq"
df2<-df
df<-rbind(df1,df2)
rownames(df)<-NULL


df$value<-round(df$value,digits = 2)
df$type<-factor(as.character(df$type),levels=c("Upregulated in RNA-seq","Upregulated in Ribo-seq"))
df$biotype<-as.character(df$biotype)
df$biotype[df$biotype=="novel_antisense_intr"]<-"novel_antisense_intronic"

df$biotype<-factor(df$biotype,levels=c("protein_coding","non-coding","DUX4_targets","novel_intergenic","novel_antisense","novel_antisense_intronic"))


gg_bar_chip <- ggplot(df, aes(x = biotype , y = value,label=value,fill=biotype)) +
  #stat_summary(fun.data = mean_se,size=1.5) +
  #stat_summary(fun.y = mean, geom="point",colour="magenta", size=2) +
  scale_fill_manual(values = colssi,"") +
  geom_bar(stat="identity") +
  ylab("% promoters overlapping DUX4 peaks") +
  theme_bw()+
  xlab("")+
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=16), 
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16),
        panel.grid.major = element_blank(),  # remove grid lines from panel
        panel.grid.minor = element_blank(),  # remove grid lines from panel
        strip.text=element_text(size=16),legend.position="none") +
  facet_wrap(~type,scales="free_x",nrow = 1) + geom_text(vjust=-.1, color="black")


df2<-df[df$type=="Upregulated in RNA-seq",]
gg_bar_chip2 <- ggplot(df2, aes(x = biotype , y = value,label=value,fill=biotype)) +
  #stat_summary(fun.data = mean_se,size=1.5) +
  #stat_summary(fun.y = mean, geom="point",colour="magenta", size=2) +
  scale_fill_manual(values = colssi,"") +
  geom_bar(stat="identity") +
  ylab("% promoters overlapping DUX4 peaks") +
  theme_bw()+
  xlab("")+
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=16), 
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16),
        panel.grid.major = element_blank(),  # remove grid lines from panel
        panel.grid.minor = element_blank(),  # remove grid lines from panel
        strip.text=element_text(size=16),legend.position="none") +
  facet_wrap(~type,scales="free_x",nrow = 1) + geom_text(vjust=-.1, color="black")


pdf(file = "Figures/SupplementaryFigure5.pdf",width = 8,height = 5)
gg_bar_chip
dev.off()

pdf(file = "Figures/Figure2_E.pdf",width = 6,height = 5)
gg_bar_chip2
dev.off()

#put the part with R1 R2 etc...

load("data/LRalso_UniqueORFq_Rannot")



multiDE_all1<-get(load("data/LRalso_04h_vs_00h_multiDE_results_summary_allgenes.RData"))
multiDE_all2<-get(load("data/LRalso_08h_vs_00h_multiDE_results_summary_allgenes.RData"))
multiDE_all3<-get(load("data/LRalso_14h_vs_00h_multiDE_results_summary_allgenes.RData"))

multiDE_all<-rbind(multiDE_all1$uniq[[1]],
                   multiDE_all2$uniq[[1]],
                   multiDE_all3$uniq[[1]])
multiDE_all$gene_biotype<-as.character(multiDE_all$gene_biotype)
multiDE_all$gene_biotype[grepl(multiDE_all$gene_id,pattern="^R1|^R2")]<-multiDE_all$biotype[grepl(multiDE_all$gene_id,pattern="^R1|^R2")]
multiDE_all$gene_biotype[multiDE_all$gene_id%in%target_genes]<-"DUX4_targets"
multiDE_all$gene_biotype[multiDE_all$gene_biotype%in%c("lncRNA","pseudogene","ncRNA")]<-"non-coding"

multiDE_all$novelty_gene<-"none"
multiDE_all$novelty_gene[grep(multiDE_all$gene_id,pattern="^R1")]<-"R1"
multiDE_all$novelty_gene[grep(multiDE_all$gene_id,pattern="^R1R2")]<-"R1R2"
multiDE_all$novelty_gene[grep(multiDE_all$gene_id,pattern="^R2")]<-"R2"

#IMPORTANT: here filter for not overlapping DE genes
multiDE_all<-multiDE_all[!multiDE_all$gene_overlapping_noDE,]

multiDE_all<-multiDE_all[multiDE_all$biotype!="artifact",]

multiDE_all$experiment<-gsub(multiDE_all$experiment,pattern="_vs_00h",replacement="")

df<-melt(table(multiDE_all$RNA_padj<.05,multiDE_all$experiment,multiDE_all$gene_biotype,multiDE_all$RNA_log2FC>0,multiDE_all$novelty_gene))
df<-df[df$Var1,]
df<-df[df$Var4,]
df<-df[df$Var5!="none",]
df$type<-"RNA-seq"
df1<-df
df<-melt(table(multiDE_all$Ribo_padj<.05,multiDE_all$experiment,multiDE_all$gene_biotype,multiDE_all$Ribo_log2FC>0,multiDE_all$novelty_gene))
df<-df[df$Var1,]
df<-df[df$Var4,]
df<-df[df$Var5!="none",]
df$type<-"Ribo-seq"
df2<-df
df<-rbind(df1,df2)
df$Var2<-gsub(as.character(df$Var2),pattern="_vs_00h",replacement="")
df$Var3<-as.character(df$Var3)
df$Var3[df$Var3=="novel_antisense_intr"]<-"novel_antisense_intronic"

df$Var3<-factor(df$Var3,levels=c("protein_coding","non-coding","DUX4_targets","novel_intergenic","novel_antisense","novel_antisense_intronic"))
df$type<-factor(as.character(df$type),levels=c("RNA-seq","Ribo-seq"))
df$Var5<-factor(as.character(df$Var5),levels=c("R1","R2","R1R2"))
df<-df[df$value>0,]
colssi<-c("#0072B2", "#56B4E9", "#009E73", "grey30","grey", "goldenrod1")
names(colssi)<-c("novel_antisense", "novel_antisense_intronic", "novel_intergenic","protein_coding", "non-coding", "DUX4_targets")

gg_barDELR <- ggplot(df, aes(x = Var2 , y = value,label=value,fill=Var3)) +
  #stat_summary(fun.data = mean_se,size=1.5) +
  #stat_summary(fun.y = mean, geom="point",colour="magenta", size=2) +
  scale_fill_manual(values = colssi,"") +
  geom_bar(stat="identity",position=position_dodge(width = .9)) +
  #scale_y_log10()+
  ylab("# upregulated genes") +
  theme_bw()+
  xlab("")+
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=16), 
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16),
        panel.grid.major = element_blank(),  # remove grid lines from panel
        panel.grid.minor = element_blank(),  # remove grid lines from panel
        strip.text=element_text(size=16),legend.position="top",legend.text=element_text(size=14),
        legend.title=element_text(size=16)) +
  facet_wrap(type~Var5,scales="fixed",nrow = 2) + geom_text(vjust=-.1, color="black",position=position_dodge(width = .9))
# # pdf(file = "Figures/Fig5DE_genes.# pdf",width = 9,height = 6)
gg_barDELR + coord_cartesian(ylim=c(1,max(df$value)*1.1))
# # dev.off()

dfall<-unique(multiDE_all[,c("gene_id","RNA_padj","RNA_log2FC","Ribo_padj","Ribo_log2FC","gene_biotype","novelty_gene")])
df<-melt(table(dfall$RNA_padj<.05,dfall$gene_biotype,dfall$RNA_log2FC>0,dfall$novelty_gene))
df<-df[df$Var1,]
df<-df[df$Var3,]
df<-df[df$Var4!="none",]
df$type<-"RNA-seq"
df1<-df
df<-melt(table(dfall$Ribo_padj<.05,dfall$gene_biotype,dfall$Ribo_log2FC>0,dfall$novelty_gene))
df<-df[df$Var1,]
df<-df[df$Var3,]
df<-df[df$Var4!="none",]
df$type<-"Ribo-seq"
df2<-df
df<-rbind(df1,df2)
df$Var2<-as.character(df$Var2)
df$Var2[df$Var2=="novel_antisense_intr"]<-"novel_antisense_intronic"

df$Var2<-factor(df$Var2,levels=c("protein_coding","non-coding","DUX4_targets","novel_intergenic","novel_antisense","novel_antisense_intronic"))
df$type<-factor(as.character(df$type),levels=c("RNA-seq","Ribo-seq"))
df$Var4<-factor(as.character(df$Var4),levels=c("R1","R2","R1R2"))
df<-df[df$value>0,]
colssi<-c("#0072B2", "#56B4E9", "#009E73", "grey30","grey", "goldenrod1")
names(colssi)<-c("novel_antisense", "novel_antisense_intronic", "novel_intergenic","protein_coding", "non-coding", "DUX4_targets")

gg_barDELRunq <- ggplot(df, aes(x = Var4 , y = value,label=value,fill=Var2)) +
  #stat_summary(fun.data = mean_se,size=1.5) +
  #stat_summary(fun.y = mean, geom="point",colour="magenta", size=2) +
  scale_fill_manual(values = colssi,"") +
  geom_bar(stat="identity",position=position_dodge(width = .9)) +
  #scale_y_log10()+
  ylab("# upregulated genes") +
  theme_bw()+
  xlab("")+
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=16), 
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16),
        panel.grid.major = element_blank(),  # remove grid lines from panel
        panel.grid.minor = element_blank(),  # remove grid lines from panel
        strip.text=element_text(size=16),legend.position="top",legend.text=element_text(size=14),
        legend.title=element_text(size=16)) +
  facet_wrap(~type,scales="free_y",nrow = 2) + geom_text(vjust=-.1, color="black",position=position_dodge(width = .9))
pdf(file = "Figures/Figure5_A.pdf",width = 7,height = 6)
#gg_barDELRunq + coord_cartesian(ylim=c(1,max(df$value)*1.1))
gg_barDELRunq
dev.off()




multiDE_all<-rbind(yao_tubes$uniq[[1]],
                   wang_G1$uniq[[1]],
                   wang_G2$uniq[[1]],
                   wang_G3$uniq[[1]],
                   wang_G4$uniq[[1]])

multiDE_all$gene_biotype<-as.character(multiDE_all$gene_biotype)
multiDE_all$gene_biotype[grepl(multiDE_all$gene_id,pattern="^R1")]<-multiDE_all$biotype[grepl(multiDE_all$gene_id,pattern="^R1")]
multiDE_all$gene_biotype[multiDE_all$gene_id%in%target_genes]<-"DUX4_targets"
multiDE_all$gene_biotype[multiDE_all$gene_biotype%in%c("lncRNA","pseudogene","ncRNA")]<-"non-coding"

#IMPORTANT: here filter for not overlapping DE genes
multiDE_all<-multiDE_all[!multiDE_all$gene_overlapping_noDE,]
multiDE_all<-multiDE_all[multiDE_all$biotype!="artifact",]
multiDE_all_patients<-multiDE_all

supp6<-multiDE_all_patients
supp6<-supp6[,!grepl(colnames(supp6),pattern="ovl|xy|len|Intron|dist|GCpc|n_exons")]
write.table(supp6, file="/group/calviello/Roberto/R2T2P_scripts_for_Figures_with_container/R2T2P_Figures_GitHub_repo/Supplementary_Tables/Supplementary_Table6.tsv",sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE)


df<-melt(table(multiDE_all$RNA_padj<.05,multiDE_all$experiment,multiDE_all$gene_biotype,multiDE_all$RNA_log2FC>0))
df<-df[df$Var1,]
df<-df[df$Var4,]
df$Var2<-gsub(df$Var2,pattern="_FSHD_vs_G1_Control",replacement="_vs_control")
df$Var2<-gsub(df$Var2,pattern="FSHD_vs_Control",replacement="FSHD_vs_control")

df$type<-"Patient biopsies"
df$type[grep(df$Var2,pattern="FSHD_vs_control")]<-"Patient-derived cells"

df$Var3<-as.character(df$Var3)
df$Var3[df$Var3=="novel_antisense_intr"]<-"novel_antisense_intronic"

df$Var3<-factor(df$Var3,levels=c("protein_coding","non-coding","DUX4_targets","novel_intergenic","novel_antisense","novel_antisense_intronic"))
df$type<-factor(as.character(df$type),levels=c("Patient-derived cells","Patient biopsies"))

df$Var2<-gsub(df$Var2,pattern="_vs_control",replacement="")

colssi<-c("#0072B2", "#56B4E9", "#009E73", "grey30","grey", "goldenrod1")
names(colssi)<-c("novel_antisense", "novel_antisense_intronic", "novel_intergenic","protein_coding", "non-coding", "DUX4_targets")

gg_barDE <- ggplot(df, aes(x = Var2 , y = value+1,label=value,fill=Var3)) +
  #stat_summary(fun.data = mean_se,size=1.5) +
  #stat_summary(fun.y = mean, geom="point",colour="magenta", size=2) +
  scale_fill_manual(values = colssi,"") +
  geom_bar(stat="identity",position=position_dodge(width = .9)) +
  scale_y_log10()+
  ylab("# upregulated genes") +
  theme_bw()+
  xlab("")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16),
        panel.grid.major = element_blank(),  # remove grid lines from panel
        panel.grid.minor = element_blank(),  # remove grid lines from panel
        strip.text=element_text(size=16),legend.position="top",legend.text=element_text(size=14),
        legend.title=element_text(size=16)) +
  facet_grid(~type,scales="free",space="free") + geom_text(vjust=-.1, color="black",position=position_dodge(width = .9))
pdf(file = "Figures/SupplementaryFigure13.pdf",width = 11,height = 6)
gg_barDE+ scale_x_discrete(expand = c(0.1, 0.4)) + coord_cartesian(ylim=c(1,max(df$value)*1.1))
dev.off()

df2<-df[!df$Var3%in%c("protein_coding","non-coding"),]
gg_barDE2 <- ggplot(df2, aes(x = Var2 , y = value+1,label=value,fill=Var3)) +
  #stat_summary(fun.data = mean_se,size=1.5) +
  #stat_summary(fun.y = mean, geom="point",colour="magenta", size=2) +
  scale_fill_manual(values = colssi,"") +
  geom_bar(stat="identity",position=position_dodge(width = .9)) +
  #scale_y_log10()+
  ylab("# upregulated genes") +
  theme_bw()+
  xlab("")+
  theme(axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16),
        panel.grid.major = element_blank(),  # remove grid lines from panel
        panel.grid.minor = element_blank(),  # remove grid lines from panel
        strip.text=element_text(size=16),legend.position="top",legend.text=element_text(size=13),
        legend.title=element_text(size=16)) +
  facet_grid(~type,scales="free",space="free") + geom_text(vjust=-.1, color="black",position=position_dodge(width = .9))
pdf(file = "Figures/Figure6_A.pdf",width = 8,height = 6)
gg_barDE2 + scale_x_discrete(expand = c(0.1, 0.4)) + coord_cartesian(ylim=c(1,max(df2$value)*1.1))
dev.off()

gns<-multiDE_all$gene_id[which(multiDE_all$RNA_padj<.05 & multiDE_all$RNA_log2FC>0)]
#why this?
df<-multiDE_all[multiDE_all$gene_id%in%gns,]

# this is better

df<-multiDE_all[which(multiDE_all$RNA_padj<.05 & multiDE_all$RNA_log2FC>0),]

df$experiment<-gsub(df$experiment,pattern="_FSHD_vs_G1_Control",replacement="_severity_vs_control")
df$experiment<-gsub(df$experiment,pattern="fshd_vs_ctrl",replacement="FSHD_vs_control")

df$type<-"Patient biopsies"
df$type[grep(df$experiment,pattern="FSHD_vs_control")]<-"Patient-derived cells"

df$gene_biotype<-as.character(df$gene_biotype)
df$gene_biotype[df$gene_biotype=="novel_antisense_intr"]<-"novel_antisense_intronic"

df$gene_biotype<-factor(df$gene_biotype,levels=c("protein_coding","non-coding","DUX4_targets","novel_intergenic","novel_antisense","novel_antisense_intronic"))
df$type<-factor(as.character(df$type),levels=c("Patient-derived cells","Patient biopsies"))

colssi<-c("#0072B2", "#56B4E9", "#009E73", "grey30","grey", "goldenrod1")
names(colssi)<-c("novel_antisense", "novel_antisense_intronic", "novel_intergenic","protein_coding", "non-coding", "DUX4_targets")


#change for testing

gg_boxfc <- ggplot(df, aes(x = experiment , y = RNA_log2FC,color=gene_biotype)) +
  stat_summary(fun.data = mean_se,size=1.5,position = position_dodge(width = .9)) +
  stat_summary(fun.data = mean_se,size=1.5,position = position_dodge(width = .9),geom="errorbar",width=.2) +
  
  scale_color_manual(values = alpha(colssi,.8),"") +
  #geom_bar(stat="identity",position=position_dodge(width = .9)) +
  #scale_y_log10()+
  ylab("RNA-seq log2FC") +
  theme_bw()+
  xlab("")+
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=16), 
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16),
        panel.grid.major = element_blank(),  # remove grid lines from panel
        panel.grid.minor = element_blank(),  # remove grid lines from panel
        strip.text=element_text(size=16),legend.position="top",legend.text=element_text(size=14),
        legend.title=element_text(size=16)) +
  facet_grid(~type,scales="free",space="free") 

# # pdf(file = "Figures/Fig6DE_genesFCs.# pdf",width = 11,height = 6)
gg_boxfc+ scale_x_discrete(expand = c(0.1, 0.4)) + coord_cartesian(ylim=c(0,10))
# # dev.off()

#make a plot with both genes and FCs, but with testing



gg_boxfc <- ggplot(df, aes(x = gene_biotype , y = RNA_log2FC,color=gene_biotype)) +
  stat_summary(fun.data = mean_se,size=1.5,position = position_dodge(width = .9)) +
  stat_summary(fun.data = mean_se,size=1.5,position = position_dodge(width = .9),geom="errorbar",width=.2) +
  
  scale_color_manual(values = alpha(colssi,.8),"") +
  #geom_bar(stat="identity",position=position_dodge(width = .9)) +
  #scale_y_log10()+
  ylab("RNA-seq log2FC") +
  theme_bw()+
  xlab("")+
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=16), 
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16),
        panel.grid.major = element_blank(),  # remove grid lines from panel
        panel.grid.minor = element_blank(),  # remove grid lines from panel
        strip.text=element_text(size=16),legend.position="top",legend.text=element_text(size=14),
        legend.title=element_text(size=16)) +
  facet_grid(~experiment,scales="free",space="free") 
# # pdf(file = "Figures/Fig6DE_genesFCs_compare.# pdf",width = 11,height = 6)
gg_boxfc + stat_compare_means(ref.group = "protein_coding",label="p.signif",method.args=list(alternative="g"))
# # dev.off()


