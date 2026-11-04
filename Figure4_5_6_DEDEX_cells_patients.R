library(cowplot)
library(zoo)
library(reshape2)
library(limma)
library(UpSetR)
library(readxl)

suppressPackageStartupMessages(library("RiboseQC"))
suppressPackageStartupMessages(library("ggfortify"))
trans <- function(x) {reso=-log(x, 10);reso[is.na(reso)]=0;return(reso)}
inv <- function(x) 10^(-x)

setwd("/group/calviello/Roberto/R2T2P_scripts_for_Figures_with_container/R2T2P_Figures_GitHub_repo")

source("data/new_riboseqc_all_HT_oct23_26.R")

#orfs<-get(load("/scratch/lorenzo.calviello/lopit_R2T2P/ORFquant/combined_forORFquant.RData_final_ORFquant_results"))
load("data/gencode.v47.annotation.gtf.gz_Rannot")

annot_files_SR<-get(load("data/annot_files_SR"))
annot_files_LR<-get(load("data/annot_files_LR"))
annot_files_LR2<-get(load("data/annot_files_LR2"))

multiDE_coding<-get(load("data/14h_vs_00h_multiDE_results_summary.RData"))
multiDE_all<-get(load("data/14h_vs_00h_multiDE_results_summary_allgenes.RData"))
multiDE_bins<-get(load("data/14h_vs_00h_multiDE_results_summary_subgene.RData"))

target_genes_file <- read_excel("data/Table_S1_DUX4_targets.xls", skip=1, col_names=T)
target_genes <- unique(GTF_annotation$trann$gene_id[sapply(strsplit(GTF_annotation$trann$gene_id,"[.]"),"[[",1) %in% sapply(strsplit(target_genes_file$Ensembl,"[.]"),"[[",1)])
target_genes2 <- unique(GTF_annotation$trann$gene_id[match(target_genes_file$gene.name,GTF_annotation$trann$gene_name)])
target_genes<-unique(c(target_genes,target_genes2))

load("data/SRonly_UniqueORFq_Rannot")

wang_G1<-get(load("data/G1_FSHD_vs_G1_Control_multiDE_results_summary.RData"))
wang_G2<-get(load("data/G2_FSHD_vs_G1_Control_multiDE_results_summary.RData"))
wang_G3<-get(load("data/G3_FSHD_vs_G1_Control_multiDE_results_summary.RData"))
wang_G4<-get(load("data/G4_FSHD_vs_G1_Control_multiDE_results_summary.RData"))
yao_tubes<-get(load("data/FSHD_vs_Control_multiDE_results_summary.RData"))

#df<-aa$uniq[[1]]

multiDE_bins1<-get(load("data/04h_vs_00h_multiDE_results_summary_subgene.RData"))
multiDE_bins2<-get(load("data/08h_vs_00h_multiDE_results_summary_subgene.RData"))
multiDE_bins3<-get(load("data/14h_vs_00h_multiDE_results_summary_subgene.RData"))
multiDE_bins1<-multiDE_bins1$uniq[[1]]
multiDE_bins2<-multiDE_bins2$uniq[[1]]
multiDE_bins3<-multiDE_bins3$uniq[[1]]
mcols(multiDE_bins1)<-cbind(mcols(multiDE_bins1),mcols(annot_files_SR$dex_bins_todex))
mcols(multiDE_bins2)<-cbind(mcols(multiDE_bins2),mcols(annot_files_SR$dex_bins_todex))
mcols(multiDE_bins3)<-cbind(mcols(multiDE_bins3),mcols(annot_files_SR$dex_bins_todex))

supp5<-mcols(multiDE_bins3)
supp5$coords<-as.character(multiDE_bins3)
supp5<-supp5[,!grepl(colnames(supp5),pattern="exonic_part|exon_name|exon_rank|exon_id|transcript_biotype|tx_id|original_column|_DE")]
supp5<-supp5[supp5$type%in%c("E"),]
supp5<-supp5[,-c(8,9,10,11,12)]
supp5$tx_name<-paste(unname(supp5$tx_name),collapse=";")
supp5<-supp5[which(supp5$RNA_padj<.05 | supp5$Ribo_padj<.05 ),]
write.table(supp5, file="/group/calviello/Roberto/R2T2P_scripts_for_Figures_with_container/R2T2P_Figures_GitHub_repo/Supplementary_Tables/Supplementary_Table5.tsv",sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE)


multiDE_coding1<-get(load("data/04h_vs_00h_multiDE_results_summary.RData"))
multiDE_coding2<-get(load("data/08h_vs_00h_multiDE_results_summary.RData"))
multiDE_coding3<-get(load("data/14h_vs_00h_multiDE_results_summary.RData"))
multiDE_coding1<-multiDE_coding1$uniq[[1]]
multiDE_coding2<-multiDE_coding2$uniq[[1]]
multiDE_coding3<-multiDE_coding3$uniq[[1]]

dex_evrna<-c(sum(multiDE_bins1$RNA_padj<.05,na.rm = T),sum(multiDE_bins2$RNA_padj<.05,na.rm = T),sum(multiDE_bins3$RNA_padj<.05,na.rm = T))
dex_evribo<-c(sum(multiDE_bins1$Ribo_padj<.05,na.rm = T),sum(multiDE_bins2$Ribo_padj<.05,na.rm = T),sum(multiDE_bins3$Ribo_padj<.05,na.rm = T))

dex_evrnag<-c(sum(multiDE_coding1$binsdexRNA_E<.05,na.rm = T),sum(multiDE_coding2$binsdexRNA_E<.05,na.rm = T),sum(multiDE_coding3$binsdexRNA_E<.05,na.rm = T))
dex_evribog<-c(sum(multiDE_coding1$binsdexRibo_E<.05,na.rm = T),sum(multiDE_coding2$binsdexRibo_E<.05,na.rm = T),sum(multiDE_coding3$binsdexRibo_E<.05,na.rm = T))

de_riborna_up<-c(sum(multiDE_coding1$RiboRNA_padj<.05 & multiDE_coding1$RiboRNA_log2FC>0,na.rm = T),sum(multiDE_coding2$RiboRNA_padj<.05 & multiDE_coding2$RiboRNA_log2FC>0,na.rm = T),sum(multiDE_coding3$RiboRNA_padj<.05 & multiDE_coding3$RiboRNA_log2FC>0,na.rm = T))
de_riborna_down<-c(sum(multiDE_coding1$RiboRNA_padj<.05 & multiDE_coding1$RiboRNA_log2FC<0,na.rm = T),sum(multiDE_coding2$RiboRNA_padj<.05 & multiDE_coding2$RiboRNA_log2FC<0,na.rm = T),sum(multiDE_coding3$RiboRNA_padj<.05 & multiDE_coding3$RiboRNA_log2FC<0,na.rm = T))

de_trascrrna_up<-c(sum(multiDE_coding1$RNA_padj<.05 & multiDE_coding1$RNA_log2FC>0,na.rm = T),sum(multiDE_coding2$RNA_padj<.05 & multiDE_coding2$RNA_log2FC>0,na.rm = T),sum(multiDE_coding3$RNA_padj<.05 & multiDE_coding3$RNA_log2FC>0,na.rm = T))
de_trascrrna_down<-c(sum(multiDE_coding1$RNA_padj<.05 & multiDE_coding1$RNA_log2FC<0,na.rm = T),sum(multiDE_coding2$RNA_padj<.05 & multiDE_coding2$RNA_log2FC<0,na.rm = T),sum(multiDE_coding3$RNA_padj<.05 & multiDE_coding3$RNA_log2FC<0,na.rm = T))



multiDE_bins_idux4<-multiDE_bins
bins_DE<-multiDE_bins$uniq[[1]]
mcols(bins_DE)<-cbind(mcols(bins_DE),mcols(annot_files_SR$dex_bins_todex))
bins_DE<-bins_DE[!bins_DE$gene_overlapping_noDE]
#bins_DE$RNA_padj[is.na(bins_DE$RNA_padj)]<-1
#bins_DE$Ribo_padj[is.na(bins_DE$Ribo_padj)]<-1

tocheck <- rev(c("cds","fiveutrs","threeutrs",
                 "ncIsof","ncRNAs","introns"))
bins_DE$region<-NA
regions <- list(reduce(unlist(GTF_annotation$cds_genes)), GTF_annotation$fiveutrs, GTF_annotation$threeutrs,
                GTF_annotation$ncIsof, GTF_annotation$ncRNAs, GTF_annotation$introns, GTF_annotation$intergenicRegions)
names(regions) <- c("cds","fiveutrs","threeutrs",
                    "ncIsof","ncRNAs","introns","intergenic")

for(i in rev(names(regions))){
  bins_DE$region[bins_DE%over%regions[[i]]]<-i
}

#TODO
#specify novel isoforms in coding genes and novel genes, and color them a little

bins_DE$region[bins_DE$type=="J"]<-"junction"
bins_DE$novel<-sum(!grepl(bins_DE$transcript_biotype,pattern="novel"))==0
bins_DE$ribo_sig<-bins_DE$Ribo_padj<.05
bins_DE$ribo_sig[is.na(bins_DE$ribo_sig)]<-F

bins_DE$rna_sig<-bins_DE$RNA_padj<.05
bins_DE$rna_sig[is.na(bins_DE$rna_sig)]<-F

vallli<-c("fiveutrs","cds","threeutrs","ncIsof","ncRNAs","novel_R1","novel_R2","novel_R1R2")

bins_ok<-bins_DE

bins_ok<-bins_ok[which(bins_ok$RNA_padj<.05 | bins_ok$Ribo_padj<.05 )]
bins_ok<-bins_ok[grep(bins_ok$gene_id,pattern="ENS"),]

bins_ok<-bins_ok[bins_ok$type%in%c("E"),]
bins_ok$region[bins_ok$novel]<-"novel_R1R2"
bins_ok$novel_R1<-sum(!grepl(bins_ok$tx_name,pattern="^R1."))==0
#R2
bins_ok$novel_R2<-sum(!grepl(bins_ok$tx_name,pattern="^R2."))==0

bins_ok$region[bins_ok$novel]<-"novel_R1R2"
bins_ok$region[bins_ok$novel & bins_ok$novel_R1]<-"novel_R1"
bins_ok$region[bins_ok$novel & bins_ok$novel_R2]<-"novel_R2"
bins_ok<-bins_ok[grep(bins_ok$gene_id,pattern="ENS"),]

#df_dex<-melt(table(bins_ok$region,bins_ok$ribo_sig,bins_ok$rna_sig,bins_ok$novel))
#colnames(df_dex)<-c("region","ribo_sig","rna_sig","novel","value")

bins_ok_g<-bins_ok
bins_ok_g<-mcols(bins_ok_g)[,c("gene_id","region","ribo_sig","rna_sig","novel")]
df_dex_ev<-melt(table(bins_ok_g$region,bins_ok_g$ribo_sig,bins_ok_g$rna_sig))
colnames(df_dex_ev)<-c("region","ribo_sig","rna_sig","value")

df_dex_ev$region<-as.character(df_dex_ev$region)
df_dex_ev<-df_dex_ev[which(df_dex_ev$ribo_sig | df_dex_ev$rna_sig),]
df_dex_evb<-df_dex_ev[which(df_dex_ev$ribo_sig & df_dex_ev$rna_sig),]
df_dex_ev<-df_dex_ev[-which(df_dex_ev$ribo_sig & df_dex_ev$rna_sig),]
df_dex_ev$value<-df_dex_ev$value+rep(df_dex_evb$value,2)
df_dex_ev$ribo_sig<-ifelse(df_dex_ev$ribo_sig,"Ribo-seq","RNA-seq")
df_dex_ev$ribo_sig<-factor(df_dex_ev$ribo_sig,levels = c("RNA-seq","Ribo-seq"))
df_dex_ev$region<-factor(df_dex_ev$region,levels = vallli)

#gene

bins_ok_g<-bins_ok
bins_ok_g<-mcols(bins_ok_g)[,c("gene_id","region","ribo_sig","rna_sig","novel")]
bins_ok_g<-unique(bins_ok_g)

df_dex_g_ok<-melt(table(bins_ok_g$region,bins_ok_g$ribo_sig,bins_ok_g$rna_sig))
colnames(df_dex_g_ok)<-c("region","ribo_sig","rna_sig","value")

df_dex_g_ok$region<-as.character(df_dex_g_ok$region)
df_dex_g_ok<-df_dex_g_ok[which(df_dex_g_ok$ribo_sig | df_dex_g_ok$rna_sig),]
df_dex_g_okb<-df_dex_g_ok[which(df_dex_g_ok$ribo_sig & df_dex_g_ok$rna_sig),]
df_dex_g_ok<-df_dex_g_ok[-which(df_dex_g_ok$ribo_sig & df_dex_g_ok$rna_sig),]
df_dex_g_ok$value<-df_dex_g_ok$value+rep(df_dex_g_okb$value,2)
df_dex_g_ok$ribo_sig<-ifelse(df_dex_g_ok$ribo_sig,"Ribo-seq","RNA-seq")
df_dex_g_ok$ribo_sig<-factor(df_dex_g_ok$ribo_sig,levels = c("RNA-seq","Ribo-seq"))
df_dex_g_ok$region<-factor(df_dex_g_ok$region,levels = vallli)

gg_bar <- ggplot(df_dex_g_ok, aes(x = region , y = value,label=value,fill=region)) +
  #stat_summary(fun.data = mean_se,size=1.5) +
  #stat_summary(fun.y = mean, geom="point",colour="magenta", size=2) +
  scale_fill_manual(values = c(rep("darkgrey",5),"red1","red3","darkred")) +
  geom_bar(stat="identity") +
  ylab("n genes with differential exonic bin") +
  theme_bw()+
  xlab("")+
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=16), 
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16),
        panel.grid.major = element_blank(),  # remove grid lines from panel
        panel.grid.minor = element_blank(),  # remove grid lines from panel
        strip.text=element_text(size=16),legend.position="none") +
  facet_wrap(~ribo_sig,scales="free_y",nrow = 2) + geom_text(vjust=-.1, color="black")



gg_bar_ev <- ggplot(df_dex_ev, aes(x = region , y = value,label=value,fill=region)) +
  #stat_summary(fun.data = mean_se,size=1.5) +
  #stat_summary(fun.y = mean, geom="point",colour="magenta", size=2) +
  scale_fill_manual(values = c(rep("darkgrey",5),"red1","red3","darkred")) +
  geom_bar(stat="identity") +
  ylab("n differential exonic bins") +
  theme_bw()+
  xlab("")+
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=16), 
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16),
        panel.grid.major = element_blank(),  # remove grid lines from panel
        panel.grid.minor = element_blank(),  # remove grid lines from panel
        strip.text=element_text(size=16),legend.position="none") +
  facet_wrap(~ribo_sig,scales="free_y",nrow = 2) + geom_text(vjust=-.1, color="black")


pdf(file = "Figures/Figure4_A.pdf",width = 7,height = 6)
gg_bar
dev.off()

# pdf(file = "Figures/Fig4_barDEX_events.pdf",width = 7,height = 6)
gg_bar_ev
# dev.off()




df_suja<-multiDE_coding$uniq[[1]]
df_suja<-df_suja[!df_suja$gene_overlapping_noDE,]
df_suja$Ribo_DEX_padj<-apply(cbind(df_suja$binsdexRibo_E,df_suja$binsdexRibo_J),1,function(x){min(x,na.rm=T)})
df_suja$Ribo_DEX_padj[which(df_suja$Ribo_DEX_padj==-Inf)]<-NA
df_suja$Ribo_DEX_padj[which(df_suja$Ribo_DEX_padj==Inf)]<-NA

df_suja$RNA_DEX_padj<-apply(cbind(df_suja$binsdexRNA_E,df_suja$binsdexRNA_J),1,function(x){min(x,na.rm=T)})
df_suja$RNA_DEX_padj[which(df_suja$RNA_DEX_padj==-Inf)]<-NA
df_suja$RNA_DEX_padj[which(df_suja$RNA_DEX_padj==Inf)]<-NA

suja_dexribo<-df_suja$gene_id[which(df_suja$Ribo_DEX_padj<.05)]
suja_dex<-df_suja$gene_id[which(df_suja$RNA_DEX_padj<.05)]

df<-multiDE_coding$uniq[[1]]
df<-df[!df$gene_overlapping_noDE,]

df$Ribo_DEX_padj<-apply(cbind(df$binsdexRibo_E,df$binsdexRibo_J),1,function(x){min(x,na.rm=T)})
df$Ribo_DEX_padj[which(df$Ribo_DEX_padj==-Inf)]<-NA
df$Ribo_DEX_padj[which(df$Ribo_DEX_padj==Inf)]<-NA
df$RNA_DEX_padj<-apply(cbind(df$binsdexRNA_E,df$binsdexRNA_J),1,function(x){min(x,na.rm=T)})
df$RNA_DEX_padj[which(df$RNA_DEX_padj==-Inf)]<-NA
df$RNA_DEX_padj[which(df$RNA_DEX_padj==Inf)]<-NA
df$Ribo_padj[is.na(df$Ribo_padj)]<-1
df$Ribo_DEX_padj[is.na(df$Ribo_DEX_padj)]<-1
df$RNA_padj[is.na(df$RNA_padj)]<-1
df$RNA_DEX_padj[is.na(df$RNA_DEX_padj)]<-1

df$category<-"Other"
#df$category[df$gene_id%in%suja_dex]<-"DTU_iDUX4"
df$category[df$gene_id%in%target_genes]<-"DUX4_target"
df$category<-factor(df$category,levels=c("DUX4_target","Other"))
df_plot<-data.frame(df)
dfok<-df[df$category=="DUX4_target",]

gdex_ribo <- ggplot(data = df_plot, aes(x=-log10(Ribo_padj), y=-log10(Ribo_DEX_padj), color=category)) + 
  #geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
  geom_point(size=2.5) + 
  #labs(x="-log10(RNA_DE pval)", y="-log10(RNA_DEX pval)") +
  labs(x="Adj. P-value\ngene level -  Ribo-seq", y="Adj. P-value\nsub-gene level - Ribo-seq") +
  scale_color_manual(values = alpha(colour = c("gold","grey"),c(.7,.7)),"") +
  geom_point(size=2.5,data=dfok) + 
  scale_y_continuous(breaks = c(0,20,40,60,80),labels = c(1,1e-20,1e-40,1e-60,1e-80))+
  scale_x_continuous(breaks = c(0,80,160,240,300),labels = c(1,1e-80,1e-160,1e-240,1e-300))+
  theme_bw() +
  theme(plot.margin = ggplot2::margin(.5, .5, .5, .5, "cm"),
        panel.background=element_blank(), 
        axis.ticks=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=15), 
        axis.text.y=element_text(size=15,angle=45),
        axis.title=element_text(size=15),legend.text=element_text(size=15),legend.position="top") 



gdex_rna <- ggplot(data = df_plot, aes(x=-log10(RNA_padj), y=-log10(RNA_DEX_padj), color=category)) + 
  #geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
  geom_point(size=2.5) + 
  #labs(x="-log10(RNA_DE pval)", y="-log10(RNA_DEX pval)") +
  labs(x="Adj. P-value\ngene level - RNA-seq", y="Adj. P-value\nsub-gene level - RNA-seq") +
  scale_color_manual(values = alpha(colour = c("gold","grey"),c(.7,.7)),"") +
  scale_y_continuous(breaks = c(0,50,100,150),labels = c(1,1e-50,1e-100,1e-150))+
  scale_x_continuous(breaks = c(0,80,160,240,300),labels = c(1,1e-80,1e-160,1e-240,1e-300))+
  geom_point(size=2.5,data=dfok) + 
  
  theme_bw() +
  theme(plot.margin = ggplot2::margin(.5, .5, .5, .5, "cm"),
        panel.background=element_blank(), 
        axis.ticks=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=15), 
        axis.text.y=element_text(size=15,angle=45),
        axis.title=element_text(size=15),legend.text=element_text(size=15),legend.position="top") 

pdf(file = "Figures/Figure4_C_and_SupplementaryFigure9_left.pdf",width = 7,height = 6)
gdex_ribo
dev.off()

pdf(file = "Figures/SupplementaryFigure9_right.pdf",width = 7,height = 6)
gdex_rna
dev.off()


df<-wang_G4$uniq[[1]]
df<-df[!df$gene_overlapping_noDE,]
df$RNA_DEX_padj<-apply(cbind(df$binsdexRNA_E,df$binsdexRNA_J),1,function(x){min(x,na.rm=T)})
df$RNA_DEX_padj[which(df$RNA_DEX_padj==-Inf)]<-NA
df$RNA_DEX_padj[which(df$RNA_DEX_padj==Inf)]<-NA

df$RNA_padj[is.na(df$RNA_padj)]<-1
df$RNA_DEX_padj[is.na(df$RNA_DEX_padj)]<-1
df2<-df

df$category<-"Other"
df$category[df$gene_id%in%suja_dexribo]<-"Sub-gene DE iDUX4 Ribo-seq"
df$category[df$gene_id%in%target_genes]<-"DUX4_target"
df$category<-factor(df$category,levels=c("DUX4_target","Sub-gene DE iDUX4 Ribo-seq","Other"))

df_plot<-data.frame(df)
dfok<-df[df$category%in%c("DUX4_target","Sub-gene DE iDUX4 Ribo-seq"),]

gdex_rna <- ggplot(data = df_plot, aes(x=-log10(RNA_padj), y=-log10(RNA_DEX_padj), color=category)) + 
  #geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
  geom_point(size=2.5) + 
  #labs(x="-log10(RNA_DE pval)", y="-log10(RNA_DEX pval)") +
  labs(x="Adj. P-value\ngene level - RNA-seq", y="Adj. P-value\nsub-gene level - RNA-seq") +
  scale_color_manual(values = alpha(colour = c("gold","cornflowerblue","grey"),c(.7,.7,.3)),"") +
  geom_point(size=2.5,data=dfok) + 
  scale_y_continuous(breaks = c(0,20,40,60,80),labels = c(1,1e-20,1e-40,1e-60,1e-80))+
  scale_x_continuous(breaks = c(0,20,40,60,80),labels = c(1,1e-20,1e-40,1e-60,1e-80))+
  theme_bw() +
  theme(plot.margin = ggplot2::margin(.5, .5, .5, .5, "cm"),
        panel.background=element_blank(), 
        axis.ticks=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=15), 
        axis.text.y=element_text(size=15,angle=45),
        axis.title=element_text(size=15),legend.text=element_text(size=15),legend.position="top") 

# pdf(file = "Figures/Fig6_dedex_WangRNA.pdf",width = 7,height = 6)
gdex_rna
# dev.off()


df2$category<-"Other"
df2$category[df2$gene_id%in%target_genes]<-"DUX4_target"
df2$category<-factor(df2$category,levels=c("DUX4_target","Other"))

df_plot2<-data.frame(df2)
dfok2<-df2[df2$category%in%c("DUX4_target"),]

gdex_rna2 <- ggplot(data = df_plot2, aes(x=-log10(RNA_padj), y=-log10(RNA_DEX_padj), color=category)) + 
  #geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
  geom_point(size=2.5) + 
  #labs(x="-log10(RNA_DE pval)", y="-log10(RNA_DEX pval)") +
  labs(x="Adj. P-value\ngene level - RNA-seq", y="Adj. P-value\nsub-gene level - RNA-seq") +
  scale_color_manual(values = alpha(colour = c("gold","grey"),c(.7,.7)),"") +
  geom_point(size=2.5,data=dfok2) + 
  scale_y_continuous(breaks = c(0,20,40,60,80),labels = c(1,1e-20,1e-40,1e-60,1e-80))+
  scale_x_continuous(breaks = c(0,20,40,60,80),labels = c(1,1e-20,1e-40,1e-60,1e-80))+
  theme_bw() +
  theme(plot.margin = ggplot2::margin(.5, .5, .5, .5, "cm"),
        panel.background=element_blank(), 
        axis.ticks=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=15), 
        axis.text.y=element_text(size=15,angle=45),
        axis.title=element_text(size=15),legend.text=element_text(size=15),legend.position="top") 

# pdf(file = "Figures/Fig6_dedex_WangRNA2.pdf",width = 7,height = 6)
gdex_rna2
# dev.off()


df<-yao_tubes$uniq[[1]]
df<-df[!df$gene_overlapping_noDE,]
df$RNA_DEX_padj<-apply(cbind(df$binsdexRNA_E,df$binsdexRNA_J),1,function(x){min(x,na.rm=T)})
df$RNA_DEX_padj[which(df$RNA_DEX_padj==-Inf)]<-NA
df$RNA_DEX_padj[which(df$RNA_DEX_padj==Inf)]<-NA

df$RNA_padj[is.na(df$RNA_padj)]<-1
df$RNA_DEX_padj[is.na(df$RNA_DEX_padj)]<-1
df2<-df

df1<-df
df1$dataset<-"Patient-derived myotubes"

df<-wang_G4$uniq[[1]]
df<-df[!df$gene_overlapping_noDE,]
df$RNA_DEX_padj<-apply(cbind(df$binsdexRNA_E,df$binsdexRNA_J),1,function(x){min(x,na.rm=T)})
df$RNA_DEX_padj[which(df$RNA_DEX_padj==-Inf)]<-NA
df$RNA_DEX_padj[which(df$RNA_DEX_padj==Inf)]<-NA

df$RNA_padj[is.na(df$RNA_padj)]<-1
df$RNA_DEX_padj[is.na(df$RNA_DEX_padj)]<-1
df2<-df
df2$dataset<-"Muscle biopsies"

df<-rbind(df1,df2)
df$dataset<-factor(df$dataset,levels=c("Patient-derived myotubes","Muscle biopsies"))
df$category<-"Other"
#df$category[df$gene_id%in%suja_dexribo]<-"Sub-gene DE iDUX4 Ribo-seq"
df$category[df$gene_id%in%target_genes]<-"DUX4_target"
df$category<-factor(df$category,levels=c("DUX4_target","Other"))


df_plot<-data.frame(df)
dfok<-df[df$category%in%c("DUX4_target"),]

gdex_rna_patients <- ggplot(data = df_plot, aes(x=-log10(RNA_padj), y=-log10(RNA_DEX_padj), color=category)) + 
  #geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
  geom_point(size=2.5) + 
  #labs(x="-log10(RNA_DE pval)", y="-log10(RNA_DEX pval)") +
  labs(x="Adj. P-value\ngene level - RNA-seq", y="Adj. P-value\nsub-gene level - RNA-seq") +
  scale_color_manual(values = alpha(colour = c("gold","grey"),c(.7,.7)),"") +
  geom_point(size=2.5,data=dfok) + 
  scale_y_continuous(breaks = c(0,20,40,60,80),labels = c(1,1e-20,1e-40,1e-60,1e-80))+
  scale_x_continuous(breaks = c(0,20,40,60,80),labels = c(1,1e-20,1e-40,1e-60,1e-80))+
  theme_bw() +
  theme(plot.margin = ggplot2::margin(.5, .5, .5, .5, "cm"),
        panel.background=element_blank(), 
        axis.ticks=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=15), 
        axis.text.y=element_text(size=15,angle=45),
        axis.title=element_text(size=15),legend.text=element_text(size=15),legend.position="top") 

pdf(file = "Figures/Figure6_C.pdf",width = 7,height = 8)
gdex_rna_patients + facet_wrap(~dataset,scales="free_y",nrow = 2)
dev.off()



list_all<-list(multiDE_coding,yao_tubes,wang_G1,wang_G2,wang_G3,wang_G4)
names(list_all)<-c("iDUX4","Yao_tubes","Wang_G1","Wang_G2","Wang_G3","Wang_G4")
list_all<-lapply(list_all,function(x){
  x<-x$uniq[[1]]
  x<-x[!x$gene_overlapping_noDE,]
  x
})
list_dex<-lapply(list_all,function(x){
  x$RNA_DEX_padj<-apply(cbind(x$binsdexRNA_E,x$binsdexRNA_J),1,function(x){min(x,na.rm=T)})
  x$RNA_DEX_padj[which(x$RNA_DEX_padj==-Inf)]<-NA
  x$RNA_DEX_padj[which(x$RNA_DEX_padj==Inf)]<-NA
  x$gene_id[which(x$RNA_DEX_padj<.01)]
})

list_dex<-list_dex[c("iDUX4","Wang_G4","Wang_G3","Wang_G2","Wang_G1","Yao_tubes")]
#list_dex$iDUX4<-suja_dex

# pdf(file = "Figures/Fig6_upset_all.pdf",width = 6,height = 5,onefile =F)
upset(fromList(list_dex), order.by = "freq",nsets = 6,keep.order = F,nintersects = 19)
# dev.off()
list_dex2<-list_dex[c("iDUX4","Wang_G4","Yao_tubes")]
#save(list_dex2,file = "/group/calviello/Lorenzo/R2T2P_paper/list_DEXgenes_patients.RData")
pdf(file = "Figures/Figure6_D.pdf",width = 5,height = 3.5,onefile =F)
upset(fromList(list_dex2), order.by = "freq",nsets = 6,keep.order = F,nintersects = 19)
dev.off()

#TODO

wang_G1_dex<-get(load("data/G1_FSHD_vs_G1_Control_multiDE_results_summary_subgene.RData"))
wang_G2_dex<-get(load("data/G2_FSHD_vs_G1_Control_multiDE_results_summary_subgene.RData"))
wang_G3_dex<-get(load("data/G3_FSHD_vs_G1_Control_multiDE_results_summary_subgene.RData"))
wang_G4_dex<-get(load("data/G4_FSHD_vs_G1_Control_multiDE_results_summary_subgene.RData"))
yao_tubes_dex<-get(load("data/FSHD_vs_Control_multiDE_results_summary_subgene.RData"))
multiDE_bins<-get(load("data/14h_vs_00h_multiDE_results_summary_subgene.RData"))


yao_bn<-yao_tubes_dex$uniq[[1]]
mcols(yao_bn)<-cbind(mcols(yao_bn),mcols(annot_files_SR$dex_bins_todex))

wang_bn<-wang_G4_dex$uniq[[1]]
mcols(wang_bn)<-cbind(mcols(wang_bn),mcols(annot_files_SR$dex_bins_todex))

supp7<-mcols(yao_bn)
supp7$coords<-as.character(yao_bn)
supp7<-supp7[,!grepl(colnames(supp7),pattern="exonic_part|exon_name|exon_rank|exon_id|transcript_biotype|tx_id|original_column|_DE")]
supp7<-supp7[supp7$type%in%c("E"),]
supp7<-supp7[,-c(8,9,10,11,12)]
supp7$tx_name<-paste(unname(supp7$tx_name),collapse=";")

supp7_w<-mcols(wang_bn)
supp7_w$coords<-as.character(wang_bn)

supp7_w<-supp7_w[,!grepl(colnames(supp7_w),pattern="exonic_part|exon_name|exon_rank|exon_id|transcript_biotype|tx_id|original_column|_DE")]
supp7_w<-supp7_w[supp7_w$type%in%c("E"),]
supp7_w<-supp7_w[,-c(8,9,10,11,12)]
supp7_w$tx_name<-paste(unname(supp7_w$tx_name),collapse=";")
supp7<-rbind(supp7,supp7_w)
supp7<-supp7[which(supp7$RNA_padj<.05 | supp7$Ribo_padj<.05 ),]
supp7$condition[supp7$condition=="FSHD"]<-"Yao_FSHD"
supp7$condition[supp7$condition=="G4_FSHD"]<-"Wang_FSHD"

write.table(supp7, file="/group/calviello/Roberto/R2T2P_scripts_for_Figures_with_container/R2T2P_Figures_GitHub_repo/Supplementary_Tables/Supplementary_Table7.tsv",sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE)


list_all_dex<-list(multiDE_bins,yao_tubes_dex,wang_G1_dex,wang_G2_dex,wang_G3_dex,wang_G4_dex)
names(list_all_dex)<-c("iDUX4","Yao_tubes","Wang_G1","Wang_G2","Wang_G3","Wang_G4")
list_all_dex<-lapply(list_all_dex,function(x){
  x<-x$uniq[[1]]
  x<-x[!x$gene_overlapping_noDE]
  x
})
warning("need a unique SR gtf for boh idux4 data and patients, redo!")
list_dex_events<-lapply(list_all_dex,function(x){
  x<-x[x%over%GTF_annotation$exons_txs[names(GTF_annotation$cds_txs)]]
  as.character(x[which(x$RNA_padj<.01)])
})

list_dex_events<-list_dex_events[c("iDUX4","Wang_G4","Wang_G3","Wang_G2","Wang_G1","Yao_tubes")]


# pdf(file = "Figures/Fig6_upset_events_all.pdf",width = 6,height = 5,onefile =F)
upset(fromList(list_dex_events), order.by = "freq",nsets = 6,keep.order = F,nintersects = 19)
# dev.off()
list_dex_events2<-list_dex_events[c("iDUX4","Wang_G4","Yao_tubes")]
#save(list_dex_events2,file = "/group/calviello/Lorenzo/R2T2P_paper/list_DEXevents_patients.RData")
# pdf(file = "Figures/Fig6_upset_events_few.pdf",width = 6,height = 5,onefile =F)
upset(fromList(list_dex_events2), order.by = "freq",nsets = 6,keep.order = F,nintersects = 19)
# dev.off()



cnt_regions_LR<-get(load("data/cnt_regions_LR2"))
multiDE_coding<-get(load("data/LRalso_14h_vs_00h_multiDE_results_summary.RData"))
multiDE_all<-get(load("data/LRalso_14h_vs_00h_multiDE_results_summary_allgenes.RData"))
multiDE_bins<-get(load("data/LRalso_14h_vs_00h_multiDE_results_summary_subgene.RData"))


bins_DE<-multiDE_bins$uniq[[1]]
mcols(bins_DE)<-cbind(mcols(bins_DE),mcols(annot_files_LR2$dex_bins_todex))
bins_DE<-bins_DE[!bins_DE$gene_overlapping_noDE]

tocheck <- rev(c("cds","fiveutrs","threeutrs",
                 "ncIsof","ncRNAs","introns"))
bins_DE$region<-NA
regions <- list(reduce(unlist(GTF_annotation$cds_genes)), GTF_annotation$fiveutrs, GTF_annotation$threeutrs,
                GTF_annotation$ncIsof, GTF_annotation$ncRNAs, GTF_annotation$introns, GTF_annotation$intergenicRegions)
names(regions) <- c("cds","fiveutrs","threeutrs",
                    "ncIsof","ncRNAs","introns","intergenic")

for(i in rev(names(regions))){
  bins_DE$region[bins_DE%over%regions[[i]]]<-i
}
bins_okall<-bins_DE
#TODO
#specify novel isoforms in coding genes and novel genes, and color them a little


bins_DE$region[bins_DE$type=="J"]<-"junction"
bins_DE$novel<-sum(!grepl(bins_DE$transcript_biotype,pattern="novel"))==0
bins_DE$ribo_sig<-bins_DE$Ribo_padj<.05
bins_DE$ribo_sig[is.na(bins_DE$ribo_sig)]<-F

bins_DE$rna_sig<-bins_DE$RNA_padj<.05
bins_DE$rna_sig[is.na(bins_DE$rna_sig)]<-F

bins_DE<-bins_DE[bins_DE$type%in%c("E"),]


bins_ok<-bins_DE

bins_ok<-bins_ok[which(bins_ok$RNA_padj<.05 | bins_ok$Ribo_padj<.05 )]
bins_ok<-bins_ok[grep(bins_ok$gene_id,pattern="ENS"),]

bins_ok<-bins_ok[bins_ok$type%in%c("E"),]
bins_ok$region[bins_ok$novel]<-"novel_R1R2"
bins_ok$novel_R1<-sum(!grepl(bins_ok$tx_name,pattern="^R1."))==0
#R2
bins_ok$novel_R2<-sum(!grepl(bins_ok$tx_name,pattern="^R2."))==0

bins_ok$region[bins_ok$novel]<-"novel_R1R2"
bins_ok$region[bins_ok$novel & bins_ok$novel_R1]<-"novel_R1"
bins_ok$region[bins_ok$novel & bins_ok$novel_R2]<-"novel_R2"
bins_ok<-bins_ok[grep(bins_ok$gene_id,pattern="ENS"),]

bins_ok_g<-bins_ok


bins_ok_g<-mcols(bins_ok_g)[,c("gene_id","region","ribo_sig","rna_sig","novel")]
df_dex_ev<-melt(table(bins_ok_g$region,bins_ok_g$ribo_sig,bins_ok_g$rna_sig))
colnames(df_dex_ev)<-c("region","ribo_sig","rna_sig","value")

df_dex_ev$region<-as.character(df_dex_ev$region)
df_dex_ev<-df_dex_ev[which(df_dex_ev$ribo_sig | df_dex_ev$rna_sig),]
df_dex_evb<-df_dex_ev[which(df_dex_ev$ribo_sig & df_dex_ev$rna_sig),]
df_dex_ev<-df_dex_ev[-which(df_dex_ev$ribo_sig & df_dex_ev$rna_sig),]
df_dex_ev$value<-df_dex_ev$value+rep(df_dex_evb$value,2)
df_dex_ev$ribo_sig<-ifelse(df_dex_ev$ribo_sig,"Ribo-seq","RNA-seq")
df_dex_ev$ribo_sig<-factor(df_dex_ev$ribo_sig,levels = c("RNA-seq","Ribo-seq"))
df_dex_ev$region<-factor(df_dex_ev$region,levels = vallli)

#gene

bins_ok_g<-bins_ok
bins_ok_g<-mcols(bins_ok_g)[,c("gene_id","region","ribo_sig","rna_sig","novel")]
bins_ok_g<-unique(bins_ok_g)

df_dex_g_ok<-melt(table(bins_ok_g$region,bins_ok_g$ribo_sig,bins_ok_g$rna_sig))
colnames(df_dex_g_ok)<-c("region","ribo_sig","rna_sig","value")

df_dex_g_ok$region<-as.character(df_dex_g_ok$region)
df_dex_g_ok<-df_dex_g_ok[which(df_dex_g_ok$ribo_sig | df_dex_g_ok$rna_sig),]
df_dex_g_okb<-df_dex_g_ok[which(df_dex_g_ok$ribo_sig & df_dex_g_ok$rna_sig),]
df_dex_g_ok<-df_dex_g_ok[-which(df_dex_g_ok$ribo_sig & df_dex_g_ok$rna_sig),]
df_dex_g_ok$value<-df_dex_g_ok$value+rep(df_dex_g_okb$value,2)
df_dex_g_ok$ribo_sig<-ifelse(df_dex_g_ok$ribo_sig,"Ribo-seq","RNA-seq")
df_dex_g_ok$ribo_sig<-factor(df_dex_g_ok$ribo_sig,levels = c("RNA-seq","Ribo-seq"))
df_dex_g_ok$region<-factor(df_dex_g_ok$region,levels = vallli)

gg_bar <- ggplot(df_dex_g_ok, aes(x = region , y = value,label=value,fill=region)) +
  #stat_summary(fun.data = mean_se,size=1.5) +
  #stat_summary(fun.y = mean, geom="point",colour="magenta", size=2) +
  scale_fill_manual(values = c(rep("darkgrey",5),"red1","red3","darkred")) +
  geom_bar(stat="identity") +
  ylab("n genes with differential exonic bin") +
  theme_bw()+
  xlab("")+
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=16), 
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16),
        panel.grid.major = element_blank(),  # remove grid lines from panel
        panel.grid.minor = element_blank(),  # remove grid lines from panel
        strip.text=element_text(size=16),legend.position="none") +
  facet_wrap(~ribo_sig,scales="free_y",nrow = 2) + geom_text(vjust=-.1, color="black")
# pdf(file = "Figures/Fig5_barDEX.pdf",width = 9,height = 7)
gg_bar
# dev.off()


gg_bar_ev <- ggplot(df_dex_ev, aes(x = region , y = value,label=value,fill=region)) +
  #stat_summary(fun.data = mean_se,size=1.5) +
  #stat_summary(fun.y = mean, geom="point",colour="magenta", size=2) +
  scale_fill_manual(values = c(rep("darkgrey",5),"red1","red3","darkred")) +
  geom_bar(stat="identity") +
  ylab("n of differential exonic bins") +
  theme_bw()+
  xlab("")+
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=16), 
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16),
        panel.grid.major = element_blank(),  # remove grid lines from panel
        panel.grid.minor = element_blank(),  # remove grid lines from panel
        strip.text=element_text(size=16),legend.position="none") +
  facet_wrap(~ribo_sig,scales="free_y",nrow = 2) + geom_text(vjust=-.1, color="black")

# pdf(file = "Figures/Fig5_barDEX_events.pdf",width = 9,height = 7)
gg_bar_ev
# dev.off()


df_dex_g_ok_few<-df_dex_g_ok[grep(df_dex_g_ok$region,pattern="novel"),]
gg_bar <- ggplot(df_dex_g_ok_few, aes(x = region , y = value,label=value,fill=region)) +
  #stat_summary(fun.data = mean_se,size=1.5) +
  #stat_summary(fun.y = mean, geom="point",colour="magenta", size=2) +
  scale_fill_manual(values = c("red1","red3","darkred")) +
  geom_bar(stat="identity") +
  ylab("n genes with differential exonic bin") +
  theme_bw()+
  xlab("")+
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=16), 
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16),
        panel.grid.major = element_blank(),  # remove grid lines from panel
        panel.grid.minor = element_blank(),  # remove grid lines from panel
        strip.text=element_text(size=16),legend.position="none") +
  facet_wrap(~ribo_sig,scales="free_y",nrow = 2) + geom_text(vjust=-.1, color="black")
pdf(file = "Figures/Figure5_C.pdf",width = 5,height = 5)
gg_bar
dev.off()

df_dex_ev_few<-df_dex_ev[grep(df_dex_ev$region,pattern="novel"),]

gg_bar_ev <- ggplot(df_dex_ev_few, aes(x = region , y = value,label=value,fill=region)) +
  #stat_summary(fun.data = mean_se,size=1.5) +
  #stat_summary(fun.y = mean, geom="point",colour="magenta", size=2) +
  scale_fill_manual(values = c("red1","red3","darkred")) +
  geom_bar(stat="identity") +
  ylab("n of differential exonic bins") +
  theme_bw()+
  xlab("")+
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=16), 
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16),
        panel.grid.major = element_blank(),  # remove grid lines from panel
        panel.grid.minor = element_blank(),  # remove grid lines from panel
        strip.text=element_text(size=16),legend.position="none") +
  facet_wrap(~ribo_sig,scales="free_y",nrow = 2) + geom_text(vjust=-.1, color="black")

# pdf(file = "Figures/Fig5_barDEX_events_few.pdf",width = 5,height = 5)
gg_bar_ev 
# dev.off()
