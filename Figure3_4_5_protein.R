library(cowplot)
library(zoo)
library(reshape2)
library(limma)
library(limpa)
library(UpSetR)
library(ggVennDiagram)
suppressMessages(library("readxl"))

suppressPackageStartupMessages(library("ggfortify"))
suppressPackageStartupMessages(library("ggpubr"))

trans <- function(x) {reso=-log(x, 10);reso[is.na(reso)]=0;return(reso)}
inv <- function(x) 10^(-x)

setwd("/group/calviello/Roberto/R2T2P_scripts_for_Figures_with_container/R2T2P_Figures_GitHub_repo")

source("data/new_riboseqc_all_HT_oct23_26.R")
load("data/SRonly_UniqueORFq_Rannot")

target_genes_file <- read_excel("data/Table_S1_DUX4_targets.xls", skip=1, col_names=T)
target_genes <- unique(GTF_annotation$trann$gene_id[sapply(strsplit(GTF_annotation$trann$gene_id,"[.]"),"[[",1) %in% sapply(strsplit(target_genes_file$Ensembl,"[.]"),"[[",1)])
target_genes2 <- unique(GTF_annotation$trann$gene_id[match(target_genes_file$gene.name,GTF_annotation$trann$gene_name)])
target_genes<-unique(c(target_genes,target_genes2))

annot_files<-get(load("data/annot_files_SR"))

multiDE_bins<-get(load("data/14h_vs_00h_multiDE_results_summary_subgene.RData"))


list_DE<-get(load("data/SRonly_UniqueORFq_combined_database_search_proteomics_multiDE.RData"))
list_intensities<-get(load("data/SRonly_UniqueORFq_combined_database_search_proteomics_all_intensities.RData"))
peptide_data<-get(load("data/SRonly_UniqueORFq_combined_database_search_peptide_results_annotated.RData"))

supp4<-peptide_data$metadata
supp4<-supp4[,grepl(colnames(supp4),pattern="Proteins_collapsed|novelty_category|gene_name|gene_id")]
supp4$Sequence<-rownames(supp4)
rownames(supp4)<-NULL
supp4<-data.frame(supp4)
write.table(supp4, file="/group/calviello/Roberto/R2T2P_scripts_for_Figures_with_container/R2T2P_Figures_GitHub_repo/Supplementary_Tables/Supplementary_Table4.tsv",sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE)

df_meta<-peptide_data$metadata

multi<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))>1,]
df_meta<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))==1,]

df_meta$gene_biotype[grepl(df_meta$gene_name,pattern="^IG") & df_meta$gene_biotype=="other"]<-"protein_coding"
df_meta$gene_biotype[df_meta$gene_biotype!="protein_coding"]<-"non-coding"
df_meta$DUX4_target<-df_meta$gene_id%in%target_genes
df_meta$gene_biotype[df_meta$DUX4_target]<-"DUX4_target"

df_meta_unq<-unique(df_meta[,c("gene_id","gene_biotype","novelty_category")])


df<-melt(table(df_meta$gene_biotype,df_meta$novelty_category))
df<-df[df$value>0,]
df$id_type<-"Peptides"
df$Var1<-factor(df$Var1,levels=names(sort(table(df_meta$gene_biotype),decreasing = T)))
df2<-melt(table(df_meta_unq$gene_biotype,df_meta_unq$novelty_category))
df2<-df2[df2$value>0,]
df2$Var1<-factor(df2$Var1,levels=names(sort(table(df_meta$gene_biotype),decreasing = T)))
df2$id_type<-"Genes"
df<-rbind(df,df2)

gg_ids <- ggplot(data = df, aes(x=Var1, y=value+1, fill=Var2,label=value)) + 
  #geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
  geom_bar(stat="identity", position=position_dodge(width = .9)) + 
  #geom_col(stat="identity", position=position_dodge(width = .9),preserve = c("total", "single")) + 
  
  labs(x="", y="# identifications") +
  scale_fill_manual(values = c("grey","orange","red"),"peptide_type") +
  theme_bw() +
  scale_y_log10() +
  theme(plot.margin = ggplot2::margin(.5, .5, .5, .5, "cm"),
        panel.background=element_blank(), 
        axis.ticks=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=14), 
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20),
        #legend.position="bottom",
        legend.text=element_text(size=18),
        legend.title=element_text(size=20),
        panel.grid.major.y = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.x = element_line(linewidth =.1, color="dark grey" ,linetype=2)) +
  geom_text(position=position_dodge(width = .9),vjust = -0.25) + facet_wrap(~id_type,nrow = 2,scales = "free_y")

# pdf(file = "Figures/Fig3_ids.pdf",width = 9,height = 6)
# gg_ids
# dev.off()


df_meta<-peptide_data$metadata

multi<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))>1,]
df_meta<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))==1,]

df_meta$gene_biotype[grepl(df_meta$gene_name,pattern="^IG") & df_meta$gene_biotype=="other"]<-"protein_coding"
df_meta$gene_biotype[df_meta$gene_biotype!="protein_coding"]<-"non-coding"
df_meta$DUX4_target<-df_meta$gene_id%in%target_genes
df_meta$gene_biotype[df_meta$DUX4_target]<-"DUX4_target"

df_meta_unq<-unique(df_meta[,c("gene_id","gene_biotype","novelty_category")])

DE_res<-list_DE$raw$DE

pepsig<-DE_res$Peptides[["DUX4_DMSO_14h-Control_DMSO_14h"]]
pepsig<-rownames(pepsig)[which(pepsig$padj<.05 & pepsig$logFC>0)]
ginsig<-DE_res$gene_id_novelty[["DUX4_DMSO_14h-Control_DMSO_14h"]]
ginsig<-rownames(ginsig)[which(ginsig$padj<.05 & ginsig$logFC>0)]

df_meta$sig_up<-rownames(df_meta)%in%pepsig
df_meta_unq$sig_up<-paste(df_meta_unq$gene_id,df_meta_unq$novelty_category,sep=";")%in%ginsig

df<-melt(table(df_meta$gene_biotype,df_meta$novelty_category,df_meta$sig_up))

df$id_type<-"Peptides"
df$Var1<-factor(df$Var1,levels=names(sort(table(df_meta$gene_biotype),decreasing = T)))

df2<-melt(table(df_meta_unq$gene_biotype,df_meta_unq$novelty_category,df_meta_unq$sig_up))

df2$Var1<-factor(df2$Var1,levels=names(sort(table(df_meta$gene_biotype),decreasing = T)))
df2$id_type<-"Genes"

df<-rbind(df,df2)

df$value2<-df$value
df$value2[!df$Var3]<-df$value[!df$Var3]+df$value[df$Var3]

df<-df[df$value>0,]
df$sigup<-"all"
df$sigup[df$Var3]<-"Upregulated"
df$sigup <- factor(df$sigup, levels = c("all", "Upregulated"))

gg_idsig <- ggplot(data = df, aes(x=Var2, y=value+1, fill=Var2,label=value2,alpha=sigup)) + 
  #geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
  geom_bar(stat="identity", position="stack") + 
  #geom_col(stat="identity", position=position_dodge(width = .9),preserve = c("total", "single")) + 
  
  labs(x="", y="# identifications") +
  scale_fill_manual(values = c("grey","orange","red"),"peptide_type") +
  theme_bw() +
  scale_y_log10() +
  scale_alpha_manual(values = c(.5,1),"")+
  theme(plot.margin = ggplot2::margin(.5, .5, .5, .5, "cm"),
        panel.background=element_blank(), 
        axis.ticks=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=14), 
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20),
        #legend.position="bottom",
        legend.text=element_text(size=18),
        legend.title=element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_text(aes(x=Var2, y=(value+1)*1.1,label=value2),position="stack",check_overlap = F,show.legend = FALSE) + facet_grid(id_type~Var1,scales = "free_x",space = "free_x")
pdf(file = "Figures/Figure3_B.pdf",width = 9,height = 6)
gg_idsig
dev.off()

# # here summarize n of DE events
DE_res<-list_DE$raw$DE
# 
# df<-DE_res$gene_id_novelty[["DUX4_DMSO_04h-Control_DMSO_04h"]]
# 
# gns<-sapply(strsplit(rownames(df),";"),"[[",1)
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# #df$novelty<-df_meta$novelty_category[match(gns,df_meta$gene_id)]
# df$novelty<-sapply(strsplit(rownames(df),";"),"[[",2)
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Genes"
# df$time<-"04h"
# df1<-df
# 
# df<-DE_res$gene_id_novelty[["DUX4_DMSO_08h-Control_DMSO_08h"]]
# 
# gns<-sapply(strsplit(rownames(df),";"),"[[",1)
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# #df$novelty<-df_meta$novelty_category[match(gns,df_meta$gene_id)]
# df$novelty<-sapply(strsplit(rownames(df),";"),"[[",2)
# 
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Genes"
# df$time<-"08h"
# df2<-df
# 
df<-DE_res$gene_id_novelty[["DUX4_DMSO_14h-Control_DMSO_14h"]]

gns<-sapply(strsplit(rownames(df),";"),"[[",1)
df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
df$gene_id<-gns
df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
df$novelty<-sapply(strsplit(rownames(df),";"),"[[",2)
df$category<-paste(df$gene_biotype,df$novelty,sep=";")
df$sig<-F
df$sig[df$padj<.05]<-T
df<-df[order(df$padj,decreasing = F),]
df<-df[which(df$padj<.05),]
df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
df$id_type="Genes"
df$time<-"14h"
df3<-df

# dfg<-rbind(df1,df2,df3)
dfg <- df3

# df<-DE_res$Peptides[["DUX4_DMSO_04h-Control_DMSO_04h"]]
# gns<-df_meta$gene_id
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# df$novelty<-df_meta$novelty_category
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Peptides"
# df$time<-"04h"
# df1<-df
# 
# df<-DE_res$Peptides[["DUX4_DMSO_08h-Control_DMSO_08h"]]
# 
# gns<-df_meta$gene_id
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# df$novelty<-df_meta$novelty_category
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Peptides"
# df$time<-"08h"
# df2<-df
# 
df<-DE_res$Peptides[["DUX4_DMSO_14h-Control_DMSO_14h"]]

gns<-df_meta$gene_id
df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
df$gene_id<-gns
df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
df$novelty<-df_meta$novelty_category
df$category<-paste(df$gene_biotype,df$novelty,sep=";")
df$sig<-F
df$sig[df$padj<.05]<-T
df<-df[order(df$padj,decreasing = F),]
df<-df[which(df$padj<.05),]
df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
df$id_type="Peptides"
df$time<-"14h"
df3<-df
# dfp<-rbind(df1,df2,df3)
dfp <- df3

df<-rbind(dfg,dfp)
dfall<-df
df<-df[df$Var3,]

df<-df[df$time=="14h",]
#df<-df[df$value>0,]
df$value[df$value==0]<-NA

df$Var1<-factor(as.character(df$Var1),levels=c("protein_coding","DUX4_target","non-coding"))
gg_deres <- ggplot(data = df, aes(x=Var1, y=value+1, fill=Var2,label=value)) +
  #geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
  geom_bar(stat="identity", position=position_dodge(width = .9)) +
  labs(x="", y="# significantly upregulated") +
  scale_fill_manual(values = c("grey","orange","red"),"peptide_type") +
  theme_bw() +
  scale_y_log10() +
  theme(plot.margin = ggplot2::margin(.5, .5, .5, .5, "cm"),
        panel.background=element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=14),
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20),
        legend.text=element_text(size=14),
        panel.grid.major.y = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.x = element_line( size=.1, color="dark grey" ,linetype=2)) +
  geom_text(position=position_dodge(width = .9),vjust = -0.25) + facet_wrap(id_type~time,nrow = 2,scales = "fixed")

# pdf(file = "Figures/Fig3_deres.pdf",width = 8,height = 6)
# gg_deres
# dev.off()


# df<-dfall
# df<-df[df$Var3,]
# 
# #df<-df[df$time=="14h",]
# #df<-df[df$value>0,]
# df$value[df$value==0]<-NA
# 
# df$Var1<-factor(as.character(df$Var1),levels=c("protein_coding","DUX4_target","non-coding"))
# gg_deresalltime <- ggplot(data = df, aes(x=Var1, y=value+1, fill=Var2,label=value)) +
#   #geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
#   geom_bar(stat="identity", position=position_dodge(width = .9)) +
#   labs(x="", y="# significantly upregulated") +
#   scale_fill_manual(values = c("grey","orange","red"),"peptide_type") +
#   theme_bw() +
#   scale_y_log10() +
#   theme(plot.margin = ggplot2::margin(.5, .5, .5, .5, "cm"),
#         panel.background=element_blank(),
#         axis.ticks=element_blank(),
#         axis.text.x=element_text(angle=45, hjust=1, size=14),
#         axis.text.y=element_text(size=20),
#         axis.title=element_text(size=20),
#         legend.text=element_text(size=14),
#         panel.grid.major.y = element_blank() ,
#         # explicitly set the horizontal lines (or they will disappear too)
#         panel.grid.major.x = element_line( size=.1, color="dark grey" ,linetype=2)) +
#   geom_text(position=position_dodge(width = .9),vjust = -0.25) + facet_wrap(id_type~time,nrow = 2,scales = "fixed")
# 
# # pdf(file = "Figures/Fig3_deresalltime.pdf",width = 10,height = 6)
# gg_deresalltime
# # dev.off()




#make MA plots

df<-DE_res$gene_id_novelty[["DUX4_DMSO_14h-Control_DMSO_14h"]]
df$avg_expr<-apply(list_intensities$raw$gene_id_novelty[,grep(colnames(list_intensities$raw$gene_id_novelty),pattern="Control_DMSO")],1,function(x){
  mean(x[!x%in%c(NA,-Inf,Inf)])
})

gns<-sapply(strsplit(rownames(df),";"),"[[",1)
df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
df$gene_id<-gns
df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]

#ERROR, match only assigns 1, one gene can be in multiple novelty categories
#df$novelty<-df_meta$novelty_category[match(gns,df_meta$gene_id)]

df$novelty<-sapply(strsplit(rownames(df),";"),"[[",2)

gns<-sapply(strsplit(rownames(df),";"),"[[",2)
df$category<-paste(df$gene_biotype,df$novelty,sep=";")
df$sig<-F
df$sig[df$padj<.05]<-T
df<-df[order(df$padj,decreasing = F),]
df$type="Genes"
dfg<-df

df<-DE_res$Peptides[["DUX4_DMSO_14h-Control_DMSO_14h"]]
df$avg_expr<-apply(list_intensities$raw$Peptides[,grep(colnames(list_intensities$raw$Peptides),pattern="Control_DMSO")],1,function(x){
  mean(x[!x%in%c(NA,-Inf,Inf)])
})
df$gene_id<-df_meta$gene_id
gns<-df$gene_id
df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
df$gene_id<-gns

df$gene_biotype<-df_meta$gene_biotype

df$novelty<-df_meta$novelty_category

df$category<-paste(df$gene_biotype,df$novelty,sep=";")
df$sig<-.1
df$sig[df$padj<.05]<-.5
df<-df[order(df$padj,decreasing = F),]
df$type<-"Peptides"

df<-rbind(dfg,df)
df$sig<-.1
df$sig[df$padj<.05]<-.5
df$gene_biotype<-factor(as.character(df$gene_biotype),levels = c("protein_coding","DUX4_target","non-coding"))
colsi<-alpha(c("grey","gold","cornflowerblue"),.8)

names(colsi)<-c("protein_coding","DUX4_target","non-coding")


colsi2<-alpha(c("grey","orange","red"),.75)

names(colsi2)<-c("in_ref","novel","novel_tx_R1")
df2<-df[which(df$padj<.05),]

ma<-ggplot(df2,aes(x=avg_expr,y=logFC,size=padj,color=gene_biotype,shape=novelty)) + geom_point()
ma<-ma + theme_bw() +
  ylab("log2FC\n14h vs DMSO") +
  xlab("Average Intensity") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18),legend.text=element_text(size=12))  +
  scale_color_manual(values = colsi,"category") +
  scale_size_continuous("Adj. p-value",trans = scales::trans_new(name = "test",transform = trans,inverse = inv, 
                                                                 domain = c(1e-100, Inf)),breaks = c(10e-51,10e-11,10e-5,10e-3,10e-2))
#ma<-ma + coord_cartesian(xlim=c(12,35))
maplot1<-ma + facet_wrap(~type,ncol = 1,scales="free_y") + coord_cartesian(xlim=c(13,32))
df2<-df[which(df$padj<.05),]
ma<-ggplot(df2,aes(x=avg_expr,y=logFC,size=padj,color=novelty,shape=gene_biotype)) + geom_point()
ma<-ma + theme_bw() +
  ylab("log2FC\n14h vs DMSO") +
  xlab("Average Intensity(log2)") +
  theme(axis.title.x = element_text(size=22),axis.text.x  = element_text(angle=45, vjust=0.5, size=18)) +
  theme(axis.title.y = element_text(size=24),axis.text.y  = element_text(angle=45, vjust=0.5, size=18),legend.text=element_text(size=12))  +
  scale_color_manual(values = colsi2,"category") +
  scale_size_continuous("Adj. p-value",trans = scales::trans_new(name = "test",transform = trans,inverse = inv, 
                                                                 domain = c(1e-100, Inf)),breaks = c(10e-51,10e-11,10e-5,10e-3,10e-2))
#ma<-ma + coord_cartesian(xlim=c(12,35))
maplot2<-ma + facet_wrap(~type,ncol = 1,scales="free_y") + coord_cartesian(xlim=c(13,32))


# pdf(file = "Figures/Fig3_maplot1.pdf",width = 9,height = 6)
# maplot1
# dev.off()

pdf(file = "Figures/Figure3_C.pdf",width = 9,height = 6)
maplot2
dev.off()



bins_DE<-multiDE_bins$uniq[[1]]
mcols(bins_DE)<-cbind(mcols(bins_DE),mcols(annot_files$dex_bins_todex))
bins_DE <- bins_DE[!bins_DE$gene_overlapping_noDE]
bins_DE<-bins_DE[bins_DE$type=="E"]

diffpep<-list_DE$raw$diffsplice$Peptides[["DUX4_DMSO_14h-Control_DMSO_14h"]]
diffpep<-diffpep[!is.na(diffpep$group_id),]
df_meta<-peptide_data$metadata[rownames(diffpep),]
peps_gen<-df_meta$coords

gns_ok<-unique(bins_DE$gene_id[which(bins_DE$RNA_padj<.05)])
diffp<-diffpep[diffpep$group_id%in%gns_ok,]
peps<-peps_gen[diffpep$group_id%in%gns_ok]
diffp$type<-"ns"
diffp$type[peps%over%bins_DE[which(bins_DE$RNA_padj<.05 & bins_DE$RNA_log2FC>0)]]<-"up"
diffp$type[peps%over%bins_DE[which(bins_DE$RNA_padj<.05 & bins_DE$RNA_log2FC<(0))]]<-"down"
diffp$assay<-"RNA-seq"
diffprna<-diffp

gns_ok<-unique(bins_DE$gene_id[which(bins_DE$Ribo_padj<.05)])
diffp<-diffpep[diffpep$group_id%in%gns_ok,]
peps<-peps_gen[diffpep$group_id%in%gns_ok]
diffp$type<-"ns"
diffp$type[peps%over%bins_DE[which(bins_DE$Ribo_padj<.05 & bins_DE$Ribo_log2FC>0)]]<-"up"
diffp$type[peps%over%bins_DE[which(bins_DE$Ribo_padj<.05 & bins_DE$Ribo_log2FC<(0))]]<-"down"
diffp$assay<-"Ribo-seq"
diffpribo<-diffp
df_diff<-rbind(diffprna,diffpribo)
df_diff$assay<-factor(as.character(df_diff$assay),levels=c("RNA-seq","Ribo-seq"))
comps<-list(c("down","ns"),c("ns","up"))



gg <- ggplot(df_diff, aes(x = type, y = delta_log2FC, color=type, group=type)) +
  stat_summary(fun.data = mean_se,size=1.5) +
  #stat_summary(fun.y = mean, geom="point",colour="magenta", size=2) +
  xlab("\nDEXseq significance at 14h") +
  scale_color_manual(values =alpha(c("blue","grey","red"),.8)) +
  scale_y_continuous(name="Peptide-level\nDelta log2FC") +
  theme(legend.position = "none") +
  theme_bw()+
  guides(fill=guide_legend(title="")) +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=0.5, size=16), 
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16),
        panel.grid.major = element_blank(),  # remove grid lines from panel
        panel.grid.minor = element_blank(),  # remove grid lines from panel
        strip.text=element_text(size=16))+ geom_hline(yintercept = 0,size=.5,lwd=1,lty=2) +
  stat_compare_means(label="p.format",label.y=c(-.31,-.24),method = "wilcox.test",comparisons = comps,method.args = list(alternative=c("l")),tip.length=.002)

gg_diffpep<-gg + facet_wrap(~assay,ncol = 2)
pdf(file = "Figures/Figure4_D.pdf",width = 9,height = 5)
gg_diffpep
dev.off()
#wilcox.test(df_diff$delta_log2FC[df_diff$type=="up" & df_diff$assay=="RNA-seq"],df_diff$delta_log2FC[df_diff$type=="up" & df_diff$assay=="Ribo-seq"],alternative = "l")



bins_DE<-multiDE_bins$uniq[[1]]
mcols(bins_DE)<-cbind(mcols(bins_DE),mcols(annot_files$dex_bins_todex))
bins_DE <- bins_DE[!bins_DE$gene_overlapping_noDE]
bins_DE<-bins_DE[bins_DE$type=="E"]

diffpep<-list_DE$raw$diffsplice$Peptides[["DUX4_DMSO_14h-Control_DMSO_14h"]]
diffpep<-diffpep[!is.na(diffpep$group_id),]
df_meta<-peptide_data$metadata[rownames(diffpep),]
peps_gen<-df_meta$coords

gns_ok<-unique(bins_DE$gene_id[which(bins_DE$RNA_padj<.05)])
binsok<-bins_DE[order(bins_DE$RNA_padj,decreasing = F),]
diffp<-diffpep[diffpep$group_id%in%gns_ok,]
peps<-peps_gen[diffpep$group_id%in%gns_ok]
diffp$type<-"ns"
diffp$type[peps%over%binsok[which(binsok$RNA_log2FC>0)][1:200]]<-"up"
diffp$type[peps%over%binsok[which(binsok$RNA_log2FC<(0))][1:200]]<-"down"
diffp$assay<-"RNA-seq"
diffprna<-diffp

gns_ok<-unique(bins_DE$gene_id[which(bins_DE$Ribo_padj<.05)])
diffp<-diffpep[diffpep$group_id%in%gns_ok,]
peps<-peps_gen[diffpep$group_id%in%gns_ok]
binsok<-bins_DE[order(bins_DE$Ribo_padj,decreasing = F),]

diffp$type<-"ns"
diffp$type[peps%over%binsok[which(binsok$Ribo_log2FC>0)][1:200]]<-"up"
diffp$type[peps%over%binsok[which(binsok$Ribo_log2FC<(0))][1:200]]<-"down"
diffp$assay<-"Ribo-seq"
diffpribo<-diffp
df_diff<-rbind(diffprna,diffpribo)
df_diff$assay<-factor(as.character(df_diff$assay),levels=c("RNA-seq","Ribo-seq"))
comps<-list(c("down","ns"),c("ns","up"))

gg <- ggplot(df_diff, aes(x = type, y = delta_log2FC, color=type, group=type)) +
  stat_summary(fun.data = mean_se,size=1.5) +
  #stat_summary(fun.y = mean, geom="point",colour="magenta", size=2) +
  xlab("\nDEXseq significance at 14h\ntop200 bins") +
  scale_color_manual(values =alpha(c("blue","grey","red"),.8)) +
  scale_y_continuous(name="Peptide-level\nDelta log2FC") +
  theme(legend.position = "none") +
  theme_bw()+
  guides(fill=guide_legend(title="")) +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=0.5, size=16), 
        axis.text.y=element_text(size=16),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=16),
        panel.grid.major = element_blank(),  # remove grid lines from panel
        panel.grid.minor = element_blank(),  # remove grid lines from panel
        strip.text=element_text(size=16))+ geom_hline(yintercept = 0,size=.5,lwd=1,lty=2) +
  stat_compare_means(label="p.format",label.y=c(-.31,-.14),method = "wilcox.test",comparisons = comps,method.args = list(alternative=c("l")),tip.length=.002)

gg_diffpep<-gg + facet_wrap(~assay,ncol = 2)
pdf(file = "Figures/SupplementaryFigure10.pdf",width = 9,height = 5)
gg_diffpep
dev.off()
#wilcox.test(df_diff$delta_log2FC[df_diff$type=="up" & df_diff$assay=="RNA-seq"],df_diff$delta_log2FC[df_diff$type=="up" & df_diff$assay=="Ribo-seq"],alternative = "l")


#meta<-peptide_data$metadata
#neworfs<-which(sum(meta$ORF_category_Tx_compatible=="novel",na.rm=T)==elementNROWS(meta$ORF_category_Tx_compatible))
#metanew<-meta[neworfs,]
#metanew<-metanew[sum(metanew$ORF_category_Gen=="novel")==0,]
#metanew<-metanew[which(sum(metanew$adj_unique_features_reads,na.rm = T)>10),]
#metanew<-metanew[elementNROWS(metanew$Mapped_Proteins)==1,]
#metanew$Mapped_Proteins<-unlist(metanew$Mapped_Proteins)
#metanew[which(ORFs$ORFs_tx$NMD_candidate[match(gsub(metanew$Mapped_Proteins,pattern="R2T2P_",replacement=""),ORFs$ORFs_tx$ORF_id_tr)]),]


peptide_data1<-get(load("data/SRonly_UniqueORFq_combined_database_search_peptide_results_annotated.RData"))
peptide_data2<-get(load("data/SRonly_AllORFq_combined_database_search_peptide_results_annotated.RData"))

list_peps<-list(rownames(peptide_data1$intensities_raw),rownames(peptide_data2$intensities_raw))
names(list_peps)<-c("Unique\nreads","All\nreads")
#upset(fromList(list_peps), order.by = "freq")
# pdf(file = "Figures/Fig3_ids_unqallvenn.pdf",width = 7,height = 6)
# ggVennDiagram(list_peps) + scale_fill_gradient(low="grey90",high = "darkred")
# dev.off()

df_meta<-peptide_data1$metadata
multi<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))>1,]
df_meta<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))==1,]
df_meta$gene_biotype[grepl(df_meta$gene_name,pattern="^IG") & df_meta$gene_biotype=="other"]<-"protein_coding"
df_meta$gene_biotype[df_meta$gene_biotype!="protein_coding"]<-"non-coding"
df_meta$DUX4_target<-df_meta$gene_id%in%target_genes
df_meta$gene_biotype[df_meta$DUX4_target]<-"DUX4_target"
df_meta_unq<-unique(df_meta[,c("gene_id","gene_biotype","novelty_category")])
df<-melt(table(df_meta$gene_biotype,df_meta$novelty_category))
df<-df[df$value>0,]
df$id_type<-"Peptides"
df$Var1<-factor(df$Var1,levels=names(sort(table(df_meta$gene_biotype),decreasing = T)))
df2<-melt(table(df_meta_unq$gene_biotype,df_meta_unq$novelty_category))
df2<-df2[df2$value>0,]
#df2$Var1<-factor(df2$Var1,levels=c("protein_coding","DUX4_target","non-coding"))
df2$id_type<-"Genes"
df<-rbind(df,df2)
df$Var1<-factor(df$Var1,levels=c("protein_coding","DUX4_target","non-coding"))
df_unq<-df
df_unq$mapping_reads<-"unique"


df_meta<-peptide_data2$metadata
multi<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))>1,]
df_meta<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))==1,]
df_meta$gene_biotype[grepl(df_meta$gene_name,pattern="^IG") & df_meta$gene_biotype=="other"]<-"protein_coding"
df_meta$gene_biotype[df_meta$gene_biotype!="protein_coding"]<-"non-coding"
df_meta$DUX4_target<-df_meta$gene_id%in%target_genes
df_meta$gene_biotype[df_meta$DUX4_target]<-"DUX4_target"
df_meta_unq<-unique(df_meta[,c("gene_id","gene_biotype","novelty_category")])
df<-melt(table(df_meta$gene_biotype,df_meta$novelty_category))
df<-df[df$value>0,]
df$id_type<-"Peptides"
df$Var1<-factor(df$Var1,levels=names(sort(table(df_meta$gene_biotype),decreasing = T)))
df2<-melt(table(df_meta_unq$gene_biotype,df_meta_unq$novelty_category))
df2<-df2[df2$value>0,]
#df2$Var1<-factor(df2$Var1,levels=c("protein_coding","DUX4_target","non-coding"))
df2$id_type<-"Genes"
df<-rbind(df,df2)
df$Var1<-factor(df$Var1,levels=c("protein_coding","DUX4_target","non-coding"))
df_all<-df
df_all$mapping_reads<-"all"
df<-rbind(df_unq,df_all)

gg_idsreads <- ggplot(data = df, aes(x=Var2, y=value+1, fill=mapping_reads,label=value)) + 
  #geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
  geom_bar(stat="identity", position=position_dodge(width = .9)) + 
  #geom_col(stat="identity", position=position_dodge(width = .9),preserve = c("total", "single")) + 
  
  labs(x="", y="# identifications") +
  scale_fill_manual(values = c("grey78","grey9"),"Ribo-seq mapping") +
  theme_bw() +
  scale_y_log10() +
  theme(plot.margin = ggplot2::margin(.5, .5, .5, .5, "cm"),
        panel.background=element_blank(), 
        axis.ticks=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=14), 
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20),
        #legend.position="bottom",
        legend.text=element_text(size=18),
        legend.title=element_text(size=20),
        panel.grid.major.y = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.x = element_line(linewidth=.1, color="dark grey" ,linetype=2)) +
  geom_text(position=position_dodge(width = .9),vjust = -0.25) + facet_wrap(id_type~Var1,nrow = 2,scales = "free")

# pdf(file = "Figures/Fig3_ids_mapping.pdf",width = 11,height = 6)
# gg_idsreads
# dev.off()



peptide_data1<-get(load("data/SRonly_UniqueORFq_combined_database_search_peptide_results_annotated.RData"))
peptide_data2<-get(load("data/LRalso_UniqueORFq_combined_database_search_peptide_results_annotated.RData"))

list_peps<-list(rownames(peptide_data1$intensities_raw),rownames(peptide_data2$intensities_raw))
names(list_peps)<-c("R1","R1R2")
#upset(fromList(list_peps), order.by = "freq")
# pdf(file = "Figures/Fig5_ids_r1r2.pdf",width = 6,height = 6)
# ggVennDiagram(list_peps) + scale_fill_gradient(low="grey90",high = "darkred")
# dev.off()

df_meta<-peptide_data1$metadata
multi<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))>1,]
df_meta<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))==1,]
df_meta$gene_biotype[grepl(df_meta$gene_name,pattern="^IG") & df_meta$gene_biotype=="other"]<-"protein_coding"
df_meta$gene_biotype[df_meta$gene_biotype!="protein_coding"]<-"non-coding"
df_meta$DUX4_target<-df_meta$gene_id%in%target_genes
df_meta$gene_biotype[df_meta$DUX4_target]<-"DUX4_target"
df_meta_unq<-unique(df_meta[,c("gene_id","gene_biotype","novelty_category")])
df<-melt(table(df_meta$gene_biotype,df_meta$novelty_category))
df<-df[df$value>0,]
df$id_type<-"Peptides"
df$Var1<-factor(df$Var1,levels=names(sort(table(df_meta$gene_biotype),decreasing = T)))
df2<-melt(table(df_meta_unq$gene_biotype,df_meta_unq$novelty_category))
df2<-df2[df2$value>0,]
#df2$Var1<-factor(df2$Var1,levels=c("protein_coding","DUX4_target","non-coding"))
df2$id_type<-"Genes"
df<-rbind(df,df2)
df$Var1<-factor(df$Var1,levels=c("protein_coding","DUX4_target","non-coding"))
df_r1<-df
df_r1$R1R2<-"R1"

df_meta<-peptide_data2$metadata
multi<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))>1,]
df_meta<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))==1,]
df_meta$gene_biotype[grepl(df_meta$gene_name,pattern="^IG") & df_meta$gene_biotype=="other"]<-"protein_coding"
df_meta$gene_biotype[df_meta$gene_biotype!="protein_coding"]<-"non-coding"
df_meta$DUX4_target<-df_meta$gene_id%in%target_genes
df_meta$gene_biotype[df_meta$DUX4_target]<-"DUX4_target"
df_meta_unq<-unique(df_meta[,c("gene_id","gene_biotype","novelty_category")])
df<-melt(table(df_meta$gene_biotype,df_meta$novelty_category))
df<-df[df$value>0,]
df$id_type<-"Peptides"
df$Var1<-factor(df$Var1,levels=names(sort(table(df_meta$gene_biotype),decreasing = T)))
df2<-melt(table(df_meta_unq$gene_biotype,df_meta_unq$novelty_category))
df2<-df2[df2$value>0,]
#df2$Var1<-factor(df2$Var1,levels=c("protein_coding","DUX4_target","non-coding"))
df2$id_type<-"Genes"
df<-rbind(df,df2)
df$Var1<-factor(df$Var1,levels=c("protein_coding","DUX4_target","non-coding"))
df_r1r2<-df
df_r1r2$R1R2<-"R1R2"
df<-rbind(df_r1,df_r1r2)

gg_idsreads <- ggplot(data = df, aes(x=Var2, y=value+1, fill=R1R2,label=value)) + 
  #geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
  geom_bar(stat="identity", position=position_dodge(width = .9)) + 
  #geom_col(stat="identity", position=position_dodge(width = .9),preserve = c("total", "single")) + 
  
  labs(x="", y="# identifications") +
  scale_fill_manual(values = c("grey78","grey9"),"Transcriptome") +
  theme_bw() +
  scale_y_log10() +
  theme(plot.margin = ggplot2::margin(.5, .5, .5, .5, "cm"),
        panel.background=element_blank(), 
        axis.ticks=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=14), 
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20),
        #legend.position="bottom",
        legend.text=element_text(size=18),
        legend.title=element_text(size=20),
        panel.grid.major.y = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.x = element_line( size=.1, color="dark grey" ,linetype=2)) +
  geom_text(position=position_dodge(width = .9),vjust = -0.1) + facet_wrap(id_type~Var1,nrow = 2,scales = "fixed")

# pdf(file = "Figures/Fig5_ids_mapping.pdf",width = 11,height = 6)
# gg_idsreads
# dev.off()



peptide_data1<-get(load("data/SRonly_AllORFq_combined_database_search_peptide_results_annotated.RData"))
peptide_data2<-get(load("data/LRalso_AllORFq_combined_database_search_peptide_results_annotated.RData"))

df_meta<-peptide_data1$metadata
multi<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))>1,]
df_meta<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))==1,]
df_meta$gene_biotype[grepl(df_meta$gene_name,pattern="^IG") & df_meta$gene_biotype=="other"]<-"protein_coding"
df_meta$gene_biotype[df_meta$gene_biotype!="protein_coding"]<-"non-coding"
df_meta$DUX4_target<-df_meta$gene_id%in%target_genes
df_meta$gene_biotype[df_meta$DUX4_target]<-"DUX4_target"
df_meta_unq<-unique(df_meta[,c("gene_id","gene_biotype","novelty_category")])
df<-melt(table(df_meta$gene_biotype,df_meta$novelty_category))
df<-df[df$value>0,]
df$id_type<-"Peptides"
df$Var1<-factor(df$Var1,levels=names(sort(table(df_meta$gene_biotype),decreasing = T)))
df2<-melt(table(df_meta_unq$gene_biotype,df_meta_unq$novelty_category))
df2<-df2[df2$value>0,]
#df2$Var1<-factor(df2$Var1,levels=c("protein_coding","DUX4_target","non-coding"))
df2$id_type<-"Genes"
df<-rbind(df,df2)
df$Var1<-factor(df$Var1,levels=c("protein_coding","DUX4_target","non-coding"))
df_r1<-df
df_r1$R1R2<-"R1"

df_meta<-peptide_data2$metadata
multi<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))>1,]
df_meta<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))==1,]
df_meta$gene_biotype[grepl(df_meta$gene_name,pattern="^IG") & df_meta$gene_biotype=="other"]<-"protein_coding"
df_meta$gene_biotype[df_meta$gene_biotype!="protein_coding"]<-"non-coding"
df_meta$DUX4_target<-df_meta$gene_id%in%target_genes
df_meta$gene_biotype[df_meta$DUX4_target]<-"DUX4_target"
df_meta_unq<-unique(df_meta[,c("gene_id","gene_biotype","novelty_category")])
df<-melt(table(df_meta$gene_biotype,df_meta$novelty_category))
df<-df[df$value>0,]
df$id_type<-"Peptides"
df$Var1<-factor(df$Var1,levels=names(sort(table(df_meta$gene_biotype),decreasing = T)))
df2<-melt(table(df_meta_unq$gene_biotype,df_meta_unq$novelty_category))
df2<-df2[df2$value>0,]
#df2$Var1<-factor(df2$Var1,levels=c("protein_coding","DUX4_target","non-coding"))
df2$id_type<-"Genes"
df<-rbind(df,df2)
df$Var1<-factor(df$Var1,levels=c("protein_coding","DUX4_target","non-coding"))
df_r1r2<-df
df_r1r2$R1R2<-"R1R2"
df<-rbind(df_r1,df_r1r2)

gg_idsreads <- ggplot(data = df, aes(x=Var2, y=value, fill=R1R2,label=value)) + 
  #geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
  geom_bar(stat="identity", position=position_dodge(width = .9)) + 
  #geom_col(stat="identity", position=position_dodge(width = .9),preserve = c("total", "single")) + 
  
  labs(x="", y="# identifications") +
  scale_fill_manual(values = c("grey78","grey9"),"Transcriptome") +
  theme_bw() +
  scale_y_log10() +
  theme(plot.margin = ggplot2::margin(.5, .5, .5, .5, "cm"),
        panel.background=element_blank(), 
        axis.ticks=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=14), 
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20),
        #legend.position="bottom",
        legend.text=element_text(size=18),
        legend.title=element_text(size=20),
        panel.grid.major.y = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.x = element_line( size=.1, color="dark grey" ,linetype=2)) +
  geom_text(position=position_dodge(width = .9),vjust = -0.1) + facet_wrap(id_type~Var1,nrow = 2,scales = "fixed")

# pdf(file = "Figures/Fig5_ids_mapping_all.pdf",width = 11,height = 6)
# gg_idsreads
# dev.off()


### Here DE

list_DE<-get(load("data/SRonly_UniqueORFq_combined_database_search_proteomics_multiDE.RData"))
peptide_data<-get(load("data/SRonly_UniqueORFq_combined_database_search_peptide_results_annotated.RData"))

DE_res<-list_DE$raw$DE

df_meta<-peptide_data$metadata
df_meta<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))==1,]
df_meta$gene_biotype[grepl(df_meta$gene_name,pattern="^IG") & df_meta$gene_biotype=="other"]<-"protein_coding"
df_meta$gene_biotype[df_meta$gene_biotype!="protein_coding"]<-"non-coding"
df_meta$DUX4_target<-df_meta$gene_id%in%target_genes
df_meta$gene_biotype[df_meta$DUX4_target]<-"DUX4_target"


# df<-DE_res$gene_id_novelty[["DUX4_DMSO_04h-Control_DMSO_04h"]]
# 
# gns<-sapply(strsplit(rownames(df),";"),"[[",1)
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# #df$novelty<-df_meta$novelty_category[match(gns,df_meta$gene_id)]
# df$novelty<-sapply(strsplit(rownames(df),";"),"[[",2)
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Genes"
# df$time<-"04h"
# df1<-df
# 
# df<-DE_res$gene_id_novelty[["DUX4_DMSO_08h-Control_DMSO_08h"]]
# 
# gns<-sapply(strsplit(rownames(df),";"),"[[",1)
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# #df$novelty<-df_meta$novelty_category[match(gns,df_meta$gene_id)]
# df$novelty<-sapply(strsplit(rownames(df),";"),"[[",2)
# 
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Genes"
# df$time<-"08h"
# df2<-df
# 
df<-DE_res$gene_id_novelty[["DUX4_DMSO_14h-Control_DMSO_14h"]]

gns<-sapply(strsplit(rownames(df),";"),"[[",1)
df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
df$gene_id<-gns
df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
df$novelty<-sapply(strsplit(rownames(df),";"),"[[",2)
df$category<-paste(df$gene_biotype,df$novelty,sep=";")
df$sig<-F
df$sig[df$padj<.05]<-T
df<-df[order(df$padj,decreasing = F),]
df<-df[which(df$padj<.05),]
df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
df$id_type="Genes"
df$time<-"14h"
df3<-df

# dfg<-rbind(df1,df2,df3)
dfg <- df3

# df<-DE_res$Peptides[["DUX4_DMSO_04h-Control_DMSO_04h"]]
# gns<-df_meta$gene_id
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# df$novelty<-df_meta$novelty_category
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Peptides"
# df$time<-"04h"
# df1<-df
# 
# df<-DE_res$Peptides[["DUX4_DMSO_08h-Control_DMSO_08h"]]
# 
# gns<-df_meta$gene_id
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# df$novelty<-df_meta$novelty_category
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Peptides"
# df$time<-"08h"
# df2<-df
# 
df<-DE_res$Peptides[["DUX4_DMSO_14h-Control_DMSO_14h"]]

gns<-df_meta$gene_id
df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
df$gene_id<-gns
df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
df$novelty<-df_meta$novelty_category
df$category<-paste(df$gene_biotype,df$novelty,sep=";")
df$sig<-F
df$sig[df$padj<.05]<-T
df<-df[order(df$padj,decreasing = F),]
df<-df[which(df$padj<.05),]
df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
df$id_type="Peptides"
df$time<-"14h"
df3<-df
# dfp<-rbind(df1,df2,df3)
dfp <- df3

df<-rbind(dfg,dfp)

df<-df[df$Var3,]
 
df<-df[df$time!="04h",]
#df<-df[df$value>0,]
df$value[df$value==0]<-NA

df$Var1<-factor(as.character(df$Var1),levels=c("protein_coding","DUX4_target","non-coding"))
df$search<-"unique"
df_unq<-df


list_DE<-get(load("data/SRonly_AllORFq_combined_database_search_proteomics_multiDE.RData"))
peptide_data<-get(load("data/SRonly_AllORFq_combined_database_search_peptide_results_annotated.RData"))

DE_res<-list_DE$raw$DE

df_meta<-peptide_data$metadata
df_meta<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))==1,]
df_meta$gene_biotype[grepl(df_meta$gene_name,pattern="^IG") & df_meta$gene_biotype=="other"]<-"protein_coding"
df_meta$gene_biotype[df_meta$gene_biotype!="protein_coding"]<-"non-coding"
df_meta$DUX4_target<-df_meta$gene_id%in%target_genes
df_meta$gene_biotype[df_meta$DUX4_target]<-"DUX4_target"


# df<-DE_res$gene_id_novelty[["DUX4_DMSO_04h-Control_DMSO_04h"]]
# 
# gns<-sapply(strsplit(rownames(df),";"),"[[",1)
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# #df$novelty<-df_meta$novelty_category[match(gns,df_meta$gene_id)]
# df$novelty<-sapply(strsplit(rownames(df),";"),"[[",2)
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Genes"
# df$time<-"04h"
# df1<-df
# 
# df<-DE_res$gene_id_novelty[["DUX4_DMSO_08h-Control_DMSO_08h"]]
# 
# gns<-sapply(strsplit(rownames(df),";"),"[[",1)
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# #df$novelty<-df_meta$novelty_category[match(gns,df_meta$gene_id)]
# df$novelty<-sapply(strsplit(rownames(df),";"),"[[",2)
# 
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Genes"
# df$time<-"08h"
# df2<-df

df<-DE_res$gene_id_novelty[["DUX4_DMSO_14h-Control_DMSO_14h"]]

gns<-sapply(strsplit(rownames(df),";"),"[[",1)
df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
df$gene_id<-gns
df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
df$novelty<-sapply(strsplit(rownames(df),";"),"[[",2)
df$category<-paste(df$gene_biotype,df$novelty,sep=";")
df$sig<-F
df$sig[df$padj<.05]<-T
df<-df[order(df$padj,decreasing = F),]
df<-df[which(df$padj<.05),]
df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
df$id_type="Genes"
df$time<-"14h"
df3<-df

# dfg<-rbind(df1,df2,df3)
dfg<-df3

# df<-DE_res$Peptides[["DUX4_DMSO_04h-Control_DMSO_04h"]]
# gns<-df_meta$gene_id
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# df$novelty<-df_meta$novelty_category
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Peptides"
# df$time<-"04h"
# df1<-df
# 
# df<-DE_res$Peptides[["DUX4_DMSO_08h-Control_DMSO_08h"]]
# 
# gns<-df_meta$gene_id
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# df$novelty<-df_meta$novelty_category
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Peptides"
# df$time<-"08h"
# df2<-df

df<-DE_res$Peptides[["DUX4_DMSO_14h-Control_DMSO_14h"]]

gns<-df_meta$gene_id
df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
df$gene_id<-gns
df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
df$novelty<-df_meta$novelty_category
df$category<-paste(df$gene_biotype,df$novelty,sep=";")
df$sig<-F
df$sig[df$padj<.05]<-T
df<-df[order(df$padj,decreasing = F),]
df<-df[which(df$padj<.05),]
df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
df$id_type="Peptides"
df$time<-"14h"
df3<-df
# dfp<-rbind(df1,df2,df3)
dfp <- df3

df<-rbind(dfg,dfp)

df<-df[df$Var3,]

df<-df[df$time!="04h",]
#df<-df[df$value>0,]
df$value[df$value==0]<-NA

df$Var1<-factor(as.character(df$Var1),levels=c("protein_coding","DUX4_target","non-coding"))
df$search<-"all"
df_all<-df
df<-rbind(df_unq,df_all)

dfok<-df[df$time=="14h",]
gg_deres_unqall <- ggplot(data = dfok, aes(x=Var2, y=value+1, fill=search,label=value)) + 
  #geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
  geom_bar(stat="identity", position=position_dodge(width = .9)) + 
  labs(x="", y="# significantly upregulated (14h)") +
  scale_fill_manual(values = c("grey78","grey9"),"Ribo-seq mapping") +
  theme_bw() +
  scale_y_log10() +
  theme(plot.margin = ggplot2::margin(.5, .5, .5, .5, "cm"),
        panel.background=element_blank(), 
        axis.ticks=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=14), 
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20),
        legend.text=element_text(size=14),
        panel.grid.major.y = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.x = element_line( size=.1, color="dark grey" ,linetype=2)) +
  geom_text(position=position_dodge(width = .9),vjust = -0.25) + facet_wrap(id_type~Var1,nrow = 2,scales = "fixed")

# pdf(file = "Figures/Fig3_deres_mapping.pdf",width = 11,height = 6)
# gg_deres_unqall
# dev.off()


list_DE<-get(load("data/SRonly_UniqueORFq_combined_database_search_proteomics_multiDE.RData"))
peptide_data<-get(load("data/SRonly_UniqueORFq_combined_database_search_peptide_results_annotated.RData"))

DE_res<-list_DE$raw$DE

df_meta<-peptide_data$metadata
df_meta<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))==1,]
df_meta$gene_biotype[grepl(df_meta$gene_name,pattern="^IG") & df_meta$gene_biotype=="other"]<-"protein_coding"
df_meta$gene_biotype[df_meta$gene_biotype!="protein_coding"]<-"non-coding"
df_meta$DUX4_target<-df_meta$gene_id%in%target_genes
df_meta$gene_biotype[df_meta$DUX4_target]<-"DUX4_target"


# df<-DE_res$gene_id_novelty[["DUX4_DMSO_04h-Control_DMSO_04h"]]
# 
# gns<-sapply(strsplit(rownames(df),";"),"[[",1)
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# #df$novelty<-df_meta$novelty_category[match(gns,df_meta$gene_id)]
# df$novelty<-sapply(strsplit(rownames(df),";"),"[[",2)
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Genes"
# df$time<-"04h"
# df1<-df
# 
# df<-DE_res$gene_id_novelty[["DUX4_DMSO_08h-Control_DMSO_08h"]]
# 
# gns<-sapply(strsplit(rownames(df),";"),"[[",1)
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# #df$novelty<-df_meta$novelty_category[match(gns,df_meta$gene_id)]
# df$novelty<-sapply(strsplit(rownames(df),";"),"[[",2)
# 
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Genes"
# df$time<-"08h"
# df2<-df

df<-DE_res$gene_id_novelty[["DUX4_DMSO_14h-Control_DMSO_14h"]]

gns<-sapply(strsplit(rownames(df),";"),"[[",1)
df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
df$gene_id<-gns
df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
df$novelty<-sapply(strsplit(rownames(df),";"),"[[",2)
df$category<-paste(df$gene_biotype,df$novelty,sep=";")
df$sig<-F
df$sig[df$padj<.05]<-T
df<-df[order(df$padj,decreasing = F),]
df<-df[which(df$padj<.05),]
df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
df$id_type="Genes"
df$time<-"14h"
df3<-df

# dfg<-rbind(df1,df2,df3)
dfg <- df3

# df<-DE_res$Peptides[["DUX4_DMSO_04h-Control_DMSO_04h"]]
# gns<-df_meta$gene_id
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# df$novelty<-df_meta$novelty_category
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Peptides"
# df$time<-"04h"
# df1<-df
# 
# df<-DE_res$Peptides[["DUX4_DMSO_08h-Control_DMSO_08h"]]
# 
# gns<-df_meta$gene_id
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# df$novelty<-df_meta$novelty_category
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Peptides"
# df$time<-"08h"
# df2<-df

df<-DE_res$Peptides[["DUX4_DMSO_14h-Control_DMSO_14h"]]

gns<-df_meta$gene_id
df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
df$gene_id<-gns
df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
df$novelty<-df_meta$novelty_category
df$category<-paste(df$gene_biotype,df$novelty,sep=";")
df$sig<-F
df$sig[df$padj<.05]<-T
df<-df[order(df$padj,decreasing = F),]
df<-df[which(df$padj<.05),]
df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
df$id_type="Peptides"
df$time<-"14h"
df3<-df
# dfp<-rbind(df1,df2,df3)
dfp <- df3

df<-rbind(dfg,dfp)

df<-df[df$Var3,]

df<-df[df$time!="04h",]
#df<-df[df$value>0,]
df$value[df$value==0]<-NA

df$Var1<-factor(as.character(df$Var1),levels=c("protein_coding","DUX4_target","non-coding"))
df_r1<-df
df_r1$R1R2<-"R1"


list_DE<-get(load("data/LRalso_UniqueORFq_combined_database_search_proteomics_multiDE.RData"))
peptide_data<-get(load("data/LRalso_UniqueORFq_combined_database_search_peptide_results_annotated.RData"))


DE_res<-list_DE$raw$DE

df_meta<-peptide_data$metadata
df_meta<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))==1,]
df_meta$gene_biotype[grepl(df_meta$gene_name,pattern="^IG") & df_meta$gene_biotype=="other"]<-"protein_coding"
df_meta$gene_biotype[df_meta$gene_biotype!="protein_coding"]<-"non-coding"
df_meta$DUX4_target<-df_meta$gene_id%in%target_genes
df_meta$gene_biotype[df_meta$DUX4_target]<-"DUX4_target"


# df<-DE_res$gene_id_novelty[["DUX4_DMSO_04h-Control_DMSO_04h"]]
# 
# gns<-sapply(strsplit(rownames(df),";"),"[[",1)
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# #df$novelty<-df_meta$novelty_category[match(gns,df_meta$gene_id)]
# df$novelty<-sapply(strsplit(rownames(df),";"),"[[",2)
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Genes"
# df$time<-"04h"
# df1<-df
# 
# df<-DE_res$gene_id_novelty[["DUX4_DMSO_08h-Control_DMSO_08h"]]
# 
# gns<-sapply(strsplit(rownames(df),";"),"[[",1)
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# #df$novelty<-df_meta$novelty_category[match(gns,df_meta$gene_id)]
# df$novelty<-sapply(strsplit(rownames(df),";"),"[[",2)
# 
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Genes"
# df$time<-"08h"
# df2<-df

df<-DE_res$gene_id_novelty[["DUX4_DMSO_14h-Control_DMSO_14h"]]

gns<-sapply(strsplit(rownames(df),";"),"[[",1)
df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
df$gene_id<-gns
df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
df$novelty<-sapply(strsplit(rownames(df),";"),"[[",2)
df$category<-paste(df$gene_biotype,df$novelty,sep=";")
df$sig<-F
df$sig[df$padj<.05]<-T
df<-df[order(df$padj,decreasing = F),]
df<-df[which(df$padj<.05),]
df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
df$id_type="Genes"
df$time<-"14h"
df3<-df

# dfg<-rbind(df1,df2,df3)
dfg <- df3

# df<-DE_res$Peptides[["DUX4_DMSO_04h-Control_DMSO_04h"]]
# gns<-df_meta$gene_id
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# df$novelty<-df_meta$novelty_category
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Peptides"
# df$time<-"04h"
# df1<-df
# 
# df<-DE_res$Peptides[["DUX4_DMSO_08h-Control_DMSO_08h"]]
# 
# gns<-df_meta$gene_id
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# df$novelty<-df_meta$novelty_category
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Peptides"
# df$time<-"08h"
# df2<-df

df<-DE_res$Peptides[["DUX4_DMSO_14h-Control_DMSO_14h"]]

gns<-df_meta$gene_id
df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
df$gene_id<-gns
df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
df$novelty<-df_meta$novelty_category
df$category<-paste(df$gene_biotype,df$novelty,sep=";")
df$sig<-F
df$sig[df$padj<.05]<-T
df<-df[order(df$padj,decreasing = F),]
df<-df[which(df$padj<.05),]
df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
df$id_type="Peptides"
df$time<-"14h"
df3<-df
# dfp<-rbind(df1,df2,df3)
dfp <- df3

df<-rbind(dfg,dfp)

df<-df[df$Var3,]

df<-df[df$time!="04h",]
#df<-df[df$value>0,]
df$value[df$value==0]<-NA

df$Var1<-factor(as.character(df$Var1),levels=c("protein_coding","DUX4_target","non-coding"))
df_r1r2<-df
df_r1r2$R1R2<-"R1R2"

df<-rbind(df_r1,df_r1r2)

dfok<-df[df$time=="14h",]
gg_deres_r1r2 <- ggplot(data = dfok, aes(x=Var2, y=value+1, fill=R1R2,label=value)) + 
  #geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
  geom_bar(stat="identity", position=position_dodge(width = .9)) + 
  labs(x="", y="# significantly upregulated (14h)") +
  scale_fill_manual(values = c("grey78","grey9"),"Transcriptome") +
  theme_bw() +
  scale_y_log10() +
  theme(plot.margin = ggplot2::margin(.5, .5, .5, .5, "cm"),
        panel.background=element_blank(), 
        axis.ticks=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=14), 
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20),
        legend.text=element_text(size=14),
        panel.grid.major.y = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.x = element_line( size=.1, color="dark grey" ,linetype=2)) +
  geom_text(position=position_dodge(width = .9),vjust = -0.25) + facet_wrap(id_type~Var1,nrow = 2,scales = "fixed")

# pdf(file = "Figures/Fig5_deres_mapping.pdf",width = 11,height = 6)
# gg_deres_r1r2
# dev.off()


list_DE<-get(load("data/SRonly_AllORFq_combined_database_search_proteomics_multiDE.RData"))
peptide_data<-get(load("data/SRonly_AllORFq_combined_database_search_peptide_results_annotated.RData"))

DE_res<-list_DE$raw$DE

df_meta<-peptide_data$metadata
df_meta<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))==1,]
df_meta$gene_biotype[grepl(df_meta$gene_name,pattern="^IG") & df_meta$gene_biotype=="other"]<-"protein_coding"
df_meta$gene_biotype[df_meta$gene_biotype!="protein_coding"]<-"non-coding"
df_meta$DUX4_target<-df_meta$gene_id%in%target_genes
df_meta$gene_biotype[df_meta$DUX4_target]<-"DUX4_target"


# df<-DE_res$gene_id_novelty[["DUX4_DMSO_04h-Control_DMSO_04h"]]
# 
# gns<-sapply(strsplit(rownames(df),";"),"[[",1)
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# #df$novelty<-df_meta$novelty_category[match(gns,df_meta$gene_id)]
# df$novelty<-sapply(strsplit(rownames(df),";"),"[[",2)
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Genes"
# df$time<-"04h"
# df1<-df
# 
# df<-DE_res$gene_id_novelty[["DUX4_DMSO_08h-Control_DMSO_08h"]]
# 
# gns<-sapply(strsplit(rownames(df),";"),"[[",1)
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# #df$novelty<-df_meta$novelty_category[match(gns,df_meta$gene_id)]
# df$novelty<-sapply(strsplit(rownames(df),";"),"[[",2)
# 
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Genes"
# df$time<-"08h"
# df2<-df

df<-DE_res$gene_id_novelty[["DUX4_DMSO_14h-Control_DMSO_14h"]]

gns<-sapply(strsplit(rownames(df),";"),"[[",1)
df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
df$gene_id<-gns
df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
df$novelty<-sapply(strsplit(rownames(df),";"),"[[",2)
df$category<-paste(df$gene_biotype,df$novelty,sep=";")
df$sig<-F
df$sig[df$padj<.05]<-T
df<-df[order(df$padj,decreasing = F),]
df<-df[which(df$padj<.05),]
df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
df$id_type="Genes"
df$time<-"14h"
df3<-df

# dfg<-rbind(df1,df2,df3)
dfg <- df3

# df<-DE_res$Peptides[["DUX4_DMSO_04h-Control_DMSO_04h"]]
# gns<-df_meta$gene_id
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# df$novelty<-df_meta$novelty_category
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Peptides"
# df$time<-"04h"
# df1<-df
# 
# df<-DE_res$Peptides[["DUX4_DMSO_08h-Control_DMSO_08h"]]
# 
# gns<-df_meta$gene_id
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# df$novelty<-df_meta$novelty_category
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Peptides"
# df$time<-"08h"
# df2<-df

df<-DE_res$Peptides[["DUX4_DMSO_14h-Control_DMSO_14h"]]

gns<-df_meta$gene_id
df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
df$gene_id<-gns
df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
df$novelty<-df_meta$novelty_category
df$category<-paste(df$gene_biotype,df$novelty,sep=";")
df$sig<-F
df$sig[df$padj<.05]<-T
df<-df[order(df$padj,decreasing = F),]
df<-df[which(df$padj<.05),]
df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
df$id_type="Peptides"
df$time<-"14h"
df3<-df
# dfp<-rbind(df1,df2,df3)
dfp <- df3

df<-rbind(dfg,dfp)

df<-df[df$Var3,]

df<-df[df$time!="04h",]
#df<-df[df$value>0,]
df$value[df$value==0]<-NA

df$Var1<-factor(as.character(df$Var1),levels=c("protein_coding","DUX4_target","non-coding"))
df_r1<-df
df_r1$R1R2<-"R1"


list_DE<-get(load("data/LRalso_AllORFq_combined_database_search_proteomics_multiDE.RData"))
peptide_data<-get(load("data/LRalso_AllORFq_combined_database_search_peptide_results_annotated.RData"))

DE_res<-list_DE$raw$DE

df_meta<-peptide_data$metadata
df_meta<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))==1,]
df_meta$gene_biotype[grepl(df_meta$gene_name,pattern="^IG") & df_meta$gene_biotype=="other"]<-"protein_coding"
df_meta$gene_biotype[df_meta$gene_biotype!="protein_coding"]<-"non-coding"
df_meta$DUX4_target<-df_meta$gene_id%in%target_genes
df_meta$gene_biotype[df_meta$DUX4_target]<-"DUX4_target"


# df<-DE_res$gene_id_novelty[["DUX4_DMSO_04h-Control_DMSO_04h"]]
# 
# gns<-sapply(strsplit(rownames(df),";"),"[[",1)
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# #df$novelty<-df_meta$novelty_category[match(gns,df_meta$gene_id)]
# df$novelty<-sapply(strsplit(rownames(df),";"),"[[",2)
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Genes"
# df$time<-"04h"
# df1<-df
# 
# df<-DE_res$gene_id_novelty[["DUX4_DMSO_08h-Control_DMSO_08h"]]
# 
# gns<-sapply(strsplit(rownames(df),";"),"[[",1)
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# #df$novelty<-df_meta$novelty_category[match(gns,df_meta$gene_id)]
# df$novelty<-sapply(strsplit(rownames(df),";"),"[[",2)
# 
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Genes"
# df$time<-"08h"
# df2<-df

df<-DE_res$gene_id_novelty[["DUX4_DMSO_14h-Control_DMSO_14h"]]

gns<-sapply(strsplit(rownames(df),";"),"[[",1)
df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
df$gene_id<-gns
df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
df$novelty<-sapply(strsplit(rownames(df),";"),"[[",2)
df$category<-paste(df$gene_biotype,df$novelty,sep=";")
df$sig<-F
df$sig[df$padj<.05]<-T
df<-df[order(df$padj,decreasing = F),]
df<-df[which(df$padj<.05),]
df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
df$id_type="Genes"
df$time<-"14h"
df3<-df

# dfg<-rbind(df1,df2,df3)
dfg <- df3

# df<-DE_res$Peptides[["DUX4_DMSO_04h-Control_DMSO_04h"]]
# gns<-df_meta$gene_id
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# df$novelty<-df_meta$novelty_category
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Peptides"
# df$time<-"04h"
# df1<-df
# 
# df<-DE_res$Peptides[["DUX4_DMSO_08h-Control_DMSO_08h"]]
# 
# gns<-df_meta$gene_id
# df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
# df$gene_id<-gns
# df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
# df$novelty<-df_meta$novelty_category
# df$category<-paste(df$gene_biotype,df$novelty,sep=";")
# df$sig<-F
# df$sig[df$padj<.05]<-T
# df<-df[order(df$padj,decreasing = F),]
# df<-df[which(df$padj<.05),]
# df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
# df$id_type="Peptides"
# df$time<-"08h"
# df2<-df

df<-DE_res$Peptides[["DUX4_DMSO_14h-Control_DMSO_14h"]]

gns<-df_meta$gene_id
df$gene_name<-df_meta$gene_name[match(gns,df_meta$gene_id)]
df$gene_id<-gns
df$gene_biotype<-df_meta$gene_biotype[match(gns,df_meta$gene_id)]
df$novelty<-df_meta$novelty_category
df$category<-paste(df$gene_biotype,df$novelty,sep=";")
df$sig<-F
df$sig[df$padj<.05]<-T
df<-df[order(df$padj,decreasing = F),]
df<-df[which(df$padj<.05),]
df<-melt(table(df$gene_biotype,df$novelty,df$logFC>0))
df$id_type="Peptides"
df$time<-"14h"
df3<-df
# dfp<-rbind(df1,df2,df3)
dfp <- df3

df<-rbind(dfg,dfp)

df<-df[df$Var3,]

df<-df[df$time!="04h",]
#df<-df[df$value>0,]
df$value[df$value==0]<-NA

df$Var1<-factor(as.character(df$Var1),levels=c("protein_coding","DUX4_target","non-coding"))
df_r1r2<-df
df_r1r2$R1R2<-"R1R2"

df<-rbind(df_r1,df_r1r2)

dfok<-df[df$time=="14h",]
gg_deres_r1r2 <- ggplot(data = dfok, aes(x=Var2, y=value+1, fill=R1R2,label=value)) + 
  #geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
  geom_bar(stat="identity", position=position_dodge(width = .9)) + 
  labs(x="", y="# significantly upregulated (14h)") +
  scale_fill_manual(values = c("grey78","grey9"),"Transcriptome") +
  theme_bw() +
  scale_y_log10() +
  theme(plot.margin = ggplot2::margin(.5, .5, .5, .5, "cm"),
        panel.background=element_blank(), 
        axis.ticks=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=14), 
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20),
        legend.text=element_text(size=14),
        panel.grid.major.y = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.x = element_line( size=.1, color="dark grey" ,linetype=2)) +
  geom_text(position=position_dodge(width = .9),vjust = -0.25) + facet_wrap(id_type~Var1,nrow = 2,scales = "fixed")

# pdf(file = "Figures/Fig5_deres_mapping_all.pdf",width = 11,height = 6)
# gg_deres_r1r2
# dev.off()


#diff db


peptide_data1<-get(load("data/SRonly_UniqueORFq_custom_database_search_peptide_results_annotated.RData"))
peptide_data2<-get(load("data/SRonly_UniqueORFq_ref_database_search_peptide_results_annotated.RData"))
peptide_data3<-get(load("data/SRonly_UniqueORFq_combined_database_search_peptide_results_annotated.RData"))


list_peps<-list(rownames(peptide_data1$intensities_raw),rownames(peptide_data2$intensities_raw),rownames(peptide_data3$intensities_raw))
names(list_peps)<-c("ORFquant","GENCODE","combined")
#upset(fromList(list_peps), order.by = "freq")
pdf(file = "Figures/SupplementaryFigure8_left.pdf",width = 6,height = 6)
ggVennDiagram(list_peps) + scale_fill_gradient(low="grey90",high = "darkred")
dev.off()

list_gns<-list(unique(peptide_data1$metadata$gene_id[elementNROWS(peptide_data1$metadata$Mapped_Genes)==1]),
               unique(peptide_data2$metadata$gene_id[elementNROWS(peptide_data2$metadata$Mapped_Genes)==1]),
               unique(peptide_data3$metadata$gene_id[elementNROWS(peptide_data3$metadata$Mapped_Genes)==1]))
names(list_gns)<-c("ORFquant","GENCODE","combined")
#upset(fromList(list_peps), order.by = "freq")
pdf(file = "Figures/SupplementaryFigure8_right.pdf",width = 6,height = 6)
ggVennDiagram(list_gns) + scale_fill_gradient(low="grey90",high = "darkred")
dev.off()



#TODO: summarize findings by ORFs (e.g. uORFs etc.. and calculate if they are DE in the ORFquant-only search)

list_DE1<-get(load("data/SRonly_UniqueORFq_ref_database_search_proteomics_multiDE.RData"))
list_DE2<-get(load("data/SRonly_UniqueORFq_custom_database_search_proteomics_multiDE.RData"))
list_DE3<-get(load("data/SRonly_UniqueORFq_combined_database_search_proteomics_multiDE.RData"))

list_DE1ok<-list_DE1$raw$DE$gene_id$`DUX4_DMSO_14h-Control_DMSO_14h`
list_DE2ok<-list_DE2$raw$DE$gene_id$`DUX4_DMSO_14h-Control_DMSO_14h`
list_DE3ok<-list_DE3$raw$DE$gene_id$`DUX4_DMSO_14h-Control_DMSO_14h`
gns1<-sapply(strsplit(rownames(list_DE3ok[which(list_DE3ok$padj<.05),]),";"),"[[",1)
gns2<-sapply(strsplit(rownames(list_DE2ok[which(list_DE2ok$padj<.05),]),";"),"[[",1)
list_gns<-list(gns1,gns2)
names(list_gns)<-c("combined","ORFquant")
ggVennDiagram(list_gns) + scale_fill_gradient(low="grey90",high = "darkred")
#GTF_annotation$trann$gene_biotype[match(gns1[!gns1%in%gns2],GTF_annotation$trann$gene_id)]
#GTF_annotation$trann$gene_biotype[match(gns2[!gns2%in%gns1],GTF_annotation$trann$gene_id)]

multiDE_all<-get(load("data/14h_vs_00h_multiDE_results_summary_allgenes.RData"))
multiDE_all<-multiDE_all$uniq[[1]]
multiDE_all <- multiDE_all[!multiDE_all$gene_overlapping_noDE,]

#table(multiDE_all[multiDE_all$gene_id%in%gns1[!gns1%in%gns2],]$Ribo_padj<.05 & multiDE_all[multiDE_all$gene_id%in%gns1[!gns1%in%gns2],]$Ribo_log2FC>0 )
#table(multiDE_all[multiDE_all$gene_id%in%gns2[!gns2%in%gns1],]$Ribo_padj<.05 & multiDE_all[multiDE_all$gene_id%in%gns2[!gns2%in%gns1],]$Ribo_log2FC>0 )


#####divide by ORF


peptide_data<-get(load("data/SRonly_UniqueORFq_combined_database_search_peptide_results_annotated.RData"))
df_meta<-peptide_data$metadata
df_meta<-df_meta[elementNROWS(df_meta$Mapped_Genes)==1,]

cat_tx<-unlist(df_meta$ORF_category_Tx_compatible)
cat_tx[cat_tx=="novel"]<-"noCDS_info"
cat_tx[cat_tx=="novel"]<-"noCDS_info"
cat_tx[cat_tx=="overl_uORF"]<-"uORF"
cat_tx[cat_tx=="overl_dORF"]<-"other"
cat_tx[cat_tx=="dORF"]<-"other"
cat_tx[cat_tx%in%c("C_extension","C_truncation","NC_extension","N_extension","nested_ORF","overl_dORF")]<-"other"
cat_tx[cat_tx=="N_truncation"]<-"ORF_annotated"
cat_tx[grep(unlist(df_meta$compatible_biotype),pattern="^novel_tx_")]<-unlist(df_meta$compatible_biotype)[grep(unlist(df_meta$compatible_biotype),pattern="^novel_tx_")]


cat_biot<-unlist(df_meta$compatible_biotype)
cat_biot_gene<-rep(df_meta$gene_biotype,elementNROWS(df_meta$ORF_category_Tx_compatible))
cat_biot[grepl(cat_biot, pattern="^novel_tx_") & cat_biot_gene=="protein_coding"] <- "novel isoform"
cat_biot[!cat_biot %in% c("protein_coding", "novel isoform") & cat_biot_gene=="protein_coding" & !is.na(cat_biot)]<-"non-coding isoform"
cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
cat_biot[!cat_biot%in%c("protein_coding","pseudogene","non-coding isoform", "novel isoform") & !is.na(cat_biot)]<-"non-coding RNA"
cat_biot[grepl(cat_biot_gene, pattern="^novel")] <- "novel gene"


repss<-unlist(apply(cbind(1:dim(df_meta)[1],elementNROWS(df_meta$ORF_category_Tx_compatible)),1,function(x){rep(x[1],x[2])}))
df_meta$ORF_category_Tx_compatible<-CharacterList(split(cat_tx,repss))
df_meta$compatible_biotype<-CharacterList(split(cat_biot,repss))
df_metaok<-df_meta


levvs<-c("ORF_annotated","uORF","other","noCDS_info","novel_tx_R1","novel_tx_R2","novel_tx_R1R2")

df_metaok<-df_meta
df_metaok$ORF_category_Tx_compatible<-unique(df_metaok$ORF_category_Tx_compatible)
df_metaok$compatible_biotype<-unique(df_metaok$compatible_biotype)
if(sum(is.na(unlist(df_metaok$ORF_category_Tx_compatible)))>0){
  df_metaok$ORF_category_Tx_compatible[sum(is.na(df_metaok$ORF_category_Tx_compatible) | df_metaok$ORF_category_Tx_compatible=="ORF_annotated")==elementNROWS(df_metaok$ORF_category_Tx_compatible)]<-"ORF_annotated/in_ref"
  df_metaok$compatible_biotype[sum(is.na(df_metaok$compatible_biotype) | df_metaok$compatible_biotype=="protein_coding")==elementNROWS(df_metaok$compatible_biotype)]<-"protein_coding"
  
  levvs[1]<-"ORF_annotated/in_ref"
}
df_metaunq<-df_metaok[elementNROWS(unique(df_metaok$compatible_biotype))==1 & elementNROWS(unique(df_metaok$ORF_category_Tx_compatible))==1,]

cat_tx<-unlist(df_metaunq$ORF_category_Tx_compatible)
cat_biot<-unlist(df_metaunq$compatible_biotype)

df<-melt(table(cat_tx,cat_biot))

df$value[df$value==0]<-NA
df$cat_tx<-factor(df$cat_tx,levels=levvs)
df$cat_biot<-factor(df$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","non-coding RNA", "novel isoform", "novel gene"))
colli <- c("#4E79A7",  "#F28E2B", "#E15759", "#76B7B2", "#59A14F","#EDC948")  


df<-df[!is.na(df$value),]

df1<-df
df1$search<-"combined"


peptide_data<-get(load("data/SRonly_UniqueORFq_custom_database_search_peptide_results_annotated.RData"))
df_meta<-peptide_data$metadata
df_meta<-df_meta[elementNROWS(df_meta$Mapped_Genes)==1,]

cat_tx<-unlist(df_meta$ORF_category_Tx_compatible)
cat_tx[cat_tx=="novel"]<-"noCDS_info"
cat_tx[cat_tx=="novel"]<-"noCDS_info"
cat_tx[cat_tx=="overl_uORF"]<-"uORF"
cat_tx[cat_tx=="overl_dORF"]<-"other"
cat_tx[cat_tx=="dORF"]<-"other"
cat_tx[cat_tx%in%c("C_extension","C_truncation","NC_extension","N_extension","nested_ORF","overl_dORF")]<-"other"
cat_tx[cat_tx=="N_truncation"]<-"ORF_annotated"
cat_tx[grep(unlist(df_meta$compatible_biotype),pattern="^novel_tx_")]<-unlist(df_meta$compatible_biotype)[grep(unlist(df_meta$compatible_biotype),pattern="^novel_tx_")]


cat_biot<-unlist(df_meta$compatible_biotype)
cat_biot_gene<-rep(df_meta$gene_biotype,elementNROWS(df_meta$ORF_category_Tx_compatible))
cat_biot[grepl(cat_biot, pattern="^novel_tx_") & cat_biot_gene=="protein_coding"] <- "novel isoform"
cat_biot[!cat_biot %in% c("protein_coding", "novel isoform") & cat_biot_gene=="protein_coding" & !is.na(cat_biot)]<-"non-coding isoform"
cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
cat_biot[!cat_biot%in%c("protein_coding","pseudogene","non-coding isoform", "novel isoform") & !is.na(cat_biot)]<-"non-coding RNA"
cat_biot[grepl(cat_biot_gene, pattern="^novel")] <- "novel gene"


repss<-unlist(apply(cbind(1:dim(df_meta)[1],elementNROWS(df_meta$ORF_category_Tx_compatible)),1,function(x){rep(x[1],x[2])}))
df_meta$ORF_category_Tx_compatible<-CharacterList(split(cat_tx,repss))
df_meta$compatible_biotype<-CharacterList(split(cat_biot,repss))
df_metaok<-df_meta


levvs<-c("ORF_annotated","uORF","other","noCDS_info","novel_tx_R1","novel_tx_R2","novel_tx_R1R2")

df_metaok<-df_meta
df_metaok$ORF_category_Tx_compatible<-unique(df_metaok$ORF_category_Tx_compatible)
df_metaok$compatible_biotype<-unique(df_metaok$compatible_biotype)
if(sum(is.na(unlist(df_metaok$ORF_category_Tx_compatible)))>0){
  df_metaok$ORF_category_Tx_compatible[sum(is.na(df_metaok$ORF_category_Tx_compatible) | df_metaok$ORF_category_Tx_compatible=="ORF_annotated")==elementNROWS(df_metaok$ORF_category_Tx_compatible)]<-"ORF_annotated/in_ref"
  df_metaok$compatible_biotype[sum(is.na(df_metaok$compatible_biotype) | df_metaok$compatible_biotype=="protein_coding")==elementNROWS(df_metaok$compatible_biotype)]<-"protein_coding"
  
  levvs[1]<-"ORF_annotated/in_ref"
}
df_metaunq<-df_metaok[elementNROWS(unique(df_metaok$compatible_biotype))==1 & elementNROWS(unique(df_metaok$ORF_category_Tx_compatible))==1,]

cat_tx<-unlist(df_metaunq$ORF_category_Tx_compatible)
cat_biot<-unlist(df_metaunq$compatible_biotype)

df<-melt(table(cat_tx,cat_biot))

df$value[df$value==0]<-NA
df$cat_tx<-factor(df$cat_tx,levels=levvs)
df$cat_biot<-factor(df$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","non-coding RNA", "novel isoform", "novel gene"))
colli <- c("#4E79A7",  "#F28E2B", "#E15759", "#76B7B2", "#59A14F","#EDC948")  


df<-df[!is.na(df$value),]

df2<-df
df2$search<-"ORFquant"
df<-rbind(df1,df2)

a<-ggplot(df,aes(x=cat_biot,y=value+1,fill=cat_biot))
a<-a + geom_bar(stat="identity",position = "dodge")
a<-a + facet_grid(search ~ cat_tx ,drop = T,scales = "free_x",space = "free_x")
a<-a + geom_text(aes(x=cat_biot, y=value+1, hjust=.5,vjust=-.2,label=value),check_overlap = TRUE,colour="black",position = position_dodge(width = .9),size=5)
a<-a + theme_bw()
a<-a + ylab(expression("Peptides uniquely assigned to category"))
a<-a + xlab("")
a<-a + scale_fill_manual(values = colli,"")
a<-a + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
a<-a + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=20))
peps_found_orf<-a + theme(strip.text.x = element_text(size=15, face="bold"),strip.text.y = element_text(size=20), legend.text=element_text(size=20), legend.title=element_text(size=23)) + scale_y_log10(limits = c(1, max(df$value, na.rm = TRUE) * 1.1))

#pdf(file = "Figures/Fig3_ids_ORFs_unq_combined.pdf",width = 13,height = 6)
# peps_found_orf
#dev.off()



peptide_data<-get(load("data/SRonly_UniqueORFq_combined_database_search_peptide_results_annotated.RData"))
df_meta<-peptide_data$metadata
df_meta<-df_meta[elementNROWS(df_meta$Mapped_Genes)==1,]

cat_tx<-unlist(df_meta$ORF_category_Tx_compatible)
cat_tx[cat_tx=="novel"]<-"noCDS_info"
cat_tx[cat_tx=="novel"]<-"noCDS_info"
cat_tx[cat_tx=="overl_uORF"]<-"uORF"
cat_tx[cat_tx=="overl_dORF"]<-"other"
cat_tx[cat_tx=="dORF"]<-"other"
cat_tx[cat_tx%in%c("C_extension","C_truncation","NC_extension","N_extension","nested_ORF","overl_dORF")]<-"other"
cat_tx[cat_tx=="N_truncation"]<-"ORF_annotated"
cat_tx[grep(unlist(df_meta$compatible_biotype),pattern="^novel_tx_")]<-unlist(df_meta$compatible_biotype)[grep(unlist(df_meta$compatible_biotype),pattern="^novel_tx_")]


cat_biot<-unlist(df_meta$compatible_biotype)
cat_biot_gene<-rep(df_meta$gene_biotype,elementNROWS(df_meta$ORF_category_Tx_compatible))
cat_biot[grepl(cat_biot, pattern="^novel_tx_") & cat_biot_gene=="protein_coding"] <- "novel isoform"
cat_biot[!cat_biot %in% c("protein_coding", "novel isoform") & cat_biot_gene=="protein_coding" & !is.na(cat_biot)]<-"non-coding isoform"
cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
cat_biot[!cat_biot%in%c("protein_coding","pseudogene","non-coding isoform", "novel isoform") & !is.na(cat_biot)]<-"non-coding RNA"
cat_biot[grepl(cat_biot_gene, pattern="^novel")] <- "novel gene"


repss<-unlist(apply(cbind(1:dim(df_meta)[1],elementNROWS(df_meta$ORF_category_Tx_compatible)),1,function(x){rep(x[1],x[2])}))
df_meta$ORF_category_Tx_compatible<-CharacterList(split(cat_tx,repss))
df_meta$compatible_biotype<-CharacterList(split(cat_biot,repss))
df_metaok<-df_meta


levvs<-c("ORF_annotated","uORF","other","noCDS_info","novel_tx_R1","novel_tx_R2","novel_tx_R1R2")

df_metaok<-df_meta
df_metaok$ORF_category_Tx_compatible<-unique(df_metaok$ORF_category_Tx_compatible)
df_metaok$compatible_biotype<-unique(df_metaok$compatible_biotype)
if(sum(is.na(unlist(df_metaok$ORF_category_Tx_compatible)))>0){
  df_metaok$ORF_category_Tx_compatible[sum(is.na(df_metaok$ORF_category_Tx_compatible) | df_metaok$ORF_category_Tx_compatible=="ORF_annotated")==elementNROWS(df_metaok$ORF_category_Tx_compatible)]<-"ORF_annotated/in_ref"
  df_metaok$compatible_biotype[sum(is.na(df_metaok$compatible_biotype) | df_metaok$compatible_biotype=="protein_coding")==elementNROWS(df_metaok$compatible_biotype)]<-"protein_coding"
  
  levvs[1]<-"ORF_annotated/in_ref"
}
df_metaunq<-df_metaok[elementNROWS(unique(df_metaok$compatible_biotype))==1 & elementNROWS(unique(df_metaok$ORF_category_Tx_compatible))==1,]

cat_tx<-unlist(df_metaunq$ORF_category_Tx_compatible)
cat_biot<-unlist(df_metaunq$compatible_biotype)

df<-melt(table(cat_tx,cat_biot))

df$value[df$value==0]<-NA
df$cat_tx<-factor(df$cat_tx,levels=levvs)
df$cat_biot<-factor(df$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","non-coding RNA", "novel isoform", "novel gene"))
colli <- c("#4E79A7",  "#F28E2B", "#E15759", "#76B7B2", "#59A14F","#EDC948")  


df<-df[!is.na(df$value),]

df3<-df
df3$search<-"R1_combined"

peptide_data<-get(load("data/LRalso_UniqueORFq_combined_database_search_peptide_results_annotated.RData"))
df_meta<-peptide_data$metadata
df_meta<-df_meta[elementNROWS(df_meta$Mapped_Genes)==1,]

cat_tx<-unlist(df_meta$ORF_category_Tx_compatible)
cat_tx[cat_tx=="novel"]<-"noCDS_info"
cat_tx[cat_tx=="novel"]<-"noCDS_info"
cat_tx[cat_tx=="overl_uORF"]<-"uORF"
cat_tx[cat_tx=="overl_dORF"]<-"other"
cat_tx[cat_tx=="dORF"]<-"other"
cat_tx[cat_tx%in%c("C_extension","C_truncation","NC_extension","N_extension","nested_ORF","overl_dORF")]<-"other"
cat_tx[cat_tx=="N_truncation"]<-"ORF_annotated"
cat_tx[grep(unlist(df_meta$compatible_biotype),pattern="^novel_tx_")]<-unlist(df_meta$compatible_biotype)[grep(unlist(df_meta$compatible_biotype),pattern="^novel_tx_")]


cat_biot<-unlist(df_meta$compatible_biotype)
cat_biot_gene<-rep(df_meta$gene_biotype,elementNROWS(df_meta$ORF_category_Tx_compatible))
cat_biot[grepl(cat_biot, pattern="^novel_tx_") & cat_biot_gene=="protein_coding"] <- "novel isoform"
cat_biot[!cat_biot %in% c("protein_coding", "novel isoform") & cat_biot_gene=="protein_coding" & !is.na(cat_biot)]<-"non-coding isoform"
cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
cat_biot[!cat_biot%in%c("protein_coding","pseudogene","non-coding isoform", "novel isoform") & !is.na(cat_biot)]<-"non-coding RNA"
cat_biot[grepl(cat_biot_gene, pattern="^novel")] <- "novel gene"


repss<-unlist(apply(cbind(1:dim(df_meta)[1],elementNROWS(df_meta$ORF_category_Tx_compatible)),1,function(x){rep(x[1],x[2])}))
df_meta$ORF_category_Tx_compatible<-CharacterList(split(cat_tx,repss))
df_meta$compatible_biotype<-CharacterList(split(cat_biot,repss))
df_metaok<-df_meta


levvs<-c("ORF_annotated","uORF","other","noCDS_info","novel_tx_R1","novel_tx_R2","novel_tx_R1R2")

df_metaok<-df_meta
df_metaok$ORF_category_Tx_compatible<-unique(df_metaok$ORF_category_Tx_compatible)
df_metaok$compatible_biotype<-unique(df_metaok$compatible_biotype)
if(sum(is.na(unlist(df_metaok$ORF_category_Tx_compatible)))>0){
  df_metaok$ORF_category_Tx_compatible[sum(is.na(df_metaok$ORF_category_Tx_compatible) | df_metaok$ORF_category_Tx_compatible=="ORF_annotated")==elementNROWS(df_metaok$ORF_category_Tx_compatible)]<-"ORF_annotated/in_ref"
  df_metaok$compatible_biotype[sum(is.na(df_metaok$compatible_biotype) | df_metaok$compatible_biotype=="protein_coding")==elementNROWS(df_metaok$compatible_biotype)]<-"protein_coding"
  
  levvs[1]<-"ORF_annotated/in_ref"
}
df_metaunq<-df_metaok[elementNROWS(unique(df_metaok$compatible_biotype))==1 & elementNROWS(unique(df_metaok$ORF_category_Tx_compatible))==1,]

cat_tx<-unlist(df_metaunq$ORF_category_Tx_compatible)
cat_biot<-unlist(df_metaunq$compatible_biotype)

df<-melt(table(cat_tx,cat_biot))

df$value[df$value==0]<-NA
df$cat_tx<-factor(df$cat_tx,levels=levvs)
df$cat_biot<-factor(df$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","non-coding RNA", "novel isoform", "novel gene"))
colli <- c("#4E79A7",  "#F28E2B", "#E15759", "#76B7B2", "#59A14F","#EDC948")  


df<-df[!is.na(df$value),]

df4<-df
df4$search<-"R1R2_combined"

df<-rbind(df3,df4)

a<-ggplot(df,aes(x=cat_biot,y=value+1,fill=cat_biot))
a<-a + geom_bar(stat="identity",position = "dodge")
a<-a + facet_grid(search ~ cat_tx ,drop = T,scales = "free_x",space = "free_x")
a<-a + geom_text(aes(x=cat_biot, y=value+1, hjust=.5,vjust=-.2,label=value),check_overlap = TRUE,colour="black",position = position_dodge(width = .9),size=5)
a<-a + theme_bw()
a<-a + ylab(expression("Peptides uniquely assigned to category"))
a<-a + xlab("")
a<-a + scale_fill_manual(values = colli,"")
a<-a + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
a<-a + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=20))
peps_found_orf<-a + theme(strip.text.x = element_text(size=15, face="bold"),strip.text.y = element_text(size=20), legend.text=element_text(size=20), legend.title=element_text(size=23)) + scale_y_log10(limits = c(1, max(df$value, na.rm = TRUE) * 1.1))

# pdf(file = "Figures/Fig5_ids_ORFs_R1R2unq_ref_combined.pdf",width = 13,height = 6)
# peps_found_orf
# dev.off()


peptide_data<-get(load("data/SRonly_UniqueORFq_combined_database_search_peptide_results_annotated.RData"))
de_data<-get(load("data/SRonly_UniqueORFq_combined_database_search_proteomics_multiDE.RData"))

#de_data<-get(load(gsub("/group/calviello/Lorenzo/rob_runs/SRonly_UniqueORFq/R2T2P_output_dir/annot_and_orfquant_prot_results_peptide_results_annotated.RData",pattern="peptide_results_annotated.RData",replacement="proteomics_multiDE.RData.RData")))
df_meta<-peptide_data$metadata
df_meta<-df_meta[elementNROWS(df_meta$Mapped_Genes)==1,]
df_meta$up<-F
df_meta$up[which(de_data$raw$DE$Peptides$`DUX4_DMSO_14h-Control_DMSO_14h`$padj<.05 & de_data$raw$DE$Peptides$`DUX4_DMSO_14h-Control_DMSO_14h`$logFC>0)]<-T



cat_tx<-unlist(df_meta$ORF_category_Tx_compatible)
cat_tx[cat_tx=="novel"]<-"noCDS_info"
cat_tx[cat_tx=="novel"]<-"noCDS_info"
cat_tx[cat_tx=="overl_uORF"]<-"uORF"
cat_tx[cat_tx=="overl_dORF"]<-"other"
cat_tx[cat_tx=="dORF"]<-"other"
cat_tx[cat_tx%in%c("C_extension","C_truncation","NC_extension","N_extension","nested_ORF","overl_dORF")]<-"other"
cat_tx[cat_tx=="N_truncation"]<-"ORF_annotated"
cat_tx[grep(unlist(df_meta$compatible_biotype),pattern="^novel_tx_")]<-unlist(df_meta$compatible_biotype)[grep(unlist(df_meta$compatible_biotype),pattern="^novel_tx_")]


cat_biot<-unlist(df_meta$compatible_biotype)
cat_biot_gene<-rep(df_meta$gene_biotype,elementNROWS(df_meta$ORF_category_Tx_compatible))
cat_biot[grepl(cat_biot, pattern="^novel_tx_") & cat_biot_gene=="protein_coding"] <- "novel isoform"
cat_biot[!cat_biot %in% c("protein_coding", "novel isoform") & cat_biot_gene=="protein_coding" & !is.na(cat_biot)]<-"non-coding isoform"
cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
cat_biot[!cat_biot%in%c("protein_coding","pseudogene","non-coding isoform", "novel isoform") & !is.na(cat_biot)]<-"non-coding RNA"
cat_biot[grepl(cat_biot_gene, pattern="^novel")] <- "novel gene"


repss<-unlist(apply(cbind(1:dim(df_meta)[1],elementNROWS(df_meta$ORF_category_Tx_compatible)),1,function(x){rep(x[1],x[2])}))
df_meta$ORF_category_Tx_compatible<-CharacterList(split(cat_tx,repss))
df_meta$compatible_biotype<-CharacterList(split(cat_biot,repss))
df_metaok<-df_meta


levvs<-c("ORF_annotated","uORF","other","noCDS_info","novel_tx_R1","novel_tx_R2","novel_tx_R1R2")

df_metaok<-df_meta
df_metaok$ORF_category_Tx_compatible<-unique(df_metaok$ORF_category_Tx_compatible)
df_metaok$compatible_biotype<-unique(df_metaok$compatible_biotype)
if(sum(is.na(unlist(df_metaok$ORF_category_Tx_compatible)))>0){
  df_metaok$ORF_category_Tx_compatible[sum(is.na(df_metaok$ORF_category_Tx_compatible) | df_metaok$ORF_category_Tx_compatible=="ORF_annotated")==elementNROWS(df_metaok$ORF_category_Tx_compatible)]<-"ORF_annotated/in_ref"
  df_metaok$compatible_biotype[sum(is.na(df_metaok$compatible_biotype) | df_metaok$compatible_biotype=="protein_coding")==elementNROWS(df_metaok$compatible_biotype)]<-"protein_coding"
  
  levvs[1]<-"ORF_annotated/in_ref"
}
df_metaunq<-df_metaok[elementNROWS(unique(df_metaok$compatible_biotype))==1 & elementNROWS(unique(df_metaok$ORF_category_Tx_compatible))==1,]

cat_tx<-unlist(df_metaunq$ORF_category_Tx_compatible)
cat_biot<-unlist(df_metaunq$compatible_biotype)

df<-melt(table(cat_tx,cat_biot,df_metaunq$up))

df$value[df$value==0]<-NA
df$cat_tx<-factor(df$cat_tx,levels=levvs)
df$cat_biot<-factor(df$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","non-coding RNA", "novel isoform", "novel gene"))
colli <- c("#4E79A7",  "#F28E2B", "#E15759", "#76B7B2", "#59A14F","#EDC948")  


df<-df[!is.na(df$value),]

df3<-df
df3$search<-"R1_combined"


#peptide_data<-get(load("/group/calviello/Lorenzo/rob_runs/LRalso_UniqueORFq/R2T2P_output_dir/annot_and_orfquant_prot_results_peptide_results_annotated.RData"))
#de_data<-get(load(gsub("/group/calviello/Lorenzo/rob_runs/LRalso_UniqueORFq/R2T2P_output_dir/annot_and_orfquant_prot_results_peptide_results_annotated.RData",pattern="peptide_results_annotated.RData",replacement="proteomics_multiDE.RData.RData")))

peptide_data<-get(load("data/LRalso_UniqueORFq_combined_database_search_peptide_results_annotated.RData"))
de_data<-get(load("data/LRalso_UniqueORFq_combined_database_search_proteomics_multiDE.RData"))

df_meta<-peptide_data$metadata
df_meta<-df_meta[elementNROWS(df_meta$Mapped_Genes)==1,]
df_meta$up<-F
df_meta$up[which(de_data$raw$DE$Peptides$`DUX4_DMSO_14h-Control_DMSO_14h`$padj<.05 & de_data$raw$DE$Peptides$`DUX4_DMSO_14h-Control_DMSO_14h`$logFC>0)]<-T


cat_tx<-unlist(df_meta$ORF_category_Tx_compatible)
cat_tx[cat_tx=="novel"]<-"noCDS_info"
cat_tx[cat_tx=="novel"]<-"noCDS_info"
cat_tx[cat_tx=="overl_uORF"]<-"uORF"
cat_tx[cat_tx=="overl_dORF"]<-"other"
cat_tx[cat_tx=="dORF"]<-"other"
cat_tx[cat_tx%in%c("C_extension","C_truncation","NC_extension","N_extension","nested_ORF","overl_dORF")]<-"other"
cat_tx[cat_tx=="N_truncation"]<-"ORF_annotated"
cat_tx[grep(unlist(df_meta$compatible_biotype),pattern="^novel_tx_")]<-unlist(df_meta$compatible_biotype)[grep(unlist(df_meta$compatible_biotype),pattern="^novel_tx_")]


cat_biot<-unlist(df_meta$compatible_biotype)
cat_biot_gene<-rep(df_meta$gene_biotype,elementNROWS(df_meta$ORF_category_Tx_compatible))
cat_biot[grepl(cat_biot, pattern="^novel_tx_") & cat_biot_gene=="protein_coding"] <- "novel isoform"
cat_biot[!cat_biot %in% c("protein_coding", "novel isoform") & cat_biot_gene=="protein_coding" & !is.na(cat_biot)]<-"non-coding isoform"
cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
cat_biot[!cat_biot%in%c("protein_coding","pseudogene","non-coding isoform", "novel isoform") & !is.na(cat_biot)]<-"non-coding RNA"
cat_biot[grepl(cat_biot_gene, pattern="^novel")] <- "novel gene"


repss<-unlist(apply(cbind(1:dim(df_meta)[1],elementNROWS(df_meta$ORF_category_Tx_compatible)),1,function(x){rep(x[1],x[2])}))
df_meta$ORF_category_Tx_compatible<-CharacterList(split(cat_tx,repss))
df_meta$compatible_biotype<-CharacterList(split(cat_biot,repss))
df_metaok<-df_meta


levvs<-c("ORF_annotated","uORF","other","noCDS_info","novel_tx_R1","novel_tx_R2","novel_tx_R1R2")

df_metaok<-df_meta
df_metaok$ORF_category_Tx_compatible<-unique(df_metaok$ORF_category_Tx_compatible)
df_metaok$compatible_biotype<-unique(df_metaok$compatible_biotype)
if(sum(is.na(unlist(df_metaok$ORF_category_Tx_compatible)))>0){
  df_metaok$ORF_category_Tx_compatible[sum(is.na(df_metaok$ORF_category_Tx_compatible) | df_metaok$ORF_category_Tx_compatible=="ORF_annotated")==elementNROWS(df_metaok$ORF_category_Tx_compatible)]<-"ORF_annotated/in_ref"
  df_metaok$compatible_biotype[sum(is.na(df_metaok$compatible_biotype) | df_metaok$compatible_biotype=="protein_coding")==elementNROWS(df_metaok$compatible_biotype)]<-"protein_coding"
  
  levvs[1]<-"ORF_annotated/in_ref"
}
df_metaunq<-df_metaok[elementNROWS(unique(df_metaok$compatible_biotype))==1 & elementNROWS(unique(df_metaok$ORF_category_Tx_compatible))==1,]

cat_tx<-unlist(df_metaunq$ORF_category_Tx_compatible)
cat_biot<-unlist(df_metaunq$compatible_biotype)

df<-melt(table(cat_tx,cat_biot,df_metaunq$up))

df$value[df$value==0]<-NA
df$cat_tx<-factor(df$cat_tx,levels=levvs)
df$cat_biot<-factor(df$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","non-coding RNA", "novel isoform", "novel gene"))
colli <- c("#4E79A7",  "#F28E2B", "#E15759", "#76B7B2", "#59A14F","#EDC948")  


df<-df[!is.na(df$value),]

df4<-df
df4$search<-"R1R2_combined"

df<-rbind(df3,df4)


a<-ggplot(df,aes(x=cat_biot,y=value+1,fill=cat_biot,alpha=Var3))
a<-a + geom_bar(stat="identity",position = "stack")
a<-a + facet_grid(search ~ cat_tx ,drop = T,scales = "free_x",space = "free_x")
a<-a + geom_text(aes(x=cat_biot, y=value+1, hjust=.5,vjust=-.2,label=value),check_overlap = F,colour="black",position = "stack",size=5)
a<-a + theme_bw()
a<-a + ylab(expression("Peptides uniquely assigned to category"))
a<-a + xlab("")
a<-a + scale_fill_manual(values = colli,"")
a<-a + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
a<-a + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=20))
peps_found_orf<-a + theme(strip.text.x = element_text(size=15, face="bold"),strip.text.y = element_text(size=20), legend.text=element_text(size=20), legend.title=element_text(size=23)) + scale_y_log10()

# pdf(file = "Figures/Fig5_ids_ORFs_R1R2unq_ref_combined_sig.pdf",width = 13,height = 6)
# peps_found_orf + scale_alpha_manual(values = c(.5,1),"Upregulated")
# dev.off()

list_DE<-get(load("data/LRalso_UniqueORFq_combined_database_search_proteomics_multiDE.RData"))
peptide_data<-get(load("data/LRalso_UniqueORFq_combined_database_search_peptide_results_annotated.RData"))

#list_DE<-get(load("/group/calviello/Lorenzo/rob_runs/LRalso_UniqueORFq/R2T2P_output_dir/annot_and_orfquant_prot_results_proteomics_multiDE.RData.RData"))
#peptide_data<-get(load("/group/calviello/Lorenzo/rob_runs/LRalso_UniqueORFq/R2T2P_output_dir/annot_and_orfquant_prot_results_peptide_results_annotated.RData"))

df_meta<-peptide_data$metadata
DE_res<-list_DE$raw$DE
multi<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))>1,]
df_meta<-df_meta[elementNROWS(unique(df_meta$Mapped_Genes))==1,]

df_meta$gene_biotype[grepl(df_meta$gene_name,pattern="^IG") & df_meta$gene_biotype=="other"]<-"protein_coding"
df_meta$gene_biotype[df_meta$gene_biotype!="protein_coding"]<-"non-coding"
df_meta$DUX4_target<-df_meta$gene_id%in%target_genes
df_meta$gene_biotype[df_meta$DUX4_target]<-"DUX4_target"


df_meta_unq<-unique(df_meta[,c("gene_id","gene_biotype","novelty_category")])

pepsig<-DE_res$Peptides[["DUX4_DMSO_14h-Control_DMSO_14h"]]
pepsig<-rownames(pepsig)[which(pepsig$padj<.05 & pepsig$logFC>0)]
ginsig<-DE_res$gene_id_novelty[["DUX4_DMSO_14h-Control_DMSO_14h"]]
ginsig<-rownames(ginsig)[which(ginsig$padj<.05 & ginsig$logFC>0)]

df_meta$sig_up<-rownames(df_meta)%in%pepsig
df_meta_unq$sig_up<-paste(df_meta_unq$gene_id,df_meta_unq$novelty_category,sep=";")%in%ginsig

df<-melt(table(df_meta$gene_biotype,df_meta$novelty_category,df_meta$sig_up))

df$id_type<-"Peptides"
df$Var1<-factor(df$Var1,levels=names(sort(table(df_meta$gene_biotype),decreasing = T)))

df2<-melt(table(df_meta_unq$gene_biotype,df_meta_unq$novelty_category,df_meta_unq$sig_up))

df2$Var1<-factor(df2$Var1,levels=names(sort(table(df_meta$gene_biotype),decreasing = T)))
df2$id_type<-"Genes"

df<-rbind(df,df2)

df$value2<-df$value
df$value2[!df$Var3]<-df$value[!df$Var3]+df$value[df$Var3]

df<-df[df$value>0,]
df$Var2<-factor(as.character(df$Var2),levels=c("in_ref","novel","novel_tx_R1","novel_tx_R2","novel_tx_R1R2"))

df$sigup<-"all"
df$sigup[df$Var3]<-"Upregulated"
df$sigup <- factor(df$sigup, levels = c("all", "Upregulated"))

gg_idsig <- ggplot(data = df, aes(x=Var2, y=value+1, fill=Var2,label=value2,alpha=sigup)) + 
  #geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
  geom_bar(stat="identity", position="stack") + 
  #geom_col(stat="identity", position=position_dodge(width = .9),preserve = c("total", "single")) + 
  
  labs(x="", y="# identifications") +
  scale_fill_manual(values = c("grey","orange","red1","red3","darkred"),"peptide_type") +
  theme_bw() +
  scale_y_log10() +
  scale_alpha_manual(values = c(.5,1),"")+
  theme(plot.margin = ggplot2::margin(.5, .5, .5, .5, "cm"),
        panel.background=element_blank(), 
        axis.ticks=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=14), 
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20),
        #legend.position="bottom",
        legend.text=element_text(size=18),
        legend.title=element_text(size=20),
        panel.grid.major.y = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.x = element_line( size=.1, color="dark grey" ,linetype=2)) +
  geom_text(aes(x=Var2, y=(value+1)*1.1,label=value2),position="stack",check_overlap = F,show.legend = FALSE) + facet_wrap(id_type~Var1,nrow = 2,scales = "fixed")
# pdf(file = "Figures/Fig5_idsSIG.pdf",width = 9,height = 6)
# gg_idsig
# dev.off()


df2<-df[df$Var2!="in_ref",]
gg_idsig2 <- ggplot(data = df2, aes(x=Var2, y=value+1, fill=Var2,label=value2,alpha=sigup)) + 
  #geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
  geom_bar(stat="identity", position="stack") + 
  #geom_col(stat="identity", position=position_dodge(width = .9),preserve = c("total", "single")) + 
  
  labs(x="", y="# identifications") +
  scale_fill_manual(values = c("orange","red1","red3","darkred"),"peptide_type") +
  theme_bw() +
  scale_y_log10() +
  scale_alpha_manual(values = c(.5,1),"")+
  theme(plot.margin = ggplot2::margin(.5, .5, .5, .5, "cm"),
        panel.background=element_blank(), 
        axis.ticks=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, size=14), 
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=20),
        #legend.position="bottom",
        legend.text=element_text(size=18),
        legend.title=element_text(size=20),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_text(aes(x=Var2, y=(value+1)*1.1,label=value2),position="stack",check_overlap = F,show.legend = FALSE) + facet_grid(id_type~Var1,scales = "free_x",space = "free_x")
pdf(file = "Figures/Figure5_E.pdf",width = 9,height = 6)
gg_idsig2
dev.off()
