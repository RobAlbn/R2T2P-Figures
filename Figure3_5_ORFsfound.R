library(cowplot)
library(zoo)
library(reshape2)
library(limma)
library(UpSetR)
library(readxl)
library(ggfortify)

suppressPackageStartupMessages(library("RiboseQC"))
suppressPackageStartupMessages(library("ggfortify"))
trans <- function(x) {reso=-log(x, 10);reso[is.na(reso)]=0;return(reso)}
inv <- function(x) 10^(-x)

setwd("/group/calviello/Roberto/R2T2P_scripts_for_Figures_with_container/R2T2P_Figures_GitHub_repo")
source("data/new_riboseqc_all_HT_oct23_26.R")

GTF_all<-get(load("data/LRalso_UniqueORFq_Rannot"))
load("data/SRonly_UniqueORFq_Rannot")

# Loading ORFquant results  

orfs_sr_all<-get(load("data/SRonly_AllORFq_ORFquant_results"))
orfs_sr_unq<-get(load("data/SRonly_UniqueORFq_ORFquant_results"))
orfs_srlr_all<-get(load("data/LRalso_AllORFq_ORFquant_results"))
orfs_srlr_unq<-get(load("data/LRalso_UniqueORFq_ORFquant_results"))


ORFquant_results<-orfs_sr_unq
ORFs_tx<-ORFquant_results$ORFs_tx
ORFs_gen<-ORFquant_results$ORFs_gen
ORFs_gen<-split(ORFs_gen,names(ORFs_gen))


#make R2 instead of PB when problem is solved
ORFs_tx$tx_novel<-sum(!grepl(ORFs_tx$compatible_with,pattern="^R1|^R2"))==0
ORFs_tx$tx_novel_R1R2<-"none"
ORFs_tx$tx_novel_R1R2[ORFs_tx$tx_novel]<-"R1R2"
ORFs_tx$tx_novel_R1R2[sum(!grepl(ORFs_tx$compatible_with,pattern="^R1"))==0]<-"R1"
ORFs_tx$tx_novel_R1R2[sum(!grepl(ORFs_tx$compatible_with,pattern="^R2"))==0]<-"R2"

ORFs_tx$compatible_biotype[ORFs_tx$tx_novel]<-paste("novel_tx",ORFs_tx$tx_novel_R1R2[ORFs_tx$tx_novel],sep="_")
supp3<-mcols(ORFs_tx)
supp3<-supp3[,grepl(colnames(supp3),pattern="ORF_id_tr|pval|Protein|gene_id|gene_name|ORF_pct_P_sites|ORF_pM|tx_novel_R1R2|compatible_biotype|ORF_category_Tx_compatible|ORF_category_Gen")]
supp3$ORF_pct_P_sites_pN<-NULL
supp3$pval_uniq<-NULL
supp3<-supp3[,c(2,1,3:11)]
rownames(supp3)<-NULL
supp3<-data.frame(supp3)
write.table(supp3, file="/group/calviello/Roberto/R2T2P_scripts_for_Figures_with_container/R2T2P_Figures_GitHub_repo/Supplementary_Tables/Supplementary_Table3.tsv",sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE)




gids<-ORFs_tx$gene_id
cat_biot<-ORFs_tx$gene_biotype
cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
cat_biot[!cat_biot%in%c("protein_coding","pseudogene")]<-"non-coding RNA"
tbid<-table(gids)
tbid[tbid>3]<-">3"
#boxplot(log(tpms)~tbid)
tb_bio<-cat_biot[match(names(tbid),gids)]


cat_tx<-ORFs_tx$ORF_category_Tx_compatible
cat_tx[cat_tx=="novel"]<-"noCDS_info"
cat_tx[cat_tx=="overl_uORF"]<-"uORF"
cat_tx[cat_tx=="overl_dORF"]<-"other"
cat_tx[cat_tx=="dORF"]<-"other"
cat_tx[grep(ORFs_tx$compatible_biotype,pattern="^novel_tx_")]<-ORFs_tx$compatible_biotype[grep(ORFs_tx$compatible_biotype,pattern="^novel_tx_")]

cat_gen<-ORFs_tx$ORF_category_Gen
cat_biot<-ORFs_tx$compatible_biotype
cat_biot_gene<-ORFs_tx$gene_biotype
cat_biot[grepl(cat_biot, pattern="^novel_tx_") & cat_biot_gene=="protein_coding"] <- "novel isoform"
cat_biot[!cat_biot %in% c("protein_coding", "novel isoform") & cat_biot_gene=="protein_coding"]<-"non-coding isoform"
cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
cat_biot[!cat_biot%in%c("protein_coding","pseudogene","non-coding isoform", "novel isoform")]<-"non-coding RNA"
cat_biot[grepl(cat_biot_gene, pattern="^novel")] <- "novel gene"

#cat_tx[intersect(grep(ORFs_tx$gene_biotype,pattern = "pseudogene"),which(cat_tx=="novel"))]<-"novel_pseudogene"
cat_tx[cat_tx%in%c("C_extension","C_truncation","NC_extension","N_extension","nested_ORF","overl_dORF")]<-"other"
cat_tx[cat_tx=="N_truncation"]<-"ORF_annotated"

cat_gen[grep(cat_gen,pattern = "novel")]<-"non-overlapping\nCDS regions"
cat_gen[grep(cat_gen,pattern = "non-overlapping\nCDS regions",invert = T)]<-"overlapping\nCDS regions"
levvs<-c("ORF_annotated","uORF","other","noCDS_info","novel_tx_R1","novel_tx_R2","novel_tx_R1R2")

df<-melt(table(cat_tx,cat_biot))

df$value[df$value==0]<-NA
df$cat_tx<-factor(df$cat_tx,levels=levvs)
df$cat_biot<-factor(df$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","non-coding RNA", "novel isoform", "novel gene"))


colli <- c("#4E79A7",  "#F28E2B", "#E15759", "#76B7B2", "#59A14F","#EDC948")  

dfiso<-data.frame(cat_tx,cat_biot,ORFs_tx$ORF_pct_P_sites)
dfiso$cat_tx<-factor(dfiso$cat_tx,levels=levvs)
dfiso$cat_biot<-factor(dfiso$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","non-coding RNA", "novel isoform", "novel gene"))

dfpnpm<-data.frame(cat_tx,cat_biot,ORFs_tx$ORFs_pM)
#dfpnpm$value[dfpnpm$value==0]<-NA
dfpnpm$cat_tx<-factor(dfpnpm$cat_tx,levels=levvs)
dfpnpm$cat_biot<-factor(dfpnpm$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","non-coding RNA", "novel isoform", "novel gene"))

dfwid<-data.frame(cat_tx,cat_biot,width(ORFs_tx))
dfwid$cat_tx<-factor(dfwid$cat_tx,levels=levvs)
dfwid$cat_biot<-factor(dfwid$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","non-coding RNA", "novel isoform", "novel gene"))

df<-df[!is.na(df$value),]


a<-ggplot(df,aes(x=cat_biot,y=value,fill=cat_biot))
a<-a + geom_bar(stat="identity",position = "dodge")
a<-a + facet_grid(. ~ cat_tx ,drop = T,scales = "free_x", space = "free_x")
a<-a + geom_text(aes(x=cat_biot, y=value, hjust=.5,vjust=-.2,label=value),check_overlap = TRUE,colour="black",position = position_dodge(width = .9),size=5)
a<-a + theme_bw()
a<-a + ylab(expression("Number of ORFs"))
a<-a + xlab("")
a<-a + scale_fill_manual(values = colli,"Biotype")
a<-a + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
a<-a + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=20))
orfs_found_n<-a + theme(strip.text.x = element_text(size=15, face="bold"),strip.text.y = element_text(size=20), legend.text=element_text(size=20), legend.title=element_text(size=23), panel.grid = element_blank()) + scale_y_log10(limits = c(1, max(df$value, na.rm = TRUE) * 1.1))



d<-ggplot(dfwid,aes(x=cat_biot,y=width.ORFs_tx.,fill=cat_biot))
d<-d + geom_violin(scale="width",draw_quantiles=.5,adjust = 1)
#c<-c + geom_jitter(aes(x=cat_tx,y=ORFs_tx.ORFs_pM),position = position_jitterdodge(jitter.width = .5), alpha = 0.2)
d<-d + scale_y_log10(breaks=c(10,100,1000,10000),labels=c(10,100,1000,10000),limits=c(min(dfwid$width.ORFs_tx.),max(dfwid$width.ORFs_tx.)*1.1))
#c<-c + geom_text(aes(x=cat_tx, y=ORFs_tx.ORFs_pM, hjust="top",label=ORFs_tx.ORFs_pM),colour="black",position = position_dodge(width = 1),size=5)
d<-d + theme_bw()
d<-d + facet_grid(. ~ cat_tx ,drop = T,scales = "free_x", space = "free_x")
d<-d + xlab("")
d<-d + ylab("ORF length (nt)")
d<-d + theme(axis.title.x = element_blank(),axis.text.x = element_blank())
#c<-c + theme(legend.position="none")
d<-d + scale_fill_manual(values = colli,"Biotype")
d<-d + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
d<-d + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=20))
orfs_found_len<-d + theme(strip.text.x = element_text(size=15, face="bold"),strip.text.y = element_text(size=20), legend.text=element_text(size=20), legend.title=element_text(size=23))

#orfs_found_n<-orfs_found_n+theme(legend.position = "none")


#colli <- c("#4E79A7",  "#F28E2B", "#E15759", "#76B7B2", "#59A14F","#EDC948","#76B7B2", "#59A14F","#EDC948")  

d<-ggplot(dfiso,aes(x=cat_biot,y=ORFs_tx.ORF_pct_P_sites,fill=cat_biot))
d<-d + geom_violin(scale="width",draw_quantiles=.5,adjust = 1)
#c<-c + geom_jitter(aes(x=cat_tx,y=ORFs_tx.ORFs_pM),position = position_jitterdodge(jitter.width = .5), alpha = 0.2)
#d<-d + scale_y_log10(breaks=c(10,100,1000,10000),labels=c(10,100,1000,10000),limits=c(min(dfwid$width.ORFs_tx.),max(dfwid$width.ORFs_tx.)*1.1))
#c<-c + geom_text(aes(x=cat_tx, y=ORFs_tx.ORFs_pM, hjust="top",label=ORFs_tx.ORFs_pM),colour="black",position = position_dodge(width = 1),size=5)
d<-d + theme_bw()
d<-d + facet_grid(. ~ cat_tx ,drop = T,scales = "free_x", space = "free_x")
d<-d + xlab("")
d<-d + ylab("ORF_pct_Psites")
d<-d + theme(axis.title.x = element_blank(),axis.text.x = element_blank())
#c<-c + theme(legend.position="none")
d<-d + scale_fill_manual(values = colli,"Biotype")
d<-d + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
d<-d + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=20))
orfs_pct<-d + theme(strip.text.x = element_text(size=13, face="bold"),strip.text.y = element_text(size=20), legend.text=element_text(size=20), legend.title=element_text(size=23))



# pdf(file = "Figures/Fig3_ORFs.pdf",width = 14,height = 4)
# orfs_found_n
# dev.off()

pdf(file = "Figures/Figure3_A.pdf",width = 14,height = 4)
orfs_found_n + theme(legend.position='top')
dev.off()

pdf(file = "Figures/SupplementaryFigure6.pdf",width = 14,height = 8)
plot_grid(orfs_found_n,orfs_found_len,orfs_pct,nrow = 3,rel_widths = c(1,1,1),align = "hv")
dev.off()

colsi2<-c("red", "blue", "purple")
colsi2<-c("red1","red3","darkred")


ORFquant_results<-orfs_srlr_unq
ORFs_tx<-ORFquant_results$ORFs_tx
ORFs_gen<-ORFquant_results$ORFs_gen
ORFs_gen<-split(ORFs_gen,names(ORFs_gen))

#make R2 instead of PB when problem is solved
ORFs_tx$tx_novel<-sum(!grepl(ORFs_tx$compatible_with,pattern="^R1|^R2"))==0
ORFs_tx$tx_novel_R1R2<-"none"
ORFs_tx$tx_novel_R1R2[ORFs_tx$tx_novel]<-"R1R2"
ORFs_tx$tx_novel_R1R2[sum(!grepl(ORFs_tx$compatible_with,pattern="^R1"))==0]<-"R1"
ORFs_tx$tx_novel_R1R2[sum(!grepl(ORFs_tx$compatible_with,pattern="^R2"))==0]<-"R2"

ORFs_tx$compatible_biotype[ORFs_tx$tx_novel]<-paste("novel_tx",ORFs_tx$tx_novel_R1R2[ORFs_tx$tx_novel],sep="_")

gids<-ORFs_tx$gene_id
cat_biot<-ORFs_tx$gene_biotype
cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
cat_biot[!cat_biot%in%c("protein_coding","pseudogene")]<-"non-coding RNA"
tbid<-table(gids)
tbid[tbid>3]<-">3"
#boxplot(log(tpms)~tbid)
tb_bio<-cat_biot[match(names(tbid),gids)]


cat_tx<-ORFs_tx$ORF_category_Tx_compatible
cat_tx[cat_tx=="novel"]<-"no_CDS_info"
cat_tx[cat_tx=="overl_uORF"]<-"uORF"
cat_tx[cat_tx=="overl_dORF"]<-"other"
cat_tx[cat_tx=="dORF"]<-"other"
cat_tx[grep(ORFs_tx$compatible_biotype,pattern="^novel_tx_")]<-ORFs_tx$compatible_biotype[grep(ORFs_tx$compatible_biotype,pattern="^novel_tx_")]

cat_gen<-ORFs_tx$ORF_category_Gen
cat_biot<-ORFs_tx$compatible_biotype
cat_biot_gene<-ORFs_tx$gene_biotype
cat_biot[grepl(cat_biot, pattern="^novel_tx_") & cat_biot_gene=="protein_coding"] <- "novel isoform"
cat_biot[!cat_biot %in% c("protein_coding", "novel isoform") & cat_biot_gene=="protein_coding"]<-"non-coding isoform"
cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
cat_biot[!cat_biot%in%c("protein_coding","pseudogene","non-coding isoform", "novel isoform")]<-"non-coding RNA"
cat_biot[grepl(cat_biot_gene, pattern="^novel")] <- "novel gene"

#cat_tx[intersect(grep(ORFs_tx$gene_biotype,pattern = "pseudogene"),which(cat_tx=="novel"))]<-"novel_pseudogene"
cat_tx[cat_tx%in%c("C_extension","C_truncation","NC_extension","N_extension","nested_ORF","overl_dORF")]<-"other"

cat_gen[grep(cat_gen,pattern = "novel")]<-"non-overlapping\nCDS regions"
cat_gen[grep(cat_gen,pattern = "non-overlapping\nCDS regions",invert = T)]<-"overlapping\nCDS regions"
levvs<-c("ORF_annotated","N_truncation","uORF","other","no_CDS_info","novel_tx_R1","novel_tx_R2","novel_tx_R1R2")

df<-melt(table(cat_tx,cat_biot))

df$value[df$value==0]<-NA
df$cat_tx<-factor(df$cat_tx,levels=levvs)
df$cat_biot<-factor(df$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","novel isoform", "non-coding RNA", "novel gene"))


colli <- c("#4E79A7",  "#F28E2B", "#E15759", "#76B7B2", "#59A14F","#EDC948")  

dfiso<-data.frame(cat_tx,cat_biot,ORFs_tx$ORF_pct_P_sites)
dfiso$cat_tx<-factor(dfiso$cat_tx,levels=levvs)
dfiso$cat_biot<-factor(dfiso$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","novel isoform","non-coding RNA", "novel gene"))
dfiso<-dfiso[grep(dfiso$cat_tx,pattern="novel_tx"),]

dfpnpm<-data.frame(cat_tx,cat_biot,ORFs_tx$ORFs_pM)
#dfpnpm$value[dfpnpm$value==0]<-NA
dfpnpm$cat_tx<-factor(dfpnpm$cat_tx,levels=levvs)
dfpnpm$cat_biot<-factor(dfpnpm$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","novel isoform","non-coding RNA", "novel gene"))

dfwid<-data.frame(cat_tx,cat_biot,width(ORFs_tx))
dfwid$cat_tx<-factor(dfwid$cat_tx,levels=levvs)
dfwid$cat_biot<-factor(dfwid$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","novel isoform","non-coding RNA", "novel gene"))
dfwid<-dfwid[grep(dfwid$cat_tx,pattern="novel_tx"),]

dfwid$cat_tx<-factor(as.character(dfwid$cat_tx),levels=c("novel_tx_R1","novel_tx_R2","novel_tx_R1R2"))

df<-df[!is.na(df$value),]
df<-df[grep(df$cat_tx,pattern="novel_tx"),]
df$cat_tx<-factor(as.character(df$cat_tx),levels=c("novel_tx_R1","novel_tx_R2","novel_tx_R1R2"))
a<-ggplot(df,aes(x=cat_tx,y=value,fill=cat_tx,group=cat_tx))
a<-a + geom_bar(stat="identity",position = "dodge")
a<-a + facet_grid(. ~ cat_biot ,drop = T,scales = "free_x", space = "free_x")
a<-a + geom_text(aes(x=cat_tx, y=value, hjust=.5,vjust=-.2,label=value),check_overlap = TRUE,colour="black",position = position_dodge(width = .9),size=5)
a<-a + theme_bw()
a<-a + ylab(expression("Number of ORFs"))
a<-a + xlab("")
a<-a + scale_fill_manual(values = colsi2,"Annotation")
a<-a + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
a<-a + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=20))
orfs_found_n<-a + theme(strip.text.x = element_text(size=13, face="bold"),strip.text.y = element_text(size=20), legend.text=element_text(size=20), legend.title=element_text(size=23), panel.grid = element_blank()) + scale_y_log10(limits = c(1, max(df$value, na.rm = TRUE) * 1.1))


d<-ggplot(dfwid,aes(x=cat_tx,y=width.ORFs_tx.,fill=cat_tx))
d<-d + geom_violin(scale="width",draw_quantiles=.5,adjust = 1)
#c<-c + geom_jitter(aes(x=cat_tx,y=ORFs_tx.ORFs_pM),position = position_jitterdodge(jitter.width = .5), alpha = 0.2)
d<-d + scale_y_log10(breaks=c(10,100,1000,10000),labels=c(10,100,1000,10000),limits=c(min(dfwid$width.ORFs_tx.),max(dfwid$width.ORFs_tx.)*1.1))
#c<-c + geom_text(aes(x=cat_tx, y=ORFs_tx.ORFs_pM, hjust="top",label=ORFs_tx.ORFs_pM),colour="black",position = position_dodge(width = 1),size=5)
d<-d + theme_bw()
d<-d + facet_grid(. ~ cat_biot ,drop = T,scales = "free_x", space = "free_x")
d<-d + xlab("")
d<-d + ylab("ORF length (nt)")
d<-d + theme(axis.title.x = element_blank(),axis.text.x = element_blank())
#c<-c + theme(legend.position="none")
d<-d + scale_fill_manual(values = colsi2,"Annotation")
d<-d + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
d<-d + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=20))
orfs_found_len<-d + theme(strip.text.x = element_text(size=13, face="bold"),strip.text.y = element_text(size=20), legend.text=element_text(size=20), legend.title=element_text(size=23))

#orfs_found_len<-orfs_found_len + stat_compare_means(ref.group = "novel_tx_R1",label="p.signif",method.args=list(alternative="g"))


d<-ggplot(dfiso,aes(x=cat_tx,y=ORFs_tx.ORF_pct_P_sites,fill=cat_tx))
d<-d + geom_violin(scale="width",draw_quantiles=.5,adjust = 1)
#c<-c + geom_jitter(aes(x=cat_tx,y=ORFs_tx.ORFs_pM),position = position_jitterdodge(jitter.width = .5), alpha = 0.2)
#d<-d + scale_y_log10(breaks=c(10,100,1000,10000),labels=c(10,100,1000,10000),limits=c(min(dfwid$width.ORFs_tx.),max(dfwid$width.ORFs_tx.)*1.1))
#c<-c + geom_text(aes(x=cat_tx, y=ORFs_tx.ORFs_pM, hjust="top",label=ORFs_tx.ORFs_pM),colour="black",position = position_dodge(width = 1),size=5)
d<-d + theme_bw()
d<-d + facet_grid(. ~ cat_biot ,drop = T,scales = "free_x", space = "free_x")
d<-d + xlab("")
d<-d + ylab("ORF_pct_Psites")
d<-d + theme(axis.title.x = element_blank(),axis.text.x = element_blank())
#c<-c + theme(legend.position="none")
d<-d + scale_fill_manual(values = colsi2,"Annotation")
d<-d + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
d<-d + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=20))
orfs_pct<-d + theme(strip.text.x = element_text(size=13, face="bold"),strip.text.y = element_text(size=20), legend.text=element_text(size=20), legend.title=element_text(size=23))

#orfs_pct<-orfs_pct + stat_compare_means(ref.group = "novel_tx_R1",label="p.signif",method.args=list(alternative="g"))


# pdf(file = "Figures/Fig5_ORFs.pdf",width = 10,height = 4)
# orfs_found_n
# dev.off()

pdf(file = "Figures/SupplementaryFigure12.pdf",width = 10,height = 7)
plot_grid(orfs_found_n,orfs_found_len,orfs_pct,nrow = 3,rel_widths = c(1,1),align = "hv")
dev.off()

colsi3<-c("#76B7B2", "#59A14F","#EDC948")
a<-ggplot(df,aes(x=cat_biot ,y=value,fill=cat_biot ,group=cat_tx))
a<-a + geom_bar(stat="identity",position = "dodge")
a<-a + facet_grid(. ~ cat_tx ,drop = T,scales = "free_x", space = "free_x")
a<-a + geom_text(aes(x=cat_biot, y=value, hjust=.5,vjust=-.2,label=value),check_overlap = TRUE,colour="black",position = position_dodge(width = .9),size=5)
a<-a + theme_bw()
a<-a + ylab(expression("Number of ORFs"))
a<-a + xlab("")
a<-a + scale_fill_manual(values = colsi3,"")
a<-a + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
a<-a + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=20))
orfs_found_n<-a + theme(strip.text.x = element_text(size=13, face="bold"),strip.text.y = element_text(size=20), legend.text=element_text(size=20), legend.title=element_text(size=23),legend.position="top", panel.grid = element_blank()) + scale_y_log10(limits = c(1, max(df$value, na.rm = TRUE) * 1.1))


d<-ggplot(dfwid,aes(x=cat_biot ,y=width.ORFs_tx.,fill=cat_biot ))
d<-d + geom_violin(scale="width",draw_quantiles=.5,adjust = 1)
#c<-c + geom_jitter(aes(x=cat_tx,y=ORFs_tx.ORFs_pM),position = position_jitterdodge(jitter.width = .5), alpha = 0.2)
d<-d + scale_y_log10(breaks=c(10,100,1000,10000),labels=c(10,100,1000,10000),limits=c(min(dfwid$width.ORFs_tx.),max(dfwid$width.ORFs_tx.)*1.1))
#c<-c + geom_text(aes(x=cat_tx, y=ORFs_tx.ORFs_pM, hjust="top",label=ORFs_tx.ORFs_pM),colour="black",position = position_dodge(width = 1),size=5)
d<-d + theme_bw()
d<-d + facet_grid(. ~ cat_tx ,drop = T,scales = "free_x", space = "free_x")
d<-d + xlab("")
d<-d + ylab("ORF length (nt)")
d<-d + theme(axis.title.x = element_blank(),axis.text.x = element_blank())
#c<-c + theme(legend.position="none")
d<-d + scale_fill_manual(values = colsi3,"")
d<-d + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
d<-d + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=20))
orfs_found_len<-d + theme(strip.text.x = element_text(size=13, face="bold"),strip.text.y = element_text(size=20), legend.text=element_text(size=20), legend.title=element_text(size=23),legend.position="top")


d<-ggplot(dfiso,aes(x=cat_biot ,y=ORFs_tx.ORF_pct_P_sites,fill=cat_biot ))
d<-d + geom_violin(scale="width",draw_quantiles=.5,adjust = 1)
#c<-c + geom_jitter(aes(x=cat_tx,y=ORFs_tx.ORFs_pM),position = position_jitterdodge(jitter.width = .5), alpha = 0.2)
#d<-d + scale_y_log10(breaks=c(10,100,1000,10000),labels=c(10,100,1000,10000),limits=c(min(dfwid$width.ORFs_tx.),max(dfwid$width.ORFs_tx.)*1.1))
#c<-c + geom_text(aes(x=cat_tx, y=ORFs_tx.ORFs_pM, hjust="top",label=ORFs_tx.ORFs_pM),colour="black",position = position_dodge(width = 1),size=5)
d<-d + theme_bw()
d<-d + facet_grid(. ~ cat_tx ,drop = T,scales = "free_x", space = "free_x")
d<-d + xlab("")
d<-d + ylab("ORF_pct_Psites")
d<-d + theme(axis.title.x = element_blank(),axis.text.x = element_blank())
#c<-c + theme(legend.position="none")
d<-d + scale_fill_manual(values = colsi3,"")
d<-d + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
d<-d + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=20))
orfs_pct<-d + theme(strip.text.x = element_text(size=13, face="bold"),strip.text.y = element_text(size=20), legend.text=element_text(size=20), legend.title=element_text(size=23),legend.position="top")

#orfs_found_len<-orfs_found_len + stat_compare_means(ref.group = "novel_tx_R1",label="p.signif",method.args=list(alternative="g"))

# pdf(file = "Figures/Fig5_ORFs2.pdf",width = 10,height = 4)
# orfs_found_n
# dev.off()

# pdf(file = "Figures/Fig5_ORFs2_more.pdf",width = 10,height = 7)
# plot_grid(orfs_found_n,orfs_found_len,orfs_pct,nrow = 3,rel_widths = c(1,1,1),align = "hv")
# dev.off()



ORFquant_results<-orfs_srlr_all
ORFs_tx<-ORFquant_results$ORFs_tx
ORFs_gen<-ORFquant_results$ORFs_gen
ORFs_gen<-split(ORFs_gen,names(ORFs_gen))

#make R2 instead of PB when problem is solved
ORFs_tx$tx_novel<-sum(!grepl(ORFs_tx$compatible_with,pattern="^R1|^R2"))==0
ORFs_tx$tx_novel_R1R2<-"none"
ORFs_tx$tx_novel_R1R2[ORFs_tx$tx_novel]<-"R1R2"
ORFs_tx$tx_novel_R1R2[sum(!grepl(ORFs_tx$compatible_with,pattern="^R1"))==0]<-"R1"
ORFs_tx$tx_novel_R1R2[sum(!grepl(ORFs_tx$compatible_with,pattern="^R2"))==0]<-"R2"

ORFs_tx$compatible_biotype[ORFs_tx$tx_novel]<-paste("novel_tx",ORFs_tx$tx_novel_R1R2[ORFs_tx$tx_novel],sep="_")

gids<-ORFs_tx$gene_id
cat_biot<-ORFs_tx$gene_biotype
cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
cat_biot[!cat_biot%in%c("protein_coding","pseudogene")]<-"non-coding RNA"
tbid<-table(gids)
tbid[tbid>3]<-">3"
#boxplot(log(tpms)~tbid)
tb_bio<-cat_biot[match(names(tbid),gids)]


cat_tx<-ORFs_tx$ORF_category_Tx_compatible
cat_tx[cat_tx=="novel"]<-"no_CDS_info"
cat_tx[cat_tx=="overl_uORF"]<-"uORF"
cat_tx[cat_tx=="overl_dORF"]<-"other"
cat_tx[cat_tx=="dORF"]<-"other"
cat_tx[grep(ORFs_tx$compatible_biotype,pattern="^novel_tx_")]<-ORFs_tx$compatible_biotype[grep(ORFs_tx$compatible_biotype,pattern="^novel_tx_")]

cat_gen<-ORFs_tx$ORF_category_Gen
cat_biot<-ORFs_tx$compatible_biotype
cat_biot_gene<-ORFs_tx$gene_biotype
cat_biot[grepl(cat_biot, pattern="^novel_tx_") & cat_biot_gene=="protein_coding"] <- "novel isoform"
cat_biot[!cat_biot %in% c("protein_coding", "novel isoform") & cat_biot_gene=="protein_coding"]<-"non-coding isoform"
cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
cat_biot[!cat_biot%in%c("protein_coding","pseudogene","non-coding isoform", "novel isoform")]<-"non-coding RNA"
cat_biot[grepl(cat_biot_gene, pattern="^novel")] <- "novel gene"

#cat_tx[intersect(grep(ORFs_tx$gene_biotype,pattern = "pseudogene"),which(cat_tx=="novel"))]<-"novel_pseudogene"
cat_tx[cat_tx%in%c("C_extension","C_truncation","NC_extension","N_extension","nested_ORF","overl_dORF")]<-"other"

cat_gen[grep(cat_gen,pattern = "novel")]<-"non-overlapping\nCDS regions"
cat_gen[grep(cat_gen,pattern = "non-overlapping\nCDS regions",invert = T)]<-"overlapping\nCDS regions"
levvs<-c("ORF_annotated","N_truncation","uORF","other","no_CDS_info","novel_tx_R1","novel_tx_R2","novel_tx_R1R2")

df<-melt(table(cat_tx,cat_biot))

df$value[df$value==0]<-NA
df$cat_tx<-factor(df$cat_tx,levels=levvs)
df$cat_biot<-factor(df$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","non-coding RNA", "novel isoform", "novel gene"))


colli <- c("#4E79A7",  "#F28E2B", "#E15759", "#76B7B2", "#59A14F","#EDC948")  

df<-df[!is.na(df$value),]
#df<-df[grep(df$cat_tx,pattern="novel_tx"),]
#df$cat_tx<-factor(as.character(df$cat_tx),levels=c("novel_tx_R1","novel_tx_R2","novel_tx_R1R2"))
df$reads<-"all"
df_all<-df


ORFquant_results<-orfs_srlr_unq
ORFs_tx<-ORFquant_results$ORFs_tx
ORFs_gen<-ORFquant_results$ORFs_gen
ORFs_gen<-split(ORFs_gen,names(ORFs_gen))

#make R2 instead of PB when problem is solved
ORFs_tx$tx_novel<-sum(!grepl(ORFs_tx$compatible_with,pattern="^R1|^R2"))==0
ORFs_tx$tx_novel_R1R2<-"none"
ORFs_tx$tx_novel_R1R2[ORFs_tx$tx_novel]<-"R1R2"
ORFs_tx$tx_novel_R1R2[sum(!grepl(ORFs_tx$compatible_with,pattern="^R1"))==0]<-"R1"
ORFs_tx$tx_novel_R1R2[sum(!grepl(ORFs_tx$compatible_with,pattern="^R2"))==0]<-"R2"

ORFs_tx$compatible_biotype[ORFs_tx$tx_novel]<-paste("novel_tx",ORFs_tx$tx_novel_R1R2[ORFs_tx$tx_novel],sep="_")

gids<-ORFs_tx$gene_id
cat_biot<-ORFs_tx$gene_biotype
cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
cat_biot[!cat_biot%in%c("protein_coding","pseudogene")]<-"non-coding RNA"
tbid<-table(gids)
tbid[tbid>3]<-">3"
#boxplot(log(tpms)~tbid)
tb_bio<-cat_biot[match(names(tbid),gids)]


cat_tx<-ORFs_tx$ORF_category_Tx_compatible
cat_tx[cat_tx=="novel"]<-"no_CDS_info"
cat_tx[cat_tx=="overl_uORF"]<-"uORF"
cat_tx[cat_tx=="overl_dORF"]<-"other"
cat_tx[cat_tx=="dORF"]<-"other"
cat_tx[grep(ORFs_tx$compatible_biotype,pattern="^novel_tx_")]<-ORFs_tx$compatible_biotype[grep(ORFs_tx$compatible_biotype,pattern="^novel_tx_")]

cat_gen<-ORFs_tx$ORF_category_Gen
cat_biot<-ORFs_tx$compatible_biotype
cat_biot_gene<-ORFs_tx$gene_biotype
cat_biot[grepl(cat_biot, pattern="^novel_tx_") & cat_biot_gene=="protein_coding"] <- "novel isoform"
cat_biot[!cat_biot %in% c("protein_coding", "novel isoform") & cat_biot_gene=="protein_coding"]<-"non-coding isoform"
cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
cat_biot[!cat_biot%in%c("protein_coding","pseudogene","non-coding isoform", "novel isoform")]<-"non-coding RNA"
cat_biot[grepl(cat_biot_gene, pattern="^novel")] <- "novel gene"

#cat_tx[intersect(grep(ORFs_tx$gene_biotype,pattern = "pseudogene"),which(cat_tx=="novel"))]<-"novel_pseudogene"
cat_tx[cat_tx%in%c("C_extension","C_truncation","NC_extension","N_extension","nested_ORF","overl_dORF")]<-"other"

cat_gen[grep(cat_gen,pattern = "novel")]<-"non-overlapping\nCDS regions"
cat_gen[grep(cat_gen,pattern = "non-overlapping\nCDS regions",invert = T)]<-"overlapping\nCDS regions"
levvs<-c("ORF_annotated","N_truncation","uORF","other","no_CDS_info","novel_tx_R1","novel_tx_R2","novel_tx_R1R2")

df<-melt(table(cat_tx,cat_biot))

df$value[df$value==0]<-NA
df$cat_tx<-factor(df$cat_tx,levels=levvs)
df$cat_biot<-factor(df$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","non-coding RNA", "novel isoform", "novel gene"))


colli <- c("#4E79A7",  "#F28E2B", "#E15759", "#76B7B2", "#59A14F","#EDC948")  

df<-df[!is.na(df$value),]
#df<-df[grep(df$cat_tx,pattern="novel_tx"),]
#df$cat_tx<-factor(as.character(df$cat_tx),levels=c("novel_tx_R1","novel_tx_R2","novel_tx_R1R2"))
df$reads<-"unique"
df_unq<-df
df<-rbind(df_unq,df_all)

df<-df[grep(df$cat_tx,pattern="novel_tx"),]
df$cat_tx<-factor(as.character(df$cat_tx),levels=c("novel_tx_R1","novel_tx_R2","novel_tx_R1R2"))


colsi3<-c("#76B7B2", "#59A14F","#EDC948")
a<-ggplot(df,aes(x=cat_biot ,y=value,fill=cat_biot ,group=cat_tx))
a<-a + geom_bar(stat="identity",position = "dodge")
a<-a + facet_grid(reads ~ cat_tx ,drop = T,scales = "free_x", space = "free_x")
a<-a + geom_text(aes(x=cat_biot, y=value, hjust=.5,vjust=-.2,label=value),check_overlap = TRUE,colour="black",position = position_dodge(width = .9),size=5)
a<-a + theme_bw()
a<-a + ylab(expression("Number of ORFs"))
a<-a + xlab("")
a<-a + scale_fill_manual(values = colsi3,"")
a<-a + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
a<-a + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=20))
orfs_found_n<-a + theme(strip.text.x = element_text(size=13, face="bold"),strip.text.y = element_text(size=20), legend.text=element_text(size=20), legend.title=element_text(size=23), panel.grid = element_blank()) + scale_y_log10(limits = c(1, max(df$value, na.rm = TRUE) * 1.1))

# pdf(file = "Figures/Fig5_ORFs_unqall.pdf",width = 12,height = 7)
# orfs_found_n
# dev.off()


ORFquant_results<-orfs_sr_all
ORFs_tx<-ORFquant_results$ORFs_tx
ORFs_gen<-ORFquant_results$ORFs_gen
ORFs_gen<-split(ORFs_gen,names(ORFs_gen))

#make R2 instead of PB when problem is solved
ORFs_tx$tx_novel<-sum(!grepl(ORFs_tx$compatible_with,pattern="^R1|^R2"))==0
ORFs_tx$tx_novel_R1R2<-"none"
ORFs_tx$tx_novel_R1R2[ORFs_tx$tx_novel]<-"R1R2"
ORFs_tx$tx_novel_R1R2[sum(!grepl(ORFs_tx$compatible_with,pattern="^R1"))==0]<-"R1"
ORFs_tx$tx_novel_R1R2[sum(!grepl(ORFs_tx$compatible_with,pattern="^R2"))==0]<-"R2"

ORFs_tx$compatible_biotype[ORFs_tx$tx_novel]<-paste("novel_tx",ORFs_tx$tx_novel_R1R2[ORFs_tx$tx_novel],sep="_")

gids<-ORFs_tx$gene_id
cat_biot<-ORFs_tx$gene_biotype
cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
cat_biot[!cat_biot%in%c("protein_coding","pseudogene")]<-"non-coding RNA"
tbid<-table(gids)
tbid[tbid>3]<-">3"
#boxplot(log(tpms)~tbid)
tb_bio<-cat_biot[match(names(tbid),gids)]


cat_tx<-ORFs_tx$ORF_category_Tx_compatible
cat_tx[cat_tx=="novel"]<-"no_CDS_info"
cat_tx[cat_tx=="overl_uORF"]<-"uORF"
cat_tx[cat_tx=="overl_dORF"]<-"other"
cat_tx[cat_tx=="dORF"]<-"other"
cat_tx[grep(ORFs_tx$compatible_biotype,pattern="^novel_tx_")]<-ORFs_tx$compatible_biotype[grep(ORFs_tx$compatible_biotype,pattern="^novel_tx_")]

cat_gen<-ORFs_tx$ORF_category_Gen
cat_biot<-ORFs_tx$compatible_biotype
cat_biot_gene<-ORFs_tx$gene_biotype
cat_biot[grepl(cat_biot, pattern="^novel_tx_") & cat_biot_gene=="protein_coding"] <- "novel isoform"
cat_biot[!cat_biot %in% c("protein_coding", "novel isoform") & cat_biot_gene=="protein_coding"]<-"non-coding isoform"
cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
cat_biot[!cat_biot%in%c("protein_coding","pseudogene","non-coding isoform", "novel isoform")]<-"non-coding RNA"
cat_biot[grepl(cat_biot_gene, pattern="^novel")] <- "novel gene"

#cat_tx[intersect(grep(ORFs_tx$gene_biotype,pattern = "pseudogene"),which(cat_tx=="novel"))]<-"novel_pseudogene"
cat_tx[cat_tx%in%c("C_extension","C_truncation","NC_extension","N_extension","nested_ORF","overl_dORF")]<-"other"

cat_gen[grep(cat_gen,pattern = "novel")]<-"non-overlapping\nCDS regions"
cat_gen[grep(cat_gen,pattern = "non-overlapping\nCDS regions",invert = T)]<-"overlapping\nCDS regions"
levvs<-c("ORF_annotated","N_truncation","uORF","other","no_CDS_info","novel_tx_R1","novel_tx_R2","novel_tx_R1R2")

df<-melt(table(cat_tx,cat_biot))

df$value[df$value==0]<-NA
df$cat_tx<-factor(df$cat_tx,levels=levvs)
df$cat_biot<-factor(df$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","non-coding RNA", "novel isoform", "novel gene"))


colli <- c("#4E79A7",  "#F28E2B", "#E15759", "#76B7B2", "#59A14F","#EDC948")  

df<-df[!is.na(df$value),]
#df<-df[grep(df$cat_tx,pattern="novel_tx"),]
#df$cat_tx<-factor(as.character(df$cat_tx),levels=c("novel_tx_R1","novel_tx_R2","novel_tx_R1R2"))
df$reads<-"all"
df_all<-df


ORFquant_results<-orfs_sr_unq
ORFs_tx<-ORFquant_results$ORFs_tx
ORFs_gen<-ORFquant_results$ORFs_gen
ORFs_gen<-split(ORFs_gen,names(ORFs_gen))

#make R2 instead of PB when problem is solved
ORFs_tx$tx_novel<-sum(!grepl(ORFs_tx$compatible_with,pattern="^R1|^R2"))==0
ORFs_tx$tx_novel_R1R2<-"none"
ORFs_tx$tx_novel_R1R2[ORFs_tx$tx_novel]<-"R1R2"
ORFs_tx$tx_novel_R1R2[sum(!grepl(ORFs_tx$compatible_with,pattern="^R1"))==0]<-"R1"
ORFs_tx$tx_novel_R1R2[sum(!grepl(ORFs_tx$compatible_with,pattern="^R2"))==0]<-"R2"

ORFs_tx$compatible_biotype[ORFs_tx$tx_novel]<-paste("novel_tx",ORFs_tx$tx_novel_R1R2[ORFs_tx$tx_novel],sep="_")

gids<-ORFs_tx$gene_id
cat_biot<-ORFs_tx$gene_biotype
cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
cat_biot[!cat_biot%in%c("protein_coding","pseudogene")]<-"non-coding RNA"
tbid<-table(gids)
tbid[tbid>3]<-">3"
#boxplot(log(tpms)~tbid)
tb_bio<-cat_biot[match(names(tbid),gids)]


cat_tx<-ORFs_tx$ORF_category_Tx_compatible
cat_tx[cat_tx=="novel"]<-"no_CDS_info"
cat_tx[cat_tx=="overl_uORF"]<-"uORF"
cat_tx[cat_tx=="overl_dORF"]<-"other"
cat_tx[cat_tx=="dORF"]<-"other"
cat_tx[grep(ORFs_tx$compatible_biotype,pattern="^novel_tx_")]<-ORFs_tx$compatible_biotype[grep(ORFs_tx$compatible_biotype,pattern="^novel_tx_")]

cat_gen<-ORFs_tx$ORF_category_Gen
cat_biot<-ORFs_tx$compatible_biotype
cat_biot_gene<-ORFs_tx$gene_biotype
cat_biot[grepl(cat_biot, pattern="^novel_tx_") & cat_biot_gene=="protein_coding"] <- "novel isoform"
cat_biot[!cat_biot %in% c("protein_coding", "novel isoform") & cat_biot_gene=="protein_coding"]<-"non-coding isoform"
cat_biot[grep(cat_biot,pattern = "pseudo")]<-"pseudogene"
cat_biot[!cat_biot%in%c("protein_coding","pseudogene","non-coding isoform", "novel isoform")]<-"non-coding RNA"
cat_biot[grepl(cat_biot_gene, pattern="^novel")] <- "novel gene"

#cat_tx[intersect(grep(ORFs_tx$gene_biotype,pattern = "pseudogene"),which(cat_tx=="novel"))]<-"novel_pseudogene"
cat_tx[cat_tx%in%c("C_extension","C_truncation","NC_extension","N_extension","nested_ORF","overl_dORF")]<-"other"

cat_gen[grep(cat_gen,pattern = "novel")]<-"non-overlapping\nCDS regions"
cat_gen[grep(cat_gen,pattern = "non-overlapping\nCDS regions",invert = T)]<-"overlapping\nCDS regions"
levvs<-c("ORF_annotated","N_truncation","uORF","other","no_CDS_info","novel_tx_R1","novel_tx_R2","novel_tx_R1R2")

df<-melt(table(cat_tx,cat_biot))

df$value[df$value==0]<-NA
df$cat_tx<-factor(df$cat_tx,levels=levvs)
df$cat_biot<-factor(df$cat_biot,levels=c("protein_coding","non-coding isoform","pseudogene","non-coding RNA", "novel isoform", "novel gene"))


colli <- c("#4E79A7",  "#F28E2B", "#E15759", "#76B7B2", "#59A14F","#EDC948")  

df<-df[!is.na(df$value),]
#df<-df[grep(df$cat_tx,pattern="novel_tx"),]
#df$cat_tx<-factor(as.character(df$cat_tx),levels=c("novel_tx_R1","novel_tx_R2","novel_tx_R1R2"))
df$reads<-"unique"
df_unq<-df
df<-rbind(df_unq,df_all)

#df$cat_tx<-factor(as.character(df$cat_tx),levels=c("novel_tx_R1","novel_tx_R2","novel_tx_R1R2"))


colsi3<-c("#76B7B2", "#59A14F","#EDC948")
colli <- c("#4E79A7",  "#F28E2B", "#E15759", "#76B7B2", "#59A14F","#EDC948")
a<-ggplot(df,aes(x=cat_biot ,y=value,fill=cat_biot ,group=cat_tx))
a<-a + geom_bar(stat="identity",position = "dodge")
a<-a + facet_grid(reads ~ cat_tx ,drop = T,scales = "free_x", space = "free_x")
a<-a + geom_text(aes(x=cat_biot, y=value, hjust=.5,vjust=-.2,label=value),check_overlap = TRUE,colour="black",position = position_dodge(width = .9),size=5)
a<-a + theme_bw()
a<-a + ylab(expression("Number of ORFs"))
a<-a + xlab("")
a<-a + scale_fill_manual(values = colli,"")
a<-a + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
a<-a + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=20))
orfs_found_n<-a + theme(strip.text.x = element_text(size=13, face="bold"),strip.text.y = element_text(size=20), legend.text=element_text(size=20), legend.title=element_text(size=23), panel.grid = element_blank()) + scale_y_log10(limits = c(1, max(df$value, na.rm = TRUE) * 1.1))

pdf(file = "Figures/SupplementaryFigure7.pdf",width = 16,height = 7)
orfs_found_n
dev.off()

