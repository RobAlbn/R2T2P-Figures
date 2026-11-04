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

lens_tx<-data.frame(sum(width(GTF_annotation$exons_txs)))

trann_all<-GTF_annotation$trann
trann_all$lentx<-lens_tx[match(trann_all$transcript_id,rownames(lens_tx)),1]


cat_biot_gene<-trann_all$gene_biotype
cat_biot_gene[grep(cat_biot_gene,pattern = "pseudo")]<-"pseudogene"
#cat_biot_gene[grep(cat_biot_gene,pattern = "^novel")]<-"novel gene"
cat_biot_gene[!cat_biot_gene%in%c("protein_coding","pseudogene","lncRNA","novel_antisense","novel_antisense_intr","novel_intergenic")]<-"non-coding_RNA"

notok<-cat_biot_gene=="non-coding_RNA"
trann_all<-trann_all[!notok,]
cat_biot_gene<-cat_biot_gene[!notok]

cat_biot<-trann_all$transcript_biotype
cat_biot[grepl(cat_biot,patter="protein_coding")]<-"protein_coding"
cat_biot[grepl(cat_biot, pattern="^novel") & cat_biot_gene=="protein_coding"] <- "novel_isoform"

cat_biot[grepl(cat_biot, pattern="^novel") & !cat_biot_gene%in%c("protein_coding","novel_antisense","novel_antisense_intr","novel_intergenic")] <- "novel_non-coding"
cat_biot[!grepl(cat_biot, pattern="^novel") & !cat_biot_gene%in%c("protein_coding","novel_antisense","novel_antisense_intr","novel_intergenic")] <- "non-coding_RNA"
cat_biot[!cat_biot %in% c("protein_coding", "novel_isoform") & cat_biot_gene=="protein_coding"]<-"non-coding_isoform"
cat_biot[cat_biot %in% c("novel_antisense","novel_antisense_intr","novel_intergenic")]<-"novel_tx"
cat_biot[cat_biot=="non-coding_RNA" & cat_biot_gene=="pseudogene"]<-"pseudogene"
cat_biot[cat_biot=="non-coding_RNA" & cat_biot_gene=="lncRNA"]<-"lncRNA"

#add R1 and R2 and make the two figures

catunq<-unique(trann_all[,c("gene_id","gene_biotype")])
catunq<-catunq[grep(catunq$gene_biotype,pattern="^novel"),"gene_biotype"]

dfgn<-data.frame(melt(table(catunq)))
coli<-c("#0072B2", "#56B4E9", "#009E73")
names(coli)<-c("novel_antisense","novel_antisense_intronic","novel_intergenic")
dfgn$catunq<-as.character(dfgn$catunq)
#dfgn$catunq<-gsub(dfgn$catunq,pattern="novel_",replacement="")
#dfgn$catunq<-gsub(dfgn$catunq,pattern="_",replacement="\n")
dfgn$catunq<-gsub(dfgn$catunq,pattern="intr",replacement="intronic")
dfgn$catunq<-factor(dfgn$catunq,levels=c("novel_intergenic","novel_antisense","novel_antisense_intronic"))

a<-ggplot(dfgn,aes(x=catunq,y=value,fill=catunq))
a<-a + geom_bar(stat="identity",position = "dodge")
a<-a + geom_text(aes(x=catunq, y=value, hjust=.5,vjust=-.2,label=value),check_overlap = TRUE,colour="black",position = position_dodge(width = .9),size=5)
a<-a + theme_bw()
a<-a + ylab(expression("Novel genes"))
a<-a + xlab("")
a<-a + scale_fill_manual(values = coli,"Biotype")
a<-a + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
#a<-a + theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(angle=45, vjust=0.5, size=20))
a<-a + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=20))
gn_found_n<-a + theme(strip.text.x = element_text(size=15, face="bold"),strip.text.y = element_text(size=20), legend.text=element_text(size=14), legend.title=element_text(size=15),legend.position="top") 
  
# pdf(file = "Figures/Fig2_genes.pdf",width = 7,height = 5)
# gn_found_n
# dev.off()

dfgn$catunq<-as.character(dfgn$catunq)
dfgn$catunq<-gsub(dfgn$catunq,pattern="novel_",replacement="")
dfgn$catunq<-gsub(dfgn$catunq,pattern="_",replacement="\n")
dfgn$catunq<-factor(dfgn$catunq,levels=c("intergenic","antisense","antisense\nintronic"))
coli22<-c("#0072B2", "#56B4E9", "#009E73")
names(coli22)<-c("antisense","antisense\nintronic","intergenic")

a<-ggplot(dfgn,aes(x=catunq,y=value,fill=catunq))
a<-a + geom_bar(stat="identity",position = "dodge")
a<-a + geom_text(aes(x=catunq, y=value, hjust=.5,vjust=-.2,label=value),check_overlap = TRUE,colour="black",position = position_dodge(width = .9),size=5)
a<-a + theme_bw()
a<-a + ylab(expression("Novel genes"))
a<-a + xlab("")
a<-a + scale_fill_manual(values = coli22,"")
a<-a + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
#a<-a + theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(angle=45, vjust=0.5, size=20))
a<-a + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=20))
gn_found_n<-a + theme(strip.text.x = element_text(size=15, face="bold"),strip.text.y = element_text(size=20), legend.text=element_text(size=14), legend.title=element_text(size=15), legend.position="top", panel.grid = element_blank()) 

# On November 3rd at 10:46 PM I made this plot shorter and wider; at 10:58 I specified "top" as legend position (previous line)
# On November 4th at 8:54 PM I set width to 5 and height to 4
pdf(file = "Figures/Figure2_B.pdf",width = 5,height = 4)
gn_found_n
dev.off()

cat_biot_orig<-cat_biot
cat_biot_gene_orig<-cat_biot_gene

cat_biot_gene<-cat_biot_gene_orig
cat_biot_gene[grep(cat_biot_gene,pattern="novel")]<-"novel_gene"

cat_biot<-cat_biot_orig
cat_biot[cat_biot=="protein_coding"]<-"annotated_coding"
cat_biot[grep(cat_biot,pattern="^novel")]<-"novel_tx"
cat_biot[!cat_biot%in%c("annotated_coding","novel_tx")]<-"annotated_non-coding"

dftx<-data.frame(melt(table(cat_biot,cat_biot_gene)))

dftx$cat_biot<-factor(dftx$cat_biot,levels=c("annotated_coding","annotated_non-coding","novel_tx"))
dftx$cat_biot_gene<-factor(dftx$cat_biot_gene,levels=c("protein_coding","lncRNA","pseudogene","non-coding_RNA","novel_gene"))
dftx<-dftx[dftx$value>0,]

colli <- c("black","grey","red")  
a<-ggplot(dftx,aes(x=cat_biot,y=value,fill=cat_biot))
a<-a + geom_bar(stat="identity",position = "dodge")

a<-a + facet_grid(. ~ cat_biot_gene,drop = T,scales = "free_x", space = "free_x")
a<-a + geom_text(aes(x=cat_biot, y=value, hjust=.5,vjust=-.2,label=value),check_overlap = TRUE,colour="black",position = position_dodge(width = .9),size=5)
a<-a + theme_bw()
a<-a + ylab(expression("Number of txs"))
a<-a + xlab("")
a<-a + scale_fill_manual(values = colli,"")
a<-a + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
a<-a + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=20))
tx_found_n<-a + theme(strip.text.x = element_text(size=15, face="bold"),strip.text.y = element_text(size=20), legend.text=element_text(size=20), legend.title=element_text(size=23), panel.grid = element_blank()) + scale_y_log10(limits = c(1, max(dftx$value, na.rm = TRUE) * 1.1))
pdf(file = "Figures/SupplementaryFigure2_top.pdf",width = 15,height = 5)
tx_found_n
dev.off()

dfwid<-data.frame(cat_biot,cat_biot_gene,trann_all$lentx)
dfwid$cat_biot<-factor(dfwid$cat_biot,levels=c("annotated_coding","annotated_non-coding","novel_tx"))
dfwid$cat_biot_gene<-factor(dfwid$cat_biot_gene,levels=c("protein_coding","lncRNA","pseudogene","non-coding_RNA","novel_gene"))

dfwid$len<-dfwid$trann_all.lentx

d<-ggplot(dfwid,aes(x=cat_biot,y=len,fill=cat_biot))
d<-d + geom_violin(scale="width",draw_quantiles=.5,adjust = 1)
d<-d + scale_y_log10(breaks=c(10,100,1000,10000,100000),labels=c(10,100,1000,10000,100000),limits=c(min(dfwid$len),max(dfwid$len)*1.1))
d<-d + theme_bw()
d<-d + facet_grid(. ~ cat_biot_gene ,drop = T,scales = "free_x",space = "free_x")
d<-d + xlab("")
d<-d + ylab("Tx length (nt)")
d<-d + theme(axis.title.x = element_blank(),axis.text.x = element_blank())
#c<-c + theme(legend.position="none")
d<-d + scale_fill_manual(values = colli,"")
d<-d + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
d<-d + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=20))
tx_found_len<-d + theme(strip.text.x = element_text(size=13, face="bold"),strip.text.y = element_text(size=20), legend.text=element_text(size=20), legend.title=element_text(size=23))

pdf(file = "Figures/SupplementaryFigure2_bottom.pdf",width = 15,height = 5)
tx_found_len
dev.off()


lens_tx<-data.frame(sum(width(GTF_all$exons_txs)))

trann_all<-GTF_all$trann
trann_all$lentx<-lens_tx[match(trann_all$transcript_id,rownames(lens_tx)),1]

trann_all$novelty_gene<-"none"
trann_all$novelty_gene[grep(trann_all$gene_id,pattern="^R1")]<-"R1"
trann_all$novelty_gene[grep(trann_all$gene_id,pattern="^R1R2")]<-"R1R2"
trann_all$novelty_gene[grep(trann_all$gene_id,pattern="^R2")]<-"R2"

trann_all$novelty_tx<-"none"
trann_all$novelty_tx[grep(trann_all$transcript_id,pattern="^R1")]<-"R1"
trann_all$novelty_tx[grep(trann_all$transcript_id,pattern="^R2")]<-"R2"

cat_biot_gene<-trann_all$gene_biotype
cat_biot_gene[grep(cat_biot_gene,pattern = "pseudo")]<-"pseudogene"
cat_biot_gene[!cat_biot_gene%in%c("protein_coding","pseudogene","lncRNA","novel_antisense","novel_antisense_intr","novel_intergenic")]<-"non-coding_RNA"

cat_biot<-trann_all$transcript_biotype
cat_biot[grepl(cat_biot,patter="protein_coding")]<-"protein_coding"
cat_biot[grepl(cat_biot, pattern="^novel") & cat_biot_gene=="protein_coding"] <- "novel_isoform"

cat_biot[grepl(cat_biot, pattern="^novel") & !cat_biot_gene%in%c("protein_coding","novel_antisense","novel_antisense_intr","novel_intergenic")] <- "novel_non-coding"
cat_biot[!grepl(cat_biot, pattern="^novel") & !cat_biot_gene%in%c("protein_coding","novel_antisense","novel_antisense_intr","novel_intergenic")] <- "non-coding_RNA"
cat_biot[!cat_biot %in% c("protein_coding", "novel_isoform") & cat_biot_gene=="protein_coding"]<-"non-coding_isoform"
cat_biot[cat_biot %in% c("novel_antisense","novel_antisense_intr","novel_intergenic")]<-"novel_tx"
cat_biot[cat_biot=="non-coding_RNA" & cat_biot_gene=="pseudogene"]<-"pseudogene"
cat_biot[cat_biot=="non-coding_RNA" & cat_biot_gene=="lncRNA"]<-"lncRNA"

#add R1 and R2 and make the two figures


notok<-cat_biot_gene=="non-coding_RNA"
trann_all<-trann_all[!notok,]
cat_biot_gene<-cat_biot_gene[!notok]
cat_biot<-cat_biot[!notok]

catunq<-unique(trann_all[,c("gene_id","gene_biotype","novelty_gene","novelty_tx")])
catunq<-catunq[grep(catunq$gene_biotype,pattern="^novel"),]
#catunq<-catunq[grep(catunq$gene_biotype,pattern="^novel"),"gene_biotype"]
dfgn<-data.frame(melt(table(catunq$gene_biotype,catunq$novelty_gene)))
coli<-c("#0072B2", "#56B4E9", "#009E73")
dfgn$catunq<-dfgn$Var1
dfgn$catunq<-as.character(dfgn$catunq)
dfgn$catunq<-gsub(dfgn$catunq,pattern="novel_",replacement="")
#dfgn$catunq<-gsub(dfgn$catunq,pattern="_",replacement="\n")
dfgn$catunq<-gsub(dfgn$catunq,pattern="intr",replacement="intronic")
dfgn$Var2<-factor(as.character(dfgn$Var2),levels=c("R1","R2","R1R2"))

a<-ggplot(dfgn,aes(x=catunq,y=value,fill=catunq))
a<-a + geom_bar(stat="identity",position = "dodge")
a<-a + geom_text(aes(x=catunq, y=value, hjust=.5,vjust=-.2,label=value),check_overlap = TRUE,colour="black",position = position_dodge(width = .9),size=5)
a<-a + theme_bw()
a<-a + ylab(expression("Novel genes"))
a<-a + xlab("")
a<-a + scale_fill_manual(values = coli,"Biotype")
a<-a + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
#a<-a + theme(axis.title.x = element_text(size=20),axis.text.x  = element_text(angle=45, vjust=0.5, size=20))
a<-a + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=20))
gn_found_n<-a + theme(strip.text.x = element_text(size=15, face="bold"),strip.text.y = element_text(size=20), legend.text=element_text(size=20), legend.title=element_text(size=23)) 

# pdf(file = "Figures/Fig5_genes.pdf",width = 10,height = 5)
# gn_found_n + facet_grid(. ~ Var2,drop = T,scales = "free_x", space = "free_x")
# dev.off()

cat_biot_orig<-cat_biot
cat_biot_gene_orig<-cat_biot_gene

cat_biot_gene<-cat_biot_gene_orig
cat_biot_gene[grep(cat_biot_gene,pattern="novel")]<-"novel_gene"

cat_biot<-cat_biot_orig
cat_biot[cat_biot=="protein_coding"]<-"annotated_coding"
cat_biot[grep(cat_biot,pattern="^novel")]<-"novel_tx"
cat_biot[!cat_biot%in%c("annotated_coding","novel_tx")]<-"annotated_non-coding"
cat_biot[cat_biot=="novel_tx"]<-paste(cat_biot[cat_biot=="novel_tx"],trann_all$novelty_tx[cat_biot=="novel_tx"],sep="_")

dftx<-data.frame(melt(table(cat_biot,cat_biot_gene)))

dftx$cat_biot<-factor(dftx$cat_biot,levels=c("annotated_coding","annotated_non-coding","novel_tx_R1","novel_tx_R2"))
dftx$cat_biot_gene<-factor(dftx$cat_biot_gene,levels=c("protein_coding","lncRNA","pseudogene","non-coding_RNA","novel_gene"))
dftx<-dftx[dftx$value>0,]

colli <- c("black","grey","red","red3")  
a<-ggplot(dftx,aes(x=cat_biot,y=value,fill=cat_biot))
a<-a + geom_bar(stat="identity",position = "dodge")

a<-a + facet_grid(. ~ cat_biot_gene,drop = T,scales = "free_x", space = "free_x")
a<-a + geom_text(aes(x=cat_biot, y=value, hjust=.5,vjust=-.2,label=value),check_overlap = TRUE,colour="black",position = position_dodge(width = .9),size=4)
a<-a + theme_bw()
a<-a + ylab(expression("Number of txs"))
a<-a + xlab("")
a<-a + scale_fill_manual(values = colli,"")
a<-a + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
a<-a + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=20))
tx_found_n<-a + theme(strip.text.x = element_text(size=15, face="bold"),strip.text.y = element_text(size=20), legend.text=element_text(size=20), legend.title=element_text(size=23), panel.grid = element_blank()) + scale_y_log10(limits = c(1, max(dftx$value, na.rm = TRUE) * 1.1))
pdf(file = "Figures/SupplementaryFigure11_top.pdf",width = 15,height = 5)
tx_found_n
dev.off()

dfwid<-data.frame(cat_biot,cat_biot_gene,trann_all$lentx)
dfwid$cat_biot<-factor(dfwid$cat_biot,levels=c("annotated_coding","annotated_non-coding","novel_tx_R1","novel_tx_R2"))
dfwid$cat_biot_gene<-factor(dfwid$cat_biot_gene,levels=c("protein_coding","lncRNA","pseudogene","non-coding_RNA","novel_gene"))

dfwid$len<-dfwid$trann_all.lentx

d<-ggplot(dfwid,aes(x=cat_biot,y=len,fill=cat_biot))
d<-d + geom_violin(scale="width",draw_quantiles=.5,adjust = 1)
d<-d + scale_y_log10(breaks=c(10,100,1000,10000,100000),labels=c(10,100,1000,10000,100000),limits=c(min(dfwid$len),max(dfwid$len)*1.1))
d<-d + theme_bw()
d<-d + facet_grid(. ~ cat_biot_gene ,drop = T,scales = "free_x",space = "free_x")
d<-d + xlab("")
d<-d + ylab("Tx length (nt)")
d<-d + theme(axis.title.x = element_blank(),axis.text.x = element_blank())
#c<-c + theme(legend.position="none")
d<-d + scale_fill_manual(values = colli,"")
d<-d + theme(axis.title.x = element_blank(),axis.text.x  = element_blank())
d<-d + theme(axis.title.y = element_text(size=20),axis.text.y  = element_text(angle=45, vjust=0.5, size=20))
tx_found_len<-d + theme(strip.text.x = element_text(size=13, face="bold"),strip.text.y = element_text(size=20), legend.text=element_text(size=20), legend.title=element_text(size=23))

pdf(file = "Figures/SupplementaryFigure11_bottom.pdf",width = 15,height = 5)
tx_found_len
dev.off()

