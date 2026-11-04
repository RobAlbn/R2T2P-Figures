### 
### Supplementary Table 1 --- Getting numbers of reads (raw reads, cleaned reads, mapped reads, and uniquely mapped reads)
### for the different datasets (Campbell et al. (R1 only as well as R1+R2), Wang et al., and Yao et al.)
###

# Attaching the R package ShortRead
library("ShortRead")

### Campbell et al.

# Reading the file allfiles_final
allfiles_final <- read.table("/home/roberto.albanese/suja_pipeline/alignments_pipeline/allfiles_final", sep="\t")

# Getting numbers of raw reads
fastq_files <- list.files("/group/calviello/Lorenzo/data/Suja/fromgeo", pattern="fastq.gz$", full.names=T)
nums_reads <- sapply(fastq_files, function(x){countFastq(x)$records})
names(nums_reads) <- allfiles_final$V3[match(fastq_files, allfiles_final$V1)]

### --- R1 only 

# Getting number of cleaned reads, uniquely mapped reads, and mapped reads
log_files <- list.files("/scratch/roberto.albanese/R2T2P_paper_runs/SRonly_UniqueORFquant/R2T2P_output_dir/final_alignment", pattern="Log.final.out$", full.names=T, recursive=T)
identical(names(nums_reads),gsub(basename(log_files), pattern="_RNA\\.Log.*|_Ribo\\.Log.*", replacement=""))
log_files <- sapply(log_files, readLines)
cleaned_reads <- as.numeric(sapply(strsplit(log_files[grepl(log_files, pattern="Number of input reads")], split="\t"), "[[", 2))
unique_reads <-  as.numeric(sapply(strsplit(log_files[grepl(log_files, pattern="Uniquely mapped reads number")], split="\t"), "[[", 2))
multim_reads <-  as.numeric(sapply(strsplit(log_files[grepl(log_files, pattern="Number of reads mapped to multiple loci")], split="\t"), "[[", 2))
mapped_reads <- unique_reads + multim_reads

# Creating a dataframe with numbers of raw reads, cleaned reads, mapped reads, and uniquely mapped reads
reads_mat <- cbind(nums_reads, cleaned_reads, mapped_reads, unique_reads)
colnames(reads_mat) <- c("Raw reads", "Cleaned reads", "Mapped reads", "Uniquely mapped reads")
reads_df <- as.data.frame(reads_mat)
reads_df <- rbind(reads_df[grepl(pattern="RNA", rownames(reads_df)),], 
                  reads_df[grepl(pattern="Ribo", rownames(reads_df)),])
reads_df[["Sample name"]] <- rownames(reads_df)
reads_df$Dataset <- "Campbell_R1"
rownames(reads_df) <- NULL
reads_df_Campbell_R1 <- reads_df

### --- R1 + R2

# Getting number of cleaned reads, uniquely mapped reads, and mapped reads
log_files <- list.files("/scratch/roberto.albanese/R2T2P_paper_runs/LRalso_UniqueORFquant/R2T2P_output_dir/final_alignment", pattern="Log.final.out$", full.names=T, recursive=T)
identical(names(nums_reads),gsub(basename(log_files), pattern="_RNA\\.Log.*|_Ribo\\.Log.*", replacement=""))
log_files <- sapply(log_files, readLines)
cleaned_reads <- as.numeric(sapply(strsplit(log_files[grepl(log_files, pattern="Number of input reads")], split="\t"), "[[", 2))
unique_reads <-  as.numeric(sapply(strsplit(log_files[grepl(log_files, pattern="Uniquely mapped reads number")], split="\t"), "[[", 2))
multim_reads <-  as.numeric(sapply(strsplit(log_files[grepl(log_files, pattern="Number of reads mapped to multiple loci")], split="\t"), "[[", 2))
mapped_reads <- unique_reads + multim_reads

# Creating a dataframe with numbers of raw reads, cleaned reads, mapped reads, and uniquely mapped reads
reads_mat <- cbind(nums_reads, cleaned_reads, mapped_reads, unique_reads)
colnames(reads_mat) <- c("Raw reads", "Cleaned reads", "Mapped reads", "Uniquely mapped reads")
reads_df <- as.data.frame(reads_mat)
reads_df <- rbind(reads_df[grepl(pattern="RNA", rownames(reads_df)),], 
                  reads_df[grepl(pattern="Ribo", rownames(reads_df)),])
reads_df[["Sample name"]] <- rownames(reads_df)
reads_df$Dataset <- "Campbell_R1R2"
rownames(reads_df) <- NULL
reads_df_Campbell_R1R2 <- reads_df

### Wang dataset

# Reading the file with samples and groups
multiDE_wang <- read.table("/group/calviello/Lorenzo/fshd_samples/multiDE_wang.tsv", sep="\t", header=F)
multiDE_wang$sample_name <- basename(dirname(dirname(multiDE_wang$V1)))

# Numbers of cleaned reads, mapped reads, and uniquely mapped reads
log_files_paths <- list.files("/scratch/roberto.albanese/R2T2P_paper_runs/Wang_biopsies/R2T2P_output_dir/final_alignment", pattern="Log.final.out$", full.names=T, recursive=T)
log_files_paths <- log_files_paths[basename(dirname(log_files_paths)) %in% multiDE_wang$sample_name]
log_files <- log_files_paths[match(multiDE_wang$sample_name,basename(dirname(log_files_paths)))]
log_files <- sapply(log_files, readLines)
cleaned_reads <- as.numeric(sapply(strsplit(log_files[grepl(log_files, pattern="Number of input reads")], split="\t"), "[[", 2))
unique_reads <-  as.numeric(sapply(strsplit(log_files[grepl(log_files, pattern="Uniquely mapped reads number")], split="\t"), "[[", 2))
multim_reads <-  as.numeric(sapply(strsplit(log_files[grepl(log_files, pattern="Number of reads mapped to multiple loci")], split="\t"), "[[", 2))
mapped_reads <- unique_reads + multim_reads

# Creating a dataframe with numbers of raw reads, cleaned reads, mapped reads, and uniquely mapped reads
reads_mat <- cbind(cleaned_reads, cleaned_reads, mapped_reads, unique_reads)
colnames(reads_mat) <- c("Raw reads", "Cleaned reads", "Mapped reads", "Uniquely mapped reads")
reads_df <- as.data.frame(reads_mat)
reads_df[["Sample name"]] <- multiDE_wang$sample_name
reads_df$Dataset <- paste0("Wang_", multiDE_wang$V3)
reads_df_Wang <- reads_df

### Yao dataset

# Reading the file with samples
multiDE_yao <- read.table("/group/calviello/Lorenzo/fshd_samples/multiDE_yao.tsv", sep="\t", header=F)
multiDE_yao_tubes <- multiDE_yao[grepl(pattern="tubes", multiDE_yao$V3),]
# Adding groups
multiDE_yao_tubes$Group <- ifelse(grepl(pattern="^C", multiDE_yao_tubes$V3), "Control", "FSHD")

# Numbers of cleaned reads, mapped reads, and uniquely mapped reads
log_files_paths <- list.files("/scratch/roberto.albanese/R2T2P_paper_runs/Yao_tubes/R2T2P_output_dir/final_alignment", pattern="Log.final.out$", full.names=T, recursive=T)
log_files_paths <- log_files_paths[basename(dirname(log_files_paths)) %in% multiDE_yao_tubes$V3]
log_files <- log_files_paths[match(multiDE_yao_tubes$V3,basename(dirname(log_files_paths)))]
log_files <- sapply(log_files, readLines)
cleaned_reads <- as.numeric(sapply(strsplit(log_files[grepl(log_files, pattern="Number of input reads")], split="\t"), "[[", 2))
unique_reads <-  as.numeric(sapply(strsplit(log_files[grepl(log_files, pattern="Uniquely mapped reads number")], split="\t"), "[[", 2))
multim_reads <-  as.numeric(sapply(strsplit(log_files[grepl(log_files, pattern="Number of reads mapped to multiple loci")], split="\t"), "[[", 2))
mapped_reads <- unique_reads + multim_reads

# Creating a dataframe with numbers of raw reads, cleaned reads, mapped reads, and uniquely mapped reads
reads_mat <- cbind(cleaned_reads, cleaned_reads, mapped_reads, unique_reads)
colnames(reads_mat) <- c("Raw reads", "Cleaned reads", "Mapped reads", "Uniquely mapped reads")
reads_df <- as.data.frame(reads_mat)
reads_df[["Sample name"]] <- multiDE_yao_tubes$V3
reads_df$Dataset <- paste0("Yao_", multiDE_yao_tubes$Group)
reads_df_Yao <- reads_df

### Concatenating all the dataframes
supp1 <- rbind(reads_df_Campbell_R1,
               reads_df_Campbell_R1R2,
               reads_df_Wang,
               reads_df_Yao)

### Writing Supplementary Table 1 to a TSV file
write.table(supp1, file="/group/calviello/Roberto/R2T2P_scripts_for_Figures_with_container/Supplementary_Tables/Supplementary_Table1.tsv",sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE)

