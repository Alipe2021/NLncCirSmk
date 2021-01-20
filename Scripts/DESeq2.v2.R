
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager",update=F)
if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2", update=F) 
if (!requireNamespace("tibble", quietly = TRUE))
  BiocManager::install("tibble", update=F) 
if (!requireNamespace("BiocParallel", quietly = TRUE))
  BiocManager::install("BiocParallel",update=F) 
if (!requireNamespace("optparse", quietly = TRUE))
  install.packages("optparse",update=F)
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr",update=F)
library(BiocParallel)
library(optparse)
library(tibble)
library(DESeq2)
library(dplyr)

rm(list = ls())
#==============================================================
option_list <- list(
  make_option(opt_str = c("-m", "--matrix"), action = "stroe", type = "character", 
              default = FALSE, help = "Gene or transcript count matrix."),
  make_option(opt_str = c("-g", "--group"), action = "stroe", type = "character",
              default = FALSE, help = "group file for each sample"),
  make_option(opt_str = c("-c", "--compare"), action = "stroe", type = "character",
              default = FALSE, help = "Comapare information for each contral and treatment."),
  make_option(opt_str = c("-o", "--outdir"), action = "stroe", type = "character",
              default = "./", help = "Comapare information for each contral and treatment. [default: ./]"),
  make_option(opt_str = c("-f", "--foldchange"), action = "stroe", type = "double",
              default = 2.0, help = "filter by fold change. [default: 2.0] "),
  make_option(opt_str = c("-p", "--fdr"), action = "stroe", type = "double",
              default = 0.05, help = "filter by padjust FDR. [default: 0.05] "),
  make_option(opt_str = c("-t", "--threads"), action = "stroe", type = "integer",
              default = 4, help = "Threads for compute. [defualt: 4] ")
  make_option(opt_str = c("-h", "--help"), type = "logical", action = "store_TRUE",
              default = FALSE, help = "Show this help message and exit ")
)
#==============================================================
opt = parse_args(OptionParser(option_list = option_list, usage = "This Script is a DGE analysis by DESeq2."))

threads = opt$threads
count_matrix_file = opt$matrix
group_file = opt$group
compare_file = opt$compare
outdir = opt$outdir
foldchange_cutoff = opt$foldchange
fdr_cutoff = opt$fdr

if (!dir.exists(outdir)) {
  dir.create(outdir)
}
setwd(outdir); getwd()

# register(MulticoreParam(8)) ## for Linux
register(SnowParam(threads))  ## For Windows

## 
mrna_expr <- read.csv(count_matrix_file, header = T, row.names = 1)
mrna_expr[1:5,1:5];dim(mrna_expr)
## filter median = 0
mrna_expr_median <- apply(X = as.matrix(mrna_expr), MARGIN = 1,FUN = median)
mrna_expr_filter1 <- mrna_expr[mrna_expr_median > 0,]
## filter var < 1
mrna_expr_variance <- apply(X = as.matrix(mrna_expr_filter1), MARGIN = 1,FUN = var)
mrna_expr_filter2 <- mrna_expr_filter1[mrna_expr_variance > 1,]
## filter sum < 1
mrna_expr_final <- mrna_expr_filter2[rowSums(mrna_expr_filter2)>=10, ]
##
group <- read.table(group_file, header = T, sep = "\t")
compare <- read.table(compare_file, header = T, sep = "\t")
head(group);head(compare)

## build dds
dds <- DESeqDataSetFromMatrix(mrna_expr_final, colData = group, design = ~ Group)
dds <- DESeq(dds)

## DEGs
for (i in 1:nrow(compare)) {
  res <- results(
    object = dds, 
    contrast = c("Group",compare[i,1],compare[i,2]),   # FC=expA / expB
    parallel = TRUE
  )  
  
  ## output temp
  out=select(as.data.frame(res),'baseMean','log2FoldChange','pvalue','padj')
  out=rownames_to_column(out,'GeneID')
  out=filter(out,!is.na(padj))
  ## filter by FDR
  FCcut=2
  FDRcut=0.05
  diff <- filter(out, padj<FDRcut, abs(log2FoldChange)>=abs(log2(FCcut)))
  ## output
  file_name <- paste0(compare[i,1],"-vs-", compare[i,2])
  ##
  write.table(out, paste0("DESeq2.", file_name, ".All.xls"), sep = "\t", row.names = F, quote=F)
  write.table(diff, paste0("DESeq2.", file_name, ".DEG.xls"), sep = "\t", row.names = F, quote=F)
  write.table(diff$GeneID, paste0("DESeq2.", file_name, ".sig.txt"), sep = "\t", row.names = F, quote=F, col.names=F)
}

