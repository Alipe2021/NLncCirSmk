# BiocManager::install("DESeq2")
# BiocManager::install("BiocParallel")
library(BiocParallel)
library(tibble)
library(DESeq2)
library(dplyr)

# register(MulticoreParam(8)) ## for Linux
register(SnowParam(8))  ## For Windows

rm(list = ls())
##
setwd("G:/01-Project/04.LncRNA/Analysis/Expression/")
## 
mrna_expr <- read.csv("transcript_count_matrix.csv", header = T, row.names = 1)
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
group <- read.table("../Pipeline/group.txt", header = T, sep = "\t")
compare <- read.table("../Pipeline/compare.txt", header = T, sep = "\t")
# head(group);head(compare)
## build dds
dds <- DESeqDataSetFromMatrix(mrna_expr_final, colData = group, design = ~ Group)
dds <- DESeq(dds)

## DEGs
for (i in 1:nrow(compare)) {
  # FC=expA / expB
  res <- results(
    object = dds, 
    contrast = c("Group",compare[i,1],compare[i,2]), 
    # alpha = 0.05, lfcThreshold = 2, altHypothesis = "greaterAbs",
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
  file_name <- paste0(abc[i,1],"-vs-", abc[i,2])
  write.table(out, paste0("DESeq2.", file_name, ".All.xls"), sep = "\t", row.names = F, quote=F)
  write.table(diff, paste0("DESeq2.", file_name, ".DEGs.xls"), sep = "\t", row.names = F, quote=F)
  write.table(diff$ID, paste0("DESeq2.", file_name, ".sig.txt"), sep = "\t", row.names = F, quote=F, col.names=F)
}

