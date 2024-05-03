setwd("RNA_Seq_QS/")

library(tximport)
library(DESeq2)
library(dplyr)
library(apeglm)
library(tidyverse)
library(pheatmap)

tx2gene <- read.csv("tx2gene.csv") ## required by tximport
lowly_expressed <- c('K562_RS35590','K562_RS15205','K562_RS35800','K562_RS33735','K562_RS32565',
                     'K562_RS27690','K562_RS25690','K562_RS32835','K562_RS24985','K562_RS00435',
                     'K562_RS33730','K562_RS25525','K562_RS01855','K562_RS35275','K562_RS25095',
                     'K562_RS05900','K562_RS35840','K562_RS35270','K562_RS25970','K562_RS22285',
                     'K562_RS08655','K562_RS27430','K562_RS20155','K562_RS35650','K562_RS26940',
                     'K562_RS26775','K562_RS17600','K562_RS33265','K562_RS25825','K562_RS09380',
                     'K562_RS35695','K562_RS35625','K562_RS35395','K562_RS19675','K562_RS09800',
                     'K562_RS06240','K562_RS35445','K562_RS33540','K562_RS25135','K562_RS35525',
                     'K562_RS35475','K562_RS10715','K562_RS34860','K562_RS35575','K562_RS35345',
                     'K562_RS35225','K562_RS28120','K562_RS26040','K562_RS17860','K562_RS04625',
                     'K562_RS35810','K562_RS35795','K562_RS35430','K562_RS35415','K562_RS34190',
                     'K562_RS30545','K562_RS29140','K562_RS18355','K562_RS06640','K562_RS06450',
                     'K562_RS05645','K562_RS00480','K562_RS18405','K562_RS18410')

#Lowly expressed genes where the NumReads row sum < 20. Also added CciI and CciR
tx2gene <- tx2gene[!(tx2gene$GENEID %in% lowly_expressed) & !(tx2gene$TXNAME %in% lowly_expressed), ]

## import quant.sf files as produced by salmon
files <- list.files("Salmon results/")
txi <- tximport(paste("Salmon results/", files, sep=""),
                type="salmon",
                tx2gene = tx2gene,
                countsFromAbundance = "no")

columns <- gsub(pattern = ".sf", replacement = "", x=files)
colnames(txi$counts) <- columns
coldata <- data.frame(sample=colnames(txi$counts), 
                      strain=as.factor(c(rep(c("wt", "wt", "wt", "mut", "mut", "mut"), 3))), 
                      hour=as.factor(c(rep("20", 6),rep("25", 6),rep("30", 6))))

## Strain + hour + combination
dds <- DESeqDataSetFromTximport(txi=txi, colData=coldata, design = ~ strain + hour + strain:hour)
dds <- DESeq(dds, test="LRT", reduced = ~ strain + hour)
resultsNames(dds)
res <- lfcShrink(dds, coef="strain_wt_vs_mut", type="apeglm") # I apply the LFCshrink because I get a lot of genes with large LFC but not significant!
res <- na.omit(res)



vsd <- vst(dds, blind=TRUE)
pcaData <- plotPCA(vsd, intgroup = c("sample", "strain"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = sample, shape = strain)) + 
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + 
  geom_text_repel(aes(label = pcaData$name), box.padding = 0.2) +
  theme_light() +
  ggtitle("PCA on variance stabilized expression levels")


## Plot LFC over counts
#plotMA(res, alpha=0.01, cex=0.6, ylim = c(-3, 3), main='α = 0.01')
#plotMA(res, alpha=0.05, cex=0.6, ylim = c(-3, 3), main='α = 0.05')

## Get top genes
length(res[res$padj<0.01,]$log2FoldChange)
res[res$padj<0.01,]
select <- order(abs(res[res$padj<0.01,]$log2FoldChange), decreasing = TRUE)
top_36 <- res[res$padj<0.01,][select,]
ids <- rownames(top_36)
write.csv(as.data.frame(res), file="lfc_lrt_results.csv")

## The full heatmap
ntd <- normTransform(dds) # which is log2(n + 1)
df <- as.data.frame(colData(dds)[,c("hour","strain")])
pheatmap(assay(ntd), cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


## Short heatmap
mat <- counts(dds, normalized = T)[rownames(top_36),]
mat.z <- t(apply(mat, 1, scale))
colnames(mat.z) <- coldata$sample
pheatmap(mat.z, 
         cluster_rows=FALSE, 
         show_rownames=TRUE,
         cluster_cols=FALSE,
         fontsize_row=6,
         fontsize_col=7,
         annotation_col=df)
#

## Volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                title = '',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 2.5,
                legendPosition = 'top',
                legendLabSize = 7,
                legendIconSize = 3.0)
