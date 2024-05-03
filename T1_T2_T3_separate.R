setwd("RNA_Seq_QS/")

library(tximport)
library(DESeq2)
library(dplyr)
library(apeglm)
library(tidyverse)
library(EnhancedVolcano)
library(RColorBrewer)
library(VennDiagram)

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

"
sample strain hour
1  1549668     wt   20
2  1549669     wt   20
3  1549670     wt   20
4  1549671    mut   20
5  1549672    mut   20
6  1549673    mut   20
7  1549674     wt   25
8  1549675     wt   25
9  1549676     wt   25
10 1549677    mut   25
11 1549678    mut   25
12 1549679    mut   25
13 1549680     wt   30
14 1549681     wt   30
15 1549682     wt   30
16 1549683    mut   30
17 1549684    mut   30
18 1549685    mut   30
"

## Mut vs WT in T1
files <- c("1549671.sf", "1549672.sf", "1549673.sf", "1549668.sf", "1549669.sf", "1549670.sf")
txi <- tximport(paste("Salmon results/", files, sep=""),
                type="salmon",
                tx2gene = tx2gene,
                countsFromAbundance = "no")

columns <- gsub(pattern = ".sf", replacement = "", x=files)
colnames(txi$counts) <- columns
coldata <- data.frame(sample=colnames(txi$counts), 
                      strain=as.factor(c("mut", "mut", "mut", "wt", "wt", "wt")), 
                      hour=as.factor(rep("20", 6)))

dds <- DESeqDataSetFromTximport(txi=txi, colData=coldata, design = ~ strain)
dds <- DESeq(dds, test="Wald")
resultsNames(dds)
res <- results(dds, contrast=c("strain", "mut", "wt")) # I don't want "strain_wt_vs_mut" so I explicitly request mut vs wt
res <- na.omit(res)
nrow(res[res$padj<0.01,]) # 337
nrow(res[res$padj<0.05,]) # 611
nrow(res[abs(res$log2FoldChange)>1,]) # 48
nrow(res[which(res$padj<0.05 & abs(res$log2FoldChange)>1),]) # 45

#T1_upreg <- rownames(res[which(res$padj<0.01 & res$log2FoldChange>0),])
#T1_downreg <- rownames(res[which(res$padj<0.01 & res$log2FoldChange<=0),])
sig_threshold <- 0.05

T1_upreg <- data.frame(res[which(res$padj<sig_threshold & res$log2FoldChange>1),])
T1_downreg <- data.frame(res[which(res$padj<sig_threshold & res$log2FoldChange<(-1)),])

T1_upreg <- rownames(res[which(res$padj<sig_threshold & res$log2FoldChange>1),])
T1_downreg <- rownames(res[which(res$padj<sig_threshold & res$log2FoldChange<(-1)),])

T1_res <-  res
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Wild type vs Mutant at 20h',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0,
                labSize = 2.5,
                legendPosition = 'top',
                legendLabSize = 7,
                legendIconSize = 3.0)

vsd <- vst(dds, blind=TRUE)
pcaData <- plotPCA(vsd, intgroup = c("sample", "strain"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = sample, shape = strain)) + 
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + 
  geom_text_repel(aes(label = pcaData$name), box.padding = 0.2)

## Mut vs WT in T2
files <- c("1549677.sf", "1549678.sf", "1549679.sf", "1549674.sf", "1549675.sf", "1549676.sf")
txi <- tximport(paste("Salmon results/", files, sep=""),
                type="salmon",
                tx2gene = tx2gene,
                countsFromAbundance = "no")

columns <- gsub(pattern = ".sf", replacement = "", x=files)
colnames(txi$counts) <- columns
coldata <- data.frame(sample=colnames(txi$counts), 
                      strain=as.factor(c("mut", "mut", "mut", "wt", "wt", "wt")), 
                      hour=as.factor(rep("25", 6)))

dds <- DESeqDataSetFromTximport(txi=txi, colData=coldata, design = ~ strain)
dds <- DESeq(dds, test="Wald")
resultsNames(dds)
res <- results(dds, contrast=c("strain", "mut", "wt"))
res <- na.omit(res)
nrow(res[res$padj<0.01,]) # 66
nrow(res[res$padj<0.05,]) # 104
nrow(res[abs(res$log2FoldChange)>1,]) # 53
nrow(res[which(res$padj<0.05 & abs(res$log2FoldChange)>1),]) #33

#T2_upreg <- rownames(res[which(res$padj<0.01 & res$log2FoldChange>0),])
#T2_downreg <- rownames(res[which(res$padj<0.01 & res$log2FoldChange<=0),])
sig_threshold <- 0.05

T2_upreg <- data.frame(res[which(res$padj<sig_threshold & res$log2FoldChange>1),])
T2_downreg <- data.frame(res[which(res$padj<sig_threshold & res$log2FoldChange<(-1)),])

T2_upreg <- rownames(res[which(res$padj<sig_threshold & res$log2FoldChange>1),])
T2_downreg <- rownames(res[which(res$padj<sig_threshold & res$log2FoldChange<(-1)),])

T2_res <- res
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Wild type vs Mutant at 25h',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0,
                labSize = 2.5,
                legendPosition = 'top',
                legendLabSize = 7,
                legendIconSize = 3.0)

vsd <- vst(dds, blind=TRUE)
pcaData <- plotPCA(vsd, intgroup = c("sample", "strain"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = sample, shape = strain)) + 
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + 
  geom_text_repel(aes(label = pcaData$name), box.padding = 0.2)

## Mut vs WT in T3
files <- c("1549683.sf", "1549684.sf", "1549685.sf", "1549680.sf", "1549681.sf", "1549682.sf")
txi <- tximport(paste("Salmon results/", files, sep=""),
                type="salmon",
                tx2gene = tx2gene,
                countsFromAbundance = "no")

columns <- gsub(pattern = ".sf", replacement = "", x=files)
colnames(txi$counts) <- columns
coldata <- data.frame(sample=colnames(txi$counts), 
                      strain=as.factor(c("mut", "mut", "mut", "wt", "wt", "wt")), 
                      hour=as.factor(rep("30", 6)))

dds <- DESeqDataSetFromTximport(txi=txi, colData=coldata, design = ~ strain)
dds <- DESeq(dds, test="Wald")
resultsNames(dds)
res <- results(dds, contrast=c("strain", "mut", "wt"))
res <- na.omit(res)
nrow(res[res$padj<0.01,]) # 71
nrow(res[res$padj<0.05,]) # 139
nrow(res[abs(res$log2FoldChange)>1,]) # 50
nrow(res[which(res$padj<0.05 & abs(res$log2FoldChange)>1),]) # 44

#T3_upreg <- rownames(res[which(res$padj<0.01 & res$log2FoldChange>=0),])
#T3_downreg <- rownames(res[which(res$padj<0.01 & res$log2FoldChange<=0),])
sig_threshold <- 0.05

T3_upreg <- data.frame(res[which(res$padj<sig_threshold & res$log2FoldChange>1),])
T3_downreg <- data.frame(res[which(res$padj<sig_threshold & res$log2FoldChange<(-1)),])

T3_upreg <- rownames(res[which(res$padj<sig_threshold & res$log2FoldChange>1),])
T3_downreg <- rownames(res[which(res$padj<sig_threshold & res$log2FoldChange<(-1)),])

T3_res <- res

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Wild type vs Mutant at 30h',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0,
                labSize = 2.5,
                legendPosition = 'top',
                legendLabSize = 7,
                legendIconSize = 3.0)

vsd <- vst(dds, blind=TRUE)
pcaData <- plotPCA(vsd, intgroup = c("sample", "strain"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = sample, shape = strain)) + 
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() + 
  geom_text_repel(aes(label = pcaData$name), box.padding = 0.2)

## Log fold change histograms -2 2
ggplot() +
  geom_density(aes(log2FoldChange, fill="20h"), 
               alpha=0.5,
               color='darkgrey',
               data=as.data.frame(T1_res)) +
  geom_density(aes(log2FoldChange, fill="25h"), 
               alpha=0.2, 
               color='darkgrey',
               data=as.data.frame(T2_res)) +
  geom_density(aes(log2FoldChange, fill="30h"), 
               alpha=0.2,
               color='black',
               data=as.data.frame(T3_res)) +
  theme_minimal() +
  scale_fill_manual(values=c('#FFB575', '#6EB5FF', '#FFFFFF')) +
  xlim(-2,2)

intersection <-  intersect(intersect(c(T1_upreg, T1_downreg), c(T2_upreg, T2_downreg)), c(T3_upreg, T3_downreg))
write.table(intersection, file='DEGs_intersection.txt', quote=F, row.names = F, col.names = F)
venn_colors <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(c(T1_upreg, T1_downreg), 
           c(T2_upreg, T2_downreg), 
           c(T3_upreg, T3_downreg)),
  disable.logging = TRUE,
  category.names = c("20h" , "25h" , "30h"),
  filename = 'Shared.png',
  # Output features
  imagetype="png" ,
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = venn_colors,
  # Numbers
  fontface = "bold",
  fontfamily = "sans",
  cat.fontfamily = "sans",resolution = 800)

## Write tables to files
write.csv(as.data.frame(T1_res), file="DEGs_T1.csv")
write.csv(as.data.frame(T2_res), file="DEGs_T2.csv")
write.csv(as.data.frame(T3_res), file="DEGs_T3.csv")
