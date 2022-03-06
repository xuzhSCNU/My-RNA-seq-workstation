# Using featureCounts file to Analysis
# DEseq2 -> MA plot/volcano plot/heatmap -> GO/KEGG

# ##########################################################
# DEseq2 -> MAplot/volcano plot/heatmap -> GO/KEGG/GESA
# ##########################################################

# -------------------------->DEseq2<------------------------
library(DESeq2)
library(ggplot2)

# load table and set datafram and names col and row
raw_df<- read.table(file = "F:\\研究生学习文件\\xu_SCNU\\生物信息学分析\\RNA-seq\\My RNA-seq\\293T-RNASeq-Ctrl_vs_KD.STAR.hg38.featureCounts.tsv",header = T,skip = 1,sep = "\t")
count_df <- raw_df[,c(7:10)]
rownames(count_df) <- as.character(raw_df[,1])
colnames(count_df) <- c("Ctrl_rep1","Ctrl_rep2","METTL3_KD_rep1","METTL3_KD_rep2")

# -------------------------------------------------------->>>>>>>>>>
# make obj 
# -------------------------------------------------------->>>>>>>>>>

# count table 

# filter 
colnames(count_df)
count_df.filter <- count_df[rowSums(count_df) > 20 & apply(count_df,1,function(x){ all(x > 0) }),]

# condition table
sample_df <- data.frame(
  condition = c(rep("ctrl",2), rep("KD",2)),
  cell_line = "293T"
)
rownames(sample_df) <- colnames(count_df.filter)

deseq2.obj <- DESeqDataSetFromMatrix(countData = count_df.filter, colData = sample_df, design = ~condition)
deseq2.obj

# -------------------------------------------------------->>>>>>>>>>
# directly get test result 
# -------------------------------------------------------->>>>>>>>>>

# test 
deseq2.obj <- DESeq(deseq2.obj)

# get result
deseq2.obj.res <- results(deseq2.obj)
deseq2.obj.res.df <- as.data.frame(deseq2.obj.res)

# -------------------------------------------------------->>>>>>>>>>
# step by step get test result 
# -------------------------------------------------------->>>>>>>>>>
deseq2.obj <- DESeqDataSetFromMatrix(countData = count_df, colData = sample_df, design = ~condition)
deseq2.obj

# normalization 
deseq2.obj <- estimateSizeFactors(deseq2.obj)
sizeFactors(deseq2.obj)

# dispersion
deseq2.obj <- estimateDispersions(deseq2.obj)
dispersions(deseq2.obj)

# plot dispersion
plotDispEsts(deseq2.obj, ymin = 1e-4)

# test 
deseq2.obj <- nbinomWaldTest(deseq2.obj)
deseq2.obj.res <- results(deseq2.obj)

# -------------------------->MA plot<------------------------
DESeq2::plotMA(deseq2.obj.res, alpha=0.001)

df <- deseq2.obj.res.df
df$sig[df$padj<= 0.0001] <- c("yes")
df$sig[df$padj > 0.0001] <- c("no")

# -------------------------->volcano plot<------------------------
ggplot(df, aes(x=log2FoldChange, y=-log10(pvalue))) + 
  geom_point(aes(color=sig)) + 
  geom_abline(intercept = 0,slope = 0,linetype="dashed") + 
  theme_bw() + 
  ylim(-1,50) + 
  xlim(-5,5) + 
  scale_color_manual(values=c("gray", "red"))

# -------------------------->GO<------------------------
colnames(deseq2.obj.res.df)
colnames(deseq2.obj.res.df)[2]="log2FC"

# real sign
real_sign = rep("no",nrow(deseq2.obj.res.df))
real_sign
# 1. log2FC > 1 , log2FC < -1
select.log2FC <- abs(deseq2.obj.res.df$log2FC) > 1
table(select.log2FC)
# 2. adj.pVal < 0.05
select.qval <- (deseq2.obj.res.df$padj < 0.001)
table(select.qval)
# 3. count real data
real_sign[select.log2FC & select.qval] <- "yes"
table(real_sign)
# 4 .select sign DEGs
deseq2.filter.DEG <- deseq2.obj.res.df[select.log2FC & select.qval,]

# org.Hs.eg.db是对应物种的参考数据库
# Bioconductor中的orgdb
library(org.Hs.eg.db)
library(clusterProfiler)
library(tidyverse)
# GO 
# 做基因list
DEG.gene_symbol = as.character(rownames(deseq2.filter.DEG))

# 1. 对BP做富集
erich.go.BP = enrichGO(gene = DEG.gene_symbol,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
# 2. 对CC做富集
erich.go.CC = enrichGO(gene = DEG.gene_symbol,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "CC",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
# 3. 对MF做富集
erich.go.MF = enrichGO(gene = DEG.gene_symbol,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "MF",
                       pvalueCutoff = 0.5,
                       qvalueCutoff = 0.5)
# gene symbol : WDR55,TP53....
# gene id : 
# entrez id : 14225
# ensembl id : ENSG00000222564
# ucsc id : ucsc1.2f
# refseq id : NR_xxx,NM_xxx
# UNIPROT : P04637

# barplot(enrich.go.CC)

barplot(erich.go.CC)
barplot(erich.go.BP)
barplot(erich.go.MF)

# dotplot
dotplot(erich.go.CC)
dotplot(erich.go.BP)
dotplot(erich.go.MF)
# adj.pVal一般小于0.01，10^-5,-4

# -------------------------->KEGG<------------------------
# org.Hs.eg.db是对应物种的参考数据库
# Bioconductor中的orgdb
library(org.Hs.eg.db)

# 做基因list
DEG.gene_symbol = as.character(rownames(deseq2.filter.DEG))

# convert id ,KEGG的gene只能使用entrezID，所以要进行ID转换
columns(org.Hs.eg.db)
DEG.entrez_id = mapIds(x = org.Hs.eg.db,
                       keys = DEG.gene_symbol,
                       keytype = "SYMBOL",
                       column = "ENTREZID",
                       multiVals = "first")

# 关于organism，其为KEGG中三字母的缩写
erich.kegg.res <- enrichKEGG(gene = DEG.entrez_id,
                             organism = "hsa",
                             keyType = "kegg")

barplot(erich.kegg.res)
dotplot(erich.kegg.res)

# -------------------------->GSEA<------------------------
library(fgsea)
library(clusterProfiler)

#---------------------------------------------------------------------------->>>>>>>
# 1. load gene set -- value info
#---------------------------------------------------------------------------->>>>>>>
# 下载基因集————MSigDB公用数据库（可以直接百度）
# 使用从MSigDB下载的gmt文件（基因set）直接作为输入，使用read.gmt() 即可

#----------------------------------------------------------->>>>>>>
# 1.1 load gmt file
#----------------------------------------------------------->>>>>>>
gmt.CELL_MIG = clusterProfiler::read.gmt(gmtfile = "F:\\研究生学习文件\\xu_SCNU\\生物信息学分析\\RNA-seq\\My RNA-seq\\CELL_MIGRATION.geneset.gmt")
gmt.CELL_PRO = clusterProfiler::read.gmt(gmtfile = "F:\\研究生学习文件\\xu_SCNU\\生物信息学分析\\RNA-seq\\My RNA-seq\\CELL_PROLIFERATION_GO_0008283.geneset.gmt")

# make gene list（利用fgsea包）
gmt.CELL_MIG.list = list()
# gmt.CELL_MIG.list[[levels(gmt.CELL_MIG$ont)]] = gmt.CELL_MIG$gene
gmt.CELL_MIG.list[["CELL_MIGRATION"]] <- as.character(gmt.CELL_MIG$gene)

#---------------------------------------------------------------------------->>>>>>>
# 2. run GSEA
#---------------------------------------------------------------------------->>>>>>>
# 2.1 get data
DGE_log2FC = as.data.frame(deseq2.filter.DEG$log2FC,row.names = as.character(rownames(deseq2.filter.DEG)))
colnames(DGE_log2FC)[1] <- "log2FC"

# 2.2 get log2FC
gene_log2FC = as.numeric(DGE_log2FC$log2FC)

# 2.3 get rownames(gene_id)
names(gene_log2FC) = as.character(rownames(DGE_log2FC))

# 2.4 filter Inf value
DGE_log2FC.filter = gene_log2FC[! is.infinite(gene_log2FC)]
as.data.frame(DGE_log2FC.filter)

# run 
fgseaRes <- fgsea(pathways = gmt.CELL_MIG.list, 
                  stats    = DGE_log2FC.filter,  #fgsea.value.filter
                  minSize  = 1,
                  maxSize  = 500,        
                  nperm = 1 #
)

plotEnrichment(gmt.CELL_MIG.list[["CELL_MIGRATION"]],
               gene_log2FC   #此处geneList就是value值
)
