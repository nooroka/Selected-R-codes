library(data.table)
library(DESeq2)
library(limma)
library(fgsea)
library(ggplot2)
library(ggrepel)
library("AnnotationDbi")
library(devtools)
library(dplyr)
library("org.Mm.eg.db")
library(phantasus)
library(msigdbr)
library(clusterProfiler)
library(VennDiagram)
library(pheatmap)
library(dplyr)
library(remotes)
library(rUtils)
table          <- read.table("LMNA_AD_raw_counts.txt", header = T)
cond           <- read.csv("conditions.tsv", header=T, sep = '\t', row.names = 1)

# Change names of columns - delete "X" letter
tags           <- sapply(colnames(table), function (x) sub("X", "", basename(x)))
colnames(table) <- tags

# subtract HS treatment from conditions table
cond <- cond[which(cond$Treatment == "control" | cond$Treatment=="A"),]

# the same order of conditions rows and table columns
table[,-1] <- table[,rownames(cond)] #СЌС‚Рѕ РЅРµ СЃСЂР°Р±РѕС‚Р°РµС‚


# Basic statistics - number of reads
col_sums <- colSums(table[,-1])
col_sums
plot(col_sums, xlab="sample")
mean(col_sums)
#СЃСЂРµРґРЅРµРµ РєРѕР»РёС‡РµСЃС‚РІРѕ СЂРёРґРѕРІ РЅР° РѕР±СЂР°Р·РµС†:  10745505
#РѕРґРёРЅ РѕР±СЂР°Р·РµС† РІС‹Р±РёРІР°РµС‚СЃСЏ РЅР° СЂРёСЃСѓРЅРєРµ
#(РїР»СЋСЃ РµС‰С‘ С‚СЂРё СЃСЂРµРґРЅРёРµ С‚РѕР¶Рµ РєР°Рє-С‚Рѕ РЅРµ РѕС‡РµРЅСЊ РІС‹РіР»СЏРґСЏС‚).
col_sums[col_sums < 3000000]
#wt_A_term_2 # Р·Р°РїРёС€РµРј, С‡С‚Рѕ РјР°Р»РµРЅСЊРєРѕРµ РїРѕРєСЂС‹С‚РёРµ
#335 
boxplot(table[,-1])
# Annotation of rows
#table$Gene_id <- sub("\\..*", "", rownames(table))         #for removing .1 ending in ENSEMBL ID if it is
table <- unique(table, by = "Gene_id")      # make unique geneENSMUSG00000102693s
table <- as.data.table(table)
table[, Symbol:= mapIds(org.Mm.eg.db, keys=Gene_id, 
                        keytype="ENSEMBL", column="SYMBOL")] # Create new annotation - gene symbol 
# based on ENSEMBL ID
es <- ExpressionSet(as.matrix(table[,!c("Gene_id", "Symbol")]))
head(exprs(es),20)
fData(es) <- table[,c("Gene_id", "Symbol")]
head(fData(es))
rownames(es) <- fData(es)$Gene_id
pData(es) <- cond
head(pData(es))
write.gct(es,"es_raw.gct")
# filtering (NA exclusion, deduplication)
table <- table[!is.na(table$Symbol)]   # remove rows with NA Symbol
length(unique(table$Symbol))   # Check if we have duplicated Symbol ids
# table <- unique(table, by="Symbol") # simple way to exclude duplicates
#31950
table_dedup <- table[,!c("Gene_id")] %>% group_by(Symbol) %>% 
  summarise_all(median) %>% as.data.table()  # deduplication with summarising reads - remove duplicated SYMBOL
table_dedup[, Gene_id := table[match(table_dedup$Symbol, table$Symbol), "Gene_id"] ]  # Add Gene_id ENSEMBL annotation because it was deleted
table_dedup[, Entrez := mapIds(org.Mm.eg.db, keys=Symbol, keytype="SYMBOL", column="ENTREZID")] # Add Entrez gene annotation
table_dedup$Entrez <- as.character(table_dedup$Entrez) %>% replace(list = table_dedup$Entrez == "NULL", values = NA)  # Replace NULL to NA
table_dedup <- table_dedup[!is.na(table_dedup$Entrez)]
# # Top genes filtering using mean
table_dedup$mean <- rowMeans(table_dedup[,!c("Gene_id", "Symbol", "Entrez")])
table_dedup <- table_dedup[mean > 1,]
table_dedup$mean <- NULL # Delete column "mean"
### 15324 unique (by Gene_id, Entrez, Symbol) genes, mean_expression > 1, not NA Symbols and Entrez
# Save
es <- ExpressionSet(as.matrix(table_dedup[, !c("Gene_id", "Symbol", "Entrez")]))
fData(es) <- table_dedup[, c("Gene_id", "Symbol", "Entrez")]
rownames(es) <- fData(es)$Gene_id
pData(es) <- cond
write.gct(es,"es_filter.gct")
# Normalization and filtering top 12000 genes

es.qnorm.top12K <- es       # create new object of es
exprs(es.qnorm.top12K) <- normalizeBetweenArrays(log2(exprs(es.qnorm.top12K) + 1), method="quantile")  # log2 and quantile normalization
fData(es.qnorm.top12K)$mean <- apply(exprs(es.qnorm.top12K), 1, mean)  # mean calculation
es.qnorm.top12K <- es.qnorm.top12K[order(fData(es.qnorm.top12K)$mean, decreasing = TRUE), ] # ordering of table by mean 
head(exprs(es.qnorm.top12K))
es.qnorm.top12K <- es.qnorm.top12K[1:12000,]  #top 12000 genes substruction
write.gct(es.qnorm.top12K, file="./es.qnorm.top12k.gct")  #save

boxplot(exprs(es.qnorm.top12K), )
## PCA

# Plot PCA on filtered normalized data

# 1 way
pca <- prcomp(t(exprs(es.qnorm.top12K)))
plot(pca$x[, 1:2])

# 3 way
library(affycoretools)
plotPCA(es.qnorm.top12K)
# Clustering: hierarchical and k-means
library(pheatmap)
# z-score
es.qnorm.top12K.z <- es.qnorm.top12K
exprs(es.qnorm.top12K.z) <- t(apply(exprs(es.qnorm.top12K.z),1,scale))
phmap_z <- pheatmap(exprs(es.qnorm.top12K.z), 
                    kmeans_k=16, 
                    clustering_distance_cols = "correlation", 
                    clustering_method = "average",
                    annotation_col = cond,
                    color=colorRampPalette(c("navy", "white", "red"))(50))
head(phmap_z$kmeans$cluster, 20)
outliers <- c("wt_A_term_2")
es.qnorm.top12K.no_out <- es.qnorm.top12K[, !(colnames(es.qnorm.top12K) %in% outliers)]
colnames(es.qnorm.top12K.no_out)
## Lets remove outliers from raw data and update ExpressionSet
# Remove outliers from "table_dedup" (unnormalized table!): two outliers defined
table_dedup <- as.data.frame(table_dedup)
table_dedup <- table_dedup[,!(names(table_dedup) %in% outliers)] %>% as.data.table()
cond <- cond[!(rownames(cond) %in% outliers),]

# New Expression Set
es <- ExpressionSet(as.matrix(table_dedup[, !c("Gene_id", "Symbol", "Entrez")]))
fData(es) <- table_dedup[, c("Gene_id", "Symbol", "Entrez")]
rownames(es) <- fData(es)$Gene_id
pData(es) <- cond
write.gct(es, "es_filt_no_out.gct")

es.design <- model.matrix(~0+Condition, data=pData(es.qnorm.top12K.no_out)) # Choose what condition column will be used for comparison
es.design

fit       <- lmFit(es.qnorm.top12K.no_out, es.design)
fit1      <- contrasts.fit(fit, makeContrasts2(c("Condition", "232_A_term", "wt_A_term"),
                                               levels=es.design))
fit1      <- eBayes(fit1)
de        <- topTable(fit1, adjust.method="BH", number=Inf)
de        <- as.data.table(de, keep.rownames=TRUE)
de        <- de[order(t, decreasing = T)]
de_filt   <- de[adj.P.Val<0.05 & (logFC > 1 | logFC < -1)] # Just filter our DEGs #log,kratny izmeneniu expressii, smotrim
## differential expression analysis with DESeq2
## The input data for the DESeq is the table - filtered and unnormalized - like "es" object for now
## Now we take table with ~15000 genes (it can be cropper to the 12000 genes if you want)


library(DESeq2)
exprs(es) <- apply(exprs(es), 2, as.integer) #make integers type

# Make a normalization with dds
dds <- DESeqDataSetFromMatrix(countData = exprs(es), 
                              colData = pData(es),
                              design=~Condition)     # Condition is a variable in cond table, for which you a going to make DE analysis
dds
dds <- DESeq(dds)
dds

# Now table is normalized with DESeq
# That is more accurate normalization
# Now use it for more accurate PCA plot 
vst <- varianceStabilizingTransformation(dds) # get data from dds
plotPCA(vst,intgroup=c("Condition")) + coord_fixed(ratio = 4)
plotPCA(vst,intgroup=c("Day")) + coord_fixed(ratio = 4)
plotPCA(vst,intgroup=c("Cell_type")) + coord_fixed(ratio = 4)

dists <-  dist(t(assay(vst)))
plot(hclust(dists)) # hclust for columns

## differential expression analysis with DESeq2
## The input data for the DESeq is the table - filtered and unnormalized - like "es" object for now
## Now we take table with ~15000 genes (it can be cropper to the 12000 genes if you want)






# Example of differential expression (DE) calculation for one comparison: wt_control vs 232_control
dir.create("./de/", showWarnings = F)  # Create an output directory
unique(dds$Condition)
# Pay attention to the order of compared conditions
# log2FC = 232_control -  wt_control   
# If  log2FC>0 => gene is upregulated in 232_control; if log2FC<0 => gene is downregulated

de <- results(dds, contrast = c("Condition", "232_A_term", "wt_A_term"), cooksCutoff = F)
head(de)
de <- data.table(ID=rownames(de), as.data.table(de))  
head(de)
tail(de)
head(fData(es))
de <- cbind(de, fData(es)) # add annotation to the DE table
de <- de[order(stat, decreasing = T), ]#sort by t
head(de)
# Look at the specific gene
de[Symbol == "Tnc"]

# Save table
fwrite(de, file="./de/wt_control.vs.232_control.de.tsv", sep="\t")

#Some plots
# Volcanoplot
cols <- densCols(de$log2FoldChange, -log10(de$pvalue))
plot(de$log2FoldChange, -log10(de$padj), col=cols, panel.first=grid(),
     main="Volcano plot", 
     xlab="log2(fold-change)", ylab="-log10(p-value)",
     pch=20, cex=0.6)
abline(v=0, col="grey")
abline(v=c(-1.5,1.5), col="darkgrey")
abline(h=-log10(0.01), col="darkgrey")

gn.selected <- abs(de$log2FoldChange) > 1.5 & de$padj < 0.01
text(de$log2FoldChange[gn.selected],
     -log10(de$padj)[gn.selected],
     lab=de$Symbol[gn.selected ], cex=0.4)
## Pathway analysis

## Lets analyse DEGs (Differentially Expressed Genes) and pathways they belong to
## Take pathways from GO BP (Gene Ontology - Biological Processes) database
# NB! Analyze UP and DOWN genes SEPARATELY
library(clusterProfiler)
de_filt <- de[padj<0.01 & (log2FoldChange > 0.5 | log2FoldChange < -0.5)] # Just filter our DEGs
genes <- de_filt[log2FoldChange   >0.6]$Symbol # Lets take a vector of Entrez names of upregulated genes
go_deg <- enrichGO(genes, 'org.Mm.eg.db', ont="BP", pvalueCutoff=0.05, keyType = "SYMBOL")
head(go_deg)
head(go_deg$Description)
barplot(go_deg, showCategory=20) # first 20 pathways with FDR=0.05 
plotGOgraph(go_deg, firstSigNodes = 10) # plot Graph

# Also you can copy genes and go to MSigDB of another online tool
paste0(genes, collapse = " ")
## GSEA - Gene Set Enrichment Analysis
# FDR = 0.05

# First lets do it using clusterProfiler package :
library(clusterProfiler)
genes <- as.vector(de$stat) # take stat variable from de table - that is t-statistics from DE analysis that is used as ORDER of genes
names(genes) <- as.vector(de$Symbol) #name genes with their symbol
genes <- sort(genes, decreasing = TRUE)
gse <- gseGO(geneList     = genes,
             OrgDb        = org.Mm.eg.db,
             ont          = "BP",
             keyType = "SYMBOL",
             nPerm=10000)
# Dot plot of GSEA
dotplot(gse, showCategory=20, split=".sign",font.size = 7) + facet_grid(.~.sign) #+ ggtitle("dotplot for GSEA")
gseaplot(gse, geneSetID = 5, by = "runningScore", title = gse$Description[5])
gseaplot(gse, geneSetID = 1, by = "runningScore", title = gse$Description[1])

# Now lets do in manually using fgsea package
# Lets download pathways from database
library(msigdbr)
# Gene Ontology biological processes pathways from MSigDB (https://www.gsea-msigdb.org/gsea/msigdb/index.jsp)
m_df <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")#C5 category Gene Ontology
m_df
pathways <- split(m_df$entrez_gene, m_df$gs_name)
# t statistics - just a column from de-table
stats <- de[, setNames(stat, Entrez)]

# Lets get a table of molecular pathways with information about p-value, NES etc
library(fgsea)
fr <- fgsea(pathways, stats, nperm = 100000, nproc=4, minSize=15, maxSize=500)
fr[order(padj)]     # Look at this
fr_res <- fr[order(NES)][padj < 0.01]#filtration
# Let's exclude pathways that are similar to each other, and remain only the most important (padj<0.01)
collapsedPathways <- collapsePathways(fr[order(pval)][padj < 0.01], pathways, stats)
str(collapsedPathways)#there were 103 pathways--> 58 pathways
# The list of remaining pathways is in the collapsedPathways$mainPathways

# Let's order these pathways according to the NES and p-value
mainPathways <- fr[pathway %in% collapsedPathways$mainPathways][
  order(sign(ES)*log(pval)), pathway]

frMain <- fr[match(mainPathways, pathway)] # From the table fr select only mainPathways
# Change ENTREZID to the SYMBOL in the column "leadingEdge"
frMain[, leadingEdge := lapply(leadingEdge, mapIds, 
                               x=org.Mm.eg.db, keytype="ENTREZID", column="SYMBOL")]
# Save table of interesting pathways
pdf("wt_A.vs.232_A.pdf", width=12, height=2 + length(mainPathways) * 0.25)
plotGseaTable(pathways = pathways[mainPathways][1:25], stats = stats, fgseaRes=frMain, gseaParam = 0.5)#first 25 pathways
dev.off()
dir.create("gsea")
fwrite(frMain, file="gsea/wt_control.vs.232_control.filtered.tsv", sep="\t", sep2=c("", " ", ""))
# Draw GSEA plots and save
pdf("wt_A.vs.232_A.pdf", width=12, height=2 + length(mainPathways) * 0.25)
plotGseaTable(pathways = pathways[mainPathways][1:25], stats = stats, fgseaRes=frMain, gseaParam = 0.5)#first 25 pathways
dev.off()


library(msigdbr)
# Gene Ontology biological processes pathways from MSigDB (https://www.gsea-msigdb.org/gsea/msigdb/index.jsp)
hall_df <- msigdbr(species = "Mus musculus", category = "H")
hall_df
pathways <- split(hall_df$entrez_gene, hall_df$gs_name)

# t statistics - just a column from de-table
stats <- de[, setNames(stat, Entrez)]

# Lets get a table of molecular pathways with information about p-value, NES etc
library(fgsea)
fr <- fgsea(pathways, stats, nperm = 100000, nproc=4, minSize=15, maxSize=500)
fr[order(padj)]     # Look at this
fr_res <- fr[order(NES)][padj < 0.01]

# Let's exclude pathways that are similar to each other, and remain only the most important (padj<0.01)
collapsedPathways <- collapsePathways(fr[order(pval)][padj < 0.01], pathways, stats)
str(collapsedPathways)

# The list of remaining pathways is in the collapsedPathways$mainPathways

# Let's order these pathways according to the NES and p-value
mainPathways <- fr[pathway %in% collapsedPathways$mainPathways][
  order(sign(ES)*log(pval)), pathway]

frMain <- fr[match(mainPathways, pathway)] # From the table fr select only mainPathways
# Change ENTREZID to the SYMBOL in the column "leadingEdge"
frMain[, leadingEdge := lapply(leadingEdge, mapIds, 
                               x=org.Mm.eg.db, keytype="ENTREZID", column="SYMBOL")]

# Save table of interesting pathways
#dir.create("gsea")
fwrite(frMain, file="gsea/wt_control.vs.232_control_hall.filtered.tsv", sep="\t", sep2=c("", " ", ""))

# Draw GSEA plots and save
pdf("gsea/wt_control.vs.232_control_hall.pdf", width=12, height=2 + length(mainPathways) * 0.25)
plotGseaTable(pathways = pathways[mainPathways], stats = stats, fgseaRes=frMain, gseaParam = 0.5)
dev.off()
















