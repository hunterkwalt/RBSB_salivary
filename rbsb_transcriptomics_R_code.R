library(tximportData)
library(tximport)
library(DESeq2)
library(SummarizedExperiment)
library(knitr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library("pheatmap")
library("RColorBrewer")
library(ggrepel)
library(ggplotify)


setwd("/home/hunter/Desktop/RBSB_transcriptomics/")

############### PREPARE DATA FOR DESEQ2 ####################

#make sample files and tx2gene. The tx2gene is just the first two columns of 
##RSEM's isoform level results
samples= read.table("quant_files_genes/samples.txt", header = T)
tx2gene <- read.delim("new_quant_files_isoform/tx2gene.tsv")

#read in files
files <- file.path("new_quant_files_isoform", paste0(samples$sample, ".isoforms.results"))

#name the files by sample
names(files) <- paste0(samples$sample)
txi.rsem <- tximport(files, type = "rsem", txIn = T, txOut = F, tx2gene = tx2gene)
head(txi.rsem$counts)

#isoform level
txi.rsem.iso <- tximport(files, type = "rsem", txIn = T, txOut = T)
head(txi.rsem.iso$counts)

#add metadata as factors
sampleTable <- data.frame(condition = factor(rep(c("ASD", "ROB","PSG","ASG"), each = 4)))
rownames(sampleTable) <- colnames(txi.rsem$counts)
sampleTable$sex <-factor(rep(c("female", "male"), times = 8))

###################### DESEQ2 ANALYSIS ###############################


#make deseq2 dataset
dds <- DESeqDataSetFromTximport(txi.rsem, sampleTable, ~condition + sex)
levels(dds$sex)

#make the design multifactored w/r/t sex and tissue type
design(dds) <- formula( ~ sex + condition)
dds <- DESeq(dds)


#re-level dataset making the carcass data the control
dds$condition <- relevel(dds$condition, ref = "ROB")
dds <- DESeq(dds)



#this keeps counts >= 10 across all four biological replicates. The ">=4" isn't really necessary. This would change if I had different sized groups
## and 4 would be the size of the smallest group. 
keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep,]


#create results objects. 
res_ASG <- results(dds, name="condition_ASG_vs_ROB", alpha = 0.05)
res_PSG <- results(dds, name="condition_PSG_vs_ROB", alpha = 0.05)
res_ASD <- results(dds, name="condition_ASD_vs_ROB",  alpha = 0.05)

#check differences in sex
res_sex <- results(dds, name="sex_male_vs_female",  alpha = 0.05)
res_sex_sig <- subset(res_sex, padj <= 0.05)
#write.csv(as.data.frame(res_sex_sig), file = "male_vs_female.csv")

#subset significant salivary DEGs
res_ASG_sig <- subset(res_ASG, padj <=0.05)
res_PSG_sig <- subset(res_PSG, padj <=0.05)
res_ASD_sig <- subset(res_ASD, padj <=0.05)

#write.csv(as.data.frame(res_ASG_sig), file = "ASG_vs_ROB.csv")

#write.csv(as.data.frame(res_PSG_sig), file = "PSG_vs_ROB.csv")

#write.csv(as.data.frame(res_AsD_sig), file = "ASD_vs_ROB.csv")

#various data transformations
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
ntd <- normTransform(dds)


#some heatmaps to explore data
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
df <- as.data.frame(colData(dds)[,c("condition", "sex")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=F,
         cluster_cols=T, annotation_col=df)

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=T,
         cluster_cols=F, annotation_col=df)

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=F,
         cluster_cols=T, annotation_col=df)

#correlations/clustering
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition,vsd$sex, sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
distance <- pheatmap(sampleDistMatrix,
                     clustering_distance_rows=sampleDists,
                     clustering_distance_cols=sampleDists,
                     col=colors)
#ggsave("sample_distance_matrix.png", distance)

#PCA
pca_samples <- plotPCA(vsd, intgroup=c("condition", "sex"))
#ggsave("pca_samples.png", pca_samples)
#make supp fig.
supp_fig <- ggarrange(as.grob(distance), pca_samples, labels = c("A", "B"), ncol = 2)

## make individual ma plots
##I made a list of the genes that are significantly differentially expressed in salivary tissues and duplicated (according to OMA).

psg_duplicated <- read.table("/home/hunter/Desktop/oma_07_17/duplicated/significant_annotated/psg_for_R", header = F)
psg_duplicated <- psg_duplicated$V1

asd_duplicated <- read.table("/home/hunter/Desktop/oma_07_17/duplicated/significant_annotated/asd_for_R", header = F)
asd_duplicated <- asd_duplicated$V1

asg_duplicated <- read.table("/home/hunter/Desktop/oma_07_17/duplicated/significant_annotated/asg_for_R", header = F)
asg_duplicated <- asg_duplicated$V1


asd_e_subset <- res_ASD[rownames(res_ASD) %in% asd_duplicated,]
psg_e_subset <- res_PSG[rownames(res_PSG) %in% psg_duplicated,]
asg_e_subset <- res_ASG[rownames(res_ASG) %in% asg_duplicated,]

#get gained (according to OMA) DE genes

psg_gained <- read.table("/home/hunter/Desktop/oma_07_17/gained/psg_significant_gained", header = F)
psg_gained <- psg_gained$V1

asd_gained <- read.table("/home/hunter/Desktop/oma_07_17/gained/asd_significant_gained", header = F)
asd_gained <- asd_gained$V1

asg_gained <- read.table("/home/hunter/Desktop/oma_07_17/gained/asg_significant_gained", header = F)
asg_gained <- asg_gained$V1


asd_g_subset <- res_ASD[rownames(res_ASD) %in% asd_gained,]
psg_g_subset <- res_PSG[rownames(res_PSG) %in% psg_gained,]
asg_g_subset <- res_ASG[rownames(res_ASG) %in% asg_gained,]

asgma <- ggmaplot(res_ASG, main =expression("Accessory Salivary Gland vs. Rest of Body"),
                  fdr = 0.05, fc = 2, size = 2,
                  palette = c("gold2", "dodgerblue3", "darkgray"),
                  genenames = as.vector(rownames(res_ASG)),
                  font.label = c("bold", 11),
                  font.legend = c("bold",12),
                  font.main = c("bold", 15),
                  #label.select = asg_duplicated,
                  legend = "right", top = 0,
                  ggtheme = ggplot2::theme_minimal())

asdma <- ggmaplot(res_ASD, main =expression("Accessory Salivary Duct vs. Rest of Body"),
                  fdr = 0.05, fc = 2, size = 2,
                  palette = c("gold2", "dodgerblue3", "darkgray"),
                  genenames = as.vector(rownames(res_ASD)),
                  font.label = c("bold", 11),
                  font.legend = c("bold",12),
                  font.main = c("bold", 15),
                  #label.select = asd_duplicated, 
                  legend = "right", top = 0,
                  #scale_fill_manual(values = alpha("darkgray", .3)),
                  ggtheme = ggplot2::theme_minimal())

psgma <- ggmaplot(res_PSG, main =expression("Principal Salivary Gland vs. Rest of Body"),
                  fdr = 0.05, fc = 2, size = 2,
                  palette = c("gold2", "dodgerblue3", "darkgray"),
                  genenames = as.vector(rownames(res_PSG)),
                  font.label = c("bold", 11),
                  font.legend = c("bold",12),
                  font.main = c("bold", 15),
                  #label.select =  psg_duplicated,
                  legend = "right", top = 0,
                  ggtheme = ggplot2::theme_minimal())

#make final ma plots with duplicated and gained genes mapped
asdma = asdma + geom_point(data = data.frame(asd_e_subset), aes(y = log2FoldChange, x = log(baseMean, 2)), color = "red", size = 2) +
  geom_point(data = data.frame(asd_g_subset), aes(y = log2FoldChange, x = log(baseMean, 2)), color = "green", size = 2,alpha = 0.7)

psgma = psgma + geom_point(data = data.frame(psg_e_subset), aes(y = log2FoldChange, x = log(baseMean, 2)), color = "red", size = 2) +
  geom_point(data = data.frame(psg_g_subset), aes(y = log2FoldChange, x = log(baseMean, 2)), color = "green", size = 2,alpha = 0.7)

asgma = asgma + geom_point(data = data.frame(asg_e_subset), aes(y = log2FoldChange, x = log(baseMean, 2)), color = "red", size = 2) +
  geom_point(data = data.frame(asg_g_subset), aes(y = log2FoldChange, x = log(baseMean, 2)), color = "green", size = 2,alpha = 0.7)


################################# GO Analysis ########################################

############################## GO Analysis ##############################
library(topGO)

#this is the gene list. I mapped every protein returned from genome annotation against hmmer2go to annotate their GO terms. The file is essentially a tabulated list of
## gene IDs matched to their corresponding GO terms
geneID2GO_full <- readMappings(file = "/home/hunter/Desktop/oma_07_17/GO_RBSB_map.tsv") #GO map from eggnog on the most recent RBSB proteome (07/20/2023)

length(geneID2GO_full)

#subset the expressed genes for GO enrichment within DE.
geneID2GO <- geneID2GO_full[which(rownames(dds) %in% names(geneID2GO_full))] 

head(geneID2GO)
geneNames <- names(geneID2GO)


#function to get gene universe. The score will be the adjusted p-value
getGeneUniverse <- function(deseqRes){
  vgs = as.numeric(deseqRes$padj)
  names(vgs) = rownames(deseqRes)
  vgs
}

#subset the relevant data. In this case, it is all genes that are differentially upregulated. log2FC > 1 = Fold change of 2. 
res_PSG_upregulated <- subset(res_PSG, log2FoldChange >= 1)
res_ASG_upregulated <- subset(res_ASG, log2FoldChange >= 1)
res_ASD_upregulated <- subset(res_ASD, log2FoldChange >= 1)

# write.csv(res_PSG_upregulated, file = "/home/hunter/Desktop/RBSB_transcriptomics/new_analysis/DE/PSG_upregulated.tsv")
# write.csv(res_ASG_upregulated, file = "/home/hunter/Desktop/RBSB_transcriptomics/new_analysis/DE/ASG_upregulated.tsv")
# write.csv(res_ASD_upregulated, file = "/home/hunter/Desktop/RBSB_transcriptomics/new_analysis/DE/ASD_upregulated.tsv")


#downregulated
res_PSG_downregulated <- subset(res_PSG, log2FoldChange <= -1)
res_ASG_downregulated <- subset(res_ASG, log2FoldChange <= -1)
res_ASD_downregulated <- subset(res_ASD, log2FoldChange <= -1)

# write.csv(res_PSG_downregulated, file = "/home/hunter/Desktop/RBSB_transcriptomics/new_analysis/DE/PSG_downregulated.tsv")
# write.csv(res_ASG_downregulated, file = "/home/hunter/Desktop/RBSB_transcriptomics/new_analysis/DE/ASG_downregulated.tsv")
# write.csv(res_ASD_downregulated, file = "/home/hunter/Desktop/RBSB_transcriptomics/new_analysis/DE/ASD_downregulated.tsv")


sel_pval <- function(allScore){ return(allScore < 0.05)}
sel_fc <- function(allScore){ return(allScore > 1)}

#use the above functions to get the gene universe for the GO analysis for each tissue. The gene universe is the set of genes that I am interested in. I also turn all NA's to 1.
asg_gu <- getGeneUniverse(res_ASG_upregulated)
asg_gu[is.na(asg_gu)] <- 1
asd_gu <- getGeneUniverse(res_ASD_upregulated)
asd_gu[is.na(asd_gu)] <- 1
psg_gu <- getGeneUniverse(res_PSG_upregulated)
psg_gu[is.na(psg_gu)] <- 1

#get gene universe downregulated
psg_gu_d <- getGeneUniverse(res_PSG_downregulated)
psg_gu_d[is.na(psg_gu_d)] <- 1

asg_gu_d <- getGeneUniverse(res_ASG_downregulated)
asg_gu_d[is.na(asg_gu_d)] <- 1

asd_gu_d <- getGeneUniverse(res_ASD_downregulated)
asd_gu_d[is.na(asd_gu_d)] <- 1



#this command creates a topGOdata object with whatever ontology you want: Biological Process = "BP", Molecular Function = "MF", Cellular Component = "CC"
## The gene selection in this scenario is based off pvalue, and the gene universe I'm using is the genes that are upregulated in the salivary tissue
#get GO result for Accessory salivary gland
GO_data_asg <- new("topGOdata", ontology = "MF", allGenes = asg_gu, geneSel = sel_pval,
                   annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
allGO = usedGO(object = GO_data_asg)
asg_result_ks.elim <- runTest(GO_data_asg, statistic = "ks", algorithm = "elim")
asg_result_ks <- runTest(GO_data_asg, statistic = "ks", algorithm = "classic")
asg_fisher <- runTest(GO_data_asg, statistic = "fisher", algorithm = "classic")

#allres_asg <- GenTable(GO_data_asg, pval = asg_result, topNodes = length(allGO))
allres_asg <- GenTable(GO_data_asg, classicFisher = asg_fisher, classicKS = asg_result_ks, elimKS = asg_result_ks.elim, 
                       orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 30)


#get go Result for PSG

#this command creates a topGOdata object with whatever ontology you want: Biological Process = "BP", Molecular Function = "MF", Cellular Component = "CC"
## The gene selection in this scenario is based off pvalue, and the gene universe I'm using is the genes that are upregulated in the salivary tissue
GO_data_psg <- new("topGOdata", ontology = "MF", allGenes = psg_gu, geneSel = sel_pval,
                   annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
allGO = usedGO(object = GO_data_psg)
psg_result_ks.elim <- runTest(GO_data_psg, statistic = "ks", algorithm = "classic")
psg_result_ks <- runTest(GO_data_psg, statistic = "ks", algorithm = "classic")
psg_fisher <- runTest(GO_data_psg, statistic = "fisher", algorithm = "classic")

allres_psg_fisher <- GenTable(GO_data_psg, pval = psg_fisher, topNodes = length(allGO))
allres_psg <- GenTable(GO_data_psg, classicFisher = psg_fisher, classicKS = psg_result_ks, elimKS = psg_result_ks.elim, 
                       orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 30)


#get GO result for ASD

#this command creates a topGOdata object with whatever ontology you want: Biological Process = "BP", Molecular Function = "MF", Cellular Component = "CC"
## The gene selection in this scenario is based off pvalue, and the gene universe I'm using is the genes that are upregulated in the salivary tissue
GO_data_asd <- new("topGOdata", ontology = "MF", allGenes = asd_gu, geneSel = sel_pval,
                   annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
allGO = usedGO(object = GO_data_asd)
asd_result_ks.elim <- runTest(GO_data_asd, statistic = "ks", algorithm = "elim")
asd_result_ks <- runTest(GO_data_asd, statistic = "ks", algorithm = "classic")
asd_fisher <- runTest(GO_data_asd, statistic = "fisher", algorithm = "classic")

allres_asd <- GenTable(GO_data_asd, pval = asd_fisher, topNodes = length(allGO))

allres_asd <- GenTable(GO_data_asd, classicFisher = asd_fisher, classicKS = asd_result_ks, elimKS = asd_result_ks.elim, 
                       orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 30)



#### upregulated in ROB vs. PSG ######
GO_data_psg_d <- new("topGOdata", ontology = "MF", allGenes = psg_gu_d, geneSel = sel_pval,
                     annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
allGO = usedGO(object = GO_data_psg_d)
psg_result_d_ks.elim <- runTest(GO_data_psg_d, statistic = "ks", algorithm = "classic")
psg_result_d_ks <- runTest(GO_data_psg_d, statistic = "ks", algorithm = "classic")
psg_fisher_d <- runTest(GO_data_psg_d, statistic = "fisher", algorithm = "classic")

allres_psg_fisher_d <- GenTable(GO_data_psg_d, pval = psg_fisher_d, topNodes = length(allGO))
allres_psg_d <- GenTable(GO_data_psg_d, classicFisher = psg_fisher_d, classicKS = psg_result_d_ks, elimKS = psg_result_d_ks.elim, 
                         orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 100)
kable(allres_psg_d[1:50,])


#### upregulated in ROB vs. ASD ######

GO_data_asg_d <- new("topGOdata", ontology = "MF", allGenes = asg_gu_d, geneSel = sel_pval,
                     annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
allGO = usedGO(object = GO_data_asg_d)
asg_result_d_ks.elim <- runTest(GO_data_asg_d, statistic = "ks", algorithm = "classic")
asg_result_d_ks <- runTest(GO_data_asg_d, statistic = "ks", algorithm = "classic")
asg_fisher_d <- runTest(GO_data_asg_d, statistic = "fisher", algorithm = "classic")

allres_asg_fisher_d <- GenTable(GO_data_asg_d, pval = asg_fisher_d, topNodes = length(allGO))
allres_asg_d <- GenTable(GO_data_asg_d, classicFisher = asg_fisher_d, classicKS = asg_result_d_ks, elimKS = asg_result_d_ks.elim, 
                         orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 100)

kable(allres_asg_d[1:50,])

#### upregulated in ROB vs. ASD ######3

GO_data_asd_d <- new("topGOdata", ontology = "MF", allGenes = asd_gu_d, geneSel = sel_pval,
                     annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
allGO = usedGO(object = GO_data_asd_d)
asd_result_d_ks.elim <- runTest(GO_data_asd_d, statistic = "ks", algorithm = "classic")
asd_result_d_ks <- runTest(GO_data_asd_d, statistic = "ks", algorithm = "classic")
asd_fisher_d <- runTest(GO_data_asd_d, statistic = "fisher", algorithm = "classic")

allres_asd_fisher_d <- GenTable(GO_data_asd_d, pval = asd_fisher_d, topNodes = length(allGO))
allres_asd_d <- GenTable(GO_data_asd_d, classicFisher = asd_fisher_d, classicKS = asd_result_d_ks, elimKS = asd_result_d_ks.elim, 
                         orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 100)
kable(allres_asd_d[1:50,])


##view results
kable(allres_asg[1:30,])
#write.table(allres_asg[1:30,], sep = "\t", quote = F, row.names = F, file = "/home/hunter/Desktop/RBSB_transcriptomics/tables/asg_DE_upregulated_GO.tsv")

kable(allres_psg[1:30,])
#write.table(allres_psg[1:30,], sep = "\t", quote = F, row.names = F, file = "/home/hunter/Desktop/RBSB_transcriptomics/tables/psg_DE_upregulated_GO.tsv")


kable(allres_asd[1:30,])
#write.table(allres_asd[1:30,], quote = F, sep = "\t", row.names = F, file = "/home/hunter/Desktop/RBSB_transcriptomics/tables/asd_DE_upregulated_GO.tsv")

##################### MAKE VENN DIAGRAM #############################
library(VennDiagram)
#subset important data: significant and upregulated.
asg_sig_up <- subset(res_ASG_upregulated, padj <= 0.05)
psg_sig_up <- subset(res_PSG_upregulated, padj <= 0.05)
asd_sig_up <- subset(res_ASD_upregulated, padj <= 0.05)

#get lists of upregulated genes
asg_venn <- row.names(asg_sig_up)
psg_venn <- row.names(psg_sig_up)
asd_venn <- row.names(asd_sig_up)

#make venn diagram. I didn't end up using this one because I wanted a proportional one.
venn.diagram(
  x = list(asg_venn, psg_venn, asd_venn),
  category.names = c("ASG", "PSG", "ASD"),
  filename = "Sailvary_overlap_venn.png",
  output = T,
  col=c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3))
)

#I used the numbers generated from the previous figure to make a proportional venn diagram.
library(eulerr)

fit1 <- euler(c("ASG" = 135, "PSG" = 678, "ASD" = 480, "ASG&PSG" = 87, "PSG&ASD" = 81,
                "ASG&ASD" = 94, "ASG&PSG&ASD" = 50))
par(cex.main = 0.5)
proportional_venn <- plot(fit1,shape = "ellipse", quantities = list(cex = 1), labels = list(font = 2, cex = 1), fill = c("dodgerblue", "gold2", "gray"))


#make figure 2: MA plots with venn diagram.
four_paneled <- plot_grid(psgma, asgma, asdma, proportional_venn, labels = c("A","B","C","D"))

################## DIFFERENTIALLY RETAINED GENE EXPRESSION ########################

#I previously made a list of genes that were lost in other stink bug lineages but retained and expressed in rbsb salivary tissues
## this data came from OMA
ehl <- read.table("/home/hunter/Desktop/oma_07_17/lost/lost_but_expressed/eh", header = F)
ehl <- ehl$V1
nvl <- read.table("/home/hunter/Desktop/oma_07_17/lost/lost_but_expressed/nv", header = F)
nvl <- nvl$V1
hhl <- read.table("/home/hunter/Desktop/oma_07_17/lost/lost_but_expressed/hh", header = F)
hhl <- hhl$V1

ehl_expr <- psg_sig_up[ehl,]
nvl_expr <-psg_sig_up[nvl,]
hhl_expr <-psg_sig_up[hhl,]


#make volcano plots
ggplot(psg_sig_up, aes(x = log2FoldChange, y = -log(padj, 10))) + geom_point() + geom_point(data = ehl_expr, aes(x = log2FoldChange, y = -log(padj, 10), color = "red"))
ggplot(psg_sig_up, aes(x = log2FoldChange, y = -log(padj, 10))) + geom_point() + geom_point(data = nvl_expr, aes(x = log2FoldChange, y = -log(padj, 10), color = "red"))
ggplot(psg_sig_up, aes(x = log2FoldChange, y = -log(padj, 10))) + geom_point() + geom_point(data = hhl_expr, aes(x = log2FoldChange, y = -log(padj, 10), color = "red"))

### GO Enrichment of retained and expressed genes
ehl_gu <- getGeneUniverse(ehl_expr)
ehl_gu[is.na(ehl_gu)] <- 1

nvl_gu <- getGeneUniverse(nvl_expr)
nvl_gu[is.na(nvl_gu)] <- 1

hhl_gu <- getGeneUniverse(hhl_expr)
hhl_gu[is.na(hhl_gu)] <- 1


########### NBSB ################
GO_data_ehl <- new("topGOdata", ontology = "BP", allGenes = ehl_gu, geneSel = sel_pval,
                   annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
allGO = usedGO(object = GO_data_ehl)
ehl_result_ks.elim <- runTest(GO_data_ehl, statistic = "ks", algorithm = "elim")
ehl_result_ks <- runTest(GO_data_ehl, statistic = "ks", algorithm = "classic")
ehl_fisher <- runTest(GO_data_ehl, statistic = "fisher", algorithm = "classic")

#allres_ehl <- GenTable(GO_data_ehl, pval = ehl_result, topNodes = length(allGO))
allres_ehl <- GenTable(GO_data_ehl, classicFisher = ehl_fisher, classicKS = ehl_result_ks, elimKS = ehl_result_ks.elim, 
                       orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 30)

########### BMSB ################
GO_data_hhl <- new("topGOdata", ontology = "BP", allGenes = hhl_gu, geneSel = sel_pval,
                   annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
allGO = usedGO(object = GO_data_hhl)
hhl_result_ks.elim <- runTest(GO_data_hhl, statistic = "ks", algorithm = "elim")
hhl_result_ks <- runTest(GO_data_hhl, statistic = "ks", algorithm = "classic")
hhl_fisher <- runTest(GO_data_hhl, statistic = "fisher", algorithm = "classic")

#allres_hhl <- GenTable(GO_data_hhl, pval = hhl_result, topNodes = length(allGO))
allres_hhl <- GenTable(GO_data_hhl, classicFisher = hhl_fisher, classicKS = hhl_result_ks, elimKS = hhl_result_ks.elim, 
                       orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 30)


########### SGSB ################
GO_data_nvl <- new("topGOdata", ontology = "BP", allGenes = nvl_gu, geneSel = sel_pval,
                   annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 5)
allGO = usedGO(object = GO_data_nvl)
nvl_result_ks.elim <- runTest(GO_data_nvl, statistic = "ks", algorithm = "elim")
nvl_result_ks <- runTest(GO_data_nvl, statistic = "ks", algorithm = "classic")
nvl_fisher <- runTest(GO_data_nvl, statistic = "fisher", algorithm = "classic")

#allres_nvl <- GenTable(GO_data_nvl, pval = nvl_result, topNodes = length(allGO))
allres_nvl <- GenTable(GO_data_nvl, classicFisher = nvl_fisher, classicKS = nvl_result_ks, elimKS = nvl_result_ks.elim, 
                       orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 30)

######## View tables #####

#Brown marmorated results
kable(allres_hhl)
#neotropical brown results
kable(allres_ehl)
#southern green results
kable(allres_nvl)