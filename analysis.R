#------------------------------------------------------------
# Installations
#------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
BiocManager::install("genefilter")
BiocManager::install("multiClust")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
BiocManager::install("ConsensusClusterPlus")
BiocManager::install("preprocessCore")
BiocManager::install("biotmle")
BiocManager::install("biotmleData")
BiocManager::install("SuperLearner")
BiocManager::install("survtype")
BiocManager::install("MultiAssayExperiment")

#------------------------------------------------------------

#------------------------------------------------------------
# Reinstalls if errors
#------------------------------------------------------------
withr::with_envvar(c(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true"), 
                   remotes::install_github('BioinformaticsFMRP/TCGAbiolinks')
)
withr::with_envvar(c(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true"), 
                   remotes::install_github('grimbough/biomaRt')
)



#------------------------------------------------------------
# Requirements
#------------------------------------------------------------
require(TCGAbiolinks)
require(genefilter)
require(multiClust)
require(preprocessCore)
require(SummarizedExperiment)
require(SuperLearner)
require(biotmle)
require(edgeR)
require(rtracklayer) #(Lawrence, Carey, and Gentleman 2019) provides functions to help parse GTF files.
require(plyr)
require(dplyr)
require(dendextend) #clustering
require(ggplot2)
require(factoextra)
require(ConsensusClusterPlus)
require(MultiAssayExperiment)


#------------------------------------------------------------
# for the next files, we have to download the repository code from
# https://github.com/frangam/wearable-sensor-ml-pipeline
# and include it into our project folder.
# This is for working with machine learning argorithms to make classifications
#------------------------------------------------------------
source("utils/utils.R")
source("utils/predict_fun.R")
source("utils/multiplot.R")
source("utils/sensor_utils.R")
#------------------------------------------------------------



suppressPackageStartupMessages({
  library(survival)
  library(survminer)
})

query <- GDCquery(project="TCGA-SKCM",
                  data.category="Transcriptome Profiling",
                  data.type="Gene Expression Quantification",
                  workflow.type="HTSeq - Counts",
                  experimental.strategy = "RNA-Seq"
                  # file.type = "normalized_results"
)
GDCdownload(query)

# clinical_query <- GDCquery(project="TCGA-SKCM",
#                   data.category="Clinical",
#                   data.type="Clinical Supplement"
# )
# GDCdownload(clinical_query)
# clinical.data <- GDCprepare(clinical_query)
# WriteMatrixToFile(tmpMatrix=clinical.data, tmpFileName="Datathon2020.clinical.txt",
#                   blnRowNames=TRUE, blnColNames=TRUE)

# query <- GDCquery(project="TCGA-SKCM",
#                   data.category="Gene expression",
#                   data.type="Gene expression quantification",
#                   experimental.strategy = "RNA-Seq",
#                   platform = "Illumina HiSeq",
#                   file.type = "normalized_results",
#                   legacy = T
# )

# GDCdownload(query)

# RnaseqSE <- GDCprepare(query)
# Rnaseq_CorOutliers <- TCGAanalyze_Preprocessing(RnaseqSE)
# dataNorm <- TCGAanalyze_Normalization(tabDF = RnaseqSE, geneInfo = geneInfo)
# dataFilt <- TCGAanalyze_Filtering(tabDF = dataNorm, method = "quantile", qnt.cut =  0.25)
# samplesNT <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), typesample = c("NT"))
# samplesTP <- TCGAquery_SampleTypes(barcode = colnames(dataFilt), typesample = c("TP"))

data <- GDCprepare(query)
# clinical_SKCM <- GDCquery_clinic("TCGA-SKCM", "clinical")
# clinical_SKCM$bcr_patient_barcode <- as.character(clinical_SKCM$bcr_patient_barcode)
# clinical_SKCM$bcr_patient_barcode[1:5]

# Rnaseq_CorOutliers <- TCGAanalyze_Preprocessing(data)
class(data)
dim(assay(data))
dim(colData(data))
head(assay(data))
head(colData(data))
# original_data$barcode[1:10]
# original_data$patient[1:10]

datatable(as.data.frame(colData(data)), 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)
datatable(assay(data)[1:100,], 
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = TRUE)
a<-assay(data)[1:100,]

original_data <- as.data.frame(colData(data))
write.csv(colnames(original_data), "exported_data/variables_clinicas.csv", row.names = F)
original_data$age_at_index

# table(original_data$definition)
# table(original_data$sample_type)
# sum(table(original_data$sample_type))
# colnames(original_data)
# rn_od<-rownames(original_data)
# rn_od_uq <- unique(rn_od)
# total_na(original_data)

data_counts <- as.data.frame(assay(data))
# data_counts<- data_counts[, clinical_SKCM$bcr_patient_barcode]
colnames(data_counts)
write.csv(original_data, "exported_data/original_data.csv", row.names = T)
# write.csv(data_dge$counts, "exported_data/data_counts.csv", row.names = T)
write.csv(rowData(data), "exported_data/genes.csv", row.names=T)

#------------------------------------------------------------
# Filter out genes with low counts
#------------------------------------------------------------
original_data_aux <- read.csv("exported_data/original_data.csv", na.strings = "")
data_counts <- read.csv("exported_data/data_counts.csv", row.names = "X")
genes <- read.csv("exported_data/genes.csv", row.names = "X")
# data_dge = DGEList(counts = assay(data), genes = rowData(data))
data_dge = DGEList(counts = data_counts, genes = genes)
# head(data)


# data_counts_matrix<- as.matrix(data_counts)
# colnames(data_counts_matrix) <- colnames(data_counts)
# rownames(data_counts_matrix) <- rownames(data_counts)

#median centered
# median.centered.df = DelayedArray::sweep(data_dge$counts, 1, apply(data_dge$counts, 1, median, na.rm=T))

# Quantile normalization of the dataset
data.norm <- preprocessCore::normalize.quantiles(data_dge$counts, copy=FALSE)
head(data.norm)


# shift data before log scaling to prevent errors from log scaling
# negative numbers
if (min(data.norm)<=0) {
  mindata.norm=abs(min(data.norm)) + .001
  data.norm=data.norm + mindata.norm  
}


# Log2 scaling of the dataset
data.log <- t(apply(data.norm, 1, log2))
min(data.log)
total_na(data.log)



# Write the gene expression and clinical data to text files
WriteMatrixToFile(tmpMatrix=data.log,
                  tmpFileName="Datathon2020.normalized.expression.txt",
                  blnRowNames=TRUE, blnColNames=TRUE)
write.csv(data.log, "exported_data/normalized.expression.csv", row.names = T)

# Obtain gene expression matrix
data.exprs <- read.csv("exported_data/normalized.expression.csv", row.names = "X")
dim(data.exprs)
rownames(data.exprs)
colnames(data.exprs)



# exp_file <- system.file("extdata", "Datathon2020.normalized.expression.txt", package= "multiClust")
# # Load the gene expression matrix 
# data.exprs <- input_file(input=exp_file)

ranked.exprs <- multiClust::probe_ranking(input="Datathon2020.normalized.expression.txt",
                              probe_number=1500, 
                              probe_num_selection="Fixed_Probe_Num",
                              data.exp=data.exprs, 
                              method="SD_Rank")
ranked.exprs[1:4,1:4]
rownames(ranked.exprs) <- rownames(data.exprs)[as.integer(rownames(ranked.exprs))]
write.csv(ranked.exprs,  "exported_data/ranked_normalized.expression.genes.csv", row.names = T)

ranked.exprs <- read.csv("exported_data/ranked_normalized.expression.genes.csv", row.names = "X")
# colnames(ranked.exprs) <- colnames(data.exprs)
#**********************
#nombres de los genes
#**********************


##
## Return a data.frame with gene ids and names
get_genes_df <- function(completedata, mydata){
  gene_ids<-c()
  gene_names<-c()
  rn <- rownames(mydata)
  # print(rn[1:5])
  for(i in 1:length(rn)){
    n <- as.character(completedata[completedata$ensembl_gene_id == rn[i],]$external_gene_name)
    gene_ids[length(gene_ids)+1] <- rn[i]
    gene_names[length(gene_names)+1] <- n
  }
  selected_genes <- as.data.frame(cbind(gene_ids, gene_names))
  return(selected_genes)
}

selected_genes <- get_genes_df(data_dge$genes, ranked.exprs)

head(selected_genes)
# write.csv(selected_genes, "exported_data/TCGA-SKCM_transcriptomic_most_expressed_genes.csv", row.names = F)
selected_genes <- read.csv("exported_data/TCGA-SKCM_transcriptomic_most_expressed_genes.csv")
selected_genes[1:5,]
rownames(ranked.exprs)[1:5]

#**********************
#cambiamos el nombre de las columnas (id genes) por el nombre de los genes
rownames(ranked.exprs) <- selected_genes$gene_names
head(ranked.exprs)
rownames(ranked.exprs)[1480:1500]

write.csv(ranked.exprs, "exported_data/ranked.expression.csv", row.names = T)

ranked.exprs <- read.csv("exported_data/ranked.expression.csv", row.names = "X")
head(ranked.exprs)
# cluster_num <- multiClust::number_clusters(data.exp=data.exprs, Fixed=NULL, gap_statistic=T)

# Call the cluster_analysis function
cluster_num <- 3

hclust_analysis <- multiClust::cluster_analysis(sel.exp=ranked.exprs,
                                    cluster_type="HClust",
                                    distance="euclidean", linkage_type="ward.D2", 
                                    gene_distance="correlation",
                                    num_clusters=cluster_num, data_name="SKCM", 
                                    probe_rank="SD_Rank", probe_num_selection="Fixed_Probe_Num",
                                    cluster_num_selection="Fixed_Clust_Num")

head(hclust_analysis)
table(hclust_analysis)
hclust_analysis <- read.csv(paste("exported_data/SKCM HClust euclidean ward.D2 SD_Rank Fixed_Probe_Num Fixed_Clust_Num Samples.Clusters_", cluster_num,"clusters.csv", sep =""))
hclust_analysis_names <- hclust_analysis$X
hclust_analysis <- hclust_analysis$x
names(hclust_analysis) <- hclust_analysis_names
head(hclust_analysis)
table(hclust_analysis)

#--------------------------------------------------------------------------------------------
# Obtaining the Average Expression for Each Gene/Probe in Each Cluster
#--------------------------------------------------------------------------------------------
# Call the avg_probe_exp function
avg_matrix <- avg_probe_exp(sel.exp=ranked.exprs,
                            samp_cluster=hclust_analysis,
                            data_name="SKCM", cluster_type="HClust", distance="euclidean",
                            linkage_type="ward.D2", probe_rank="SD_Rank",
                            probe_num_selection="Fixed_Probe_Num",
                            cluster_num_selection="Fixed_Clust_Num")
head(avg_matrix)
avg_df <- as.data.frame(avg_matrix)


cluster_1 <- rownames(avg_df[rev(order(avg_df$`Cluster 1`)),])
cluster_2 <- rownames(avg_df[rev(order(avg_df$`Cluster 2`)),])
cluster_3 <- rownames(avg_df[rev(order(avg_df$`Cluster 3`)),])
cluster_4 <- rownames(avg_df[rev(order(avg_df$`Cluster 4`)),])
head(cluster_2)


top10_genes <- data.frame(cluster1=cluster_1, cluster2=cluster_2)
top10_genes <- data.frame(cluster1=cluster_1, cluster2=cluster_2, cluster3=cluster_3)
top10_genes <- data.frame(cluster1=cluster_1, cluster2=cluster_2, cluster3=cluster_3, cluster4=cluster_4)
top10_genes[1:20,]

write.csv(top10_genes, paste("exported_data/genes_expressed_by_cluster_",cluster_num,"clusters.csv", sep=""), row.names = F)
top10_genes <- read.csv(paste("exported_data/genes_expressed_by_cluster_",cluster_num,"clusters.csv", sep=""))
#functions of top genes by each cluster
# cluster1_gen_functions <- queryMany(cluster_1[1:30], scopes='symbol', fields=c('entrezgene', 'go'), species='human')
# cluster2_gen_functions <- queryMany(cluster_2[1:30], scopes='symbol', fields=c('entrezgene', 'go'), species='human')
# cluster3_gen_functions <- queryMany(cluster_3[1:30], scopes='symbol', fields=c('entrezgene', 'go'), species='human')



#PCA
pca_t <- t(ranked.exprs)
head(pca_t)
res.pca <- prcomp(pca_t, scale = TRUE)
names(res.pca)
head(res.pca$x)
res.pca.df <- as.data.frame(res.pca$x[,1:50])
head(res.pca.df)
rownames(res.pca.df)[1:10]
names(hclust_analysis)[1:10]

clusters_values<-c()
for(r in rownames(res.pca.df)){
  clusters_values[length(clusters_values)+1] <- hclust_analysis[r]
}
res.pca.df$cluster <- factor(clusters_values)
var_explained <- res.pca$sdev^2/sum(res.pca$sdev^2)
factoextra::fviz_eig(res.pca)

ggplot(res.pca.df, aes(x=PC1,y=PC2) )+
  geom_point(aes(color=cluster, shape=cluster)) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))+
  theme(legend.position="top")
ggsave(paste("exported_data/pc1_pc2_tcars_",cluster_num,"clusters.png", sep = ""))

head(ranked.exprs)


original_data2 <- original_data
rownames(original_data2) <- gsub("-", ".", rownames(original_data2))
original_data2$barcode <- gsub("-", ".", original_data2$barcode)
df_transpose_clusters <- res.pca.df#as.data.frame(t(ranked.exprs))
rownames(df_transpose_clusters)
clusters_values_2<-c()
sample_type<-c()
tumor_stage<-c()
initial_weight<-c()
year_of_birth<-c()
year_of_death <- c()
gender <- c()
vital_status <- c()


for(r in rownames(df_transpose_clusters)){
  clusters_values_2[length(clusters_values_2)+1] <- hclust_analysis[r]
  sample_type[length(sample_type)+1]<-original_data2[original_data2$barcode == r,]$sample_type
  tumor_stage[length(tumor_stage)+1]<-original_data2[original_data2$barcode == r,]$tumor_stage
  initial_weight[length(initial_weight)+1]<-original_data2[original_data2$barcode == r,]$initial_weight
  year_of_birth[length(year_of_birth)+1]<-original_data2[original_data2$barcode == r,]$year_of_birth
  year_of_death[length(year_of_death)+1]<-original_data2[original_data2$barcode == r,]$year_of_death
  gender[length(gender)+1] <- original_data2[original_data2$barcode == r,]$gender
  vital_status[length(vital_status)+1] <- original_data2[original_data2$barcode == r,]$vital_status
}

# original_data2[original_data2$barcode == "TCGA.EE.A2MC.06A.12R.A18S.07",]$sample_type
# total_na(original_data)

df_transpose_clusters$cluster <- clusters_values_2
df_transpose_clusters$sample_type <- sample_type
df_transpose_clusters$tumor_stage <- tumor_stage
df_transpose_clusters$initial_weight <- initial_weight
df_transpose_clusters$year_of_birth <- year_of_birth
df_transpose_clusters$year_of_death <- year_of_death
df_transpose_clusters$vital_status <- vital_status
df_transpose_clusters$gender <- gender


# pred_df$sample_type <- sample_type
# pred_df$tumor_stage <- tumor_stage
# pred_df$initial_weight <- initial_weight
# pred_df$year_of_birth <- year_of_birth
# pred_df$year_of_death <- year_of_death
# pred_df$vital_status <- vital_status
# pred_df$gender <- gender

ncol(df_transpose_clusters)


#ordenar columnas, poner las dos ultimas al principio
c1 <- ncol(df_transpose_clusters)-7 #cluster
c2 <- ncol(df_transpose_clusters)-6 #sample_type
c3 <- ncol(df_transpose_clusters)-5 #tumor_stage
c4 <- ncol(df_transpose_clusters)-4 #initial_weight
c5 <- ncol(df_transpose_clusters)-3 #year_of_birth
c6 <- ncol(df_transpose_clusters)-2 #year_of_death
c7 <- ncol(df_transpose_clusters)-1 #vital_status
c8 <- ncol(df_transpose_clusters)   #gender
start_cols <- c(c1,c2,c3,c4,c5,c6,c7,c8)
col_rest <- 1:(ncol(df_transpose_clusters)-length(start_cols))
df_transpose_clusters <- df_transpose_clusters[, c(start_cols, col_rest)] 
write.csv(df_transpose_clusters, paste("exported_data/TCGA-SKCM_transcriptomic_most_expressed_clusters_",cluster_num,"clusters.csv", sep = ""), row.names = T)
head(df_transpose_clusters)


df_transpose_clusters <- read.csv(paste("exported_data/TCGA-SKCM_transcriptomic_most_expressed_clusters_",cluster_num,"clusters.csv", sep = ""), row.names = "X")
head(df_transpose_clusters)


df_transpose_clusters$cluster <- as.factor(df_transpose_clusters$cluster)
ggplot(df_transpose_clusters, aes(x=PC1, y = PC2, color = cluster, shape=cluster)) + geom_point() + facet_wrap(~sample_type+vital_status , scales = "free") 
ggsave(paste("exported_data/pc1_pc2_sample_type_vital_status_",cluster_num,"clusters.png", sep = ""))

ggplot(df_transpose_clusters, aes(x=PC1, y = PC2, color = cluster, shape=cluster)) + geom_point() + facet_wrap(~gender+sample_type+vital_status , scales = "free") 
ggsave(paste("exported_data/pc1_pc2_sample_type_vital_status_gender_",cluster_num,"clusters.png", sep = ""))


df_transpose_clusters$KRAS_cna <- pred_df_biomarkers$KRAS_cna
  
ggplot(df_transpose_clusters, aes(x=PC1, y = PC2, color = cluster, shape=cluster)) + geom_point() + facet_wrap(~KRAS_cna+vital_status , scales = "free") 
ggsave(paste("exported_data/pc1_pc2_krascna_vital_status_",cluster_num,"clusters.png", sep = ""))




#prediction
#PCA2
top_n <- 50
top_genes_pred <- unique(c(as.character(top10_genes$cluster1[1:top_n]), as.character(top10_genes$cluster2[1:top_n]), as.character(top10_genes$cluster3[top_n])))
length(top_genes_pred)
ranked.exprs2 <- ranked.exprs[top_genes_pred,]
rownames(ranked.exprs2)
pca_t2 <- t(ranked.exprs2)
# res.pca.pred.df <- as.data.frame(pca_t2)

head(pca_t2)
res.pca2 <- prcomp(pca_t2, scale = TRUE)
names(res.pca2)
head(res.pca2$x)
# factoextra::fviz_eig(res.pca2)
# res.pca2$rotation
# biplot(x = res.pca2, scale = 0, cex = 0.6, col = c("blue4", "brown3"))
pca_num<-2
res.pca.pred.df <- as.data.frame(res.pca2$x[,1:pca_num])
head(res.pca.pred.df)

clusters_values<-c()
for(r in rownames(res.pca.pred.df)){
  clusters_values[length(clusters_values)+1] <- hclust_analysis[r]
}
pred_df2 <- as.data.frame(t(ranked.exprs2))
pred_df2$y <-clusters_values
res.pca.pred.df$y <- factor(clusters_values)
pred_df <- res.pca.pred.df
# pred_df <- res.pca.df
# colnames(pred_df) <- c(colnames(pred_df)[1:top_n], "y")
pred_df$y <- as.integer(pred_df$y)
head(pred_df)
rownames(pred_df)[1:10]
rownames(df_transpose_clusters)[1:10]


#--------------------------------------------------------------------------------------------
# Biomarkers
#--------------------------------------------------------------------------------------------
# clusters_values<-c()
# for(r in colnames(ranked.exprs2)){
#   clusters_values[length(clusters_values)+1] <- hclust_analysis[r]
# }
# biomarkers_df<-as.data.frame(t(ranked.exprs2))
biomarkers_df<-as.data.frame(t(ranked.exprs))
# biomarkers_df$y <- factor(clusters_values)

colnames(df_transpose_clusters)[1:10]



se_col <- df_transpose_clusters[, colnames(df_transpose_clusters) %in% c("sample_type", "tumor_stage", "initial_weight", "vital_status", "gender", "cluster")]
se_col$sample_type <- as.integer(se_col$sample_type)
se_col$tumor_stage <- as.integer(se_col$tumor_stage)
se_col$gender <- as.integer(se_col$gender)
se_col$vital_status <- as.integer(se_col$vital_status)
se_col$initial_weight <- as.integer(se_col$initial_weight)
se_col$cluster <- as.integer(se_col$cluster)
se_col$y <- se_col$cluster 
se_col$cluster <-NULL
se_col <- as.data.frame(cbind(pred_df2[,seq(30)], se_col)) #integramos variables clinicas
total_na(se_col)
se_col[is.na(se_col)] <- 0
total_na(se_col)

# se_col<-replace_missing_values(se_col, check_y_col = F)
se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = DataFrame(t(biomarkers_df))),
                                                 colData = DataFrame(se_col))#DataFrame(pred_df))

require(biotmleData)
# data(illuminaData)
# colData(illuminaData)
# class(illuminaData)
# dim(assay(illuminaData))
# dim(colData(illuminaData))
# head(assay(illuminaData))
colData(se)
dim(assay(se))
idx <- which(names(colData(se)) %in% "cluster")
# total_na(colData(se))
# index_na(colData(se))
# colData(se)[index_na(colData(se)),]
total_na(colData(se))
total_na(assay(se))
biotmle_out2 <- biomarkertmle(se = se,
                             varInt = idx,
                             g_lib = c("SL.mean", "SL.glm"),
                             # Q_lib = c("SL.xgboost", "SL.randomForest", "SL.glmnet", "SL.nnet", "SL.ksvm",
                             #           "SL.bartMachine", "SL.kernelKnn", "SL.rpartPrune", "SL.lm", "SL.mean"),
                             Q_lib = c("SL.xgboost"),
                             cv_folds = 2,
                             parallel = FALSE
)
modtmle_out2 <- modtest_ic(biotmle = biotmle_out2)





# top_biomarkers <- biotmle::toptable(modtmle_out)
# top_biomarkers <- top_biomarkers[top_biomarkers$adj.P.Val < 0.05, ]
# top_biomarkers <- top_biomarkers[order(top_biomarkers$adj.P.Val),]

top_number <- 50
topbiomarkersFDR <- modtmle_out2@topTable %>% subset(adj.P.Val < 
                    1) %>% dplyr::arrange(adj.P.Val) %>% dplyr::slice(seq_len(top_number))
topbiomarkersFDR$ID
rownames(topbiomarkersFDR)
topbiomarkersFDR_names <- rownames(modtmle_out2)[as.integer(topbiomarkersFDR$ID)]



# topbiomarkersFDR_heat <- data.frame()
# topbiomarkersFDR_heat

# require(pheatmap)
# df_se <- as.data.frame(assay(se))
# df_se <- df_se[rownames(df_se) %in% topbiomarkersFDR_names,]
# clusters <- pred_df$y
# df_se <- as.data.frame(t(df_se))
# rownames(t(df_se))[1:10]
# pred_df$y
# rownames(pred_df)[1:10]
# hclust_analysis

# 
# hclust_analysis
# 
# df_col_se<-as.data.frame(colData(se))
# df_col_se <- df_col_se$y
# head()
# topbiomarkersFDR_heat_subset <- df_se[rownames(df_se) %in% ,]
# clusters <- as.data.frame(clusters)
# rownames(clusters) <- colnames(df_se)
# pheatmap(df_se, annotation_row = my_gene_col, annotation_col = my_sample_col)
  # pheatmap(df_se, annotation_col = clusters, cutree_cols = 3)



# top_biomarkers$ID <- rownames(modtmle_out)[as.integer(top_biomarkers$ID)]
# modtmle_out_bionames <- modtmle_out
# rownames(modtmle_out_bionames@tmleOut)
# modtmle_out_bionames@$ID <- top_biomarkers$ID
# rownames(modtmle_out)[as.integer(top_biomarkers$ID[1:10])]
write.csv(topbiomarkersFDR_names, "exported_data/top10_biomarkers.csv", row.names = F)
top_biomarkers <- read.csv("exported_data/top10_biomarkers.csv", stringsAsFactors = F)$x
top_biomarkers_aux <- top_biomarkers


fc_bound <- 1.5
pval_bound <- 0.2
top_biomarkers_aux$volcano = ifelse((top_biomarkers_aux$logFC > fc_bound) & (top_biomarkers_aux$adj.P.Val < pval_bound), 1, 
               ifelse((top_biomarkers_aux$logFC < -fc_bound) & (top_biomarkers_aux$adj.P.Val < pval_bound), -1, 0)) 
top_biomarkers_aux
table(top_biomarkers_aux$volcano)


plot(x = modtmle_out2, type = "pvals_adj")
plot(x = modtmle_out2, type = "pvals_raw")
# design <- as.numeric(designVar == max(designVar))
# heatmap_ic(x = modtmle_out, design = design, FDRcutoff = 0.05, top = 10)
volcano_ic(biotmle = modtmle_out2, fc_bound = 1.5)
ggsave(paste("exported_data/volcano_top10biomarkers_",cluster_num,"clusters.png", sep = ""))


designVar <- as.data.frame(colData(se))[, idx]
designVar <- as.numeric(designVar == 2) # as.numeric(designVar == max(designVar))
# modtmle_out@tmleOut

# build heatmap with top 10 biomarkers
h<-heatmap_ic(x = modtmle_out, left.label = "variable", scale = TRUE,
           clustering.method = "hierarchical", row.dendrogram = TRUE,
           design = designVar, FDRcutoff = 1, top = 10)
topbiomarkersFDR_names <- topbiomarkersFDR_names[c(7,6,5,3,9,8,2,1,10,4)]
ggsave(paste("exported_data/heatmap_top10biomarkers2_",cluster_num,"clusters.png", sep = ""))

# topbiomarkersFDR_names <- rownames(modtmle_out)[as.integer(topbiomarkersFDR$ID)]
# topbiomarkersFDR_names <- topbiomarkersFDR_names[c(29,1,23,19,14,21,28,22,20,12,6,13,26,25,11,3,15,5,7,2,18,9,10,8,4,27,16,24,30,17)]
# topbiomarkersFDR_names <- topbiomarkersFDR_names[c(5,4,3,6,2,13,1,14,9,7,11,12,10,15,8)]

topbiomarkersFDR_names


# ngs_cases <- t(biomarkers_df[biomarkers_df$y==1 | biomarkers_df$y==3,])
# ngs_controls <- t(biomarkers_df[biomarkers_df$y==2,])
# 
# ngs_data <- as.data.frame(cbind(ngs_cases, ngs_controls))
# dim(ngs_data)
# n<-dim(ngs_data)[2]
# exp_var <- c(rep(1, n))
# design <- as.data.frame(exp_var)
# dim(design)
# 
# head(ngs_data[, 1:7])
# 
# 
# 
# ## ----data_proc-------------------------------------------------------------
# se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = DataFrame(ngs_data)),
#                                                  colData = DataFrame(design))
# se
# 
# ## ----biomarkertmle, eval=FALSE---------------------------------------------
# rnaseqTMLEout <-biomarkertmle::biomarkertmle(se = se,
#                                              varInt = 1,
#                                              ngscounts = TRUE,
#                                              parallel = TRUE,
#                                              family = "gaussian",
#                                              g_lib = c("SL.mean", "SL.glm",
#                                                        "SL.randomForest"),
#                                              Q_lib = c("SL.mean", "SL.glm",
#                                                        "SL.randomForest", "SL.nnet")
# )
head(eif(rnaseqTMLEout)$E[, seq_len(6)])


#**********************
#nombres de los genes
#**********************
# gene_ids<-c()
# gene_names<-c()
# rn <- colnames(pred_df)[1:ncol(pred_df)-1]
# for(i in 1:length(rn)){
#   gene_ids[i] <- data_dge$genes[data_dge$genes$external_gene_name == rn[i],]$ensembl_gene_id
#   gene_names[i] <- rn[i]
# }
# selected_genes <- as.data.frame(cbind(gene_ids, gene_names))
# head(selected_genes)
# colnames(pred_df) <- c(selected_genes$gene_ids, "y")
# pred_df$y <- clusters_values_2
colnames(pred_df) <- gsub("-", ".", colnames(pred_df))


pred_df_biomarkers <- se_col

# pred_df2 <- biomarkers_df
# pred_df2$y <- clusters_values
# pred_df2 <- pred_df2[,colnames(pred_df2) %in% c(top_biomarkers, "y")]
colnames(pred_df_biomarkers) <- gsub("-", ".", colnames(pred_df_biomarkers))
pred_df_biomarkers

Sys.setlocale("LC_ALL", "en_US.UTF-8")
table(pred_df_biomarkers$y)
set.seed(33)
parts <- createDataPartition(pred_df_biomarkers$y, p = 0.8, list = F)
validation_data <- pred_df_biomarkers[-parts, ]
train_data <- pred_df_biomarkers[parts, ] #rbind(data.clean_noise[parts, ], data.noise)
table(validation_data$y)
table(train_data$y)



# pred_df$`NKX2-5`



#------------------------------------------------------
# Preprocessing
#------------------------------------------------------
set.seed(33)
if (Sys.getenv("JAVA_HOME")!="") Sys.setenv(JAVA_HOME="")

prep_models <- preprocess_data_several_methods(pred_df_biomarkers, validation_data, fs_method="rf", model="knn", metric = "f1",
                                               prep_methods=list(na=0, enn=c(5)), #rose=0, enn=c(11)),  #tomek=0, smote=c(520, 180)
                                               folds=5, reps=1, verbose=T)

table(prep_models$data$y)
# f1_3_3 <- search_features_selection(prep_models$data, method="rf", model="knn", minfeatTune=1, maxFeatTune = 30, metric="f1", search_method="rfe", folds=5, reps=1, verboseOnlyFeatSel=T, verbose=F)


pred_df_biomarkers$y <- as.integer(pred_df_biomarkers$y)
ggplot(data = pred_df_biomarkers, aes(x = y, y = vital_status, color = y)) +geom_bar() +theme_bw()
ggplot(pred_df_biomarkers, aes(vital_status)) + geom_histogram(aes(y=..density..),binwidth = 3, colour="black", fill="white") + labs(y="Density") + geom_density(alpha=.2, fill="#487C7A", color="#19535F") + facet_wrap(~y) 

# require(car)
# scatterplotMatrix(~pred_df_biomarkers$y+pred_df_biomarkers$vital_status, var.labels=c("Compactness", "Elongatedness"))


# 
# #Como el dataset tiene m??s de 50 samples, usamos el test Kolmogorov-Smirnov con la correcci??n de Lilliefors
# # si tuviera menos de 50 usar??amos test Shapiro-Wilk
# require(nortest)
# by(data = pred_df_biomarkers, INDICES = pred_df_biomarkers$y,FUN = function(x){ nortest::lillie.test(x$initial_weight)})
# #p-value < 2.2e-16 < 0.05 para los 3 clusters
# 
# fligner.test(initial_weight ~ y, pred_df_biomarkers)

plot_normality <- function(data, col="initial_weight"){
  ggplot(data = data, aes_string(x=col)) +
    geom_histogram(aes(y = ..density.., fill = ..count..)) +
    scale_fill_gradient(low = "#DCDCDC", high = "#7C7C7C") +
    stat_function(fun = dnorm, colour = "firebrick",
                  args = list(mean = mean(data[,col]),
                              sd = sd(data[,col]))) +
    ggtitle("Histograma + curva normal te??rica") + facet_wrap(~y) +
    theme_bw()
}




# pheatmap(ranked.exprs, annotation_col = clusters, cutree_cols = 3)
original_data$exposure_id
my_biomarkers_df <- biomarkers_df
my_biomarkers_df$initial_weight <- as.numeric(original_data$initial_weight)
my_biomarkers_df$race <- as.factor(ifelse(original_data$race=="not reported", NA, as.factor(original_data$race)))
my_biomarkers_df$pigment.score <- original_data$paper_PIGMENT.SCORE
my_biomarkers_df$inmune.score <- original_data$paper_LYMPHOCYTE.SCORE
my_biomarkers_df$mutation.subtypes <- as.factor(ifelse(original_data$paper_MUTATIONSUBTYPES=="-", NA, as.character(original_data$paper_MUTATIONSUBTYPES)))
my_biomarkers_df$vital_status <- as.factor(ifelse(original_data$vital_status=="Not Reported", NA, as.character(original_data$vital_status)))
my_biomarkers_df$gender <- original_data$gender
my_biomarkers_df$age <- as.numeric(original_data$age_at_index) # as.factor(ifelse(original_data$age=="-", NA, as.character(original_data$paper_MUTATIONSUBTYPES)))
my_biomarkers_df$uv_rate <- as.numeric(original_data$paper_UV.RATE) # as.factor(ifelse(original_data$age=="-", NA, as.character(original_data$paper_MUTATIONSUBTYPES)))
# my_biomarkers_df$ <- original_data$paper_LYMPHOCYTE.SCORE # as.factor(ifelse(original_data$age=="-", NA, as.character(original_data$paper_MUTATIONSUBTYPES)))
my_biomarkers_df$KRAS_cna <- as.factor(ifelse(original_data$paper_KRAS_cna=="-", NA, ifelse(original_data$paper_KRAS_cna==0, 0, ifelse(original_data$paper_KRAS_cna==-1, -1, ifelse(original_data$paper_KRAS_cna==1, 1, 2)))))
table(original_data$paper_KRAS_cna)
table(my_biomarkers_df$KRAS_cna )
table(my_biomarkers_df$inmune.score,my_biomarkers_df$y)

my_biomarkers_df$sample_type <- original_data$sample_type
my_biomarkers_df$y <- factor(pred_df_biomarkers$y)
write.csv(my_biomarkers_df, "exported_data/all_genes_heatmap.csv", row.names = T)
my_biomarkers_df <- read.csv("exported_data/all_genes_heatmap.csv", row.names = "X")
my_biomarkers_df$y <- factor(my_biomarkers_df$y)
my_biomarkers_df$pigment.score <- factor(my_biomarkers_df$pigment.score)
my_biomarkers_df$inmune.score <- factor(my_biomarkers_df$inmune.score)



require(pheatmap)
os.makedirs("results/heatmaps/", exist_ok=True) 
cf_p <- 200
cf_g <- 1500
top_biomarkers_selected <- ranked.exprs[top_biomarkers,]
top_biomarkers_selected <- ranked.exprs
cf_g <- nrow(top_biomarkers_selected)

plot_heatmap(ranked.exprs[top_biomarkers,],show_rownames = T, my_biomarkers_df, width = 10, height = 8, fontsize_row = 12, filename = "results/heatmaps/heatmap_top10.png")
plot_heatmap(ranked.exprs, my_biomarkers_df, width = 7, height = 7, filename = "results/heatmaps/heatmap_all_reduced.png")

# plot_heatmap(ranked.exprs, my_biomarkers_df, width = 12, height = 11, filename = "exported_data/heatmaps/heatmap_all.png")
plot_heatmap(ranked.exprs, my_biomarkers_df, width = 15, height = 11, filename = "results/heatmaps/heatmap_all.png")

# plot_heatmap(ranked.exprs[top_biomarkers,], my_biomarkers_df, top_biomarkers_selected, width = 15, height = 11, filename = "exported_data/heatmaps/heatmap_top10.png")


#heatmap top biomarkers
top_n <- 30
top_genes_pred <- unique(c(as.character(top10_genes$cluster1[1:top_n]), as.character(top10_genes$cluster2[1:top_n]), as.character(top10_genes$cluster3[top_n])))
length(top_genes_pred)
ranked.exprs2 <- ranked.exprs[top_genes_pred,]
plot_heatmap(ranked.exprs2, my_biomarkers_df, show_rownames = T, width = 12, height = 9.5, fontsize_row = 12, scale="none", filename = "results/heatmaps/heatmap_top10.png")


#------------------------------------------------------
# Survival Analysis
#------------------------------------------------------

os.makedirs("results/survival/", exist_ok=True) 


# 
# cutoff_genes <- nrow(ranked.exprs[top_biomarkers,])
# cutoff_patients <-  ncol(ranked.exprs[top_biomarkers,])
# clusters_hp <-data.frame(my_biomarkers_df$y[1:cutoff_patients], my_biomarkers_df$sample_type[1:cutoff_patients]
#                          , my_biomarkers_df$inmune.score[1:cutoff_patients], my_biomarkers_df$pigment.score[1:cutoff_patients]
#                          , my_biomarkers_df$mutation.subtypes[1:cutoff_patients]
#                          , my_biomarkers_df$initial_weight[1:cutoff_patients], my_biomarkers_df$uv_rate[1:cutoff_patients]
#                          , my_biomarkers_df$KRAS_cna[1:cutoff_patients])
# head(clusters_hp)
# rownames(clusters_hp) <- rownames(my_biomarkers_df[,1:cutoff_patients])
# colnames(clusters_hp) <- c("Cluster", "Tissue Origin", "Immune Score", "Pigment Score", "Mutation Subtypes", "Initial Weight", "UV Rate", "KRAS CNA")
# head(ranked.exprs[top_biomarkers,][1:cutoff_genes,1:cutoff_patients])
# pheatmap(ranked.exprs[top_biomarkers,][1:cutoff_genes,1:cutoff_patients], color = colorRampPalette(c("green", "black", "red"))(7)
#          , cluster_rows = T,  clustering_distance_rows = "correlation",fontsize_row=12
#          , clustering_distance_cols = "euclidean", scale="column", width=15, height = 10
#          , show_rownames = T, show_colnames = F, annotation_col = clusters_hp, filename = "exported_data/heatmaps/heatmap_top10.png")
# 
# pheatmap(ranked.exprs[top_biomarkers,][1:cutoff_genes,1:cutoff_patients], color = colorRampPalette(c("green", "black", "red"))(7)
#          , cluster_rows = T,  clustering_distance_rows = "correlation",fontsize_row=12
#          , clustering_distance_cols = "euclidean", scale="none", width=15, height = 10
#          , show_rownames = T, show_colnames = F, filename = "exported_data/heatmaps/heatmap_top10.png")
# 
# 


dim(ranked.exprs)
dim(my_biomarkers_df)

my_biomarkers_df$sample_type[1:5]

plot_heatmap <- function(exp.df, my_biomarkers_df,top_biomarkers_selected, cutoff_patients=ncol(exp.df), cutoff_genes=nrow(exp.df), width=18.75, height = 15.625, cellwidth = NA, cellheight = NA, fontsize_row=4.6, filename = "exported_data/heatmap.png",
                         show_rownames = FALSE, show_colnames = FALSE, scale="column"){
  clusters_hp <-data.frame(my_biomarkers_df$y[1:cutoff_patients], my_biomarkers_df$sample_type[1:cutoff_patients]
                           , my_biomarkers_df$inmune.score[1:cutoff_patients], my_biomarkers_df$pigment.score[1:cutoff_patients]
                           , my_biomarkers_df$mutation.subtypes[1:cutoff_patients]
                           , my_biomarkers_df$initial_weight[1:cutoff_patients], my_biomarkers_df$uv_rate[1:cutoff_patients]
                           , my_biomarkers_df$KRAS_cna[1:cutoff_patients])
  rownames(clusters_hp) <- rownames(my_biomarkers_df[,1:cutoff_patients])
  colnames(clusters_hp) <- c("Cluster", "Tissue Origin", "Immune Score", "Pigment Score", "Mutation Subtypes", "Initial Weight", "UV Rate", "KRAS CNA")
  pheatmap(exp.df[1:cutoff_genes,1:cutoff_patients], color = colorRampPalette(c("green", "black", "red"))(7)
           , cluster_rows = T,  clustering_distance_rows = "correlation",fontsize_row=fontsize_row
           , clustering_distance_cols = "euclidean", scale=scale, width=width, height = height, cellwidth = cellwidth, cellheight = cellheight
           , show_rownames = show_rownames, show_colnames = show_colnames, annotation_col = clusters_hp, filename = filename)
  
}



colnames(ranked.exprs)[1:5]

#---------------------------
# Tests
#---------------------------
pred_df_biomarkers

colnames(pred_df_biomarkers)

original_data$sample_type
original_data$paper_PIGMENT.SCORE
original_data$paper_LYMPHOCYTE.SCORE
original_data$paper_MUTATIONSUBTYPES
rownames(original_data)[1:10]
rownames(pred_df_biomarkers)[1:10]


pred_df_biomarkers$pigment.score <- ifelse(is.na(original_data$paper_PIGMENT.SCORE), -1,
                                           as.integer(original_data$paper_PIGMENT.SCORE)-2)
pred_df_biomarkers$inmune.score <- ifelse(is.na(original_data$paper_LYMPHOCYTE.SCORE), -1, 
                                          ifelse(original_data$paper_LYMPHOCYTE.SCORE==0, 0,
                                                 as.integer(original_data$paper_LYMPHOCYTE.SCORE)-1))
pred_df_biomarkers$mutation.subtypes <- ifelse(is.na(original_data$paper_MUTATIONSUBTYPES), -1,
                                               as.integer(original_data$paper_MUTATIONSUBTYPES))
pred_df_biomarkers$age <- ifelse(is.na(original_data$age_at_index), -1,
                                               as.integer(original_data$age_at_index))
pred_df_biomarkers$race <-ifelse(original_data$race=="not reported", -1, as.factor(original_data$race))
pred_df_biomarkers$uv_rate <-ifelse(is.na(original_data$paper_UV.RATE), -1, as.factor(original_data$paper_UV.RATE))

pred_df_biomarkers$KRAS_cna<- as.integer(ifelse(is.na(original_data$paper_KRAS_cna), -2, ifelse(original_data$paper_KRAS_cna==0, 0, ifelse(original_data$paper_KRAS_cna==-1, 3, ifelse(original_data$paper_KRAS_cna==1, 1, 2)))))

pred_df_biomarkers[is.na(pred_df_biomarkers)] <- 0
table(original_data$paper_MUTATIONSUBTYPES)


#si no vienen de distribución normal mejor usar el test de flinger o el de levene
fligner.test(initial_weight ~ y, pred_df_biomarkers)
fligner.test(tumor_stage ~ y, pred_df_biomarkers) #No
fligner.test(sample_type ~ y, pred_df_biomarkers) #No
fligner.test(vital_status ~ y, pred_df_biomarkers)
fligner.test(gender ~ y, pred_df_biomarkers)
fligner.test(pigment.score ~ y, pred_df_biomarkers[pred_df_biomarkers$pigment.score!=-1,])
fligner.test(inmune.score ~ y, pred_df_biomarkers[pred_df_biomarkers$mutation.subtypes!=-1,])
fligner.test(mutation.subtypes ~ y, pred_df_biomarkers[pred_df_biomarkers$mutation.subtypes!=-1,])
fligner.test(age ~ y, pred_df_biomarkers[pred_df_biomarkers$age!=-1,])
fligner.test(race ~ y, pred_df_biomarkers[pred_df_biomarkers$race!=-1,])
fligner.test(uv_rate ~ y, pred_df_biomarkers[pred_df_biomarkers$uv_rate!=-1,])
fligner.test(KRAS_cna ~ y, pred_df_biomarkers[pred_df_biomarkers$KRAS_cna!=-2,])

pred_df_biomarkers$uv_rate


plot_normality(pred_df_biomarkers, "initial_weight")
plot_normality(pred_df_biomarkers, "vital_status")
plot_normality(pred_df_biomarkers, "sample_type")
plot_normality(pred_df_biomarkers, "tumor_stage")
plot_normality(pred_df_biomarkers, "gender")
plot_normality(pred_df_biomarkers, "inmune.score")
plot_normality(pred_df_biomarkers, "pigment.score")
plot_normality(pred_df_biomarkers, "mutation.subtypes")





anova <- aov(pred_df_biomarkers$vital_status ~pred_df_biomarkers$y)
summary(anova)
par(mfrow = c(2,2))
plot(anova)
par(mfrow = c(1,1))
par(mfrow = c(1,3))
qqnorm(pred_df_biomarkers[pred_df_biomarkers$y == 1, "vital_status"], main = "Cluster 1")
qqline(pred_df_biomarkers[pred_df_biomarkers$y == 1,"vital_status"])
qqnorm(pred_df_biomarkers[pred_df_biomarkers$y == 2, "vital_status"], main = "Cluster 2")
qqline(pred_df_biomarkers[pred_df_biomarkers$y == 2,"vital_status"])
qqnorm(pred_df_biomarkers[pred_df_biomarkers$y == 3, "vital_status"], main = "Cluster 3")
qqline(pred_df_biomarkers[pred_df_biomarkers$y == 3,"vital_status"])
par(mfrow = c(1,1))

#Test de normalidad
by(data = pred_df_biomarkers, INDICES = pred_df_biomarkers$y,FUN = function(x){ nortest::lillie.test(x$initial_weight)})
by(data = pred_df_biomarkers, INDICES = pred_df_biomarkers$y,FUN = function(x){ nortest::lillie.test(x$vital_status)})
by(data = pred_df_biomarkers, INDICES = pred_df_biomarkers$y,FUN = function(x){ nortest::lillie.test(x$sample_type)})
by(data = pred_df_biomarkers, INDICES = pred_df_biomarkers$y,FUN = function(x){ nortest::lillie.test(x$tumor_stage)})
by(data = pred_df_biomarkers, INDICES = pred_df_biomarkers$y,FUN = function(x){ nortest::lillie.test(x$gender)})
by(data = pred_df_biomarkers[pred_df_biomarkers$pigment.score!=-1,], INDICES = pred_df_biomarkers[pred_df_biomarkers$pigment.score!=-1,]$y,FUN = function(x){ nortest::lillie.test(x$pigment.score)})
by(data = pred_df_biomarkers[pred_df_biomarkers$inmune.score!=-1,], INDICES = pred_df_biomarkers[pred_df_biomarkers$inmune.score!=-1,]$y,FUN = function(x){ nortest::lillie.test(x$inmune.score)})


#Test de homogeneidad
plot(pigment.score ~ y, data = pred_df_biomarkers)
plot(inmune.score ~ y, data = pred_df_biomarkers)
plot(mutation.subtypes ~ y, data = pred_df_biomarkers)
plot(sample_type ~ y, data = pred_df_biomarkers)
plot(vital_status ~ y, data = pred_df_biomarkers)
plot(vital_status ~ y, data = pred_df_biomarkers)

#bartlett requiere que sean distribución normal
# bartlett.test(initial_weight ~ y, data=pred_df_biomarkers)
# bartlett.test(tumor_stage ~ y, data=pred_df_biomarkers)
# bartlett.test(sample_type ~ y, data=pred_df_biomarkers)
# bartlett.test(vital_status ~ y, data=pred_df_biomarkers)
# bartlett.test(gender ~ y, data=pred_df_biomarkers)
# bartlett.test(pigment.score ~ y, data=pred_df_biomarkers[pred_df_biomarkers$pigment.score!=-1,])
# bartlett.test(inmune.score ~ y, data=pred_df_biomarkers[pred_df_biomarkers$inmune.score!=-1,])
# bartlett.test(mutation.subtypes ~ y, data=pred_df_biomarkers[pred_df_biomarkers$mutation.subtypes!=-1,])




leveneTest(pigment.score ~ y, data = pred_df_biomarkers)
leveneTest(inmune.score ~ y, data = pred_df_biomarkers)
leveneTest(mutation.subtypes ~ y, data = pred_df_biomarkers)
leveneTest(sample_type ~ y, data = pred_df_biomarkers)
leveneTest(vital_status ~ y, data = pred_df_biomarkers)
leveneTest(gender ~ y, data = pred_df_biomarkers)



ggplot(pred_df_biomarkers, aes(y,vital_status)) + 
  geom_point() +
  theme_bw()+ geom_smooth(method=lm)

# fisher.test(initial_weight ~ y, data=pred_df_biomarkers)


require(moments)
# moments::agostino.test(pred_df_biomarkers$initial_weight ~ pred_df_biomarkers$y)


pred_df_biomarkers$tumor_stage
# require(survtype)
# SKCM.survtype <- survtype(clinical_SKCM, time = "days_to_last_follow_up",
#                                 status = "vital_status", assay(data),
#                                 num.genes = 100, scale = "row",
#                                 gene.sel = FALSE, clustering_method = "ward.D2",
#                                 show_colnames = FALSE)
coldat <- as.data.frame(colData(data))
coldat$y <- survival::Surv(data$days_to_death, as.integer(data$vital_status))
colData(data) <- DataFrame(coldat)
data_aux <- data[complete.cases(coldat$y), ]
coldat <- as(colData(data_aux), "data.frame")
fit <- survfit(y ~ tumor_stage, data = coldat)
ggsurvplot(fit, data = coldat, risk.table = TRUE)

original_data$state
data("ovarian")


require(survtype)
surv.df <- my_biomarkers_df
surv.df$days_to_death <- as.integer(original_data$days_to_death)
surv.df$days_to_last_follow_up <-  as.integer(original_data$days_to_last_follow_up)
surv.df$days_to_collection <- as.integer(ifelse(is.na(original_data$days_to_collection), 0, original_data$days_to_collection))
surv.df$death_event <- ifelse(original_data$vital_status=="Not Reported", NA, 
                               ifelse(original_data$vital_status=="Alive", 0, 1))

ydata$time
surv.df$patientID <- rownames(surv.df)
original_data$days_to_death
original_data$tum

# time = max(surv.df$Days.to.date.of.Death, 
#            surv.df$Days.to.Last.Contact, na.rm = TRUE)


ydata.raw <- surv.df %>% as.data.frame %>% 
  # Keep only data relative to survival or samples
  select(patientID, death_event, days_to_death, days_to_last_follow_up, days_to_collection, y
     ) %>% 
  # # Convert days to integer
  # mutate(Days.to.date.of.Death = as.integer(Days.to.date.of.Death)) %>%
  # mutate(Days.to.Last.Contact  = as.integer(Days.to.Date.of.Last.Contact)) %>%
  # Find max time between all days (ignoring missings)
  rowwise %>%
  mutate(time = ifelse(is.na(days_to_death), days_to_last_follow_up, days_to_death)) %>%
  # mutate(time = max(days_to_last_follow_up, days_to_death, na.rm = TRUE)) %>%

  
  # mutate(time = max(days_to_death, na.rm = TRUE)) %>%
  
  # Keep only survival variables and codes
  select(patientID, status = death_event, time, y) %>% 
  # Discard individuals with survival time less or equal to 0
  filter(!is.na(time) & time > 0 & !is.na(status)) %>% as.data.frame

# Set index as the patientID
rownames(ydata.raw) <- ydata.raw$patientID

# Get matches between survival and assay data
xdata.raw <- surv.df[rownames(surv.df) %in% rownames(ydata.raw),]
# xdata.raw <- xdata.raw[,1:ncol(xdata.raw)] %>%
# { (apply(., 2, sd) != 0) } %>%
# { xdata.raw[, .] } %>%
#   scale

# Order ydata the same as assay
ydata.raw <- ydata.raw[rownames(xdata.raw), ]
ydata <- ydata.raw %>% select(time, status, y)

# ydata$time <- ydata$time/365.25
# ydata <- ydata[ydata$time<=10,]
# 
# skcm.survtype <- survtype::Surv.survtype(ydata, time="time", status = "status")
# plot.survtype(skcm.survtype, pval = TRUE)
# ggsave("results/survival/survival.png")
# survtype::Sur
ydata$time
ydata$x<-ydata$y
ydata$y <- NULL
levels(ydata$x) <- c('Cluster 1: Immune','Cluster 2: Keratin','Cluster 3: MAGE')
skcm.surv.fit <- survival::survfit(survival::Surv(time,status) ~ x, data = ydata)
surv_pvalue(skcm.surv.fit)
ggsurvplot(
  skcm.surv.fit,
  conf.int = FALSE,
  # surv.median.line = c('hv'), 
  data = ydata, 
  pval = TRUE,
  pval.method = TRUE,
  risk.table = FALSE, xlab="Year", xscale="d_y", #transformar días a años
  pval.method.coord=c(0,0), pval.coord=c(365.25*5,0),
  font.legend=12
  , break.x.by=365.25*2
  # , xlim=365.25*10
  )

table(my_biomarkers_df$sample_type, my_biomarkers_df$y)

ggsave("results/survival/survival.png")
(original_data$days_to_collection - original_data$days_to_birth) / 365.25

plot.survtype(skcm.surv, pval = TRUE, pval.method.coord=c(5,0), pval.coord=c(5,-1), font.legend=12)
ggsave("results/survival/survival.png")
autoplot(skcm.surv, conf.int = F) +
  ggtitle('KM survival curve for BLCA cohort')

which(clinical_SKCM$days_to_last_follow_up == -Inf)
rownames(clinical_SKCM)
colnames(clinical_SKCM)
original_data2$barcode
colnames(original_data2)
clinical_SKCM$bcr_patient_barcode[1:5]
original_data2$patient[1:5]

#------------------------------------------------------
# Classification
#------------------------------------------------------
# features <- search_features_selection(prep_models$data, method="rf", model="rf", minfeatTune=5, maxFeatTune = 15, metric="f1", search_method="rfe", folds=5, reps=1, verboseOnlyFeatSel=T, verbose=F)
# features <- search_features_selection(data, method="rf", model="svm", minfeatTune=5, maxFeatTune = 15, metric="f1", search_method="rfe", folds=5, reps=1, verboseOnlyFeatSel=T, verbose=F)

dim(pred_df)

f1_3_1 <- search_features_selection(pred_df_biomarkers, method="rf", model="rf", minfeatTune=1, metric="f1", search_method="rfe", folds=3, reps=1, verboseOnlyFeatSel=T, verbose=F)
f1_3_2 <- search_features_selection(pred_df_biomarkers, method="rf", model="svm", minfeatTune=1, metric="f1", search_method="rfe", folds=3, reps=1, verboseOnlyFeatSel=T, verbose=F)
f1_3_3 <- search_features_selection(pred_df_biomarkers, method="rf", model="knn", minfeatTune=1, metric="f1", search_method="rfe", folds=3, reps=1, verboseOnlyFeatSel=T, verbose=F)
f1_3_4 <- search_features_selection(pred_df_biomarkers, method="rf", model="ann", minfeatTune=1, metric="f1", search_method="rfe", folds=3, reps=1, verboseOnlyFeatSel=T, verbose=F)
f1_3_5 <- search_features_selection(pred_df_biomarkers, method="rf", model="gbm", minfeatTune=1, metric="f1", search_method="rfe", folds=3, reps=1, verboseOnlyFeatSel=T, verbose=F)

plot_features_vs_acc_list_features(list(f1_3_1$all, f1_3_2$all, f1_3_3$all, f1_3_5$all), metric="f1", tech=c("RF", "SVM", "k-NN", "GBM"), lineBest=F, interv = 8)
plot_features_vs_acc_list_features(list(f1_3_1$all, f1_3_2$all, f1_3_3$all, f1_3_5$all), metric="acc", tech=c("RF", "SVM", "k-NN", "GBM"), lineBest=F, interv = 8)

f1_3_2$total_features
f1_3_2$f1score
f1_3_2$accuracy
f1_3_2$all[[19]]

f1_3_3$all[[6]]

f1_3_1$all[[5]]

f1_3_5$total_features
f1_3_5$all[[4]]



ggplot(df_transpose_clusters, aes(x=PC1, y = PC2, color = factor(cluster))) + geom_point() + facet_wrap(~sample_type, scales = "free") 
ggplot(df_transpose_clusters, aes(x=PC1, y = PC2, color = factor(cluster))) + geom_point() + facet_wrap(~vital_status, scales = "free") 
ggplot(df_transpose_clusters, aes(x=PC1, y = PC2, color = factor(cluster))) + geom_point() + facet_wrap(~sample_type+vital_status , scales = "free") 

# res.pca.df_t <- as.data.frame(t(as.matrix(res.pca.df)))
# head(res.pca.df_t)
# ggplot(res.pca.df, aes(x=PC1, y = PC2, color = cluster)) + geom_point()




#--------------------------------------------------------------------------------------------
# Clinical Analysis of Selected Gene Probes and Samples
#--------------------------------------------------------------------------------------------



#--------------------------------------------------------------------------------------------


sd_df_filetered <- data_dge$counts[apply( data_dge$counts, 1 , function(x) sd(x) != 0 && all(x!=0)), ] #seleccionamos genes que tenga variacion
total_na(sd_df_filetered)
dim(sd_df_filetered)
any(sd_df_filetered==0)
log2_df <- apply(sd_df_filetered, 1:2, function(x) log(x, 2))
any(log2_df==0)
total_na(log2_df)
sd_df_filetered <- log2_df[apply( log2_df, 1 , function(x) sd(x) != 0 && all(x!=0)), ] #seleccionamos genes que tenga variacion
total_na(sd_df_filetered)
any(sd_df_filetered==0)
dim(sd_df_filetered)
# log2_df_no_zero_values <- log2_df[apply(log2_df, 1, function(x) all(x!=0)), ]
# total_na(log2_df_no_zero_values)
# sd_df_filetered <- na.omit(sd_df_filetered)
# data.quantileAll <- apply(log2_df, 2, function(x){quantile(x, 0.75)})
# data.variably<-apply(data_dge$counts, 1, sd)
# data.quantileExpressed <- apply(data_dge$counts, 1, function(x){quantile(x[x>0], 0.75)});

most_informative_genes<-genefilter::varFilter(sd_df_filetered)[1:1500,]
most_informative_genes <- as.data.frame(sweep(most_informative_genes,1, apply(most_informative_genes,1,median,na.rm=T)))
rownames(most_informative_genes)

#**********************
# gene names
#**********************
gene_ids<-c()
gene_names<-c()
rn <- rownames(most_informative_genes)
for(i in 1:length(rn)){
  gene_ids[i] <- rn[i]
  gene_names[i] <- data_dge$genes[data_dge$genes$ensembl_gene_id == rn[i],]$external_gene_name
}
selected_genes <- as.data.frame(cbind(gene_ids, gene_names))
head(selected_genes)
write.csv(selected_genes, "exported_data/TCGA-SKCM_transcriptomic_most_expressed_genes.csv", row.names = F)
#**********************
#cambiamos el nombre de las columnas (id genes) por el nombre de los genes
rownames(most_informative_genes) <- selected_genes$gene_names
head(most_informative_genes)
rownames(most_informative_genes)[1480:1500]


# 
# 
# prop_expressed = rowMeans(cpm(data_dge) > 1)
# keep = prop_expressed > 0.5
# op = par(no.readonly = TRUE)
# par(mfrow = c(1, 2))
# hist(cpm(data_dge, log = TRUE), main = 'Unfiltered', xlab = 'logCPM')
# abline(v = log(1), lty = 2, col = 2)
# hist(cpm(data_dge[keep, ], log = TRUE), main = 'Filtered', xlab = 'logCPM')
# abline(v = log(1), lty = 2, col = 2)
# par(op)
# #subset the data
# data_dge = data_dge[keep, , keep.lib.sizes = FALSE]
# data_se = data[keep, ]


# #------------------------------------------------------------
# # Transformation to FPKM values and normalisation
# #------------------------------------------------------------
# #download v22 of the GENCODE annotation
# gencode_file = 'gencode.v22.annotation.gtf.gz'
# gencode_link = paste(
#   'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22',
#   gencode_file,
#   sep = '/'
# )
# download.file(gencode_link, gencode_file, method = 'libcurl')
# gtf = import.gff(gencode_file, format = 'gtf', genome = 'GRCm38.71', feature.type = 'exon')
# #split records by gene to group exons of the same gene
# grl = reduce(split(gtf, elementMetadata(gtf)$gene_id))
# gene_lengths = ldply(grl, function(x) {
#   #sum up the length of individual exons
#   return(c('gene_length' = sum(width(x))))
# }, .id = 'ensembl_gene_id')
# 
# # Genes are also annotated with their biotype for further analysis. The annotation file uses Ensembl IDs with versions as keys to records, which then need to be converted to Ensembl IDs. This is simply achieved by truncating the trailing version number.
# #extract information on gene biotype
# genetype = unique(elementMetadata(gtf)[, c('gene_id', 'gene_type')])
# colnames(genetype)[1] = 'ensembl_gene_id'
# gene_lengths = merge(genetype, gene_lengths)
# 
# #remove ENSEMBL ID version numbers
# gene_lengths$ensembl_gene_id = gsub('\\.[0-9]*', '', gene_lengths$ensembl_gene_id)
# saveRDS(gene_lengths, file = 'gene_lengths_HTSeq_gencodev22.rds')
# 
# #allocate rownames for ease of indexing
# rownames(gene_lengths) = gene_lengths$ensembl_gene_id
# SummarizedExperiment::rowData(data_se)$gene_length = gene_lengths[rownames(data_se), 'gene_length']
# SummarizedExperiment::rowData(data_se)$gene_biotype = gene_lengths[rownames(data_se), 'gene_type']
# 
# #annotate gene lengths for the DGE object
# data_dge$genes$length = gene_lengths[rownames(data_dge), 'gene_length']
# 
# data_dge_tmm <- edgeR::calcNormFactors(data_dge, method = 'TMM')
# 
# #compute FPKM values and append to assays
# SummarizedExperiment::assay(data_se, 'logFPKM_TMM') = rpkm(data_dge_tmm, log = TRUE)
# data_se

#------------------------------------------------------------
# Selected Data
#------------------------------------------------------------

# data.quantileAll <- apply(SummarizedExperiment::assay(data_se), 2, function(x){quantile(x, 0.75)});
# data.quantileExpressed <- apply(SummarizedExperiment::assay(data_se), 2, function(x){quantile(x[x>0], 0.75)});
# data.norm <- as.data.frame(t(t(SummarizedExperiment::assay(data_se)) / data.quantileExpressed))
# 
# data.norm <- apply(data.norm, 1, mad)
# 
# 
# 
# # most_expressed_genes <- SummarizedExperiment::rowData(data_se)
# # most_expressed_genes<-as.data.frame(most_expressed_genes[1:1500,])
# # most_expressed_genes_ids <- most_expressed_genes$ensembl_gene_id
# most_expressed_genes_ids <- names(data.norm)[1:1500]


# #filtramos los datos por los genes mas expresados
# df<-SummarizedExperiment::assay(data)
# df <- as.data.frame(df[rownames(df) %in% most_expressed_genes_ids,])
# head(df)



head(most_informative_genes)
df_transpose <- as.data.frame(t(as.matrix(most_informative_genes))) #pacientes (rows) x genes (cols)
head(df_transpose)
write.csv(df_transpose, "exported_data/TCGA-SKCM_transcriptomic_most_expressed.csv", row.names = T)
df_aux <- read.csv("exported_data/TCGA-SKCM_transcriptomic_most_expressed.csv", row.names = "X")
# View(df_aux)



res.pca <- prcomp(df_transpose, scale = TRUE)
names(res.pca)
res.pca$x[1:5,1:5]
res.pca.df <- as.data.frame(res.pca$x[,1:5])
head(res.pca.df)
res.pca.df_t <- as.data.frame(t(as.matrix(res.pca.df)))
head(res.pca.df_t)


# factoextra::fviz_eig(res.pca)
# factoextra::fviz_pca_ind(res.pca,
#              col.ind = "cos2", # Color by the quality of representation
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#              repel = TRUE     # Avoid text overlapping
# )

# # Results for Variables
# res.var <- factoextra::get_pca_var(res.pca)
# res.var$coord          # Coordinates
# res.var$contrib        # Contributions to the PCs
# res.var$cos2           # Quality of representation 



#------------------------------------------------------------
# Clustering
#------------------------------------------------------------


# df_transpose_log2 <- as.matrix(log(df_transpose, 2))
# total_na(df_transpose_log2)

#median center genes
# dc = DelayedArray::sweep(df, 1, apply(df, 1, median, na.rm=T))
# total_na(dc)
# title=getwd()
# results = ConsensusClusterPlus(dc,maxK=5,reps=1000,pItem=0.8,pFeature=1,title=title,clusterAlg="hc",innerLinkage="average",
#                                finalLinkage="average",
#                                distance="pearson",seed=1262118388.71279,plot="png")

# same as above but with pre-computed distance matrix, useful for large datasets (>1,000's of items)
# dt = as.dist(1-cor(dc,method="pearson"))
# total_na(dt)

#--------------------------------------------------------
# requirement: genes in the rows ; patients in the columns
#--------------------------------------------------------
#data to be clustered; either a data matrix where columns=items/samples and
# rows are features. For example, a gene expression matrix of genes in rows and
# microarrays in columns, or ExpressionSet object, or a distance object (only for
#                                                                        cases of no feature resampling)
title=getwd()
rcc = ConsensusClusterPlus(as.matrix(res.pca.df_t),maxK=6,reps=200,pItem=0.8,pFeature=1,title=title,distance="pearson",clusterAlg="hc", seed=33, plot="png")
icl = calcICL(rcc,title=title,plot="png")

clusters_number=3
hclust_avg <- rcc[[clusters_number]]$consensusTree
clusters_class <- rcc[[clusters_number]]$consensusClass
table(clusters_class)
table(rcc[[3]]$consensusClass)
table(rcc[[4]]$consensusClass)
table(rcc[[5]]$consensusClass)
table(rcc[[2]]$consensusClass)


# set.seed(33)
# clusters_number <- 3
# seeds_df_sc <- as.data.frame(scale(df_transpose_log2))
# dist_mat <- dist(seeds_df_sc, method = 'euclidean')
# hclust_avg <- hclust(dist_mat, method = 'average')
# plot(hclust_avg)
# 
# cut_avg <- cutree(hclust_avg, k = clusters_number)
# plot(hclust_avg)
# rect.hclust(hclust_avg , k = clusters_number, border = 2:6)
# abline(h = clusters_number, col = 'red')
# 
# avg_dend_obj <- as.dendrogram(hclust_avg)
# avg_col_dend <- dendextend::color_branches(avg_dend_obj, h = clusters_number)
# plot(avg_col_dend)



factoextra::fviz_dend(x = hclust_avg, k = clusters_number, cex = 0.6) +
  geom_hline(yintercept = 5.5, linetype = "dashed") +
  labs(title = "Hierarchical clustering",
       subtitle = "pearson")


factoextra::fviz_cluster(object = list(data=df_transpose_log2, cluster=cutree(hclust_avg, k=clusters_number)),
             ellipse.type = "convex", repel = TRUE, show.clust.cent = FALSE,
             labelsize = 8)  +
  labs(title = "Hierarchical clustering",
       subtitle = "pearson") +
  theme_bw() +
  theme(legend.position = "bottom")



# 
# 
# res.pca <- prcomp(df_transpose, scale = TRUE)
# names(res.pca)
# res.pca$x[1:5,1:5]
# res.pca.df <- as.data.frame(res.pca$x[,1:5])
# head(res.pca.df)
# res.pca.df_t <- as.data.frame(t(as.matrix(res.pca.df)))

# adding cluster column
df_transpose_clusters <-res.pca.df
df_transpose_clusters$cluster <- clusters_class
df_transpose_clusters<- df_transpose_clusters[order(df_transpose_clusters$cluster), ]
head(df_transpose_clusters)

most_informative_genes_clusters <- as.data.frame(t(as.matrix(most_informative_genes)))
most_informative_genes_clusters$cluster <- clusters_class
head(most_informative_genes_clusters)


ggplot(df_transpose_clusters, aes(x=PC1, y = PC2, color = factor(cluster))) + geom_point()



# adding extra columns with patient information
sample_type<-c()
tumor_stage<-c()
initial_weight<-c()
year_of_birth<-c()
year_of_death <- c()
gender <- c()
vital_status <- c()
rn<-rownames(df_transpose_clusters)
for( i in 1:nrow(df_transpose_clusters)){
  sample_type[i]<-original_data[original_data$barcode == rn[i],]$sample_type
  tumor_stage[i]<-original_data[original_data$barcode == rn[i],]$tumor_stage
  initial_weight[i]<-original_data[original_data$barcode == rn[i],]$initial_weight
  year_of_birth[i]<-original_data[original_data$barcode == rn[i],]$year_of_birth
  year_of_death[i]<-original_data[original_data$barcode == rn[i],]$year_of_death
  gender[i] <- original_data[original_data$barcode == rn[i],]$gender
  vital_status[i] <- original_data[original_data$barcode == rn[i],]$vital_status
}
df_transpose_clusters$sample_type <- sample_type
df_transpose_clusters$tumor_stage <- tumor_stage
df_transpose_clusters$initial_weight <- initial_weight
df_transpose_clusters$year_of_birth <- year_of_birth
df_transpose_clusters$year_of_death <- year_of_death
# df_transpose_clusters$is_dead <- ifelse(!is.na(df_transpose_clusters$year_of_death), "dead", "alive")
df_transpose_clusters$vital_status <- vital_status
df_transpose_clusters$gender <- gender




#sort columns, put the last ones at the first position
c1 <- ncol(df_transpose_clusters)-7 #cluster
c2 <- ncol(df_transpose_clusters)-6 #sample_type
c3 <- ncol(df_transpose_clusters)-5 #tumor_stage
c4 <- ncol(df_transpose_clusters)-4 #initial_weight
c5 <- ncol(df_transpose_clusters)-3 #year_of_birth
c6 <- ncol(df_transpose_clusters)-2 #year_of_death
c7 <- ncol(df_transpose_clusters)-1 #vital_status
c8 <- ncol(df_transpose_clusters)   #gender
start_cols <- c(c1,c2,c3,c4,c5,c6,c7,c8)
col_rest <- 1:(ncol(df_transpose_clusters)-length(start_cols))
df_transpose_clusters <- df_transpose_clusters[, c(start_cols, col_rest)] 
write.csv(df_transpose_clusters, "TCGA-SKCM_transcriptomic_most_expressed_clusters.csv", row.names = T)


ggplot(df_transpose_clusters, aes(x=PC1, y = PC2, color = factor(cluster))) + geom_point() + facet_wrap(~tumor_stage, scales = "free") 



ggplot(df_transpose_clusters, aes(x=tumor_stage, y = CACNA2D2, color = factor(cluster))) + geom_point() + facet_wrap(~vital_status, scales = "free") 

ggplot(df_transpose_clusters, aes(x=tumor_stage, y = XPO1, color = factor(cluster))) + geom_point() + facet_wrap(~vital_status, scales = "free") 
ggplot(df_transpose_clusters, aes(x=tumor_stage, y = XPO1, color = factor(cluster))) + geom_point()
ggplot(df_transpose_clusters, aes(x=tumor_stage, y = CELF2, color = factor(cluster))) + geom_point()
ggplot(df_transpose_clusters, aes(x=tumor_stage, y = POLR2B, color = factor(cluster))) + geom_point()
ggplot(df_transpose_clusters, aes(x=tumor_stage, y = EPS15, color = factor(cluster))) + geom_point()



#--------------------------------------------------------
# Top Genes expresados por cada cluster
#--------------------------------------------------------
top_genes_by_cluster <- most_informative_genes_clusters
# cluster1 <- top_genes_by_cluster[top_genes_by_cluster$cluster == 1, c((length(start_cols)+1):(ncol(top_genes_by_cluster)))]
# cluster1_t <- as.data.frame(t(as.matrix(cluster1)))
# cluster2 <- top_genes_by_cluster[top_genes_by_cluster$cluster == 2, c((length(start_cols)+1):(ncol(top_genes_by_cluster)))]
# cluster2_t <- as.data.frame(t(as.matrix(cluster2)))
# cluster3 <- top_genes_by_cluster[top_genes_by_cluster$cluster == 3, c((length(start_cols)+1):(ncol(top_genes_by_cluster)))]
# cluster3_t <- as.data.frame(t(as.matrix(cluster3)))

cluster1 <- top_genes_by_cluster[top_genes_by_cluster$cluster == 1, 1:(ncol(top_genes_by_cluster)-1)]
cluster1_t <- as.data.frame(t(as.matrix(cluster1)))
cluster2 <- top_genes_by_cluster[top_genes_by_cluster$cluster == 2, 1:(ncol(top_genes_by_cluster)-1)]
cluster2_t <- as.data.frame(t(as.matrix(cluster2)))
cluster3 <- top_genes_by_cluster[top_genes_by_cluster$cluster == 3, 1:(ncol(top_genes_by_cluster)-1)]
cluster3_t <- as.data.frame(t(as.matrix(cluster3)))


# cluster_1_s <- apply(cluster1_t, 1, sum)
# cluster_2_s <- apply(cluster2_t, 1, sum)
# cluster_3_s <- apply(cluster3_t, 1, sum)
# cluster1_top10=cluster_1_s[rev(order(cluster_1_s))[1:30]]
# cluster2_top10=cluster_2_s[rev(order(cluster_2_s))[1:30]]
# cluster3_top10=cluster_3_s[rev(order(cluster_3_s))[1:30]]
# cluster_1_top_names <- names(cluster1_top10)
# cluster_2_top_names <- names(cluster2_top10)
# cluster_3_top_names <- names(cluster3_top10)

#seleccion por measured by median absolute deviation
cluster1_mads=apply(cluster1_t,1,mad)
cluster1_top10=cluster1_t[rev(order(cluster1_mads))[1:30],]
cluster2_mads=apply(cluster2_t,1,mad)
cluster2_top10=cluster2_t[rev(order(cluster2_mads))[1:30],]
cluster3_mads=apply(cluster3_t,1,mad)
cluster3_top10=cluster3_t[rev(order(cluster3_mads))[1:30],]
cluster_1_top_names <- rownames(cluster1_top10)
cluster_2_top_names <- rownames(cluster2_top10)
cluster_3_top_names <- rownames(cluster3_top10)


top10_genes <- data.frame(cluster1=cluster_1_top_names, cluster2=cluster_2_top_names, cluster3=cluster_3_top_names)
top10_genes
write.csv(top10_genes, "exported_data/TCGA-SKCM_top10genes_by_cluster.csv", row.names = F)
rownames(cluster3_top10)
rownames(cluster3_top10)


ggplot(df_transpose_clusters, aes(x=tumor_stage, y = CACNA2D2, color = factor(cluster))) + geom_point() + facet_wrap(~is_dead, scales = "free") + facet_wrap(~gender, scales = "free") 


table(df_transpose_clusters$cluster)

#nombres de los genes
gene_ids<-c()
gene_names<-c()
cn <- colnames(df_transpose_clusters)[length(start_cols)+1:length(colnames(df_transpose_clusters))]
for(i in 1:length(cn)){
  gene_ids[i] <- cn[i]
  gene_names[i] <- data_dge$genes[data_dge$genes$ensembl_gene_id == cn[i],]$external_gene_name
}
selected_genes <- as.data.frame(cbind(gene_ids, gene_names))
head(selected_genes)
write.csv(selected_genes, "exported_data/TCGA-SKCM_transcriptomic_most_expressed_genes.csv", row.names = F)

#------------------------------------------------------------


seeds_df_cl <- mutate(df_transpose, cluster = cut_avg)
# count(seeds_df_cl, 3)
ggplot(seeds_df_cl, aes(x=area, y = perimeter, color = factor(cluster))) + geom_point()


set.seed(786)

class(my_data)

dim(df)

# data <- GDCprepare(query, summarizedExperiment = F)
write.csv(df, "TCGA-SKCM_transcriptomic_most_expressed.csv", row.names = F)

df <- read.csv("TCGA-SKCM_transcriptomic.csv")

hclust <- TCGAanalyze_Clustering(df, "hclust", methodHC = "ward.D2")
clustconsensus <- TCGAanalyze_Clustering(df, "consensus", methodHC = "ward.D2")

plot(hclust)
