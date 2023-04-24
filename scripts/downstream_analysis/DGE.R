#Load libraries
library(edgeR)
library(readr)
library(data.table)
library(tximport)
library(stringr)
library(readr)
library(biomaRt)


#Load data
txi <- readRDS('/path/2021-07-21-RNA-seq_transcript_level_quant_final.rds')
meta <- readRDS( '/path/2021-07-01_meta_final_pca_w_2942genes_50pcs_LDA_w_48PC.rds')

#Annotate raw count
mart <- useMart("ensembl")
datasets <- listDatasets(mart)
mart <- useDataset("hsapiens_gene_ensembl",mart)

genes_hgnc <- getBM(attributes = c("ensembl_gene_id",'ensembl_transcript_id', 
                                   "external_gene_name", "hgnc_id"),    
                      mart = mart)

#Raw count matrix
RNAmatrix <- txi$counts

transcript_id_short <- lapply(strsplit(rownames(RNAmatrix), "[.]"), "[[", 1) %>% as.character()

rownames(RNAmatrix) <- transcript_id_short

RNAmatrix_annot <- merge(RNAmatrix, genes_hgnc[, c('ensembl_gene_id', 'external_gene_name', 'ensembl_transcript_id')], 
                         by.x = 0, by.y = 'ensembl_transcript_id', all.x = TRUE )
RNAmatrix_annot <- unique(RNAmatrix_annot)

qc_gene_id <- lapply(strsplit(rownames(RNAmatrixLM), "[.]"), "[[", 1) %>% as.character()

RNAmatrix_qc <- RNAmatrix_annot[ RNAmatrix_annot$ensembl_gene_id %in% qc_gene_id, 
                                c(colnames(RNAmatrixLM), 'external_gene_name')]


test <- aggregate(. ~ external_gene_name, data=RNAmatrix_qc, FUN=sum)
df_gene_count <- test[-1,]
df_gene_count2 <- df_gene_count[,-1]
rownames(df_gene_count2) <- df_gene_count[,1]


#DGE analysis

##Rearrange matrix for this comparison

deg_mat <- df_gene_count2

biop_pred <-  meta$SPITS_class

dgList <- DGEList(counts=deg_mat, genes=rownames(deg_mat))

#Keep the genes that have counts >6 in at least 75 samples
keep <- rowSums(cpm(dgList) > 1 ) >= 75

dgList <- dgList[keep, keep.lib.sizes = FALSE]

designMat <- model.matrix(~ biop_pred)

#Normalize the data
dgList <- calcNormFactors(dgList, method = 'TMM')


dgList <- estimateGLMCommonDisp(dgList, design=designMat)
dgList <- estimateGLMTrendedDisp(dgList, design=designMat)
dgList <- estimateGLMTagwiseDisp(dgList, design=designMat)

fit <- glmFit(dgList, designMat)


up_index <- which(lrt$table$FDR < 0.05 & lrt$table$logFC > 1.5)
down_index <- which(lrt$table$FDR < 0.05 & lrt$table$logFC < -1.5 )
