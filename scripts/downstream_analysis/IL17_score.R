#Load libraries
library("readxl")
library(msigdbr)

#Load data
tpm <- readRDS('../data/2022-05-26-tpm_by_GeneName.rds') #PAUSE gene TPM
subject <- read_excel("../data/external/SHUNSUKE-KC-IL17andIL36.MA.KC.BigTable.DEGs.FCH1.5.p05.xlsx", sheet = 1)#KC IL-17 genes


#KC IL17 score
IL17 <- subject$SYMBOL[which( abs( subject$FCH_IL17A.vs.Control ) > 5 )]
uni_IL17 <- unique(IL17)

uni_IL17 <- c(uni_IL17, 'IL17A')

rownames(tpm)[which(rownames(tpm) == 'VNN3P' )] <- 'VNN3'


df.IL17 <- tpm[uni_IL17,]

df.IL17 <- df.IL17[ complete.cases(df.IL17),]

IL17_score <- colSums( df.IL17[-23,] )

IL17_score <- scale(IL17_score)


#MSigDB IL17 score
g17 <- read.table(file = paste0( '../data/external/WP_IL17_SIGNALING_PATHWAY.v7.5.1.tsv'), sep = '\t', header = TRUE)
g17 <- unlist( strsplit( g17[g17$STANDARD_NAME == 'MAPPED_SYMBOLS','WP_IL17_SIGNALING_PATHWAY'], ',') )

df.g17 <- tpm[g17,]

df.g17 <- df.g17[ complete.cases(df.g17),]

g17_score <- scale( colSums( df.g17 ))

