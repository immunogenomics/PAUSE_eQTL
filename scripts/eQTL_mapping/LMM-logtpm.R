#!/usr/bin/env Rscript

all <- commandArgs(trailingOnly = TRUE)

#Specify the chromosome number and number of RNA-PCs to include in the model
chr <- all[1]
PC_num <- all[2]

#Load libraries
library(tidyverse)
library(lme4)
library(lmerTest)
library(parallel)
library(future.apply)
library(foreach)
library(doParallel)

#Load gene expression, metadata, genotype dosages data and PCs
RNAmatrixLM <- readRDS(paste0("log2tpm_final_chr", chr,".rds")) # scaled and QN RNA
MetaTable <- readRDS("meta_final_pca_w_2942genes_50pcs-2-10-21.rds")
GenePairs <- readRDS(paste0("GenePairs_noMHC_chr",chr,".rds"))
dosagefiltered_DR <- readRDS(paste0("dosage_filterednoMHC_chr",chr,".rds"))
covariatesPCsLM <- readRDS(paste0("covariate-files/covariate_",PC_num,"ePCs3gPCsSite-2-11-21.rds"))

#Define PCs covariates
coFixedLM <- paste(paste("covariatesPCsLM$", colnames(covariatesPCsLM), sep=""), collapse=" + ")
donor <- as.factor(MetaTable$USUBJID)

#Run analysis
res <- apply(GenePairs,1,function(pair){

    gene <- as.numeric(RNAmatrixLM[pair[1],]) #gene expression vector
#     print(length(gene))
    snp <- as.numeric(dosagefiltered_DR[pair[2],]) # snp genotypes as dosage
#     print(length(snp))
    form1 <- as.formula(paste("gene","~","snp", "+", coFixedLM, "+ (1|donor)")) #Linear model. Full

    lm1 <- lmer(form1)

    test <- summary(lm1)
    test$coefficients[2,]

    })

colnames(res) <- GenePairs$ID
res <- res %>% t() %>% as_tibble(rownames = "ID")
res$Gene <- GenePairs$gene
res <- res %>% arrange(`Pr(>|t|)`)

saveRDS(res, paste0("LMM-output-logtpm-chr",chr,'_', PC_num, "PC.rds"))
