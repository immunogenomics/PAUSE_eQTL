#Load libraries
suppressPackageStartupMessages({
    library(tidyverse)
    library(lme4)
    library(vcfR)
    library(ggthemes)
    library(preprocessCore)
    library(lmerTest)
    library(foreach)
    library(doParallel)
    library(future.apply)
    library(janitor)
    })

#Load data
RNAmatrixLM <- readRDS("/path/log2tpm_full_final_w_27100genes.rds") # scaled and QN RNA
MetaTable <- readRDS('/path/2022-05-26_meta_final_w_IL17score.rds')
sigeQTLall <- readRDS("/path/sigeGenesallCHROMlogtpmpasslowexp.rds")
dosagefiltered_DR <- readRDS("/path/dosage_filterednoMHC_full.rds")
covariatesPCsLM <- readRDS("/path/covariate_20ePCs3gPCsSite-2-11-21.rds")


GenePairs <- dplyr::select(sigeQTLall, c("Gene", "ID"))

#Define PCs covariates
coFixedLM <- paste(paste("covariatesPCsLM$", colnames(covariatesPCsLM), sep=""), collapse=" + ")

donor <- as.factor(MetaTable$USUBJID)
type <- as.factor(MetaTable$SPITS_class)#or any other variable we want to test interactions

#Run models
suppressMessages({res <- apply(GenePairs,1,function(pair){

    gene <- as.numeric(RNAmatrixLM[pair[1],]) #gene expression vector
#     print(length(gene))
    snp <- as.numeric(dosagefiltered_DR[pair[2],]) # snp genotypes as dosage
#     print(length(snp))
    form0 <- as.formula(paste("gene","~","snp", "+", coFixedLM, "+ (1|donor) + type")) #Linear model. Null
    
    lm0 <- lmer(form0)
    
    form1 <- as.formula(paste("gene","~","snp", "+", coFixedLM, "+ (1|donor) + type + type*snp")) #Linear model. Full
    
    lm1 <- lmer(form1)
    test <- anova(lm0,lm1, refit = F)
    test[2,]

    })

                     })
res_final <- data.frame(matrix(unlist(res), nrow=length(res), byrow=T))

colnames(res_final) <- c("npar", "AIC", "BIC", "logLik", "deviance", "Chisq", "Df", "Pr(>Chisq)")
res_final$ID <- GenePairs$ID
res_final$Gene <- GenePairs$Gene
res_final <- res_final %>% arrange(`Pr(>Chisq)`)

res_final %<>% mutate("FDR" = p.adjust(`Pr(>Chisq)`)) %>% arrange(FDR)

#Run models
 suppressMessages({res <- apply(GenePairs,1,function(pair){

    gene <- as.numeric(RNAmatrixLM[pair[1],]) #gene expression vector

    snp <- as.numeric(dosagefiltered_DR[pair[2],]) # snp genotypes as dosage

    form0 <- as.formula(paste("gene","~","snp", "+", coFixedLM, "+ (1|donor) + type")) #Linear model. Null

    
    lm0 <- lmer(form0)
    
    form1 <- as.formula(paste("gene","~","snp", "+", coFixedLM, "+ (1|donor) + type + type*snp")) #Linear model. Full
    
    lm1 <- lmer(form1)
    test <- summary(lm1) 
    test$coefficients[c(2,35,36),]
    

    })

                     })

res <- res %>% t() %>% as.data.frame()

colnames(res) <- c(paste0( c('eQTL_', 'Gene_', 'Interaction_'),'Beta'),
                   paste0( c('eQTL_', 'Gene_', 'Interaction_'),'SE'),
                   paste0( c('eQTL_', 'Gene_', 'Interaction_'),'df'),
                   paste0( c('eQTL_', 'Gene_', 'Interaction_'),'tValue'),
                   paste0( c('eQTL_', 'Gene_', 'Interaction_'),'PValue'))
res$ID <- GenePairs$ID

res$Gene <- GenePairs$Gene

res <- res %>% arrange(Interaction_PValue)
