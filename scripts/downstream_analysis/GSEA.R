#Load libraries
library(dplyr)
library(msigdbr)
library(fgsea)
library(tibble)

#Load data
res_final_all <- readRDS( '/path/2021-07-19_LMM-lowexp-eQTL-LDATypeint-Output-allCHROMFwHGNCandclassification.rds')
res <- readRDS( "/path/2021-07-18-LMM-lowexp-eQTL-LDAType-lm-classification-allCHROMF.rds")

#Score calculated by estimate*-log10(FDR)
res$Rank_score <- abs(res$Estimate) * (-log10(res$`Pr(>|t|)`))
res$gene.short <- lapply(strsplit(res$Gene, "[.]"), "[[", 1)

##format for gsea 
FC.vec <- abs( res$Estimate )
names(FC.vec) <- res$gene.short

scoreType <- "pos"

#Curated gene sets(or can be other gene sets, such as H for hallmark)
c2_gene_sets = msigdbr(species = "human", category = 'C2')

C2.ensembl.ls <- c2_gene_sets %>% 
  select(gs_name, ensembl_gene) %>% 
  group_by(gs_name) %>% 
  summarise(all.genes = list(unique(ensembl_gene))) %>% 
  deframe()

#Run GSEA
gsea.c2 <- fgseaSimple(pathways = C2.ensembl.ls,
                      stats = FC.vec,
                      scoreType = scoreType,
                      nperm=10000)
