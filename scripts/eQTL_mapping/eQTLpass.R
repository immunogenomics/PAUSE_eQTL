

#Load LMM results from each chromosome
chr1 <- readRDS("../LMM-output-logtpm-chr1-2-11-21.rds")
chr1$CHROM <- "chr1"
chr2 <- readRDS("../LMM-output-logtpm-chr2-2-11-21.rds")
chr2$CHROM <- "chr2"
chr3 <- readRDS("../LMM-output-logtpm-chr3-2-11-21.rds")
chr3$CHROM <- "chr3"
chr4 <- readRDS("../LMM-output-logtpm-chr4-2-11-21.rds")
chr4$CHROM <- "chr4"
chr5 <- readRDS("../LMM-output-logtpm-chr5-2-11-21.rds")
chr5$CHROM <- "chr5"
chr6 <- readRDS("../LMM-output-logtpm-chr6-2-11-21.rds")
chr6$CHROM <- "chr6"
chr7 <- readRDS("../LMM-output-logtpm-chr7-2-11-21.rds")
chr7$CHROM <- "chr7"
chr8 <- readRDS("../LMM-output-logtpm-chr8-2-11-21.rds")
chr8$CHROM <- "chr8"
chr9 <- readRDS("../LMM-output-logtpm-chr9-2-11-21.rds")
chr9$CHROM <- "chr9"
chr10 <- readRDS("../LMM-output-logtpm-chr10-2-11-21.rds")
chr10$CHROM <- "chr10"
chr11 <- readRDS("../LMM-output-logtpm-chr11-2-11-21.rds")
chr11$CHROM <- "chr11"
chr12 <- readRDS("../LMM-output-logtpm-chr12-2-11-21.rds")
chr12$CHROM <- "chr12"
chr13 <- readRDS("../LMM-output-logtpm-chr13-2-11-21.rds")
chr13$CHROM <- "chr13"
chr14 <- readRDS("../LMM-output-logtpm-chr14-2-11-21.rds")
chr14$CHROM <- "chr14"
chr15 <- readRDS("../LMM-output-logtpm-chr15-2-11-21.rds")
chr15$CHROM <- "chr15"
chr16 <- readRDS("../LMM-output-logtpm-chr16-2-11-21.rds")
chr16$CHROM <- "chr16"
chr17 <- readRDS("../LMM-output-logtpm-chr17-2-11-21.rds")
chr17$CHROM <- "chr17"
chr18 <- readRDS("../LMM-output-logtpm-chr18-2-11-21.rds")
chr18$CHROM <- "chr18"
chr19 <- readRDS("../LMM-output-logtpm-chr19-2-11-21.rds")
chr19$CHROM <- "chr19"
chr20 <- readRDS("../LMM-output-logtpm-chr20-2-11-21.rds")
chr20$CHROM <- "chr20"
chr21 <- readRDS("../LMM-output-logtpm-chr21-2-11-21.rds")
chr21$CHROM <- "chr21"
chr22 <- readRDS("../LMM-output-logtpm-chr22-2-11-21.rds")
chr22$CHROM <- "chr22"

eQTLall <- rbind(chr1, chr2, chr3, chr4,
                 chr5, chr6, chr7, chr8,
                 chr9, chr10, chr11, chr12,
                 chr13, chr14, chr15, chr16,
                 chr17, chr18, chr19, chr20,
                 chr21, chr22)

eQTLall$Gene_ID <- lapply(strsplit(eQTLall$Gene, "[.]"), "[[", 1)

#Deduplicate the data by Ensembl ID, picking the lead SNP-gene pair with the largest estimate
dedup <- eQTLall[!duplicated(eQTLall$Gene),]

#Bonferroni Threshold
pvalue <- 0.05/6305752

#Record whether a SNP-Gene pair is significant
dedup$sig <- "no"

dedup[dedup$`Pr(>|t|)` < pvalue,]$sig <- "yes"

saveRDS(dedup, "../cispassresults-lowexp_allchrom_dedup.rds")

sig <- dedup[dedup$sig == "yes",]

saveRDS(sig, "../sigeGenesallCHROMlogtpmpasslowexp.rds")
