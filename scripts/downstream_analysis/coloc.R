#Load eQTL data

eqtl_summ <- readRDS('/path/cispassresults-lowexp_allchrom.rds')
genotype <- readRDS("/path/filtered_snps2-final-7-23-21.rds")

our.data <- left_join(eqtl_summ, genotype, by = 'ID')
our.data <- our.data %>%
mutate(id2 = paste0(CHROM.y, '_', POS, '_', ALT, '_', REF ))

#Load GWAS data
##eczema 
ez_gwas <- read.table( gzfile('/data/srlab2/qxiao/ITN/data/external/eczema_build37_summstats.tsv.gz'), sep = '\t', header = TRUE)

ez_gwas <- ez_gwas %>% 
mutate(ID = paste0(chromosome, '_', base_pair_location, '_', effect_allele, '_', other_allele)) 

#Get eQTL loci
end.data <- read.csv('/path/Table_S2.csv')
Loci.of.Interest <- end.data$Ensembl_ID

#Prepare data for coloc
our.data <- our.data %>% separate(Gene, into = c('Gene_short', 'Transcript'), sep = '[.]', remove = FALSE)
our.data.loi <- our.data %>% filter(Gene_short %in% end.data$Ensembl_ID)

#Run coloc
gwas_summstats <- ez_gwas


nsnps <- c()
loi <- c()
h3 <- c()
h4 <- c()


suppressWarnings({



for (i in 1:length(Loci.of.Interest)){
    
    our.data.onelocus <- our.data.loi %>% filter( Gene_short == Loci.of.Interest[i] )
    
    test1 <- inner_join(our.data.onelocus, gwas_summstats, by = 'ID')
    test2 <- inner_join(our.data.onelocus, gwas_summstats, by = c('id2'='ID') )
    
  
    test_combined <- rbind(test1, test2)
    
    if (nrow(test_combined) > 0) {
        result.onegene <- coloc.abf(
                              
                                dataset1=list(type="cc", beta=test_combined$beta, varbeta=(test_combined$standard_error)**2), 
                                dataset2=list(pvalues=test_combined$`Pr(>|t|)`, type="quant",beta=test_combined$Estimate, varbeta=(test_combined$`Std. Error`)**2,  N=375), 
                                MAF=test_combined$MAF) 
    
       
        loi <- append(loi, Loci.of.Interest[i])
        nsnps <- append( nsnps, result.onegene$summary['nsnps'])
        h3 <- append( h3, result.onegene$summary['PP.H3.abf'] )
        h4 <- append( h4, result.onegene$summary['PP.H4.abf'])
        
        }else{
        
        loi <- append(loi, NA)
        nsnps <- append( nsnps, NA)
        h3 <- append( h3, NA)
        h4 <- append( h4, NA)
    
    }
    
    
    
    }
    
})

#Write Table

df.coloc.ez <- data.frame(
    Locus = loi,
    nSNPs = nsnps,
    PP.H3 = h3,
    PP.H4 = h4,
    GWAS.Trait = 'Eczema'
)

df.coloc.ez <- left_join( df.coloc.ez, end.data[,c('Ensembl_ID', 'Symbol')], by = c('Locus'='Ensembl_ID') )

df.coloc.ez2 <- df.coloc.ed %>% filter(PP.H4 > 0.75) #Filter HH4 > 0.75

