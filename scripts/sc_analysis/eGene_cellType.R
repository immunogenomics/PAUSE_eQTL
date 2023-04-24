#Load libraries
library(edgeR)
library(readr)
library(data.table)
library(tximport)
library(stringr)
library(readr)
library(biomaRt)


#Load sc data
cpm_M <- readRDS('/path/2021-10-06_scMatrix_intx_genes_CPM_CellTypeNames_wGSTT2wAC.rds')

#Categorize cell type
Category = c()
cell_type = unlist(colnames(cpm_M))

for (i in seq(1:length(cell_type)) ){
    
    if(cell_type[i] %in% c('MigDC', 'DC1', 'DC2', 'moDC_1',  'moDC_2', 'moDC_3') ){
        Category <- append(Category, 'DC')
    }else if(cell_type[i] %in% c('Tc17_Th17','Treg', 'Th', 'Tc') ){
        Category <- append(Category, 'T_cell')
    }else if(cell_type[i] %in% c('F1', 'F2', 'F3')){
        Category <- append(Category, 'Fibroblast')
    }else if(cell_type[i] %in% c('Differentiated_KC*', 'Differentiated_KC', 'Undifferentiated_KC', 'Proliferating_KC')){
        Category <- append(Category, 'KC') 
    }else if(cell_type[i] %in% c('LC_1', 'LC_2', 'LC_3', 'LC_4')){
        Category <- append(Category, 'LC') 
    }else if(cell_type[i] %in% c('Macro_1', 'Macro_2', 'Inf_mac', 'Mono_mac')){
        Category <- append(Category, 'Macrophage')
    }else if(cell_type[i] %in% c('ILC1_NK', 'ILC1_3', 'ILC2', 'NK')){
        Category <- append(Category, 'ILC_NK')
    }else if(cell_type[i] %in% c('VE1', 'VE2', 'VE3')){
        Category <- append(Category, 'VE')
    }else if(cell_type[i] %in% c('Pericyte_1', 'Pericyte_2')){
        Category <- append(Category, 'Pericyte')
    }else if(cell_type[i] %in% c('Melanocyte')){
        Category <- append(Category, 'Melanocyte')
    
    }else if(cell_type[i] %in% c('VE1', 'VE2', 'VE3')){
        Category <- append(Category, 'VE')
        
    }else if(cell_type[i] %in% c('LE1', 'LE2')){
        Category <- append(Category, 'LE')
   
    }else if(cell_type[i] %in% c('Plasma')){
        Category <- append(Category, 'Plasma')
        
    }else if(cell_type[i] %in% c('Mast_cell')){
        Category <- append(Category, 'Mast_cell')
    
    }else if(cell_type[i] %in% c('Schwann_1', 'Schwann_2')){
        Category <- append(Category, 'Schwann')
    }else if(cell_type[i] %in% c('Gene')){
        Category <- append(Category, 'Genes')}
    

        }
        
        
#Normalize so that rows add up to 1
egene_norm = sweep(cpm_M, 1, rowSums(cpm_M), FUN = '/')
        
df_egene_norm <- as.data.frame(as.matrix(egene_norm)) 
        
        
colnames(df_egene_norm) <- Category
        

heat_egene <- as.data.frame( 
    sapply(unique(names(df_egene_norm)), # for each unique column name
       function(col) rowMeans(df_egene_norm[names(df_egene_norm) == col]) # calculate row means
    )
  )
        
