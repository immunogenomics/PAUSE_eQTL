#Load libraries
library(MASS)
library(ggplot2)

#Load data
meta <- readRDS('/path/meta_final_pca_w_2942genes_50pcs_fix_relapse-4-28-21.rds')

pca_obj <- readRDS('/path/pca_var_labs_w_2942_genes_0.7_0.7-2-10-21.rds')


#Determine the number of PCs to use for LDA
df_scree <- data.frame(
             PC = paste('PC', seq(375), sep = ''),
             var_explained = ls_perc,
             cum_sum = cumsum(ls_perc))

which(df_scree$cum_sum > 90)[1] #48 PCs

temp <- meta

#Scale the predictors(PC1-48) 
temp[,49:96] <- scale(temp[,49:96])

#Split train and test sets based on visit numbers
first_visit <- temp$Visit.Number == 0

train <- temp[ first_visit,c(5,49:96)]
test <- temp[ !first_visit,c(5,49:96)]

#Reset the outcome variable as factor
train$Biopsy.Type <- factor(train$Biopsy.Type, levels = c('Non-Lesional', 'Lesional'))

#LDA fitting
model <- lda(Biopsy.Type ~ ., data=train)

#Predict the data 
predicted <- predict(model, test)

SPITS <- predicted$x

SPITS_class <- ifelse( predicted$x > 0, 'SPITS positive', 'SPITS negative') 


