library(WGCNA)
library(dplyr)

set.seed(123)

####constructing the COE network using the "supplement_table1_expression_matrix.csv"
##Add column headings to row names for ease of operation
exprmatr <- read.csv("supplement_table1_expression_matrix.csv",header = F,sep = "\t")
exprmatr[1:5,1:5]
###Clear the line, import exprmatr, and assign the value of Exprmat_select with the line name
exprmatr_select <- exprmatr[-1,-1]
exprmatr_select[1:5,1:5]
exprmatr <- read.csv("supplement_table1_expression_matrix.csv",header = T,sep = "\t")
colnames(exprmatr_select) <- colnames(exprmatr)
exprmatr_select$Sample <- rownames(exprmatr)
exprmatr_select <- select(exprmatr_select,941,1:940)
exprmatr_select[1:5,1:5]
exprmatr <- exprmatr_select
genes_exprmatr <- substring(exprmatr$Sample,1)
# Transpose the matrix
expr_matrix_transposed <- t(exprmatr)
# Set the column names to the gene names
colnames(expr_matrix_transposed) <- exprmatr$Sample
# View the transposed matrix
expr_matrix_transposed[1:10,1:10]
expr_matrix_transposed <- expr_matrix_transposed[-1,]
#Save the results for later
#write.table(expr_matrix_transposed,file = "expr_matrix_transposed.txt",sep = "\t",quote = F)
write.table(exprmatr,"matrix_before_PCA.txt",col.names = T,row.names = F,quote = F,sep = "\t")
#COE_network construction
exprmatr_matrix <- expr_matrix_transposed
set.seed(123)
TF_gene <-read.table("filter_genes_plantTFDB_filtered.txt",header = T,sep = "\t") 
genes <- substring(TF_gene$TF_name,1)
positions <- ifelse(genes_exprmatr %in% genes, TRUE, FALSE)
adj <- adjacency(exprmatr_matrix,selectCols = positions,power = 6)
##write.table(adj,"adj_pearson.txt",quote = F,sep = "\t")
#Filter the matrix
#Delete a row in the matrix whose values are all NA
delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]
}
filter_adj <- delete.na(adj,2310)
#Delete a column in the matrix whose values are all NA
x_RowsAllNa_removed =  filter_adj[apply(filter_adj, 1, function(y) any(!is.na(y))),]
#replace NA to 0
filter_adj[is.na(filter_adj)]=0
#The first thousand co-expressed genes were selected as targets
# Initialize list to hold dataframes
df_list <- list()
# Loop over columns and extract top thousand rows
for (i in 1:ncol(filter_adj)) {
  df_list[[i]] <- data.frame(filter_adj[order(filter_adj[,i], decreasing = TRUE)[1:1000], i])
  rownames(df_list[[i]]) <- rownames(filter_adj)[order(filter_adj[,i], decreasing = TRUE)[1:1000]]
  colnames(df_list[[i]]) <- colnames(filter_adj)[i]
}
# Initialize list to hold transformed dataframes
transformed_df_list <- list()
# Loop over dataframes and transform each one
for (i in 1:length(df_list)) {
  transformed_df <- data.frame(colnames(df_list[[i]]), rownames(df_list[[i]]), select(df_list[[i]], 1))
  rownames(transformed_df) <- NULL
  colnames(transformed_df) <- c("TF_name","TARGET","scores")
  transformed_df_list[[i]] <- transformed_df
}
# Merge all dataframes by row
newdata <- bind_rows(transformed_df_list)
#Get the subset where scores is not zero
pearson <- subset(newdata,scores !=0)
#They are listed in descending order by scores
pearson_filter <- pearson[order(pearson$scores,decreasing = T),]
#save
#Remove the column where TF and target are the same gene
pearson_filter <- subset(pearson_filter,pearson_filter$TF_name != pearson_filter$TARGET)
write.table(pearson_filter,"COE_network.txt",sep = "\t",quote = F,row.names = F)








