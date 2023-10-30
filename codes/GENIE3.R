library(GENIE3)
set.seed(123)
exprDataframe <- read.table("supplement_table1_expression_matrix.csv",header = T,row.names = 1,sep = "\t")
exprDataframe <- exprDataframe[which(rowSums(exprDataframe)>0),]
exprMatr <- as.matrix(exprDataframe)
dim(exprMatr)
gene_name <- read.table("TF_list_from_plantTFDB.txt",header = T,sep = "\t")
head(gene_name)
regulators <- Reduce(intersect,list(gene_name$TF_name,rownames(exprDataframe)))
head(regulators)
weightMat <- GENIE3(exprMatr,regulators = regulators,nCores=16)
linklist <- getLinkList(weightMat)
filter_linklist <- subset(linklist,linklist$weight>0.005)
write.table(filter_linklist,file = "GENIE3_default.txt",row.names = F,quote = F)






