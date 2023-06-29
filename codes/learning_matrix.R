###########The codes provide an example build the train matrix(lines 3 to lines 25) and test matrix(lines 30 to lines51).
######Attention :the example codes to build the train matrix using four features but codes to build test matrix only using three features.
######You need to correct the code according to the number of features you have 
library(dplyr)
gold <- read.table("Golden_standard_iGRN.txt",header = T,sep = "\t")
length(unique(gold$target))
genie3net <- read.table("GENIE3_iGRN.txt",header = T,sep = "\t")
COEnet <- read.table("COE_iGRN.txt",header = T,sep = "\t")
pwmnet <- read.table("CB_motif_iGRN.txt",header = T,sep = "\t")
CLUSTERnet <- read.table("CB_cluster_iGRN.txt",header = T,sep = "\t")
# create example data frames
TF <- as.data.frame(unique(gold$TF))
colnames(TF) <- c("TF_name")
Target <- as.data.frame(unique(gold$target))
colnames(Target) <- c("Target_name")
# use expand.grid to get all possible combinations of TF and Target
interactions <- expand.grid(TF_name = TF$TF_name, Target_name = Target$Target_name)
interactions$COE <- ifelse(paste0(interactions$TF_name, "\t", interactions$Target_name) %in% paste0(COEnet$TF, "\t", COEnet$target), COEnet$score[match(paste0(interactions$TF_name, "\t", interactions$Target_name), paste0(COEnet$TF, "\t", COEnet$target))], 0)
interactions$pwm <- ifelse(paste0(interactions$TF_name, "\t", interactions$Target_name) %in% paste0(pwmnet$TF, "\t", pwmnet$target), pwmnet$score[match(paste0(interactions$TF_name, "\t", interactions$Target_name), paste0(pwmnet$TF, "\t", pwmnet$target))], 0)
interactions$CLUSTER <- ifelse(paste0(interactions$TF_name, "\t", interactions$Target_name) %in% paste0(CLUSTERnet$TF, "\t", CLUSTERnet$target), CLUSTERnet$score[match(paste0(interactions$TF_name, "\t", interactions$Target_name), paste0(CLUSTERnet$TF, "\t", CLUSTERnet$target))], 0)
interactions$genie3 <- ifelse(paste0(interactions$TF_name, "\t", interactions$Target_name) %in% paste0(genie3net$TF, "\t", genie3net$target), genie3net$score[match(paste0(interactions$TF_name, "\t", interactions$Target_name), paste0(genie3net$TF, "\t", genie3net$target))], 0)
interactions$class <- ifelse(paste0(interactions$TF_name, "\t", interactions$Target_name) %in% paste0(gold$TF, "\t", gold$target), 1, 0)
interactions$interaction <- paste(interactions$TF_name,interactions$Target_name,sep = "_")
learning_matrix <- dplyr::select(interactions,8,3,4,5,6,7)
colnames(learning_matrix) <- c("Interaction","COE","PWM","CLUSTER","GENIE3","Class")
write.table(learning_matrix,"interactions.txt",sep = "\t",quote = F,row.names = F)




####Construct a feature matrix for prediction
genie3net <- read.table("GENIE3_iGRN.txt",header = T,sep = "\t")
COEnet <- read.table("COE_iGRN.txt",header = T,sep = "\t")
pwmnet <- read.table("CB_motif_iGRN.txt",header = T,sep = "\t")
CLUSTERnet <- read.table("CB_cluster_iGRN.txt",header = T,sep = "\t")
#Extract TF and target from all networks and select all possible relationships,expand.grid()
TF <- c(CLUSTERnet$TF,COEnet$TF,genie3net$TF,pwmnet$TF)
TF <- as.data.frame(unique(TF))
Target <- c(CLUSTERnet$target,COEnet$target,genie3net$target,pwmnet$target)
Target <- as.data.frame(unique(Target))
# generate all possible combinations
interactions <- expand.grid(TF_name = TF$`unique(TF)`, Target_name = Target$`unique(Target)`)

interactions$COE <- ifelse(paste0(interactions$TF_name, "\t", interactions$Target_name) %in% paste0(COEnet$TF, "\t", COEnet$target), COEnet$score[match(paste0(interactions$TF_name, "\t", interactions$Target_name), paste0(COEnet$TF, "\t", COEnet$target))], 0)
#interactions$pwm <- ifelse(paste0(interactions$TF_name, "\t", interactions$Target_name) %in% paste0(pwmnet$TF, "\t", pwmnet$target), pwmnet$score[match(paste0(interactions$TF_name, "\t", interactions$Target_name), paste0(pwmnet$TF, "\t", pwmnet$target))], 0)
interactions$CLUSTER <- ifelse(paste0(interactions$TF_name, "\t", interactions$Target_name) %in% paste0(CLUSTERnet$TF, "\t", CLUSTERnet$target), CLUSTERnet$score[match(paste0(interactions$TF_name, "\t", interactions$Target_name), paste0(CLUSTERnet$TF, "\t", CLUSTERnet$target))], 0)
interactions$genie3 <- ifelse(paste0(interactions$TF_name, "\t", interactions$Target_name) %in% paste0(genie3net$TF, "\t", genie3net$target), genie3net$score[match(paste0(interactions$TF_name, "\t", interactions$Target_name), paste0(genie3net$TF, "\t", genie3net$target))], 0)
interactions$interaction <- paste(interactions$TF_name,interactions$Target_name,sep = "_")
test_matrix <- dplyr::select(interactions,6,3,4,5)
colnames(test_matrix) <- c("Interaction","COE","CLUSTER","GENIE3")
head(test_matrix)
write.table(test_matrix,file = "test_matrix.txt",sep = "\t",quote = F,row.names = F,col.names = T)


