###And then the R code shows how the graph is drawn,the "Zma_TF_list.txt" is from PlantTFDB database
###Seed germination bar plot
sub_network <- read.table("sub_network_without_annotation.txt",header = T,sep = "\t")
TF_list <- read.table("Zma_TF_list.txt",header = T,sep = "\t")
sub_network$Family_name <- TF_list$Family[match(sub_network$TF,TF_list$Gene_ID)]
table_seed_germination <- data.frame(table(sub_network$Family_name))
table_seed_germination <- table_seed_germination[order(table_seed_germination$Freq,decreasing = F),]
head(table_seed_germination)
table_seed_germination$Freq <- as.numeric(table_seed_germination$Freq)
table_seed_germination$Var1 <- factor(table_seed_germination$Var1,levels = table_seed_germination$Var1)
ggplot(table_seed_germination, aes(x = Freq, y = Var1)) +
  geom_bar(stat = "identity") +
  xlab("Family_Frequency") +
  ylab("Transcription Factor Family") +
  ggtitle("Seed Germination TF family Bar Graph")+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y.left = element_text(hjust = 1),
        )+
  scale_x_continuous(expand = c(0,0.5))

####four network venn
cluster <- read.table("CB_cluster_iGRN.txt",header = T,sep = "\t")
genie3 <- read.table("GENIE3_iGRN.txt",header = T,sep = "\t")
COE <- read.table("COE_iGRN.txt",header = T,sep = "\t")
IGRN <- read.table("IGRN_seperate.txt",header = T,sep = "\t")

combined_data <- list(CLUSTER=paste(cluster$TF,sep = "_",cluster$target),GENIE3=paste(genie3$TF,sep = "_",genie3$target),COE=paste(COE$TF,sep = "_",COE$target),IGRN=IGRN$Interaction)
p1=ggvenn(combined_data)
p1
####Bubble map of all transcription factor families in front of four networks
cluster <- read.table("CB_cluster_iGRN.txt",header = T,sep = "\t")
genie3 <- read.table("GENIE3_iGRN.txt",header = T,sep = "\t")
COE <- read.table("COE_iGRN.txt",header = T,sep = "\t")
IGRN <- read.table("IGRN_seperate.txt",header = T,sep = "\t")
TF_list <- read.table("Zma_TF_list.txt",header = T,sep = "\t")
###Find the family name corresponding to TF_ID in each network
cluster$family <- TF_list$Family[match(cluster$TF,TF_list$Gene_ID)]
genie3$family <- TF_list$Family[match(genie3$TF,TF_list$Gene_ID)]
COE$family <- TF_list$Family[match(COE$TF,TF_list$Gene_ID)]
IGRN$family <- TF_list$Family[match(IGRN$TF,TF_list$Gene_ID)]
####The all transcription factor families  in each network were extracted
table_cluster <- data.frame(table(cluster$family))
table_cluster_sub <- table_cluster[order(table_cluster$Freq,decreasing = T),]
table_genie3 <- data.frame(table(genie3$family))
table_genie3_sub <- table_genie3[order(table_genie3$Freq,decreasing = T),]
table_COE <- data.frame(table(COE$family))
table_COE_sub <- table_COE[order(table_COE$Freq,decreasing = T),]
table_IGRN <- data.frame(table(IGRN$family))
table_IGRN_sub <- table_IGRN[order(table_IGRN$Freq,decreasing = T),]
colnames(table_cluster_sub) <- c("Var1","cluster")
colnames(table_COE_sub) <- c("Var1","COE")
colnames(table_genie3_sub) <- c("Var1","genie3")
colnames(table_IGRN_sub) <- c("Var1","IGRN")
TF_family <- list(table_cluster_sub,table_COE_sub,table_genie3_sub,table_IGRN_sub)

# Merge the dataframes based on the 'Var' column
dataframe1 <- Reduce(function(x, y) merge(x, y, by = "Var1", all = TRUE), TF_family)
# Replace missing values with 0
dataframe1[is.na(dataframe1)] <- 0
###data_normalization
dataframe2 <- as.data.frame(scale(dataframe1[2:5]))
dataframe2$Family_name <- dataframe1$Var1
dataframe3 <- select(dataframe2,5,1,2,3,4)
colnames(dataframe1) <- c("Family_name","CLUSTER","COE","GENIE3","IGRN")
###draw the bubble plot
library(reshape2)
library(ggplot2)
data_melt<-melt (dataframe3)
names(data_melt) = c('Family_name', 'variable', 'value')

p2 <- ggplot(data_melt, aes(x = variable, y = Family_name, size = value, color=variable)) + geom_point()+ 
  theme(panel.background = element_blank(),panel.grid.major = element_line(colour = "gray"),
        panel.border = element_rect(colour="black",fill=NA))
p2
####Bubble maps of the top 10 transcription factor families in four networks
cluster <- read.table("CB_cluster_iGRN.txt",header = T,sep = "\t")
genie3 <- read.table("GENIE3_iGRN.txt",header = T,sep = "\t")
COE <- read.table("COE_iGRN.txt",header = T,sep = "\t")
IGRN <- read.table("IGRN_seperate.txt",header = T,sep = "\t")
TF_list <- read.table("Zma_TF_list.txt",header = T,sep = "\t")
###Find the family name corresponding to TF_ID in each network
cluster$family <- TF_list$Family[match(cluster$TF,TF_list$Gene_ID)]
genie3$family <- TF_list$Family[match(genie3$TF,TF_list$Gene_ID)]
COE$family <- TF_list$Family[match(COE$TF,TF_list$Gene_ID)]
IGRN$family <- TF_list$Family[match(IGRN$TF,TF_list$Gene_ID)]
####The top ten transcription factor families with the largest number in each network were extracted
table_cluster <- data.frame(table(cluster$family))
table_cluster_sub <- table_cluster[order(table_cluster$Freq,decreasing = T),][1:10,]
table_genie3 <- data.frame(table(genie3$family))
table_genie3_sub <- table_genie3[order(table_genie3$Freq,decreasing = T),][1:10,]
table_COE <- data.frame(table(COE$family))
table_COE_sub <- table_COE[order(table_COE$Freq,decreasing = T),][1:10,]
table_IGRN <- data.frame(table(IGRN$family))
table_IGRN_sub <- table_IGRN[order(table_IGRN$Freq,decreasing = T),][1:10,]
colnames(table_cluster_sub) <- c("Var1","cluster")
colnames(table_COE_sub) <- c("Var1","COE")
colnames(table_genie3_sub) <- c("Var1","genie3")
colnames(table_IGRN_sub) <- c("Var1","IGRN")
TF_family <- list(table_cluster_sub,table_COE_sub,table_genie3_sub,table_IGRN_sub)

# Merge the dataframes based on the 'Var' column
dataframe1 <- Reduce(function(x, y) merge(x, y, by = "Var1", all = TRUE), TF_family)
# Replace missing values with 0
dataframe1[is.na(dataframe1)] <- 0
###data_normalization
dataframe2 <- as.data.frame(scale(dataframe1[2:5]))
dataframe2$Family_name <- dataframe1$Var1
dataframe3 <- select(dataframe2,5,1,2,3,4)
colnames(dataframe3) <- c("Family_name","CLUSTER","COE","GENIE3","IGRN")
###draw the bubble plot
library(reshape2)
library(ggplot2)
data_melt<-melt (dataframe3)
names(data_melt) = c('Family_name', 'variable', 'value')

p3 <- ggplot(data_melt, aes(x = variable, y = Family_name, size = value, color=variable)) + geom_point()+ 
  theme(panel.background = element_blank(),panel.grid.major = element_line(colour = "gray"),
        panel.border = element_rect(colour="black",fill=NA))
p3


