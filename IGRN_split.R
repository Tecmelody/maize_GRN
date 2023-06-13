###导入最终网络
IGRN <- read.table("iGRN_network.txt",header = T,sep = "\t")
AC <- IGRN[grep("^AC",IGRN$Interaction),]
AC$TF <- substr(AC$Interaction,1,16)
AC$target <- substring(AC$Interaction,18)
GRMZM <- IGRN[grep("^GRMZM",IGRN$Interaction),]
GRMZM$TF <- substr(GRMZM$Interaction,1,13)
GRMZM$target <- substring(GRMZM$Interaction,15)
EF <- IGRN[grep("^EF",IGRN$Interaction),]
EF$TF <- substr(EF$Interaction,1,16)
EF$target <- substring(EF$Interaction,18)
AF <- IGRN[grep("^AF",IGRN$Interaction),]
AF$TF <- substr(AF$Interaction,1,16)
AF$target <- substring(AF$Interaction,18)


IGRN_seperate <- rbind(AC,AF,EF,GRMZM)
write.table(IGRN_seperate,"IGRN_seperate.txt",row.names = F,col.names = T,quote = F,sep = "\t")
target <- as.data.frame(unique(IGRN_seperate$target))
colnames(target) <- "target"
write.table(target,"IGRN_target.txt",col.names = T,row.names = F,quote = F,sep = "\t")
