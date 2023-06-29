data_dir <- "/home/melody/pwms_all_motifs/"

pwm_file_list <- as.list(paste0(data_dir,list.files(data_dir)))

pwm_file_list <- lapply(pwm_file_list, function(x){
  t <- read.table(x,header = T,row.names = 1)
  apply(t,2,as.character)
})
c <- c("letter-probability matrix:alength=4 E=0","","","")
for (i in 1:length(pwm_file_list)) {
  name <- gsub(pattern = "*.txt","",list.files(data_dir))[i]
  colnames(pwm_file_list[[i]]) <- c("MOTIF",name,"","")
  pwm_file_list[[i]] <- rbind(colnames(pwm_file_list[[i]]),c("letter-probability matrix:alength=4 E=0","","",""),pwm_file_list[[i]])
  colnames(pwm_file_list[[i]]) <- NULL
}

merged_pwm <- do.call(rbind,pwm_file_list)

write.table(merged_pwm,file = paste0(data_dir,"merged_pwm.txt"),sep = "\t",
            row.names = F,col.names = F,quote = F)
