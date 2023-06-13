###"IGRN_seperate.txt" is the result of "IGRN_split.R",Split the TF and target columns from the interaction column to facilitate subsequent operations
### "entrze.txt" contains the interaction between target gene in IGRN_network and its entrze ID.it can be acquied from maizeGDB database using ID convert tools
entrze <- read.table("entrze.txt",header = T,sep = "\t")
colnames(entrze) <- c("target","entrze_id")
IGRN_seperate <- read.table("IGRN_seperate.txt",header = T,sep = "\t")

order <- IGRN_seperate[order(IGRN_seperate$y,decreasing = T),]
net <- select(IGRN_seperate,8,9,6)
net$entrze_id <- merge(net,entrze,by = "target")

# Create a new column 'entrze_id' in 'net' and initialize it with NA
net$entrze_id <- NA
# Use match() function to find the indices of matching 'target' values in 'entrze'
matching_indices <- match(net$target, entrze$target)
# Find the indices where a match was found
valid_indices <- !is.na(matching_indices)
# Assign the corresponding 'entrze_id' values to the matching rows in 'net'
net$entrze_id[valid_indices] <- entrze$entrze_id[matching_indices[valid_indices]]
###supplement_table9_annotation_result.xlsx
annotation_information <- read.table("annotation_result.txt",header = T,sep="\t")
# Create a new column 'entrze_id' in 'net' and initialize it with NA
net$annotation_information <- NA
# Use match() function to find the indices of matching 'target' values in 'entrze'
matching_indices <- match(net$entrze_id, annotation_information$ID)
# Find the indices where a match was found
valid_indices <- !is.na(matching_indices)
# Assign the corresponding 'annotation_information' values to the matching rows in 'net'
net$annotation_information[valid_indices] <- annotation_information$GOTERM_BP_DIRECT[matching_indices[valid_indices]]
write.table(net,"TF_target_annotation.txt",col.names = T,row.names = F,sep = "\t")



IGRn <- read.table("IGRN_seperate.txt",header=T,sep = "\t")
annotation_information <- read.table("TF_target_annotation.txt",header = T,sep = "\t")
gene_seed_germination <- c("100285725","542310","103631821","103638756","100281175","100284718","103647149","103651483"
                           ,"103632346","100192497","100193516","100274557","100283778","100383462","778433",	
                           "100282389","100217313","100273980","541642","100191553","100281525")
gene <- annotation_information[annotation_information$entrze_id %in% gene_seed_germination, ]

TF_germination <- unique(gene$TF)
target_germination <- unique(gene$target)

germination_subnetwork <- gene[order(gene$y,decreasing = T),]
write.table(germination_subnetwork,"germination_subnetwork.txt",col.names = T,row.names = F,quote=F,sep = "\t")


sub_netwok <- read.table("germination_subnetwork.txt",header = T,sep = "\t")

sub_netwok <- select(sub_netwok,1:3)
write.table(sub_netwok,"sub_network_without_annotation.txt",col.names = T,quote = F,sep = "\t",row.names = F)
