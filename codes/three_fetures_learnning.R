###three feature evaluate
library(data.table)
library(dplyr)
library(randomForest)
library(caret)
library(ROCR)
library(gbm)

set.seed(1)
#read in data
full <- fread("./interactions_three_features.txt", header=T)

#normalize
sc <- names(full)[2:4]
full <- copy(full)[ , (sc) := as.data.table(scale(.SD)), .SDcols = sc]

splitdf <- function(dataframe, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  index <- 1:nrow(dataframe)
  trainindex <- sample(index, trunc(length(index)/5))
  trainset <- dataframe[-trainindex, ]
  testset <- dataframe[trainindex, ]
  list(trainset=trainset,testset=testset)
}

list_all <- splitdf(full)
list_pos <- splitdf(filter(full, Class==1))

train_all <- list_all$trainset
neg <- sample_n(filter(train_all, Class==0), 3*nrow(list_pos$trainset))
trainset <- as.data.table(rbind(list_pos$trainset,neg))

test_all <- list_all$testset
neg <- sample_n(filter(test_all, Class==0), 3*nrow(list_pos$testset))
testset <- as.data.table(rbind(list_pos$testset,neg))


shuffled <- trainset[sample(nrow(trainset)),]
cross <- 10
#CV testing
accs <- rep(0,cross)
prec <- rep(0,cross)
rec <- rep(0,cross)


for (i in 1:cross) {
  indices <- (((i-1) * round((1/cross)*nrow(shuffled))) + 1):((i*round((1/cross) * nrow(shuffled))))
  train <- shuffled[-indices,]
  test <- shuffled[indices,]
  
  tree_gbm <- gbm(Class ~ COE+CLUSTER+GENIE3, data=train, distribution="bernoulli", n.trees=1000, shrinkage=0.01, interaction.depth=3)
  pred <- predict(tree_gbm, test, n.trees=1000,type="response")
  
  pred$class <- ifelse(pred > 0.4, 1, 0)
  conf <- table(test$Class,pred$class)
  
  accs[i] <- sum(diag(conf))/sum(conf)
  prec[i] <- conf[2,2]/sum(conf[,2])
  rec[i] <- conf[2,2]/sum(conf[2,])
}
F1 <- 2*mean(prec)*mean(rec)/(mean(prec)+mean(rec))
mean(accs)
mean(prec)
mean(rec)
F1 
tree_gbm <- gbm(Class ~ COE+CLUSTER+GENIE3, data=trainset, distribution="bernoulli", n.trees=1000, shrinkage=0.01, interaction.depth=3)
pred <- predict(tree_gbm, testset, n.trees=1000,type="response")

pred$class <- ifelse(pred > 0.25, 1, 0)
conf <- table(testset$Class,pred$class)

probs <- predict(tree_gbm, testset, type = 'response', n.trees=1000)
pred <- prediction(probs, testset$Class)
perf <- performance(pred, "tpr", "fpr")
plot(perf)
perf <- performance(pred, "auc")
perf@y.values[[1]]

accs <- sum(diag(conf))/sum(conf)
prec <- conf[2,2]/sum(conf[,2])
rec <- conf[2,2]/sum(conf[2,])
F1 <- 2*prec*rec/(prec+rec)
accs
prec
rec
F1


####The next section shows how the final network is generated(the test_matrix.txt are the result of learning_matrix.R)
predictor <- fread("test_matrix.txt", header=T)
sc <- names(predictor)[2:5]
predictor <- copy(predictor)[ , (sc) := as.data.table(scale(.SD)), .SDcols = sc]

tree_gbm <- gbm(Class ~ COE+PWM+CLUSTER+GENIE3, data=trainset, distribution="bernoulli", n.trees=1000, shrinkage=0.01, interaction.depth=3)
pred <- predict(tree_gbm, predictor, n.trees=1000,type="response")

predictor <- as.data.frame(predictor)
merged <- merge(predictor, pred, by="row.names")
merged_f <- filter(merged, y>0.25)
merged_f <- mutate(merged_f,rank=rank(-y, ties.method="random"))
merged_f <- arrange(merged_f, rank)
write.table(merged_f, file="IGRN_netowk.txt", quote = F, sep = "\t")
