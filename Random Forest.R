###################
## Random Forest ##
###################

library(tidyverse)
library(caret)
library(randomForest)
library(data.table)

setwd("C:/Liang Wang/gene level data")
data <- fread("Qianxia_raw_counts_by_Genes_FINAL.csv")
Final <- read_excel("temporal samples.xlsx", sheet = "Final")

#####
# Time 1
#####

## Set up data ----

# Select Time 1 data
DESeq.data.1 <- data[,c(1,7:61)]

# Set up metadata/design matrix
metaData.1 <- data.frame(Final[1:56,c(2,4)])
metaData.1 <- metaData.1[-41,]  # A46 at row 41
colnames(metaData.1)[2] <- "Pheno"

# Create DESeq Data format
dds.1 <- DESeqDataSetFromMatrix(countData=DESeq.data.1, 
                                colData=metaData.1, 
                                design=~Pheno, tidy = TRUE)

# Run DESeq
dds.1 <- DESeq(dds.1)

# Normalize data to write out
normal_RC.1 <- counts(dds.1, normalized = T)

t_normal_RC.1 <- t(normal_RC.1)
colnames(t_normal_RC.1) <- DESeq.data.1$Geneid


## Random Forest ----

# Coefficient of variation
cv <- apply(log(normal_RC.1+1), 1, function(x) sd(x, na.rm=T) / mean(x, na.rm=T) )

hist(cv, nclass=50)
head(sort(cv, decreasing=T), 100)

sum(cv>0.75)

# Filter genes with highest CV
top <- as.data.frame(t_normal_RC.1[,order(cv)[1:5000]])
top <- cbind(as.factor(metaData.1$Pheno), top)
colnames(top)[1] <- "Pheno"

names(top) <- make.names(names(top))


# Cross-validation
n <- 100
acc <- rep(NA,n)
spec <- rep(NA,n)
sens <- rep(NA,n)

for(i in 1:n){
  ind <- sample(1:nrow(top), size=6)
  top.train <- top[-ind,]
  top.test <- top[ind,]
    
  rf.model <- randomForest(Pheno ~ .,
                           data=top.train,
                           ntree=100,
                           mtry=2,
                           importance=TRUE)
  
  pred <- predict(rf.model, top.test)
  
  conf.mat <- table(pred, top.test$Pheno)
  acc[i] <- mean(pred == top.test$Pheno)
  spec[i] <- specificity(conf.mat)
  sens[i] <- sensitivity(conf.mat)
}

mean(acc)
mean(spec, na.rm=T)
mean(sens)

varImpPlot(rf.model)





# Train model
# tc <- trainControl(method='repeatedcv',
#                    number=5,
#                    repeats=3,
#                    classProbs=T)
# 
# model <- train(Pheno ~ .,
#                data = top,
#                method = 'ranger',
#                importance = 'impurity',
#                trControl = tc,
#                tuneLength = 6)
# beepr::beep(2)
# 
# model
# 
# #
# pred <- predict(model, top)
# prob <- predict(model, top, type="prob")
# prob$class <- top$Pheno
# prob$pred <- pred
# 
# confusionMatrix(top$Pheno, pred)
# 
# # Variable Importance
# imp <- varImp(model)
# plot(imp, top=20)
