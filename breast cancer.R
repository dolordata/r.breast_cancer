rm(list=ls())

library(e1071)
library(caret)
library(kernlab)
library(corrplot)

setwd("~/Desktop/CSC 529/case study 2")
bcancer <- read.table("breast-cancer-wisconsin.data", sep = ",", header = FALSE)


#data cleaning
#add variables names
colnames(bcancer)[1] <- "code.number"
colnames(bcancer)[2] <- "clump.thickness"
colnames(bcancer)[3] <- "cell.size"
colnames(bcancer)[4] <- "cell.shape"
colnames(bcancer)[5] <- "marginal.adhesion"
colnames(bcancer)[6] <- "epithelial.size"
colnames(bcancer)[7] <- "bare.nuclei"
colnames(bcancer)[8] <- "bland.chromatin"
colnames(bcancer)[9] <- "normal.nucleoli"
colnames(bcancer)[10] <- "mitoses"
colnames(bcancer)[11] <- "class"

#remove id column
bcancer$code.number <- NULL

#remove missing values
na.omit(bcancer)

#convert factor value to numeric for linear regression
bcancer$bare.nuclei <- as.numeric(bcancer$bare.nuclei)
str(bcancer)

#exploratory analysis

str(bcancer)
summary(bcancer)
prop.table(table(bcancer$class))

#correlation plot
M <- cor(bcancer)
corrplot(M, method="circle")


#fit linear regression model
fit <- glm(class~., data=bcancer)
summary(fit)

#change classification lable to factor
bcancer$class[bcancer$class=="2"] <- "benign"
bcancer$class[bcancer$class=="4"] <- "malignant"
bcancer$class <- factor(bcancer$class)
bcancer$bare.nuclei <- factor(bcancer$bare.nuclei)

bcancer$class <- factor(bcancer$class, levels=c("malignant","benign"))


#prepare training and test dataset
ind <- sample(2, nrow(bcancer), replace=TRUE, prob=c(0.8, 0.2))
train <- bcancer[ind==1,]
test <- bcancer[ind==2,]


#build the original svm model
svm_model <- svm(class ~ ., data=train, cross=10)
summary(svm_model)

pred <- predict(svm_model, test)
t <- table(pred, test$class)
confusionMatrix(t)


#svm linear using cv
svm_model.1 <- svm(class ~ ., kernel="linear", data=train, cross=10)
summary(svm_model.1)
pred.1 <- predict(svm_model.1, test)
t1 <- table(pred.1,test$class)
confusionMatrix(t1)

#svm polynomial using cv
svm_model.2 <- svm(class ~ ., kernel="polynomial", data=train, cross=10)
summary(svm_model.2)
pred.2 <- predict(svm_model.2,test)
t2 <- table(pred.2,test$class)
confusionMatrix(t2)

#svm radial using cv
svm_model.3 <- svm(class ~ ., kernel='radial', data=train, cross=10)
summary(svm_model.3)
pred.3 <- predict(svm_model.3,test)
t3 <- table(pred.3,test$class)
confusionMatrix(t3)

#svm sigmoid using cv
svm_model.4 <- svm(class ~ ., kernel="sigmoid", data=train, cross=10)
summary(svm_model.4)
pred.4 <- predict(svm_model.4,test)
t4 <- table(pred.4,test$class)
confusionMatrix(t4)

#tuned kernel svm polynomial using cv
svm_model.poly.tuned <- tune.svm(class ~ ., kernel = "polynomial", data = train, coef0 = (-1:4), degree = (1:4), tunecontrol=tune.control(sampling="cross", cross=10))
summary(svm_model.poly.tuned)
plot(svm_model.poly.tuned,xlab="degree", ylab="coef0")

svm_model.2.tuned <- svm(class ~ ., data=train, kernel = "polynomial", coef0 = 1, degree = 2, tunecontrol=tune.control(sampling="cross", cross=10))
summary(svm_model.2.tuned)
pred.2.tuned <- predict(svm_model.2.tuned, test)
t2.tuned <- table(pred.2.tuned, test$class)
confusionMatrix(t2.tuned)

#tuned radial using cv
svm_model.radial.tuned <- tune.svm(class ~ ., kernel = "radial", data = train, cost=(1:4), gamma=(1:4), tunecontrol=tune.control(sampling="cross", cross=10))
summary(svm_model.radial.tuned)
plot(svm_model.radial.tuned,xlab="gamma", ylab="cost")

svm_model.3.tuned <- svm(class ~ ., data=train, kernel = "radial", cost=2, gamma=1, tunecontrol=tune.control(sampling="cross", cross=10))
summary(svm_model.3.tuned)
pred.3.tuned <- predict(svm_model.3.tuned, test)
t3.tuned <- table(pred.3.tuned, test$class)
confusionMatrix(t3.tuned)


#tuned linear using cv
svm_model.linear.tuned <- tune.svm(class ~ ., kernel = "linear", data = train, cost=(1:4), tunecontrol=tune.control(sampling="cross", cross=10))
summary(svm_model.linear.tuned)
plot(svm_model.linear.tuned)

svm_model.1.tuned <- svm(class ~ ., data=train, kernel = "linear", cost=1, tunecontrol=tune.control(sampling="cross", cross=10))
summary(svm_model.1.tuned)
pred.1.tuned <- predict(svm_model.1.tuned, test)
t1.tuned <- table(pred.1.tuned, test$class)
confusionMatrix(t1.tuned)


