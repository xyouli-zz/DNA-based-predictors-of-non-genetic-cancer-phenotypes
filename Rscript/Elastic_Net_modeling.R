library(sampling)
library(caret)
library(glmnet)
library(ROCR)
library(doMC)
registerDoMC(cores = 8)
set.seed(1992)
source("~/script/helper.R")

# load data, BRCA_signatue for example
# predictor: CN_score
# response: signature_score/rppa/mutation/clinical variables
# for balanced stratification: balancing_variables
load('~/data/BRCA_signature_elastic_net_data.rda')

# choose a phenotype to model
# for example, RB-LOH signature
pheno <- "UNC_RB_LOH_Median_Breast.Cancer.Res.2008_PMID.18782450"
score <- unlist(signature_score[pheno,])
score_bi <- ifelse(score >= quantile(score,0.67),1,0)

# split sample into 70% training set and 30% test set
strata <- score_bi
pik <- rep(7/10,times=length(strata))
balanced_split <- balancedstratification(balancing_variables,strata,pik,comment = F)
train <- which(balanced_split==1)
test <- which(balanced_split==0)
trainX <- CN_score[train,]
testX <- CN_score[test,]
trainY <- score_bi[train]
testY <- score_bi[test]

glmnet_obj <- caret_wrap(trainX,trainY,testX,testY,bi = T)

# look at model performance: receiving operating curve and precision recall curve
pred_train <- predict(glmnet_obj,newdata = trainX,type = 'prob')
pred_test <- predict(glmnet_obj,newdata = testX,type = 'prob')
pred_train <- prediction(pred_train$pos,labels = trainY)
pred_test <- prediction(pred_test$pos,labels = testY)
perf_train <- performance(pred_train,measure = 'tpr',x.measure = 'fpr')
perf_test <- performance(pred_test,measure = 'tpr',x.measure = 'fpr')
auc_train <- signif(performance(pred_train,measure = 'auc')@y.values[[1]][1],2)
auc_test <- signif(performance(pred_test,measure = 'auc')@y.values[[1]][1],2)

plot_ROC(perf_train,perf_test,auc_train,auc_test,pheno)

# look at feature landscape
beta <- as.matrix(coef(glmnet_obj$finalModel,glmnet_obj$bestTune$lambda))
beta <- beta[-1,1]
plot_seg_ss(beta,pheno)






