library(readxl)
library(caret)
library(pROC)
library(dplyr)
library(MLmetrics)
library(boot)
library(gmodels)

Test_slide <- read.csv(paste("~/documents/GBM/Results/1029/S1is_gcimp/out/Test_slide.csv", sep=''))
Test_tile <- read.csv(paste("~/documents/GBM/Results/1029/S1is_gcimp/out/Test_tile.csv", sep=''))
pos = "is_gcimp"

# per patient level
answers <- factor(Test_slide$True_label)
results <- factor(Test_slide$Prediction)
# statistical metrics
CMP = confusionMatrix(data=results, reference=answers, positive = pos)
# ROC
roc =  roc(answers, Test_slide$POS_score, levels=c('negative', pos))
rocdf = t(data.frame(ci.auc(roc)))
colnames(rocdf) = c('ROC.95.CI_lower', 'ROC', 'ROC.95.CI_upper')
# PRC
SprcR = PRAUC(Test_slide$POS_score, factor(Test_slide$True_label))
Sprls = list()
for (j in 1:100){
  sampleddf = Test_slide[sample(nrow(Test_slide), round(nrow(Test_slide)*0.8)),]
  Sprc = PRAUC(sampleddf$POS_score, factor(sampleddf$True_label))
  Sprls[j] = Sprc
}
Sprcci = ci(as.numeric(Sprls))
Sprcdf = data.frame('PRC.95.CI_lower' = Sprcci[2], 'PRC' = SprcR, 'PRC.95.CI_upper' = Sprcci[3])
# Combine and add prefix
soverall = cbind(rocdf, Sprcdf, data.frame(t(CMP$overall)), data.frame(t(CMP$byClass)))
colnames(soverall) = paste('Patient', colnames(soverall), sep='_')