##############################################################################################################################
################################################ Predictive models ###########################################################
##############################################################################################################################


##################################  

#1. Merge the minimum set of growth parameters for TS cases, GHD cases and referents (sexes and ages combined)

#2. Create age-specific predictive models 

#3. Retain the risk thresholds that offered a pre-defined specificity of 98% 

##################################  




# Library
# -------
library(caret)
library(dplyr)
library(pROC)
library(gd)

path="our path" 


# Data set
# --------
# traindata includes : 
## id : subject
## age : age of measurement (days)
## sex : sex of subject
## x : TS cases (=1), GHD cases (=2) or referents (=0)
## x2 : TS cases (="MTS"), GHD cases (="MGHD") or referents (="NM")
## A
## B
## c
## D
## E 

load(file=paste0(path,"/traindata.rda"))




## selection of last set of growth parameters by age range
selection <- function(training_data, age_min, age_max) {
  traindata_age <- subset(training_data, age>=age_min & age <age_max)
  traindata_age <- traindata_age %>% group_by(id) %>% top_n(1,age)
  traindata_age[['rowIndex']]<- as.numeric(rownames(traindata_age))
  return(traindata_age)
}
traindata_1_12 <- selection(traindata,365,4381)
traindata_1_2  <- selection(traindata,365,730)
traindata_2_3  <- selection(traindata,730,1095)
traindata_3_5  <- selection(traindata,1095,1825)
traindata_5_8  <- selection(traindata,1825,2920)
traindata_8_12 <- selection(traindata,2920,4381)


# Predictive models
# -----------------
##### Multinomial logistic regression 
ML <- function(training_data, meth, metr, tunegrid, trcontrol) {
  training_data$x2 <- factor(ifelse(training_data$ts==1, "MTS", ifelse(training_data$ghd==1, "MGHD", ifelse(training_data$referents==1, "NM", NA))))
  training_data$x2 <-relevel(training_data$x2, ref = "NM")
  
  if (max(training_data[["age"]]) <2920) {
    fit <- train(x2~A+B+C+D+sex,   data= training_data, method=meth, family = "multinomial",   metric=metr,  tuneGrid=tunegrid, trControl=trcontrol)
  }
  if (max(training_data[["age"]])>=2920) {
    fit <- train(x2~A+B+C+D+E+sex, data= training_data, method=meth, family = "multinomial",  metric=metr,  tuneGrid=tunegrid, trControl=trcontrol)
  }
  return(fit)
}
fit.rlm.mf.1_12 <- ML(traindata_1_12, "multinom",  "logLoss", NULL,  trainControl(method="repeatedcv",    number=5,  repeats = 1,   allowParallel=TRUE, summaryFunction = multiClassSummary, classProbs = TRUE, savePredictions="final"))
fit.rlm.mf.1_2  <- ML(traindata_1_2,  "multinom",  "logLoss", NULL,  trainControl(method="repeatedcv",    number=5,  repeats = 1,   allowParallel=TRUE, summaryFunction = multiClassSummary, classProbs = TRUE, savePredictions="final")) 
fit.rlm.mf.2_3  <- ML(traindata_2_3,  "multinom",  "logLoss", NULL,  trainControl(method="repeatedcv",    number=5,  repeats = 1,   allowParallel=TRUE, summaryFunction = multiClassSummary, classProbs = TRUE, savePredictions="final")) 
fit.rlm.mf.3_5  <- ML(traindata_3_5,  "multinom",  "logLoss", NULL,  trainControl(method="repeatedcv",    number=5,  repeats = 1,   allowParallel=TRUE, summaryFunction = multiClassSummary, classProbs = TRUE, savePredictions="final"))
fit.rlm.mf.5_8  <- ML(traindata_5_8,  "multinom",  "logLoss", NULL,  trainControl(method="repeatedcv",    number=5,  repeats = 1,   allowParallel=TRUE, summaryFunction = multiClassSummary, classProbs = TRUE, savePredictions="final")) 
fit.rlm.mf.8_12 <- ML(traindata_8_12, "multinom",  "logLoss", NULL,  trainControl(method="repeatedcv",    number=5,  repeats = 1,   allowParallel=TRUE, summaryFunction = multiClassSummary, classProbs = TRUE, savePredictions="final")) 


##### Performances (varying thresholds including thresholds=0.5)
thresholder_rlm <- function(algorithme_training, final, mean, threshold, k) {
  statistics <- c("Sensitivity", "Specificity", "AUC")
  
  
  # Expand the predicted values with the candidate values of the threshold
  expand_preds <- function(df, th, excl = NULL) {
    th <- unique(th)
    nth <- length(th)
    ndf <- nrow(df)
    if (!is.null(excl))
      df <- df[, !(names(df) %in% excl), drop = FALSE]
    df <- df[rep(1:nrow(df), times = nth),]
    df$prob_threshold <- rep(th, each = ndf)
    df
  }
  
  # Modify the threshold for varying the specificity
  recode <- function(dat) {
    lvl <- levels(dat$obs)
    
    dat$MTSPHEs  <-rowSums(data.frame(dat[, lvl[2]],dat[, lvl[3]]), na.rm=FALSE)
    dat$MTSPHEm  <- apply(data.frame(dat[, lvl[2]],dat[, lvl[3]]), 1, which.max)
    
    # en 3 classes
    dat$obs2 <- dat$obs;  dat$obs2  <- factor(dat$obs2, levels=lvl)
    dat$pred2 <- ifelse(dat$MTSPHEs >= dat$prob_threshold,colnames(data.frame(dat[, c("MTS","MGHD")]))[dat$MTSPHEm],"NM");dat$pred2 <- factor(dat$pred2, levels=lvl)
    
    # en 2 classes
    dat$obs1  <- ifelse(dat$obs=="NM", "NM", ifelse(dat$obs=="MTS" | dat$obs=="MGHD", "M", NA));  dat$obs1  <- factor(dat$obs1, levels=c("NM","M"))
    dat$pred1  <- ifelse(dat$MTSPHEs >= dat$prob_threshold, "M", "NM");  dat$pred1  <- factor(dat$pred1, levels=c("NM","M"))
    
    dat
  }
  
  # Compute statistics per threshold and tuning parameters
  stats<- function(dat) {
    tab <- caret::confusionMatrix(dat$pred1, dat$obs1, positive = "M")
    
    # Performances
    sp=tab$table[1]/(tab$table[1] + tab$table[2])
    se=tab$table[4]/(tab$table[3] + tab$table[4])
    IClowspe    <- round(binom.test(tab$table[1],tab$table[1] + tab$table[2], alternative = c("two.sided"), conf.level = 0.95)$conf.int[1]*100,2)
    IChighspe   <- round(binom.test(tab$table[1],tab$table[1] + tab$table[2], alternative = c("two.sided"), conf.level = 0.95)$conf.int[2]*100,2)
    IClowsens   <- round(binom.test(tab$table[4],tab$table[3] + tab$table[4], alternative = c("two.sided"), conf.level = 0.95)$conf.int[1]*100,2)
    IChighsens  <- round(binom.test(tab$table[4],tab$table[3] + tab$table[4], alternative = c("two.sided"), conf.level = 0.95)$conf.int[2]*100,2)
    se_CI <- paste0(IClowsens, "-",IChighsens)
    sp_CI <- paste0(IClowspe, "-",IChighspe)
    
    # AUC and CI
    AUC    <- round(auc(dat$obs1, dat$MTSPHEs),3)
    AUC_CI <- paste0(round(ci.auc(dat$obs1, dat$MTSPHEs)[1],3),"-",round(ci.auc(dat$obs1, dat$MTSPHEs)[3],3))
    
    res <- c(Sensitivity=as.numeric(se), se_CI=se_CI, Specificity=as.numeric(sp), sp_CI=sp_CI, AUC=as.numeric(AUC), AUC_CI=AUC_CI)
    return(res)
  }
  
  # Summarize over resamples (by mean)
  summ_stats <- function(x, cols, r) {
    colSd <- function (x, na.rm=FALSE) {apply(X=x, MARGIN=2, FUN=sd, na.rm=na.rm)}
    
    na_cols <- apply(x, 2, function(x) any(is.na(x)))
    na_col_names <- colnames(x)[na_cols]
    relevant_col_names <- intersect(na_col_names, cols)
    if (length(relevant_col_names) > 0) warning("The following columns have missing values (NA), which have been ", "removed: '", paste0(relevant_col_names, collapse = "', '"), "'.\n")
    if (r==1) warning("The calculation of CI is impossible")
    
    x[cols]<-lapply(x[cols], as.numeric)
    
    # performance
    resmeans <- colMeans(x[, cols, drop = FALSE], na.rm = TRUE)
    ressd    <- colSd(x[,cols, drop = FALSE], na.rm=TRUE)
    
    # CI
    IClow   <- resmeans  - 1.96*ressd/sqrt(r)
    IChigh  <- resmeans  + 1.96*ressd/sqrt(r)
    
    
    res <- c(se=resmeans[[1]],
             se_CI=paste0(round(IClow[1]*100,2),"-",round(IChigh[1]*100,2)), 
             sp=resmeans[[2]],
             sp_CI=paste0(round(IClow[2]*100,2),"-",round(IChigh[2]*100,2)),
             AUC=resmeans[[3]],
             AUC_CI=paste0(round(IClow[3],3),"-",round(IChigh[3],3)))
    
    
    return(res)
  }
  
  
  pred_dat <- expand_preds(df=if(final==F) algorithme_training$pred 
                           else  cbind.data.frame(decay=algorithme_training$bestTune, 
                                                  obs=algorithme_training$trainingData$.outcome, 
                                                  pred=predict(algorithme_training, algorithme_training$trainingData, type="raw"),
                                                  predict(algorithme_training,      algorithme_training$trainingData, type="prob")),
                           th=threshold)
  
  
  ## Based on the threshold, recode the predicted classes
  if (mean==F) {pred_dat <- ddply(pred_dat, .variables = c("prob_threshold"),             recode) }
  if (mean==T) {pred_dat <- ddply(pred_dat, .variables = c("prob_threshold", "Resample"), recode) }
  
  ## Compute statistics per threshold
  if (mean==F) {pred_stats <- ddply(pred_dat, .variables = c("prob_threshold"),             stats)}
  if (mean==T) {pred_stats <- ddply(pred_dat, .variables = c("prob_threshold", "Resample"), stats)}
  
  ## Summarize statistics / resample
  if (mean==T) {pred_stats <- ddply(pred_stats, .variables = c("prob_threshold"), summ_stats,statistics, r=k)}
  
  return(pred_stats)
}
fit.rlm.mf.th_1_12 <-thresholder_rlm(algorithme_training=fit.rlm.mf.1_12,  final=T, mean=F, threshold = seq(0, 1, by = 0.001), k=5)
fit.rlm.mf.th_1_2  <-thresholder_rlm(algorithme_training=fit.rlm.mf.1_2,   final=T, mean=F, threshold = seq(0, 1, by = 0.001), k=5)
fit.rlm.mf.th_2_3  <-thresholder_rlm(algorithme_training=fit.rlm.mf.2_3,   final=T, mean=F, threshold = seq(0, 1, by = 0.001), k=5)
fit.rlm.mf.th_3_5  <-thresholder_rlm(algorithme_training=fit.rlm.mf.3_5,   final=T, mean=F, threshold = seq(0, 1, by = 0.001), k=5)
fit.rlm.mf.th_5_8  <-thresholder_rlm(algorithme_training=fit.rlm.mf.5_8,   final=T, mean=F, threshold = seq(0, 1, by = 0.001), k=5)
fit.rlm.mf.th_8_12 <-thresholder_rlm(algorithme_training=fit.rlm.mf.8_12,  final=T, mean=F, threshold = seq(0, 1, by = 0.001), k=5)


##### Thresholds
thresholds <- function (algorithme_th,thres_spe_min) {
s <-  cbind(threshold=min(subset(cbind.data.frame(th=algorithme_th$prob_threshold, ind=findInterval(algorithme_th$Specificity, thres_spe_min)), ind==1)$th, na.rm=T),  
            specificity=algorithme_th$Specificity[algorithme_th$prob_threshold==min(subset(cbind.data.frame(th=algorithme_th$prob_threshold, ind=findInterval(algorithme_th$Specificity, thres_spe_min)), ind==1)$th, na.rm=T)], 
            sensitivity=algorithme_th$Sensitivity[algorithme_th$prob_threshold==min(subset(cbind.data.frame(th=algorithme_th$prob_threshold, ind=findInterval(algorithme_th$Specificity, thres_spe_min)), ind==1)$th, na.rm=T)])

return(s)
}
s98_rlm.mf_1_12 <-thresholds(algorithme_th=fit.rlm.mf.th_1_12, thres_spe_min=0.980)     
s98_rlm.mf_1_2  <-thresholds(algorithme_th=fit.rlm.mf.th_1_2,  thres_spe_min=0.980)     
s98_rlm.mf_2_3  <-thresholds(algorithme_th=fit.rlm.mf.th_2_3,  thres_spe_min=0.980)     
s98_rlm.mf_3_5  <-thresholds(algorithme_th=fit.rlm.mf.th_3_5,  thres_spe_min=0.980)     
s98_rlm.mf_5_8  <-thresholds(algorithme_th=fit.rlm.mf.th_5_8,  thres_spe_min=0.980)    
s98_rlm.mf_8_12 <-thresholds(algorithme_th=fit.rlm.mf.th_8_12, thres_spe_min=0.980) 


##### Prediction from newdata
valdata=cbind.data.frame(id=c(1,2,3,4,5), 
                         sex=c(1,2,2,1,1),
                         A=c(3.88, 3.90, 3.99, 3.87, 3.77), 
                         B=c(-3.51,-3.56,-2.61, -3.49, -4.30),
                         C=c(2.76, 2.76, 3.71, 2.86,  2.78),
                         D=c(-4.30,-4.37, -5.66, -3.67, -4.89),
                         E=c(-19.0,-19.4,-19.1,-20.1,-18.1))

prediction <- function (algorithme_training, newdata, thres){
  pred <- cbind.data.frame(pred =predict(algorithme_training,  newdata,  type="prob"),
                           predc=rowSums(predict(algorithme_training,  newdata,  type="prob")[,c("MTS", "MGHD")])>=thres)
  
return(pred)
}
prediction(algorithme_training=fit.rlm.mf.1_12,newdata=valdata, thres=s98_rlm.mf_1_12[[1]])
prediction(algorithme_training=fit.rlm.mf.1_2, newdata=valdata, thres=s98_rlm.mf_1_2[[1]])
prediction(algorithme_training=fit.rlm.mf.2_3, newdata=valdata, thres=s98_rlm.mf_2_3[[1]])
prediction(algorithme_training=fit.rlm.mf.3_5, newdata=valdata, thres=s98_rlm.mf_3_5[[1]])
prediction(algorithme_training=fit.rlm.mf.5_8, newdata=valdata, thres=s98_rlm.mf_5_8[[1]])
prediction(algorithme_training=fit.rlm.mf.8_12,newdata=valdata, thres=s98_rlm.mf_8_12[[1]])



