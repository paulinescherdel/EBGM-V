##############################################################################################################################
################################################ Jenss-bayley modeling #######################################################
##############################################################################################################################


##################################  

##Step A

#1. Create a first data set, which includes the first two measurements of a given subject (case or referent), and merge it with the modelling sample.

#2. Run the Jenss–Bayley model on this data set (according to age).

#3. Store the minimum set of growth parameters obtained for this subject in a results data set.


##Step B

#1. Create a second data set that includes the first three measurements of the same subject (case or referent), and merge it with the modelling sample.

#2. Run the Jenss–Bayley model on this second data set (according to age).

#3. Merge the new minimum set of growth parameters for this subject with the previous results data set.


##Next steps

# Repeat these steps each new measurement for this subject, and then apply the full procedure to all subjects.

##################################  



# Library
# -------
library(saemix)


path="our path" 

# Data set
# --------
# tab includes : 
## id : subject
## age : age of measurement (days)
## sex : sex of subject (girl=2, boy=1)
## height : height (cm)
## x : cases or referents (=1) versus modeling sample (=0)

load(file=paste0(path,"/sample_G.rda"))
load(file=paste0(path,"/sample_B.rda"))


# Height trajectory modeling according to sex
# -------------------------------------------
## girls
result_data_G <- data.frame()
for (i in 1:length(c(unique(subset(sample_G, x==1)$id)))){
  
  # Selection of first subject
  # --------------------------
  listind = c(unique(subset(sample_G, x==1)$id))[[i]]
  
  subjet <- subset(sample_G, x==1 & id==listind)
  subjet <- merge (subjet, as.data.frame(table(subjet$id)), by.x = "id", by.y="Var1")
  subjet <- subset(subjet, Freq>=2)
  subjet <- droplevels(subjet)
  
  # Selection of modeling sample
  # ----------------------------
  tabind <- subset(sample_G, x==0)
  tabind <- merge (tabind, as.data.frame(table(tabind$id)), by.x = "id", by.y="Var1")
  tabind <- subset(tabind)
  tabind <- droplevels(tabind)
  
  
  if (nrow(subjet)>0  & nrow(tabind)>0){
    
    # Number of counter ranging from 1 to j measurements of subject
    # ------------------------------------------------------------
    subjet          <- subjet[order(subjet[['id']], subjet[['age']]),]
    subjet[['cpt']] <- unlist(lapply(table(subjet[['id']]), function(x){1:x}))
    
    
    # Number of counter ranging from 1 to j measurements for each children included in modeling sample
    # ------------------------------------------------------------------------------------------------
    tabind          <- tabind[order(tabind[['id']], tabind[['age']]),]
    tabind[['cpt']] <- unlist(lapply(table(tabind[['id']]), function(x){1:x}))
    
    
    # Maximum number of measurement of subject
    # ---------------------------------------
    nbmesmax <- max(subjet[['Freq']])
    
    
    if(nbmesmax>0){
      
      for (j in 2:max(subjet$Freq)){
        
        # Selection of j measurements of subject
        # --------------------------------------
        subj <- subset(subjet, Freq>=j & cpt<=j)
        
        
        # Selection of measurements of children included in modeling sample / age range
        # --------------------------------------------------------
        if (subset(subj, cpt==j, select=c('age'))[['age']] >= 0    & subset(subj, cpt==j, select=c('age'))[['age']] < 730)   {tabind <- subset(tabind, age < 730)}
        if (subset(subj, cpt==j, select=c('age'))[['age']] >= 730  & subset(subj, cpt==j, select=c('age'))[['age']] < 1095)  {tabind <- subset(tabind, age < 1095)}
        if (subset(subj, cpt==j, select=c('age'))[['age']] >= 1095 & subset(subj, cpt==j, select=c('age'))[['age']] < 1460)  {tabind <- subset(tabind, age < 1460)}
        if (subset(subj, cpt==j, select=c('age'))[['age']] >= 1460 & subset(subj, cpt==j, select=c('age'))[['age']] < 1825)  {tabind <- subset(tabind, age < 1825)}
        if (subset(subj, cpt==j, select=c('age'))[['age']] >= 1825 & subset(subj, cpt==j, select=c('age'))[['age']] < 2920)  {tabind <- subset(tabind, age < 2920)}
        if (subset(subj, cpt==j, select=c('age'))[['age']] >= 2920 & subset(subj, cpt==j, select=c('age'))[['age']] <= 4380) {tabind <- subset(tabind, age <= 4380)}
        
        
        # Merge of data of cases and children included in modeling sample 
        # ----------------------------------------------------------------
        tab_m  <- rbind.data.frame(subj, tabind)
        
        
        # Height growth modeling
        # ---------------------
        
        # data : format SAEMIX
        saemix.data <- saemixData(name.data=tab_m, 
                                  header=TRUE, 
                                  sep=" ", 
                                  na=NA, 
                                  name.group=c("id"), 
                                  name.predictors=c("age"), 
                                  name.response=c("height"), 
                                  units=list(x="days",y="cm"),name.X = "age")
        
        
        
        
        # models : age from 0 to 8 years
        modelT.jens <-function(psi,id,xidep)  {
          age<-xidep[,1]
          A<-psi[id,1]
          B<-psi[id,2]
          C<-psi[id,3]
          D<-psi[id,4]
          Tpredit<-exp(A) + exp(B)*age + exp(C)*(1-exp(-exp(D)*age))
          return(Tpredit)
        }
        saemix.model.jens_0_8<-saemixModel(model=modelT.jens,
                                           description="height trajectory",
                                           psi0=matrix(c(3.9,-4,3.2,-5.5),ncol=4, byrow=TRUE, dimnames=list(NULL, c("A","B","C","D"))),
                                           transform.par=c(0,0,0,0), 
                                           fixed.estim=c(1,1,1,1),
                                           covariance.model=matrix(rep(1,16),ncol=4,byrow=TRUE),
                                           omega.init=matrix(c(0.003,	0.001,	0.001,	0.004, 0.001,	0.017,	0.001,	0.015, 0.001,	0.001,	0.024,	0.024, 0.004,	0.015,	0.024,	0.077) ,ncol=4,byrow=TRUE),
                                           error.model="proportional") 
        
        
        
        # models : age from 8 to 12 years
        modelT.jensA<-function(psi,id,xidep) {
          age<-xidep[,1]
          A<-psi[id,1]
          B<-psi[id,2]
          C<-psi[id,3]
          D<-psi[id,4]
          E<-psi[id,5]
          Tpredit<-exp(A) + exp(B)*age + exp(C)*(1-exp(-exp(D)*age)) + exp(E)*age^2
          return(Tpredit)
        }
        saemix.model.jens.8_12<-saemixModel(model=modelT.jensA,
                                            description="height trajectory",
                                            psi0=matrix(c(3.9,-4,3.2,-5.5,-19.3),ncol=5, byrow=TRUE, dimnames=list(NULL, c("A","B","C","D","E"))),
                                            transform.par=c(0,0,0,0,0), 
                                            fixed.estim=c(1,1,1,1,1),
                                            covariance.model=matrix(rep(1,25),ncol=5,byrow=TRUE),
                                            omega.init=matrix(c(0.003,	0.001,	0.001,	0.004,	0.002, 0.001,	0.017,	0.001,	0.015,	0.001, 0.001,	0.001,	0.024,	0.024,	0.002, 0.004,	0.015,	0.024,	0.077,	0.037, 0.002,	0.001,	0.002,	0.037,	0.051) ,ncol=5,byrow=TRUE),
                                            error.model="proportional") 
        
        
        
        # options 
        saemix.options   <- list(seed=19082,save=TRUE,save.graphs=TRUE,directory="../path/", nb.chains=8, nbiter.saemix=c(500,300))
        
        
        # Application
        ## if age of given subject is from 0 to 8 years
        if (subset(subj, cpt==j, select=c('age'))[['age']] >= 0       & subset(subj, cpt==j, select=c('age'))[['age']] <  2920) {
          saemix.fit.jens    <- saemix(model=saemix.model.jens_0_8, data=saemix.data, control=saemix.options)
          pred_ind.jens <- cbind(id = unique(saemix.fit.jens["data"]["data"][, saemix.fit.jens["data"]["name.group"]]), saemix.fit.jens["results"]["map.psi"])
          colnames(pred_ind.jens)[1] <- saemix.fit.jens["data"]["name.group"]
          result_data   <- data.frame(subset(pred_ind.jens, id==listind, select=c(id,A,B,C,D)),E=0)
        }
        
        ## if age of given subject is from 8 to 12 years
        if (subset(subj, cpt==j, select=c('age'))[['age']] >= 2920    & subset(subj, cpt==j, select=c('age'))[['age']] <= 4380) {
          saemix.fit.jens    <- saemix(model=saemix.model.jens.8_12, data=saemix.data, control=saemix.options)
          pred_ind.jens <- cbind(id = unique(saemix.fit.jens["data"]["data"][, saemix.fit.jens["data"]["name.group"]]), saemix.fit.jens["results"]["map.psi"])
          colnames(pred_ind.jens)[1] <- saemix.fit.jens["data"]["name.group"]
          result_data   <-  data.frame(subset(pred_ind.jens, id==listind, select=c(id,A,B,C,D,E)))
        }
        
        ## results
        result_data_G <- rbind(result_data_G,result_data)

      }
    }
  }
}
tabG         <- anti_join(subset(sample_G, x==1), subset(sample_G, x==1) %>% group_by(id) %>% top_n(-1, age))
resultsG      <- cbind.data.frame(result_data_G, x=tabG[,c("x")], sex=tabG[,c("sex")], age=tabG[,c("age")], 
                                  ts=tabG[,c("ts")], ghd=tabG[,c("ghd")], referents=tabG[,c("referents")], controls=tabG[,c("controls")])

# boys
result_data_B <- data.frame()
for (i in 1:length(c(unique(subset(sample_B, x==1)$id)))){
  
# Selection of first subject
# --------------------------
  listind = c(unique(subset(sample_B, x==1)$id))[[i]]
  
  subjet <- subset(sample_B, x==1 & id==listind)
  subjet <- merge (subjet, as.data.frame(table(subjet$id)), by.x = "id", by.y="Var1")
  subjet <- subset(subjet, Freq>=2)
  subjet <- droplevels(subjet)
  
# Selection of modeling sample
# ----------------------------
  tabind <- subset(sample_B, x==0)
  tabind <- merge (tabind, as.data.frame(table(tabind$id)), by.x = "id", by.y="Var1")
  tabind <- subset(tabind)
  tabind <- droplevels(tabind)
  
  
  if (nrow(subjet)>0  & nrow(tabind)>0){
      
      # Number of counter ranging from 1 to j measurements of subject
      # ------------------------------------------------------------
      subjet          <- subjet[order(subjet[['id']], subjet[['age']]),]
      subjet[['cpt']] <- unlist(lapply(table(subjet[['id']]), function(x){1:x}))
      
      
      # Number of counter ranging from 1 to j measurements for each children included in modeling sample
      # ------------------------------------------------------------------------------------------------
      tabind          <- tabind[order(tabind[['id']], tabind[['age']]),]
      tabind[['cpt']] <- unlist(lapply(table(tabind[['id']]), function(x){1:x}))
      
      
      # Maximum number of measurement of subject
      # ---------------------------------------
      nbmesmax <- max(subjet[['Freq']])
      
      
      if(nbmesmax>0){
        
        for (j in 2:max(subjet$Freq)){
          
          # Selection of j measurements of subject
          # --------------------------------------
          subj <- subset(subjet, Freq>=j & cpt<=j)
          
       
          # Selection of measurements of children included in modeling sample / age range
          # --------------------------------------------------------
          if (subset(subj, cpt==j, select=c('age'))[['age']] >= 0    & subset(subj, cpt==j, select=c('age'))[['age']] < 730)   {tabind <- subset(tabind, age < 730)}
          if (subset(subj, cpt==j, select=c('age'))[['age']] >= 730  & subset(subj, cpt==j, select=c('age'))[['age']] < 1095)  {tabind <- subset(tabind, age < 1095)}
          if (subset(subj, cpt==j, select=c('age'))[['age']] >= 1095 & subset(subj, cpt==j, select=c('age'))[['age']] < 1460)  {tabind <- subset(tabind, age < 1460)}
          if (subset(subj, cpt==j, select=c('age'))[['age']] >= 1460 & subset(subj, cpt==j, select=c('age'))[['age']] < 1825)  {tabind <- subset(tabind, age < 1825)}
          if (subset(subj, cpt==j, select=c('age'))[['age']] >= 1825 & subset(subj, cpt==j, select=c('age'))[['age']] < 2920)  {tabind <- subset(tabind, age < 2920)}
          if (subset(subj, cpt==j, select=c('age'))[['age']] >= 2920 & subset(subj, cpt==j, select=c('age'))[['age']] <= 4380) {tabind <- subset(tabind, age <= 4380)}
          
          
          # Merge of data of cases and children included in modeling sample 
          # ----------------------------------------------------------------
          tab_m  <- rbind.data.frame(subj, tabind)
  

          # Height growth modeling
          # ---------------------
          
            # data : format SAEMIX
            saemix.data <- saemixData(name.data=tab_m, 
                                       header=TRUE, 
                                       sep=" ", 
                                       na=NA, 
                                       name.group=c("id"), 
                                       name.predictors=c("age"), 
                                       name.response=c("height"), 
                                       units=list(x="days",y="cm"),name.X = "age")
            
            
            
            
            # models : age from 0 to 8 years
            modelT.jens <-function(psi,id,xidep)  {
              age<-xidep[,1]
              A<-psi[id,1]
              B<-psi[id,2]
              C<-psi[id,3]
              D<-psi[id,4]
              Tpredit<-exp(A) + exp(B)*age + exp(C)*(1-exp(-exp(D)*age))
              return(Tpredit)
            }
            saemix.model.jens_0_8<-saemixModel(model=modelT.jens,
                                            description="height trajectory",
                                            psi0=matrix(c(3.9,-4,3.2,-5.5),ncol=4, byrow=TRUE, dimnames=list(NULL, c("A","B","C","D"))),
                                            transform.par=c(0,0,0,0), 
                                            fixed.estim=c(1,1,1,1),
                                            covariance.model=matrix(rep(1,16),ncol=4,byrow=TRUE),
                                            omega.init=matrix(c(0.003,	0.001,	0.001,	0.004, 0.001,	0.017,	0.001,	0.015, 0.001,	0.001,	0.024,	0.024, 0.004,	0.015,	0.024,	0.077) ,ncol=4,byrow=TRUE),
                                            error.model="proportional") 
            
            
            
            # models : age from 8 to 12 years
            modelT.jensA<-function(psi,id,xidep) {
              age<-xidep[,1]
              A<-psi[id,1]
              B<-psi[id,2]
              C<-psi[id,3]
              D<-psi[id,4]
              E<-psi[id,5]
              Tpredit<-exp(A) + exp(B)*age + exp(C)*(1-exp(-exp(D)*age)) + exp(E)*age^2
              return(Tpredit)
            }
            saemix.model.jens.8_12<-saemixModel(model=modelT.jensA,
                                             description="height trajectory",
                                             psi0=matrix(c(3.9,-4,3.2,-5.5,-19.3),ncol=5, byrow=TRUE, dimnames=list(NULL, c("A","B","C","D","E"))),
                                             transform.par=c(0,0,0,0,0), 
                                             fixed.estim=c(1,1,1,1,1),
                                             covariance.model=matrix(rep(1,25),ncol=5,byrow=TRUE),
                                             omega.init=matrix(c(0.003,	0.001,	0.001,	0.004,	0.002, 0.001,	0.017,	0.001,	0.015,	0.001, 0.001,	0.001,	0.024,	0.024,	0.002, 0.004,	0.015,	0.024,	0.077,	0.037, 0.002,	0.001,	0.002,	0.037,	0.051) ,ncol=5,byrow=TRUE),
                                             error.model="proportional") 
            
            
            
            # options 
            saemix.options   <- list(seed=19082,save=TRUE,save.graphs=TRUE,directory="../path/", nb.chains=8, nbiter.saemix=c(500,300))
            
            
            # Application
              ## if age of given subject is from 0 to 8 years
              if (subset(subj, cpt==j, select=c('age'))[['age']] >= 0       & subset(subj, cpt==j, select=c('age'))[['age']] <  2920) {
                saemix.fit.jens    <- saemix(model=saemix.model.jens_0_8, data=saemix.data, control=saemix.options)
                pred_ind.jens <- cbind(id = unique(saemix.fit.jens["data"]["data"][, saemix.fit.jens["data"]["name.group"]]), saemix.fit.jens["results"]["map.psi"])
                colnames(pred_ind.jens)[1] <- saemix.fit.jens["data"]["name.group"]
                result_data   <- data.frame(subset(pred_ind.jens, id==listind, select=c(id,A,B,C,D)),E=0)
              }
              
              ## if age of given subject is from 8 to 12 years
              if (subset(subj, cpt==j, select=c('age'))[['age']] >= 2920    & subset(subj, cpt==j, select=c('age'))[['age']] <= 4380) {
                saemix.fit.jens    <- saemix(model=saemix.model.jens.8_12, data=saemix.data, control=saemix.options)
                pred_ind.jens <- cbind(id = unique(saemix.fit.jens["data"]["data"][, saemix.fit.jens["data"]["name.group"]]), saemix.fit.jens["results"]["map.psi"])
                colnames(pred_ind.jens)[1] <- saemix.fit.jens["data"]["name.group"]
                result_data   <- data.frame(subset(pred_ind.jens, id==listind, select=c(id,A,B,C,D,E)))
              }
            
          ## results
          result_data_B <- rbind(result_data_B,result_data) 
          
        }
         
    }
  } 
}
tabB         <- anti_join(subset(sample_B, x==1), subset(sample_B, x==1) %>% group_by(id) %>% top_n(-1, age))
resultsB      <- cbind.data.frame(result_data_B, x=tabB[,c("x")], sex=tabB[,c("sex")], age=tabB[,c("age")], 
                                  ts=tabB[,c("ts")], ghd=tabB[,c("ghd")], referents=tabB[,c("referents")], controls=tabB[,c("controls")])

traindata <- rbind.data.frame(resultsB, resultsG)

# Save results
# ------------
save(resultsG, file=paste0(path,"/resultsG.rda"))
save(resultsB, file=paste0(path,"/resultsB.rda"))
save(traindata,file=paste0(path,"/traindata.rda"))






  
  
