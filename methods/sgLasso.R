#This script is used for running parallelized calculations of sgLasso regression model

##################################################################
#             Fonction Main resampling                           #
##################################################################

#Important: First, the database should be shaped according to the SGL package recommendations

resampling_function = function(database,nb_iterations){

set.seed(1)  #We fix a seed for reproducibility matters
  
library(foreach)
library(doParallel)
library(doRNG)             
  
#Initialization
vec_list_var = c()
vec_nb_var = c()
index = c()
cl = makeCluster(detectCores())                    
registerDoParallel(cl, cores = detectCores())       
  
#################################################################
#  Repeated double cross-validation function : result_sampling 
#################################################################
  
result_sampling  =  foreach(i=1:nb_iterations) %dorng% {
  
                library(SGL) 
                library(R.utils)
                library(psych)
                
              
                #Initialization
                nb_var_vec = c()
                name_var = c()
                M=10
                n = nrow(database)
                MSEP.vec = Ypred =  rep(0,n)
                folds = split(sample(1:n), rep(1:M, length = n))      # 10 fold CV samples 
                alphas = seq(0.05, 0.95, by=.1)
                grp = c(rep(1,length(grep("_EC", colnames(database)))),rep(2,length(grep("_PC", colnames(database)))),rep(3,length(grep("_pl", colnames(database)))),rep(4,length(grep("_GR", colnames(database)))),rep(5,length(grep("_lcms", colnames(database)))))
                
                #Begin outerloop
                
                for(j in 1:M){
                  
                            #Split the data into training and test sets
                            omit = folds[[j]]
                            data.train = database[-omit,]
                            data.test = database[omit,]
                            x.train = data.train[,-1]
                            y.train = data.train[,1]
                            x.test = data.test[,-1]
                            y.test = data.test[,1]
                            
                            # Begin inner loop
                            #Training the model to obtain the best tuning parameters
                            mses = numeric(length(alphas))
                            lambda.min = numeric(length(alphas))
                            SGL_cv = list()
                            length(SGL_cv) = length(alphas)
                            dataset = list(x=x.train,y=y.train)
                            
                            for(k in 1:length(alphas)) {
                              
                                                      SGL_cv[[k]] = cvSGL(data=dataset,index=grp,type="linear",alpha=alphas[k],nfold=n) 
                                                      mses[k] = min(SGL_cv[[k]]$lldiff) 
                                                      lambda.min[k] = SGL_cv[[k]]$fit$lambdas[which.min(SGL_cv[[k]]$lldiff)]
                                                      
                                                        }
                            
                            #Saving the outputs: names and number of predictors associated to the model with the optimal tuning parameters
                            idx =  which.min(SGL_cv[[which.min(mses)]]$lldiff)
                            coefs_sglasso = SGL_cv[[which.min(mses)]]$fit$beta[,idx]
                            names(coefs_sglasso) = colnames(x.train) 
                            coefs_sglasso = coefs_sglasso[coefs_sglasso!=0] 
                            nb_var_vec = c(nb_var_vec,length(coefs_sglasso))
                            name_var = c(name_var,names(coefs_sglasso))
                            
                            # End inner loop
                            
                            #saving prediction error (MSEP) calculated on the validation set 
                            sglasso = SGL(data=dataset,index=grp,type="linear",alpha=alphas[which.min(mses)],lambdas=SGL_cv[[which.min(mses)]]$fit$lambdas)   
                            y.hat = predictSGL(x=sglasso,newX=x.test,lam=idx)
                            Ypred[omit] = y.hat
                            MSEP.vec[omit] = (y.test - y.hat)^2
                            
                          }
                
                #End outerloop
                
                #Prediction indices calculation and saving
                RMSEP =  sqrt(mean(MSEP.vec))
                R2 = cor(database[,1],Ypred)^2
                
                return(list(RMSEP=RMSEP,R2=R2,Nb_var=nb_var_vec,Num_ind=folds,Var_name=name_var))
                                                    }
  
#Formatting the outputs
for(l in 1:nb_iterations){
    vec_list_var = c(vec_list_var,result_sampling[[l]]$Var_name)
    vec_nb_var = c(vec_nb_var,result_sampling[[l]]$Nb_var)
    index = cbind(index,c(round(result_sampling[[l]]$RMSEP,2),round(result_sampling[[l]]$R2,2)))
                          }
  
table_index = as.data.frame(index)
rownames(table_index) = c("RMSEP","R2")
colnames(table_index) = paste("iter",1:ncol(table_index),sep="")
  
return(list(list_var=vec_list_var,nb_var=vec_nb_var,indexes=table_index))
}

#Calling resampling_funtion
sgLasso = resampling_function(database,100)

# Obtaining variables selection frequency with a threshold = 0.6
sort(table(sgLasso$list_var)/1000,decreasing=T)
fqcy_60 = length(table(sgLasso$list_var)[table(sgLasso$list_var)/1000>=0.6])
list_60 =  names(table(sgLasso$list_var)[table(sgLasso$list_var)/1000>=0.6])

#Plot of the variable selection frequency 
source("barplot_fct.R")
barplot_sgLasso = barplot_fct(sgLasso)
