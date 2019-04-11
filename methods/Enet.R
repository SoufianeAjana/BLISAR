
#This script is used for running parallelized calculations of Elastic net regression model

##################################################################
#             Fonction Main resampling                           #
##################################################################

#Important: First, the database should be shaped according to the glmnet package recommendations

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

                library(glmnet)  # Package to fit lasso/elastic net models
                library(R.utils)
                    
                #Initialization
                nb_var_vec = c()
                name_var = c()
                alpha_vec = c()
                M=10
                n = nrow(database)
                MSEP.vec = Ypred =  rep(0,n)
                folds = split(sample(1:n), rep(1:M, length = n))    
                alphas = seq(0.05, 0.95, by=.1)
                    
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
                            enet_cv = list()
                            length(enet_cv) = length(alphas)
                      
                            for(k in 1:length(alphas)) {
                              
                                                      enet_cv[[k]] = cv.glmnet(x.train, y.train,alpha=alphas[k],parallel=T,nfolds=n) 
                                                      mses[k] = enet_cv[[k]]$cvm[which(enet_cv[[k]]$lambda==enet_cv[[k]]$lambda.min)] 
                                                       
                                                        }
                      
                       
                            #Saving the outputs: names and number of predictors associated to the model with the optimal tuning parameters
                            coefs_enet = coef(enet_cv[[which.min(mses)]], s = "lambda.min")[which(coef(enet_cv[[which.min(mses)]], s = "lambda.min")!=0)][-1]  
                            names(coefs_enet) = rownames(coef(enet_cv[[which.min(mses)]], s = "lambda.min"))[which(coef(enet_cv[[which.min(mses)]], s = "lambda.min")!=0)][-1]
                            nb_var_vec = c(nb_var_vec,length(coefs_enet))
                            name_var = c(name_var,names(coefs_enet))
                      
                            # End inner loop
                                           
                      
                            #saving prediction error (MSEP) calculated on the validation set
                            y.hat = predict(enet_cv[[which.min(mses)]],newx=x.test,s="lambda.min")
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
Enet = resampling_function(database,100)


# Obtaining variables selection frequency with a threshold = 0.6
sort(table(Enet$list_var)/1000,decreasing=T)
fqcy_60 = length(table(Enet$list_var)[table(Enet$list_var)/1000>=0.6])
list_60 =  names(table(Enet$list_var)[table(Enet$list_var)/1000>=0.6])

#Plot of the variable selection frequency 
source("barplot_fct.R")
barplot_Enet = barplot_fct(Enet)
