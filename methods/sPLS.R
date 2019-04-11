
#This script is used for running parallelized calculations of sPLS regression model

##################################################################
#             Fonction Main resampling                           #
##################################################################

#Important: First, the database should be shaped according to the sgPLS package recommendations

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

                library(sgPLS)
                library(R.utils)
                library(caret)
                
                #Initialization
                nb_var_vec = c()
                name_var = c()
                M=10
                n = nrow(database)
                MSEP.vec = Ypred =  rep(0,n)
                folds = split(sample(1:n), rep(1:M, length = n))     
                l = c(1:10,15,20,50,100)
                
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
                              msep.all = NULL
                              for(k in l){
                                        result = sPLS(x.train,y.train,ncomp = 1,keepX = k , mode = "regression")           
                                        loo = perf(result, criterion = "MSEP", validation = 'loo')
                                        msep.all = c(msep.all,loo$MSEP)
                                        }
                  
                              #Saving the outputs: names and number of predictors associated to the model with the optimal tuning parameters
                              res.spls =sPLS(x.train,y.train, ncomp = 1, keepX =l[which.min(msep.all)], mode = "regression")
                              CP1_var_name = names(select.sgpls(res.spls)$select.X[[1]])
                              nb_var_vec = c(nb_var_vec,length(CP1_var_name))
                              name_var = c(name_var,CP1_var_name)
                              
                              # End inner loop
                              
                              #prediction error calculated on the independant test set
                              y.hat = predict(res.spls ,newdata=apply(x.test,2,as.numeric))$predict
                              Ypred[omit] = y.hat
                              MSEP.vec[omit] = (y.test - y.hat)^2
                              
                              }
                
                #End outerloop
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
Spls = resampling_function(database,100)

# Obtaining variables selection frequency with a threshold = 0.6
sort(table(Spls$list_var)/1000,decreasing=T)
fqcy_60 = length(table(Spls$list_var)[table(Spls$list_var)/1000>=0.6])
list_60 =  names(table(Spls$list_var)[table(Spls$list_var)/1000>=0.6])

#Plot of the variable selection frequency
source("barplot_fct.R") 
barplot_Spls = barplot_fct(Spls)













