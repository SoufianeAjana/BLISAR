library(ggplot2)

#Load the outputs of all the methods

load("sPLS.RDATA")
load("gPLS.RDATA")
load("sgPLS.RDATA")
load("Lasso.RDATA")
load("gLasso.RDATA")
load("sgLasso.RDATA")
load("Enet.RDATA")

#RMSEP

delta_GPLS <- as.numeric(Gpls$indexes[1,])-as.numeric(SGpls$indexes[1,])
delta_SPLS <- as.numeric(Spls$indexes[1,])-as.numeric(SGpls$indexes[1,])
delta_Lasso <- as.numeric(Lasso$indexes[1,])-as.numeric(SGpls$indexes[1,])
delta_gLasso <- as.numeric(gLasso$indexes[1,])-as.numeric(SGpls$indexes[1,])
delta_sgLasso <- as.numeric(sgLasso$indexes[1,])-as.numeric(SGpls$indexes[1,])
delta_Enet <- as.numeric(Enet$indexes[1,])-as.numeric(SGpls$indexes[1,])

delta_RMSEP <- c(delta_Lasso,delta_gLasso,delta_sgLasso,delta_Enet,delta_SPLS,delta_GPLS)

method <- as.factor(c(rep("Lasso",100),rep("gLasso",100),rep("sgLasso",100),rep("Elastic net",100),rep("sPLS",100),rep("gPLS",100)))

dataset <- data.frame(delta_RMSEP=delta_RMSEP,method=method)

p10 <- ggplot(dataset, aes(x = method, y = delta_RMSEP)) +
        geom_boxplot()+ geom_hline(yintercept = 0) 

#change y axis label

p10 <- p10 + scale_x_discrete(name = "Method") +
        scale_y_continuous(name = "RMSEP")
        
        
        
#R²             

delta_GPLS <- as.numeric(Gpls$indexes[2,])-as.numeric(SGpls$indexes[2,])
delta_SPLS <- as.numeric(Spls$indexes[2,])-as.numeric(SGpls$indexes[2,])
delta_Lasso <- as.numeric(Lasso$indexes[2,])-as.numeric(SGpls$indexes[2,])
delta_gLasso <- as.numeric(gLasso$indexes[2,])-as.numeric(SGpls$indexes[2,])
delta_sgLasso <- as.numeric(sgLasso$indexes[2,])-as.numeric(SGpls$indexes[2,])
delta_Enet <- as.numeric(Enet$indexes[2,])-as.numeric(SGpls$indexes[2,])

delta_RMSEP <- c(delta_Lasso,delta_gLasso,delta_sgLasso,delta_Enet,delta_SPLS,delta_GPLS)

method <- as.factor(c(rep("Lasso",100),rep("gLasso",100),rep("sgLasso",100),rep("Elastic net",100),rep("sPLS",100),rep("gPLS",100)))

dataset <- data.frame(delta_RMSEP=delta_RMSEP,method=method)

p10 <- ggplot(dataset, aes(x = method, y = delta_RMSEP)) +
        geom_boxplot()+ geom_hline(yintercept = 0) 

#change y axis label

p10 <- p10 + scale_x_discrete(name = "Method") +
        scale_y_continuous(name = "R²")
        