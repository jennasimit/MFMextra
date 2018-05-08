
#' @title Generate model matrix of all possible models with maximum specified size
#' @param mT1 maximum model size
#' @param msnps SNPs for model matrix
#' @return model matrix for SNPs specified in msnps; columns are SNPs, each row corresponds to a model 
T1modsmat.fn <- function(mT1,msnps) {
 myLetters <- letters[1:26]
 s <- length(msnps)

 modT1 <- vector("list",mT1)
  for(i in 1:mT1) {
   mci <- combn(letters[1:s],i,simplify=TRUE)  # all combinations of i msnps
   nci <- dim(mci)[2]
   modT1[[i]] <- matrix(0,nrow=nci,ncol=s)
   for(j in 1:nci)  modT1[[i]][j,match(mci[,j],myLetters)] <- 1 
  }
 modsT1 <- rep(0,s)
 for(i in 1:mT1) modsT1 <- rbind(modsT1,modT1[[i]])
 T1mods <- modsT1
return(T1mods)
}

#' @title Approximate Bayes' factors (ABFs) for one case-control study
#' @param sim output from phen.gen.fn
#' @export
#' @return list of data.frames for each trait; column 1 gives case-control status (1/0) and remaining columns form a SnpMatrix object
abfcalc.format <- function(sim) { 
 Gm <- new("SnpMatrix",(sim$G+1)) # snp cols, indivs rows
 c0 <- grep("control.",rownames(Gm))
 c1 <- grep("case1.",rownames(Gm))
 c2 <- grep("case2.",rownames(Gm))

 G1 <- Gm[c(c0,c1),]
 G2 <- Gm[c(c0,c2),]

 data1.1 <- data.frame(Y=sim$y[c(c0,c1)],sim$G[c(c0,c1),])
 data1.2 <- data.frame(Y=sim$y[c(c0,c2)],sim$G[c(c0,c2),])
 data1.2$Y[data1.2$Y==2] <- 1
return(list(data1.1,data1.2))
}


#' @title Approximate Bayes' factors (ABFs) for one case-control study
#' @param data1 data.frame that has case-control status (1-0) in column 1 and remaining columns are the genotype scores (rows are individuals)
#' @param mT1 maximum number of causal variants
#' @param msnps vector of SNPs to consider in models
#' @export
#' @return data.frame of model ABFs (each row is a model): column 1 has ABFs, remaining columns specify inclusion/excusion of SNPs in each model
abfT1.fn <- function(data1,mT1=3,msnps) {

mods <- T1modsmat.fn(mT1,msnps)
colnames(mods) <- msnps

mod1 <- BMA::glib(x=data1[,-1],y=data1[,1],error="binomial", link="logit",models=mods)

logABF <- mod1$bf$twologB10[,1]*0.5

out <- data.frame(logABF=logABF,M=mod1$models,row.names=NULL)

cnames <- c("logABF",names(data1)[-1])
names(out) <- cnames
return(out)
}


#' @title Joint Approximate Bayes' factors (ABFs) for two case-control studies with shared controls
#' @param sim output from phen.gen.fn
#' @param msnps vector of SNPs to consider in models
#' @param mT1 maximum number of causal variants for trait 1
#' @param mT2 maximum number of causal variants for trait 2
#' @export
#' @return data.frame of joint model ABFs (each row is a model): column 1 has joint ABFs, remaining columns specify inclusion/excusion of SNPs in each joint model
abf.fn <- function(sim,msnps,mT1,mT2) {

data1 <- data.frame(Y=sim$y,sim$G)
s <- length(msnps) # number of snps considered in models for the 2 traits
data.m1 <- data1[,c("Y",msnps)]
m1 <- mlogit2logit(Y ~ 1|. -Y,data.m1,choices=0:2,base.choice=1)

T1mods <- T1modsmat.fn(mT1,msnps)
if(mT2==mT1) {T2mods <- T1mods} else {T2mods <- T1modsmat.fn(mT2,msnps)}

nT1 <- dim(T1mods)[1]
nT2 <- dim(T2mods)[1]

T1modsrep <- matrix(rep(t(T1mods),nT2),ncol=ncol(T1mods),byrow=TRUE)
T2modsrep <- matrix(rep(T2mods,each=nT1),ncol=ncol(T2mods),byrow=FALSE)
T1T2mods <- cbind(T1modsrep,T2modsrep)

#' add column of 1's to the models for the trait2*effect variable
T1T2mods1 <- cbind(T1T2mods,1)
T1T2mods <- T1T2mods1

mod1 <- BMA::glib(x=m1$data[,(4+s):(4+3*s)],y=m1$data$Y.star,error="binomial", link="logit",models=T1T2mods)

logABF <- mod1$bf$twologB10[,1]*0.5

out <- data.frame(logABF=logABF,M=mod1$models)

cnames <- c("logABF",names(m1$data[,(4+s):(4+3*s)]))
names(out) <- cnames

out <- out[,-dim(out)[2]] # rm last column, which is 1 for z_2

return(out) # logABF,mod
}


#' @title Joint Approximate Bayes' factors (ABFs) for two case-control studies with shared controls
#' @param sim output from phen.gen.fn
#' @param msnps vector of SNPs to consider in models
#' @param mT1 maximum number of causal variants for trait 1
#' @param mT2 maximum number of causal variants for trait 2
#' @export
#' @return data.frame of multinomial logABFs and logistic logABFs
bf.compare.fn <- function(sim,msnps,mT1,mT2) {

#' find marginal ABFs for each trait
data12 <- abfcalc.format(sim)
bft1 <- abfT1.fn(data12[[1]],mT1=mT1,msnps)
bft2 <- abfT1.fn(data12[[2]],mT1=mT2,msnps)

#' find joint ABFs for traits 1 and 2
bft1t2 <- abf.fn(sim,msnps=msnps,mT1=mT1,mT2=mT2)

t1mods <- as.matrix(bft1[,-1])
t2mods <- as.matrix(bft2[,-1])
s <- dim(t1mods)[2]
t1t2mods <- as.matrix(bft1t2[,-1])
colnames(t1mods) <- NULL
colnames(t2mods) <- NULL
colnames(t1t2mods) <- NULL

bfall <- c()

for(k in 1:dim(t1mods)[1]) {
 for(j in 1:dim(t2mods)[1]) {
ind <- which(apply(t1t2mods, 1, identical, c(t1mods[k,],t2mods[j,])))
tmp <- cbind(bft1t2[ind,],logBF1=bft1[k,1],logBF2=bft2[j,1])
bfall <- rbind(bfall, tmp)
}
}
return(data.frame(BF12=bfall[,1],BF1=bfall[,"logBF1"],BF2=bfall[,"logBF2"])) 
} 


#' @title Relationship between multinomial ABFs and logistic ABFs
#' @param BFs data.frame of multinomial logABFs BF12 and logistic logABFs
#' @param Bplot logical, if TRUE plot ABF12 against (log(ABF1)+log(ABF2)
#' @export
#' @return fitted regression model R2 and coefficient estimate summaries from log(jointABF) ~ (log(ABF1)+log(ABF2))
bf.relations.fn <- function(BFs,Bplot=FALSE) {
b1.b2 <- BFs$BF1 + BFs$BF2
out <- lm(BFs$BF12~b1.b2)
R2 <- summary(out)$r.squared
betas <- summary(out)$coefficients[,1:2]
output <- c(R2=R2,beta0=betas[1,],beta1=betas[2,])
if(Bplot) plot(BFs$BF12~b1.b2,pch=20,xlab="log(BF1)+log(BF2)",ylab="log(BF12)",main=paste("N = ",N))
return(output)
}

