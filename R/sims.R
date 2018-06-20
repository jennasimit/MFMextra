#' @title inverse logit function
#' @param x value to evaluate inverse logit function at
#' @return value of inverse logit function at x
inv.logit.fn <-function(x) return(exp(x)/(1+exp(x)))

#' Convert genotype calls at a SNP to genotype score
#' @param i SNP index
#' @param gcalls data.frame of genotype calls where rows are SNPs and columns are individuals
#' @return genotype score of SNP
 gen.fn <- function(i,gcalls) {
  N <- dim(gcalls)[2]
  G <- numeric(N/3)
  c1 <- seq(1,N,by=3) #AA
  c2 <- seq(2,N,by=3) #AB
  c3 <- seq(3,N,by=3) #BB
  G[which(gcalls[i,c1]==1)] <- 0
  G[which(gcalls[i,c2]==1)] <- 1
  G[which(gcalls[i,c3]==1)] <- 2
  return(G)
  }


#' @title Convert genotype calls, as output from hapgen2, to a genotype score matrix
#' @param gcalls data.frame of genotype calls where rows are SNPs and columns are individuals
#' @export
#' @return genotype score matrix where rows are SNPs and columns are individuals
convert.fn <- function(gcalls) {
 N <- dim(gcalls)[2]-5 # number of individuals*3; first 5 cols are snp info
 m <- dim(gcalls)[1] # number of snps
 snames <- gcalls[,2]
 gcalls <- gcalls[,-(1:5)]
 G <- matrix(0,nrow=m,ncol=N/3)
 rownames(G) <- snames
 colnames(G) <- paste("indiv",1:(N/3),sep="")

 Gmat <- apply(matrix(1:m,ncol=1),1,gen.fn,gcalls)
 G <- t(Gmat)
 rownames(G) <- snames
 colnames(G) <- paste("indiv",1:(N/3),sep="")

  return(G)
                                }


#' @title Generate case-control data with shared controls for traits 1 and 2 
#' @param beta1 vector of model parameters for trait 1: (log(prevalence),causal variant effects)
#' @param beta2 vector of model parameters for trait 2: (log(prevalence),causal variant effects)
#' @param snpG matrix of genotype scores for null data (rows=SNPs, columns=individuals)
#' @param N0 number of shared controls
#' @param N1 number of cases for trait 1
#' @param N2 number of cases for trait 2
#' @param causals1.ind indices of trait 1 causal variants with respect to snpG 
#' @param causals2.ind indices of trait 2 causal variants with respect to snpG 
#' @export
#' @return list consisting of G=genotype matrix (rows=indiv, cols=snps), y=vector of case-control status (0=control,1=trait 1 case, 2=trait 2 case)
phen.gen.fn <-function(beta1=c(-2.3,.2,.2),beta2=c(-2.3,.2,.2),snpG,N0=100,N1=100,N2=100,causals1.ind,causals2.ind) {
#' beta[1] = log(prev)
 N <-dim(snpG)[2] # number of indivs
 n0=0
 n1=0
 n2=0
 G0 <-NULL
 G1 <- NULL
 G2 <- NULL
 j=1
 while( (n0<=N0 | n1 <=N1 | n2 <=N2) & (j <= N) ) {
  p1 <-inv.logit.fn(beta1[1]+sum(beta1[-1]*snpG[causals1.ind,j]))
  p2 <-inv.logit.fn(beta2[1]+sum(beta2[-1]*snpG[causals2.ind,j]))
  u<-runif(2)
   if(u[1]>p1 & u[1]>p2) {
     n0<-n0 + 1
     G0 <-rbind(G0,snpG[,j])
   } else if(u[1]<=p1 & u[1]>p2) {
    n1<-n1 + 1
     G1 <-rbind(G1,snpG[,j])
   } else if(u[1]<=p2 & u[1]>p1) {
     n2<-n2 + 1
     G2 <-rbind(G2,snpG[,j])
    } else {
     if(u[2]<0.5 & n1<N1 ) { n1<-n1 + 1; G1 <-rbind(G1,snpG[,j]) 
     } else { n2<-n2 + 1; G2 <-rbind(G2,snpG[,j]) }
    }
j <- j+1
	}
G<-rbind(G0[1:N0,],G1[1:N1,],G2[1:N2,])
rownames(G)<-c(paste("control.",1:N0,sep=""),paste("case1.",1:N1,sep=""),paste("case2.",1:N2,sep=""))
colnames(G) <- rownames(snpG)
 return(list(G=G,y=c(rep(0,N0),rep(1,N1),rep(2,N2)) ))
  	  	     }


#' @title Generate case-control data with shared controls for traits 1 and 2, where trait 2 has no associations in this region
#' @param beta1 vector of model parameters for trait 1: (log(prevalence),causal variant effects)
#' @param snpG matrix of genotype scores for null data (rows=SNPs, columns=individuals)
#' @param N0 number of shared controls
#' @param N1 number of cases for trait 1
#' @param N2 number of cases for trait 2
#' @param causals1.ind indices of trait 1 causal variants with respect to snpG
#' @export
#' @return list consisting of G=genotype matrix (rows=indiv, cols=snps), y=vector of case-control status (0=control,1=trait 1 case, 2=trait 2 case)
phen.gen.t2null.fn <- function (beta1 = c(-2.3, 0.2, 0.2), snpG, N0 = 100, N1 = 100, N2 = 100, causals1.ind) {    
    N <- dim(snpG)[2]
    n0 = 0
    n1 = 0
    G0 <- NULL
    G1 <- NULL
    G2 <- NULL       
    j = 1
    while ((n0 <= (N0+N2) | n1 <= N1 ) & (j <= N)) {
        p1 <- inv.logit.fn(beta1[1] + sum(beta1[-1] * snpG[causals1.ind, j]))       
        u <- runif(1)
        if (u[1] > p1) {
            n0 <- n0 + 1
            G0 <- rbind(G0, snpG[, j])
        }
        else  {
            n1 <- n1 + 1
            G1 <- rbind(G1, snpG[, j])
        }                											
        j <- j + 1    							}												
    G <- rbind(G0[1:N0, ], G1[1:N1, ], G0[(N0+1):(N0+N2), ])
    rownames(G) <- c(paste("control.", 1:N0, sep = ""), paste("case1.",1:N1, sep = ""), paste("case2.", 1:N2, sep = ""))
    colnames(G) <- rownames(snpG)
    return(list(G = G, y = c(rep(0, N0), rep(1, N1), rep(2, N2))))
}
