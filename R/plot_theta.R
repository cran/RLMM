plot_theta <- function(genotypefile="Xba.rlmm",thetafile="Xba.theta",Pick.Obj="FALSE",plotfile="plots.ps",snpsfile="snps.lst")
{

  if(Pick.Obj==FALSE){
    dat <- read.table(thetafile,as.is=T)
    geno <- read.table(genotypefile,as.is=T)
  }

NUMSAMP<- (ncol(geno)-1)/2
NUMSNP<- nrow(geno)

gn.idx <- seq(2,(2*NUMSAMP+1),2) #-- for Xba.rlmm

if(plotfile!=""){
postscript(plotfile)
}

par(pty="s")
par(mfrow=c(1,1))

babe <- read.table(snpsfile,as.is=T)

match(babe[,1],geno[,1])->idx
print("Matching at Index")
print(idx)

geno<- geno[,gn.idx]	#Extract only the genotype columns

if ( sum(is.na(idx)) == length(idx) ) #none of the requested SNPs are found
{
  print("Could not find the theta values for the SNPs requested for plotting ")
}
if ( sum(is.na(idx)) != length(idx) ) #some of the requested SNPs are found
{
 idx <- idx[!is.na(idx)]
 dat<-dat[idx,]		#Extract only the matched SNP rows
 geno<-geno[idx,]

 THETA.A <- dat[,2:(NUMSAMP+1)]
 THETA.B <- dat[,(NUMSAMP+2):(2*NUMSAMP+1)]

Ac <- array(-1,ncol(geno))
ec <- array(-1,ncol(geno))

for (i in 1:nrow(THETA.A))	#for each SNP in the list
{
g <- as.character(geno[i,])
ind.AA <- g=="AA"
ind.AB <- g=="AB"
ind.BB <- g=="BB"

num.AA <- sum(ind.AA)
num.AB <- sum(ind.AB)
num.BB <- sum(ind.BB)
T.A <- as.numeric(THETA.A[i,])
T.B <- as.numeric(THETA.B[i,])

#DO THE RLMM PLOT
#nc <- round(sum(n=="NC")/NUMSAMP * 100,2)
Ac[g=="AA"] <- 4
Ac[g=="AB"] <- 2
Ac[g=="BB"] <- 3
Ac[g=="NC"] <- 1
ec[g=="AA"] <-24
ec[g=="AB"] <-23
ec[g=="BB"] <-25
ec[g=="NC"] <- 1
label.A <- paste("allele A")
label.B <- paste("allele B")
snpid <- dat[i,1]
plot(T.A,T.B,col=Ac,pch=ec,main=paste(snpid," Allele Summary Plot",sep=" "),xlab=label.A,ylab=label.B,xlim=c(8,15),ylim=c(8,15))

legend (13,15,legend=c("AA","AB","BB","NC"),pch=c(24,23,25,1),col=c(4,2,3,1))
} # end for loop
} # end some of the requested SNPs are found

if(plotfile!=""){
   dev.off()
   }

}
