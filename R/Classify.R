Classify <- function(genotypefile="",regionsfile="",thetafile="",callrate=100)
{

md.classify <- function(theta.A,theta.B,m.AA,m.AB,m.BB,s.AA,s.AB,s.BB,callrate,threshold,debug)
{

theta <- cbind(theta.A,theta.B)
md.AA <- mahalanobis(theta,m.AA,s.AA)
md.AB <- mahalanobis(theta,m.AB,s.AB)
md.BB <- mahalanobis(theta,m.BB,s.BB)
geno <- array(-1,length(md.AA))
mind<- array(-1,length(md.AA))
res <- ""
for (i in 1:length(md.AA))
{
### Always make a call
mind[i]<- min (md.AA[i],md.AB[i],md.BB[i])
if (mind[i]==md.AA[i]) geno[i] <- "AA"
if (mind[i]==md.AB[i]) geno[i] <- "AB"
if (mind[i]==md.BB[i]) geno[i] <- "BB"
if ((callrate <100) & (mind[i]>threshold)) {geno[i] <-"NC"}
res <- paste(res,geno[i],round(mind[i],4),sep=" ")
} #for i loop

return(res)
} # end function md.classify

calc.meancov <-function(num.AA,num.AB,num.BB,M.AA,M.AB,M.BB,sigma.AA,sigma.AB,sigma.BB,B,al,BnB,Bal,CB,Cal, T.A,T.B, SAA, SAB, SBB,rate,thresh)
{

if (num.AA <= 5)
{
if (num.AB>1 & num.BB>1) {M.AA <- B%*%   t(t(c(M.AB,M.BB))) + t(t(al))}
}
if (num.AB <= 5)
{
if (num.AA>1 & num.BB>1) {M.AB <-BnB%*% t(t(c(M.AA,M.BB))) + t(t(Bal))}
}
if (num.BB <= 5)
{
if (num.AA>1 & num.AB>1) {M.BB <-CB%*%  t(t(c(M.AA,M.AB))) + t(t(Cal))}
}
if (num.BB <= 5 | num.AB <= 5 | num.AA <= 5 )
{S.AA <- SAA
 S.AB <- SAB
 S.BB <- SBB
}
else 
{
if (num.AA>5)
{
S.AA <- sigma.AA + .008		#Add constant
}
if (num.AB>5)
{
S.AB <- sigma.AB + .008
}
if (num.BB>5)
{
S.BB <- sigma.BB + .008
}
} # end else

res <- md.classify(T.A,T.B,M.AA,M.AB,M.BB,S.AA,S.AB,S.BB,rate,thresh,dflag)
return(res)
} #end function


###### START THE MAIN PROGRAM ###########

if(thetafile=="")
   stop("Please specify a filename for the thetafile parameter")

if(genotypefile=="")
   stop("Please specify a filename for the genotypefile parameter")

if(regionsfile=="")
   stop("Please specify a filename for the regionsfile parameter")


#Read the SNP names in the thetafile
snpfn <- read.table(thetafile,as.is=T)
snpfn <- snpfn[,1]

#Process call rates
rates <- seq(80,100,2)
m <- is.na(match(callrate,rates))
if (m ==TRUE ) {print("Setting Call Rate to 80%"); rateindex<-1; callrate<-80} 
else {status<-paste("Setting Call Rate to",callrate); print(status);which(callrate==rates)->rateindex}
thresholds<-c(2.606, 2.787,2.996,3.236,3.523,3.873,4.318,4.929,5.874,7.811,298.484)
thresh <- thresholds[rateindex]
#END - Process call rates

#BEGIN - Now read in regression coefficients 
#Do some processing with the grand covariance matrix - (from grandcov.out)

M <- read.table(regionsfile,nrows=1,as.is=T)
M <- as.numeric(M)
dat <- read.table(regionsfile,nrows=9,as.is=T,skip=2)
S <- matrix(NA,9,9)
for (i in 1:9) { S[i,] <- as.numeric(dat[i,]) }
SAA <- matrix(NA,2,2)
SAB <- matrix(NA,2,2)
SBB <- matrix(NA,2,2)
SAA[1,1] <- M[1]
SAA[2,2] <- M[2]
SAA[1,2] <- SAA[2,1] <- (M[3])*sqrt(SAA[1,1])*sqrt(SAA[2,2])
SAB[1,1] <- M[4]
SAB[2,2] <- M[5]
SAB[1,2] <- SAB[2,1] <- (M[6])*sqrt(SAB[1,1])*sqrt(SAB[2,2])
SBB[1,1] <- M[7]
SBB[2,2] <- M[8]
SBB[1,2] <- SBB[2,1] <- (M[9])*sqrt(SBB[1,1])*sqrt(SBB[2,2])


#Do some processing with the grandmean matrix - (from grandmean.out)
dat <- read.table(regionsfile,as.is=T,skip=18,nrows=14)
al1 <- dat[1,]
al2 <- dat[2,]

b11 <- dat[3,]
b12 <- dat[4,]
b13 <- dat[5,]
b14 <- dat[6,]
b21 <- dat[7,]
b22 <- dat[8,]
b23 <- dat[9,]
b24 <- dat[10,]

g11 <- dat[11,]
g12 <- dat[12,]
g21 <- dat[13,]
g22 <- dat[14,]

s1g2 <- matrix(c(g11,g12,g21,g22),byrow=T,nrow=2)
B <- matrix(c(b11,b12,b13,b14,b21,b22,b23,b24),byrow=T,nrow=2)
al <- c(al1,al2)

#Do some processing with the 2nd grandmean matrix - (from grandmean2.out)
Bdat <- read.table(regionsfile,as.is=T,skip=33,nrows=14)
Bal1 <- Bdat[1,]
Bal2 <- Bdat[2,]
Bb11 <- Bdat[3,]
Bb12 <- Bdat[4,]
Bb13 <- Bdat[5,]
Bb14 <- Bdat[6,]
Bb21 <- Bdat[7,]
Bb22 <- Bdat[8,]
Bb23 <- Bdat[9,]
Bb24 <- Bdat[10,]

Bg11 <- Bdat[11,]
Bg12 <- Bdat[12,]
Bg21 <- Bdat[13,]
Bg22 <- Bdat[14,]

Bs1g2 <- matrix(c(Bg11,Bg12,Bg21,Bg22),byrow=T,nrow=2)
BnB <- matrix(c(Bb11,Bb12,Bb13,Bb14,Bb21,Bb22,Bb23,Bb24),byrow=T,nrow=2)
Bal <- c(Bal1,Bal2)

#Do some processing with the 3rd grandmean matrix - (from grandmean3.out)
Cdat <- read.table(regionsfile,as.is=T,skip=48,nrows=14)
Cal1 <- Cdat[1,]
Cal2 <- Cdat[2,]

Cb11 <- Cdat[3,]
Cb12 <- Cdat[4,]
Cb13 <- Cdat[5,]
Cb14 <- Cdat[6,]
Cb21 <- Cdat[7,]
Cb22 <- Cdat[8,]
Cb23 <- Cdat[9,]
Cb24 <- Cdat[10,]
Cg11 <- Cdat[11,]
Cg12 <- Cdat[12,]
Cg21 <- Cdat[13,]
Cg22 <- Cdat[14,]

Cs1g2 <- matrix(c(Cg11,Cg12,Cg21,Cg22),byrow=T,nrow=2)
CB <- matrix(c(Cb11,Cb12,Cb13,Cb14,Cb21,Cb22,Cb23,Cb24),byrow=T,nrow=2)
Cal <- c(Cal1,Cal2)
#END - Now read in regression coefficients 

#BEGIN - READ THE FILES
regionfile.obj <- read.table(regionsfile,as.is=T,skip=63)
thetafile.obj <- read.table(thetafile,as.is=T)
print(dim(thetafile.obj))
#END - READ THE FILES

NUMSAMP <- (ncol(thetafile.obj)-1)/2	# Number of chips
NUMSNP <- nrow(regionfile.obj)		# Number of SNPs in regionsfile


# Match the SNPs in the thetafile to the SNPs in the regions file
snpid.thetafile<-thetafile.obj[,1]
print(snpid.thetafile)
snpid.regionfile<-regionfile.obj[,1]
m <- match(snpid.thetafile,snpid.regionfile)

# END - Match the SNPs in the thetafile to the SNPs in the regions file

#BEGIN - MAIN LOOP
for (r in 1:nrow(thetafile.obj)) # start looping over the SNPs in thetafile
{ 
discontinue <- FALSE

snpid <- thetafile.obj[r,1]
n <- m[r]	#corresponding index in the regions file
if (is.na(n)) {discontinue<-TRUE; error<-paste("Can't find ",snpid," in the regions file"); print(error)}
else {		#If match obtained for SNP in the regions file
num.AA <- regionfile.obj[n,2]
num.AB <- regionfile.obj[n,3]
num.BB <- regionfile.obj[n,4]
M.AA <- as.numeric(regionfile.obj[n,5:6])
M.AB <- as.numeric(regionfile.obj[n,7:8])
M.BB <- as.numeric(regionfile.obj[n,9:10])
S.AA <- as.numeric(regionfile.obj[n,11:13])
S.AB <- as.numeric(regionfile.obj[n,14:16])
S.BB <- as.numeric(regionfile.obj[n,17:19])
S.AA <- matrix(c(S.AA[1],S.AA[3],S.AA[3],S.AA[2]),nrow=2,byrow=T)
S.AB <- matrix(c(S.AB[1],S.AB[3],S.AB[3],S.AB[2]),nrow=2,byrow=T)
S.BB <- matrix(c(S.BB[1],S.BB[3],S.BB[3],S.BB[2]),nrow=2,byrow=T)


  #Exclude singletons or monomorphic SNPs
  zero.AA <- num.AA==0
  zero.AB <- num.AB==0
  zero.BB <- num.BB==0
  one.AA <- num.AA==1
  one.AB <- num.AB==1
  one.BB <- num.BB==1
 
 if ((zero.AA+zero.AB==2) | (zero.AB+zero.BB==2) | (zero.AA+zero.BB==2)) {discontinue <- TRUE}
 if ((zero.AA+one.AB==2) | (zero.AB+one.BB==2) | (zero.AA+one.BB==2)) {discontinue <- TRUE}
 if ((one.AA+zero.AB==2) | (one.AB+zero.BB==2) | (one.AA+zero.BB==2)) {discontinue <- TRUE}
  #END - Exclude singletons or monomorphic SNPs

} # END - If match obtained for SNP in the regions file

if (discontinue==0)
{
g <- array("-1",NUMSAMP)

THETA.A <- as.numeric(thetafile.obj[r,2:(NUMSAMP+1)])
THETA.B <- as.numeric(thetafile.obj[r,(NUMSAMP+2):(2*NUMSAMP+1)])


res<- calc.meancov(num.AA,num.AB,num.BB,M.AA,M.AB,M.BB,S.AA,S.AB,S.BB,B,al,BnB,Bal,CB,Cal, THETA.A,THETA.B, SAA, SAB, SBB,callrate,thresh)

mdist <- res
rlmm <- c(as.character(snpid),res)
write.table(t(rlmm),genotypefile,quote=F,sep=" ",row.names=F,col.names=F,append=T)
status <- paste("Processed",snpid)
print(status)
} #end if NOT discontinue

if (discontinue==1)	#when this SNP cannot be genotyped
{
res <- c("NC","NA")
res <- rep(res,NUMSAMP)
rlmm <- c(as.character(snpid),res)
write.table(t(rlmm),genotypefile,quote=F,sep=" ",row.names=F,col.names=F,append=T)
} #end if discontinue

} # end r loop of all the HapMap SNPs
#END  - MAIN LOOP

}

attr(Classify, "source") <- "See documentation or e-mail Nusrat Rabbee at nrabbee@post.harvard.edu"
