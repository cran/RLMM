create_Thetafile <- function(probefiledir = getwd(),start=1,end=-1,thetafile=""){

   ## This program creates allele summaries (theta_A and theta_B values) for
   ## each chip , each SNP. It uses Robust Linear Model with constraints 
   ## for the A probes and B probes separately, in order to produce the
   ## the theta_A and theta_B values
   ## --------- 
   ## Functions within create_Thetafile: RMA & write.theta
   ##
   ## BY: Nusrat Rabbee & Gary Wong
   #########################################################################


####### DEFINE THE FUNCTIONS #########
RMA <- function(A,B,flag=FALSE) 
{

 new.y.A <- as.vector(A)
 new.y.B <- as.vector(B)
 num.probes <- nrow(A)
 num.chips <- ncol(A)
 
 temp<- seq(1,num.chips,1)
 temp2<- rep(1,num.probes)
 C2 <- kronecker(temp,temp2)
 temp <- rep(1,num.chips)
 temp2 <- seq(1,num.probes,1)
 P2 <- kronecker(temp,temp2)
       
 aa.A <- rlm(new.y.A ~ as.factor(C2) + C(as.factor(P2),"contr.sum"),maxit=100)
 aa.B <- rlm(new.y.B ~ as.factor(C2) + C(as.factor(P2),"contr.sum"),maxit=100)
 int.A <- coefficients(aa.A)[1]
 int.B <- coefficients(aa.B)[1] 
 chips <- coefficients(aa.A)[2:(num.chips)]
 THETA.A <- c(int.A,int.A+chips)
 chips <- coefficients(aa.B)[2:(num.chips)]
 THETA.B <- c(int.B,int.B+chips)

 THETA.A <- as.vector(THETA.A)
 THETA.B <- as.vector(THETA.B)
 return(c(THETA.A,THETA.B))
}

write.theta <- function(babe,bits,num.samp,num.probes,thetafile)
{

snpid <- snpnames[babe]
print(snpid)
training.size <- length(bits)

idx <- babe

num.samp-> NUMSAMP
x40 <- matrix(20,NUMSAMP)
st <- (idx-1)*num.probes + 1
fn <- (idx-1)*num.probes + num.probes 
x40 <- PROBEDAT[st:fn,]			#extract the probe data; used to be PROBEDAT


y<- x40[1:num.probes,]
y.A <- y[1:10,bits]
y.B <- y[11:20,bits]

aa <- RMA(y.A,y.B)	#Obtain theta for all the samples

THETA.A <- aa[1:(training.size)]
THETA.B <- aa[(1+training.size):(2*training.size)]

THETA.A <- round(THETA.A,4)
THETA.B <- round(THETA.B,4)

res <- c(as.character(snpid),as.character(THETA.A), as.character(THETA.B))
write.table(t(res),thetafile,row.names=F,col.names=F,append=T,sep=" ",quote=F)
return(res)
} # end function

# ==================================================================================
# Main
library(MASS)

if(thetafile=="")
   stop("Please specify a name for the thetafile parameter")

currentdir <- getwd()
print(probefiledir)
list.files(probefiledir)
setwd(probefiledir)
dirFiles <- list.files()
normfiles <- dirFiles[grep(".norm",dirFiles)]
rawfiles <- dirFiles[grep(".raw",dirFiles)]
arawfile <- rawfiles[1]
fn <- read.table(rawfiles[1],as.is=T)
snpnames <- fn[,1]

NUMSAMP <- length(normfiles)
print(paste("Total number of .norm files found in probefiledir:",NUMSAMP))

if(length(normfiles)==F)
   stop("The correct number of *.norm files is not in the directory")

if(grep(".norm",normfiles[1])==F)
   stop("There are no *.norm files in this directory")

fn <- read.table(normfiles[1],as.is=T)
fn <- fn[,1]
fn <- matrix(data=fn)

#creating PROBEDAT

NUMSNP <- length(snpnames)

print(paste("NUMSNP is", NUMSNP))
print(paste("NUMSAMP is", NUMSAMP))

if(start<1 || start>NUMSNP)
   start <- 1
 
if(end==-1 || end>NUMSNP)
   end <- NUMSNP

if (start>end)
   stop("start value is greater than the end value, please change")

PROBEDAT <- matrix(-1,NUMSNP*20,NUMSAMP)
for (j in 1:NUMSAMP) # of normfiles
{
d <- paste(normfiles[j],sep="")
d <- as.character(normfiles[j])
print(paste("Processing ",d))
temp <- read.table(d,as.is=T)
PROBEDAT[,j] <- temp[,1]
}

babe <- seq(start,end,1)

tr.bits <- seq(1,NUMSAMP,1)	#Get the theta's for all the samples together
setwd(currentdir)
RES <- apply(t(babe),2,write.theta,bits=tr.bits,num.samp=NUMSAMP,num.probes=20,thetafile=thetafile)

rm(PROBEDAT)
#END WRITE_THETA function
}

attr(create_Thetafile, "source") <- "See documentation or e-mail Nusrat Rabbee at nrabbee@post.harvard.edu"
