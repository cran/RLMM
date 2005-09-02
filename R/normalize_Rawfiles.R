#--------------------------------------------------------------------#
#-- Creates *.norm files from *.raw files which                    --#
#-- has normalized intensity values for PMA and PMB using a pseudo --#
#-- normalization method by using the CQV vector                   --#
#-- By: Nusrat Rabbee & Gary Wong                                  --#
#--------------------------------------------------------------------#

normalize_Rawfiles <- function(cqvfile="", probefiledir = getwd()){  
   
   
   # CQV is a sorted column vector of normalized values
   # Xba.cqv is a text file
   
   print("Reading CQV")
   
   CQV <- read.table(cqvfile)
   CQV <- CQV[[1]]
   
   print("Reading CQV ...Done")
   
   if(length(CQV) == 0)
       stop("Cannot locate CQV within the directory")
   
   CQV <- sort(CQV)
    
# ------
# Read all file names in dir

   currentdir <- getwd()
   setwd(probefiledir)
   
   fn <- list.files()
   fn <- fn[grep(".raw",fn)]
   
   NUMSAMP <- length(fn) #number of *.raw files
   NUMSNP <- length(readLines(fn[1])) #take the length of the first file

   if(grep(".raw",fn[1])==F)
   stop("There are no *.raw files in this directory")
   
   fn <- data.frame(fn)
   fn <- t(fn)
   pma.idx <- seq(2,41,4) # index for 10 pma values
   pmb.idx <- seq(4,41,4) # index for 10 pmb values


   
   
       if((length(CQV)/20) != NUMSNP) {
          error <- paste("File",cqvfile,"is not the correct size")
          stop(error)
	}

   for (j in 1:length(fn)) # of *.raw 
   {
   if(length(readLines(fn[j])) != NUMSNP){
      stop(paste("# of SNPs is not equal to",NUMSNP,"in", fn[j]))    
      cat("The # of SNPs in the *.raw files is correct \n")
      }

   print("Processing File")
   print(as.character(fn[j]))
   
   dat <- read.table(fn[j],as.is=T)
   new <- dat[,c(pma.idx,pmb.idx)] # extract pma & pmb values
   CQVTEMP <- as.vector(t(new))
   NORM <- CQV[floor(rank(CQVTEMP))]
   
   Norm.RawFiles<-matrix(NORM, nrow=length(new[,1]),ncol=length(new[1,]),byrow=TRUE)
   Norm.RawFiles<-as.data.frame(Norm.RawFiles)
   Norm.RawFiles<-data.frame(cbind(dat[,1],Norm.RawFiles))

   write.table(Norm.RawFiles, file = paste(strsplit(fn[j],split=".raw"),"norm",sep="."), quote = FALSE, row.names=FALSE, col.names=FALSE)
   }  # end for loop
  
   setwd(currentdir)
}

attr(normalize_Rawfiles, "source") <- "See documentation or e-mail Nusrat Rabbee at nrabbee@post.harvard.edu"
