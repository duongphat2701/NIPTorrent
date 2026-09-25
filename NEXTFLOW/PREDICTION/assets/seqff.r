##################################################################################################################################################
##################################################################################################################################################
##### This source is derived from the original work of:
#####
##### SeqFF.R
##### Sung Kim, Phd
##### Sequenom, Inc.
##### 
##################################################################################################################################################
##### Instructions
##### Command Line
##### R CMD BATCH
##### --i input directory
##### --f input file name
##### --d output directory
##### --o output file name
##### --t data type; sam file (without header) or tabulated read counts ordered by genomic coordinates found in SupplementalTable1.csv
##### SeqFF.R
###################################################################################################################################################
##### Modified by:
##### Marco Milanesio, PhD
##### marco.milanesio@univ.cotedazur.fr
#####
##### Part of the NiPTUNE framework
##### Usage: Rscript seqff.r target.file seqff.dir
##### Where:
##### target.file: a file composed by columns 3,4 of the original bam/sam (e.g., ['chr1 10001\n', ...]
##### seqff.dir: location of the original seqff (i.e., where *.Rdata and *.csv reside)
###################################################################################################################################################

library(MASS) #to write a matrix
suppressMessages(library(Rsamtools))
suppressMessages(library(stringi))

ff.pred <- function(gc.norm.bc.61927, B, mu, parameter.1, parameter.2){
  gc.norm.bc.61927[is.na(gc.norm.bc.61927)] <- 0
  gc.norm.bc <- gc.norm.bc.61927[grepl('chr[0-9]', names(gc.norm.bc.61927))]
  gc.norm.bc.c <- gc.norm.bc - mu
  y.hat <- matrix(c(1, gc.norm.bc.c), nrow = 1) %*% B
  y.hat.rep <- sum(y.hat, na.rm = T) / sum(gc.norm.bc)
  ff.hat <- (y.hat.rep + parameter.1) * parameter.2
  return(ff.hat)
  
}

##################################################################################################################################################
##################################################################################################################################################
##### COMMAND LINE ARGUMENTS
args <- commandArgs(trailingOnly = TRUE)
datatype = "sam"  # "counts"
file.name = args[1]
seqff.dir = args[2]
input.dir = dirname(file.name)
# output.filename = paste(file.name,"out",sep=".") 
output.filename = ""
if (!file.exists(file.name)) {
    stop(paste("Input file does not exist:", file.name))
}


load(file.path(seqff.dir, "SupplementalFile1.RData"))
bininfo = read.csv(file.path(seqff.dir, "SupplementalTable2.csv"))
colnames(bininfo)[1]="binName"
bininfo$binorder=c(1:61927)

setwd(input.dir)
##################################################################################################################################################
##### READ IN DATA

# dat <- read.table(pipe(paste("cut -f3,4",file.name,sep=" ")), header=FALSE,colClasses= c("character","integer"))
if( datatype =="sam" )
{
# dat <- read.table(pipe(paste("cut -f1,2",file.name,sep=" ")), header=FALSE,colClasses= c("character","integer"))
# colnames(dat)=c("refChr","begin")

dat <- as.data.frame(stri_list2matrix(scanBam(file.name)[[1]][c('rname','pos')]),stringsAsFactors = F)
tmp <- scanBam(file.name)[[1]][c('rname','pos')]
dat <- as.data.frame(matrix(NA, ncol=2,nrow=length(tmp$rname)))
names(dat) <- c("refChr","begin")
dat[,1] <- as.character(tmp$rname)
dat[,2] <- as.numeric(tmp$pos)
dat[,'begin'] <- as.numeric(dat[,'begin'])

dat=dat[dat$refChr!="*" & dat$refChr!="chrM" ,]
binindex = ifelse((dat$begin%%50000)==0,floor(dat$begin/50000)-1,floor(dat$begin/50000))
fastbincount = table(paste(dat$refChr,binindex, sep="_"))
newtemp=as.data.frame(matrix(NA, ncol=2, nrow=length(fastbincount)))
newtemp[,1]=paste(names(fastbincount))
newtemp[,2]=as.numeric(paste(fastbincount))
colnames(newtemp)=c("binName","counts")
rm(dat,binindex,fastbincount)

}
if( datatype =="counts" )
{
newtemp <- read.table(file.name, header=FALSE,colClasses= c("character","integer"))
colnames(newtemp)=c("binName","counts")
}

bininfo = merge(bininfo, newtemp, by="binName",all.x=T)
bininfo=bininfo[order(bininfo$binorder),]
##################################################################################################################################################
##### DATA PROCESSING
autosomebinsonly = bininfo$BinFilterFlag==1 & bininfo$CHR!="chrX" & bininfo$CHR!="chrY"
alluseablebins = bininfo$BinFilterFlag==1 
autoscaledtemp  <- bininfo$counts[autosomebinsonly]/sum(bininfo$counts[autosomebinsonly], na.rm=T)
allscaledtemp  <- bininfo$counts[alluseablebins]/sum(bininfo$counts[autosomebinsonly], na.rm=T)
# additive loess correction
mediancountpergc <- tapply(
autoscaledtemp,bininfo$GC[autosomebinsonly], function(x) median(x, na.rm=T))
## prediction 
loess.fitted  <- predict( loess(mediancountpergc ~ as.numeric(names(mediancountpergc))), bininfo$GC[alluseablebins]) 
normalizedbincount <- allscaledtemp + ( median(autoscaledtemp, na.rm=T) - loess.fitted )  

bincounts=rep(0,61927)
names(bincounts) = bininfo$binName
bincounts[alluseablebins] <- (normalizedbincount/sum(normalizedbincount, na.rm=T)) * length(normalizedbincount)

wrsc=ff.pred(bincounts,B,mu,parameter.1,parameter.2)
enet = bincounts %*% elnetbeta+elnetintercept
ff=c(((wrsc+enet)/2)*100, enet, wrsc)
names(ff)=c("seqff","Enet","WRSC")

setwd(input.dir)

# write.csv(ff, file=output.filename)
# cat(output.filename)
print(ff)



