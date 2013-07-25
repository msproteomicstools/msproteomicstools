library(lattice)

extract_chrom <- function(t, thisdata, productmz, extraction_window=0.05)
{
  this_spectrum = subset(thisdata, SEC == t)
  return(sum(subset(this_spectrum, MZ > productmz-(extraction_window/2) & MZ < productmz+(extraction_window/2))$INT))
}

graphme <- function(xxp,allmx){
  xxp <- xxp[length(xxp):1]
  allmx <- allmx[allmx$MZ > 400,]
  sum(is.element(allmx$label,xxp))
  allmx <- allmx[is.element(allmx$label,xxp),]
  print(dim(allmx))
  allmx$MZ <- as.factor(allmx$MZ)
  return(allmx)
}

irt2rt <- function(x,c=2148.68,m=33.87) {
  return(m*x+c)
}

plotgraph <- function(peptide,background,swath,productmz,irt) {
  txtfiles <- dir(patter=glob2rx(paste("*",background,"*_",swath,"_filtered.dta2d",sep="")))

  rawdata <- list()
  for(i in 1:length(txtfiles))
  {
    rawdata[[i]] <- read.csv(txtfiles[i], sep="\t")
    names(rawdata[[i]])<-c("SEC","MZ","INT")
  }

  data <- list()
  for(i in 1:length(txtfiles))
  {
    df<-data.frame()
    for(j in 1:length(productmz)) {
      dfj <- data.frame("INT" = sapply( unique(rawdata[[i]]$SEC), extract_chrom, thisdata=rawdata[[i]], productmz=productmz[j]), "SEC"=unique(rawdata[[i]]$SEC))
      dfj$MZ <- rep(productmz[j],dim(dfj)[1])
      df<-rbind(df,dfj)
    }
    
    data[[i]] = df
  }

  xx <- c("x512","x256","x128","x064","x032","x016","x008","x004","x002","x001")
  length(xx)
  allm <- NULL
  label <- NULL
  for(i in 1:10){
    allm <- rbind(allm,data[[i]])
    labelt <- rep(xx[i],dim(data[[i]])[1])
    label <- c(label, labelt)
  }
  allm <- cbind(label, allm)
  allm <- data.frame(as.factor(allm$label), as.numeric(allm$SEC), as.numeric(allm$MZ), as.numeric(allm$INT))

  colnames(allm) <- c("label","SEC","MZ","INT")
  colnames(allm)
  allm$label[1:10]

  xxs <- c("x512","x256","x128","x064","x032","x016","x008","x004","x002","x001")
  allmx <- allm

  pdf(file=paste(peptide,"_",background,".pdf",sep=""),width=6, height=length(xxs)*1.5)
  print(xyplot(INT ~  SEC | label ,data=allmx,type="l",xlim=c(irt2rt(irt)-100,irt2rt(irt)+100),scales=list(y=list(relation="free", cex=0.7,rot=45)),groups=MZ,layout=c(1,length(xxs)),xlab="RT [s]", ylab="INT",as.table=TRUE))
  dev.off()
}

setwd("/IMSB/users/georger/html/osw_peptides/data")

background<-list("water"="no_background","yeast"="yeast","human"="human")
assay_swath<-list("VGDTVLYGK"=3,"IADIQLEGLR"=6,"TGGDEFDEAIIK"=10,"LITVEGPDGAGK"=7,"LVDEEGNDVTPEK"=13)
assay_precursormz<-list("VGDTVLYGK"=480.268,"IADIQLEGLR"=569.329,"TGGDEFDEAIIK"=651.819,"LITVEGPDGAGK"=582.821,"LVDEEGNDVTPEK"=726.851)
assay_productmz<-list("VGDTVLYGK"=c(688.412,375.212,803.439,860.46),"IADIQLEGLR"=c(597.359,838.502,725.418),"TGGDEFDEAIIK"=c(972.513,1144.56,696.402,843.47),"LITVEGPDGAGK"=c(609.308,938.467,738.351,837.419),"LVDEEGNDVTPEK"=c(1125.51,996.472,1240.54,867.43))
assay_irt<-list("VGDTVLYGK"=3.7,"IADIQLEGLR"=49.4,"TGGDEFDEAIIK"=40.8,"LITVEGPDGAGK"=10.9,"LVDEEGNDVTPEK"=-5.1)

for(i in 1:length(assay_swath)) {
  for(j in 1:length(background)) {
    plotgraph(names(assay_swath)[[i]],background[[j]],assay_swath[[i]],assay_productmz[[i]],assay_irt[[i]])
  }
}
