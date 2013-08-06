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

plotgraph <- function(assay_irt,background,rt_extraction_window=180) {
  txtfiles <- dir(pattern=glob2rx(paste("*",background,"*","._chrom.mzML.dta2d",sep="")))

  rawdata <- list()
  for(i in 1:length(txtfiles))
  {
    rawdata[[i]] <- read.csv(txtfiles[i], sep="\t")
    names(rawdata[[i]])<-c("SEC","MZ","INT")
  }

  # use this code to extract chromatograms
  # data <- list()
  # for(i in 1:length(txtfiles))
  # {
  #   df<-data.frame()
  #   for(j in 1:length(productmz)) {
  #     dfj <- data.frame("INT" = sapply( unique(rawdata[[i]]$SEC), extract_chrom, thisdata=rawdata[[i]], productmz=productmz[j]), "SEC"=unique(rawdata[[i]]$SEC))
  #     dfj$MZ <- rep(productmz[j],dim(dfj)[1])
  #     df<-rbind(df,dfj)
  #   }
  #   data[[i]] = df
  # }

  data<-rawdata

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

  if (background=="human") {
    irt<-irt2rt(assay_irt[[1]],1687.64,33.61)
  }
  else if (background=="yeast") {
    irt<-irt2rt(assay_irt[[1]],2105.2,34.27)
  }
  else if (background=="no_background") {
    irt<-irt2rt(assay_irt[[1]],2150.32,35.05)
  }

  pdf(file=paste(names(assay_irt)[[1]],"_",background,".pdf",sep=""),width=6, height=length(xxs)*1.5)
  print(xyplot(INT ~  SEC | label ,data=subset(allmx,SEC >= irt-rt_extraction_window & SEC <= irt+rt_extraction_window),type="l",xlim=c(irt-rt_extraction_window,irt+rt_extraction_window),scales=list(y=list(relation="free", cex=0.7,rot=45)),groups=MZ,layout=c(1,length(xxs)),xlab="RT [s]", ylab="INT",as.table=TRUE))
  dev.off()
}

background<-list("water"="no_background","yeast"="yeast","human"="human")
assays<-list("VGDTVLYGK"=3.7,"IADIQLEGLR"=49.4,"TGGDEFDEAIIK"=40.8,"LITVEGPDGAGK"=10.9,"LVDEEGNDVTPEK"=-5.1)
assay_irt<-assays[tail(strsplit(getwd(),"/")[[1]],n=1)]

for(j in 1:length(background)) {
  plotgraph(assay_irt,background[[j]])
}
