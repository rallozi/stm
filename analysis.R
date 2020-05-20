#####################################################
####Combining tab delimited text files by dataset####
#####################################################
dat.list<-function(mypath) {
  all.files<-list.files(mypath)
  fmri.dat<-list()
  for(i.file in 1:length(all.files)) {
    myfile<-all.files[i.file]
    fmri<-read.table(file.path(mypath, myfile), sep="\t", header=TRUE, comment.char="^")
    fmri.dat[[i.file]]<-fmri
  }
  names(fmri.dat)<-all.files
  return(fmri.dat)
}

sites<-c("Caltech", 
         "CMU", 
         "KKI", 
         "Leuven 1",
         "Leuven 2",
         "MaxMun",
         "NYU" , 
         "OHSU",
         "Olin", 
         "Pitt", 
         "SBL", 
         "SDSU", 
         "Stanford", 
         "Trinity", 
         "UCLA 1",
         "UCLA 2", 
         "UM 1",
         "UM 2",
         "USM", 
         "YALE")

#abide.source.path is the path to the folder that contains folders by dataset of downloaded ABIDE rsfMRI data 
#dat.path is the path that the Rdata containing fMRI data by site outputs to
for(isite in sites) {
  print(isite)
  fmri<-dat.list(file.path(abide.source.path, paste0(isite)))
  print(length(fmri))
  save(fmri, file=file.path(dat.path, paste0(isite, ".Rdata")))
}


###################################
####Preparing data for analysis####
###################################
#analysis.path contains necessary functions
source(file.path(analysis.path, "formX.R"))
source(file.path(analysis.path, "estK.R"))
source(file.path(analysis.path, "estFD.R"))
source(file.path(analysis.path, "bessel.moment.R"))

library(corpcor)
library(Matrix)

sites.p<-c("CALTECH", 
           "CMU", 
           "KKI",
           "LEUVEN_1",
           "LEUVEN_2",
           "MAX_MUN",
           "NYU" ,
           "OHSU",
           "OLIN", 
           "PITT", 
           "SBL", 
           "SDSU", 
           "STANFORD", 
           "TRINITY", 
           "UCLA_1", 
           "UCLA_2", 
           "UM_1", 
           "UM_2", 
           "USM", 
           "YALE")

sites.fmri<-c("Caltech", 
             "CMU", 
             "KKI", 
             "Leuven 1",
             "Leuven 2",
             "MaxMun",
             "NYU" , 
             "OHSU",
             "Olin", 
             "Pitt", 
             "SBL", 
             "SDSU", 
             "Stanford", 
             "Trinity", 
             "UCLA 1",
             "UCLA 2", 
             "UM 1",
             "UM 2",
             "USM", 
             "YALE")
sites<-data.frame(cbind(sites.p, sites.fmri), stringsAsFactors=FALSE)

all.phenotype<-c()
for(isite in sites.p) {
  #import phenotypic data
  X.01<-read.csv(file=file.path(abide.source.path, paste0("phenotypic_", isite, ".csv")),
                 stringsAsFactors=FALSE)
  all.phenotype<-rbind(all.phenotype, X.01)
}

#include all individuals in age range 
age.max<-14
age.min<-7
child.pheno<-all.phenotype[which(all.phenotype$AGE_AT_SCAN<age.max & all.phenotype$AGE_AT_SCAN>age.min),]
nrow(child.pheno) #486

mysites.01<-unique(child.pheno$SITE_ID)
mysites.02<-table(child.pheno$SITE_ID)
mysites.03<-names(mysites.02[mysites.02>10])
mysites.03<-mysites.03[!mysites.03=="LEUVEN_2"] #does not have FIQ

##number of time points per site
sites.fmri<-sites[which(sites$sites.p %in% mysites.03),"sites.fmri"]

ntp<-c()
for(isite in sites.fmri) {
  #import fMRI data
  load(file=file.path(dat.path, paste0(isite, ".Rdata")))
  itp<-unique(lapply(fmri, nrow))[[1]]
  ntp<-c(ntp, itp)
}
names(ntp)<-sites.fmri
sort(ntp)

train.sites<-names(ntp[which(ntp>=146)])

#remove individuals that have less than 146
train.fmri.01<-list()
for(isite in train.sites) {
  #import fMRI data created in fmri.data.prep.R
  load(file=file.path(dat.path, paste0(isite, ".Rdata")))
  train.fmri.01<-c(train.fmri.01, fmri)
}
names(train.fmri.01)<-unlist(strsplit(names(train.fmri.01), split="_rois_ho.1D"))

#rename fmri data using subject ids
get.ids.01<-unlist(strsplit(names(train.fmri.01), split="_rois_ho.1D"))
get.ids<-sapply(strsplit(get.ids.01, "\\CMU_a_00|\\CMU_b_00|\\KKI_00|\\Leuven_1_00|\\Leuven_2_00|\\NYU_00|\\Olin_00|\\Pitt_00|\\SDSU_00|\\Stanford_00|\\UCLA_1_00|\\UCLA_2_00|\\UM_1_00|\\UM_2_00|\\USM_00|\\Yale_00|\\SBL_00|\\Trinity_00| "), "[[", 2)
names(train.fmri.01)<-get.ids

rm.ids<-names(which(unlist(lapply(train.fmri.01, nrow))<146))
train.fmri.02<-train.fmri.01[-which(names(train.fmri.01) %in% rm.ids)]

##column 83 is all 0s for some and this region does not match a column description
check.83<-lapply(train.fmri.02, function(i) { #no 83 in region description
  i[,83]
})

train.fmri.03<-lapply(train.fmri.02, function(i) {
  i[-83]
})

train.sites.p<-sites[which(sites$sites.fmri %in% train.sites),"sites.p"]

train.phenotype<-c()
for(isite in train.sites.p) {
  #import phenotypic data
  X.01<-read.csv(file=file.path(abide.source.path, paste0("phenotypic_", isite, ".csv")),
                 stringsAsFactors=FALSE)
  train.phenotype<-rbind(train.phenotype, X.01)
}
train.phenotype<-train.phenotype[which(train.phenotype$AGE_AT_SCAN<age.max & train.phenotype$AGE_AT_SCAN>age.min),]
train.phenotype<-train.phenotype[-which(train.phenotype$SUB_ID %in% as.numeric(rm.ids)),]
nrow(train.phenotype)

#remove missing FIQ
any(train.phenotype$FIQ==-9999)
train.phenotype<-train.phenotype[-which(train.phenotype$FIQ==-9999),]

#data preparation
train.phenotype$DX_GROUP_CAT<-ifelse(train.phenotype$DX_GROUP==1,
                                     "ASD",
                                     "Control")
train.phenotype$SEX_CAT<-ifelse(train.phenotype$SEX==1,
                                "Male",
                                "Female")

train.phenotype$SITE_ID<-gsub(pattern="_", replacement=" ", x=train.phenotype$SITE_ID)

##function for table of descriptive statistics
cont.desc<-function(y, ylab, pheno.desc) {
  pheno.desc$y<-pheno.desc[,y]
  mymean<-aggregate(y~SITE_ID, data=pheno.desc, mean)
  colnames(mymean)[2]<-"Mean"
  
  mysd<-aggregate(y~SITE_ID, data=pheno.desc, sd)
  colnames(mysd)[2]<-"SD"
  
  out<-merge(mymean, mysd, by="SITE_ID")  
  out[,ylab]<-paste0(trimws(format(round(out$Mean, 1), nsmall=1)),
                     " (", 
                     trimws(format(round(out$SD, 1), nsmall=1)),
                     ")")
  
  return(out[,c("SITE_ID", ylab)])
}

cat.desc<-function(y, pheno.desc) {
  tab.n<-table(pheno.desc$SITE_ID, pheno.desc[,y])
  tab.perc<-trimws(format(round(prop.table(table(pheno.desc$SITE_ID, pheno.desc[,y]), 1)*100, 1), nsmall=1))
  
  out<-data.frame(SITE_ID=rownames(tab.n),
                  COL1=paste0(tab.n[,1], " (", tab.perc[,1], ")"),
                  COL2=paste0(tab.n[,2], " (", tab.perc[,2], ")"))
  colnames(out)[c(2,3)]<-paste0(colnames(tab.n))
  return(out)
}

pheno.desc<-train.phenotype
pheno.desc$SITE_ID<-"TOTAL"
pheno.desc<-rbind(pheno.desc, train.phenotype)

###########################################
####Table 3: ASD Descriptive Statistics####
###########################################
asd.desc<-pheno.desc[which(pheno.desc$DX_GROUP_CAT=="ASD"),]
age.desc<-cont.desc(y="AGE_AT_SCAN", ylab="Age, Mean (SD)", pheno.desc=asd.desc)
fiq.desc<-cont.desc(y="FIQ", ylab="Full IQ, Mean (SD)", pheno.desc=asd.desc)
sex.desc<-cat.desc(y="SEX_CAT", pheno.desc=asd.desc)
samp.size<-data.frame(table(asd.desc$SITE_ID), stringsAsFactors=FALSE)
colnames(samp.size)<-c("SITE_ID", "N")

tab.01<-merge(samp.size, sex.desc, by="SITE_ID", all.x=TRUE, all.y=TRUE)
tab.02<-merge(tab.01, age.desc, by="SITE_ID", all.x=TRUE, all.y=TRUE)
tab.03<-merge(tab.02, fiq.desc, by="SITE_ID", all.x=TRUE, all.y=TRUE)

rn.site<-tab.03[which(tab.03$SITE_ID!="TOTAL"),]
rn.site$SITE_ID<-paste0(seq(1:length(unique(train.phenotype$SITE_ID))),".", rn.site$SITE_ID)

tab<-rbind(rn.site, tab.03[which(tab.03$SITE_ID=="TOTAL"),])
rownames(tab)<-NULL
colnames(tab)[c(1)]<-c("Site")


################################################
####Table 4: Controls Descriptive Statistics####
################################################
norm.desc<-pheno.desc[which(pheno.desc$DX_GROUP_CAT=="Control"),]
age.desc<-cont.desc(y="AGE_AT_SCAN", ylab="Age, Mean (SD)", pheno.desc=norm.desc)
fiq.desc<-cont.desc(y="FIQ", ylab="Full IQ, Mean (SD)", pheno.desc=norm.desc)
sex.desc<-cat.desc(y="SEX_CAT", pheno.desc=norm.desc)
samp.size<-data.frame(table(norm.desc$SITE_ID), stringsAsFactors=FALSE)
colnames(samp.size)<-c("SITE_ID", "N")

tab.01<-merge(samp.size, sex.desc, by="SITE_ID", all.x=TRUE, all.y=TRUE)
tab.02<-merge(tab.01, age.desc, by="SITE_ID", all.x=TRUE, all.y=TRUE)
tab.03<-merge(tab.02, fiq.desc, by="SITE_ID", all.x=TRUE, all.y=TRUE)

rn.site<-tab.03[which(tab.03$SITE_ID!="TOTAL"),]
rn.site$SITE_ID<-paste0(seq(1:length(unique(train.phenotype$SITE_ID))),".", rn.site$SITE_ID)

tab<-rbind(rn.site, tab.03[which(tab.03$SITE_ID=="TOTAL"),])
rownames(tab)<-NULL
colnames(tab)[c(1)]<-c("Site")


#######################################
####Figure 1: Boxplot of age and IQ####
#######################################
library(ggpubr)
library(gridExtra)

age.p1.01<-ggboxplot(train.phenotype[which(train.phenotype$SITE_ID %in% unique(train.phenotype$SITE_ID)),],
                     x = "SITE_ID" ,
                     y =  "AGE_AT_SCAN"  ,
                     color = "SITE_ID" ,
                     palette ="jco",
                     add = "jitter",
                     ylim=c(age.min,age.max),
                     xlab="Dataset",
                     ylab="Age at Scan")
age.p1<-age.p1.01 + rremove("legend") + font("xlab", size = 24)+
  font("ylab", size = 24) + font("xy.text", size = 18)


fiq.p1.01<-ggboxplot(train.phenotype[which(train.phenotype$SITE_ID %in% unique(train.phenotype$SITE_ID)),],
                     x = "SITE_ID" ,
                     y =  "FIQ"   ,
                     color = "SITE_ID" ,
                     palette ="jco",
                     add = "jitter",
                     ylim=c(50,150),
                     xlab="Dataset",
                     ylab="Full IQ")
fiq.p1<-fiq.p1.01 + rremove("legend") + font("xlab", size = 24)+
  font("ylab", size = 24) + font("xy.text", size = 18)

grid.arrange(age.p1, fiq.p1, ncol=2, nrow=1)

#subset fmri to age group studies
train.fmri.04<-train.fmri.03[which(names(train.fmri.03) %in% train.phenotype$SUB_ID)]

#Restrict data to minimum number of time points
tps<-unique(unlist(lapply(train.fmri.04, function(x) {
  nrow(x)
})))
tps
tp<-min(tps)

train.fmri<-lapply(train.fmri.04, function(x) { #tp(x)n
  x[1:tp,]
})

n<-unique(unlist(lapply(train.fmri, function(x) {
  ncol(x)
})))
n

#subset phenotypic data
X.01<-train.phenotype
X<-X.01[,c("AGE_AT_SCAN", "FIQ", "DX_GROUP", "SEX")]
rownames(X)<-X.01$SUB_ID

#make sex 0,1
X$SEX<-ifelse(X$SEX==2, 1, 0)

#make autism status 0,1
X$DX_GROUP<-ifelse(X$DX_GROUP==2, 0, 1)

#standardize FIQ
X[,"FIQ"]<-(X[,"FIQ"]-mean(X[,"FIQ"]))/sd(X[,"FIQ"])

#datasets for ASD, CONTROL
X.ASD<-X[which(X$DX_GROUP==1),c("AGE_AT_SCAN", "FIQ", "SEX")]
X.CONT<-X[which(X$DX_GROUP==0),c("AGE_AT_SCAN", "FIQ", "SEX")]

fmri.ASD<-train.fmri[names(train.fmri) %in% rownames(X.ASD)]
fmri.CONT<-train.fmri[names(train.fmri) %in% rownames(X.CONT)]
  

################################################
####Figure 2: Temporal Autocorrelation Plots####
################################################
set.seed(22)
ind<-sort(sample(1:110, 3))
all.acf<-lapply(1:length(train.fmri), function(i) {
  acf.i<-list()
  for(j in 1:length(ind)) {
    acfobj<-acf(ts(train.fmri[[i]][,ind[j]]), plot=FALSE)
    acf.i[[j]]<-acfobj$acf[-1]
  }
  return(acf.i)
})

ind.acf<-lapply(1:length(ind), function(j) {
  acf.j<-sapply(all.acf, "[[", j)
  apply(acf.j, 1, mean)
})

plot.acf <- function(ACFest, region) {
  rr <- ACFest
  kk <- length(rr)
  nn <- tp
  plot(seq(kk),rr,type="h",lwd=1,yaxs="i",xaxs="i",
       ylim=c(floor(min(rr)),1),xlim=c(0,kk+1),
       xlab="Lag",
       ylab="Correlation",
       las=1, 
       cex.lab=1.5,
       cex.main=1.7,
       cex.axis=1.3,
       main=paste0("Region ", region))
  abline(h=0)
}


par(mfrow=c(1,3))
for(j in 1:length(ind.acf)) {
  plot.acf(ind.acf[[j]], region=ind[j])
}
par(mfrow=c(1,1))


###########################################
###Figure 3: Spatial Weighting Functions###
###########################################
library(minpack.lm)
library(gstat)
semivar<-function(FD.dat) {
  FD<-estFD(dat=FD.dat)
  
  #  Get array indices
  ind <- which(upper.tri(FD), arr.ind=TRUE) 
  #  cbind indices to values
  FD.df <- data.frame(cbind(FD[upper.tri(FD)], ind))
  colnames(FD.df)<-c("Distance", "R1", "R2")
  FD.bins<-seq(min(FD.df$Distance),max(FD.df$Distance),length.out=20)
  
  
  spv.01<-lapply(FD.dat, function(i) {

    all.v<-c()
    all.n<-c()
    for(bin in 1:(length(FD.bins)-1)) {
      mind<-FD.bins[bin]
      maxd<-FD.bins[bin+1]
      
      grpd<-FD.df[which(FD.df$Distance>=mind & FD.df$Distance<maxd),]  
      
      if(nrow(grpd)==0) {
        all.v<-c(all.v, NA)
      } else {
        bin.v<-c()
        for(l in 1:nrow(grpd)) {
          l.d<-(1/(2*nrow(i)))*(t(matrix(i[,grpd[l,"R1"]]-i[,grpd[l,"R2"]])) %*% matrix(i[,grpd[l,"R1"]]-i[,grpd[l,"R2"]]))
          bin.v<-c(bin.v, l.d)
        }
        all.v<-c(all.v,mean(bin.v))
      }
      all.n<-c(all.n,nrow(grpd))
    }
    return(list(all.v=all.v, n=all.n))
  })
  spv.02<-sapply(spv.01, "[[", 1)
  spv.03<-apply(spv.02,1,mean)
  spv<-data.frame(x=apply(cbind(FD.bins[1:(length(FD.bins)-1)], FD.bins[2:(length(FD.bins))]), 1, mean),
                  y=1-spv.03/max(spv.03, na.rm=T))
  
  spv.cc<-spv[complete.cases(spv),]
  myfit<-nlsLM(y~(1-x^(p))^(p),data=spv.cc,start=list(p=4))
  powcoef<-coef(myfit)

  return(list(spv=spv, FD=FD, powcoef=powcoef, FD.bins=FD.bins))
}
asd.semivar<-semivar(FD.dat=fmri.ASD)
cont.semivar<-semivar(FD.dat=fmri.CONT)

par(mar=c(5.1, 4.9, 5.1, 2.1))
plot(asd.semivar$spv$x,
     asd.semivar$spv$y,
     main="Empirical and Fitted Spatial Weighting Function by Disease Status",
     pch=20,
     lty=1,
     cex=2,
     cex.lab=1.5,
     cex.main=1.5,
     cex.axis=1.2,
     ylim=c(0,1),
     xlim=c(0,1),
     ylab=expression(paste("1-",hat(gamma),"(d)")),
     xlab="Functional Distance (d)",
     col="#0073C2FF",
     xaxt="n")
d<-seq(0,1,by=0.01)
w3<-(1-d^asd.semivar$powcoef)^asd.semivar$powcoef
lines(d, w3, type="l", col="#0073C2FF", lty=1)

points(cont.semivar$spv$x,
       cont.semivar$spv$y,
       pch=18,
       cex=2,
       col="#660000")
d<-seq(0,1,by=0.01)
w3<-(1-d^cont.semivar$powcoef)^cont.semivar$powcoef
lines(d, w3, type="l", col="#660000", lty=2)

axis(1, at=seq(0, 1, by=0.1), las=2)
legend("topright", 
       bty="n", 
       legend=c("ASD Empirical Estimate", 
                expression(paste("ASD Fitted: (1-",d^{2.58},")"^{2.58})), 
                "Control Empirical Estimate",
                expression(paste("Control Fitted: (1-",d^{2.39},")"^{2.39}))), 
       lty=c(NA,1, NA, 2),
       pch=c(20,NA, 18, NA),
       col=c("#0073C2FF", "#0073C2FF", "#660000", "#660000"),
       cex=1.25)

FD.ASD<-asd.semivar$FD
FD.CONT<-cont.semivar$FD

A.ASD<-(1-FD.ASD^asd.semivar$powcoef)^asd.semivar$powcoef
A.CONT<-(1-FD.CONT^cont.semivar$powcoef)^cont.semivar$powcoef
  
#Get K for each individual
#Xi contains a list of length m containing n(x)p covariate matrices
Xi.ASD<-formX(X=X.ASD, n=n, include.int=TRUE)
Xi.CONT<-formX(X=X.CONT, n=n, include.int=TRUE)

est.r<-function(Xi, A) {
  r<-c()
  for(icp in seq(0.1, 0.9, by=0.025)) {
    Ki<-estK(Xmat=Xi, Amat=A, cp=icp)
    if(class(Ki[[1]])=="numeric") {
      ir<-1
    } else {
      ir<-ncol(Ki[[1]])
    }
    r<-c(r,ir)  
  }
  return(r)
}
ASD.r<-est.r(Xi=Xi.ASD,   A=A.ASD)
CONT.r<-est.r(Xi=Xi.CONT, A=A.CONT)


############################################
###Figure 4: Singular Value Decomposition###
############################################
par(mar=c(5.1, 4.9, 4.1, 2.1))
main.text1<-expression(paste("Dimension Reduction from SVD of (I-X(X'X",")"^{-1},"X')A(I-X(X'X",")"^{-1},"X')"))
plot(seq(0.1, 0.9, by=0.025), ASD.r,
     type="b",
     xlab="Variance Explained",
     ylab="Reduced Dimension",
     col="#0073C2FF",
     xaxt="n",
     cex=2,
     cex.lab=2.5,
     cex.axis=1.75,
     pch=20)
axis(side=1, at=seq(0.1, 0.9, by=0.1), cex.axis=1.75)
mtext(main.text1,side=3, line=2,cex=2)

lines(seq(0.1, 0.9, by=0.025), CONT.r,
      type="b",
      xaxt="n",
      col="#660000",
      pch=18,
      cex=2)
legend("topleft", 
       bty="n", 
       legend=c("ASD", 
                "Control"), 
       pch=c(20, 18),
       col=c("#0073C2FF", "#660000"),
       cex=2)

fmri.ASD<-lapply(fmri.ASD, function(i) {
  as.matrix(i)
})

fmri.CONT<-lapply(fmri.CONT, function(i) {
  as.matrix(i)
})

bes.ASD<-bessel.moment(FD=FD.ASD)
bes.CONT<-bessel.moment(FD=FD.CONT)

for(icp in  seq(0.1,0.8,by=0.05)) {
  Ki.ASD<-estK(Xmat=Xi.ASD, Amat=A.ASD, cp=icp)
  Ki.CONT<-estK(Xmat=Xi.CONT, Amat=A.CONT, cp=icp)
  
  save(Xi.ASD, Ki.ASD, A.ASD, FD.ASD, fmri.ASD, bes.ASD,
       file=paste0("ASD.stmdat.cp=",icp, ".Rdata"))
  
  save(Xi.CONT, Ki.CONT, A.CONT, FD.CONT, fmri.CONT, bes.CONT,
       file=paste0("CONT.stmdat.cp=",icp, ".Rdata"))
}


####################################################################
####Spatiotemporal Modeling and Functional Connectivity Analysis####
####################################################################
#load packages 
library(MASS)
library(mvtnorm)

#load STM codes, stm.path is path to folder containg STM codes
setwd(stm.path)
file.sources = list.files(pattern="*.R")
sapply(file.sources,source,.GlobalEnv)

#####################################################################################
load(file.path(dat.path, "ASD.stmdat.cp=0.7.Rdata"))

#Initialize parameters
p<-ncol(Xi.ASD[[1]])
r<-ncol(Ki.ASD[[1]])
m<-length(fmri.ASD)
tp<-nrow(fmri.ASD[[1]])
SigmaEtaForm<-"Diagonal"

phi0<-get.phi0(z=fmri.ASD, 
               X=Xi.ASD, 
               K=Ki.ASD, 
               SigmaEtaForm=SigmaEtaForm, 
               bessel.moments=bes.ASD)

mod.it.ASD  <- Stem.Model(z=fmri.ASD,            #a list of length m containing observation matrices of dimension tp*n
                          covariates=Xi.ASD,     #a list of length m containing covariate matrices of dimension n*p
                          funcdist=FD.ASD,       #n*n matrix of functional distance
                          phi=phi0,              #starting parameters
                          K=Ki.ASD,              #a list of length m containing loading matrices of dimension n*r
                          r=r,                   #r, dimension of the latent process
                          SigmaEta=SigmaEtaForm)

my.est.ASD <- Stem.Estimation(StemModel=mod.it.ASD, 
                              precision=1e-10, 
                              max.iter=5000, 
                              SigmaEta=SigmaEtaForm, 
                              learning.rate=0.1)

betahat<-my.est.ASD$estimates$phi.hat$beta
Ghat<-matrix(my.est.ASD$estimates$phi.hat$G, r, r)
y<-my.est.ASD$estimates$yt.tmin1
X<-Xi.ASD
K<-Ki.ASD

zhat.ASD<-lapply(1:m, function(i) {
  myy<-y[[i]] #time 1 to tp
  
  zhat.t<-lapply(1:nrow(myy), function(tt) {
    if(tt==1) {
      y.tm1<-my.est.ASD$estimates$phi.hat$m0
    } else {
      y.tm1<-myy[tt-1,]
    }
    
    X[[i]] %*% betahat + K[[i]] %*% Ghat %*% y.tm1
  })
  t(do.call(cbind, zhat.t))
})

library(Matrix)
library(corpcor)
library(whitening)
spcor.i.ASD<-lapply(1:m, function(i) {
  covt<-cov(t(zhat.ASD[[i]])) #covariance between times
  
  if(!is.positive.definite(covt)) {
    covt<-as.matrix(nearPD(covt)$mat)
  }
  
  W<-whiteningMatrix(covt, method="Cholesky")
  Zstar<-tcrossprod(t(zhat[[i]]), W)
  cor(t(Zstar))
  
})

fish.asd.i<-lapply(spcor.i.ASD, function(i) {
  0.5*log((1+i)/(1-i))
})

#fishers z transformation
fish.asd<-sumMatrices(fish.asd.i)/m
colnames(fish.asd)<-rownames(fish.asd)<-colnames(fmri.ASD[[1]])
m.asd<-m


#####################################################################################
load(file.path(dat.path, "CONT.stmdat.cp=0.7.Rdata"))

#Initialize parameters
p<-ncol(Xi.CONT[[1]])
r<-ncol(Ki.CONT[[1]])
m<-length(fmri.CONT)
tp<-nrow(fmri.CONT[[1]])
SigmaEtaForm<-"Diagonal"

phi0<-get.phi0(z=fmri.CONT, 
               X=Xi.CONT, 
               K=Ki.CONT, 
               SigmaEtaForm=SigmaEtaForm, 
               bessel.moments=bes.CONT)


mod.it.CONT  <- Stem.Model(z=fmri.CONT,            #a list of length m containing observation matrices of dimension tp*n
                           covariates=Xi.CONT,     #a list of length m containing covariate matrices of dimension n*p
                           funcdist=FD.CONT,       #n*n matrix of functional distance
                           phi=phi0,              #starting parameters
                           K=Ki.CONT,              #a list of length m containing loading matrices of dimension n*r
                           r=r,                   #r, dimension of the latent process
                           SigmaEta=SigmaEtaForm)

my.est.CONT <- Stem.Estimation(StemModel=mod.it.CONT, 
                               precision=1e-10, 
                               max.iter=5000, 
                               SigmaEta=SigmaEtaForm, 
                               learning.rate=0.1)

betahat<-my.est.CONT$estimates$phi.hat$beta
Ghat<-matrix(my.est.CONT$estimates$phi.hat$G, r, r)
y<-my.est.CONT$estimates$yt.tmin1
X<-Xi.CONT
K<-Ki.CONT

zhat.CONT<-lapply(1:m, function(i) {
  myy<-y[[i]] #time 1 to tp
  
  zhat.t<-lapply(1:nrow(myy), function(tt) {
    if(tt==1) {
      y.tm1<-my.est.CONT$estimates$phi.hat$m0
    } else {
      y.tm1<-myy[tt-1,]
    }
    
    X[[i]] %*% betahat + K[[i]] %*% Ghat %*% y.tm1
  })
  t(do.call(cbind, zhat.t))
})

spcor.i.CONT<-lapply(1:m, function(i) {
  covt<-cov(t(zhat.CONT[[i]])) #covariance between times
  
  if(!is.positive.definite(covt)) {
    covt<-as.matrix(nearPD(covt)$mat)
  }
  
  W<-whiteningMatrix(covt, method="Cholesky")
  Zstar<-tcrossprod(t(zhat[[i]]), W)
  cor(t(Zstar))
  
})

fish.cont.i<-lapply(spcor.i.CONT, function(i) {
  0.5*log((1+i)/(1-i))
})

#fishers z transformation
fish.cont<-sumMatrices(fish.cont.i)/m
colnames(fish.cont)<-rownames(fish.cont)<-colnames(fmri.CONT[[1]])
m.cont<-m


#####################################################################################
mycols<-read.csv(file.path(dat.path, "fmricolnames.csv"), stringsAsFactors=FALSE)

prep<-function(res, label) {
  res[upper.tri(res)]<-NA
  diag(res)<-NA
  
  out.01<-data.frame(cbind(which(!is.na(res),arr.ind = TRUE),na.omit(as.vector(res))))
  labs<-data.frame(column=1:ncol(res), Region=sapply(strsplit(colnames(res), split="X."), "[[", 2), stringsAsFactors=FALSE)
  out.02<-merge(out.01, labs, by.x="row", by.y="column")
  colnames(out.02)[which(colnames(out.02)=="Region")]<-"Region1"
  out.03<-merge(out.02, labs, by.x="col", by.y="column")
  colnames(out.03)[which(colnames(out.03)=="Region")]<-"Region2"
  
  out.04<-out.03[,c("Region1", "Region2", "V3")]
  colnames(out.04)[which(colnames(out.04)=="V3")]<-paste0(label, ".Value")
  
  r1.desc<-r2.desc<-mycols[,-2]
  colnames(r1.desc)<-c("Region1", "Region1.Description", "Region1.Number", "Region1.DMN")
  colnames(r2.desc)<-c("Region2", "Region2.Description", "Region2.Number", "Region2.DMN")
  
  out.05<-merge(x=out.04, y=r1.desc, by="Region1")
  out<-merge(x=out.05, y=r2.desc, by="Region2")
  out$Link<-paste0(out[,"Region1.Number"],"_", out[,"Region2.Number"])
  
  return(out)
}

dat.asd<-prep(res=fish.asd, label="ASD")
dat.cont<-prep(res=fish.cont, label="CONT")
dat<-merge(dat.asd, dat.cont, by=intersect(names(dat.asd), names(dat.cont)))

dat$Zstat<-unlist(lapply(1:nrow(dat), function(l) {
  z1<-dat[l,"ASD.Value"]
  z0<-dat[l,"CONT.Value"]
  
  zstat<-(z1-z0)*(1/(m.asd*(tp-3))+1/(m.cont*(tp-3)))^(-1/2)
  return(zstat)
}))

dat$unadjp<-2*(1-pnorm(abs(dat$Zstat)))

library(locfdr)
w <- locfdr(dat$Zstat, nulltype=2)
dat$LOCFDR<-w$fdr


############################################
###Table 5: Significant Links at FDR 0.1###
############################################
sig.link.01<-dat[which(dat$LOCFDR<0.1),]
sig.link<-sig.link.01[order(sig.link.01$LOCFDR),]
rownames(sig.link)<-NULL
sig.link$Order<-1:nrow(sig.link)
tab<-sig.link[,c("Region1.Description" , "Region2.Description", "DMN")]


#####################################
###Figure 5: Bessel plots by group###
#####################################
library(SpatialExtremes)
par(mar=c(5.1, 6.1, 4.1, 2.1))
covariance(nugget = 0,sill = 1, range = s1, smooth = nu1,
           cov.mod = "bessel", xlim = c(0,1), ylim = c(-0.25, 1),
           col="#0073C2FF", 
           lty=1,
           cex.lab=2.25,
           cex.main=2.25,
           cex.axis=1.24,
           xlab="Distance (d)", ylab=bquote(rho~"(d)"),
           main=expression(paste("Estimated Bessel Functions for ASD Subjects and Controls")))
covariance(nugget = 0, sill = 1, range = s0, smooth = nu0, lty=2,
           xlim = c(0,1), cov.mod = "bessel", add = TRUE, col="#660000")
abline(h=0, lty=3)
legend("topright", c(expression(paste("ASD: ", nu, "=1.517, ", s,"=0.036")), 
                     expression(paste("Control: ", nu, "=0.837, ", s,"=0.038"))),
       col = c("#0073C2FF", "#660000"), lty = c(1,2), bty="n",cex=2)



