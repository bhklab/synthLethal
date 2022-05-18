# Predicting the response to immunotherapy via SR partners
# Only publicly available dataset are shared

Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

library(data.table)
library(Rcpp)
library(survival)
library(survminer)
library(ROCR)
library(caTools)

# set up working directory
setwd("SELECT")

source("./R/source.r",echo=T)
load("./data/prob.TCGA.ICI.RData")
load("./data/sr.ICI.RData")
pancancer=prob

drug.predict=function(pval,flag){
	res=NULL
	if (sum(!is.na(pval))>=3){
		cx=getAUC(pval,flag)
		res$auc=cx[1]
	}
	return(res)
}
colSums2<-function(X){
	if (class(X)=="numeric") y=X
	if (class(X)=="matrix") y=colSums(X,na.rm=T)
	return(y)
}
sr.score <- function(dat,sr.final1){
		score1=rep(0,ncol(dat$mRNA))	
		partners1=match(sr.final1[,2],dat$genes)
		partners1=partners1[!is.na(partners1)]
		if (length(partners1)==1) 
					score1=1*(dat$mRNA.rank2[partners1,] <0.33)
		if (length(partners1) >1) 
					score1=colSums(dat$mRNA.rank2[partners1,] <0.33,na.rm=T)/length(partners1)
		score1=1-score1			
	return(score1)
}
sr.find<-function(sr.final1){
    sl.res=NA
    FDR=0.1
    flag=1                                              
	iy=which(apply(sr.final1[,c(5,3,4)],1,min)<FDR & apply(sr.final1[,c(8,6,7)],1,min)<FDR)  # step I: underrepresented
    if (length(iy)>1) {
        flag=2
    	sr.final1=sr.final1[iy,]
    	clinical.beta=sr.final1[,c(9,11,13,15)] 											# step II: survival
    	clinical.p=sr.final1[,c(10,12,14,16)]
    	clinical.fdr = sapply(1:ncol(clinical.beta), function(tt) {
    		p = ifelse(clinical.beta[,tt] < 0 , clinical.p[,tt], 1 )
    	})
        iz = which(rowSums(clinical.fdr < FDR ) >= 1)
        if (length(iz)>1) {
            flag=3
            sr.final1=sr.final1[iz,]
            ix=which(sr.final1[,17]< 2.315712e+01)                                            # step III: phylogenetic     
            sr.final1=sr.final1[ix,]
        }
    }
    if (flag==3) sl.res=sr.final1
    return(sl.res)
}

eval.km=function(sr.final1,dat){
	rmv=NULL
	sl.ess2=sr.final1
	logrank.p=eff=NA
	
	sr.x=cbind(prob$genes[sr.final1[,1]],prob$genes[sr.final1[,2]])
	score1=sr.score(dat,sr.x)
	dat$surv.dt=data.frame(time=dat$survival[,1],status=ifelse(dat$survival[,2]==1,0,1))
		pval=score1
		eff=eff1=eff2=NA
		q2=(score1>=quantile(score1,0.33,na.rm=T))*1+(score1>=quantile(score1,0.66,na.rm=T))*1 
		dt1=data.frame(dat$surv.dt,scoreq=q2,score=score1)		
		if (length(rmv)>0) dt1=dt1[-rmv,]
		dt1=dt1[dt1$scoreq!=1,]
		if (0 %in% dt1$scoreq & 2 %in% dt1$scoreq){
			fit <- survfit(Surv(time, status) ~ scoreq, data = dt1)
			tst=survdiff(Surv(time, status) ~ scoreq,data=dt1)
			logrank.p <- 1 - pchisq(tst$chisq, length(tst$n) - 1)
			sv1=survfit(Surv(time, status)~1,dt1[dt1$scoreq==2,])
			sv2=survfit(Surv(time, status)~1,dt1[dt1$scoreq==0,])
			eff1=as.numeric(surv_median(sv1)[2])
			eff2=as.numeric(surv_median(sv2)[2])
			eff=eff1-eff2
		}
	
	dtx=NULL
	dtx$logrank.p=logrank.p
	dtx$eff=eff
	dtx$score=score1
	return(dtx)
}
sr.ranking<-function(sr.final1){
		ps=1-sr.final1[,17]											
		ii=order(ps,decreasing=T)
		return(ii)
}
eval.auc=function(sr.final,dat,ires1,iirs1){
	res=NULL
	sr.x=cbind(prob$genes[sr.final[,1]],prob$genes[sr.final[,2]])
	score1=sr.score(dat,sr.x)
	pval=score1[c(ires1,iirs1)]
	flag=c(rep(1,length(ires1)),rep(0,length(iirs1)))
	res$score=pval
	res$flag=flag
	res$sr.x=sr.x
	return(res)
}
		
get.sr.pairs<-function(sr.final,dr){
	if(dr[ifnames]=="PD1") ix=which(as.numeric(sr.final[,1]) %in% match(c("PD1/PDL1"),prob$genes))
	if(dr[ifnames]=="CTLA4") ix=which(as.numeric(sr.final[,1]) %in% match(c("CTLA4p"),prob$genes))
	if(dr[ifnames]=="PD1/CTLA4") ix=which(as.numeric(sr.final[,1]) %in% match(c("CTLA4p","PD1/PDL1"),prob$genes))
	sr.final1=sr.final[ix,]
	class(sr.final1) <- "numeric"
	return(sr.final1)
}
sr.final=NULL
thr=10
tgs=c("PD1/PDL1","CTLA4p")
for (i in seq(length(tgs))){
	ix=which(as.numeric(sr.tot[,1]) %in% match(tgs[i],prob$genes))
	FDR=0.1;sr.tot1=sr.find(sr.tot[ix,])
	rnk=sr.ranking(sr.tot1) 
	sr.final=rbind(sr.final,sr.tot1[rnk[1:thr],])
}

fnames=dir("./data/immuno/")
fnames=paste0("./data/immuno/",fnames)
# [1] "./data/immuno/dat.Chen.RData"         
# [2] "./data/immuno/dat.Cho.RData"          
# [3] "./data/immuno/dat.Gide.RData"         
# [4] "./data/immuno/dat.Huang.RData"        
# [5] "./data/immuno/dat.Hwang.RData"        
# [6] "./data/immuno/dat.Kim.RData"          
# [7] "./data/immuno/dat.Liu.RData"          
# [8] "./data/immuno/dat.Miao.RData"         
# [9] "./data/immuno/dat.Nathanson_pre.RData"
#[10] "./data/immuno/dat.Prat.RData"         
#[11] "./data/immuno/dat.Riaz.RData"         
#[12] "./data/immuno/dat.Snyder.RData"       
#[13] "./data/immuno/dat.VanAllen.RData" 
dr=c("PD1","PD1","PD1","PD1","PD1","PD1","PD1","PD1","CTLA4","PD1","PD1","PD1","CTLA4")
tp=c("SKCM","LUAD","SKCM","SKCM","LUAD","STAD","SKCM","KIRC","SKCM","SKCM","SKCM","BLCA","SKCM")

datT=list()
for (ifnames in seq(length(fnames))){		
	dat = local({load(fnames[ifnames]);environment()})
	datT[[ifnames]]=dat$dat
}

# auc: predicted AUC of ROC curves
# lgp1: logrank P-values
# eff1: median survival differenes
auc=lgp1=eff1=rep(NA,length(fnames))
for (ifnames in seq(length(fnames))){
	dat=datT[[ifnames]]
	ires1=dat$indR
	iirs1=dat$indNR
	
	sr.final1=get.sr.pairs(sr.final,dr)
	sco=NULL;sco$pval=NA
	sco=eval.auc(sr.final1,dat,ires1,iirs1)
	dtx=drug.predict(sco$score,sco$flag)
	if ("auc" %in% names(dtx)) auc[ifnames]=ceiling(dtx$auc*100)/100;
	if ("survival" %in% names(dat)){
		rmv=NULL
		dkm=eval.km(sr.final1,dat)
		kpm=c(dkm$logrank.p,dkm$eff)
		lgp1[ifnames]=ifelse(kpm[2]>0 | is.na(kpm[2]),kpm[1],1)
		eff1[ifnames]=kpm[2]
	}
}
names(auc)=names(lgp1)=names(eff1)=fnames
