# Predicting the response to cytotoxic/targeted agents
# Only publicly available dataset are shared

Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")

library(data.table)
library(Rcpp)
library(survival)
library(ROCR)
library(caTools)
library(survminer)
library(survival)

# set up working directory
setwd("SELECT")
load("./data/prob.TCGA.targeted.RData");pancancer=prob
load("./data/sl.BRAF.RData");sl.tot1=sl.braf
load("./data/sl.targeted.RData");sl.tot2=sl.targeted
source("./R/source.r",echo=T)

sl.find<-function(sl.ess1,FDR1,FDR){
	sl.res=NA;flag=0
	ix=which(p.adjust(sl.ess1[,3],"BH")<FDR1 | p.adjust(sl.ess1[,4],"BH")<FDR1)
	if (length(ix)>1) {
		sl.ess1=sl.ess1[ix,]
		flag=2
		iz1=which(sl.ess1[,5]<0)			# step II: clinical relevance
		iz2=which(sl.ess1[,9]<0)
		iz1=iz1[which(sl.ess1[iz1,7]<0  & p.adjust(sl.ess1[iz1,8],"BH")<FDR)]
		iz2=iz2[which(sl.ess1[iz2,11]<0  & p.adjust(sl.ess1[iz2,12],"BH")<FDR)]	
		iz=union(iz1,iz2)
		if (length(iz)>1) {
			flag=3
			sl.ess1=sl.ess1[iz,]		
			ix=which(sl.ess1[,14]< .5)											# step III: phylogenetic
			if (length(ix)>0){
				sl.ess1=sl.ess1[ix,]
				flag=4
			} 											# flag: marker of SL partner
		}
	}
	if (flag==4) sl.res=sl.ess1
	return(sl.res)
}
sl.find.wrapper<- function(sl.ess0){
	FDR1=0.01;FDR=0.1;sl.ess1=sl.find(sl.ess0,FDR1,FDR)
	if (length(sl.ess1)<ncol(sl.ess0)*25) {FDR1=0.05;FDR=0.1;sl.ess1=sl.find(sl.ess0,FDR1,FDR)}
    if (length(sl.ess1)<ncol(sl.ess0)*25) {FDR1=0.2;FDR=0.2;sl.ess1=sl.find(sl.ess0,FDR1,FDR)}
    return(sl.ess1)
}
sl.score <- function(dat,sl.ess4){
		thr=0.33;thr2=0.3
		score1=rep(0,ncol(dat$mRNA))	
		target1=match(unique(sl.ess4[,1]),dat$genes)				
		partners1=match(unique(sl.ess4[,2]),dat$genes)
		partners1=partners1[!is.na(partners1)]
		tscore=1
		if (!is.na(target1)){
		if (length(target1)==1) 
					tscore=(dat$mRNA.rank[target1,]>=thr2)*1
		if (length(target1) >1)
					tscore=colSums((dat$mRNA.rank[target1,]>=thr2)*1,na.rm=T)/length(target1)
				}
		if (is.na(tscore)) tscore=1
		if (sum(tscore==0,na.rm=T) > ncol(dat$mRNA)/2) tscore=1
		if (length(partners1)==1) 
					score1=1*(dat$mRNA.rank2[partners1,] <thr)*tscore
		if (length(partners1) >1) 
					score1=colSums(dat$mRNA.rank2[partners1,] <thr,na.rm=T)/length(partners1)*tscore
	return(score1)
}
sl.ranking<-function(sl.ess1){
	cln.r=1-rank.array(												
			apply(cbind(apply(cbind(ifelse(sl.ess1[,5]<0,sl.ess1[,6],1),ifelse(sl.ess1[,7]<0,sl.ess1[,8],1)),1,max,na.rm=T),
			apply(cbind(ifelse(sl.ess1[,9]<0,sl.ess1[,10],1),ifelse(sl.ess1[,11]<0,sl.ess1[,12],1)),1,max,na.rm=T)),1,min,na.rm=T))
		
	ii=order(cln.r,decreasing=T)
	return(ii)
}
eval.auc.sl=function(sl.ess1,dat){
	thr=25
	res=NULL
	if (thr<=nrow(sl.ess1)){
		rnk=sl.ranking(sl.ess1)	
		sl.ess2=sl.ess1[rnk[1:thr],]
	}else{
		sl.ess2=sl.ess1
	}
	sl.x=cbind(prob$genes[sl.ess2[,1]],prob$genes[sl.ess2[,2]])
	score1=sl.score(dat,sl.x)
	if (sum(is.na(score1))<ncol(dat$mRNA)/2){
		score=score1[c(ires1,iirs1)]
		flag=c(rep(1,length(ires1)),rep(0,length(iirs1)))
		if (sum(!is.na(score))>5){
			res$auc=ceiling(getAUC(score,flag)[1]*100)/100
			res$sl.x=cbind(dat$drugs,sl.x,sl.ess2)
			res$score=score
			res$flag=flag
		}
	}
	return(res)
}
eval.km.sl=function(sl.ess1,rmv){
	thr=25
	res=NULL
	logrank.p=eff=NA
	rnk=sl.ranking(sl.ess1)
	sl.ess2=sl.ess1[rnk[1:thr],]
	sl.x=cbind(prob$genes[sl.ess2[,1]],prob$genes[sl.ess2[,2]])
	score1=sl.score(dat,sl.x)
	if (sum(is.na(score1))<ncol(dat$mRNA)/2 & sum(score1==0,na.rm=T)<ncol(dat$mRNA)/2){
		surv=dat$surv.dt
		pval=score1
		if (!is.null(rmv)){surv=dat$surv.dt[-rmv,];pval=pval[-rmv]}
		
		if (sum(!is.na(pval))>5 & sum(!is.na(surv))>5){
			lgp=eff=eff1=eff2=NA
			q0=quantile(pval,1/3,na.rm=T)
			q1=quantile(pval,2/3,na.rm=T)
			if (q0> 0) q2=(pval>=q0)*1+(pval>=q1)*1 
			if (q0==0) q2=(pval>q0)*1+(pval>=q1)*1 
			dt1=data.frame(surv,scoreq=q2,score=pval) 		
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
		}
	}			
	res=c(logrank.p,eff)
	return(res)
}

fnames2=dir("./data/targeted")
gse.list=do.call(rbind,strsplit(fnames2,"[/]"))
gse.list=do.call(rbind,strsplit(gse.list,"[.]"))[,2]

# auc: predicted AUC of ROC curves
# lgp: logrank P-values
setwd("./data/targeted")
npar=thr=25
auc=lgp=rep(NA,length(fnames2))
for (ifnames in seq(length(fnames2))){
	dat = local({load(fnames2[ifnames]);environment()})
	dat=dat$dat
		
	gse1=gse.list[ifnames]
	if (dat$targets[1]=="BRAF"){sl.tot=sl.tot1}else{sl.tot=sl.tot2}	
	ix=which(as.numeric(sl.tot[,1]) %in% match(dat$targets,prob$genes))
	sl.ess0=sl.tot[ix,]
	class(sl.ess0) <- "numeric"

	ires1=dat$indR
	iirs1=dat$indNR
    sl.ess1=sl.find.wrapper(sl.ess0)
  
	if (length(sl.ess1)>1){
  		if (length(ires1)>0) {dtx=eval.auc.sl(sl.ess1,dat);auc[ifnames]=dtx$auc}
  		rmv=NULL;if(gse.list[ifnames]=="GSE32603") rmv=which(dat$tempo!="T1")
		if ("survival" %in% names(dat)) {kpm=eval.km.sl(sl.ess1,rmv);lgp[ifnames]=ifelse(kpm[2]>0 | is.na(kpm[2]),kpm[1],1)}
	}
}
names(auc)=names(lgp)=gse.list
