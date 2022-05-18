coxph.run = function(formula.list,data1){
	fmla <- as.formula(paste(" Surv(time,status) ~ ", paste(formula.list, collapse= "+")))
	cox.out = coxph(fmla, data=data1)
	ll = cox.out$loglik[2]
	aa  = summary(cox.out)
	out = ll 
	if("cov" %in% rownames(aa$coefficients) ){
		out = c(aa$coefficients["cov", ], ll)
	}
	out
}
coxph.robust = function(data1, f1.list, f2.list=NULL){
	tryCatch( 
		coxph.run(c(f1.list, f2.list), data1),
		error = function(e)  coxph.run(f1.list, data1)
		)
}
sl.clinical.screen <- function(pair,prob, use.mRNA=F, f2.list = NULL, f1.list = NULL)
{
        if(use.mRNA){
                g1 = prob$mRNA.norm[pair[1],]
                g2 = prob$mRNA.norm[pair[2],]
                f1 = prob$mRNAq2[pair[1],]
                f2 = prob$mRNAq2[pair[2],]
        }else{
                g1 = prob$scna.norm[pair[1],]
                g2 = prob$scna.norm[pair[2],]
                f1 = prob$scnaq2[pair[1],]
                f2 = prob$scnaq2[pair[2],]
        }
        surv.strata = prob$surv.strata
        if(is.null(f1.list)) f1.list =  "strata(type)"
        if(is.null(f2.list)) f2.list = c(  "strata(sex)",  "age", "strata(race)")
        dt1 = cbind(surv.strata, cbind(g1 , g2))
        cntrl.out = coxph.robust(f1.list = c("g1",  "g2",  f1.list), f2.list = f2.list,data1=dt1)
        ll1 = cntrl.out
        cov = ifelse(f2 == 0 & f1 == 0, 1, 0 )
        dt1$cov = qnorm.array(cov)
        cox.out = coxph.robust(f1.list = c("g1",  "g2",  f1.list, "cov"), f2.list = f2.list,data1=dt1)
        uu  = cox.out[1:5]

        dt1 = surv.strata
        cov = ifelse(f2 == 0 & f1==0, 1, 0 )
        dt1$cov = qnorm.array(cov)
        cox.out = coxph.robust(f1.list = c(  f1.list, "cov"), f2.list = f2.list,data1=dt1)
        aa  = cox.out[1:5]
        c(uu, aa)

}
sr.clinical.screen = function(pair,prob, use.mRNA=F, f2.list = NULL, f1.list=NULL)
{
	if(use.mRNA){
		g1 = prob$mRNA.norm[pair[1],]  
		g2 = prob$mRNA.norm[pair[2],]
		f1 = prob$mRNAq2[pair[1],]  
		f2 = prob$mRNAq2[pair[2],]  
	}else{
		g1 = prob$scna.norm[pair[1],]  
		g2 = prob$scna.norm[pair[2],]  
		f1 = prob$scnaq2[pair[1],]  
		f2 = prob$scnaq2[pair[2],]  
	}
	surv.strata = prob$surv.strata
	if(is.null(f2.list)) f2.list = c(  "strata(sex)",  "age", "strata(race)") 
	if(is.null(f1.list)) f1.list =  "strata(types)"
	dt1 = cbind(surv.strata, cbind(g1 , g2))
	cov = ifelse(f2 == 0 & f1>=1, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph.robust(f1.list = c("g1",  "g2",  f1.list, "cov"), f2.list = f2.list,data1=dt1) 
	aa  = cox.out[1:5]
	cov = ifelse(f2 == 0 & f1==0, 1, 0 )
	dt1$cov = qnorm.array(cov)
	cox.out = coxph.robust(f1.list = c("g1",  "g2",  f1.list, "cov"), f2.list = f2.list,data1=dt1) 
	bb  = cox.out[1:5]
	uu = c(aa, bb)
	return(uu)
} 
qnorm.array <- function(mat)     ############# correct one
{
	mat.back = mat 
	mat = mat[!is.na(mat)]
    mat = rank(mat, ties.method = "average");
    mat = qnorm(mat / (length(mat)+1));
    mat.back[!is.na(mat.back)] = mat 
    mat.back
}
rank.array <- function(mat)     ############# correct one
{
	mat.back = mat 
	mat = mat[!is.na(mat)]
    mat = rank(mat, ties.method = "average")/length(mat);
    mat.back[!is.na(mat.back)] = mat 
    mat.back
}
phylo.profile = function(sr.gene.all){
	load("./data/phylogenetic.profile.RData")
	load("./data/feature.weight.RData")
	sr.gene1 = sr.gene.all
	sr.phylo =  cbind(match(sr.gene1$rescuer, phylo$genes), match(sr.gene1$vulnerable, phylo$genes))
	featureMat = (phylo[sr.phylo[,1],-(1:3)] - phylo[sr.phylo[,2],-(1:3)])^2
	featureMat=as.matrix(featureMat)
	class(featureMat) <- "numeric"	
	featureMat %*% t(feature.weight)
}
rank.norm= function(ss) {
	out = rank(ss, na.last="keep")/max(1, sum(!is.na(ss)))
	out[is.na(out)] = 1
	out
}
genomic.instability =function(scna) {
	scna.abs = abs(scna)
	colSums(scna.abs > 1,na.rm=T)/apply(scna.abs, 2, function(tt) sum(!is.na(tt)))
}
getAUC <- function(pval,flag){
	na.inx=which(!is.na(flag) & !is.na(pval))
	pval=pval[na.inx]
	flag=flag[na.inx]
	pred <- prediction(pval, flag)
	perf <- performance(pred,"auc")
	auc=perf@y.values[[1]][1]
	
	perf1 <- performance(pred, measure="prec", x.measure="rec")
	rec=perf1@x.values[[1]]
	prec=recall=perf1@y.values[[1]]
	prr <- trapz(rec[2:length(rec)], prec[2:length(prec)])
	return(c(auc,prr))
}
