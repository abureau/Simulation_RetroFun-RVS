library(snpStats)
library(dplyr)
library(plyr)
library(ACAT)
library(ggplot2)

source("../RETROFUN_RVS.R")

compute.signedStats.by.annot = function(null.value.by.fam,aggregate.geno.by.fam,Z_annot,W){
  
  ped_agg = aggregate.geno.by.fam$ped_agg
  
  W_sub = W[aggregate.geno.by.fam$index_variants,aggregate.geno.by.fam$index_variants]
  
  Expected = null.value.by.fam[,c("FamID", "Expected")]
  
  split_G_agg_by_fam = split(ped_agg, 1:nrow(ped_agg))
  
  diff_obs_expected_by_fam = lapply(split_G_agg_by_fam, function(x) {
    pedigree = x$pedigree
    
    x_tmp = x[,-1]
    x_tmp[x_tmp==0] = NA
    doe = x_tmp - Expected[Expected$FamID==pedigree,"Expected"]
  })
  
  diff_obs_expected_all_fam = do.call("rbind", diff_obs_expected_by_fam)
  diff_obs_expected_all_fam[is.na(diff_obs_expected_all_fam)] = 0
  
  
  S_by_var = matrix(colSums(diff_obs_expected_all_fam), nrow=1)
  Vec_unit = matrix(rep(1,ncol(S_by_var)),nrow=ncol(S_by_var))
  
  #Replace NA by zero in Z
  Z_annot[is.na(Z_annot)] = 0
  
  convert.neg.to.pos = function(col){
    if(min(col)<0) col = col - min(col)
    else col
  }
  
  if(ncol(Z_annot) > 1){
    #Make negative scores positive
    Z_annot = apply(Z_annot, 2, convert.neg.to.pos)
  } else{
    if(min(Z_annot)< 0) Z_annot = Z_annot-min(Z_annot)
    else Z_annot = Z_annot
  }
  
  Z_annot_sub = Z_annot[aggregate.geno.by.fam$index_variants,]
  
  # labels_Scores = sapply(1:ncol(Z_annot_sub), function(x) paste0("Score", x))
  # list_Scores_Annot = list()
  
  Burden_Stat_by_annot = colSums(sapply(1:ncol(Z_annot_sub), function(z){
    diag_Z = diag(Z_annot_sub[,z], length(Z_annot_sub[,z]))
    Wz = W_sub%*%diag_Z
    S_by_var%*%Wz
  }))
  
  return(list("B"=Burden_Stat_by_annot))
}


sharing.1000.null = sapply(1:1000, function(x) sharing.variants.by.fam(pedfiles_null[x], sfs_null, correction = "replace.homo", sites.of.interest = c("1045","1257", "1335", "988"),causal = F))
summary(sharing.1000.null)

variance_1000.CRH = lapply(1:1000, function(X) compute.Var.by.annot(null,list_agg_null[[X]], Z, W))
variance = sapply(variance_1000.CRH,function(x) x$Score1)

stats_1000.CRH = lapply(1:1000, function(X) compute.Stats.by.annot(null,list_agg_null[[X]], Z, W))
X2score = sapply(stats_1000.CRH,function(x) x$B[1])

signedStats_1000.CRH = lapply(1:1000, function(X) compute.signedStats.by.annot(null,list_agg_null[[X]], Z, W))
Zscore = sapply(signedStats_1000.CRH,function(x) x$B[1])
summary(Zscore)
hist(Zscore)

# Vérification qu'il n'y a pas d'erreur
sum(abs(Zscore^2 - X2score))

compute.Stats.by.annot(null,list_agg_null[[1]], Z,W=diag(rep(1,nrow(Z))))


# Vérification que la somme du partage est plus élevée
null0 = null
null0$Expected = 0

compute.Stats.by.annot(null0,list_agg_null[[1]], Z,W=diag(rep(1,nrow(Z))))
compute.Stats.by.annot(null0,list_agg_null[[1]], Z,W)

plot(sharing.1000.null, X2score,ylab="Score")
lines(lowess(sharing.1000.null, X2score))

plot(sharing.1000.null, variance,ylab="Variance",ylim=c(0,30000))
lines(lowess(sharing.1000.null, X2score),lty=2)
lines(lowess(sharing.1000.null, variance))

plot(variance,X2score,xlab="Variance",ylab="Score")
lines(lowess(variance,X2score))
lm.fit = lm(X2score~variance)
summary(lm.fit)

plot(X2score,variance,ylab="Variance",xlab="Score")
lines(lowess(X2score,variance))
lm.fit = lm(variance~X2score)
summary(lm.fit)
