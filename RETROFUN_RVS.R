library(snpStats)
library(dplyr)
library(plyr)
library(ACAT)
library(ggplot2)

set.seed(1234)
setwd("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data")
load("SPAPsimpleprob.RData")

#' Compute the null value for pedigrees
#' 
#' @param pedigree.configurations is all possible configurations for people having the risk allele: A list
#' @param pedigree.probas is all probabilities associated with all pedigree configurations: A list
#' @return The expected genotype value, variance and covariance for each pedigree within a data.frame

compute.null = function(pedigree.configurations, pedigree.probas){
  #pedigree.configurations and pedigree.probas must have the same length
  if(length(pedigree.configurations) != length(pedigree.probas)) print("Please use objects with same length")
  else{
    l = length(pedigree.configurations)
    expected_values = sapply(1:l, function(x) sum(pedigree.configurations[[x]] * pedigree.probas[[x]] ))
    
    df_expected_by_fam = data.frame("FamID" = names(pedigree.configurations), "Expected" = expected_values)
    df_var_covar_by_fam = compute.var.by.fam(pedigree.configurations, pedigree.probas)
    
    df_output = merge(df_expected_by_fam, df_var_covar_by_fam, by="FamID")
    return(df_output)
  }
  
}

#' Aggregate genotypes among affected individuals for each family
#' 
#' @param pedfile is a genotype file in ped format: A .ped file 
#' @param correction is the correction brought to treat heterozygous individuals: A string
#' @param replace_ind_geno allow to remove variants on the same haplotype : A Boolean 
#' @return A list with the ped file corrected and aggregated by family, the weight matrix and index each variants observed in families

aggregate.geno.by.fam = function(pedfile, FamID){
  p = read.table(pedfile[x], header = F)
  fam = p[,1:6]
  affected = which(fam$V6==2)
  genos = p[,7:ncol(p)]
  
  genos[genos==1] = 0
  genos[genos==2] = 1
  
  df.genos = data.frame(do.call("cbind",lapply(seq(2, ncol(genos),2), function(x){
    rowSums(genos[,c(x-1,x)])
  })))
  
  df.genos.affected = df.genos[affected,]
  df.genos.affected = data.frame(t(unique(t(df.genos.affected))))
  df.genos.affected[df.genos.affected==2] = 1
  df.genos.affected = df.genos.affected[,-which(colSums(df.genos.affected)==0)]
  
  
  df.genos.affected$pedigree = fam[affected,"V1"]
  df.genos.affected$pedigree = ifelse(df.genos.affected$pedigree %in% FamID, df.genos.affected$pedigree,sub('^[A-Z]', '', df.genos.affected$pedigree))
  
  df.genos.agg.by.fam = aggregate(.~pedigree,df.genos.affected, sum)
  
  index_null_fam = which(rowSums(df.genos.agg.by.fam[,-1])==0)
  if(length(index_null_fam) > 0) df.genos.agg.by.fam = df.genos.agg.by.fam[-index_null_fam,]
  
  locus.col = as.numeric(gsub("X", "", colnames(df.genos.agg.by.fam[,-1])))
  
  
  return(list("ped_agg"=df.genos.agg.by.fam, "index_variants"=locus.col))
  
}


#'Function computing the theoritical variance and covariance for each family
#' @param config.by.fam a list with configuration by family
#' @param sharing.proba.by.fam a list with sharing probabilities by family
#' 
#' @return a dataframe with Famid, variance and covariance
compute.var.by.fam = function(config.by.fam, sharing.proba.by.fam){
  
  var = sapply(names(config.by.fam), function(famid){
    sum(config.by.fam[[famid]]^2*sharing.proba.by.fam[[famid]]) -
      sum(config.by.fam[[famid]] * sharing.proba.by.fam[[famid]])^2
  })
  
  covar=sapply(names(config.by.fam), function(famid)
  {
    joint = outer(sharing.proba.by.fam[[famid]],sharing.proba.by.fam[[famid]],pmin)/sum(outer(sharing.proba.by.fam[[famid]],sharing.proba.by.fam[[famid]],pmin))
    marginal = apply(joint,1,sum)
    moy = sum(config.by.fam[[famid]]*marginal)
    sum(outer(config.by.fam[[famid]],config.by.fam[[famid]],"*")*joint) - moy^2
  })
  
  covar=pmin(covar, var)
  
  df_var = data.frame("FamID"=names(config.by.fam), "Var"=var)
  df_covar = data.frame("FamID"=names(config.by.fam), "CoVar"=covar)
  
  return(merge(df_var, df_covar, by="FamID"))
}

#' Compute the Burden Statistic for q functional annotations
#'
#'@param null.value.by.fam is a dataframe with two colums (FamID and Expected) return by compute.null function
#'@param aggregate.geno.by.fam is the list returned by aggregate.geno.by.fam function
#'@param Z_annot is a p*q matrix of functional annotations 
#'@param W is a p*p matrix of weights
#'
#'@return a list with each Score by annotation

compute.Burden.by.Annot = function(null.value.by.fam,aggregate.geno.by.fam,Z_annot,W){
  
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
  
  S_by_var = colSums(diff_obs_expected_all_fam)
  Wz_sub = W_sub%*%Z_annot[aggregate.geno.by.fam$index_variants,]
  
  
  S_Wz = S_by_var%*%Wz_sub
  Burden_by_annot = S_Wz^2

  return(list("B"=Burden_by_annot))
  
}

#' Compute the Variance for q functional annotations
#'
#'@param null.value.by.fam is a dataframe with two colums (FamID and Expected) return by compute.null function
#'@param aggregate.geno.by.fam is the list returned by aggregate.geno.by.fam function
#'@param Z_annot is a p*q matrix of functional annotations 
#'@param W is a p*p matrix of weights
#'
#'@return a list with each variance by annotation
compute.Var.by.annot = function(null.value.by.fam,aggregate.geno.by.fam,Z_annot,W, independence=F){
  list_var_Annot = list()
  
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
  
  Wz = W%*%Z_annot
  Wz_sub = Wz[aggregate.geno.by.fam$index_variants,]
  
  for(col_A in 1:ncol(Wz_sub)){
    Score_label = paste0("Score", col_A)
    Wz_sub_annot = Wz_sub[,col_A]
    
    cov_annot = c()
    
    for(r in 1:nrow(aggregate.geno.by.fam$ped_agg)){
      ped = aggregate.geno.by.fam$ped_agg[r,"pedigree"]
      
      x = aggregate.geno.by.fam$ped_agg[r,-1]
      diff_x = unique(x[x>0])
      
      i = which(x>0)
      Wz_sub_annot_fam = Wz_sub_annot[i]
      
      cov_tmp = NA
      var_tmp = sum(Wz_sub_annot_fam^2 * null.value.by.fam[null.value.by.fam$FamID==ped,"Var"])
      
      if(independence ==T){
        cov_annot = c(cov_annot,var_tmp)
      }
      else{
        if(length(i)==1){
          cov_tmp = var_tmp } else{
            cw = combn(Wz_sub_annot_fam,2)
            if(length(diff_x) ==1){
              cov_tmp = var_tmp + 2*sum(sapply(1:ncol(cw), function(c) prod(cw[,c]))) * null.value.by.fam[null.value.by.fam$FamID==ped,"Var"]
            } else{
              cov_tmp = var_tmp + 2*sum(sapply(1:ncol(cw),function(c) prod(cw[,c]))) * null.value.by.fam[null.value.by.fam$FamID==ped,"CoVar"]
            }
          }
        
        cov_annot = c(cov_annot, cov_tmp)
      }
    }
    list_var_Annot[[Score_label]] = sum(cov_annot)
  }
  return(list_var_Annot)
}

#' Method to adjust p-values in high tails of the distribution, considering kurtosis
#' 
#' @param Q a score statistic 
#' @param re.Q a vector the resampled score statistic 
#' 
#' @return a kurtosis adjusted p-values based on Lee et al., 2012
#' 
Get.p = function(Q,re.Q){
  re.mean<-mean(re.Q)
  re.variance<-var(re.Q)
  re.kurtosis<-mean((re.Q-re.mean)^4)/re.variance^2-3
  re.df<-(re.kurtosis>0)*12/re.kurtosis+(re.kurtosis<=0)*100000
  
  re.p<-pchisq((Q-re.mean)*sqrt(2*re.df)/sqrt(re.variance)+re.df,re.df,lower.tail=F)
  return(re.p)
}

#' Compute the p-values associated with each functional annotation
#'
#'@param null.value.by.fam is a dataframe with two colums (FamID and Expected) return by compute.null function
#'@param aggregate.geno.by.fam is the list returned by aggregate.geno.by.fam function
#'@param Z is a p*q matrix of functional annotations 
#'@param W is a p*p matrix of weights 
#'@param adjust.low.p a boolean to adjust extreme low p-values based on resampling method
#'@param p.threshold a numeric is the threshold where p-values below will be adjusted
#'@param n.Boot a numeric is the number of permutations
#'@param method.adjust a string for the permutation method used
#'
#'@return a data.frame of p-values for each Score by annotation and ACAT and Fisher combined p-values

RetroFun.RVS = function(null.value.by.fam,aggregate.geno.by.fam,Z_annot,W, independence=F){
  Burden.Stat = compute.Burden.by.Annot(null.value.by.fam,aggregate.geno.by.fam,Z_annot,W)$B
  Var.Stat = unlist(compute.Var.by.annot(null.value.by.fam,aggregate.geno.by.fam,Z_annot,W,independence = independence))
  
  p = pchisq(Burden.Stat/Var.Stat,1,lower.tail = F)
  
  df.p=data.frame(p)
  colnames(df.p) = paste0("Score",1:length(p))
  
  df.p$ACAT = apply(df.p,1,function(x) ACAT(x[!is.nan(x)]))
  df.p$Fisher = apply(df.p,1, function(x){
    pchisq(-2*sum(log(x[!is.nan(x)])),2* length(x[!is.nan(x)]), lower.tail = F)
  })
  
  return(df.p)
}

null = compute.null(forsim.N.list, forsim.pattern.prob.list)

#' Function resampling genotypes by fam 
#' @param aggregate.geno.by.fam a list return by aggregate.geno.by.fam
#' @return resample data.frame  
#' 
resample.genos.by.fam = function(aggregate.geno.by.fam){
  index_non_null = apply(aggregate.geno.by.fam$ped_agg[,-1],1,function(x) which(x>0))
  n_non_null = apply(aggregate.geno.by.fam$ped_agg[,-1],1,function(x) length(which(x>0)))
  
  agg_tmp = aggregate.geno.by.fam
  agg_tmp_ped_agg = aggregate.geno.by.fam$ped_agg[,-1]
  
  for(x in 1:length(agg_tmp$ped_agg$pedigree)){
    famid = agg_tmp$ped_agg$pedigree[x]
    
    sample_geno = sample(1:length(prob_sharing_by_famid[[famid]]),n_non_null[x])#, prob = prob_sharing_by_famid[[famid]])
    
    agg_tmp_ped_agg[x,index_non_null[[x]]] = sample_geno
  }
  agg_tmp_ped_agg = data.frame("pedigree"=agg_tmp$ped_agg[,1],agg_tmp_ped_agg)
  agg_tmp$ped_agg = agg_tmp_ped_agg
  
  return(agg_tmp)
}
