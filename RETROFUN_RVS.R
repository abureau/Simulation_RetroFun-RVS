library(snpStats)
library(dplyr)
library(plyr)
library(ACAT)

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

aggregate.geno.by.fam = function(pedfile, correction=c("remove.homo","replace.homo"), replace_ind_geno = T){
  
  sample = read.pedfile(pedfile)
  MAF_file = col.summary(sample$genotype)$MAF
  fam = sample$fam
  affected = fam[fam$affected==2,"member"]
  
  genotypes_all = data.frame(as(sample$genotypes, "numeric"))
  
  #Keep only affected individuals
  genotypes_affected = genotypes_all[rownames(genotypes_all)%in%affected,]
  
  #Remove columns with only reference alleles
  non_zero_cols = which(colSums(genotypes_affected)>0)
  genotypes_affected = select(genotypes_affected, all_of(non_zero_cols))
  
  #Keep variants present in less than half of individuals
  kept_cols =  which(colSums(genotypes_affected>0) <= nrow(genotypes_affected)/2)
  genotypes_affected = select(genotypes_affected, all_of(kept_cols))
  
  if(length(genotypes_affected)==0){
    print("No variants left after cleaning, please provide an other pedfile")
  }
  
  if(correction=="remove.homo"){
    #Keep only heterozyguous variants
    hetero_variants = which(colSums(genotypes_affected==2) == 0)
    genotypes_affected_sub = select(genotypes_affected, all_of(hetero_variants))
    
  } else if(correction=="replace.homo"){
    #Replace homozygous variants by heterogygous ones
    genotypes_affected[genotypes_affected>=1] = 1 
    genotypes_affected_sub = genotypes_affected
  } else{
    print("Please choose a valid option for correction parameter")
  }
  
  #Remove individuals with only zero genotypes
  non_null_ind = which(rowSums(genotypes_affected_sub>0)>0)
  genotypes_affected_sub = slice(genotypes_affected_sub, non_null_ind)
  
  if(replace_ind_geno){
    genotypes_affected_sub = remove.duplicated.variants.by.ind(genotypes_affected_sub)
  } else {
    genotypes_affected_sub = genotypes_affected_sub
  }
  
  genotypes_affected_sub = merge(genotypes_affected_sub, fam[,c("member", "pedigree")], by="member")
  genotypes_affected_sub$member = NULL
  
  agg_by_fam = aggregate(.~pedigree,genotypes_affected_sub, sum)
  #Useful only in simulation studies since we added a letter for duplicated pedigrees
  #To comment when time gonna come
  agg_by_fam$pedigree = substring(agg_by_fam$pedigree,2)
  
  index_col = as.numeric(gsub("locus.", "", colnames(agg_by_fam)[-1]))
  MAF_file_sub = MAF_file[index_col]
  W = diag(dbeta(MAF_file_sub,1,25),length(index_col))
  
  return(list("ped_agg" = agg_by_fam, "W" = W, "index_variants" = index_col))
  
  
}

#' Remove variants with the same genotype per individual
#' 
#' @param df is dataframe where each row is an individual
#' @return a dataframe with only kept variants
remove.duplicated.variants.by.ind = function(df){
  split_df = split(df, 1:nrow(df))
  
  l = lapply(split_df,function(row){
    row_tmp = row[row>0]
    index = sapply(unique(row_tmp), function(x){
      i = which(row==x)
      if(length(i) == 1) i
      else sample(i,1)
    })
    row[index]
  })
  
  df_output = ldply(l, rbind)
  df_output$member = rownames(df)[as.numeric(df_output$.id)]
  df_output$.id = NULL
  df_output[is.na(df_output)] = 0
  
  return(df_output)
}

#'Function computing the theoritical variance and covariance for each family
#' @param config.by.fam a list with configuration by family
#' @param sharing.proba.by.fam a list with sharing probabilities by family
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
  
  covar=pmin(covar , var)
  
  df_var = data.frame("FamID"=names(config.by.fam), "Var"=var)
  df_covar = data.frame("FamID"=names(config.by.fam), "CoVar"=covar)
  
  return(merge(df_var, df_covar, by="FamID"))
}

#' Compute the Burden Statistic for q functional annotations
#'
#'@param null.value.by.fam is a dataframe with two colums (FamID and Expected) return by compute.null function
#'@param aggregate.geno.by.fam is the list returned by aggregate.geno.by.fam function
#'@param Z is a p*q matrix of functional annotations 
#'
#'@return a list with each Score by annotation

compute.Stats.by.annot = function(null.value.by.fam,aggregate.geno.by.fam,Z_annot){
  
  ped_agg = aggregate.geno.by.fam$ped_agg
  W = aggregate.geno.by.fam$W
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
  
  #Make negative scores positive
  convert.neg.to.pos = function(col){
    if(min(col)<0) col = col - min(col)
    else col
  }
  
  Z_annot = apply(Z_annot, 2, convert.neg.to.pos)
  
  labels_Scores = sapply(1:ncol(Z_annot), function(x) paste0("Score", x))
  list_Scores_Annot = list()
  
  Burden_Stat_by_annot = sapply(1:ncol(Z_annot), function(z){
    diag_Z = diag(Z_annot[,z], length(Z_annot[,z]))
    Wz = W%*%diag_Z
    S_by_var%*%Wz%*%Vec_unit%*%t(Vec_unit)%*%Wz%*%t(S_by_var)
  })
  
  # SKAT_Stat_by_annot = sapply(1:ncol(Z_annot), function(z){
  #    diag_Z = diag(Z_annot[,z], length(Z_annot[,z]))
  #    Wz = (W%*%diag_Z)^2
  #    S_by_var%*%Wz%*%t(S_by_var)
  #  })
  
  # SKAT_O_Stat_by_annot = apply(sapply(seq(0,1,0.01), function(rho){
  #    SKAT_Stat_by_annot*(1-rho) + rho*Burden_Stat_by_annot
  #  }), 1, max)
  # 
  return(list("B"=Burden_Stat_by_annot))#"D"=SKAT_Stat_by_annot))#"C"=SKAT_O_Stat_by_annot))
}

#' Compute the Variance for q functional annotations
#'
#'@param null.value.by.fam is a dataframe with two colums (FamID and Expected) return by compute.null function
#'@param aggregate.geno.by.fam is the list returned by aggregate.geno.by.fam function
#'@param Z is a p*q matrix of functional annotations 
#'
#'@return a list with each variance by annotation
compute.Var.by.annot = function(null.value.by.fam,aggregate.geno.by.fam,Z_annot){
  list_var_Annot = list()
  
  for(col_A in 1:ncol(Z_annot)){
    Score_label = paste0("Score", col_A)
    Wz_tmp = diag(Z_annot[,col_A],nrow(Z_annot))%*%aggregate.geno.by.fam$W
    cov_annot = c()
    
    for(r in 1:nrow(aggregate.geno.by.fam$ped_agg)){
      ped = aggregate.geno.by.fam$ped_agg[r,"pedigree"]
      x = aggregate.geno.by.fam$ped_agg[r,-1]
      diff_x = unique(x[x>0])
      
      i = which(x>0)
      
      if(length(i) == 1) d_Wz = Wz_tmp[i,i]
      else d_Wz = diag(Wz_tmp[i,i])
      
      cov_tmp = NA
      var_tmp = sum(d_Wz^2 * null.value.by.fam[null.value.by.fam$FamID==ped,"Var"])
      
      if(length(i)==1){
        cov_tmp = var_tmp } else{
          cw = combn(d_Wz,2)
          if(length(diff_x) ==1){
            cov_tmp = var_tmp + 2*sum(sapply(1:ncol(cw), function(c) prod(cw[,c]))) * null.value.by.fam[null.value.by.fam$FamID==ped,"Var"]
          } else{
            cov_tmp = var_tmp + 2*sum(sapply(1:ncol(cw),function(c) prod(cw[,c]))) * null.value.by.fam[null.value.by.fam$FamID==ped,"CoVar"]
          }
        }
      
      cov_annot = c(cov_annot,cov_tmp)
    }
    list_var_Annot[[Score_label]] = sum(cov_annot)
    
  }
  return(list_var_Annot)
}

Get.p<-function(Q,re.Q){
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
#'@param
#'@param
#'@param
#'@param
#'
#'@return a vector of p-values for each Score by annotation

retrofun.RVS = function(null.value.by.fam,aggregate.geno.by.fam,Z,adjust.low.p = F, p.threshold = 0.001, n.Boot=NULL, method.adjust = c("Bootstrap", "Lee.Adjust")){
  Score = compute.Stats.by.annot(null.value.by.fam,aggregate.geno.by.fam,Z)
  Var = unlist(compute.Var.by.annot(null.value.by.fam,aggregate.geno.by.fam,Z))
  
  Stats = lapply(Score, function(x){
    S = x / Var
    S[is.na(S)] = 0
    S
  })
  
  p = lapply(Stats, function(x) pchisq(x,1,lower.tail = F))
  p = lapply(p, function(x) x[x!=1])
  
  p_tmp = unlist(p)
  Score_tmp = Score$B
  
  low_p = which(p_tmp<=p.threshold)
  
  if(adjust.low.p==F){
    df_p = data.frame(do.call(rbind,p))
    
    df_p$ACAT = apply(df_p,1, ACAT)
    df_p$Fisher = apply(df_p,1, function(x){
      pchisq(-2*sum(log(x)),2* length(x), lower.tail = F)
    })
  }
  else if(length(low_p) ==0){
    df_p = data.frame(do.call(rbind,p))
    
    df_p$ACAT = apply(df_p,1, ACAT)
    df_p$Fisher = apply(df_p,1, function(x){
      pchisq(-2*sum(log(x)),2* length(x), lower.tail = F)
    })
  }
  
  else{
    if(is.null(n.Boot)) print("Please specify a number of Bootstraps")
    
    else{
    
      agg_boot = replicate(n.Boot, resample.genos.by.fam(aggregate.geno.by.fam))
      
      if(length(low_p)==1) {
        Stats_boot = apply(agg_boot,2, function(x) compute.Stats.by.annot(null.value.by.fam, x, Z[,c(1,low_p)])$B[2])
        
        if(method.adjust=="Bootstrap") {
          p_boot = sum(floor(Stats_boot)>=floor(Score_tmp[low_p]))/n.Boot
          p_tmp[low_p] = p_boot}
        else {
          p_tmp[low_p] = Get.p(floor(Score_tmp[low_p]),floor(Stats_boot))
        }
      }
      else {
        Stats_boot = apply(agg_boot,2, function(x) compute.Stats.by.annot(null.value.by.fam, x,Z[,low_p])$B)
        if(method.adjust=="Bootstrap") p_tmp[low_p] = sapply(1:nrow(Stats_boot), function(x) sum(floor(Stats_boot[x,]) >= floor(Score_tmp[low_p[x]]))/n.Boot)
        else p_tmp[low_p] = sapply(1:nrow(Stats_boot), function(x) Get.p(floor(Score_tmp[low_p[x]])),floor(Stats_boot[x,]))
      }
    }
    p$B[low_p] = p_tmp[low_p]
    df_p = data.frame(do.call(rbind,p))
    
    #df_p[df_p==0] = 1e-4
    
    df_p$ACAT = apply(df_p,1, ACAT)
    df_p$Fisher = apply(df_p,1, function(x){
      pchisq(-2*sum(log(x)),2* length(x), lower.tail = F)
    
  })
  # Stats[is.na(Stats)] = 0
  # p = pchisq(Stats,1, lower.tail = F)
  # p = p[p!=1]
  # 
  # df_p = data.frame(p)
  # df_p[nrow(df_p)+1,"p"] = ACAT(p)
  # df_p[nrow(df_p)+1,"p"] = pchisq(-2*sum(log(df_p[1:length(p),])),2* length(p), lower.tail = F)
  # rownames(df_p) = c(paste0("Score",1:length(p)), "ACAT", "Fisher")
  }
  return(df_p)
}
#Procedure test
null = compute.null(forsim.N.list, forsim.pattern.prob.list)

ped_files_null = list.files("Null_latest_NONS_MISS_SPLICE_concat\\", full.names=T)
subset_ped_files_null = ped_files_null[1:1000]

Z_CRHs = read.table("annotation_null.mat", header=F)
Z_quant = cbind(1,rchisq(nrow(Z_CRHs),1))

#' Function resampling genotypes by fam 
#' @param aggregate.geno.by.fam a list return by aggregate.geno.by.fam
#' @return a list 
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

l_CRHs = list()
l_quant = list()
c=1

for(f in subset_ped_files_null){
  print(c)
  agg = aggregate.geno.by.fam(f, correction = "remove.homo", replace_ind_geno = T )
  
  Z_CRHs_sub = Z_CRHs[agg$index_variants,]
  Z_quant_sub = Z_quant[agg$index_variants,]
  
  retrofun_CRHS = retrofun.RVS(null,agg,Z_CRHs_sub, adjust.low.p = T, p.threshold = 0.001, n.Boot=1000, method.adjust="Bootstrap")
  retrofun_quant = retrofun.RVS(null,agg,Z_quant_sub,adjust.low.p = T, p.threshold = 0.001, n.Boot=1000, method.adjust="Bootstrap")
  
  l_CRHs[[f]]=retrofun_CRHS
  l_quant[[f]]=retrofun_quant
  c = c + 1
}
saveRDS(l_CRHs, "Burden_CRHs_10000rep.RDS")
saveRDS(l_quant, "Burden_Quant_1_10000rep.RDS")

Null_Burden_CRHs = readRDS("Burden_CRHs_10000rep.RDS")
Null_Burden_quant = readRDS("Burden_Quant_1_10000rep.RDS")

ACAT_null_CRHs = sapply(Null_Burden_CRHs, function(x) x$ACAT)
low_ACAT_CRHs = which(ACAT_null_CRHs<=0.001)
files_low_ACAT = ped_files_null[]

Fisher_null_CRHs = sapply(Null_Burden_CRHs, function(x) x$Fisher)

sum(ACAT_null_CRHs<=0.05)/10000
sum(ACAT_null_CRHs<=0.01)/10000
sum(ACAT_null_CRHs<=0.005)/10000
sum(ACAT_null_CRHs<=0.001)/10000
sum(ACAT_null_CRHs<=0.00001)/10000

ACAT_null_quant = sapply(Null_Burden_quant, function(x) x$ACAT)

sum(ACAT_null_quant<=0.05)/10000
sum(ACAT_null_quant<=0.01)/10000
sum(ACAT_null_quant<=0.005)/10000
sum(ACAT_null_quant<=0.001)/10000
sum(ACAT_null_quant<=0.00001)/10000

Fisher_null_quant = sapply(Null_Burden_quant, function(x) x$Fisher)


Get.p(3026.3433,rnorm(10000,0,sqrt(242.1075)))
