#!/usr/bin/env Rscript
set.seed(6753)

options(warn=1)
args = commandArgs(trailingOnly=TRUE)

#Sharing Probabilities
load(args[2])
#Expected values for each family
expected_values_by_fam = data.frame("FamID" = names(forsim.pattern.prob.list), "ExpectedValue"= sapply(1:length(forsim.pattern.prob.list), function(x) sum(forsim.N.list[[x]]*forsim.pattern.prob.list[[x]])))

#Variance
var_alleles=do.call(rbind,lapply(expected_values_by_fam[,1], function(x) sum(forsim.N.list[[x]]^2*forsim.pattern.prob.list[[x]]))) - expected_values_by_fam[,2]^2

#Covariance
covar_alleles=do.call(rbind,lapply(expected_values_by_fam[,1], function(x) 
{
  joint = outer(forsim.pattern.prob.list[[x]],forsim.pattern.prob.list[[x]],pmin)/sum(outer(forsim.pattern.prob.list[[x]],forsim.pattern.prob.list[[x]],pmin))
  marginal = apply(joint,1,sum)
  moy = sum(forsim.N.list[[x]]*marginal)
  sum(outer(forsim.N.list[[x]],forsim.N.list[[x]],"*")*joint) - moy^2
}))

covar_alleles=pmin(covar_alleles,var_alleles)

var_alleles=cbind(expected_values_by_fam[,1], data.frame(var_alleles))
covar_alleles=cbind(expected_values_by_fam[,1], data.frame(covar_alleles))

#Defining where the aggregated ped files are
#need to be splitted with ped and MAF directories respectively
directory = args[1]

print(directory)

ped_files_agg = list.files(paste0(directory,"ped"), pattern=".ped", full.names = T)
MAF_files = list.files(paste0(directory,"MAF"), pattern=".maf", full.names = T)

print(paste0("Looping over ", length(ped_files_agg) , " ped files..."))
df = data.frame("Stats" = rep(NA, length(ped_files_agg)), "Var"=rep(NA, length(ped_files_agg)))
Stats_with_Var = for(i in 1:length(ped_files_agg)){
  
  df_agg_genos = read.table(ped_files_agg[i], header=T)
    
  MAF = read.table(MAF_files[i], header=F)
  w = dbeta(MAF[,1],1,25)
  
  locus_col = as.numeric(sapply(strsplit(colnames(df_agg_genos)[1:(ncol(df_agg_genos)-1)], ".", fixed=T), function(x) x[2]))
  w_sub = w[locus_col]	
    
  split_agg_by_fam = split(df_agg_genos, df_agg_genos$pedigree)
  
  Stats_by_fam = sapply(1:length(split_agg_by_fam), function(x){
    famid = names(split_agg_by_fam)[x]
    famid = substring(famid,2)
    split_agg_by_fam_tmp = as.matrix(split_agg_by_fam[[x]][,-which(colnames(df_agg_genos)=="pedigree")])
    
    split_agg_by_fam_tmp[split_agg_by_fam_tmp==0] = NA
    diff_obs_expect = split_agg_by_fam_tmp - expected_values_by_fam[expected_values_by_fam$FamID == famid, 2]
    diff_obs_expect[is.na(diff_obs_expect)] = 0
    
    sum(diff_obs_expect*w_sub)
})
      
  var_geno = list()
  covar_geno = list()
  
  for(famid in names(split_agg_by_fam)){
	genotypes_in_fam = as.numeric(split_agg_by_fam[[famid]][,-which(colnames(split_agg_by_fam[[famid]]) == "pedigree")])
	index_non_null_variants = which(genotypes_in_fam >0)
	non_null_genotypes_in_fam = genotypes_in_fam[index_non_null_variants]
	
	w_subset = w[index_non_null_variants]
	famid_str = substring(famid, 2)

	var_geno[[famid]] = sum(w_subset ^2 * var_alleles[var_alleles[,1]==famid_str,2])
	
	if(length(non_null_genotypes_in_fam)==1){
		covar_geno[[famid]] = var_geno[[famid]]
		
	}

	else{	
		cw = combn(w_subset,2)
		
		if(length(unique(non_null_genotypes_in_fam)) == 1){
			covar_geno[[famid]] = var_geno[[famid]] + 2*sum(sapply(1:ncol(cw), function(i) prod(cw[,i]))) * var_alleles[var_alleles[,1]==famid_str,2]
		}
		else{
			covar_geno[[famid]] = var_geno[[famid]] + 2*sum(sapply(1:ncol(cw),function(i) prod(cw[,i]))) * covar_alleles[covar_alleles[,1]==famid_str,2]
		}
	}
	
}

Stats = sum(Stats_by_fam)
Var = sum(unlist(covar_geno))

df[i,] = c(Stats, Var)

}
print("Stats were produced")
write.table(df, "Stats_with_Var_null_CRH4.mat", col.names=F, row.names=F, quote=F)

Burden_Stats = df[,1]^2/df[,2]
pvalues_init = pchisq(Burden_Stats, 1, lower.tail=F)

#Here we adjust all p-values lower than .05
index_low_pvalues = which(pvalues_init<=0.05)
print(index_low_pvalues)
prob_sharing_by_famid = lapply(1:48, function(x) tapply(forsim.pattern.prob.list[[x]],forsim.N.list[[x]],sum))
names(prob_sharing_by_famid) = names(forsim.pattern.prob.list)

ped = read.table(args[3], header=F)
n_affected_by_fam = table(ped[ped$V6==2,1])

df_geno_low_pvalues = data.frame(matrix(NA,ncol=length(index_low_pvalues), nrow=10000))

for(i in 1:length(index_low_pvalues)){
  df_agg_genos = read.table(ped_files_agg[index_low_pvalues[i]], header=T)
  print(dim(df_agg_genos))
  MAF = read.table(MAF_files[index_low_pvalues[i]], header=F)
  w = dbeta(MAF[,1],1,25)
  if(ncol(df_agg_genos)>2) {
  	n_variants_by_fam = rowSums(df_agg_genos[,-which(colnames(df_agg_genos)=="pedigree")]!=0)
  	index_observed_variants = apply(df_agg_genos[,-which(colnames(df_agg_genos)=="pedigree")], 1, function(x) which(x!=0))}
  else{
	n_variants_by_fam = sum(df_agg_genos[,-which(colnames(df_agg_genos)=="pedigree")]!=0)
        index_observed_variants = which(df_agg_genos[,-which(colnames(df_agg_genos)=="pedigree")]!=0)
	}
  
  n_affected_by_fam_sub = n_affected_by_fam[names(n_affected_by_fam)%in%df_agg_genos$pedigree]
 
  locus_col = as.numeric(sapply(strsplit(colnames(df_agg_genos)[1:(ncol(df_agg_genos)-1)], ".", fixed=T), function(x) x[2]))
  w_sub = w[locus_col]

  split_agg_by_fam = split(df_agg_genos, df_agg_genos$pedigree)
  
  Stats_by_fam_boot = replicate(10000,sum(sapply(1:length(split_agg_by_fam), function(x){
    famid = names(split_agg_by_fam)[x]
    famid = substring(famid,2)

    sample_geno_famid = sample(1:n_affected_by_fam_sub[x], n_variants_by_fam[x],prob=prob_sharing_by_famid[[famid]],replace = T)

    sum((sample_geno_famid - expected_values_by_fam[expected_values_by_fam$FamID == famid, 2])*w_sub[index_observed_variants[[x]]])})))

  df_geno_low_pvalues[,i] = Stats_by_fam_boot
}

write.table(df_geno_low_pvalues, "Genos_boot_by_fam_4.mat", col.names=F, row.names=F, quote=F)
