library(snpStats)
library(dplyr)
library(plyr)

sharing.variants.by.fam = function(pedfile, correction="remove.homo", sites.of.interest=NULL){
  sample = read.pedfile(pedfile)
  fam = sample$fam
  fam$affected[is.na(fam$affected)] = 1
  
  affected = fam[fam$affected==2,c("member","pedigree")]
  
  genotypes_all = data.frame(as(sample$genotypes, "numeric"))
  
  #Keep only affected individuals
  genotypes_affected = genotypes_all[rownames(genotypes_all)%in%affected$member,]
  
  #Remove columns with only reference alleles
  non_zero_cols = which(colSums(genotypes_affected)>0)
  genotypes_affected = select(genotypes_affected, all_of(non_zero_cols))
  
  #Keep variants present in less than half of individuals
  kept_cols =  which(colSums(genotypes_affected>0) <= nrow(genotypes_affected)/2)
  genotypes_affected = select(genotypes_affected, all_of(kept_cols))
  
  if(length(genotypes_affected)==0){
    print("No variants left after cleaning, please provide an other pedfile")
    break;
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
    stop("Please choose a valid option for correction parameter")
  }
  # Remove duplicated variants
  # tmp = remove.duplicated.variants.by.ind(genotypes_affected_sub)
  dup = duplicated(t(genotypes_affected_sub))
  genotypes_affected_sub_unique = genotypes_affected_sub[,!dup]
  
  if (!is.null(sites.of.interest))
    genotypes_affected_sub_unique = genotypes_affected_sub_unique[,sites.of.interest]

  genotypes_affected_sub_unique$pedigree = affected$pedigree
  #Number of carriers by causal variants in each family
  agg_causal_variants = aggregate(.~pedigree,genotypes_affected_sub_unique, function(x) sum(x>0))
  
  #Number of causal variants observed in families
  df_sharing_causal_by_fam = data.frame(cbind(agg_causal_variants$pedigree,rowSums(agg_causal_variants[,-1]>0)))
  colnames(df_sharing_causal_by_fam) = c("FamID", "N_Shared_Variants")
  
  df_sharing_causal_by_fam$N_Shared_Variants=as.numeric(df_sharing_causal_by_fam$N_Shared_Variants)
  df_sharing_causal_by_fam = df_sharing_causal_by_fam[df_sharing_causal_by_fam$N_Shared_Variants>0,]
  mean(df_sharing_causal_by_fam$N_Shared_Variants)
  
}

ped_files_alter_100 = list.files("../Scenario_100_2perccausal_oneCRH", full.names=T)

# Test
sharing.variants.by.fam(ped_files_alter_100[1])
sharing.variants.by.fam(ped_files_alter_100[1],sites=1:5)

sharing.global.100.02 = sapply(ped_files_alter_100, sharing.variants.by.fam)
                                                                                        
                                                                                        