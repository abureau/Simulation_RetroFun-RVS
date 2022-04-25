aggregate.geno.by.fam_tmp = function(pedfile, FamID){
  rep1 = read.pedfile(pedfile)
  
  fam = rep1$fam
  fam$affected[is.na(fam$affected)] = 1
  
  affected = fam[fam$affected==2,c("member","pedigree")]
  
  genotypes_all = data.frame(as(rep1$genotypes, "numeric"))
  
  #Keep only affected individuals
  genotypes_affected = genotypes_all[rownames(genotypes_all)%in%affected$member,]
  
  genotypes_affected_sub_unique = ldply(apply(genotypes_affected,1,function(x) x[!duplicated(x)]), "rbind")
  genotypes_affected_sub_unique[is.na(genotypes_affected_sub_unique)] = 0
  genotypes_affected_sub_unique[genotypes_affected_sub_unique==2] = 1
  
  member = genotypes_affected_sub_unique[,1]
  genotypes_affected_sub_unique[,1] = NULL
  
  aberr_cols = which(colSums(genotypes_affected_sub_unique)>nrow(genotypes_affected_sub_unique)/2)
  if(length(aberr_cols)>0)  genotypes_affected_sub_unique = genotypes_affected_sub_unique[,-aberr_cols]
  
  index_col = as.numeric(gsub("locus.", "", colnames(genotypes_affected_sub_unique)))
  
  genotypes_affected_sub_unique = genotypes_affected_sub_unique[,colnames(genotypes_affected_sub_unique)[order(index_col)]]
  
  genotypes_affected_sub_unique$member = member
  
  genotypes_affected_sub_unique = merge(genotypes_affected_sub_unique,affected, by="member")
  genotypes_affected_sub_unique$member=NULL
  
  agg_by_fam = aggregate(.~pedigree,genotypes_affected_sub_unique,sum)
  
  if(length(which(rowSums(agg_by_fam[,-1])==0))>0) agg_by_fam = agg_by_fam[-which(rowSums(agg_by_fam[,-1])==0),]
  
  agg_by_fam$pedigree = ifelse(agg_by_fam$pedigree %in% FamID, agg_by_fam$pedigree,sub('^[A-Z]', '', agg_by_fam$pedigree))
  
  return(list("ped_agg" = agg_by_fam, "index_variants" = sort(index_col)))
}
