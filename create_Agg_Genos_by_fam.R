#!/usr/bin/env Rscript
set.seed(12345)
library(plyr, quietly=T)
library(snpStats, quietly=T)
library(dplyr, quietly=T)

#options(warn=1)
args = commandArgs(trailingOnly=TRUE)

remove_duplicated_variants = function(x){
  #This function removes duplicated variants in each family
  if(length(split_agg_by_fam[[x]]) > 2){
  	X = as.numeric(split_agg_by_fam[[x]][,-1])
  	diff_variants = unique(X[X!=0])
  	
  	random_index = sapply(diff_variants, function(y) {
		i = which(X==y)
		if(length(i)==1) i +1
		else sample(i, 1,replace=F) +1
	})
  
  split_agg_by_fam[[x]][random_index]} else split_agg_by_fam[[x]][c(2,1)]
}

directory = args[1]
print(directory)

ped_files = list.files(directory, pattern=".ped", full.names = T)

print(paste0("Looping over ", length(ped_files) , " ped files..."))
print("Start of the loop...")

for(i in 1:length(ped_files)){

  	sample = read.pedfile(ped_files[i])
    	MAF = col.summary(sample$genotype)$MAF
	write.table(MAF, paste0(args[2], "MAF/matrix_",i,".maf"), col.names=F, quote=F, row.names=F)
                
  	df_genotypes = data.frame(as(sample$genotypes, "numeric"))
	
  	#Filtering on only-affected individuals
  	affected = sample$fam[sample$fam$affected ==2,]
  	df_genotypes_affected = df_genotypes[rownames(df_genotypes) %in% affected$member,]
	
	df_genotypes_affected = select(df_genotypes_affected, which(colSums(df_genotypes_affected==0)< nrow(df_genotypes_affected)))

	if(length(df_genotypes_affected) == 0){

		print(paste0(args[2],"ped/rep_agg_",i,".ped"))
		next
	}
	
	homo_variants = which(colSums(df_genotypes_affected==2)>0)#which(apply(df_genotypes_affected,2,function(x) any(x==2)))
	hetero_variants = which(apply(df_genotypes_affected,2,function(x) sum(x!=2) == nrow(df_genotypes_affected)))

	#Removing loci where we observe homozygote variants
	if(length(homo_variants)==0) df_genotypes_affected_sub = df_genotypes_affected
	else df_genotypes_affected_sub = df_genotypes_affected[,-homo_variants]  
	
	if(length(df_genotypes_affected_sub) == 0){
		 next
	}
		
	if(!is.null(ncol(df_genotypes_affected_sub))){
		df_genotypes_affected_sub$member = rownames(df_genotypes_affected)
		df_genotypes_affected_sub = merge(df_genotypes_affected_sub, affected[,c("member", "pedigree")], by="member")
		df_genotypes_affected_sub$member = NULL
	}
	else{
		df_genotypes_affected_tmp = data.frame(df_genotypes_affected_sub, rownames(df_genotypes_affected))
		colnames(df_genotypes_affected_tmp) = c(colnames(df_genotypes_affected)[hetero_variants],"member")
		df_genotypes_affected_tmp = merge(df_genotypes_affected_tmp, affected[,c("member", "pedigree")], by="member")
		df_genotypes_affected_tmp$member = NULL
		df_genotypes_affected_sub = df_genotypes_affected_tmp
	}
	
	agg_by_fam = aggregate(.~pedigree,df_genotypes_affected_sub,sum)
	
	if(ncol(agg_by_fam)>2) agg_by_fam = agg_by_fam[apply(agg_by_fam[,-1],1,function(row) any(row!=0)),]
	else agg_by_fam = agg_by_fam[agg_by_fam[,2]>0,]
	
        split_agg_by_fam = split(agg_by_fam, agg_by_fam$pedigree)
	names(split_agg_by_fam) = agg_by_fam$pedigree

        split_agg_by_fam_clean = lapply(names(split_agg_by_fam),remove_duplicated_variants)
	df_agg_by_fam = ldply(split_agg_by_fam_clean, rbind)

	df_agg_by_fam[is.na(df_agg_by_fam)] = 0
	df_agg_by_fam$pedigree = agg_by_fam$pedigree
	
	write.table(df_agg_by_fam, paste0(args[2],"ped/rep_agg_",i,".ped"), quote=F, col.names=T, row.names=F)
}


