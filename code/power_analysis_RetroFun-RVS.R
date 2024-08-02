#This code generates data for Type-I error rate 
#The data produced by this script are available in the repo
#Unprocessed data were not provided due to storage limitations
#If you want to run this script with the original data please contact Loic Mangnier at loic.mangnier@gmail.com

library(RetroFunRVS) #Refer to https://github.com/lmangnier/RetroFun-RVS for installing the package
library(ggplot2)
library(ggpubr)

#This function is adapted from agg.genos.by.fam in the RetroFun-RVS package
#The only difference is that in simulation studies we had to repeat some families 
#in order to have reasonable computation time while reaching the observed number of affected (273)
#Letters [A-Z] have been added for repeated pedigrees, function removes this letter for repeated pedigrees,
#the corresponding function is line 42 in the function below

agg.geno.by.fam.for.sim = function(pedfile, FamID=NULL){
  p = read.table(pedfile, header = F)
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
  
  df.genos.agg.by.fam = aggregate(.~pedigree,df.genos.affected, sum)
  
  index_null_fam = which(rowSums(df.genos.agg.by.fam[,-1])==0)
  if(length(index_null_fam) > 0) df.genos.agg.by.fam = df.genos.agg.by.fam[-index_null_fam,]
  
  #This part of the code removes the first character if the pedigree is not present in a user pre-specified list
  #only used for simulation, this part is ignored in real data analysis
  df.genos.agg.by.fam$pedigree = ifelse(df.genos.agg.by.fam$pedigree %in% FamID, df.genos.agg.by.fam$pedigree,sub('^[A-Z]', '', df.genos.agg.by.fam$pedigree))
  
  locus.col = as.numeric(gsub("X", "", colnames(df.genos.agg.by.fam[,-1])))
  
  return(list("ped_agg"=df.genos.agg.by.fam, "index_variants"=locus.col))
  
}

pedfiles_2causal_100_OR5 = list.files("data\\2causal_OR5\\100causal", full.names = T, recursive = T)
pedfiles_2causal_75_OR5 = list.files("data\\2causal_OR5\\75causal", full.names = T, recursive = T)
pedfiles_2causal_50_OR5 = list.files("data\\2causal_OR5\\50causal", full.names = T, recursive = T)

pedfiles_1causal_100_OR5 = list.files("data\\1causal\\100causal", full.names=T, recursive=T)
pedfiles_1causal_75_OR5 = list.files("data\\1causal\\75causal", full.names=T, recursive=T)
pedfiles_1causal_50_OR5 = list.files("data\\1causal\\50causal", full.names=T, recursive=T)

Z = read.table("data\\Annot_Z_CRHs.mat", header=F)
null_without_consanguinity = read.table("data\\null_without_consanguinity.txt", header=TRUE)

variants_by_annot_CRHs = apply(Z,2,function(x) which(x!=0))


filter.rep.n.fam = function(agg_pedfile,variants.by.annot){
  t = agg_pedfile$ped_agg[,-1]
  
  n.unique.fam = sapply(variants.by.annot, function(x){
    subset.variants = which(agg_pedfile$index_variants%in%x)
    
    length(which(rowSums(t[,subset.variants,drop=F]!=0)!=0))
    
    #length(unique(unlist(apply(t[,subset.variants,drop=F], 2, function(y) which(y>0)))))
    
  })
  
  n.unique.fam
}


agg_2causal_100_OR5 = lapply(1:1000, function(x){
  print(x)
  agg.geno.by.fam.for.sim(pedfiles_2causal_100_OR5[[x]], null_without_consanguinity$FamID)
})

agg_2causal_75_OR5 = lapply(1:1000, function(x){
  agg.geno.by.fam.for.sim(pedfiles_2causal_75_OR5[[x]], null_without_consanguinity$FamID)
})

agg_2causal_50_OR5 = lapply(1:1000, function(x){
  agg.geno.by.fam.for.sim(pedfiles_2causal_50_OR5[[x]], null_without_consanguinity$FamID)
})


power_2causal_100_OR5 = lapply(1:1000, function(x) {
  Z.tmp = Z[,-which(filter.rep.n.fam(agg_2causal_100_OR5[[x]], variants_by_annot)<=5), drop=FALSE]
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_2causal_100_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510), independence = F)}
  )


power_2causal_75_OR5 = lapply(1:1000, function(x) {
  Z.tmp = Z[,-which(filter.rep.n.fam(agg_2causal_75_OR5[[x]], variants_by_annot)<=5), drop=FALSE]
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_2causal_75_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510),independence = F)}
  )

power_2causal_50_OR5 = lapply(1:1000, function(x) {
  print(x)
  if(length(which(filter.rep.n.fam(agg_2causal_50_OR5[[x]], variants_by_annot)<=5))>0){
    Z.tmp = Z[,-which(filter.rep.n.fam(agg_2causal_50_OR5[[x]], variants_by_annot)<=5), drop=FALSE]
  } else{
    Z.tmp=Z
  }
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_2causal_50_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510),independence = F)}
  )

saveRDS(power_2causal_100_OR5,"data\\pvalues_alter_2causal_100_OR5.RDS")
saveRDS(power_2causal_75_OR5,"data\\pvalues_alter_2causal_75_OR5.RDS")
saveRDS(power_2causal_50_OR5,"data\\pvalues_alter_2causal_50_OR5.RDS")

df_power_CRHs_Combined_2causal_OR5 = data.frame("Power"=c(sum(sapply(power_2causal_100_OR5, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                          sum(sapply(power_2causal_100_OR5, function(x) x$ACAT)<=8.33e-6)/1000,
                                                          
                                                          sum(sapply(power_2causal_75_OR5, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                          sum(sapply(power_2causal_75_OR5, function(x) x$ACAT)<=8.33e-6)/1000,
                                                          
                                                          sum(sapply(power_2causal_50_OR5, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                          sum(sapply(power_2causal_50_OR5, function(x) x$ACAT)<=8.33e-6)/1000), "Prop" = c(rep(100,2), rep(75,2), rep(50,2)), "Annot" = rep(c("Burden","ACAT-Combined"),3), "Annotation"="CRHs")

df_power_CRHs_Combined_2causal_OR5$Annot = factor(df_power_CRHs_Combined_2causal_OR5$Annot, levels=c("Burden","ACAT-Combined"))
ggplot(df_power_CRHs_Combined_2causal_OR5, aes(x=Annot, y=Power, fill=factor(Prop)))+geom_bar(position = "dodge", stat="identity",color="black")+labs(fill="Proportion")+xlab("Type")+scale_fill_brewer(palette="Pastel1")+theme_bw()+theme(legend.position = "none")


Z_SW = read.table("data\\Annot_Z_SW.mat", header=F)
Z_Pairs = read.table("data\\Annot_Z_pairs.mat", header=F)
Z_Genes = read.table("data\\Annot_Z_genes.mat", header=F)


power_2causal_100_OR5_pairs = lapply(1:1000, function(x) {
  if(length(which(filter.rep.n.fam(agg_2causal_100_OR5[[x]], variants_by_annot_Pairs)<=5))>0){
    Z.tmp = Z_Pairs[,-which(filter.rep.n.fam(agg_2causal_100_OR5[[x]], variants_by_annot_Pairs)<=5), drop=FALSE]
  } else{
    Z.tmp=Z_Pairs
  }
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_2causal_100_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510), independence = F)}
)

power_2causal_75_OR5_pairs = lapply(1:1000, function(x) {
  if(length(which(filter.rep.n.fam(agg_2causal_75_OR5[[x]], variants_by_annot_Pairs)<=5))>0){
    Z.tmp = Z_Pairs[,-which(filter.rep.n.fam(agg_2causal_75_OR5[[x]], variants_by_annot_Pairs)<=5), drop=FALSE]
  } else{
    Z.tmp=Z_Pairs
  }
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_2causal_75_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510), independence = F)}
)

power_2causal_50_OR5_pairs = lapply(1:1000, function(x) {
  if(length(which(filter.rep.n.fam(agg_2causal_50_OR5[[x]], variants_by_annot_Pairs)<=5))>0){
    Z.tmp = Z_Pairs[,-which(filter.rep.n.fam(agg_2causal_50_OR5[[x]], variants_by_annot_Pairs)<=5), drop=FALSE]
  } else{
    Z.tmp=Z_Pairs
  }
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_2causal_50_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510), independence = F)}
)

df_power_CRHs_Combined_2causal_OR5_pairs = data.frame("Power"=c(sum(sapply(power_2causal_100_OR5_pairs, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                         sum(sapply(power_2causal_100_OR5_pairs, function(x) x$ACAT)<=8.33e-6)/1000,
                                                         
                                                         sum(sapply(power_2causal_75_OR5_pairs, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                         sum(sapply(power_2causal_75_OR5_pairs, function(x) x$ACAT)<=8.33e-6)/1000,
                                                         
                                                         sum(sapply(power_2causal_50_OR5_pairs, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                         sum(sapply(power_2causal_50_OR5_pairs, function(x) x$ACAT)<=8.33e-6)/1000), "Prop" = c(rep(100,2), rep(75,2), rep(50,2)), "Annot" = rep(c("Burden","ACAT-Combined"),3),"Annotation"="Pairs")

df_power_CRHs_Combined_2causal_OR5_pairs$Annot = factor(df_power_CRHs_Combined_2causal_OR5_pairs$Annot, levels=c("Burden","ACAT-Combined"))
ggplot(df_power_CRHs_Combined_2causal_OR5_pairs, aes(x=Annot, y=Power, fill=factor(Prop)))+geom_bar(position = "dodge", stat="identity",color="black")+labs(fill="Proportion")+xlab("Type")+scale_fill_brewer(palette="Pastel1")


power_2causal_100_OR5_genes = lapply(1:1000, function(x) {
  if(length(which(filter.rep.n.fam(agg_2causal_100_OR5[[x]], variants_by_annot_Genes)<=5))>0){
    Z.tmp = Z_Genes[,-which(filter.rep.n.fam(agg_2causal_100_OR5[[x]], variants_by_annot_Genes)<=5), drop=FALSE]
  } else{
    Z.tmp=Z_Genes
  }
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_2causal_100_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510), independence = F)}
)

power_2causal_75_OR5_genes = lapply(1:1000, function(x) {
  if(length(which(filter.rep.n.fam(agg_2causal_75_OR5[[x]], variants_by_annot_Genes)<=5))>0){
    Z.tmp = Z_Genes[,-which(filter.rep.n.fam(agg_2causal_75_OR5[[x]], variants_by_annot_Genes)<=5), drop=FALSE]
  } else{
    Z.tmp=Z_Genes
  }
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_2causal_75_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510), independence = F)}
)

power_2causal_50_OR5_genes = lapply(1:1000, function(x) {
  if(length(which(filter.rep.n.fam(agg_2causal_50_OR5[[x]], variants_by_annot_Genes)<=5))>0){
    Z.tmp = Z_Genes[,-which(filter.rep.n.fam(agg_2causal_50_OR5[[x]], variants_by_annot_Genes)<=5), drop=FALSE]
  } else{
    Z.tmp=Z_Genes
  }
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_2causal_50_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510), independence = F)}
)

df_power_CRHs_Combined_2causal_OR5_genes = data.frame("Power"=c(sum(sapply(power_2causal_100_OR5_genes, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                                sum(sapply(power_2causal_100_OR5_genes, function(x) x$ACAT)<=8.33e-6)/1000,
                                                                
                                                                sum(sapply(power_2causal_75_OR5_genes, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                                sum(sapply(power_2causal_75_OR5_genes, function(x) x$ACAT)<=8.33e-6)/1000,
                                                                
                                                                sum(sapply(power_2causal_50_OR5_genes, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                                sum(sapply(power_2causal_50_OR5_genes, function(x) x$ACAT)<=8.33e-6)/1000), "Prop" = c(rep(100,2), rep(75,2), rep(50,2)), "Annot" = rep(c("Burden","ACAT-Combined"),3),"Annotation"="Genes")

df_power_CRHs_Combined_2causal_OR5_genes$Annot = factor(df_power_CRHs_Combined_2causal_OR5_genes$Annot, levels=c("Burden","ACAT-Combined"))
ggplot(df_power_CRHs_Combined_2causal_OR5_genes, aes(x=Annot, y=Power, fill=factor(Prop)))+geom_bar(position = "dodge", stat="identity",color="black")+labs(fill="Proportion")+xlab("Type")+scale_fill_brewer(palette="Pastel1")


power_2causal_100_OR5_SW = lapply(1:1000, function(x) {
  if(length(which(filter.rep.n.fam(agg_2causal_100_OR5[[x]], variants_by_annot_SW)<=5))>0){
    Z.tmp = Z_SW[,-which(filter.rep.n.fam(agg_2causal_100_OR5[[x]], variants_by_annot_SW)<=5), drop=FALSE]
  } else{
    Z.tmp=Z_SW
  }
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_2causal_100_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510), independence = F)}
)

power_2causal_75_OR5_SW = lapply(1:1000, function(x) {
  if(length(which(filter.rep.n.fam(agg_2causal_75_OR5[[x]], variants_by_annot_SW)<=5))>0){
    Z.tmp = Z_SW[,-which(filter.rep.n.fam(agg_2causal_75_OR5[[x]], variants_by_annot_SW)<=5), drop=FALSE]
  } else{
    Z.tmp=Z_SW
  }
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_2causal_75_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510), independence = F)}
)

power_2causal_50_OR5_SW = lapply(1:1000, function(x) {
  if(length(which(filter.rep.n.fam(agg_2causal_50_OR5[[x]], variants_by_annot_SW)<=5))>0){
    Z.tmp = Z_SW[,-which(filter.rep.n.fam(agg_2causal_50_OR5[[x]], variants_by_annot_SW)<=5), drop=FALSE]
  } else{
    Z.tmp=Z_SW
  }
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_2causal_50_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510), independence = F)}
)

df_power_CRHs_Combined_2causal_OR5_SW = data.frame("Power"=c(sum(sapply(power_2causal_100_OR5_SW, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                                sum(sapply(power_2causal_100_OR5_SW, function(x) x$ACAT)<=8.33e-6)/1000,
                                                                
                                                                sum(sapply(power_2causal_75_OR5_SW, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                                sum(sapply(power_2causal_75_OR5_SW, function(x) x$ACAT)<=8.33e-6)/1000,
                                                                
                                                                sum(sapply(power_2causal_50_OR5_SW, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                                sum(sapply(power_2causal_50_OR5_SW, function(x) x$ACAT)<=8.33e-6)/1000), "Prop" = c(rep(100,2), rep(75,2), rep(50,2)), "Annot" = rep(c("Burden","ACAT-Combined"),3), "Annotation"="SW")

df_power_CRHs_Combined_2causal_OR5_SW$Annot = factor(df_power_CRHs_Combined_2causal_OR5_SW$Annot, levels=c("Burden","ACAT-Combined"))
ggplot(df_power_CRHs_Combined_2causal_OR5_SW, aes(x=Annot, y=Power, fill=factor(Prop)))+geom_bar(position = "dodge", stat="identity",color="black")+labs(fill="Proportion")+xlab("Type")+theme_bw()


ggplot(rbind(df_power_CRHs_Combined_2causal_OR5[df_power_CRHs_Combined_2causal_OR5$Annot=="ACAT-Combined",],df_power_CRHs_Combined_2causal_OR5_pairs[df_power_CRHs_Combined_2causal_OR5_pairs$Annot=="ACAT-Combined",],
             df_power_CRHs_Combined_2causal_OR5_genes[df_power_CRHs_Combined_2causal_OR5_genes$Annot=="ACAT-Combined",],df_power_CRHs_Combined_2causal_OR5_SW[df_power_CRHs_Combined_2causal_OR5_SW$Annot=="ACAT-Combined",]), aes(x=Annotation,y=Power,fill=Annotation))+geom_bar(stat="identity",position="dodge", colour="black")+facet_grid(.~Prop)+theme_bw()+theme(legend.position = "none")

saveRDS(list("100causal"=power_2causal_100_OR5,"75causal"=power_2causal_75_OR5, "50causal"=power_2causal_50_OR5),"data\\pvalues_alter_2causal_CRHs.RDS")
saveRDS(list("100causal"=power_2causal_100_OR5_genes,"75causal"=power_2causal_75_OR5_genes, "50causal"=power_2causal_50_OR5_genes),"data\\pvalues_alter_2causal_Genes.RDS")
saveRDS(list("100causal"=power_2causal_100_OR5_pairs,"75causal"=power_2causal_75_OR5_pairs, "50causal"=power_2causal_50_OR5_pairs),"data\\pvalues_alter_2causal_Pairs.RDS")
saveRDS(list("100causal"=power_2causal_100_OR5_SW,"75causal"=power_2causal_75_OR5_SW, "50causal"=power_2causal_50_OR5_SW),"data\\pvalues_alter_2causal_SW.RDS")

agg_1causal_100_OR5 = lapply(1:1000, function(x){
  print(x)
  agg.geno.by.fam.for.sim(pedfiles_1causal_100_OR5[[x]], null_without_consanguinity$FamID)
})

agg_1causal_75_OR5 = lapply(1:1000, function(x){
  agg.geno.by.fam.for.sim(pedfiles_1causal_75_OR5[[x]], null_without_consanguinity$FamID)
})

agg_1causal_50_OR5 = lapply(1:1000, function(x){
  agg.geno.by.fam.for.sim(pedfiles_1causal_50_OR5[[x]], null_without_consanguinity$FamID)
})


power_1causal_100_OR5 = lapply(1:1000, function(x) {
  print(x)
  if(length(which(filter.rep.n.fam(agg_1causal_100_OR5[[x]], variants_by_annot_CRHs)<=5))>0){
    Z.tmp = Z[,-which(filter.rep.n.fam(agg_1causal_100_OR5[[x]], variants_by_annot_CRHs)<=5), drop=FALSE]
  } else{
    Z.tmp=Z
  }
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_1causal_100_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510),independence = F)}
)


power_1causal_75_OR5 = lapply(1:1000, function(x) {
  print(x)
  if(length(which(filter.rep.n.fam(agg_1causal_75_OR5[[x]], variants_by_annot_CRHs)<=5))>0){
    Z.tmp = Z[,-which(filter.rep.n.fam(agg_1causal_75_OR5[[x]], variants_by_annot_CRHs)<=5), drop=FALSE]
  } else{
    Z.tmp=Z
  }
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_1causal_75_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510),independence = F)}
)

power_1causal_50_OR5 = lapply(1:1000, function(x) {
  print(x)
  if(length(which(filter.rep.n.fam(agg_1causal_50_OR5[[x]], variants_by_annot_CRHs)<=5))>0){
    Z.tmp = Z[,-which(filter.rep.n.fam(agg_1causal_50_OR5[[x]], variants_by_annot_CRHs)<=5), drop=FALSE]
  } else{
    Z.tmp=Z
  }
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_1causal_50_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510),independence = F)}
)


df_power_CRHs_Combined_1causal_OR5 = data.frame("Power"=c(sum(sapply(power_1causal_100_OR5, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                          sum(sapply(power_1causal_100_OR5, function(x) x$ACAT)<=8.33e-6)/1000,
                                                          
                                                          sum(sapply(power_1causal_75_OR5, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                          sum(sapply(power_1causal_75_OR5, function(x) x$ACAT)<=8.33e-6)/1000,
                                                          
                                                          sum(sapply(power_1causal_50_OR5, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                          sum(sapply(power_1causal_50_OR5, function(x) x$ACAT)<=8.33e-6)/1000), "Prop" = c(rep(100,2), rep(75,2), rep(50,2)), "Annot" = rep(c("Burden","ACAT-Combined"),3), "Annotation"="CRHs")

df_power_CRHs_Combined_1causal_OR5$Annot = factor(df_power_CRHs_Combined_1causal_OR5$Annot, levels=c("Burden","ACAT-Combined"))
ggplot(df_power_CRHs_Combined_1causal_OR5, aes(x=Annot, y=Power, fill=factor(Prop)))+geom_bar(position = "dodge", stat="identity",color="black")+labs(fill="Proportion")+xlab("Type")+scale_fill_brewer(palette="Pastel1")+theme_bw()


power_1causal_100_OR5_Pairs = lapply(1:1000, function(x) {
  print(x)
  if(length(which(filter.rep.n.fam(agg_1causal_100_OR5[[x]], variants_by_annot_Pairs)<=5))>0){
    Z.tmp = Z_Pairs[,-which(filter.rep.n.fam(agg_1causal_100_OR5[[x]], variants_by_annot_Pairs)<=5), drop=FALSE]
  } else{
    Z.tmp=Z_Pairs
  }
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_1causal_100_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510),independence = F)}
)


power_1causal_75_OR5_Pairs = lapply(1:1000, function(x) {
  print(x)
  if(length(which(filter.rep.n.fam(agg_1causal_75_OR5[[x]], variants_by_annot_Pairs)<=5))>0){
    Z.tmp = Z_Pairs[,-which(filter.rep.n.fam(agg_1causal_75_OR5[[x]], variants_by_annot_Pairs)<=5), drop=FALSE]
  } else{
    Z.tmp=Z_Pairs
  }
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_1causal_75_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510),independence = F)}
)

power_1causal_50_OR5_Pairs = lapply(1:1000, function(x) {
  print(x)
  if(length(which(filter.rep.n.fam(agg_1causal_50_OR5[[x]], variants_by_annot_Pairs)<=5))>0){
    Z.tmp = Z_Pairs[,-which(filter.rep.n.fam(agg_1causal_50_OR5[[x]], variants_by_annot_Pairs)<=5), drop=FALSE]
  } else{
    Z.tmp=Z_Pairs
  }
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_1causal_50_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510),independence = F)}
)


power_1causal_100_OR5_Genes = lapply(1:1000, function(x) {
  print(x)
  if(length(which(filter.rep.n.fam(agg_1causal_100_OR5[[x]], variants_by_annot_Genes)<=5))>0){
    Z.tmp = Z_Genes[,-which(filter.rep.n.fam(agg_1causal_100_OR5[[x]], variants_by_annot_Genes)<=5), drop=FALSE]
  } else{
    Z.tmp=Z_Genes
  }
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_1causal_100_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510),independence = F)}
)


power_1causal_75_OR5_Genes = lapply(1:1000, function(x) {
  print(x)
  if(length(which(filter.rep.n.fam(agg_1causal_75_OR5[[x]], variants_by_annot_Genes)<=5))>0){
    Z.tmp = Z_Genes[,-which(filter.rep.n.fam(agg_1causal_75_OR5[[x]], variants_by_annot_Genes)<=5), drop=FALSE]
  } else{
    Z.tmp=Z_Genes
  }
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_1causal_75_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510),independence = F)}
)

power_1causal_50_OR5_Genes = lapply(1:1000, function(x) {
  print(x)
  if(length(which(filter.rep.n.fam(agg_1causal_50_OR5[[x]], variants_by_annot_Genes)<=5))>0){
    Z.tmp = Z_Genes[,-which(filter.rep.n.fam(agg_1causal_50_OR5[[x]], variants_by_annot_Genes)<=5), drop=FALSE]
  } else{
    Z.tmp=Z_Genes
  }
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_1causal_50_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510),independence = F)}
)


power_1causal_100_OR5_SW = lapply(1:1000, function(x) {
  print(x)
  if(length(which(filter.rep.n.fam(agg_1causal_100_OR5[[x]], variants_by_annot_SW)<=5))>0){
    Z.tmp = Z_SW[,-which(filter.rep.n.fam(agg_1causal_100_OR5[[x]], variants_by_annot_SW)<=5), drop=FALSE]
  } else{
    Z.tmp=Z_SW
  }
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_1causal_100_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510),independence = F)}
)


power_1causal_75_OR5_SW = lapply(1:1000, function(x) {
  print(x)
  if(length(which(filter.rep.n.fam(agg_1causal_75_OR5[[x]], variants_by_annot_SW)<=5))>0){
    Z.tmp = Z_SW[,-which(filter.rep.n.fam(agg_1causal_75_OR5[[x]], variants_by_annot_SW)<=5), drop=FALSE]
  } else{
    Z.tmp=Z_SW
  }
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_1causal_75_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510),independence = F)}
)

power_1causal_50_OR5_SW = lapply(1:1000, function(x) {
  print(x)
  if(length(which(filter.rep.n.fam(agg_1causal_50_OR5[[x]], variants_by_annot_SW)<=5))>0){
    Z.tmp = Z_SW[,-which(filter.rep.n.fam(agg_1causal_50_OR5[[x]], variants_by_annot_SW)<=5), drop=FALSE]
  } else{
    Z.tmp=Z_SW
  }
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_1causal_50_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510),independence = F)}
)



df_power_CRHs_Combined_1causal_OR5_SW = data.frame("Power"=c(sum(sapply(power_1causal_100_OR5_SW, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                             sum(sapply(power_1causal_100_OR5_SW, function(x) x$ACAT)<=8.33e-6)/1000,
                                                             
                                                             sum(sapply(power_1causal_75_OR5_SW, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                             sum(sapply(power_1causal_75_OR5_SW, function(x) x$ACAT)<=8.33e-6)/1000,
                                                             
                                                             sum(sapply(power_1causal_50_OR5_SW, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                             sum(sapply(power_1causal_50_OR5_SW, function(x) x$ACAT)<=8.33e-6)/1000), "Prop" = c(rep(100,2), rep(75,2), rep(50,2)), "Annot" = rep(c("Burden","ACAT-Combined"),3), "Annotation"="SW")

df_power_CRHs_Combined_1causal_OR5_SW$Annot = factor(df_power_CRHs_Combined_1causal_OR5_SW$Annot, levels=c("Burden","ACAT-Combined"))

df_power_CRHs_Combined_1causal_OR5_Genes = data.frame("Power"=c(sum(sapply(power_1causal_100_OR5_Genes, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                             sum(sapply(power_1causal_100_OR5_Genes, function(x) x$ACAT)<=8.33e-6)/1000,
                                                             
                                                             sum(sapply(power_1causal_75_OR5_Genes, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                             sum(sapply(power_1causal_75_OR5_Genes, function(x) x$ACAT)<=8.33e-6)/1000,
                                                             
                                                             sum(sapply(power_1causal_50_OR5_Genes, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                             sum(sapply(power_1causal_50_OR5_Genes, function(x) x$ACAT)<=8.33e-6)/1000), "Prop" = c(rep(100,2), rep(75,2), rep(50,2)), "Annot" = rep(c("Burden","ACAT-Combined"),3), "Annotation"="Genes")

df_power_CRHs_Combined_1causal_OR5_Genes$Annot = factor(df_power_CRHs_Combined_1causal_OR5_Genes$Annot, levels=c("Burden","ACAT-Combined"))


df_power_CRHs_Combined_1causal_OR5_Pairs = data.frame("Power"=c(sum(sapply(power_1causal_100_OR5_Pairs, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                             sum(sapply(power_1causal_100_OR5_Pairs, function(x) x$ACAT)<=8.33e-6)/1000,
                                                             
                                                             sum(sapply(power_1causal_75_OR5_Pairs, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                             sum(sapply(power_1causal_75_OR5_Pairs, function(x) x$ACAT)<=8.33e-6)/1000,
                                                             
                                                             sum(sapply(power_1causal_50_OR5_Pairs, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                             sum(sapply(power_1causal_50_OR5_Pairs, function(x) x$ACAT)<=8.33e-6)/1000), "Prop" = c(rep(100,2), rep(75,2), rep(50,2)), "Annot" = rep(c("Burden","ACAT-Combined"),3), "Annotation"="Pairs")

df_power_CRHs_Combined_1causal_OR5_Pairs$Annot = factor(df_power_CRHs_Combined_1causal_OR5_Pairs$Annot, levels=c("Burden","ACAT-Combined"))

ggplot(rbind(df_power_CRHs_Combined_1causal_OR5[df_power_CRHs_Combined_1causal_OR5$Annot=="ACAT-Combined",],df_power_CRHs_Combined_1causal_OR5_Pairs[df_power_CRHs_Combined_1causal_OR5_Pairs$Annot=="ACAT-Combined",],
             df_power_CRHs_Combined_1causal_OR5_Genes[df_power_CRHs_Combined_1causal_OR5_Genes$Annot=="ACAT-Combined",],df_power_CRHs_Combined_1causal_OR5_SW[df_power_CRHs_Combined_1causal_OR5_SW$Annot=="ACAT-Combined",]), aes(x=Annotation,y=Power,fill=Annotation))+geom_bar(stat="identity",position="dodge", colour="black")+facet_grid(.~Prop)+theme_bw()+theme(legend.position = "none")


saveRDS(list("100causal"=power_1causal_100_OR5,"75causal"=power_1causal_75_OR5, "50causal"=power_1causal_50_OR5),"data\\pvalues_alter_1causal_CRHs.RDS")
saveRDS(list("100causal"=power_1causal_100_OR5_Genes,"75causal"=power_1causal_75_OR5_Genes, "50causal"=power_1causal_50_OR5_Genes),"data\\pvalues_alter_1causal_Genes.RDS")
saveRDS(list("100causal"=power_1causal_100_OR5_Pairs,"75causal"=power_1causal_75_OR5_Pairs, "50causal"=power_1causal_50_OR5_Pairs),"data\\pvalues_alter_1causal_Pairs.RDS")
saveRDS(list("100causal"=power_1causal_100_OR5_SW,"75causal"=power_1causal_75_OR5_SW, "50causal"=power_1causal_50_OR5_SW),"data\\pvalues_alter_1causal_SW.RDS")


pvalues_CHP_75causal = list.files("\\data\\results_RVNPL_2causal_75\\CHP", full.names = T)
pvalues_RV_75causal = list.files("\\data\\results_RVNPL_2causal_75\\RV", full.names = T)

pvalues_ACAT_CHP_75causal_pairs = c()
pvalues_ACAT_CHP_75causal_all = c()

for(i in 1:200){
  print(i)
  df_p = read.table(pvalues_CHP_75causal[i], fill=T, header=F)
  pvalues_ACAT_CHP_75causal_pairs = c(pvalues_ACAT_CHP_75causal_pairs, ACAT::ACAT(df_p[,1]))
  pvalues_ACAT_CHP_75causal_all = c(pvalues_ACAT_CHP_75causal_all, ACAT::ACAT(df_p[,2]))
  
}

sum(pvalues_ACAT_CHP_75causal_pairs<=8.333333e-06)/200
sum(pvalues_ACAT_CHP_75causal_all<=8.333333e-06)/200

pvalues_ACAT_RV_75causal_pairs = c()
pvalues_ACAT_RV_75causal_all = c()

for(i in 1:200){
  print(i)
  df_p = read.table(pvalues_RV_75causal[i], fill=T, header=F)
  pvalues_ACAT_RV_75causal_pairs = c(pvalues_ACAT_RV_75causal_pairs, ACAT::ACAT(df_p[,1]))
  pvalues_ACAT_RV_75causal_all = c(pvalues_ACAT_RV_75causal_all, ACAT::ACAT(df_p[,2]))
  
}

sum(pvalues_ACAT_RV_75causal_pairs<=8.333333e-06)/200
sum(pvalues_ACAT_RV_75causal_all<=8.333333e-06)/200


pedfiles_2causal_smallped_100_OR5 = list.files("data\\data_smallped\\power\\2causal\\100causal", full.names=T, recursive = T)
pedfiles_2causal_smallped_75_OR5 = list.files("data\\data_smallped\\power\\2causal\\75causal", full.names=T, recursive=T)
pedfiles_2causal_smallped_50_OR5 = list.files("data\\data_smallped\\power\\2causal\\50causal", full.names=T, recursive=T)


agg_2causal_smallped_100_OR5 = lapply(1:1000, function(x) agg.geno.by.fam.for.sim(pedfiles_2causal_smallped_100_OR5[x], null_without_consanguinity$FamID))
agg_2causal_smallped_75_OR5 = lapply(1:1000, function(x) agg.geno.by.fam.for.sim(pedfiles_2causal_smallped_75_OR5[x], null_without_consanguinity$FamID))
agg_2causal_smallped_50_OR5 = lapply(1:1000, function(x) agg.geno.by.fam.for.sim(pedfiles_2causal_smallped_50_OR5[x], null_without_consanguinity$FamID))

power_2causal_smallped_100_OR5 = lapply(1:1000, function(x) {
  print(x)
  if(length(which(filter.rep.n.fam(agg_2causal_smallped_100_OR5[[x]], variants_by_annot_CRHs)<=5))>0){
    Z.tmp = Z[,-which(filter.rep.n.fam(agg_2causal_smallped_100_OR5[[x]], variants_by_annot_CRHs)<=5), drop=FALSE]
  } else{
    Z.tmp=Z
  }
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_2causal_smallped_100_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510),independence = F)}
)


power_2causal_smallped_75_OR5 = lapply(1:1000, function(x) {
  print(x)
  if(length(which(filter.rep.n.fam(agg_2causal_smallped_75_OR5[[x]], variants_by_annot_CRHs)<=5))>0){
    Z.tmp = Z[,-which(filter.rep.n.fam(agg_2causal_smallped_75_OR5[[x]], variants_by_annot_CRHs)<=5), drop=FALSE]
  } else{
    Z.tmp=Z
  }
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_2causal_smallped_75_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510),independence = F)}
)

power_2causal_smallped_50_OR5 = lapply(1:1000, function(x) {
  print(x)
  if(length(which(filter.rep.n.fam(agg_2causal_smallped_50_OR5[[x]], variants_by_annot_CRHs)<=5))>0){
    Z.tmp = Z[,-which(filter.rep.n.fam(agg_2causal_smallped_50_OR5[[x]], variants_by_annot_CRHs)<=5), drop=FALSE]
  } else{
    Z.tmp=Z
  }
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity,agg_2causal_smallped_50_OR5[[x]],Z_annot = Z.tmp,W=rep(1,510),independence = F)}
)


df_power_CRHs_Combined_2causal_OR5_smallped = data.frame("Power"=c(sum(sapply(power_2causal_smallped_100_OR5, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                          sum(sapply(power_2causal_smallped_100_OR5, function(x) x$ACAT)<=8.33e-6)/1000,
                                                          
                                                          sum(sapply(power_2causal_smallped_75_OR5, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                          sum(sapply(power_2causal_smallped_75_OR5, function(x) x$ACAT)<=8.33e-6)/1000,
                                                          
                                                          sum(sapply(power_2causal_smallped_50_OR5, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                          sum(sapply(power_2causal_smallped_50_OR5, function(x) x$ACAT)<=8.33e-6)/1000), "Prop" = c(rep(100,2), rep(75,2), rep(50,2)), "Annot" = rep(c("Burden","ACAT-Combined"),3), "Annotation"="CRHs")

df_power_CRHs_Combined_2causal_OR5_smallped$Annot = factor(df_power_CRHs_Combined_2causal_OR5_smallped$Annot, levels=c("Burden","ACAT-Combined"))
ggplot(df_power_CRHs_Combined_2causal_OR5_smallped, aes(x=Annot, y=Power, fill=Annot))+geom_bar(position = "dodge", stat="identity",color="black")+labs(fill="Proportion")+xlab("Type")+theme_bw()+theme(legend.position = "none")+facet_grid(.~Prop)


saveRDS(list("100causal"=power_2causal_smallped_100_OR5,"75causal"=power_2causal_smallped_75_OR5, "50causal"=power_2causal_smallped_50_OR5),"C:\\Users\\loicm\\Desktop\\Projets\\Github\\Simulation_RL\\data\\pvalues_alter_2causal_smallped_CRHs.RDS")


