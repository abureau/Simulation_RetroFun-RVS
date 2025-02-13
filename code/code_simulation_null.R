#This code generates data for Type-I error rate 
#The data produced by this script are available in the repo
#Unprocessed data were not provided due to storage limitations
#If you want to run this script with the original data please contact Loic Mangnier at loic.mangnier@gmail.com

library(RetroFunRVS) #Refer to https://github.com/lmangnier/RetroFun-RVS for installing the package
library(ggplot2)
library(ggpubr)

#Functions for generating fancy qqplots
gg_qqplot_facet_grid = function(list_pvalues,ci = 0.95, title="") {
  #n  <- length(ps)
  
  dfs <- lapply(1:length(list_pvalues), function(i) {
    n = length(list_pvalues[[i]])
    data.frame(
      Distance_Kernel = names(list_pvalues)[i],
      observed = -log10(sort(list_pvalues[[i]])),
      expected = -log10(ppoints(n)),
      clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
      cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
    )})
  
  all.dfs = do.call("rbind", dfs)
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  
  ggplot(all.dfs) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed, color=Distance_Kernel), size = 2) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5)+xlab(log10Pe) +
    ylab(log10Po)+theme_bw()#+theme(legend.position = "none")
  
}

gg_qqplot = function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  ggplot(df) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed), size = 3, color="#3288BD") +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po) + theme_bw()
}

#pedfiles_null = list.files("data\\Last_sim_null_10000reps\\TAD", full.names = T)

load("data\\SPAPsimpleprob.RData")

#Files with expected values, variances and covariances under the null across all the families
null_without_consanguinity = read.table("data\\null_without_consanguinity.txt", header=TRUE)
null_with_consanguinity = read.table("data\\null_with_consanguinity.txt",header = TRUE)

#Binary matrices for all the functional annotations
Z = read.table("data\\Annot_Z_CRHs.mat", header=F)
Z_SW = read.table("data\\Annot_Z_SW.mat", header=F)
Z_Pairs = read.table("data\\Annot_Z_pairs.mat", header=F)
Z_Genes = read.table("data\\Annot_Z_genes.mat", header=F)

#This function is adapted from agg.genos.by.fam in the RetroFun-RVS package
#The only difference is that in simulation studies we had to repeat some families 
#in order to have reasonable computation time while reaching the observed number of affected (273)
#Letters [A-Z] have been added for repeated pedigrees, function removes this letter for repeated pedigrees,
#the corresponding function is line 116 in the function below

agg.geno.by.fam.for.sim = function(pedfile, FamID=NULL){
  p = read.table(pedfile, header = F)
  fam = p[,1:6]
  
  #Need to be changed by 1 for scenarios with consanguinity and different pedigree structures: see line 265
  affected = which(fam$V6==2)
  genos = p[,7:ncol(p)]
  
  genos[genos==1] = 0
  genos[genos==2] = 1
  
  df.genos = data.frame(do.call("cbind",lapply(seq(2, ncol(genos),2), function(x){
    rowSums(genos[,c(x-1,x)])
  })))
  
  df.genos.affected = df.genos[affected,]
  df.genos.affected = data.frame(t(unique(t(df.genos.affected))))
  # if(length(which(colSums(df.genos.affected==2)>0))>0){
  #   df.genos.affected = df.genos.affected[,-which(colSums(df.genos.affected==2)>0)]
  # } else{
  #   df.genos.affected = df.genos.affected
  # }
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

pedfiles_agg_null = readRDS("data\\pedfiles_agg_null_10000reps.RDS")

n.unique.config.by.fam = sapply(forsim.N.list, unique)
prob.sharing.by.famid = lapply(1:length(forsim.pattern.prob.list), function(x) tapply(forsim.pattern.prob.list[[x]], forsim.N.list[[x]], sum))

names(prob.sharing.by.famid) = names(forsim.pattern.prob.list)

#Filter steps for keeping only functional annotations with 5 families minimally
variants_by_annot_CRHs = apply(Z,2,function(x) which(x!=0))
variants_by_annot_Pairs = apply(Z_Pairs,2,function(x) which(x!=0))
variants_by_annot_Genes = apply(Z_Genes,2,function(x) which(x!=0))
variants_by_annot_SW = apply(Z_SW,2,function(x) which(x!=0))

filter.rep.n.fam = function(agg_pedfile,variants.by.annot){
  t = agg_pedfile$ped_agg[,-1]
  
  n.unique.fam = sapply(variants.by.annot, function(x){
    subset.variants = which(agg_pedfile$index_variants%in%x)
    
    length(which(rowSums(t[,subset.variants,drop=F]!=0)!=0))
    
    #length(unique(unlist(apply(t[,subset.variants,drop=F], 2, function(y) which(y>0)))))
    
  })
  
  n.unique.fam
}


filter_pedfiles_agg_null=lapply(1:10000, function(X){
  filter.rep.n.fam(pedfiles_agg_null[[X]], variants.by.annot)
})

df.n.fam.by.annot = do.call("rbind",filter_pedfiles_agg_null)


#Generatin p-values for main scenario
pvalues_null = lapply(1:10000, function(rep) {
  print(rep)
  Z.tmp = Z[,-which(filter.rep.n.fam(pedfiles_agg_null[[rep]], variants_by_annot_CRHs)<=5), drop=FALSE]
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity, pedfiles_agg_null[[rep]], Z.tmp,rep(1,510))})


gg_qqplot(sapply(pvalues_null, function(x) x$ACAT))


pvalues_null_indep = lapply(1:10000, function(rep) {
  print(rep)
  Z.tmp = Z[,-which(filter.rep.n.fam(pedfiles_agg_null[[rep]], variants_by_annot_CRHs)<=5), drop=FALSE]
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity, pedfiles_agg_null[[rep]], Z.tmp,rep(1,510), independence = TRUE)})


ggpubr::ggarrange(ggpubr::ggarrange(gg_qqplot(sapply(pvalues_null_indep, function(x) x$ACAT)), labels="A"),ggarrange(gg_qqplot_facet_grid(list("Dependence"=sapply(pvalues_null, function(x) x$Score_V1),
                          "Independence"=sapply(pvalues_null_indep, function(x) x$Score_V1)))+labs(colour="Variant Structure")+theme(legend.position = "bottom"),
                  gg_qqplot_facet_grid(list("Dependence"=sapply(pvalues_null, function(x) x$Score_V2),
                                            "Independence"=unlist(sapply(pvalues_null_indep, function(x) x$Score_V2))))+labs(colour="Variant Structure")+theme(legend.position = "bottom"),
                  gg_qqplot_facet_grid(list("Dependence"=unlist(sapply(pvalues_null, function(x) x$Score_V3)),
                                            "Independence"=unlist(sapply(pvalues_null_indep, function(x) x$Score_V3))))+labs(colour="Variant Structure")+theme(legend.position = "bottom"),
                  gg_qqplot_facet_grid(list("Dependence"=unlist(sapply(pvalues_null, function(x) x$Score_V5)),
                                            "Independence"=unlist(sapply(pvalues_null_indep, function(x) x$Score_V5))))+labs(colour="Variant Structure")+theme(legend.position = "bottom"), labels=c("B","C","D","E"),ncol=2, nrow=2), nrow=2)


#pedfiles.smallped = list.files("D:\\Vraisemblance_retrospective\\Simulation\\data\\Null_smallped\\TAD", full.names = TRUE)

#Generating p-values for small to moderate pedigrees
pvalues.smallped.dep = lapply(1:1000, function(rep){
  print(rep)
  Z.tmp = Z[,-which(filter.rep.n.fam(pedfiles.smallped.agg[[rep]], variants_by_annot_CRHs)<=5), drop=FALSE]
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity, pedfiles.smallped.agg[[rep]], Z.tmp,rep(1,510))
})

pvalues.smallped.inddep = lapply(1:1000, function(rep){
  print(rep)
  Z.tmp = Z[,-which(filter.rep.n.fam(pedfiles.smallped.agg[[rep]], variants_by_annot_CRHs)<=5), drop=FALSE]
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity, pedfiles.smallped.agg[[rep]], Z.tmp,rep(1,510), independence = TRUE)
})

saveRDS(pvalues.smallped.inddep, "data\\pvalues_null_smallped_independence.RDS")
saveRDS(pvalues.smallped.dep, "data\\pvalues_null_smallped_dependence.RDS")

pvalues_null_Pairs = lapply(1:10000, function(rep) {
  print(rep)
  Z.tmp = Z_Pairs[,-which(filter.rep.n.fam(pedfiles_agg_null[[rep]], variants_by_annot_Pairs)<=5), drop=FALSE]
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity, pedfiles_agg_null[[rep]], Z.tmp,rep(1,510))})

pvalues_null_Genes = lapply(1:10000, function(rep) {
  print(rep)
  Z.tmp = Z_Genes[,-which(filter.rep.n.fam(pedfiles_agg_null[[rep]], variants_by_annot_Genes)<=5), drop=FALSE]
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity, pedfiles_agg_null[[rep]], Z.tmp,rep(1,510))})

pvalues_null_SW = lapply(1:10000, function(rep) {
  print(rep)
  
  if(length(which(filter.rep.n.fam(pedfiles_agg_null[[rep]], variants_by_annot_SW)<=5))>0){
    Z.tmp = Z_SW[,-which(filter.rep.n.fam(pedfiles_agg_null[[rep]], variants_by_annot_SW)<=5), drop=FALSE]
  } else{
    Z.tmp=Z_SW
  }
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity, pedfiles_agg_null[[rep]], as.matrix(Z.tmp),rep(1,510))
})

#Sensitivity analysis removing windows with 10,20,30 variants or more

Z_SW_filtered_10 = Z_SW[,-which(colSums(Z_SW)<10)]
Z_SW_filtered_20 = Z_SW[,-which(colSums(Z_SW)<20)]
Z_SW_filtered_30 = Z_SW[,-which(colSums(Z_SW)<30)]

pvalues_null_SW_filtered_30 = lapply(1:10000, function(rep) {
  print(rep)
  
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity, pedfiles_agg_null[[rep]], Z_SW_filtered_30,rep(1,510))
})

saveRDS(pvalues_null_SW, "data\\pvalues_null_SW.RDS")
saveRDS(pvalues_null_Pairs, "data\\pvalues_null_Pairs.RDS")
saveRDS(pvalues_null_Genes, "data\\pvalues_null_Genes.RDS")
saveRDS(pvalues_null, "data\\pvalues_null_CRHs.RDS")
saveRDS(pvalues_null_indep, "data\\pvalues_null_CRHs_indep.RDS")

ggarrange(gg_qqplot(sapply(1:10000, function(x) pvalues_null_Pairs[[x]]$ACAT)),
          gg_qqplot(sapply(1:10000, function(x) pvalues_null_Genes[[x]]$ACAT)),
          gg_qqplot(sapply(1:10000, function(x) pvalues_null_SW[[x]]$ACAT)),
          labels = c("A","B","C"), ncol=3,nrow=1)

#Analysis by pedigree structure
pedfiles_null_two_affected = list.files("data\\Null_two_affected_1000reps_concat", full.names = TRUE)
pedfiles_null_three_affected = list.files("data\\Null_three_affected_1000reps_concat", full.names = TRUE)
pedfiles_null_four_affected = list.files("data\\Null_four_affected_1000reps_concat", full.names = TRUE)
pedfiles_null_five_affected = list.files("data\\Null_five_affected_1000reps_concat", full.names = TRUE)
pedfiles_null_six_affected = list.files("data\\Null_six_affected_1000reps_concat", full.names = TRUE)
pedfiles_null_seven_affected = list.files("data\\Null_seven_affected_1000reps_concat", full.names = TRUE)

#This function is adapted from agg.genos.by.fam in the RetroFun-RVS package
#The only difference is that in simulation studies we had to repeat some families 
#in order to have reasonable computation time while reaching the observed number of affected (273)
#Letters [A-Z] have been added for repeated pedigrees, function removes this letter for repeated pedigrees,
#the corresponding function is line 305 in the function below

#Morever, we repeated the same pedigree several times to ensure reaching the original number of affected
#some additional pre-processing steps were added to be able to apply generic functions from RetroFun-RVS
#For this section, disease status has to be set to 1 and 0, we provided the right function below

agg.geno.by.fam.for.sim = function(pedfile, FamID=NULL){
  p = read.table(pedfile, header = F)
  fam = p[,1:6]
  
  #Need to be changed by 1 for scenarios with consanguinity and different pedigree structures
  affected = which(fam$V6==1)
  genos = p[,7:ncol(p)]
  
  genos[genos==1] = 0
  genos[genos==2] = 1
  
  df.genos = data.frame(do.call("cbind",lapply(seq(2, ncol(genos),2), function(x){
    rowSums(genos[,c(x-1,x)])
  })))
  
  df.genos.affected = df.genos[affected,]
  df.genos.affected = data.frame(t(unique(t(df.genos.affected))))
  # if(length(which(colSums(df.genos.affected==2)>0))>0){
  #   df.genos.affected = df.genos.affected[,-which(colSums(df.genos.affected==2)>0)]
  # } else{
  #   df.genos.affected = df.genos.affected
  # }
  df.genos.affected[df.genos.affected==2] = 1
  df.genos.affected = df.genos.affected[,-which(colSums(df.genos.affected)==0)]
  
  df.genos.affected$pedigree = fam[affected,"V1"]
  
  df.genos.agg.by.fam = aggregate(.~pedigree,df.genos.affected, sum)
  
  index_null_fam = which(rowSums(df.genos.agg.by.fam[,-1])==0)
  if(length(index_null_fam) > 0) df.genos.agg.by.fam = df.genos.agg.by.fam[-index_null_fam,]
  
  df.genos.agg.by.fam$pedigree = ifelse(df.genos.agg.by.fam$pedigree %in% FamID, df.genos.agg.by.fam$pedigree,sub('^[A-Z]', '', df.genos.agg.by.fam$pedigree))
  
  locus.col = as.numeric(gsub("X", "", colnames(df.genos.agg.by.fam[,-1])))
  
  return(list("ped_agg"=df.genos.agg.by.fam, "index_variants"=locus.col))
  
}


pedfiles_agg_null_two_affected = lapply(1:1000, function(rep){
  print(rep)
  agg.geno.by.fam.for.sim(pedfiles_null_two_affected[[rep]],null_without_consanguinity$FamID)})


for(rep in 1:1000){
  print(rep)
  tmp = pedfiles_agg_null_two_affected[[rep]]
  tmp$ped_agg$pedigree = stringr::str_extract(tmp$ped_agg$pedigree,"226",)
  pedfiles_agg_null_two_affected[[rep]] = tmp
}

p.values.two.affected = lapply(1:1000, function(x) {
  Z.tmp = Z[,-which(filter.rep.n.fam(pedfiles_agg_null_two_affected[[rep]], variants_by_annot_CRHs)<=5), drop=FALSE]
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity, pedfiles_agg_null_two_affected[[x]], Z_annot = Z.tmp, W=rep(1,510), independence = FALSE)}
)


pedfiles_agg_null_three_affected = lapply(1:1000, function(rep){
  print(rep)
  agg.geno.by.fam.for.sim(pedfiles_null_three_affected[[rep]],null_without_consanguinity$FamID)})

pedfiles_agg_null_three_affected = lapply(1:1000, function(rep){
  pedfiles_agg_null_three_affected[[rep]]$ped_agg$pedigree = ifelse(pedfiles_agg_null_three_affected[[rep]]$ped_agg$pedigree%in% null_without_consanguinity$FamID, pedfiles_agg_null_three_affected[[rep]]$ped_agg$pedigree, substring(pedfiles_agg_null_three_affected[[rep]]$ped_agg$pedigree, 2 ))
  pedfiles_agg_null_three_affected[[rep]]})

p.values.three.affected = lapply(1:1000, function(x){
  Z.tmp = Z[,-which(filter.rep.n.fam(pedfiles_agg_null_three_affected[[rep]], variants_by_annot_CRHs)<=5), drop=FALSE]
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity, pedfiles_agg_null_three_affected[[x]], Z_annot = Z.tmp, W=rep(1,510), independence = FALSE)
})

pedfiles_agg_null_four_affected = lapply(1:1000, function(rep){
  print(rep)
  agg.geno.by.fam.for.sim(pedfiles_null_four_affected[[rep]],null_without_consanguinity$FamID)})

pedfiles_agg_null_four_affected = lapply(1:1000, function(rep){
  pedfiles_agg_null_four_affected[[rep]]$ped_agg$pedigree = ifelse(pedfiles_agg_null_four_affected[[rep]]$ped_agg$pedigree%in% null_without_consanguinity$FamID, pedfiles_agg_null_four_affected[[rep]]$ped_agg$pedigree, substring(pedfiles_agg_null_four_affected[[rep]]$ped_agg$pedigree, 2 ))
  pedfiles_agg_null_four_affected[[rep]]})

p.values.four.affected = lapply(1:1000, function(x) {
  Z.tmp = Z[,-which(filter.rep.n.fam(pedfiles_agg_null_four_affected[[rep]], variants_by_annot_CRHs)<=5), drop=FALSE]
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity, pedfiles_agg_null_four_affected[[x]], Z_annot = Z.tmp, W=rep(1,510), independence = FALSE)}
)

pedfiles_agg_null_five_affected = lapply(1:1000, function(rep){
  print(rep)
  agg.geno.by.fam.for.sim(pedfiles_null_five_affected[[rep]],null_without_consanguinity$FamID)})

pedfiles_agg_null_five_affected = lapply(1:1000, function(rep){
  pedfiles_agg_null_five_affected[[rep]]$ped_agg$pedigree = ifelse(pedfiles_agg_null_five_affected[[rep]]$ped_agg$pedigree%in% null_without_consanguinity$FamID, pedfiles_agg_null_five_affected[[rep]]$ped_agg$pedigree, stringr::str_sub(pedfiles_agg_null_five_affected[[rep]]$ped_agg$pedigree,-3))
  pedfiles_agg_null_five_affected[[rep]]})

p.values.five.affected = lapply(1:1000, function(x) {
  Z.tmp = Z[,-which(filter.rep.n.fam(pedfiles_agg_null_five_affected[[rep]], variants_by_annot_CRHs)<=5), drop=FALSE]
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity, pedfiles_agg_null_five_affected[[x]], Z_annot = Z.tmp, W=rep(1,510), independence = FALSE)}
  )

pedfiles_agg_null_six_affected = lapply(1:1000, function(rep){
  print(rep)
  agg.geno.by.fam.for.sim(pedfiles_null_six_affected[[rep]],null_without_consanguinity$FamID)})

pedfiles_agg_null_six_affected = lapply(1:1000, function(rep){
  pedfiles_agg_null_six_affected[[rep]]$ped_agg$pedigree = ifelse(pedfiles_agg_null_six_affected[[rep]]$ped_agg$pedigree%in% null_without_consanguinity$FamID, pedfiles_agg_null_six_affected[[rep]]$ped_agg$pedigree, stringr::str_sub(pedfiles_agg_null_six_affected[[rep]]$ped_agg$pedigree,-3))
  pedfiles_agg_null_six_affected[[rep]]})

p.values.six.affected = lapply(1:1000, function(x){
  Z.tmp = Z[,-which(filter.rep.n.fam(pedfiles_agg_null_six_affected[[rep]], variants_by_annot_CRHs)<=5), drop=FALSE]
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity, pedfiles_agg_null_six_affected[[x]], Z_annot = Z.tmp, W=rep(1,510), independence = FALSE)
})


pedfiles_agg_null_seven_affected = lapply(1:1000, function(rep){
  print(rep)
  agg.geno.by.fam.for.sim(pedfiles_null_seven_affected[[rep]],null_without_consanguinity$FamID)})

pedfiles_agg_null_seven_affected = lapply(1:1000, function(rep){
  pedfiles_agg_null_seven_affected[[rep]]$ped_agg$pedigree = ifelse(pedfiles_agg_null_seven_affected[[rep]]$ped_agg$pedigree%in% null_without_consanguinity$FamID, pedfiles_agg_null_seven_affected[[rep]]$ped_agg$pedigree, stringr::str_sub(pedfiles_agg_null_seven_affected[[rep]]$ped_agg$pedigree,-3))
  pedfiles_agg_null_seven_affected[[rep]]})

p.values.seven.affected = lapply(1:1000, function(x) {
  Z.tmp = Z[,-which(filter.rep.n.fam(pedfiles_agg_null_seven_affected[[rep]], variants_by_annot_CRHs)<=5), drop=FALSE]
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity, pedfiles_agg_null_seven_affected[[x]], Z_annot = Z.tmp, W=rep(1,510), independence = FALSE)
  })

gg_qqplot_facet_grid(list("2 affected"=sapply(p.values.two.affected, function(x) x$ACAT),
     "3 affected"=sapply(p.values.three.affected, function(x) x$ACAT),
     "4 affected"=sapply(p.values.four.affected, function(x) x$ACAT),
     "5 affected"=sapply(p.values.five.affected, function(x) x$ACAT),
     "6 affected"=sapply(p.values.six.affected, function(x) x$ACAT),
     "7 affected"=sapply(p.values.seven.affected, function(x) x$ACAT)))+labs(colour="Pedigree")




#Analysis with consanguinity
pedfiles_null_consanguinity = list.files("data\\Null_consanguinity_concat", full.names = TRUE)

#Analysis with only consanguinity 
pedfiles_null_only_consanguinity = list.files("data\\Null_ped_only_consanguinity_concat", full.names = TRUE)

pedfiles_agg_null_only_consanguinity = lapply(1:1000, function(rep) {
  print(rep)
  agg.geno.by.fam.for.sim(pedfiles_null_only_consanguinity[[rep]],null_with_consanguinity$FamID)})

pedfiles_agg_null_only_consanguinity = lapply(1:1000, function(rep){
  pedfiles_agg_null_only_consanguinity[[rep]]$ped_agg$pedigree = ifelse(pedfiles_agg_null_only_consanguinity[[rep]]$ped_agg$pedigree%in% null_without_consanguinity$FamID, pedfiles_agg_null_only_consanguinity[[rep]]$ped_agg$pedigree, stringr::str_sub(pedfiles_agg_null_only_consanguinity[[rep]]$ped_agg$pedigree,-3))
  pedfiles_agg_null_only_consanguinity[[rep]]})

pvalues.null.with.only.consanguinity.corrected = lapply(1:1000, function(x) {
  if(length(which(filter.rep.n.fam(pedfiles_agg_null_only_consanguinity[[x]], variants_by_annot)<5))>0){
    Z.tmp = Z[,-which(filter.rep.n.fam(pedfiles_agg_null_only_consanguinity[[x]], variants_by_annot)<=5), drop=FALSE]
  } else{
    Z.tmp=Z
  }
  RetroFunRVS::RetroFun.RVS(null_with_consanguinity, pedfiles_agg_null_only_consanguinity[[x]], Z.tmp,rep(1,510))}
)

pvalues.null.with.only.consanguinity.uncorrected = lapply(1:1000, function(rep) {
  if(length(which(filter.rep.n.fam(pedfiles_agg_null_only_consanguinity[[rep]], variants_by_annot)<5))>0){
    Z.tmp = Z[,-which(filter.rep.n.fam(pedfiles_agg_null_only_consanguinity[[rep]], variants_by_annot)<=5), drop=FALSE]
  } else{
    Z.tmp=Z
  }
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity, pedfiles_agg_null_only_consanguinity[[rep]], Z.tmp,rep(1,510))}
)

pedfiles_agg_null_consanguinity = lapply(1:1000, function(rep) {
  print(rep)
  agg.geno.by.fam.for.sim(pedfiles_null_consanguinity[[rep]],null_with_consanguinity$FamID)})

pedfiles_agg_null_consanguinity = lapply(1:1000, function(rep){
  pedfiles_agg_null_consanguinity[[rep]]$ped_agg$pedigree = ifelse(pedfiles_agg_null_consanguinity[[rep]]$ped_agg$pedigree%in% null_with_consanguinity$FamID, pedfiles_agg_null_consanguinity[[rep]]$ped_agg$pedigree, sub(".","",pedfiles_agg_null_consanguinity[[rep]]$ped_agg$pedigree))
  pedfiles_agg_null_consanguinity[[rep]]})


pvalues.null.with.consanguinity.corrected = lapply(1:1000, function(rep) {
  
  if(length(which(filter.rep.n.fam(pedfiles_agg_null_consanguinity[[rep]], variants_by_annot)<=5))>0){
    Z.tmp = Z[,-which(filter.rep.n.fam(pedfiles_agg_null_consanguinity[[rep]], variants_by_annot)<=5), drop=FALSE]
  } else{
    Z.tmp=Z
  }
  print(rep)
  RetroFunRVS::RetroFun.RVS(null_with_consanguinity, pedfiles_agg_null_consanguinity[[rep]], Z.tmp,rep(1,510))}
)
pvalues.null.with.consanguinity.unccorrected = lapply(1:1000, function(rep) {
  if(length(which(filter.rep.n.fam(pedfiles_agg_null_consanguinity[[rep]], variants_by_annot)<=5))>0){
    Z.tmp = Z[,-which(filter.rep.n.fam(pedfiles_agg_null_consanguinity[[rep]], variants_by_annot)<=5), drop=FALSE]
  } else{
    Z.tmp=Z
  }
  RetroFunRVS::RetroFun.RVS(null_without_consanguinity, pedfiles_agg_null_consanguinity[[rep]], Z.tmp,rep(1,510))}
)


gg_qqplot_facet_grid(list("Corrected" = sapply(pvalues.null.with.consanguinity.corrected, function(x) x$ACAT),
          "Uncorrected"=sapply(pvalues.null.with.consanguinity.unccorrected, function(x) x$ACAT)))+labs(colour="Correction")

gg_qqplot_facet_grid(list("Corrected" = sapply(pvalues.null.with.only.consanguinity.corrected, function(x) x$ACAT),
                                    "Uncorrected"=sapply(pvalues.null.with.only.consanguinity.uncorrected, function(x) x$ACAT)))+labs(colour="Correction")


saveRDS(pvalues.null.with.only.consanguinity.uncorrected, "data\\pvalues_null_only_consanguinity_uncorrected.RDS")
saveRDS(pvalues.null.with.only.consanguinity.corrected, "data\\pvalues_null_only_consanguinity_corrected.RDS")
saveRDS(pvalues.null.with.consanguinity.corrected, "data\\pvalues_null_consanguinity_corrected.RDS")
saveRDS(pvalues.null.with.consanguinity.unccorrected, "data\\pvalues_null_only_consanguinity_uncorrected.RDS")


saveRDS(p.values.two.affected, "data\\pvalues_null_two_affected.RDS")
saveRDS(p.values.three.affected, "data\\pvalues_null_three_affected.RDS")
saveRDS(p.values.four.affected, "data\\pvalues_null_four_affected.RDS")
saveRDS(p.values.five.affected, "data\\pvalues_null_five_affected.RDS")
saveRDS(p.values.six.affected, "data\\pvalues_null_six_affected.RDS")
saveRDS(p.values.seven.affected, "data\\pvalues_null_seven_affected.RDS")


#Code for producing bootstrap-based pvalues
#First as suggested in the bootRetroFunRVS documentation n.unique.config.by.fam and prob.sharing.by.fam have to be in the same order
n.unique.config.by.fam = lapply(n.unique.config.by.fam, function(fam) sort(fam))
prob.sharing.by.famid

bootstrap.samples.CRHs = bootRetroFunRVS(agg.genos.by.fam = pedfiles_agg_null[[786]] ,n.unique.config.by.fam = n.unique.config.by.fam,Z_annot = Z,W = rep(1,510),prob.sharing.by.fam = prob.sharing.by.famid,nullval = null_without_consanguinity ,nrep = 1000)
Stat = compute.Burden.by.Annot(null.value.by.fam = null_without_consanguinity ,aggregate.geno.by.fam = pedfiles_agg_null[[786]] ,Z_annot = Z,W = rep(1,510))$B
Bootstrap.p = rowMeans(apply(bootstrap.samples.CRHs, 1, function(row) row>=Stat))

#Code for computing average running time for all the functional annotations
list.time.running = list("CRHs"=c(),"Genes"=c(), "Pairs"=c(), "SW"=c())

for(i in 1:10){
  print(i)
  start_time <- Sys.time()
  bootRetroFunRVS(pedfiles_agg_null[[1000]], n.unique.config.by.fam, Z_annot = Z,W = rep(1,510),prob.sharing.by.fam = prob.sharing.by.famid,nullval = null_without_consanguinity,nrep = 100)
  end_time <- Sys.time()
  list.time.running$CRHs = c(list.time.running$CRHs,end_time-start_time)
}

for(i in 1:10){
  print(i)
  start_time <- Sys.time()
  bootRetroFunRVS(pedfiles_agg_null[[1000]], n.unique.config.by.fam, Z_annot = Z_Genes,W = rep(1,510),prob.sharing.by.fam = prob.sharing.by.famid,nullval = null_without_consanguinity,nrep = 100)
  end_time <- Sys.time()
  list.time.running$Genes = c(list.time.running$Genes,end_time-start_time)
}


for(i in 1:10){
  print(i)
  start_time <- Sys.time()
  bootRetroFunRVS(pedfiles_agg_null[[1000]], n.unique.config.by.fam, Z_annot = Z_Pairs,W = rep(1,510),prob.sharing.by.fam = prob.sharing.by.famid,nullval = null_without_consanguinity,nrep = 100)
  end_time <- Sys.time()
  list.time.running$Pairs = c(list.time.running$Pairs,end_time-start_time)
}


for(i in 1:10){
  print(i)
  start_time <- Sys.time()
  bootRetroFunRVS(pedfiles_agg_null[[1000]], n.unique.config.by.fam, Z_annot = Z_SW,W = rep(1,510),prob.sharing.by.fam = prob.sharing.by.famid,nullval = null_without_consanguinity,nrep = 100)
  end_time <- Sys.time()
  list.time.running$SW = c(list.time.running$SW,end_time-start_time)
}



sum(list.time.running$CRHs)*10
sum(list.time.running$Genes)*10
sum(list.time.running$Pairs)*10
sum(list.time.running$SW)*10
