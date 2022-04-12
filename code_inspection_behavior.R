sharing.variants.by.fam = function(pedfile, sfs.file, correction="remove.homo", sites.of.interest=NULL,causal=T){
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
  
  index_col = as.numeric(gsub("locus.", "", colnames(genotypes_affected_sub_unique)))
  index_causal_variants_in_sites = which(sfs.file$V7==1&sfs.file$V1%in%sites.of.interest)
  index_variants_in_sites = which(sfs.file$V1%in%sites.of.interest)
  
  causal_variants_in_col = which(index_col%in%index_causal_variants_in_sites)
  variants_in_col = which(index_col%in%index_variants_in_sites)
  
  if(!is.null(sites.of.interest)) {
    if(causal){
      genotypes_affected_sub_unique = genotypes_affected_sub_unique[,causal_variants_in_col]
    } else {
      genotypes_affected_sub_unique = genotypes_affected_sub_unique[,variants_in_col]
    }
    
  }
  
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

sharing.causal.100.02.CRH = sapply(1:100,function(X) sharing.variants.by.fam(ped_files_alter_2causal[X], sfs_100_02, correction = "replace.homo", sites.of.interest = "1045",causal = T))
sharing.causal.100.04.CRH = sapply(1:100,function(X) sharing.variants.by.fam(ped_files_alter_4causal[X], sfs_100_04, correction = "replace.homo", sites.of.interest = "1045",causal = T))
sharing.causal.100.06.CRH = sapply(1:100,function(X) sharing.variants.by.fam(ped_files_alter_6causal[X], sfs_100_06, correction = "replace.homo", sites.of.interest = "1045",causal = T))

sharing.100.02.CRH = sapply(1:100,function(X) sharing.variants.by.fam(ped_files_alter_2causal[X], sfs_100_02, correction = "replace.homo", sites.of.interest = "1045",causal = F))
sharing.100.04.CRH = sapply(1:100,function(X) sharing.variants.by.fam(ped_files_alter_4causal[X], sfs_100_04, correction = "replace.homo", sites.of.interest = "1045",causal = F))
sharing.100.06.CRH = sapply(1:100,function(X) sharing.variants.by.fam(ped_files_alter_6causal[X], sfs_100_06, correction = "replace.homo", sites.of.interest = "1045",causal = F))

df.sharing.causal.variants = data.frame("Mean_Number_of_causal_variants" = c(sharing.causal.100.02.CRH,
                                                                             sharing.causal.100.04.CRH,
                                                                             sharing.causal.100.06.CRH,
                                                                             sharing.100.02.CRH,
                                                                             sharing.100.04.CRH,
                                                                             sharing.100.06.CRH), "Proportion_Causal"=c(rep(c(rep(0.2,100), rep(0.4,100), rep(0.6,100)),2)), "Type" = c(rep("Causal",300), rep("All", 300)))

ggplot(df.sharing.causal.variants, aes(x=Type,y=Mean_Number_of_causal_variants, fill=factor(Proportion_Causal)))+geom_boxplot()+xlab("Type")+ylab("Mean Number of Variants by Family")+labs(fill="Proportion Causal")

aggregate_100.02.CRH = lapply(1:100, function(X) aggregate.geno.by.fam(ped_files_alter_2causal[X],correction = "replace.homo", FamID = null$FamID))
aggregate_100.04.CRH = lapply(1:100, function(X) aggregate.geno.by.fam(ped_files_alter_4causal[X],correction = "replace.homo",FamID = null$FamID))
aggregate_100.06.CRH = lapply(1:100, function(X) aggregate.geno.by.fam(ped_files_alter_6causal[X],correction = "replace.homo",FamID = null$FamID))

variance_100.02.CRH = lapply(1:100, function(X) compute.Var.by.annot(null,aggregate_100.02.CRH[[X]], Z, W))
variance_100.04.CRH = lapply(1:100, function(X) compute.Var.by.annot(null,aggregate_100.04.CRH[[X]], Z, W))
variance_100.06.CRH = lapply(1:100, function(X) compute.Var.by.annot(null,aggregate_100.06.CRH[[X]], Z, W))

stats_100.02.CRH = lapply(1:100, function(X) compute.Stats.by.annot(null,aggregate_100.02.CRH[[X]], Z, W))
stats_100.04.CRH = lapply(1:100, function(X) compute.Stats.by.annot(null,aggregate_100.04.CRH[[X]], Z, W))
stats_100.06.CRH = lapply(1:100, function(X) compute.Stats.by.annot(null,aggregate_100.06.CRH[[X]], Z, W))

plot(sharing.100.06.CRH, sapply(variance_100.06.CRH,function(x) x$Score1))
plot(sharing.100.04.CRH, sapply(variance_100.04.CRH,function(x) x$Score1))
plot(sharing.100.02.CRH, sapply(variance_100.02.CRH,function(x) x$Score1))

plot(sharing.100.06.CRH, sapply(stats_100.06.CRH,function(x) x$B[1]))
plot(sharing.100.04.CRH, sapply(stats_100.04.CRH,function(x) x$B[1]))
plot(sharing.100.02.CRH, sapply(stats_100.02.CRH,function(x) x$B[1]))



df.variance.stats.sharing = data.frame("Proportion_Causal"= c(rep(0.6,100), rep(0.4,100), rep(0.2,100)), "Variance" = c(sapply(variance_100.06.CRH,function(x) x$Score1),
                                                                                                                        sapply(variance_100.04.CRH,function(x) x$Score1),
                                                                                                                        sapply(variance_100.02.CRH,function(x) x$Score1)),
                                       "Stats"=c(sapply(stats_100.06.CRH,function(x) x$B[1]),
                                                 sapply(stats_100.04.CRH,function(x) x$B[1]),
                                                 sapply(stats_100.02.CRH,function(x) x$B[1])), "Sharing"= c(sharing.100.06.CRH,sharing.100.04.CRH,sharing.100.02.CRH))

ggplot(df.variance.stats.sharing, aes(x=Sharing,y=Variance))+geom_point(aes(color=factor(Proportion_Causal)))+geom_smooth(method="lm", formula = y~x, aes(color=factor(Proportion_Causal)))+
  xlab("Mean Number of Shared Variants")+ylab("Variance")+labs(color="Proportion Causal")

ggplot(df.variance.stats.sharing, aes(x=Sharing,y=Stats))+geom_point(aes(color=factor(Proportion_Causal)))+geom_smooth(method="lm", formula = y~x, aes(color=factor(Proportion_Causal)))+
  xlab("Mean Number of Shared Variants")+ylab("Score")+labs(color="Proportion Causal")

ggplot(df.variance.stats.sharing, aes(x=Sharing,y=Stats/Variance))+geom_point(aes(color=factor(Proportion_Causal)))+geom_smooth(method="lm", formula = y~x, aes(color=factor(Proportion_Causal)))+
  xlab("Mean Number of Shared Variants")+ylab("Original Burden Test")+labs(color="Proportion Causal")


pedfiles_null = list.files("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data\\Scenario_Null_agg", full.names = T)
sfs_null = read.table("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data\\rare_variants_null.sfs", header=F)

sharing.100.02.null = sapply(1:1000, function(x) sharing.variants.by.fam(pedfiles_null[x], sfs_null, correction = "replace.homo", sites.of.interest = c("1045","1257", "1335", "988"),causal = F))

u = unique(sfs_null$V1)
Z = matrix(0, nrow=nrow(sfs_null),ncol=length(u))
Z_tmp = sapply(1:ncol(Z),function(x){
  Z_col = Z[,x,drop=F]
  w = which(sfs_null$V1==u[x])
  Z_col[w,] = 1
  Z_col
})
Z = cbind(1,Z_tmp)

W = diag(dbeta(sfs_null$V6,1,20),266)
list_null = list()
list_agg_null = list()

for(i in 1:1000){
  print(i)
  #agg_null = aggregate.geno.by.fam(pedfiles_null[i], correction = "replace.homo",FamID = null$FamID)
  #list_agg_null[[i]] = agg_null
  list_null[[i]] = retrofun.RVS(null,list_agg_null[[i]], Z,W)
}

retrofun.RVS(null, list_agg_null[[98]], Z,W)
compute.Var.by.annot(null, list_agg_null[[98]], Z,W)

hist(sapply(list_null, function(x)x$ACAT))

gg_qqplot <- function(ps, ci = 0.95) {
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
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    # geom_line(aes(expected, cupper), linetype = 2, size = 0.5) +
    # geom_line(aes(expected, clower), linetype = 2, size = 0.5) +
    xlab(log10Pe) +
    ylab(log10Po)
}

t = ldply(list_null, "rbind")

gg_qqplot(t$Score1[!is.na(t$Score1)])
