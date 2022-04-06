library(snpStats)
library(dplyr)
library(plyr)
library(RVS)
library(kinship2)
load("../SPAPsimpleprob.RData")
load("../pedsimple.RData")

correction="remove.homo"

# Largeur de fenêtre
lf = 5

ped_files_alter_90 = list.files("../sample_alter_10pcausal_90pinCRHs", full.names=T)

# Test avec un seul réplicat
pedfile = ped_files_alter_90[1]

#ped_files_alter_75 = list.files("../CRHs_alter_20_75_agg", full.names=T)
#pedfile = ped_files_alter_75[1]

  sample = read.pedfile(pedfile)
  MAF_file = col.summary(sample$genotype)$MAF
  fam = sample$fam
  affected = fam[fam$affected==2,"member"]

  # Créer pattern.prob.list et N.list
  fam.vec = unique(fam$pedigree)
  fam.type = substring(fam.vec,2)
  ech.pattern.prob.list = forsim.pattern.prob.list[fam.type]
  ech.N.list = forsim.N.list[fam.type]
  names(ech.pattern.prob.list) = fam.vec
  names(ech.N.list) = fam.vec
  
  # Création d'un objet pedigreeList, puis élagage des familles et sauvegarde dans une liste d'objets pedigree
  ech.ped = pedigree(id=fam$member,dadid=fam$father,momid=fam$mother,sex=fam$sex,affected=fam$affected,famid=fam$pedigree)
  ech.list=list()
  for (i in 1:length(fam.vec))
  {
    avec = ech.ped[i]$id %in% affected
    unav.list = findUnavailable(ech.ped[i],avec)
    ech.list[[i]] = pedigree.trim(unav.list,ech.ped[i])
  }
  

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
  
  genotypes_affected_sub_unique = t(unique(t(genotypes_affected_sub)))
  # Statistiques des génotypes des variants retenus
  dim(genotypes_affected_sub_unique)
  apply(genotypes_affected_sub_unique,2,sum)
  genotypes_affected_RVgene = cbind(fam[fam$affected==2,],genotypes_affected_sub_unique)

  # Tests
  # test = RVgene(genotypes_affected_RVgene,ech.list,sites=1:2,pattern.prob.list = ech.pattern.prob.list,N.list = ech.N.list,type="count")
  # test = RVgene(genotypes_affected_RVgene,ech.list,sites=1:3,pattern.prob.list = ech.pattern.prob.list,N.list = ech.N.list,type="count")
  # test = RVgene(genotypes_affected_RVgene,ech.list,sites=3:4,pattern.prob.list = ech.pattern.prob.list,N.list = ech.N.list,type="count")
  # test = RVgene(genotypes_affected_RVgene,ech.list,sites=5,pattern.prob.list = ech.pattern.prob.list,N.list = ech.N.list,type="count")
  # 
  # table(genotypes_affected_RVgene[,c(1,11)]) 
  
  # Tests partial et all sur des fenêtres de longueur maximale pour le test partiel (avec borne supérieure à lf)
  j=1
  jf = p.vec = pall.vec = numeric(0)
  while (j < ncol(genotypes_affected_RVgene)-5-lf)
  {
    cat(j,"\n")
    lsf = lf
    tf = RVgene(genotypes_affected_RVgene,ech.list,sites=j:(j+lf-1),pattern.prob.list = ech.pattern.prob.list,N.list = ech.N.list,type="count")
    while (is.na(tf$p)&lsf>1)
    {
      lsf = lsf-1
      tf = RVgene(genotypes_affected_RVgene,ech.list,sites=j:(j+lsf-1),pattern.prob.list = ech.pattern.prob.list,N.list = ech.N.list,type="count")
    }
    p.vec = c(p.vec,tf$p)
    pall.vec = c(pall.vec,tf$pall)
    jf = c(jf,j)
    # Incrémente l'indice du variant par la largeur de la fenêtre effectivement testée
    j = j+lsf
  }
# Si j < ncol(genotypes_affected_RVgene)-6, faire une dernière fenêtre de j à ncol(genotypes_affected_RVgene)-6
  if (j < ncol(genotypes_affected_RVgene)-6)
  {
    tf = RVgene(genotypes_affected_RVgene,ech.list,sites=j:(ncol(genotypes_affected_RVgene)-6),pattern.prob.list = ech.pattern.prob.list,N.list = ech.N.list,type="count")
    p.vec = c(p.vec,tf$p)
    pall.vec = c(pall.vec,tf$pall)
    jf = c(jf,j)
  }
  
  # Test all sur fenêtre de largeur lf fixe
  pall_lf.vec = numeric(ceiling((ncol(genotypes_affected_RVgene)-6)/lf))
  for (j in 1: length(pall_lf.vec))
  {
    cat(j,"\n")
    tf = RVgene(genotypes_affected_RVgene,ech.list,sites=((j-1)*lf+1):min(j*lf,ncol(genotypes_affected_RVgene)-6),pattern.prob.list = ech.pattern.prob.list,N.list = ech.N.list,type="count",partial.sharing = F)
    pall_lf.vec[j] = tf$pall
  }
  
min(pall_lf.vec[1:22])*length(pall_lf.vec)  
