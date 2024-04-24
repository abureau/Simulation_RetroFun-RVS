library(dplyr)
library(plyr)
library(RVS)
library(kinship2)
load("SPAPsimpleprob.RData")
load("pedsimple.RData")

correction="replace.homo"

# Largeur de fenêtre
lf = 5


# Test avec un seul réplicat
pedfile = "rep1.ped"


  p = read.table(pedfile, header = F)
  fam = p[,1:6]
  names(fam) = c("pedigree","member","father","mother","sex","affected")
  fam$father[fam$father==0] = NA
  fam$mother[fam$mother==0] = NA
  
  affected = fam[fam$affected==2,"member"]
  genos = p[,7:ncol(p)]
  
  genos[genos==1] = 0
  genos[genos==2] = 1
  
  genotypes_all =  data.frame(do.call("cbind",lapply(seq(2, ncol(genos),2), function(x){
    rowSums(genos[,c(x-1,x)])
  })))
  
    fam$affected[is.na(fam$affected)] = 1

  # Créer pattern.prob.list et N.list
  fam.vec = as.character(unique(fam$pedigree))
  # Les familles à partir de 41 ont une lettre ajoutée comme préfixe à leur nom
  fam.type = c(fam.vec[1:40],substring(fam.vec[41:length(fam.vec)],2))
  ech.pattern.prob.list = forsim.pattern.prob.list[fam.type]
  ech.N.list = forsim.N.list[fam.type]
  names(ech.pattern.prob.list) = fam.vec
  names(ech.N.list) = fam.vec
  
  # Création d'un objet pedigreeList, puis élagage des familles et sauvegarde dans une liste d'objets pedigree
  affected01 = ifelse(is.na(fam$affected),0,fam$affected-1)
  table(affected01)
  ech.ped = pedigree(id=fam$member,dadid=fam$father,momid=fam$mother,sex=fam$sex,affected=affected01,famid=fam$pedigree)
  ech.list=list()
  for (i in 1:length(fam.vec))
  {
    avec = ech.ped[i]$id %in% affected
    unav.list = findUnavailable(ech.ped[i],avec)
    ech.list[[i]] = pedigree.trim(unav.list,ech.ped[i])
  }
  
  rownames(genotypes_all) = fam$member

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
    break;
  }
  # Remove duplicated variants
  #genotypes_affected_unique = t(unique(t(genotypes_affected)))
  dup = duplicated(t(genotypes_affected))
  genotypes_affected_unique = genotypes_affected[,!dup]
  
  if(correction=="remove.homo"){
    #Keep only heterozyguous variants
    hetero_variants = which(colSums(genotypes_affected_unique==2) == 0)
    genotypes_affected_sub_unique = select(genotypes_affected_unique, all_of(hetero_variants))
    
  } else if(correction=="replace.homo"){
    #Replace homozygous variants by heterogygous ones
    genotypes_affected_sub_unique = genotypes_affected_unique
    genotypes_affected_sub_unique[genotypes_affected_sub_unique>1] = 1 
  } else{
    stop("Please choose a valid option for correction parameter")
  }
  # Statistiques des génotypes des variants retenus
  #dim(genotypes_affected_sub_unique)
  #apply(genotypes_affected_sub_unique,2,sum)
  genotypes_affected_RVgene = cbind(fam[affected01==1,],genotypes_affected_sub_unique)

  print("Time for partial and complete sharing:")
  start = Sys.time()
  
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
  end = Sys.time()
  time.both = end-start
  print(time.both)
  
  print("Time for complete sharing on windows of 5 SNPs:")
  start = Sys.time()
  
  # Test all sur fenêtre de largeur lf fixe
  pall_lf.vec = numeric(ceiling((ncol(genotypes_affected_RVgene)-6)/lf))
  for (j in 1: length(pall_lf.vec))
  {
    cat(j,"\n")
    tf = RVgene(genotypes_affected_RVgene,ech.list,sites=((j-1)*lf+1):min(j*lf,ncol(genotypes_affected_RVgene)-6),pattern.prob.list = ech.pattern.prob.list,N.list = ech.N.list,type="count",partial.sharing = F)
    pall_lf.vec[j] = tf$pall
  }
end = Sys.time()
time.complete = end-start
print(time.complete)