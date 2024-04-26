library(snpStats)
library(dplyr)
library(plyr)
library(RVS)
library(kinship2)
load("../SPAPsimpleprob.RData")
load("../pedsimple.RData")

correction="replace.homo"

# Largeur de fenêtre
lf = 5

ped_files_OR5_100.01 = list.files("../1causal_OR5_RVS/100causal/", full.names=T, recursive=T)
ped_files_OR5_75.01 = list.files("../1causal_OR5_RVS/75causal/", full.names = T, recursive=T)
ped_files_OR5_50.01 = list.files("../1causal_OR5_RVS/50causal/", full.names = T, recursive = T)

# Test avec un seul réplicat
#pedfile = ped_files_OR5_75[1]
#pedfile = ped_files_OR5_50[1]
#pedfile = ped_files_OR5_100.01[1]

# Boucle sur les réplicats
p50.OR5.01.list = pall50.OR5.01.list = jf50.OR5.01.list = pall_lf50.OR5.01.list = list()

for (r in 1:length(ped_files_OR5_50.01))
#for (r in 1:200)
  {
  cat("r=",r,"\n")
  pedfile = ped_files_OR5_50.01[r]

  sample = read.pedfile(pedfile)
  fam = sample$fam
  fam$affected[is.na(fam$affected)] = 1
  affected = fam[fam$affected==2,"member"]

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
p50.OR5.01.list[[r]] = p.vec
pall50.OR5.01.list[[r]] = pall.vec
pall_lf50.OR5.01.list[[r]] = pall_lf.vec
jf50.OR5.01.list[[r]] = jf
#min(pall_lf.vec)*length(pall_lf.vec)  
}

#save(p50.OR5.01.list,pall50.OR5.01.list,pall_lf50.OR5.01.list,jf50.OR5.01.list,file="OR5_50.01_v2.RData")
