library(kinship2)
library(dplyr)
library(RVsharing)

library(foreach)
library(doParallel)

# Fichiers de code à charger
source("DropHaploLocus.R")
source("drophaplolocusfn.R")
source("ped2haplo.R")
source("simulation_null_RetroFunRVS_allCRHs.R")
load("pedsimple.RData")

set.seed(100)

pedData=read.table("pedigree_for_sim_clean_trim_2.ped")
sum(pedData[,6]==2)
ped.frame = pedData

TG.ped= read.table("variants_subset_TAD.ped") #509 variants. vcf_to_ped en retire un, on ne sait pourquoi. C,est le rs199753143
TG.ped[TG.ped==3]=2
TG.haplo = ped2haplo(TG.ped)
dim(TG.haplo)

# Ici on extrait les noms de familles d'un data.frame d'un fichier ped.
fams.simple = unique(as.character(pedData[,1]))
length(fams.simple)

# Création d'un objet pedigree
pedData[pedData[,3]==0,3]=NA
pedData[pedData[,4]==0,4]=NA
pedData[,6]=pedData[,6]-1
pedsimple = pedigree(pedData[,2], pedData[,3], pedData[,4], pedData[,5], pedData[,6], famid=pedData[,1])
GC.aff = pedData[pedData[,6]==1,2]
  
pedsimple.list = list()
for(j in 1:length(fams.simple)){
  avec = pedsimple[j]$id %in% GC.aff
  unav.list = findUnavailable(pedsimple[j],avec)
  pedsimple.list[[j]] = pedigree.trim(unav.list,pedsimple[j])
}
names(pedsimple.list) = fams.simple
ID.PedForAnalysis = unlist(sapply(pedsimple.list,function(obj) obj$id))

ped.frame = data.frame(ped.frame[ped.frame[,2] %in% ID.PedForAnalysis,])
# Ici il faut mettre les IDs des parents absents de la famille à 0
#ped.frame[!(ped.frame[,3]%in%ped.frame[,2]),3:4] = 0

# Vérification que les noms de familles sont les bons
all(fams.simple==unique(as.character(pedsimple$famid)))


# Converting pedigree objects into trio objects
trioforsim = list()
for (f in names(pedsimple.list))
{
tmp1 = RVsharing::ped2trio(pedsimple.list[[f]]); 
cat(f,length(tmp1$object),"\n")
trioforsim[f] = tmp1$object
}

Z = read.table("data/Annot_Z_CRHs.mat", header=F)
#On doit retirer rs199753143, le 320e variant dans le vcf initial.
Z <- Z[-320,] 
null_without_consanguinity = read.table("data/null_without_consanguinity.txt", header=TRUE)

# Cette étape ne devrait pas être requise
null_without_consanguinity_sim = rbind(null_without_consanguinity,cbind(FamID=paste0("A",c(105,115,125,182,211,220)),null_without_consanguinity[null_without_consanguinity$FamID%in%c(105,115,125,182,211,220),-1]),cbind(FamID=paste0("B",c(105,115,125,182,211,220)),null_without_consanguinity[null_without_consanguinity$FamID%in%c(105,115,125,182,211,220),-1]))

#Retirer familles en trop
null_without_consanguinity_sim$FamID = as.character(null_without_consanguinity_sim$FamID)
keep <- c("101", "105", "108", "110", "111", "115", "119", "121", "122", "123", "124", "125", "126", "129", "133", "134", "182", "207", "208", "211", "211A", "212", "212A", "217", "218", "220", "223", "226", "227", "228", "230", "233", "234", "235", "236", "238", "239", "245", "250", "A105", "A115", "A125", "A182", "A211", "A220", "B105", "B115", "B125", "B182", "B211", "B220", "S2029")
null_without_consanguinity_sim <- null_without_consanguinity_sim[null_without_consanguinity_sim$FamID %in% keep,]

#Rouler
ncores=32
doParallel::registerDoParallel(cores=ncores)
foreach(i = 1:ncores) %dopar% {
  pvalues_null = simulation_null_RetroFunRVS_all(ped.frame, trioforsim, TG.haplo, fams.simple, Z, null_without_consanguinity_sim, nrep = 156250)
  saveRDS(pvalues_null, paste0("sim_object_arcturus_", i, ".rds"))
}


