filter.rep.n.fam = function(agg_pedfile,variants.by.annot){
  t = agg_pedfile$ped_agg[,-1]
  
  n.unique.fam = sapply(variants.by.annot, function(x){
    subset.variants = which(agg_pedfile$index_variants%in%x)
    
    length(which(rowSums(t[,subset.variants,drop=F]!=0)!=0))
    
    #length(unique(unlist(apply(t[,subset.variants,drop=F], 2, function(y) which(y>0)))))
    
  })
  
  n.unique.fam
}

simulation_null_RetroFunRVS_all = function(ped.mat,ped.trio.list,ped.haplo,f,Z,null_exp_var,nrep)
{
# ped.mat : pedigrees coded as in a ped file
# pheno.ped.list : list of pedigree objects, one object for each pedigree in peds
# ped.trio.list : list of trio objects, one object for each pedigree in peds
# ped.haplo : matrix where each row defines a haplotype, from which to sample the founder haplotypes
# f : vector of family names
# Returns a vector of p-values for Burden and each CRH, the ACAT p for all CRHs and at the end the number of families per CRH
  
#	n is the number of subjects in all pedigrees
n = nrow(ped.mat)
#pedmat = array(NA,c(n,2*ncol(ped.haplo)+6))
#ped.sim.pmpmat = matrix(NA,nrep,ncol(ped.haplo))
pvalues_null = list()

variants.by.annot = apply(Z,2,function(x) which(x!=0))

for (i in 1:nrep)
{
	cat (i,"\n")
ped.sim.haplo = HaploDropSim.fn(ped.mat,ped.trio.list,f,nrow(ped.haplo))
#if (any(ped.sim.haplo == "error")) next

pedmat = data.frame(ped.mat[,1:6],array(NA,c(n,2*ncol(ped.haplo))))
pedmat[,seq(7,ncol(pedmat)-1,2)] = ped.haplo[ped.sim.haplo[,1],]
pedmat[,seq(8,ncol(pedmat),2)] = ped.haplo[ped.sim.haplo[,2],]

  pedfiles_agg_null =  RetroFunRVS::agg.genos.by.fam(pedfile=pedmat,correction="none")
  nfam.vec = filter.rep.n.fam(pedfiles_agg_null, variants.by.annot)
  # On retire les CRHs sans aucun variant
  #Z.tmp = Z[,-which(nfam.vec<=1), drop=FALSE]
  
  #pvalues_null[[i]] = c(RetroFunRVS::RetroFun.RVS(null_exp_var, pedfiles_agg_null, Z.tmp,rep(1,ncol(ped.haplo))),nfam.vec[nfam.vec>1])
  pvalues_null[[i]] = c(RetroFunRVS::RetroFun.RVS(null_exp_var, pedfiles_agg_null, Z,rep(1,ncol(ped.haplo))),nfam.vec)
  
}
pvalues_null
}
