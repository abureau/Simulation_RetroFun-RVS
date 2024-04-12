#!/usr/bin/env Rscript
library(RetroFunRVS)

load("SPAPsimpleprob.RData")
pedfiles_agg_null = readRDS("pedfiles_agg_null_10000reps.RDS")

n.unique.config.by.fam = sapply(forsim.N.list, unique)
n.unique.config.by.fam = lapply(n.unique.config.by.fam, function(x) sort(x))

prob.sharing.by.famid = lapply(1:length(forsim.pattern.prob.list), function(x) tapply(forsim.pattern.prob.list[[x]], forsim.N.list[[x]], sum))
names(prob.sharing.by.famid) = names(forsim.pattern.prob.list)
prob.sharing.by.famid = lapply(prob.sharing.by.famid, function(x) x[order(as.numeric(names(x)))])

pvalues = readRDS("pvalues_null_Genes.RDS")
Z = read.table("Annot_Z_genes.mat")


log.pvalues = -log10(sapply(1:10000, function(x) pvalues[[x]]$ACAT))

index.with.highly.inflated.signal = which(log.pvalues>=3)
index.with.inflated.signal = which(log.pvalues>=2&log.pvalues<3)

null_without_consanguinity = read.table("null_without_consanguinity.txt", header=TRUE)

print("A")
print(index.with.inflated.signal)

resample.genos.with.inflated.signal = lapply(index.with.inflated.signal, function(x) {
  print(x)
  replicate(1000, resample.genos.by.fam(pedfiles_agg_null[[x]], n.unique.config.by.fam ,prob.sharing.by.famid))
})


bootstrap.burden.with.inflated.signal = lapply(1:length(resample.genos.with.inflated.signal), function(x){
	print(x)
	apply(resample.genos.with.inflated.signal[[x]], 2, function(x){
		compute.Burden.by.Annot(null_without_consanguinity, x, Z, W = rep(1,510))
	})
})


l=list("index.inflated.signal" = index.with.inflated.signal, "bootstrap.genotypes"=resample.genos.with.inflated.signal ,"boot.Burden"=bootstrap.burden.with.inflated.signal)

saveRDS(l,"/home/loicm/scratch/bootstrap_inflated_signal_Genes.RDS")

print("B")
print(index.with.highly.inflated.signal)

resample.genos.with.highly.inflated.signal = lapply(index.with.highly.inflated.signal, function(x) {
  replicate(10000, resample.genos.by.fam(pedfiles_agg_null[[x]], n.unique.config.by.fam ,prob.sharing.by.famid))
})

bootstrap.burden.with.highly.inflated.signal = lapply(1:length(resample.genos.with.highly.inflated.signal), function(x){
        print(x)
        apply(resample.genos.with.highly.inflated.signal[[x]], 2, function(x){
                compute.Burden.by.Annot(null_without_consanguinity, x, Z, W = rep(1,510))
        })
})


l=list("index.inflated.signal" = index.with.highly.inflated.signal, "bootstrap.genotypes"=resample.genos.with.highly.inflated.signal ,"boot.Burden"=bootstrap.burden.with.highly.inflated.signal)
saveRDS(l,"/home/loicm/scratch/bootstrap_highly_inflated_signal_Genes.RDS")

