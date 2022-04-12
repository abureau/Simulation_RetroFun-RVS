#Type I error rate
l_CRHs_lee = list()
l_quant_lee = list()
l_quant_CRHs_lee = list()


ACAT_Fisher_null = foreach(f=1:length(ped_files_null))%dopar%{
  
  library(snpStats)
  library(dplyr)
  library(plyr)
  library(ACAT)
  
  agg = aggregate.geno.by.fam(ped_files_null[f], correction = "replace.homo", replace_ind_geno = T )
  
  Z_CRHs_sub = Z_CRHs[agg$index_variants,]
  Z_quant_sub = Z_quant[agg$index_variants,]
  Z_quant_CRHs_sub = Z_quant_CRHs[agg$index_variants,]
  
  retrofun_CRHS = retrofun.RVS(null,agg,Z_CRHs_sub, adjust.low.p = F)[,c("ACAT", "Fisher")]
  retrofun_quant = retrofun.RVS(null,agg,Z_quant_sub,adjust.low.p = F)[,c("ACAT","Fisher")]
  retrofun_quant_CRHs = retrofun.RVS(null,agg,Z_quant_CRHs_sub, adjust.low.p = F)[,c("ACAT", "Fisher")]
  
  retrofun_CRHS_lee = retrofun.RVS(null,agg,Z_CRHs_sub, adjust.low.p = T, p.threshold = 0.001, n.Boot=1000, method.adjust="Lee.adjust")[,c("ACAT", "Fisher")]
  retrofun_quant_lee = retrofun.RVS(null,agg,Z_quant_sub,adjust.low.p = T, p.threshold = 0.001, n.Boot=1000, method.adjust="Lee.adjust")[,c("ACAT", "Fisher")]
  retrofun_quant_CRHs_lee = retrofun.RVS(null,agg,Z_quant_CRHs_sub, adjust.low.p = T, p.threshold = 0.001, n.Boot=1000, method.adjust="Lee.adjust")[,c("ACAT", "Fisher")]
  
  retrofun_CRHS_boot = retrofun.RVS(null,agg,Z_CRHs_sub, adjust.low.p = T, p.threshold = 0.001, n.Boot=1000, method.adjust="Bootstrap")[,c("ACAT", "Fisher")]
  retrofun_quant_boot = retrofun.RVS(null,agg,Z_quant_sub,adjust.low.p = T, p.threshold = 0.001, n.Boot=1000, method.adjust="Bootstrap")[,c("ACAT", "Fisher")]
  retrofun_quant_CRHs_boot = retrofun.RVS(null,agg,Z_quant_CRHs_sub, adjust.low.p = T, p.threshold = 0.001, n.Boot=1000, method.adjust="Bootstrap")[,c("ACAT", "Fisher")]
  
  return(list("Init"=list("CHRs"=retrofun_CRHS,"Quant"=retrofun_quant,"CRHs_Quant"=retrofun_quant_CRHs),"Lee"=list("CRHs"=retrofun_CRHS_lee,"Quant"=retrofun_quant_lee,"CRHs_Quant"=retrofun_quant_CRHs_lee),"Boot"=list("CRHs"=retrofun_CRHS_boot,"Quant"=retrofun_quant_boot,"CRHs_Quant"=retrofun_quant_CRHs_boot)))
}

stopCluster(cl)
ACAT_Fisher_sample[,1]$result.1

hist(unlist(ACAT_Fisher_sample[,5]))

saveRDS(l_CRHs, "Burden_CRHs_10000rep.RDS")
saveRDS(l_quant, "Burden_Quant_1_10000rep.RDS")

ACAT_null_quant = sapply(l_quant, function(x) x$ACAT)
ACAT_null_CRHs = sapply(l_CRHs, function(x) x$ACAT)

qqunif.plot(ACAT_null_quant)
qqunif.plot(ACAT_null_CRHs)

sum(ACAT_null_quant<=0.05)/1000
sum(ACAT_null_quant<=0.01)/1000
sum(ACAT_null_quant<=0.005)/1000
sum(ACAT_null_quant<=0.001)/1000
sum(ACAT_null_quant<=0.00001)/1000

sum(ACAT_null_CRHs<=0.05)/1000
sum(ACAT_null_CRHs<=0.01)/1000
sum(ACAT_null_CRHs<=0.005)/1000
sum(ACAT_null_CRHs<=0.001)/1000
sum(ACAT_null_CRHs<=0.0001)/1000


ACAT_null_quant_lee =sapply(l_quant_lee, function(x) x$ACAT)
ACAT_null_CRHs_lee = sapply(l_CRHs_lee, function(x) x$ACAT)
ACAT_null_quant_CRHs_lee = sapply(l_quant_CRHs_lee, function(x) x$ACAT)

sum(ACAT_null_quant_lee<=0.05)/10000
sum(ACAT_null_quant_lee<=0.01)/10000
sum(ACAT_null_quant_lee<=0.005)/10000
sum(ACAT_null_quant_lee<=0.001)/10000
sum(ACAT_null_quant_lee<=0.00001)/10000

sum(ACAT_null_CRHs_lee<=0.05)/10000
sum(ACAT_null_CRHs_lee<=0.01)/10000
sum(ACAT_null_CRHs_lee<=0.005)/10000
sum(ACAT_null_CRHs_lee<=0.001)/10000
sum(ACAT_null_CRHs_lee<=0.00001)/10000

sum(ACAT_null_quant_CRHs_lee<=0.05)/10000
sum(ACAT_null_quant_CRHs_lee<=0.01)/10000
sum(ACAT_null_quant_CRHs_lee<=0.005)/10000
sum(ACAT_null_quant_CRHs_lee<=0.001)/10000
sum(ACAT_null_quant_CRHs_lee<=0.00001)/10000