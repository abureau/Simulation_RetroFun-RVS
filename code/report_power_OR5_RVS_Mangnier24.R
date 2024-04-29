# Code to compute p-values of the RVS tests from the objects containing the test results

library(ACAT)
minp_Bonf = function(vec)
{
	min(vec,na.rm=T)*length(vec)
}
ACAT.fct = function(vec)
{
  vec[is.na(vec)] = 1
  vec = ifelse(vec>0.999999,1-1/length(vec),vec)
  ACAT(vec)
}

load("data\\powerRVS\\OR5_100.02.RData")
load("data\\powerRVS\\OR5_75.02.RData")
load("data\\powerRVS\\OR5_50.02.RData")
load("data\\powerRVS\\OR5_100.01.RData")
load("data\\powerRVS\\OR5_75.01.RData")
load("data\\powerRVS\\OR5_50.01.RData")

# 100% causal variants in CRH 1

# 2% of causal variants
minp100.2pc.vec = sapply(p100.OR5.02.list,minp_Bonf)
minpall100.2pc.vec = sapply(pall100.OR5.02.list,minp_Bonf)
minpall_lf100.2pc.vec = sapply(pall_lf100.OR5.02.list,minp_Bonf)

ACATp100.2pc.vec = sapply(p100.OR5.02.list,ACAT.fct)
ACATpall100.2pc.vec = sapply(pall100.OR5.02.list,ACAT.fct)
ACATpall_lf100.2pc.vec = sapply(pall_lf100.OR5.02.list,ACAT.fct)
power100.2pc.05 = c(sum(ACATp100.2pc.vec<0.05),sum(ACATpall100.2pc.vec<0.05),sum(ACATpall_lf100.2pc.vec<0.05))/length(ACATp100.2pc.vec)
power100.2pc.6000 = c(sum(ACATp100.2pc.vec<0.05/6000),sum(ACATpall100.2pc.vec<0.05/6000),sum(ACATpall_lf100.2pc.vec<0.05/6000))/length(ACATp100.2pc.vec)

# 1% of causal variants
minp100.1pc.vec = sapply(p100.OR5.01.list,minp_Bonf)
minpall100.1pc.vec = sapply(pall100.OR5.01.list,minp_Bonf)
minpall_lf100.1pc.vec = sapply(pall_lf100.OR5.01.list,minp_Bonf)

ACATp100.1pc.vec = sapply(p100.OR5.01.list,ACAT.fct)
ACATpall100.1pc.vec = sapply(pall100.OR5.01.list,ACAT.fct)
ACATpall_lf100.1pc.vec = sapply(pall_lf100.OR5.01.list,ACAT.fct)
power100.1pc.05 = c(sum(ACATp100.1pc.vec<0.05),sum(ACATpall100.1pc.vec<0.05),sum(ACATpall_lf100.1pc.vec<0.05))/length(ACATp100.1pc.vec)
power100.1pc.6000 = c(sum(ACATp100.1pc.vec<0.05/6000),sum(ACATpall100.1pc.vec<0.05/6000),sum(ACATpall_lf100.1pc.vec<0.05/6000))/length(ACATp100.1pc.vec)

# 75% causal variants in CRH 1

# 2% of causal variants
minp75.2pc.vec = sapply(p75.OR5.02.list,minp_Bonf)
minpall75.2pc.vec = sapply(pall75.OR5.02.list,minp_Bonf)
minpall_lf75.2pc.vec = sapply(pall_lf75.OR5.02.list,minp_Bonf)

ACATp75.2pc.vec = sapply(p75.OR5.02.list,ACAT.fct)
ACATpall75.2pc.vec = sapply(pall75.OR5.02.list,ACAT.fct)
ACATpall_lf75.2pc.vec = sapply(pall_lf75.OR5.02.list,ACAT.fct)
power75.2pc.05 = c(sum(ACATp75.2pc.vec<0.05),sum(ACATpall75.2pc.vec<0.05),sum(ACATpall_lf75.2pc.vec<0.05))/length(ACATp75.2pc.vec)
power75.2pc.6000 = c(sum(ACATp75.2pc.vec<0.05/6000),sum(ACATpall75.2pc.vec<0.05/6000),sum(ACATpall_lf75.2pc.vec<0.05/3000))/length(ACATp75.2pc.vec)

# 1% of causal variants
minp75.1pc.vec = sapply(p75.OR5.01.list,minp_Bonf)
minpall75.1pc.vec = sapply(pall75.OR5.01.list,minp_Bonf)
minpall_lf75.1pc.vec = sapply(pall_lf75.OR5.01.list,minp_Bonf)

ACATp75.1pc.vec = sapply(p75.OR5.01.list,ACAT.fct)
ACATpall75.1pc.vec = sapply(pall75.OR5.01.list,ACAT.fct)
ACATpall_lf75.1pc.vec = sapply(pall_lf75.OR5.01.list,ACAT.fct)
power75.1pc.05 = c(sum(ACATp75.1pc.vec<0.05),sum(ACATpall75.1pc.vec<0.05),sum(ACATpall_lf75.1pc.vec<0.05))/length(ACATp75.1pc.vec)
power75.1pc.6000 = c(sum(ACATp75.1pc.vec<0.05/6000),sum(ACATpall75.1pc.vec<0.05/6000),sum(ACATpall_lf75.1pc.vec<0.05/6000))/length(ACATp75.1pc.vec)

# 50% causal variants in CRH 1

# 2% of causal variants
minp50.2pc.vec = sapply(p50.OR5.02.list,minp_Bonf)
minpall50.2pc.vec = sapply(pall50.OR5.02.list,minp_Bonf)
minpall_lf50.2pc.vec = sapply(pall_lf50.OR5.02.list,minp_Bonf)

ACATp50.2pc.vec = sapply(p50.OR5.02.list,ACAT.fct)
ACATpall50.2pc.vec = sapply(pall50.OR5.02.list,ACAT.fct)
ACATpall_lf50.2pc.vec = sapply(pall_lf50.OR5.02.list,ACAT.fct)
power50.2pc.05 = c(sum(ACATp50.2pc.vec<0.05),sum(ACATpall50.2pc.vec<0.05),sum(ACATpall_lf50.2pc.vec<0.05))/length(ACATp50.2pc.vec)
power50.2pc.6000 = c(sum(ACATp50.2pc.vec<0.05/6000),sum(ACATpall50.2pc.vec<0.05/6000),sum(ACATpall_lf50.2pc.vec<0.05/6000))/length(ACATp50.2pc.vec)

# 1% of causal variants
minp50.1pc.vec = sapply(p50.OR5.01.list,minp_Bonf)
minpall50.1pc.vec = sapply(pall50.OR5.01.list,minp_Bonf)
minpall_lf50.1pc.vec = sapply(pall_lf50.OR5.01.list,minp_Bonf)

ACATp50.1pc.vec = sapply(p50.OR5.01.list,ACAT.fct)
ACATpall50.1pc.vec = sapply(pall50.OR5.01.list,ACAT.fct)
ACATpall_lf50.1pc.vec = sapply(pall_lf50.OR5.01.list,ACAT.fct)
power50.1pc.05 = c(sum(ACATp50.1pc.vec<0.05),sum(ACATpall50.1pc.vec<0.05),sum(ACATpall_lf50.1pc.vec<0.05))/length(ACATp50.1pc.vec)
power50.1pc.6000 = c(sum(ACATp50.1pc.vec<0.05/6000),sum(ACATpall50.1pc.vec<0.05/6000),sum(ACATpall_lf50.1pc.vec<0.05/6000))/length(ACATp50.1pc.vec)

power.05 = cbind(power100.2pc.05,power75.2pc.05,power50.2pc.05)
power.6000 = cbind(power100.2pc.6000,power75.2pc.6000,power50.2pc.6000)
row.names(power.6000) = row.names(power.05) = c("Partial sharing","Complete sharing","Complete sharing, fixed")
power.05
power.6000

power.1pc.05 = cbind(power100.1pc.05,power75.1pc.05,power50.1pc.05)
power.1pc.6000 = cbind(power100.1pc.6000,power75.1pc.6000,power50.1pc.6000)
row.names(power.1pc.6000) = row.names(power.1pc.05) = c("Partial sharing","Complete sharing","Complete sharing, fixed")
power.1pc.05
power.1pc.6000


# Examine number of windows for variable length
l100.2pc = sapply(p100.02.list,length)
l75.2pc = sapply(p75.02.list,length)
l50.2pc = sapply(p50.02.list,length)
summary(l100.2pc)
summary(l75.2pc)
summary(l50.2pc)

