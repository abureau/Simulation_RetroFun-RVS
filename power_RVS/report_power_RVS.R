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
# 2% of causal variants

load("alter100.RData")
load("alter75.RData")
load("alter50.RData")

# Close to 100% causal variants in CRH 1
minp100.2pc.vec = sapply(p100.02.list,minp_Bonf)
minpall100.2pc.vec = sapply(pall100.02.list,minp_Bonf)
minpall_lf100.2pc.vec = sapply(pall_lf100.02.list,minp_Bonf)
#power100.2pc.05 = c(sum(minp100.2pc.vec<0.05),sum(minpall100.2pc.vec<0.05),sum(minpall_lf100.2pc.vec<0.05))
#power100.2pc.3000 = c(sum(minp100.2pc.vec<0.05/3000),sum(minpall100.2pc.vec<0.05/3000),sum(minpall_lf100.2pc.vec<0.05/3000))

ACATp100.2pc.vec = sapply(p100.02.list,ACAT.fct)
ACATpall100.2pc.vec = sapply(pall100.02.list,ACAT.fct)
ACATpall_lf100.2pc.vec = sapply(pall_lf100.02.list,ACAT.fct)
power100.2pc.05 = c(sum(ACATp100.2pc.vec<0.05),sum(ACATpall100.2pc.vec<0.05),sum(ACATpall_lf100.2pc.vec<0.05))
power100.2pc.3000 = c(sum(ACATp100.2pc.vec<0.05/3000),sum(ACATpall100.2pc.vec<0.05/3000),sum(ACATpall_lf100.2pc.vec<0.05/3000))


# 75% causal variants in CRH 1
minp75.2pc.vec = sapply(p75.02.list,minp_Bonf)
minpall75.2pc.vec = sapply(pall75.02.list,minp_Bonf)
minpall_lf75.2pc.vec = sapply(pall_lf75.02.list,minp_Bonf)
#power75.2pc.05 = c(sum(minp75.2pc.vec<0.05),sum(minpall75.2pc.vec<0.05),sum(minpall_lf75.2pc.vec<0.05))
#power75.2pc.3000 = c(sum(minp75.2pc.vec<0.05/3000),sum(minpall75.2pc.vec<0.05/3000),sum(minpall_lf75.2pc.vec<0.05/3000))

ACATp75.2pc.vec = sapply(p75.02.list,ACAT.fct)
ACATpall75.2pc.vec = sapply(pall75.02.list,ACAT.fct)
ACATpall_lf75.2pc.vec = sapply(pall_lf75.02.list,ACAT.fct)
power75.2pc.05 = c(sum(ACATp75.2pc.vec<0.05),sum(ACATpall75.2pc.vec<0.05),sum(ACATpall_lf75.2pc.vec<0.05))
power75.2pc.3000 = c(sum(ACATp75.2pc.vec<0.05/3000),sum(ACATpall75.2pc.vec<0.05/3000),sum(ACATpall_lf75.2pc.vec<0.05/3000))

# 50% causal variants in CRH 1
minp50.2pc.vec = sapply(p50.02.list,minp_Bonf)
minpall50.2pc.vec = sapply(pall50.02.list,minp_Bonf)
minpall_lf50.2pc.vec = sapply(pall_lf50.02.list,minp_Bonf)
#power50.2pc.05 = c(sum(minp50.2pc.vec<0.05),sum(minpall50.2pc.vec<0.05),sum(minpall_lf50.2pc.vec<0.05))
#power50.2pc.00005 = c(sum(minp50.2pc.vec<0.00005),sum(minpall50.2pc.vec<0.00005),sum(minpall_lf50.2pc.vec<0.00005))

ACATp50.2pc.vec = sapply(p50.02.list,ACAT.fct)
ACATpall50.2pc.vec = sapply(pall50.02.list,ACAT.fct)
ACATpall_lf50.2pc.vec = sapply(pall_lf50.02.list,ACAT.fct)
power50.2pc.05 = c(sum(ACATp50.2pc.vec<0.05),sum(ACATpall50.2pc.vec<0.05),sum(ACATpall_lf50.2pc.vec<0.05))
power50.2pc.3000 = c(sum(ACATp50.2pc.vec<0.05/3000),sum(ACATpall50.2pc.vec<0.05/3000),sum(ACATpall_lf50.2pc.vec<0.05/3000))


# Examine number of windows for variable length
l100.2pc = sapply(p100.02.list,length)
l75.2pc = sapply(p75.02.list,length)
l50.2pc = sapply(p50.02.list,length)
summary(l100.2pc)
summary(l75.2pc)
summary(l50.2pc)

# 4% of causal variants
load("alter100.04.RData")
load("alter75.04.RData")
load("alter50.04.RData")

# Close to 100% causal variants in CRH 1
minp100.4pc.vec = sapply(p100.04.list,minp_Bonf)
minpall100.4pc.vec = sapply(pall100.04.list,minp_Bonf)
minpall_lf100.4pc.vec = sapply(pall_lf100.04.list,minp_Bonf)
power100.4pc.05 = c(sum(minp100.4pc.vec<0.05),sum(minpall100.4pc.vec<0.05),sum(minpall_lf100.4pc.vec<0.05))
power100.4pc.00005 = c(sum(minp100.4pc.vec<0.00005),sum(minpall100.4pc.vec<0.00005),sum(minpall_lf100.4pc.vec<0.00005))

# 75% causal variants in CRH 1
minp75.4pc.vec = sapply(p75.04.list,minp_Bonf)
minpall75.4pc.vec = sapply(pall75.04.list,minp_Bonf)
minpall_lf75.4pc.vec = sapply(pall_lf75.04.list,minp_Bonf)
power75.4pc.05 = c(sum(minp75.4pc.vec<0.05),sum(minpall75.4pc.vec<0.05),sum(minpall_lf75.4pc.vec<0.05))
power75.4pc.00005 = c(sum(minp75.4pc.vec<0.00005),sum(minpall75.4pc.vec<0.00005),sum(minpall_lf75.4pc.vec<0.00005))

# 50% causal variants in CRH 1
minp50.4pc.vec = sapply(p50.04.list,minp_Bonf)
minpall50.4pc.vec = sapply(pall50.04.list,minp_Bonf)
minpall_lf50.4pc.vec = sapply(pall_lf50.04.list,minp_Bonf)
power50.4pc.05 = c(sum(minp50.4pc.vec<0.05),sum(minpall50.4pc.vec<0.05),sum(minpall_lf50.4pc.vec<0.05))
power50.4pc.00005 = c(sum(minp50.4pc.vec<0.00005),sum(minpall50.4pc.vec<0.00005),sum(minpall_lf50.4pc.vec<0.00005))

#power.05 = cbind(power100.2pc.05,power75.2pc.05,power50.2pc.05,power100.4pc.05,power75.4pc.05,power50.4pc.05)
#power.3000 = cbind(power100.2pc.3000,power75.2pc.3000,power50.2pc.3000,power100.4pc.3000,power75.4pc.3000,power50.4pc.00005)
power.05 = cbind(power100.2pc.05,power75.2pc.05,power50.2pc.05)
power.3000 = cbind(power100.2pc.3000,power75.2pc.3000,power50.2pc.3000)
row.names(power.3000) = row.names(power.05) = c("Partial sharing","Complete sharing","Complete sharing, fixed")
power.05
power.3000

# Examine number of windows for variable length
l100.4pc = sapply(p100.04.list,length)
l75.4pc = sapply(p75.04.list,length)
l50.4pc = sapply(p50.04.list,length)
summary(l100.4pc)
summary(l75.4pc)
summary(l50.4pc)
