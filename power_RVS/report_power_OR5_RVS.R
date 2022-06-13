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

load("OR5_100.02.RData")
load("OR5_75.02.RData")
load("OR5_50.02.RData")
load("OR5_100.01.RData")
load("OR5_75.01.RData")
load("OR5_50.01.RData")

# 100% causal variants in CRH 1

# 2% of causal variants
minp100.2pc.vec = sapply(p100.OR5.02.list,minp_Bonf)
minpall100.2pc.vec = sapply(pall100.OR5.02.list,minp_Bonf)
minpall_lf100.2pc.vec = sapply(pall_lf100.OR5.02.list,minp_Bonf)
#power100.2pc.05 = c(sum(minp100.2pc.vec<0.05),sum(minpall100.2pc.vec<0.05),sum(minpall_lf100.2pc.vec<0.05))
#power100.2pc.3000 = c(sum(minp100.2pc.vec<0.05/3000),sum(minpall100.2pc.vec<0.05/3000),sum(minpall_lf100.2pc.vec<0.05/3000))

ACATp100.2pc.vec = sapply(p100.OR5.02.list,ACAT.fct)
ACATpall100.2pc.vec = sapply(pall100.OR5.02.list,ACAT.fct)
ACATpall_lf100.2pc.vec = sapply(pall_lf100.OR5.02.list,ACAT.fct)
power100.2pc.05 = c(sum(ACATp100.2pc.vec<0.05),sum(ACATpall100.2pc.vec<0.05),sum(ACATpall_lf100.2pc.vec<0.05))/length(ACATp100.2pc.vec)
power100.2pc.3000 = c(sum(ACATp100.2pc.vec<0.05/3000),sum(ACATpall100.2pc.vec<0.05/3000),sum(ACATpall_lf100.2pc.vec<0.05/3000))/length(ACATp100.2pc.vec)
power100.2pc.6000 = c(sum(ACATp100.2pc.vec<0.05/6000),sum(ACATpall100.2pc.vec<0.05/6000),sum(ACATpall_lf100.2pc.vec<0.05/6000))/length(ACATp100.2pc.vec)

# 1% of causal variants
minp100.1pc.vec = sapply(p100.OR5.01.list,minp_Bonf)
minpall100.1pc.vec = sapply(pall100.OR5.01.list,minp_Bonf)
minpall_lf100.1pc.vec = sapply(pall_lf100.OR5.01.list,minp_Bonf)
#power100.1pc.05 = c(sum(minp100.1pc.vec<0.05),sum(minpall100.1pc.vec<0.05),sum(minpall_lf100.1pc.vec<0.05))
#power100.1pc.3000 = c(sum(minp100.1pc.vec<0.05/3000),sum(minpall100.1pc.vec<0.05/3000),sum(minpall_lf100.1pc.vec<0.05/3000))

ACATp100.1pc.vec = sapply(p100.OR5.01.list,ACAT.fct)
ACATpall100.1pc.vec = sapply(pall100.OR5.01.list,ACAT.fct)
ACATpall_lf100.1pc.vec = sapply(pall_lf100.OR5.01.list,ACAT.fct)
power100.1pc.05 = c(sum(ACATp100.1pc.vec<0.05),sum(ACATpall100.1pc.vec<0.05),sum(ACATpall_lf100.1pc.vec<0.05))/length(ACATp100.1pc.vec)
power100.1pc.3000 = c(sum(ACATp100.1pc.vec<0.05/3000),sum(ACATpall100.1pc.vec<0.05/3000),sum(ACATpall_lf100.1pc.vec<0.05/3000))/length(ACATp100.1pc.vec)
power100.1pc.6000 = c(sum(ACATp100.1pc.vec<0.05/6000),sum(ACATpall100.1pc.vec<0.05/6000),sum(ACATpall_lf100.1pc.vec<0.05/6000))/length(ACATp100.1pc.vec)

# 75% causal variants in CRH 1

# 2% of causal variants
minp75.2pc.vec = sapply(p75.OR5.02.list,minp_Bonf)
minpall75.2pc.vec = sapply(pall75.OR5.02.list,minp_Bonf)
minpall_lf75.2pc.vec = sapply(pall_lf75.OR5.02.list,minp_Bonf)
#power75.2pc.05 = c(sum(minp75.2pc.vec<0.05),sum(minpall75.2pc.vec<0.05),sum(minpall_lf75.2pc.vec<0.05))
#power75.2pc.3000 = c(sum(minp75.2pc.vec<0.05/3000),sum(minpall75.2pc.vec<0.05/3000),sum(minpall_lf75.2pc.vec<0.05/3000))

ACATp75.2pc.vec = sapply(p75.OR5.02.list,ACAT.fct)
ACATpall75.2pc.vec = sapply(pall75.OR5.02.list,ACAT.fct)
ACATpall_lf75.2pc.vec = sapply(pall_lf75.OR5.02.list,ACAT.fct)
power75.2pc.05 = c(sum(ACATp75.2pc.vec<0.05),sum(ACATpall75.2pc.vec<0.05),sum(ACATpall_lf75.2pc.vec<0.05))/length(ACATp75.2pc.vec)
power75.2pc.3000 = c(sum(ACATp75.2pc.vec<0.05/3000),sum(ACATpall75.2pc.vec<0.05/3000),sum(ACATpall_lf75.2pc.vec<0.05/3000))/length(ACATp75.2pc.vec)
power75.2pc.6000 = c(sum(ACATp75.2pc.vec<0.05/6000),sum(ACATpall75.2pc.vec<0.05/6000),sum(ACATpall_lf75.2pc.vec<0.05/3000))/length(ACATp75.2pc.vec)

# 1% of causal variants
minp75.1pc.vec = sapply(p75.OR5.01.list,minp_Bonf)
minpall75.1pc.vec = sapply(pall75.OR5.01.list,minp_Bonf)
minpall_lf75.1pc.vec = sapply(pall_lf75.OR5.01.list,minp_Bonf)
#power75.1pc.05 = c(sum(minp75.1pc.vec<0.05),sum(minpall75.1pc.vec<0.05),sum(minpall_lf75.1pc.vec<0.05))
#power75.1pc.3000 = c(sum(minp75.1pc.vec<0.05/3000),sum(minpall75.1pc.vec<0.05/3000),sum(minpall_lf75.1pc.vec<0.05/3000))

ACATp75.1pc.vec = sapply(p75.OR5.01.list,ACAT.fct)
ACATpall75.1pc.vec = sapply(pall75.OR5.01.list,ACAT.fct)
ACATpall_lf75.1pc.vec = sapply(pall_lf75.OR5.01.list,ACAT.fct)
power75.1pc.05 = c(sum(ACATp75.1pc.vec<0.05),sum(ACATpall75.1pc.vec<0.05),sum(ACATpall_lf75.1pc.vec<0.05))/length(ACATp75.1pc.vec)
power75.1pc.3000 = c(sum(ACATp75.1pc.vec<0.05/3000),sum(ACATpall75.1pc.vec<0.05/3000),sum(ACATpall_lf75.1pc.vec<0.05/6000))/length(ACATp75.1pc.vec)
power75.1pc.6000 = c(sum(ACATp75.1pc.vec<0.05/6000),sum(ACATpall75.1pc.vec<0.05/6000),sum(ACATpall_lf75.1pc.vec<0.05/6000))/length(ACATp75.1pc.vec)

# 50% causal variants in CRH 1

# 2% of causal variants
minp50.2pc.vec = sapply(p50.OR5.02.list,minp_Bonf)
minpall50.2pc.vec = sapply(pall50.OR5.02.list,minp_Bonf)
minpall_lf50.2pc.vec = sapply(pall_lf50.OR5.02.list,minp_Bonf)
#power50.2pc.05 = c(sum(minp50.2pc.vec<0.05),sum(minpall50.2pc.vec<0.05),sum(minpall_lf50.2pc.vec<0.05))
#power50.2pc.00005 = c(sum(minp50.2pc.vec<0.00005),sum(minpall50.2pc.vec<0.00005),sum(minpall_lf50.2pc.vec<0.00005))

ACATp50.2pc.vec = sapply(p50.OR5.02.list,ACAT.fct)
ACATpall50.2pc.vec = sapply(pall50.OR5.02.list,ACAT.fct)
ACATpall_lf50.2pc.vec = sapply(pall_lf50.OR5.02.list,ACAT.fct)
power50.2pc.05 = c(sum(ACATp50.2pc.vec<0.05),sum(ACATpall50.2pc.vec<0.05),sum(ACATpall_lf50.2pc.vec<0.05))/length(ACATp50.2pc.vec)
power50.2pc.3000 = c(sum(ACATp50.2pc.vec<0.05/3000),sum(ACATpall50.2pc.vec<0.05/3000),sum(ACATpall_lf50.2pc.vec<0.05/3000))/length(ACATp50.2pc.vec)
power50.2pc.6000 = c(sum(ACATp50.2pc.vec<0.05/6000),sum(ACATpall50.2pc.vec<0.05/6000),sum(ACATpall_lf50.2pc.vec<0.05/6000))/length(ACATp50.2pc.vec)

# 1% of causal variants
minp50.1pc.vec = sapply(p50.OR5.01.list,minp_Bonf)
minpall50.1pc.vec = sapply(pall50.OR5.01.list,minp_Bonf)
minpall_lf50.1pc.vec = sapply(pall_lf50.OR5.01.list,minp_Bonf)
#power50.1pc.05 = c(sum(minp50.1pc.vec<0.05),sum(minpall50.1pc.vec<0.05),sum(minpall_lf50.1pc.vec<0.05))
#power50.1pc.00005 = c(sum(minp50.1pc.vec<0.00005),sum(minpall50.1pc.vec<0.00005),sum(minpall_lf50.1pc.vec<0.00005))

ACATp50.1pc.vec = sapply(p50.OR5.01.list,ACAT.fct)
ACATpall50.1pc.vec = sapply(pall50.OR5.01.list,ACAT.fct)
ACATpall_lf50.1pc.vec = sapply(pall_lf50.OR5.01.list,ACAT.fct)
power50.1pc.05 = c(sum(ACATp50.1pc.vec<0.05),sum(ACATpall50.1pc.vec<0.05),sum(ACATpall_lf50.1pc.vec<0.05))/length(ACATp50.1pc.vec)
power50.1pc.3000 = c(sum(ACATp50.1pc.vec<0.05/3000),sum(ACATpall50.1pc.vec<0.05/3000),sum(ACATpall_lf50.1pc.vec<0.05/3000))/length(ACATp50.1pc.vec)
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

# Examine number of windows for variable length
l100.4pc = sapply(p100.04.list,length)
l75.4pc = sapply(p75.04.list,length)
l50.4pc = sapply(p50.04.list,length)
summary(l100.4pc)
summary(l75.4pc)
summary(l50.4pc)
