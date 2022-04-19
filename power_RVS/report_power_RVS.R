minp_Bonf = function(vec)
{
	min(vec,na.rm=T)*length(vec)
}
load("alter100.RData")
load("alter75.RData")
load("alter50.RData")

# Close to 100% causal variants in CRH 1
minp100.vec = sapply(p100.list,minp_Bonf)
minpall100.vec = sapply(pall100.list,minp_Bonf)
minpall_lf100.vec = sapply(pall_lf100.list,minp_Bonf)
power100.05 = c(sum(minp100.vec<0.05),sum(minpall100.vec<0.05),sum(minpall_lf100.vec<0.05))
power100.00005 = c(sum(minp100.vec<0.00005),sum(minpall100.vec<0.00005),sum(minpall_lf100.vec<0.00005))

# 75% causal variants in CRH 1
minp75.vec = sapply(p75.list,minp_Bonf)
minpall75.vec = sapply(pall75.list,minp_Bonf)
minpall_lf75.vec = sapply(pall_lf75.list,minp_Bonf)
power75.05 = c(sum(minp75.vec<0.05),sum(minpall75.vec<0.05),sum(minpall_lf75.vec<0.05))
power75.00005 = c(sum(minp75.vec<0.00005),sum(minpall75.vec<0.00005),sum(minpall_lf75.vec<0.00005))

# 50% causal variants in CRH 1
minp50.vec = sapply(p50.list,minp_Bonf)
minpall50.vec = sapply(pall50.list,minp_Bonf)
minpall_lf50.vec = sapply(pall_lf50.list,minp_Bonf)
power50.05 = c(sum(minp50.vec<0.05),sum(minpall50.vec<0.05),sum(minpall_lf50.vec<0.05))
power50.00005 = c(sum(minp50.vec<0.00005),sum(minpall50.vec<0.00005),sum(minpall_lf50.vec<0.00005))

power.05 = cbind(power100.05,power75.05,power50.05)
power.00005 = cbind(power100.00005,power75.00005,power50.00005)
row.names(power.00005) = row.names(power.05) = c("Partial sharing","Complete sharing","Complete sharing, fixed length")
power.05
power.00005

# Examine number of windows for variable length
l100 = sapply(p100.list,length)
l75 = sapply(p75.list,length)
l50 = sapply(p50.list,length)
summary(l100)
summary(l75)
summary(l50)
