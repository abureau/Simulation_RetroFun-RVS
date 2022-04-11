minp_Bonf = function(vec)
{
	min(vec,na.rm=T)*length(vec)
}

minp.vec = sapply(p.list,minp_Bonf)
minpall.vec = sapply(pall.list,minp_Bonf)
minpall_lf.vec = sapply(pall_lf.list,minp_Bonf)
sum(minp.vec<0.05)
sum(minp.vec<0.00005)
sum(minpall.vec<0.05)
sum(minpall.vec<0.00005)
sum(minpall_lf.vec<0.05)
sum(minpall_lf.vec<0.00005)


