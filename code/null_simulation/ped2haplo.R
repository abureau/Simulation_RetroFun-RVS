ped2haplo = function(ped.mat)
{
	# Converts ordered genotypes in a ped file into haplotypes
	# ped.mat : pedigrees coded as in a ped file
	haplo1 = as.matrix(ped.mat[,seq(7,ncol(ped.mat)-1,by=2)])
	haplo2 = as.matrix(ped.mat[,seq(8,ncol(ped.mat),by=2)])	
	rbind(haplo1,haplo2)
	}
	
generate.haplos.locus = function(ped.mat,trio.list,haplo.mat,fams,site,locus=1)
{
	# Samples founder haplotypes and drops them down pedigrees according to the transmission simulated
	# at the selected locus in ped.mat
	# ped.mat : pedigrees coded as in a ped file
	# trio.list : list of trio objects, one object for each pedigree in ped.mat
	# haplo.mat : matrix where each row defines a haplotype, from which to sample the founder haplotypes
	# fams : vector of family names
	# sites : vector of the variant sites for each family. The sites must be named by the family in which they are found 
	# 		in the fams vector 
	# locus : which locus to use as variant site (for conversion from ped.mat)
	
	# Extract variant site genotype from two locus genotype
	if (!(locus %in% 1:2)) stop ("Locus must be 1 or 2")
	if (locus == 2)	variant.geno = (ped.mat[,7] - 1)%%3
	else variant.geno = (ped.mat[,7] - 1)%/%3
	names(variant.geno) = ped.mat[,2]

    nhap = nrow(haplo.mat)
	haplo.obs = matrix(NA,nrow(ped.mat),2)
	dimnames(haplo.obs)[[1]] = ped.mat[,2]
	# Loop over the families
	for (f in fams)
	{	
		#print(f)
		# Identify minor allele
		tab.alleles = table(haplo.mat[,site[f]])
		minor.allele = names(tab.alleles)[which.min(tab.alleles)]
		# Identify haplotype(s) carrying variant
		variant.haplo = which(haplo.mat[,site[f]]==minor.allele)
		# print(variant.haplo)
		# Sample founder haplotypes
		founders = ped.mat[,1]==f & ped.mat[,3]==0
		nf = sum(founders)
		if (length(variant.haplo)>1)
		{
			haplo.obs[founders,1][variant.geno[founders]>0] = sample(variant.haplo,sum(variant.geno[founders]>0),replace=TRUE)
		# Second haplotype can only be variant if founder homozygous for the variant
			haplo.obs[founders,2][variant.geno[founders]==2] = sample(variant.haplo,sum(variant.geno[founders]==2),replace=TRUE)			
		}
		else 
		{
			haplo.obs[founders,1][variant.geno[founders]>0] = variant.haplo		
		# Second haplotype can only be variant if founder homozygous for the variant
			haplo.obs[founders,2][variant.geno[founders]==2] = variant.haplo		
		}
		haplo.obs[founders,1][variant.geno[founders]==0] = sample((1:nhap)[-variant.haplo],sum(variant.geno[founders]==0),replace=TRUE)
		haplo.obs[founders,2][variant.geno[founders]<2] = sample((1:nhap)[-variant.haplo],sum(variant.geno[founders]<2),replace=TRUE)
		
		for (j in 1:length(trio.list[f]))
		# This line can only work if the families are named 1,2,..., the names returned by generate.pedigrees.geno
			haplo.obs[ped.mat[,1]==f,] = DropHaploLocus(trio.list[f][[j]],haplo.obs[ped.mat[,1]==f,],variant.geno[ped.mat[,1]==f])
		#if (any(is.na(haplo.obs))) stop("DropHaploLocus returned NA genotype.")			
	}
	haplo.obs
}

HaploDropSim.fn = function(ped.mat,trio.list,fams,nhap)
{
	# Samples founder haplotypes and drops them down pedigrees according to Mendel's laws
	# ped.mat : pedigrees coded as in a ped file
	# trio.list : list of trio objects, one object for each pedigree in ped.mat
	# fams : vector of family names
	# nhap : number of haplotypes in population
	
	haplo.obs = matrix(NA,nrow(ped.mat),2)
	dimnames(haplo.obs)[[1]] = ped.mat[,2]
	# Loop over the families
	for (f in fams)
	{	
		#print(f)
		# Sample founder haplotypes
		founders = ped.mat[,1]==f & ped.mat[,3]==0
		nf = sum(founders)
		haplo.obs[founders,1] = sample((1:nhap),nf,replace=TRUE)
		haplo.obs[founders,2] = sample((1:nhap),nf,replace=TRUE)
		
		for (j in 1:length(trio.list[f]))
			haplo.obs[ped.mat[,1]==f,] = DropHaplo(trio.list[f][[j]],haplo.obs[ped.mat[,1]==f,])
	}
	haplo.obs	
}