drop.haplo.locus.fn = function(g1,g2,go,h1,h2)
{
	# g1 : number of variant alleles in paternal genotype (0,1 or 2)
	# g2 : number of variant alleles in maternal genotype (0,1 or 2)
	# go : number of variant alleles in offspring genotype (0,1 or 2)
	# h1 : vector of length 2 with the indices of the two paternal haplotypes
	# h2 : vector of length 2 with the indices of the two maternal haplotypes
	if( any(is.na(h1)) | any(is.na(h2)) )
    hoff <- rep(NA,2)
    else 
    # If both parents are homozygous, sample randomly child haplotypes from parental haplotypes
    if (g1 %in% c(0,2) & g2 %in% c(0,2)) hoff = c(sample(h1,1),sample(h2,1))
    else
    # If both parents are heterozygous
    if (g1 == 1 & g2 == 1)
    {
    	# If child is homozygous reference, transmit haplotypes carrying reference (always second)
    	if (go == 0) hoff = c(h1[2],h2[2])
    	else
    	# If child is heterozygous, randomly sample which parent transmitted variant
    	if (go == 1)
    	{
    		if (runif(1)<0.5) hoff = c(h1[1],h2[2])
    		else hoff = c(h1[2],h2[1])
       	}
       	else # child is homozygous variant, transmit haplotypes carrying variant (always first)
       		hoff = c(h1[1],h2[1])
    }
    else
    # If father is heterozygous
    if (g1 == 1 & g2 == 0)
    {
    	# If child is homozygous reference, transmit paternal haplotype carrying reference (always second)
       	if (go == 0) hoff = c(h1[2],sample(h2,1))
    	# If child is heterozygous, transmit paternal haplotype carrying variant (always first)
       	else if (go == 1) hoff = c(h1[1],sample(h2,1))
       	# If child is homozygous variant, Mendelian inconsistency
       	else hoff = "error"
    }
    else if (g1 == 1 & g2 == 2)
    {
    	# If child is homozygous reference, Mendelian inconsistency 
       	if (go == 0) hoff = "error"
    	# If child is heterozygous, transmit paternal haplotype carrying reference (always second) 
       	else if (go == 1) hoff = c(sample(h2,1),h1[2])
       	# If child is homozygous variant, transmit paternal haplotype carrying variant (always first)
       	else hoff = c(h1[1],sample(h2,1))
    }
    else
        # If mother is heterozygous
    if (g1 == 0 & g2 == 1)
    {
    	# If child is homozygous reference, transmit maternal haplotype carrying reference (always second)
       	if (go == 0) hoff = c(sample(h1,1),h2[2])
    	# If child is heterozygous, transmit maternal haplotype carrying variant (always first)
       	else if (go == 1) hoff = c(h2[1],sample(h1,1))
       	# If child is homozygous variant, Mendelian inconsistency
       	else hoff = "error"
    }
    else if (g1 == 2 & g2 == 1)
    {
    	# If child is homozygous reference, Mendelian inconsistency 
       	if (go == 0) hoff = "error"
    	# If child is heterozygous, transmit maternal haplotype carrying reference (always second) 
       	else if (go == 1) hoff = c(sample(h1,1),h2[2])
       	# If child is homozygous variant, transmit maternal haplotype carrying variant (always first)
       	else hoff = c(sample(h1,1),h2[1])
    }
    hoff 
}

drop.haplo.fn = function(h1,h2)
{
	if( any(is.na(h1)) | any(is.na(h2)) )
    hoff <- rep(NA,2)
    else 
    hoff = c(sample(h1,1),sample(h2,1))
    hoff 
}
    
