library(GenomicRanges, quietly=T)
library(rtracklayer, quietly=T)

        
simulate_variants = function(sfs_file, GRange_regions, prob_causal=NULL, proportion_within_region=1, region_prefix = "Region_", seed=1234){
        set.seed(seed)
        if(is.null(prob_causal)) print("Please provide a proportion of causal variants")
        
        else{
                variants = read.table(sfs_file, header=F, sep="\t")
                colnames(variants) = c("chr", "pos", "ref","alt", "maf")
                
                MAF = do.call("rbind",strsplit(variants$maf,",", fixed=T))
                MAF = apply(MAF,2, as.numeric)
                
                variants$maf = apply(MAF,1,min)
                variants_rare = variants[variants$maf<=0.01,]
                
                GRange_regions = import(GRange_regions, format="bed",extraCols=c(Name="character", Annotation="numeric"))
                GRange_variants = GRanges(seqnames=paste0("chr",variants_rare$chr), ranges=IRanges(start=variants_rare$pos,end=variants_rare$pos))
                
                variants_tmp = variants_rare
                
                overlaps_variants_regions = countOverlaps(GRange_variants, GRange_regions)
                overlaps_variants_regions[overlaps_variants_regions>1] = 1 
                
                index_in_regions = which(overlaps_variants_regions==1)
                index_out_regions = which(overlaps_variants_regions==0)
                
                variants_tmp[index_in_regions,"inRegions"] = TRUE
                variants_tmp[index_out_regions,"inRegions"] = FALSE
                
                
                #Causal variants are generated based on the probability
                #In the entire regions X% of variants will be causal
                
                causal_variants = sample(0:1, nrow(variants_tmp), replace=T, prob=c((1-prob_causal),prob_causal))
                n_causal_variants = length(causal_variants[causal_variants==1])
                
                n_causal_in_regions = floor(n_causal_variants*proportion_within_region)
                n_causal_out_regions = floor(n_causal_variants*(1-proportion_within_region))
                
                if((length(n_causal_in_regions) + length(n_causal_out_regions))!=n_causal_variants) n_causal_in_regions = n_causal_in_regions +1
                
                index_causal_in_regions = sample(index_in_regions, n_causal_in_regions , replace=F)
                index_causal_out_regions = sample(index_out_regions, n_causal_out_regions , replace=F)
                
                variants_tmp[c(index_causal_in_regions, index_causal_out_regions), paste0("Scenario_",proportion_within_region)] = 1
                variants_tmp[,paste0("Scenario_",proportion_within_region)][is.na(variants_tmp[,paste0("Scenario_",proportion_within_region)])] = 0
                
                split_regions_by_annot = split(GRange_regions, GRange_regions$Annotation)
                
                variants_per_region = lapply(1:length(split_regions_by_annot), function(x){
                        c = findOverlaps(GRange_variants, split_regions_by_annot[[x]])
                        q = queryHits(c)
                        
                        v = variants_rare[unique(q),]
                        v[,"Annot"] = paste0(region_prefix,x)
                        v[,"Score"]=1
                        v
                })
                
                df_variants_per_region = do.call("rbind", variants_per_region)
                df_variants = merge(variants_tmp, df_variants_per_region, by=c("chr", "pos", "alt","ref","maf"), all.x=T)
                df_variants$Annot[is.na(df_variants$Annot)] = "TAD"
                df_variants$Score = NULL
                df_variants = df_variants[order(df_variants$Annot),]
                
                GRange_variants_with_annot = GRanges(seqnames = paste0("chr",df_variants$chr), ranges=IRanges(start=df_variants$pos,end=df_variants$pos))
                matrix_annot = sapply(split_regions_by_annot, function(x){
                        countOverlaps(GRange_variants_with_annot,x)
                })
                matrix_annot[matrix_annot>1] = 1
        }
        
        Cat_where_only_zeros = names(which(table(df_variants$Annot,df_variants[,paste0("Scenario_",proportion_within_region)])[,2]==0))
        df_variants = df_variants[!df_variants$Annot%in%c(Cat_where_only_zeros),]
        df_variants = df_variants[,c(8,1,2,3,4,5,7)]
        return(list('sfs'=df_variants, "Z" = matrix_annot))
}

sfs_95 = simulate_variants("ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.chr1_24100000_24970000.classic.annot.NONS_MISS_SPLICE_withHeader.sfs", "CRHs_in_TAD_ch1_24100000_24970000_with_annotations.bed",
                           prob_causal = 0.1, proportion_within_region = 0.95, region_prefix = "CRH", seed = 12350)

sfs_90 = simulate_variants("ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.chr1_24100000_24970000.classic.annot.NONS_MISS_SPLICE_withHeader.sfs", "CRHs_in_TAD_ch1_24100000_24970000_with_annotations.bed",
                  prob_causal = 0.1, proportion_within_region = 0.90, region_prefix = "CRH", seed = 12350)

sfs_80 = simulate_variants("ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.chr1_24100000_24970000.classic.annot.NONS_MISS_SPLICE_withHeader.sfs", "CRHs_in_TAD_ch1_24100000_24970000_with_annotations.bed",
                           prob_causal = 0.1, proportion_within_region = 0.80, region_prefix = "CRH", seed = 12350)

sfs_75 = simulate_variants("ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.chr1_24100000_24970000.classic.annot.NONS_MISS_SPLICE_withHeader.sfs", "CRHs_in_TAD_ch1_24100000_24970000_with_annotations.bed",
                           prob_causal = 0.1, proportion_within_region = 0.75, region_prefix = "CRH", seed = 12350)

table(sfs_100$sfs$Annot, sfs_100$sfs$Scenario_1)
table(sfs_90$sfs$Annot,sfs_90$sfs$Scenario_0.9)
table(sfs_80$sfs$Annot,sfs_80$sfs$Scenario_0.8)
table(sfs_70$sfs$Annot,sfs_70$sfs$Scenario_0.7)

write.table(sfs_90$sfs, "rare_variants_scenario90_0.10.sfs", col.names = F, row.names = F, quote = F)
write.table(sfs_75$sfs, "rare_variants_scenario75_0.10.sfs", col.names = F, row.names = F, quote = F)

colSums(sfs_90$Z)
table(sfs_90$sfs$Annot, sfs_90$sfs$Scenario_0.9)
