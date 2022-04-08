library(GenomicRanges, quietly=T)
library(rtracklayer, quietly=T)

        
simulate_variants = function(sfs_file, GRange_regions, prob_causal=NULL, proportion_within_region=1, seed=1234,
                             target_regions=NULL){
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
                
                n_causal = floor(length(GRange_variants)*prob_causal)
                
                n_causal_in_target_region = floor(n_causal*proportion_within_region)
                n_causal_outside_target_region = n_causal-n_causal_in_target_region

                split_regions = split(GRange_regions, GRange_regions$Annotation)
                
                annot_variants = lapply(split_regions, function(x){
                        f = findOverlaps(GRange_variants, x)
                        q = queryHits(f)
                        v = variants_rare[unique(q),]
                        v$Annotation = unique(x$Annotation)
                        v
                })
                
                df_annot_variants = do.call("rbind", annot_variants)
                df_variants_all = merge(variants_rare,df_annot_variants, by=c("chr", "pos", "ref", "alt", "maf"), all.x=T)
                
                df_variants_all$Annotation[is.na(df_variants_all$Annotation)] = "Out"
               
                
                if(!is.null(target_regions)){
                        index_in_target_regions = which(df_variants_all$Annotation%in%target_regions)
                        index_out_target_regions = which(!df_variants_all$Annotation%in%target_regions)
                }
                
                else{
                    target_regions = unique(df_variants_all$Annotation[df_variants_all$Annotation!="Out"])
                    index_in_target_regions = which(df_variants_all$Annotation%in%target_regions)
                    index_out_target_regions = which(!df_variants_all$Annotation%in%target_regions)
                }
                
                random_index_causal_in_target_regions = sample(index_in_target_regions,n_causal_in_target_region, replace = F)
                random_index_causal_out_target_regions = sample(index_out_target_regions,n_causal_outside_target_region, replace = F)
                
                df_variants_all[c(random_index_causal_in_target_regions,random_index_causal_out_target_regions),"Score"] = 1
                df_variants_all$Score[is.na(df_variants_all$Score)] = 0
                
                
        }
        
        agg_Score_by_Annot = aggregate(Score~Annotation,df_variants_all,sum)
        annot_with_zero_score = agg_Score_by_Annot$Annotation[which(agg_Score_by_Annot$Score==0)]
       
        if(length(annot_with_zero_score)==0){
                df_variants_all = df_variants_all
        } else{
                for(a in annot_with_zero_score){
                        index_annot = which(df_variants_all$Annotation == a)
                        df_variants_all[sample(index_annot,1),"Score"] = 1
                }
        }
        
        df_variants_all = df_variants_all[order(df_variants_all$Annotation),c(6,1,2,3,4,5,7)]
        
        u = unique(df_variants_all$Annotation)
        Z = matrix(0, nrow=nrow(df_variants_all),ncol=length(u)-1)
        Z = sapply(1:ncol(Z),function(x){
                w = which(df_variants_all$Annotation==u[x])
                Z[w,x] = 1
                Z
        })
        Z = cbind(1,Z)
        return(list("sfs"=df_variants_all,"Z"=Z))
}

sfs_100_02 = simulate_variants("ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.chr1_24100000_24970000.classic.annot.NONS_MISS_SPLICE_withHeader.sfs", "CRHs_in_TAD_ch1_24100000_24970000_with_annotations.bed",
                           prob_causal = 0.02, proportion_within_region = 1, seed = 8765, target_regions = c("1045"))

table(sfs_100_02$sfs$Annotation,sfs_100_02$sfs$Score)
write.table(sfs_100_02$sfs, "rare_variants_scenario100_0.02.sfs", col.names = F, row.names = F,quote = F)
