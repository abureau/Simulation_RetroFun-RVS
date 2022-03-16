#!/usr/bin/env Rscript
cat("The script is launching...\n")
library(GenomicRanges, quietly=T)
library(rtracklayer, quietly=T)

#The random seed is generated randomly to add sample noise at each generation
set.seed(sample(1:1000,1))

args = commandArgs(trailingOnly=TRUE)

#args[1] => sfs file for variants within TAD
#args[2] => CRHs elements 
#args[3] => causal probabibility 
#args[4] => Null or Alter

prob_causal = as.numeric(args[3])

cat("Import variant file \n")
#Variants within the region are considered
#TAD: chr1-Chr1:24100000-24970000
variants = read.table(args[1], header=TRUE, sep="\t")
colnames(variants) = c("chr", "pos", "ref","alt", "maf")

variants = variants[variants$maf<=0.01,]
MAF = do.call("rbind",strsplit(variants$maf,",", fixed=T))

MAF = apply(MAF,2, as.numeric)
variants$maf = apply(MAF,1,min)

n_variants = nrow(variants)


cat("Import CRH file\n")
#GRange elements in bed format with score column corresponding to CRH belonging
GRange_CRHs = import(args[2], format="bed",extraCols=c(Name="character", Annotation="numeric"))
GRange_variants = GRanges(seqnames=paste0("chr",variants$chr), ranges=IRanges(start=variants$pos,end=variants$pos))

variants_tmp = variants

overlaps_variants_CRHs = countOverlaps(GRange_variants, GRange_CRHs)
overlaps_variants_CRHs[overlaps_variants_CRHs>1] = 1 

variants_in_CRHs = findOverlaps(GRange_variants, GRange_CRHs)
index_variants_in_CRHs = unique(queryHits(variants_in_CRHs))


variants_tmp[index_variants_in_CRHs,"inCRHs"] = TRUE
index_variants_out_CRHs = which(is.na(variants_tmp$inCRHs))

#Causal variants are generated based on the probability
#In the entire regions X% of variants will be causal
causal_variants = sample(0:1, nrow(variants_tmp), replace=T, prob=c((1-prob_causal),prob_causal))
n_causal_variants = length(causal_variants[causal_variants==1])

cat("Scenario 1 \n")
#Scenario 1 : 100%
#All the causal variants are located within CRH elements 
variants_tmp[sample(index_variants_in_CRHs, n_causal_variants, replace=F), "Scenario_100"] = 1

cat("Scenario 2 \n")
#Scenario 2 : 75%
#75% of causal variants are located within CRHs the rest are outside CRH elements
n_causal_in_CRHs = floor(n_causal_variants*0.75)
n_causal_out_CRHs = floor(n_causal_variants*0.25)

if((length(n_causal_in_CRHs) + length(n_causal_out_CRHs))!=n_causal_variants) n_causal_in_CRHs = n_causal_in_CRHs +1

index_causal_in_CRHs = sample(index_variants_in_CRHs, n_causal_in_CRHs , replace=F)
index_causal_out_CRHs = sample(index_variants_out_CRHs, n_causal_out_CRHs , replace=F)

variants_tmp[c(index_causal_in_CRHs, index_causal_out_CRHs), "Scenario_75"] = 1

cat("Scenario 3 \n")
#Scenario 3 : 50%
#Equal part of causal variants within and outside CRH elements
n_causal_in_CRHs = floor(n_causal_variants*0.50)
n_causal_out_CRHs = floor(n_causal_variants*0.50)

if((length(n_causal_in_CRHs) + length(n_causal_out_CRHs))!=n_causal_variants) n_causal_in_CRHs = n_causal_in_CRHs + 1

index_causal_in_CRHs = sample(index_variants_in_CRHs, n_causal_in_CRHs , replace=F)
index_causal_out_CRHs = sample(index_variants_out_CRHs, n_causal_out_CRHs , replace=F)

variants_tmp[c(index_causal_in_CRHs, index_causal_out_CRHs), "Scenario_50"] = 1


variants_tmp[is.na(variants_tmp$Scenario_100), "Scenario_100"] = 0
variants_tmp[is.na(variants_tmp$Scenario_75), "Scenario_75"] = 0
variants_tmp[is.na(variants_tmp$Scenario_50), "Scenario_50"] = 0

split_CRHs = split(GRange_CRHs, GRange_CRHs$Annotation)

variants_per_CRH = lapply(1:length(split_CRHs), function(x){
        c = findOverlaps(GRange_variants, split_CRHs[[x]])
        q = queryHits(c)

        v = variants[unique(q),]
        v[,"Annot"] = paste0("CRH",x)
        v[,"Score"]=1
        v
})

df_variants_per_CRH = do.call("rbind", variants_per_CRH)
df_variants = merge(variants_tmp, df_variants_per_CRH, by=c("chr", "pos", "alt","ref","maf"), all.x=T)

df_variants$Annot[is.na(df_variants$Annot)] = "TAD"

variants_scenario100 = df_variants[,c("Annot","chr","pos","ref","alt","maf","Scenario_100")]
variants_scenario100$Scenario_100[variants_scenario100$Annot=="CRH4"] = 1 
 
variants_scenario75 = df_variants[,c("Annot","chr","pos","ref","alt","maf","Scenario_75")]
variants_scenario75$Scenario_75[variants_scenario75$Annot=="CRH4"] = 1

variants_scenario50 = df_variants[,c("Annot","chr","pos","ref","alt","maf","Scenario_50")]
variants_scenario50$Scenario_50[variants_scenario50$Annot=="CRH4"] = 1

write.table(variants_scenario100, paste0("rare_variants_scenario100_",prob_causal,".sfs"), quote=F, col.names=F, row.names=F)
write.table(variants_scenario75, paste0("rare_variants_scenario75_",prob_causal,".sfs"), quote=F, col.names=F, row.names=F)
write.table(variants_scenario50, paste0("rare_variants_scenario50_",prob_causal,".sfs"), quote=F, col.names=F, row.names=F)


cat("Annotation matrix generation \n")

GRanges_variants_scenario100 = GRanges(paste0("chr",seqnames=variants_scenario100[,2]), ranges=IRanges(start=variants_scenario100[,3], end=variants_scenario100[,3]))
GRanges_variants_scenario75 = GRanges(paste0("chr",seqnames=variants_scenario75[,2]), ranges=IRanges(start=variants_scenario75[,3], end=variants_scenario75[,3]))
GRanges_variants_scenario50 = GRanges(paste0("chr",seqnames=variants_scenario50[,2]), ranges=IRanges(start=variants_scenario50[,3], end=variants_scenario50[,3]))


Gl_variants_scenarios = GRangesList(GRanges_variants_scenario100, GRanges_variants_scenario75, GRanges_variants_scenario50)

#Annotation matrix 
#binary matrix are created for each scenario

split_annotation = split(GRange_CRHs, GRange_CRHs$Annotation)

annotation_matrix_scenarios = lapply(Gl_variants_scenarios, function(x) {
	sapply(split_annotation, function(y){
		countOverlaps(x,y)
	})
})

mat1 = annotation_matrix_scenarios[[1]]
mat1[mat1>1] = 1
mat2 = annotation_matrix_scenarios[[2]]
mat2[mat2>1] = 1
mat3 = annotation_matrix_scenarios[[3]]
mat3[mat3>1] = 1

write.table(mat1, paste0("annotation_matrix_scenario100_", prob_causal,".mat"), quote=F, col.names=F, row.names=F)
write.table(mat2, paste0("annotation_matrix_scenario75_", prob_causal,".mat"), quote=F, col.names=F, row.names=F)
write.table(mat3, paste0("annotation_matrix_scenario50_", prob_causal,".mat"), quote=F, col.names=F, row.names=F)
