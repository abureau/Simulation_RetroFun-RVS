#This code generates figures presented in the paper
#The data used in this script are available in the repo

library(ggplot2)
library(ggpubr)

gg_qqplot_facet_grid = function(list_pvalues,ci = 0.95, title="") {
  #n  <- length(ps)
  
  dfs <- lapply(1:length(list_pvalues), function(i) {
    n = length(list_pvalues[[i]])
    data.frame(
      Distance_Kernel = names(list_pvalues)[i],
      observed = -log10(sort(list_pvalues[[i]])),
      expected = -log10(ppoints(n)),
      clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
      cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
    )})
  
  all.dfs = do.call("rbind", dfs)
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  
  ggplot(all.dfs) +
    geom_ribbon(
      mapping = aes(x = expected, ymin = clower, ymax = cupper),
      alpha = 0.1
    ) +
    geom_point(aes(expected, observed, color=Distance_Kernel), size = 2) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5)+xlab(log10Pe) +
    ylab(log10Po)+theme_bw()#+theme(legend.position = "none")
  
}

#Main results
pvalues_null_SW = readRDS("data\\pvalues_null_SW.RDS")
pvalues_null_Pairs = readRDS("data\\pvalues_null_Pairs.RDS")
pvalues_null_Genes = readRDS("data\\pvalues_null_Genes.RDS")
pvalues_null = readRDS("data\\pvalues_null_CRHs.RDS")
pvalues_null_indep = readRDS("data\\pvalues_null_CRHs_indep.RDS")

#Figure 3 + Figure S4-S5-S8: Control of Type-I error rate
gg_qqplot_facet_grid(list("Dependence"=sapply(pvalues_null, function(x) x$ACAT)))+theme(legend.position = "none")

ggarrange(gg_qqplot_facet_grid(list("Dependence"=sapply(pvalues_null, function(x) x$ACAT), 
                                              "Independence" = sapply(pvalues_null_indep, function(x) x$ACAT)))+labs(colour="Variant Structure"),
          gg_qqplot_facet_grid(list("Dependence"=sapply(pvalues_null, function(x) x$Fisher), 
                                    "Independence" = sapply(pvalues_null_indep, function(x) x$Fisher)))+labs(colour="Variant Structure"),
          gg_qqplot_facet_grid(list("Dependence"=sapply(pvalues_null, function(x) x$Score_V1),
                                              "Independence"=sapply(pvalues_null_indep, function(x) x$Score_V1)))+labs(colour="Variant Structure"),
          gg_qqplot_facet_grid(list("Dependence"=sapply(pvalues_null, function(x) x$Score_V2),
          "Independence"=unlist(sapply(pvalues_null_indep, function(x) x$Score_V2))))+labs(colour="Variant Structure"),
          gg_qqplot_facet_grid(list("Dependence"=unlist(sapply(pvalues_null, function(x) x$Score_V3)),
          "Independence"=unlist(sapply(pvalues_null_indep, function(x) x$Score_V3))))+labs(colour="Variant Structure"),
          gg_qqplot_facet_grid(list("Dependence"=unlist(sapply(pvalues_null, function(x) x$Score_V5)),
          "Independence"=unlist(sapply(pvalues_null_indep, function(x) x$Score_V5))))+labs(colour="Variant Structure"), labels=c("A","B","C","D","E", "F"),ncol=2, nrow=3)


gg_qqplot_facet_grid(list("Genes"=sapply(pvalues_null_Genes, function(x) x$ACAT),
                          "Sliding-Windows"=c(sapply(pvalues_null_SW[-which(sapply(pvalues_null_SW, length)==1)], function(x) x$ACAT), unlist(pvalues_null_SW[which(sapply(pvalues_null_SW, length)==1)])),
                          "Pairs"=sapply(pvalues_null_Pairs, function(x) x$ACAT)))+labs(colour="Functional Annotation")

#Small pedigrees
pvalues_null_smallped_indep = readRDS("data\\pvalues_null_smallped_independence.RDS")
pvalues_null_smallped_dep = readRDS("data\\pvalues_null_smallped_dependence.RDS")

gg_qqplot_facet_grid(list("Dependence"=sapply(pvalues_null_smallped_dep, function(x) x$ACAT), "Independence"=sapply(pvalues_null_smallped_indep, function(x) x$ACAT)))+labs(color="Variant Structure")

#Different pedigree structures
p.values.two.affected = readRDS("data\\pvalues_null_two_affected.RDS")
p.values.three.affected= readRDS("data\\pvalues_null_three_affected.RDS")
p.values.four.affected = readRDS("data\\pvalues_null_four_affected.RDS")
p.values.five.affected = readRDS("data\\pvalues_null_five_affected.RDS")
p.values.six.affected = readRDS("data\\pvalues_null_six_affected.RDS")
p.values.seven.affected = readRDS("data\\pvalues_null_seven_affected.RDS")

#Figure S6
gg_qqplot_facet_grid(list("2 affected"=sapply(p.values.two.affected, function(x) x$ACAT),
                          "3 affected"=sapply(p.values.three.affected, function(x) x$ACAT),
                          "4 affected"=sapply(p.values.four.affected, function(x) x$ACAT),
                          "5 affected"=sapply(p.values.five.affected, function(x) x$ACAT),
                          "6 affected"=sapply(p.values.six.affected, function(x) x$ACAT),
                          "7 affected"=sapply(p.values.seven.affected, function(x) x$ACAT)))+labs(colour="Pedigree")


#Figure S7: Consanguinity 
pvalues.null.with.only.consanguinity.uncorrected = readRDS("data\\pvalues_null_only_consanguinity_uncorrected.RDS")
pvalues.null.with.only.consanguinity.corrected = readRDS("data\\pvalues_null_only_consanguinity_corrected.RDS")
pvalues.null.with.consanguinity.corrected = readRDS("data\\pvalues_null_consanguinity_corrected.RDS")
pvalues.null.with.consanguinity.unccorrected = readRDS("data\\pvalues_null_only_consanguinity_uncorrected.RDS")

ggarrange(gg_qqplot_facet_grid(list("Corrected" = sapply(pvalues.null.with.consanguinity.corrected, function(x) x$ACAT),
                          "Uncorrected"=sapply(pvalues.null.with.consanguinity.unccorrected, function(x) x$ACAT)))+labs(colour="Correction"),

gg_qqplot_facet_grid(list("Corrected" = sapply(pvalues.null.with.only.consanguinity.corrected, function(x) x$ACAT),
                          "Uncorrected"=sapply(pvalues.null.with.only.consanguinity.uncorrected, function(x) x$ACAT)))+labs(colour="Correction"), labels = c("A","B"), ncol=2)


#Figure 4: Power at 2% causal variants
pvalues.alter.2causal.CRHs = readRDS("data\\pvalues_alter_2causal_CRHs.RDS")
pvalues.alter.2causal.Genes = readRDS("data\\pvalues_alter_2causal_Genes.RDS")
pvalues.alter.2causal.Pairs = readRDS("data\\pvalues_alter_2causal_Pairs.RDS")
pvalues.alter.2causal.SW = readRDS("data\\pvalues_alter_2causal_SW.RDS")

pvalues_CHP_75causal = list.files("data\\results_RVNPL_2causal_75\\CHP", full.names = T)
pvalues_RV_75causal = list.files("data\\results_RVNPL_2causal_75\\RV", full.names = T)

pvalues_ACAT_CHP_75causal_pairs = c()
pvalues_ACAT_CHP_75causal_all = c()

for(i in 1:200){
  print(i)
  df_p = read.table(pvalues_CHP_75causal[i], fill=T, header=F)
  pvalues_ACAT_CHP_75causal_pairs = c(pvalues_ACAT_CHP_75causal_pairs, ACAT::ACAT(df_p[,1]))
  pvalues_ACAT_CHP_75causal_all = c(pvalues_ACAT_CHP_75causal_all, ACAT::ACAT(df_p[,2]))
  
}

sum(pvalues_ACAT_CHP_75causal_pairs<=8.333333e-06)/200
sum(pvalues_ACAT_CHP_75causal_all<=8.333333e-06)/200

pvalues_ACAT_RV_75causal_pairs = c()
pvalues_ACAT_RV_75causal_all = c()

for(i in 1:200){
  print(i)
  df_p = read.table(pvalues_RV_75causal[i], fill=T, header=F)
  pvalues_ACAT_RV_75causal_pairs = c(pvalues_ACAT_RV_75causal_pairs, ACAT::ACAT(df_p[,1]))
  pvalues_ACAT_RV_75causal_all = c(pvalues_ACAT_RV_75causal_all, ACAT::ACAT(df_p[,2]))
  
}

sum(pvalues_ACAT_RV_75causal_pairs<=8.333333e-06)/200
sum(pvalues_ACAT_RV_75causal_all<=8.333333e-06)/200

df_power_Retro_RVS_RVNPL = data.frame("Power" = c(sum(sapply(pvalues.alter.2causal.CRHs$`75causal`, function(x) x$ACAT)<=8.333333e-06)/1000,
                                                  sum(sapply(pvalues.alter.2causal.Genes$`75causal`, function(x) x$ACAT)<=8.333333e-06)/1000,
                                                  power.6000[2,2],power.6000[1,2],
                                                  sum(pvalues_ACAT_RV_75causal_pairs<=8.333333e-06)/200,
                                                  sum(pvalues_ACAT_RV_75causal_all<=8.333333e-06)/200,
                                                  sum(pvalues_ACAT_CHP_75causal_pairs<=8.333333e-06)/200,
                                                  sum(pvalues_ACAT_CHP_75causal_all<=8.333333e-06)/200), "Method" = c("RetroFun-RVS", "RetroFun-RVS",
                                                                                                                      "RVS", "RVS",
                                                                                                                      "RV-NPL", "RV-NPL", "CHP-NPL", "CHP-NPL"), 
                                      "Type" = c("CRHs","Genes", "Complete\nSharing", "Partial\nSharing", "Pairs", "All", "Pairs", "All"))

df_power_Retro_RVS_RVNPL$Method = factor(df_power_Retro_RVS_RVNPL$Method, levels=c("RetroFun-RVS", "RV-NPL", "CHP-NPL","RVS"))

df_power_CRHs_2causal_OR5 = data.frame("Power" = c(sum(sapply(pvalues.alter.2causal.CRHs$`100causal`, function(x) x$Score_V1)<=8.333333e-06)/1000,
                                                   sum(sapply(pvalues.alter.2causal.CRHs$`100causal`, function(x) x$ACAT)<=8.333333e-06)/1000,
                                                   
                                                   sum(sapply(pvalues.alter.2causal.CRHs$`75causal`, function(x) x$Score_V1)<=8.333333e-06)/1000,
                                                   sum(sapply(pvalues.alter.2causal.CRHs$`75causal`, function(x) x$ACAT)<=8.333333e-06)/1000,
                                                   
                                                   sum(sapply(pvalues.alter.2causal.CRHs$`50causal`, function(x) x$Score_V1)<=8.333333e-06)/1000,
                                                   sum(sapply(pvalues.alter.2causal.CRHs$`50causal`, function(x) x$ACAT)<=8.333333e-06)/1000), "Prop"= c(rep(100,2), rep(75,2), rep(50,2)), "Type" =rep(c("Burden\nOriginal", "ACAT-\nCombined"),3))

df_power_CRHs_2causal_OR5$Type = factor(df_power_CRHs_2causal_OR5$Type, levels=c("Burden\nOriginal", "ACAT-\nCombined"))

df_power_CRHs_Combined_2causal_OR5_SW = data.frame("Power"=c(sum(sapply(pvalues.alter.2causal.SW$`100causal`, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                             sum(sapply(pvalues.alter.2causal.SW$`100causal`, function(x) x$ACAT)<=8.33e-6)/1000,
                                                             
                                                             sum(sapply(pvalues.alter.2causal.SW$`75causal`, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                             sum(sapply(pvalues.alter.2causal.SW$`75causal`, function(x) x$ACAT)<=8.33e-6)/1000,
                                                             
                                                             sum(sapply(pvalues.alter.2causal.SW$`50causal`, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                             sum(sapply(pvalues.alter.2causal.SW$`50causal`, function(x) x$ACAT)<=8.33e-6)/1000), "Prop" = c(rep(100,2), rep(75,2), rep(50,2)), "Annot" = rep(c("Burden","ACAT-Combined"),3), "Annotation"="SW")

df_power_CRHs_Combined_2causal_OR5_SW$Annot = factor(df_power_CRHs_Combined_2causal_OR5_SW$Annot, levels=c("Burden","ACAT-Combined"))



df_power_CRHs_Combined_2causal_OR5_CRHs = data.frame("Power"=c(sum(sapply(pvalues.alter.2causal.CRHs$`100causal`, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                             sum(sapply(pvalues.alter.2causal.CRHs$`100causal`, function(x) x$ACAT)<=8.33e-6)/1000,
                                                             
                                                             sum(sapply(pvalues.alter.2causal.CRHs$`75causal`, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                             sum(sapply(pvalues.alter.2causal.CRHs$`75causal`, function(x) x$ACAT)<=8.33e-6)/1000,
                                                             
                                                             sum(sapply(pvalues.alter.2causal.CRHs$`50causal`, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                             sum(sapply(pvalues.alter.2causal.CRHs$`50causal`, function(x) x$ACAT)<=8.33e-6)/1000), "Prop" = c(rep(100,2), rep(75,2), rep(50,2)), "Annot" = rep(c("Burden","ACAT-Combined"),3), "Annotation"="CRHs")

df_power_CRHs_Combined_2causal_OR5_CRHs$Annot = factor(df_power_CRHs_Combined_2causal_OR5_CRHs$Annot, levels=c("Burden","ACAT-Combined"))


df_power_CRHs_Combined_2causal_OR5_Genes = data.frame("Power"=c(sum(sapply(pvalues.alter.2causal.Genes$`100causal`, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                             sum(sapply(pvalues.alter.2causal.Genes$`100causal`, function(x) x$ACAT)<=8.33e-6)/1000,
                                                             
                                                             sum(sapply(pvalues.alter.2causal.Genes$`75causal`, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                             sum(sapply(pvalues.alter.2causal.Genes$`75causal`, function(x) x$ACAT)<=8.33e-6)/1000,
                                                             
                                                             sum(sapply(pvalues.alter.2causal.Genes$`50causal`, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                             sum(sapply(pvalues.alter.2causal.Genes$`50causal`, function(x) x$ACAT)<=8.33e-6)/1000), "Prop" = c(rep(100,2), rep(75,2), rep(50,2)), "Annot" = rep(c("Burden","ACAT-Combined"),3), "Annotation"="Genes")

df_power_CRHs_Combined_2causal_OR5_Genes$Annot = factor(df_power_CRHs_Combined_2causal_OR5_Genes$Annot, levels=c("Burden","ACAT-Combined"))

df_power_CRHs_Combined_2causal_OR5_Pairs = data.frame("Power"=c(sum(sapply(pvalues.alter.2causal.Pairs$`100causal`, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                             sum(sapply(pvalues.alter.2causal.Pairs$`100causal`, function(x) x$ACAT)<=8.33e-6)/1000,
                                                             
                                                             sum(sapply(pvalues.alter.2causal.Pairs$`75causal`, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                             sum(sapply(pvalues.alter.2causal.Pairs$`75causal`, function(x) x$ACAT)<=8.33e-6)/1000,
                                                             
                                                             sum(sapply(pvalues.alter.2causal.Pairs$`50causal`, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                             sum(sapply(pvalues.alter.2causal.Pairs$`50causal`, function(x) x$ACAT)<=8.33e-6)/1000), "Prop" = c(rep(100,2), rep(75,2), rep(50,2)), "Annot" = rep(c("Burden","ACAT-Combined"),3), "Annotation"="Pairs")

df_power_CRHs_Combined_2causal_OR5_Pairs$Annot = factor(df_power_CRHs_Combined_2causal_OR5_Pairs$Annot, levels=c("Burden","ACAT-Combined"))



panel4 = ggarrange(ggarrange(ggplot(df_power_CRHs_2causal_OR5, aes(x=Type, y=Power, fill=Type))+geom_bar(position = "dodge", stat="identity",color="black")+theme_bw()+theme(legend.position="none",text = element_text(size = 15))+ylim(c(0,1))+facet_grid(.~Prop)+xlab("Score Type")+ylab("Power"),
                             ggplot(rbind(df_power_CRHs_Combined_2causal_OR5[df_power_CRHs_Combined_2causal_OR5$Annot=="ACAT-Combined",],df_power_CRHs_Combined_2causal_OR5_pairs[df_power_CRHs_Combined_2causal_OR5_pairs$Annot=="ACAT-Combined",],
                                          df_power_CRHs_Combined_2causal_OR5_genes[df_power_CRHs_Combined_2causal_OR5_genes$Annot=="ACAT-Combined",]), aes(x=Annotation,y=Power,fill=Annotation))+geom_bar(stat="identity",position="dodge", colour="black")+facet_grid(.~Prop)+theme_bw()+theme(legend.position = "none",text = element_text(size = 15)), labels = c("A","B"), ncol=2),
                   ggarrange(ggplot(df_power_Retro_RVS_RVNPL, aes(x=Type, y=Power, fill=Type))+geom_bar(position = "dodge", stat="identity",color="black")+theme_bw()+theme(legend.position="none",text = element_text(size = 15))+ylim(c(0,1))+xlab("Type")+ylab("Power \nACAT-Combined")+facet_grid(.~Method,scales = "free", space = "free"), labels = "C"), ncol=1, nrow=2)


ggplot(df_power_CRHs_Combined_2causal_OR5_SW, aes(x=Annot, y=Power, fill=Annot))+geom_bar(position = "dodge", stat="identity",color="black")+theme_bw()+theme(legend.position="none",text = element_text(size = 15))+ylim(c(0,1))+xlab("Score Type")+ylab("Power")+facet_grid(.~Prop)

pedfiles_2causal_smallped_100_OR5 = list.files("data\\data_smallped\\power\\2causal\\100causal", full.names=T, recursive = T)
pedfiles_2causal_smallped_75_OR5 = list.files("data\\data_smallped\\power\\2causal\\75causal", full.names=T, recursive=T)
pedfiles_2causal_smallped_50_OR5 = list.files("data\\data_smallped\\power\\2causal\\50causal", full.names=T, recursive=T)

power_2causal_smallped_100_OR5 = lapply(1:1000, function(x) RetroFun.RVS(null,agg_2causal_smallped_100_OR5[[x]], Z,W))
power_2causal_smallped_75_OR5 = lapply(1:1000, function(x) RetroFun.RVS(null,agg_2causal_smallped_75_OR5[[x]], Z,W))
power_2causal_smallped_50_OR5 = lapply(1:1000, function(x) RetroFun.RVS(null,agg_2causal_smallped_50_OR5[[x]], Z,W))

#Figure S10-14: Power at 2%-1% causal variants for competitor methods or different types of functional annotations

power_1causal_100_OR5 = readRDS("data\\pvalues_alter_1causal_100_OR5.RDS")
power_1causal_75_OR5 = readRDS("data\\pvalues_alter_1causal_75_OR5.RDS")
power_1causal_50_OR5 = readRDS("data\\pvalues_alter_1causal_50_OR5.RDS")

df_power_CRHs_Combined_1causal_OR5 = data.frame("Power"=c(sum(sapply(power_1causal_100_OR5, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                          sum(sapply(power_1causal_100_OR5, function(x) x$ACAT)<=8.33e-6)/1000,
                                                          
                                                          sum(sapply(power_1causal_75_OR5, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                          sum(sapply(power_1causal_75_OR5, function(x) x$ACAT)<=8.33e-6)/1000,
                                                          
                                                          sum(sapply(power_1causal_50_OR5, function(x) x$Score_V1)<=8.33e-6)/1000,
                                                          sum(sapply(power_1causal_50_OR5, function(x) x$ACAT)<=8.33e-6)/1000), "Prop" = c(rep(100,2), rep(75,2), rep(50,2)), "Annot" = rep(c("Burden","ACAT-Combined"),3), "Annotation"="CRHs")

df_power_CRHs_Combined_1causal_OR5$Annot = factor(df_power_CRHs_Combined_1causal_OR5$Annot, levels=c("Burden","ACAT-Combined"))
ggplot(df_power_CRHs_Combined_1causal_OR5, aes(x=Annot, y=Power, fill=Annot))+geom_bar(position = "dodge", stat="identity",color="black")+xlab("Type")+theme_bw()+facet_grid(~Prop)+theme(legend.position = "none")



pvalues_CHP_75causal = list.files("data\\results_RVNPL_1causal_75\\CHP", full.names = T)
pvalues_RV_75causal = list.files("data\\results_RVNPL_1causal_75\\RV", full.names = T)

pvalues_ACAT_CHP_75causal_pairs = c()
pvalues_ACAT_CHP_75causal_all = c()

for(i in 1:200){
  print(i)
  df_p = read.table(pvalues_CHP_75causal[i], fill=T, header=F)
  pvalues_ACAT_CHP_75causal_pairs = c(pvalues_ACAT_CHP_75causal_pairs, ACAT::ACAT(df_p[,1]))
  pvalues_ACAT_CHP_75causal_all = c(pvalues_ACAT_CHP_75causal_all, ACAT::ACAT(df_p[,2]))
  
}

sum(pvalues_ACAT_CHP_75causal_pairs<=8.333333e-06)/200
sum(pvalues_ACAT_CHP_75causal_all<=8.333333e-06)/200


pvalues_ACAT_RV_75causal_pairs = c()
pvalues_ACAT_RV_75causal_all = c()

for(i in 1:200){
  df_p = read.table(pvalues_RV_75causal[i], fill=T, header=F)
  pvalues_ACAT_RV_75causal_pairs = c(pvalues_ACAT_RV_75causal_pairs, ACAT::ACAT(df_p[,1]))
  pvalues_ACAT_RV_75causal_all = c(pvalues_ACAT_RV_75causal_all, ACAT::ACAT(df_p[,2]))
  
}

sum(pvalues_ACAT_RV_75causal_pairs<=8.333333e-06)/200
sum(pvalues_ACAT_RV_75causal_all<=8.333333e-06)/200


pvalues.alter.1causal.CRHs = readRDS("data\\pvalues_alter_1causal_CRHs.RDS")
pvalues.alter.1causal.Genes = readRDS("data\\pvalues_alter_1causal_Genes.RDS")
pvalues.alter.1causal.Pairs = readRDS("data\\pvalues_alter_1causal_Pairs.RDS")
pvalues.alter.1causal.SW = readRDS("data\\pvalues_alter_1causal_SW.RDS")

df_power_Retro_RVS_RVNPL_1causal = data.frame("Power" = c(sum(sapply(pvalues.alter.1causal.CRHs$`75causal`, function(x) x$ACAT)<=8.333333e-06)/1000,
                                                  sum(sapply(pvalues.alter.1causal.Genes$`75causal`, function(x) x$ACAT)<=8.333333e-06)/1000,
                                                  power.1pc.6000[2,2],power.1pc.6000[1,2],
                                                  sum(pvalues_ACAT_RV_75causal_pairs<=8.333333e-06)/200,
                                                  sum(pvalues_ACAT_RV_75causal_all<=8.333333e-06)/200,
                                                  sum(pvalues_ACAT_CHP_75causal_pairs<=8.333333e-06)/200,
                                                  sum(pvalues_ACAT_CHP_75causal_all<=8.333333e-06)/200), "Method" = c("RetroFun-RVS", "RetroFun-RVS",
                                                                                                                      "RVS", "RVS",
                                                                                                                      "RV-NPL", "RV-NPL", "CHP-NPL", "CHP-NPL"), 
                                      "Type" = c("CRHs","Genes", "Complete\nSharing", "Partial\nSharing", "Pairs", "All", "Pairs", "All"))

df_power_Retro_RVS_RVNPL_1causal$Method = factor(df_power_Retro_RVS_RVNPL_1causal$Method, levels=c("RetroFun-RVS", "RV-NPL", "CHP-NPL","RVS"))

ggplot(df_power_Retro_RVS_RVNPL_1causal, aes(x=Type, y=Power, fill=Type))+geom_bar(position = "dodge", stat="identity",color="black")+theme_bw()+theme(legend.position="none",text = element_text(size = 15))+ylim(c(0,1))+xlab("Type")+ylab("Power \nACAT-Combined")+facet_grid(.~Method,scales = "free", space = "free")
