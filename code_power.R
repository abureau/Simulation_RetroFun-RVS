#Power Analysis
#Test alter 100% within first CRH
ped_files_alter_6causal = list.files("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data\\Scenario_100_6perccausal_oneCRH", full.names=T)
ped_files_alter_2causal = list.files("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data\\Scenario_100_2perccausal_oneCRH", full.names=T)
ped_files_alter_4causal = list.files("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data\\Scenario_100_4perccausal_oneCRH", full.names=T)

sfs_100_06 = read.table("rare_variants_scenario100_0.06.sfs", header=F)
sfs_100_04 = read.table("rare_variants_scenario100_0.04.sfs", header=F)
sfs_100_02 = read.table("rare_variants_scenario100_0.02.sfs", header=F)

table(sfs_100_06$V1,sfs_100_06$V7)
#CRH/Causal   0   1
#1045         128  30
#1257         64   1
#1335         1   1
#988          40   1
#Out          243   1

#Code to generate Z as binary matrix
u = unique(sfs_100_06$V1)
Z = matrix(0, nrow=nrow(sfs_100_06),ncol=length(u)-1)
Z_tmp = sapply(1:ncol(Z),function(x){
  Z_col = Z[,x,drop=F]
  w = which(sfs_100_06$V1==u[x])
  Z_col[w,] = 1
  Z_col
})
Z = cbind(1,Z_tmp)

W = diag(dbeta(sfs_100_06$V6,1,20),510)



list_6_100 = list()
list_4_100 = list()
list_2_100 = list()

for(i in 1:100){
  print(i)
  agg_test_06 = aggregate.geno.by.fam(ped_files_alter_6causal[i], correction = "replace.homo", FamID = null$FamID)
  agg_test_04 = aggregate.geno.by.fam(ped_files_alter_4causal[i], correction = "replace.homo", FamID = null$FamID)
  agg_test_02 = aggregate.geno.by.fam(ped_files_alter_2causal[i], correction = "replace.homo", FamID = null$FamID)
  
  list_6_100[[i]] = retrofun.RVS(null,agg_test_06,Z,W)
  list_4_100[[i]] = retrofun.RVS(null,agg_test_04,Z,W)
  list_2_100[[i]] = retrofun.RVS(null,agg_test_02,Z,W)
}

sum(sapply(list_6_100,function(x) x$Score1)<=0.05)/100
sum(sapply(list_6_100,function(x) ACAT(c(x$Score1,x$Score2)))<=0.05)/100
sum(sapply(list_6_100,function(x) ACAT(c(x$Score1,x$Score2,x$Score3)))<=0.05)/100
sum(sapply(list_6_100,function(x) x$ACAT<=0.05))/100


sum(sapply(list_4_100,function(x) x$Score1)<=0.05)/100
sum(sapply(list_4_100,function(x) ACAT(c(x$Score1,x$Score2)))<=0.05)/100
sum(sapply(list_4_100,function(x) ACAT(c(x$Score1,x$Score2,x$Score3)))<=0.05)/100
sum(sapply(list_4_100,function(x) x$ACAT<=0.05))/100

sum(sapply(list_2_100,function(x) x$Score1)<=0.05)/100
sum(sapply(list_2_100,function(x) ACAT(c(x$Score1,x$Score2)))<=0.05)/100
sum(sapply(list_2_100,function(x) ACAT(c(x$Score1,x$Score2,x$Score3)))<=0.05)/100
sum(sapply(list_2_100,function(x) x$ACAT<=0.05))/100

df_6_100 = ldply(list_6_100,"rbind")
df_6_100$Prop = 0.06
df_4_100 = ldply(list_4_100,"rbind")
df_4_100$Prop = 0.04
df_2_100 = ldply(list_2_100,"rbind")
df_2_100$Prop = 0.02

df_all_100 = rbind(df_6_100,df_4_100,df_2_100)

#Power by Annotation
aggregate(Score1~Prop,df_all_100, function(x) sum(x<0.05, na.rm = T)/100)
aggregate(Score2~Prop,df_all_100, function(x) sum(x<0.05,na.rm = T)/100)
aggregate(Score3~Prop,df_all_100, function(x) sum(x<0.05,na.rm = T)/100)
aggregate(Score4~Prop,df_all_100, function(x) sum(x<0.05,na.rm = T)/100)
aggregate(Score5~Prop,df_all_100, function(x) sum(x<0.05,na.rm = T)/100)

df_all_100$ACAT1_2 = apply(df_all_100[,c("Score1", "Score2")],1,ACAT)

df_all_100$ACAT1_2_3 = apply(df_all_100[,c("Score1", "Score2","Score3")],1,function(r){
  if(any(is.na(r))){
    ACAT(r[-which(is.na(r))])
  } else {
    ACAT(r)
  }
})

df_all_100$ACAT1_2_3_4 = apply(df_all_100[,c("Score1", "Score2","Score3", "Score4")],1,function(r){
  if(any(is.na(r))){
    ACAT(r[-which(is.na(r))])
  } else {
    ACAT(r)
  }
})

Init_power = aggregate(Score1~Prop,df_all_100, function(x) sum(x<0.05, na.rm = T)/100)
Init_power$N_Annot = 0
colnames(Init_power) = c("Prop", "Power", "N_Annot")
Informative_Annot_power = aggregate(ACAT1_2~Prop,df_all_100, function(x) sum(x<0.05,na.rm = T)/100)
Informative_Annot_power$N_Annot = 1
Uninformative_Annot1_power = aggregate(ACAT1_2_3~Prop,df_all_100, function(x) sum(x<0.05,na.rm = T)/100)
Uninformative_Annot1_power$N_Annot = 2 
Uninformative_Annot2_power = aggregate(ACAT1_2_3_4~Prop,df_all_100, function(x) sum(x<0.05,na.rm = T)/100)
Uninformative_Annot2_power$N_Annot = 3
Uninformative_Annot3_power = aggregate(ACAT~Prop,df_all_100, function(x) sum(x<0.05,na.rm = T)/100)
Uninformative_Annot3_power$N_Annot = 4

colnames(Informative_Annot_power) = c("Prop", "Power", "N_Annot")
colnames(Uninformative_Annot1_power) = c("Prop", "Power", "N_Annot")
colnames(Uninformative_Annot2_power) = c("Prop", "Power", "N_Annot")
colnames(Uninformative_Annot3_power) = c("Prop", "Power", "N_Annot")
power_add_Annot = rbind(Init_power,Informative_Annot_power,Uninformative_Annot1_power,Uninformative_Annot2_power, Uninformative_Annot3_power)

ggplot(power_add_Annot, aes(x=N_Annot, y=Power, fill=factor(Prop)))+geom_bar(stat = "identity", position = "dodge", colour="black")+xlab("# Functional Annotations")+labs(fill="% Causal")

#Compare variances for each functional annotations at different proportion of causal variants
list_var_6 = list()
list_var_4 = list()
list_var_2 = list()

for(i in 1:100){
  agg_test_02_test = aggregate.geno.by.fam(ped_files_alter_2causal[i], correction = "replace.homo", FamID = null$FamID)
  agg_test_04_test = aggregate.geno.by.fam(ped_files_alter_4causal[i], correction = "replace.homo", FamID = null$FamID)
  agg_test_06_test = aggregate.geno.by.fam(ped_files_alter_6causal[i], correction = "replace.homo", FamID = null$FamID)
  
  
  list_var_6[[i]] = unlist(compute.Var.by.annot(null,agg_test_06_test,Z, W))
  list_var_4[[i]] = unlist(compute.Var.by.annot(null,agg_test_04_test,Z, W))
  list_var_2[[i]] = unlist(compute.Var.by.annot(null,agg_test_02_test,Z, W))
}

sum(sapply(1:100, function(x) {
  list_var_6[[x]][2]>list_var_4[[x]][2]}))

ped_files_alter_6causal_2CRHs = list.files("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data\\Scenario100_2CRHs_6perc_agg", full.names=T)
ped_files_alter_4causal_2CRHs = list.files("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data\\Scenario100_2CRHs_4perc_agg", full.names=T)
ped_files_alter_2causal_2CRHs = list.files("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data\\Scenario100_2CRHs_2perc_agg", full.names=T)

list_06_100_2CRHs = list()
list_4_100_2CRHs = list()
list_2_100_2CRHs = list()

for(i in 1:100){
  print(i)
  agg_test_06 = aggregate.geno.by.fam(ped_files_alter_6causal_2CRHs[i], correction = "replace.homo", FamID = null$FamID)
  agg_test_04 = aggregate.geno.by.fam(ped_files_alter_4causal_2CRHs[i], correction = "replace.homo", FamID = null$FamID)
  agg_test_02 = aggregate.geno.by.fam(ped_files_alter_2causal_2CRHs[i], correction = "replace.homo", FamID = null$FamID)
  
  list_06_100_2CRHs[[i]] = retrofun.RVS(null,agg_test_06,Z,W)
  list_4_100_2CRHs[[i]] = retrofun.RVS(null,agg_test_04,Z,W)
  list_2_100_2CRHs[[i]] = retrofun.RVS(null,agg_test_02,Z,W)
}

df_6_100_2CRHs = ldply(list_06_100_2CRHs,"rbind")
df_6_100_2CRHs$Prop = 0.06
df_4_100_2CRHs = ldply(list_4_100_2CRHs,"rbind")
df_4_100_2CRHs$Prop = 0.04
df_2_100_2CRHs = ldply(list_2_100_2CRHs,"rbind")
df_2_100_2CRHs$Prop = 0.02
df_2_100_2CRHs$Score4 = NA
df_all_100_2CRHs = rbind(df_6_100_2CRHs,df_4_100_2CRHs,df_2_100_2CRHs)


#Power by Annotation
aggregate(Score1~Prop,df_all_100_2CRHs, function(x) sum(x<0.05, na.rm = T)/100)
aggregate(Score2~Prop,df_all_100_2CRHs, function(x) sum(x<0.05,na.rm = T)/100)
aggregate(Score3~Prop,df_all_100_2CRHs, function(x) sum(x<0.05,na.rm = T)/100)
aggregate(Score4~Prop,df_all_100_2CRHs, function(x) sum(x<0.05,na.rm = T)/100)
aggregate(Score5~Prop,df_all_100_2CRHs, function(x) sum(x<0.05,na.rm = T)/100)

df_all_100_2CRHs$ACAT1_2 = apply(df_all_100_2CRHs[,c("Score1", "Score2")],1,ACAT)
df_all_100_2CRHs$ACAT1_2_3 = apply(df_all_100_2CRHs[,c("Score1", "Score2","Score3")],1,function(r){
  if(any(is.na(r))){
    ACAT(r[-which(is.na(r))])
  } else {
    ACAT(r)
  }
})
df_all_100_2CRHs$ACAT1_2_3_4 = apply(df_all_100_2CRHs[,c("Score1", "Score2","Score3", "Score4")],1,function(r){
  if(any(is.na(r))){
    ACAT(r[-which(is.na(r))])
  } else {
    ACAT(r)
  }
})

Init_power_2CRHs = aggregate(Score1~Prop,df_all_100_2CRHs, function(x) sum(x<0.05, na.rm = T)/100)
Init_power_2CRHs$N_Annot = 0
colnames(Init_power_2CRHs) = c("Prop", "Power", "N_Annot")

Informative_Annot_power_2CRHs = aggregate(ACAT1_2~Prop,df_all_100_2CRHs, function(x) sum(x<0.05,na.rm = T)/100)
Informative_Annot_power_2CRHs$N_Annot = 1

Uninformative_Annot1_power_2CRHs = aggregate(ACAT1_2_3~Prop,df_all_100_2CRHs, function(x) sum(x<0.05,na.rm = T)/100)
Uninformative_Annot1_power_2CRHs$N_Annot = 2 

Uninformative_Annot2_power_2CRHs = aggregate(ACAT1_2_3_4~Prop,df_all_100_2CRHs, function(x) sum(x<0.05,na.rm = T)/100)
Uninformative_Annot2_power_2CRHs$N_Annot = 3

Uninformative_Annot3_power_2CRHs = aggregate(ACAT~Prop,df_all_100_2CRHs, function(x) sum(x<0.05,na.rm = T)/100)
Uninformative_Annot3_power_2CRHs$N_Annot = 4

colnames(Informative_Annot_power_2CRHs) = c("Prop", "Power", "N_Annot")
colnames(Uninformative_Annot1_power_2CRHs) = c("Prop", "Power", "N_Annot")
colnames(Uninformative_Annot2_power_2CRHs) = c("Prop", "Power", "N_Annot")
colnames(Uninformative_Annot3_power_2CRHs) = c("Prop", "Power", "N_Annot")

power_add_Annot_2CRHs = rbind(Init_power_2CRHs,Informative_Annot_power_2CRHs,Uninformative_Annot1_power_2CRHs,Uninformative_Annot2_power_2CRHs, Uninformative_Annot3_power_2CRHs)

ggplot(power_add_Annot_2CRHs, aes(x=N_Annot, y=Power, fill=factor(Prop)))+geom_bar(stat = "identity", position = "dodge", colour="black")+xlab("# Functional Annotations")+labs(fill="% Causal")

sfs_100_02_2CRHs = read.table("rare_variants_scenario100_0.02_2CRHs.sfs", header=F)
sfs_100_04_2CRHs = read.table("rare_variants_scenario100_0.04_2CRHs.sfs", header=F)
sfs_100_06_2CRHs = read.table("rare_variants_scenario100_0.06_2CRHs.sfs", header=F)

table(sfs_100_02_2CRHs$V1,sfs_100_02_2CRHs$V7)
table(sfs_100_04_2CRHs$V1,sfs_100_04_2CRHs$V7)
table(sfs_100_06_2CRHs$V1,sfs_100_06_2CRHs$V7)

test_agg_02 = aggregate.geno.by.fam(ped_files_alter_2causal[45], correction = "replace.homo", FamID = null$FamID)
test_agg_04 = aggregate.geno.by.fam(ped_files_alter_4causal[45], correction = "replace.homo", FamID = null$FamID)
test_agg_06 = aggregate.geno.by.fam(ped_files_alter_6causal[45], correction = "replace.homo", FamID = null$FamID)

test_agg_02_2CRHs = aggregate.geno.by.fam(ped_files_alter_2causal_2CRHs[45], correction = "replace.homo", FamID = null$FamID)
test_agg_04_2CRHs = aggregate.geno.by.fam(ped_files_alter_4causal_2CRHs[45], correction = "replace.homo", FamID = null$FamID)
test_agg_06_2CRHs = aggregate.geno.by.fam(ped_files_alter_6causal_2CRHs[45], correction = "replace.homo", FamID = null$FamID)

diag(W[test_agg_02_2CRHs$index_variants,test_agg_02_2CRHs$index_variants])
diag(W[test_agg_04_2CRHs$index_variants,test_agg_04_2CRHs$index_variants])
diag(W[test_agg_06_2CRHs$index_variants,test_agg_06_2CRHs$index_variants])

colSums(Z[test_agg_02_2CRHs$index_variants,])
colSums(Z[test_agg_04_2CRHs$index_variants,])
colSums(Z[test_agg_06_2CRHs$index_variants,])

table(sfs_100_02_2CRHs[test_agg_02_2CRHs$index_variants,"V1"],sfs_100_02_2CRHs[test_agg_02_2CRHs$index_variants,"V7"])
table(sfs_100_04_2CRHs[test_agg_04_2CRHs$index_variants,"V1"],sfs_100_04_2CRHs[test_agg_04_2CRHs$index_variants,"V7"])
table(sfs_100_06_2CRHs[test_agg_06_2CRHs$index_variants,"V1"],sfs_100_06_2CRHs[test_agg_06_2CRHs$index_variants,"V7"])

table(sfs_100_02[test_agg_02$index_variants,"V1"],sfs_100_02_2CRHs[test_agg_02$index_variants,"V7"])
table(sfs_100_04[test_agg_04$index_variants,"V1"],sfs_100_04[test_agg_04$index_variants,"V7"])
table(sfs_100_06[test_agg_06$index_variants,"V1"],sfs_100_06[test_agg_06$index_variants,"V7"])

ped_files_alter_6causal_75 = list.files("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data\\Scenario75_06_1CRH_agg", full.names=T)
ped_files_alter_4causal_75 = list.files("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data\\Scenario75_04_1CRH_agg",full.names=T)
ped_files_alter_2causal_75 = list.files("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data\\Scenario75_02_1CRH_agg",full.names=T)

list_6_75 = list()
list_4_75 = list()
list_2_75 = list()

for(i in 1:100){
  print(i)
  agg_test_06 = aggregate.geno.by.fam(ped_files_alter_6causal_75[i], correction = "replace.homo", FamID = null$FamID)
  agg_test_04 = aggregate.geno.by.fam(ped_files_alter_4causal_75[i], correction = "replace.homo", FamID = null$FamID)
  agg_test_02 = aggregate.geno.by.fam(ped_files_alter_2causal_75[i], correction = "replace.homo", FamID = null$FamID)
  
  list_6_75[[i]] = retrofun.RVS(null,agg_test_06,Z,W)
  list_4_75[[i]] = retrofun.RVS(null,agg_test_04,Z,W)
  list_2_75[[i]] = retrofun.RVS(null,agg_test_02,Z,W)
}

df_6_75 = ldply(list_6_75,"rbind")
df_6_75$Prop = 0.06
df_4_75 = ldply(list_4_75,"rbind")
df_4_75$Prop = 0.04
df_4_75$Score4 = NA
df_2_75 = ldply(list_2_75,"rbind")
df_2_75$Prop = 0.02

df_all_75 = rbind(df_6_75,df_4_75,df_2_75)

#Power by Annotation
aggregate(Score1~Prop,df_all_75, function(x) sum(x<0.05, na.rm = T)/100)
aggregate(Score2~Prop,df_all_75, function(x) sum(x<0.05,na.rm = T)/100)
aggregate(Score3~Prop,df_all_75, function(x) sum(x<0.05,na.rm = T)/100)
aggregate(Score4~Prop,df_all_75, function(x) sum(x<0.05,na.rm = T)/100)
aggregate(Score5~Prop,df_all_75, function(x) sum(x<0.05,na.rm = T)/100)

df_all_75$ACAT1_2 = apply(df_all_75[,c("Score1", "Score2")],1,ACAT)
df_all_75$ACAT1_2_3 = apply(df_all_75[,c("Score1", "Score2","Score3")],1,function(r){
  if(any(is.na(r))){
    ACAT(r[-which(is.na(r))])
  } else {
    ACAT(r)
  }
})
df_all_75$ACAT1_2_3_4 = apply(df_all_75[,c("Score1", "Score2","Score3", "Score4")],1,function(r){
  if(any(is.na(r))){
    ACAT(r[-which(is.na(r))])
  } else {
    ACAT(r)
  }
})

Init_power_75 = aggregate(Score1~Prop,df_all_75, function(x) sum(x<0.05, na.rm = T)/100)
Init_power_75$N_Annot = 0
colnames(Init_power_75) = c("Prop", "Power", "N_Annot")
Informative_Annot_power_75 = aggregate(ACAT1_2~Prop,df_all_75, function(x) sum(x<0.05,na.rm = T)/100)
Informative_Annot_power_75$N_Annot = 1
Uninformative_Annot1_power_75 = aggregate(ACAT1_2_3~Prop,df_all_75, function(x) sum(x<0.05,na.rm = T)/100)
Uninformative_Annot1_power_75$N_Annot = 2 
Uninformative_Annot2_power_75 = aggregate(ACAT1_2_3_4~Prop,df_all_75, function(x) sum(x<0.05,na.rm = T)/100)
Uninformative_Annot2_power_75$N_Annot = 3
Uninformative_Annot3_power_75 = aggregate(ACAT~Prop,df_all_75, function(x) sum(x<0.05,na.rm = T)/100)
Uninformative_Annot3_power_75$N_Annot = 4

colnames(Informative_Annot_power_75) = c("Prop", "Power", "N_Annot")
colnames(Uninformative_Annot1_power_75) = c("Prop", "Power", "N_Annot")
colnames(Uninformative_Annot2_power_75) = c("Prop", "Power", "N_Annot")
colnames(Uninformative_Annot3_power_75) = c("Prop", "Power", "N_Annot")
power_add_Annot_75 = rbind(Init_power_75,Informative_Annot_power_75,Uninformative_Annot1_power_75,Uninformative_Annot2_power_75, Uninformative_Annot3_power_75)

ggplot(power_add_Annot_75, aes(x=N_Annot, y=Power, fill=factor(Prop)))+geom_bar(stat = "identity", position = "dodge", colour="black")+xlab("# Functional Annotations")+labs(fill="% Causal")

sfs_75_6 = read.table("rare_variants_scenario75_0.06.sfs", header=F)
table(sfs_75_6$V1, sfs_75_6$V7)
sfs_75_4 = read.table("rare_variants_scenario75_0.04.sfs", header=F)
table(sfs_75_4$V1, sfs_75_4$V7)
sfs_75_2 = read.table("rare_variants_scenario75_0.02.sfs", header=F)
table(sfs_75_2$V1, sfs_75_2$V7)



ped_files_alter_6causal_50 = list.files("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data\\Scenario50_06_1CRH_agg", full.names=T)
ped_files_alter_4causal_50 = list.files("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data\\Scenario50_04_1CRH_agg",full.names=T)
ped_files_alter_2causal_50 = list.files("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data\\Scenario50_02_1CRH_agg",full.names=T)

list_6_50 = list()
list_4_50 = list()
list_2_50 = list()

for(i in 1:100){
  print(i)
  agg_test_06 = aggregate.geno.by.fam(ped_files_alter_6causal_50[i], correction = "replace.homo", FamID = null$FamID)
  agg_test_04 = aggregate.geno.by.fam(ped_files_alter_4causal_50[i], correction = "replace.homo", FamID = null$FamID)
  agg_test_02 = aggregate.geno.by.fam(ped_files_alter_2causal_50[i], correction = "replace.homo", FamID = null$FamID)
  
  list_6_50[[i]] = retrofun.RVS(null,agg_test_06,Z,W)
  list_4_50[[i]] = retrofun.RVS(null,agg_test_04,Z,W)
  list_2_50[[i]] = retrofun.RVS(null,agg_test_02,Z,W)
}

df_6_50 = ldply(list_6_50,"rbind")
df_6_50$Prop = 0.06
df_4_50 = ldply(list_4_50,"rbind")
df_4_50$Prop = 0.04
df_2_50 = ldply(list_2_50,"rbind")
df_2_50$Prop = 0.02

df_all_50 = rbind(df_6_50,df_4_50,df_2_50)


#Power by Annotation
aggregate(Score1~Prop,df_all_50, function(x) sum(x<0.05, na.rm = T)/100)
aggregate(Score2~Prop,df_all_50, function(x) sum(x<0.05,na.rm = T)/100)
aggregate(Score3~Prop,df_all_50, function(x) sum(x<0.05,na.rm = T)/100)
aggregate(Score4~Prop,df_all_50, function(x) sum(x<0.05,na.rm = T)/100)
aggregate(Score5~Prop,df_all_50, function(x) sum(x<0.05,na.rm = T)/100)

df_all_50$ACAT1_2 = apply(df_all_50[,c("Score1", "Score2")],1,ACAT)
df_all_50$ACAT1_2_3 = apply(df_all_50[,c("Score1", "Score2","Score3")],1,function(r){
  if(any(is.na(r))){
    ACAT(r[-which(is.na(r))])
  } else {
    ACAT(r)
  }
})
df_all_50$ACAT1_2_3_4 = apply(df_all_50[,c("Score1", "Score2","Score3", "Score4")],1,function(r){
  if(any(is.na(r))){
    ACAT(r[-which(is.na(r))])
  } else {
    ACAT(r)
  }
})

Init_power_50 = aggregate(Score1~Prop,df_all_50, function(x) sum(x<0.05, na.rm = T)/100)
Init_power_50$N_Annot = 0
colnames(Init_power_50) = c("Prop", "Power", "N_Annot")
Informative_Annot_power_50 = aggregate(ACAT1_2~Prop,df_all_50, function(x) sum(x<0.05,na.rm = T)/100)
Informative_Annot_power_50$N_Annot = 1
Uninformative_Annot1_power_50 = aggregate(ACAT1_2_3~Prop,df_all_50, function(x) sum(x<0.05,na.rm = T)/100)
Uninformative_Annot1_power_50$N_Annot = 2 
Uninformative_Annot2_power_50 = aggregate(ACAT1_2_3_4~Prop,df_all_50, function(x) sum(x<0.05,na.rm = T)/100)
Uninformative_Annot2_power_50$N_Annot = 3
Uninformative_Annot3_power_50 = aggregate(ACAT~Prop,df_all_50, function(x) sum(x<0.05,na.rm = T)/100)
Uninformative_Annot3_power_50$N_Annot = 4

colnames(Informative_Annot_power_50) = c("Prop", "Power", "N_Annot")
colnames(Uninformative_Annot1_power_50) = c("Prop", "Power", "N_Annot")
colnames(Uninformative_Annot2_power_50) = c("Prop", "Power", "N_Annot")
colnames(Uninformative_Annot3_power_50) = c("Prop", "Power", "N_Annot")
power_add_Annot_50 = rbind(Init_power_50,Informative_Annot_power_50,Uninformative_Annot1_power_50,Uninformative_Annot2_power_50, Uninformative_Annot3_power_50)

ggplot(power_add_Annot_50, aes(x=N_Annot, y=Power, fill=factor(Prop)))+geom_bar(stat = "identity", position = "dodge", colour="black")+xlab("# Functional Annotations")+labs(fill="% Causal")

sfs_50_06 = read.table("rare_variants_scenario50_0.06.sfs", header=F)
sfs_50_04 = read.table("rare_variants_scenario50_0.04.sfs", header=F)
sfs_50_02 = read.table("rare_variants_scenario50_0.02.sfs", header=F)

list_check_50_06 = list()
list_check_50_04 = list()
list_check_50_02 = list()

for(i in 1:100){
  
  agg_test_06 = aggregate.geno.by.fam(ped_files_alter_6causal_50[i], correction = "replace.homo", FamID = null$FamID)
  agg_test_04 = aggregate.geno.by.fam(ped_files_alter_4causal_50[i], correction = "replace.homo", FamID = null$FamID)
  agg_test_02 = aggregate.geno.by.fam(ped_files_alter_2causal_50[i], correction = "replace.homo", FamID = null$FamID)
  
  
  
  
}

sharing.causal.global.100.02 = sapply(1:100, function(x) sharing.causal.variants.by.fam(ped_files_alter_2causal[[x]],sfs_100_02,c("1045", "1257", "1335", "988","Out")))
sharing.causal.global.100.04 = sapply(1:100, function(x) sharing.causal.variants.by.fam(ped_files_alter_4causal[[x]],sfs_100_04,c("1045", "1257", "1335", "988","Out")))
sharing.causal.global.100.06 = sapply(1:100, function(x) sharing.causal.variants.by.fam(ped_files_alter_6causal[[x]],sfs_100_06,c("1045", "1257", "1335", "988","Out")))

sharing.causal.CRH.100.02 = sapply(1:100, function(x) sharing.causal.variants.by.fam(ped_files_alter_2causal[[x]],sfs_100_02,c("1045")))
sharing.causal.CRH.100.04 = sapply(1:100, function(x) sharing.causal.variants.by.fam(ped_files_alter_4causal[[x]],sfs_100_04,c("1045")))
sharing.causal.CRH.100.06 = sapply(1:100, function(x) sharing.causal.variants.by.fam(ped_files_alter_6causal[[x]],sfs_100_06,c("1045")))



sfs_50_06 = read.table("rare_variants_scenario50_0.06.sfs", header=F)
sfs_50_04 = read.table("rare_variants_scenario50_0.04.sfs", header=F)
sfs_50_02 = read.table("rare_variants_scenario50_0.02.sfs", header=F)

table(sfs_50_06$V1, sfs_50_06$V7)
table(sfs_50_04$V1, sfs_50_04$V7)
table(sfs_50_02$V1, sfs_50_02$V7)

prop_sharing_4causal_50 = mean(sapply(1:100, function(x) sharing.causal.variants.by.fam(ped_files_alter_4causal_50[x], sfs_50_04, "1045")))
prop_sharing_2causal_50 = mean(sapply(1:100, function(x) sharing.causal.variants.by.fam(ped_files_alter_2causal_50[x], sfs_50_02, "1045")))
prop_sharing_6causal_50 = mean(sapply(1:100, function(x) sharing.causal.variants.by.fam(ped_files_alter_6causal_50[x], sfs_50_06, "1045")))

sfs_75_06 = read.table("rare_variants_scenario75_0.06.sfs", header=F)
sfs_75_04 = read.table("rare_variants_scenario75_0.04.sfs", header=F)
sfs_75_02 = read.table("rare_variants_scenario75_0.02.sfs", header=F)

table(sfs_75_06$V1, sfs_75_06$V7)
table(sfs_75_04$V1, sfs_75_04$V7)
table(sfs_75_02$V1, sfs_75_02$V7)

mean(sapply(1:100, function(x)sharing.causal.variants.by.fam(ped_files_alter_4causal_75[x], sfs_75_04, "1045")))
mean(sapply(1:100, function(x)sharing.causal.variants.by.fam(ped_files_alter_2causal_75[x], sfs_75_02, "1045")))
mean(sapply(1:100, function(x)sharing.causal.variants.by.fam(ped_files_alter_6causal_75[x], sfs_75_06, "1045")))

sfs_100_06 = read.table("rare_variants_scenario100_0.06.sfs", header=F)
sfs_100_04 = read.table("rare_variants_scenario100_0.04.sfs", header=F)
sfs_100_02 = read.table("rare_variants_scenario100_0.02.sfs", header=F)

table(sfs_100_06$V1, sfs_100_06$V7)
table(sfs_100_04$V1, sfs_100_04$V7)
table(sfs_100_02$V1, sfs_100_02$V7)

mean(sapply(1:100, function(x)sharing.causal.variants.by.fam(ped_files_alter_4causal[x], sfs_100_04, "1045")))
mean(sapply(1:100, function(x)sharing.causal.variants.by.fam(ped_files_alter_2causal[x], sfs_100_02, "1045")))
mean(sapply(1:100, function(x)sharing.causal.variants.by.fam(ped_files_alter_6causal[x], sfs_100_06, "1045")))

n.causal.variant.by.annot = function(pedfile, sfs.file, regions.of.interest=NULL){
  agg = aggregate.geno.by.fam(pedfile, correction = "replace.homo")
  
  sfs_sub = sfs.file[agg$index_variants,]
  table(sfs_sub$V1,sfs_sub$V7)
}

prop_in_annot_2causal_50 = lapply(1:100,function(x) prop.table(n.causal.variant.by.annot(ped_files_alter_2causal_50[x], sfs_50_02),margin = 2))
prop_in_annot_4causal_50 = lapply(1:100,function(x) prop.table(n.causal.variant.by.annot(ped_files_alter_4causal_50[x], sfs_50_04),margin = 2))
prop_in_annot_6causal_50 = lapply(1:100,function(x) prop.table(n.causal.variant.by.annot(ped_files_alter_6causal_50[x], sfs_50_06),margin = 2))



Scenario75_0.06_resample = list.files("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data\\Scenario75_0.06_1CRH_randomsampling_agg", full.names=T)
Scenario75_0.04_resample = list.files("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data\\Scenario75_0.04_1CRH_randomsampling_agg", full.names=T)
Scenario75_0.02_resample = list.files("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data\\Scenario75_0.02_1CRH_randomsampling_agg", full.names = T)

list_6_75_resample = list()
list_4_75_resample = list()
list_2_75_resample = list()


for(i in 1:100){
  agg_06 = aggregate.geno.by.fam(Scenario75_0.06_resample[i], correction = "replace.homo", FamID = null$FamID)
  agg_04 = aggregate.geno.by.fam(Scenario75_0.04_resample[i], correction = "replace.homo", FamID = null$FamID)
  agg_02 = aggregate.geno.by.fam(Scenario75_0.02_resample[i], correction = "replace.homo", FamID = null$FamID)
  
  list_6_75_resample[[i]] = retrofun.RVS(null,agg_06,Z,W )
  list_4_75_resample[[i]] = retrofun.RVS(null,agg_04,Z,W )
  list_2_75_resample[[i]] = retrofun.RVS(null,agg_02,Z,W )
  
}

df_6_75_resample = ldply(list_6_75_resample,"rbind")
df_6_75_resample$Prop = 0.06
df_6_75_resample$Score4 = NA
df_4_75_resample = ldply(list_4_75_resample,"rbind")
df_4_75_resample$Prop = 0.04
df_4_75_resample$Score4 = NA
df_2_75_resample = ldply(list_2_75_resample,"rbind")
df_2_75_resample$Prop = 0.02

df_all_75_resample = rbind(df_6_75_resample,df_4_75_resample,df_2_75_resample)

#Power by Annotation
aggregate(Score1~Prop,df_all_75_resample, function(x) sum(x<0.05, na.rm = T)/100)
aggregate(Score2~Prop,df_all_75_resample, function(x) sum(x<0.05,na.rm = T)/100)
aggregate(Score3~Prop,df_all_75_resample, function(x) sum(x<0.05,na.rm = T)/100)
aggregate(Score4~Prop,df_all_75_resample, function(x) sum(x<0.05,na.rm = T)/100)
aggregate(Score5~Prop,df_all_75_resample, function(x) sum(x<0.05,na.rm = T)/100)


df_all_75_resample$ACAT1_2 = apply(df_all_75_resample[,c("Score1", "Score2")],1,ACAT)
df_all_75_resample$ACAT1_2_3 = apply(df_all_75_resample[,c("Score1", "Score2","Score3")],1,function(r){
  if(any(is.na(r))){
    ACAT(r[-which(is.na(r))])
  } else {
    ACAT(r)
  }
})
df_all_75_resample$ACAT1_2_3_4 = apply(df_all_75_resample[,c("Score1", "Score2","Score3", "Score4")],1,function(r){
  if(any(is.na(r))){
    ACAT(r[-which(is.na(r))])
  } else {
    ACAT(r)
  }
})

Init_power_75_resample = aggregate(Score1~Prop,df_all_75_resample, function(x) sum(x<0.05, na.rm = T)/100)
Init_power_75_resample$N_Annot = 0
colnames(Init_power_75_resample) = c("Prop", "Power", "N_Annot")
Informative_Annot_power_75_resample = aggregate(ACAT1_2~Prop,df_all_75_resample, function(x) sum(x<0.05,na.rm = T)/100)
Informative_Annot_power_75_resample$N_Annot = 1
Uninformative_Annot1_power_75_resample = aggregate(ACAT1_2_3~Prop,df_all_75_resample, function(x) sum(x<0.05,na.rm = T)/100)
Uninformative_Annot1_power_75_resample$N_Annot = 2 
Uninformative_Annot2_power_75_resample = aggregate(ACAT1_2_3_4~Prop,df_all_75_resample, function(x) sum(x<0.05,na.rm = T)/100)
Uninformative_Annot2_power_75_resample$N_Annot = 3
Uninformative_Annot3_power_75_resample = aggregate(ACAT~Prop,df_all_75_resample, function(x) sum(x<0.05,na.rm = T)/100)
Uninformative_Annot3_power_75_resample$N_Annot = 4

colnames(Informative_Annot_power_75_resample) = c("Prop", "Power", "N_Annot")
colnames(Uninformative_Annot1_power_75_resample) = c("Prop", "Power", "N_Annot")
colnames(Uninformative_Annot2_power_75_resample) = c("Prop", "Power", "N_Annot")
colnames(Uninformative_Annot3_power_75_resample) = c("Prop", "Power", "N_Annot")
power_add_Annot_75_resample = rbind(Init_power_75_resample,Informative_Annot_power_75_resample,Uninformative_Annot1_power_75_resample,Uninformative_Annot2_power_75_resample, Uninformative_Annot3_power_75_resample)

ggplot(power_add_Annot_75_resample, aes(x=N_Annot, y=Power, fill=factor(Prop)))+geom_bar(stat = "identity", position = "dodge", colour="black")+xlab("# Functional Annotations")+labs(fill="% Causal")



Scenario50_0.06_resample = list.files("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data\\Scenario50_0.06_1CRH_randomsampling_agg", full.names=T)
Scenario50_0.04_resample = list.files("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data\\Scenario50_0.04_1CRH_randomsampling_agg", full.names=T)
Scenario50_0.02_resample = list.files("C:\\Users\\loicm\\Documents\\recherche\\Vraisemblance_retrospective\\Simulation\\data\\Scenario50_0.02_1CRH_randomsampling_agg", full.names = T)

list_6_50_resample = list()
list_4_50_resample = list()
list_2_50_resample = list()


for(i in 1:100){
  print(i)
  agg_06 = aggregate.geno.by.fam(Scenario50_0.06_resample[i], correction = "replace.homo", FamID = null$FamID)
  agg_04 = aggregate.geno.by.fam(Scenario50_0.04_resample[i], correction = "replace.homo", FamID = null$FamID)
  agg_02 = aggregate.geno.by.fam(Scenario50_0.02_resample[i], correction = "replace.homo", FamID = null$FamID)
  
  list_6_50_resample[[i]] = retrofun.RVS(null,agg_06,Z,W )
  list_4_50_resample[[i]] = retrofun.RVS(null,agg_04,Z,W )
  list_2_50_resample[[i]] = retrofun.RVS(null,agg_02,Z,W )
  
}


df_6_50_resample = ldply(list_6_50_resample,"rbind")
df_6_50_resample$Prop = 0.06
df_6_50_resample$Score4 = NA
df_4_50_resample = ldply(list_4_50_resample,"rbind")
df_4_50_resample$Prop = 0.04
df_4_50_resample$Score4 = NA
df_2_50_resample = ldply(list_2_50_resample,"rbind")
df_2_50_resample$Prop = 0.02

df_all_50_resample = rbind(df_6_50_resample,df_4_50_resample,df_2_50_resample)

#Power by Annotation
aggregate(Score1~Prop,df_all_50_resample, function(x) sum(x<0.05, na.rm = T)/100)
aggregate(Score2~Prop,df_all_50_resample, function(x) sum(x<0.05,na.rm = T)/100)
aggregate(Score3~Prop,df_all_50_resample, function(x) sum(x<0.05,na.rm = T)/100)
aggregate(Score4~Prop,df_all_50_resample, function(x) sum(x<0.05,na.rm = T)/100)
aggregate(Score5~Prop,df_all_50_resample, function(x) sum(x<0.05,na.rm = T)/100)


df_all_50_resample$ACAT1_2 = apply(df_all_50_resample[,c("Score1", "Score2")],1,ACAT)
df_all_50_resample$ACAT1_2_3 = apply(df_all_50_resample[,c("Score1", "Score2","Score3")],1,function(r){
  if(any(is.na(r))){
    ACAT(r[-which(is.na(r))])
  } else {
    ACAT(r)
  }
})
df_all_50_resample$ACAT1_2_3_4 = apply(df_all_50_resample[,c("Score1", "Score2","Score3", "Score4")],1,function(r){
  if(any(is.na(r))){
    ACAT(r[-which(is.na(r))])
  } else {
    ACAT(r)
  }
})

Init_power_50_resample = aggregate(Score1~Prop,df_all_50_resample, function(x) sum(x<0.05, na.rm = T)/100)
Init_power_50_resample$N_Annot = 0
colnames(Init_power_50_resample) = c("Prop", "Power", "N_Annot")
Informative_Annot_power_50_resample = aggregate(ACAT1_2~Prop,df_all_50_resample, function(x) sum(x<0.05,na.rm = T)/100)
Informative_Annot_power_50_resample$N_Annot = 1
Uninformative_Annot1_power_50_resample = aggregate(ACAT1_2_3~Prop,df_all_50_resample, function(x) sum(x<0.05,na.rm = T)/100)
Uninformative_Annot1_power_50_resample$N_Annot = 2 
Uninformative_Annot2_power_50_resample = aggregate(ACAT1_2_3_4~Prop,df_all_50_resample, function(x) sum(x<0.05,na.rm = T)/100)
Uninformative_Annot2_power_50_resample$N_Annot = 3
Uninformative_Annot3_power_50_resample = aggregate(ACAT~Prop,df_all_50_resample, function(x) sum(x<0.05,na.rm = T)/100)
Uninformative_Annot3_power_50_resample$N_Annot = 4

colnames(Informative_Annot_power_50_resample) = c("Prop", "Power", "N_Annot")
colnames(Uninformative_Annot1_power_50_resample) = c("Prop", "Power", "N_Annot")
colnames(Uninformative_Annot2_power_50_resample) = c("Prop", "Power", "N_Annot")
colnames(Uninformative_Annot3_power_50_resample) = c("Prop", "Power", "N_Annot")
power_add_Annot_50_resample = rbind(Init_power_50_resample,Informative_Annot_power_50_resample,Uninformative_Annot1_power_50_resample,Uninformative_Annot2_power_50_resample, Uninformative_Annot3_power_50_resample)

ggplot(power_add_Annot_50_resample, aes(x=N_Annot, y=Power, fill=factor(Prop)))+geom_bar(stat = "identity", position = "dodge", colour="black")+xlab("# Functional Annotations")+labs(fill="% Causal")