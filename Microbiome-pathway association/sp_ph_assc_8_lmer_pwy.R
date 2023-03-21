
library(lmerTest)

st<-read.csv('bristol.csv',row.names = 1)
st<-as.matrix.data.frame(st)
#phen_all_use<-read.csv('phen_all_use.csv',row.names = 1,header=T)
phen_all_use<-read.csv('phen_all_use6_all.csv',row.names = 1,header=T)
#data_transformed<-read.csv('data_transformed.csv',row.names = 1,header=T)
data_transformed<-read.csv('pwy_all_use.csv',row.names = 1,header=T)
phen_all_rct_use<-read.csv('phen_all_rct_use22.csv',row.names = 1,header=T)
data_transformed<-as.matrix.data.frame(data_transformed)

for(i in 1:dim(phen_all_use)[1]){
  if(is.na(phen_all_use[i,35])==F){
    if(phen_all_use[i,35]=="Dont_Know"){
      phen_all_use[i,35]<-NA
    }
  }
  if(is.na(phen_all_use[i,36])==F){
    if(phen_all_use[i,36]=="Dont_Know"){
      phen_all_use[i,36]<-NA
    }
  }
  if(is.na(phen_all_use[i,37])==F){
    if(phen_all_use[i,37]=="Dont_Know"){
      phen_all_use[i,37]<-NA
    }
  }
  if(is.na(phen_all_use[i,71])==F){
    if(phen_all_use[i,71]=="Dont_Know"){
      phen_all_use[i,71]<-NA
    }
  }
  if(is.na(phen_all_use[i,72])==F){
    if(phen_all_use[i,72]=="Dont_Know"){
      phen_all_use[i,72]<-NA
    }
  }
  if(is.na(phen_all_use[i,73])==F){
    if(phen_all_use[i,73]=="Dont_Know"){
      phen_all_use[i,73]<-NA
    }
  }
  if(is.na(phen_all_use[i,74])==F){
    if(phen_all_use[i,74]=="Dont_Know"){
      phen_all_use[i,74]<-NA
    }
  }
}

colnames(phen_all_use)[298:304]

amx<-phen_all_use[,c(298:304)]

for(i in 1:dim(amx)[1]){
  if(is.na(phen_all_use[i,32])==F){
  if(as.numeric(as.character(phen_all_use[i,32]))<5.7){
    amx[i,1]<-"Normal"
  }else if((as.numeric(as.character(phen_all_use[i,32]))>=5.7)&(as.numeric(as.character(phen_all_use[i,32]))<=6.4)){
    amx[i,1]<-"A1c(5.7-6.4)"
  }else if((as.numeric(as.character(phen_all_use[i,32]))>=6.5)&(as.numeric(as.character(phen_all_use[i,32]))<=7)){
    amx[i,1]<-"A1c(6.5-7)"
  }else if(as.numeric(as.character(phen_all_use[i,32]))>7){
    amx[i,1]<-"A1c(>7)"
  }
  }else{
    amx[i,1]<-NA
  }
}
  
for(i in 1:dim(amx)[1]){
  if((is.na(phen_all_use[i,35])==F)&(is.na(phen_all_use[i,36])==F)){
    if((as.numeric(as.character(phen_all_use[i,35]))<=119)&(as.numeric(as.character(phen_all_use[i,36]))<=79)){
      amx[i,2]<-"Normal"
    }else if((as.numeric(as.character(phen_all_use[i,35]))<90)&(as.numeric(as.character(phen_all_use[i,36]))<60)){
      amx[i,2]<-"Systolic(<90)&Diastolic(<60)"
    }else if(as.numeric(as.character(phen_all_use[i,35]))>129){
      amx[i,2]<-"Systolic(>129)"
    }else if(as.numeric(as.character(phen_all_use[i,36]))>89){
      amx[i,2]<-"Diastolic>89"
    }else if(((as.numeric(as.character(phen_all_use[i,35]))>119)&(as.numeric(as.character(phen_all_use[i,35]))<=129))&(as.numeric(as.character(phen_all_use[i,36]))>79)){
      amx[i,2]<-"Systolic(>119)&Diastolic(>79)"
    }else if(((as.numeric(as.character(phen_all_use[i,35]))>119)&(as.numeric(as.character(phen_all_use[i,35]))<=129))&(as.numeric(as.character(phen_all_use[i,36]))>89)){
      amx[i,2]<-"Systolic(>119)&Diastolic(>89)"
    }else if((as.numeric(as.character(phen_all_use[i,35]))>129)&(as.numeric(as.character(phen_all_use[i,36]))>89)){
      amx[i,2]<-"Systolic(>129)&Diastolic(>89)"
    }else if((as.numeric(as.character(phen_all_use[i,35]))>119)&(as.numeric(as.character(phen_all_use[i,36]))<=79)){
      amx[i,2]<-"Systolic(>119)&Diastolic(<=79)"
    }else{
      amx[i,2]<-NA
    }
  }else{
    amx[i,2]<-NA
  }
}


for(i in 1:dim(amx)[1]){
  if(is.na(phen_all_use[i,40])==F){
    if((as.numeric(as.character(phen_all_use[i,40]))<=25)&(as.numeric(as.character(phen_all_use[i,40]))>18)){
      amx[i,3]<-"Normal"
    }else if(as.numeric(as.character(phen_all_use[i,40]))<18){
      amx[i,3]<-"BMI(<18)"
    }else if((as.numeric(as.character(phen_all_use[i,40]))>=25)&(as.numeric(as.character(phen_all_use[i,40]))<=30)){
      amx[i,3]<-"BMI(25-30)"
    }else if((as.numeric(as.character(phen_all_use[i,40]))>=30)&(as.numeric(as.character(phen_all_use[i,40]))<=35)){
      amx[i,3]<-"BMI(30-35)"
    }else if(as.numeric(as.character(phen_all_use[i,40]))>=35){
      amx[i,3]<-"BMI(>35)"
    }
  }else{
    amx[i,3]<-NA
  }
}


for(i in 1:dim(amx)[1]){
  if(is.na(phen_all_use[i,71])==F){
    if(as.numeric(as.character(phen_all_use[i,71]))<=100){
      amx[i,4]<-"Normal"
    }else if(as.numeric(as.character(phen_all_use[i,71]))>100){
      amx[i,4]<-"HR(>100)"
    }
  }else{
    amx[i,4]<-NA
  }
}

for(i in 1:dim(amx)[1]){
  if(is.na(phen_all_use[i,72])==F){
    if(as.numeric(as.character(phen_all_use[i,72]))>=0.4){
      amx[i,5]<-"Normal"
    }else if(as.numeric(as.character(phen_all_use[i,72]))<0.4){
      amx[i,5]<-"Perfusion index(<0.4)"
    }
  }else{
    amx[i,5]<-NA
  }
}

for(i in 1:dim(amx)[1]){
  if(is.na(phen_all_use[i,74])==F){
    if(as.numeric(as.character(phen_all_use[i,74]))>=97){
      amx[i,6]<-"Normal"
    }else if(as.numeric(as.character(phen_all_use[i,74]))<97){
      amx[i,6]<-"O2 saturation(<97)"
    }
  }else{
    amx[i,6]<-NA
  }
}

for(i in 1:dim(amx)[1]){
  if(is.na(phen_all_use[i,37])==F){
    if(phen_all_use$gender[i]==0){#for males
      if((as.numeric(as.character(phen_all_use[i,37]))>=13.8)&(as.numeric(as.character(phen_all_use[i,37]))<=17.2)){
        amx[i,7]<-"Normal"
      }else if(as.numeric(as.character(phen_all_use[i,37]))<13.8){
        amx[i,7]<-"Hb_total(low)"
      }else if(as.numeric(as.character(phen_all_use[i,37]))>17.2){
        amx[i,7]<-"Hb_total(elevated)"
      }
    }
    if(phen_all_use$gender[i]==1){#for females
      if((as.numeric(as.character(phen_all_use[i,37]))>=12.1)&(as.numeric(as.character(phen_all_use[i,37]))<=15.1)){
        amx[i,7]<-"Normal"
      }else if(as.numeric(as.character(phen_all_use[i,37]))<12.1){
        amx[i,7]<-"Hb_total(low)"
      }else if(as.numeric(as.character(phen_all_use[i,37]))>15.1){
        amx[i,7]<-"Hb_total(elevated)"
      }
    }
  }else{
    amx[i,7]<-NA
  }
}


phen_all_use[,c(298:304)]<-amx

##A1c factors

esz<-array(0,dim=c(dim(data_transformed)[1],3))
pv_sp<-array(0,dim=c(dim(data_transformed)[1],3)) #p value 
kl<-298

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],8))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,7]<-as.numeric(sp[k,1])
    cleaned_data[k,8]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","sp","village")
  ind_temp<-0
  for(j in 1:dim(cleaned_data)[1]){
    na_score<-0
    for(k in 1:dim(cleaned_data)[2]){
      if(is.na(cleaned_data[j,k])==T){
        na_score<-na_score+1
      }
    }
    if(na_score==0){
      ind_temp<-rbind(ind_temp,j)
    }
  }
  cleaned_data2<-cleaned_data[ind_temp,]
  #ind_mls<-match(phen_all_rct_use$respondent_master_id,phen_all_use$respondent_master_id)
  cleaned_data2<-cbind(cleaned_data2,"health_dangerous"=factor(phen_all_use[ind_temp,kl])) #Phenotype)
  cleaned_data2 <- within(cleaned_data2, health_dangerous <- relevel(health_dangerous, ref = "Normal"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+health_dangerous+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[8:10,5]
  
  effect.size<-coef(mss)[8:10,1]
  pval2<-coef(mss)[8:10,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size1","Effect_size2","Effect_size3","Pvalue1","Pvalue2","Pvalue3")
results_ftrs$FDR1 <- p.adjust(results_ftrs$Pvalue1,method = "BH")
results_ftrs$FDR2 <- p.adjust(results_ftrs$Pvalue2,method = "BH")
results_ftrs$FDR3 <- p.adjust(results_ftrs$Pvalue3,method = "BH")

#length(which(results_ftrs$Pvalue<0.05))
#length(which(results_ftrs$FDR<0.05))


effect_size_phen<-results_ftrs[,1:3]
pval_phen<-results_ftrs[,4:6]
fdr_phen<-results_ftrs[,7:9]

colnames(effect_size_phen)<-rownames(coef(mss))[8:10]

lim<-dim(effect_size_phen)[2]


##BP factors

esz<-array(0,dim=c(dim(data_transformed)[1],4))
pv_sp<-array(0,dim=c(dim(data_transformed)[1],4)) #p value 
kl<-299

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],8))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,7]<-as.numeric(sp[k,1])
    cleaned_data[k,8]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","sp","village")
  ind_temp<-0
  for(j in 1:dim(cleaned_data)[1]){
    na_score<-0
    for(k in 1:dim(cleaned_data)[2]){
      if(is.na(cleaned_data[j,k])==T){
        na_score<-na_score+1
      }
    }
    if(na_score==0){
      ind_temp<-rbind(ind_temp,j)
    }
  }
  cleaned_data2<-cleaned_data[ind_temp,]
  #ind_mls<-match(phen_all_rct_use$respondent_master_id,phen_all_use$respondent_master_id)
  cleaned_data2<-cbind(cleaned_data2,"health_dangerous"=factor(phen_all_use[ind_temp,kl])) #Phenotype)
  cleaned_data2 <- within(cleaned_data2, health_dangerous <- relevel(health_dangerous, ref = "Normal"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+health_dangerous+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[8:11,5]
  
  effect.size<-coef(mss)[8:11,1]
  pval2<-coef(mss)[8:11,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size1","Effect_size2","Effect_size3","Effect_size4","Pvalue1","Pvalue2","Pvalue3","Pvalue4")
results_ftrs$FDR1 <- p.adjust(results_ftrs$Pvalue1,method = "BH")
results_ftrs$FDR2 <- p.adjust(results_ftrs$Pvalue2,method = "BH")
results_ftrs$FDR3 <- p.adjust(results_ftrs$Pvalue3,method = "BH")
results_ftrs$FDR4 <- p.adjust(results_ftrs$Pvalue4,method = "BH")

#length(which(results_ftrs$Pvalue<0.05))
#length(which(results_ftrs$FDR<0.05))


effect_size_phen<-cbind(effect_size_phen,results_ftrs[,1:4])
pval_phen<-cbind(pval_phen,results_ftrs[,5:8])
fdr_phen<-cbind(fdr_phen,results_ftrs[,9:12])

colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-rownames(coef(mss))[8:11]

lim<-dim(effect_size_phen)[2]

##BMI factors

esz<-array(0,dim=c(dim(data_transformed)[1],4))
pv_sp<-array(0,dim=c(dim(data_transformed)[1],4)) #p value 
kl<-300

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],8))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,7]<-as.numeric(sp[k,1])
    cleaned_data[k,8]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","sp","village")
  ind_temp<-0
  for(j in 1:dim(cleaned_data)[1]){
    na_score<-0
    for(k in 1:dim(cleaned_data)[2]){
      if(is.na(cleaned_data[j,k])==T){
        na_score<-na_score+1
      }
    }
    if(na_score==0){
      ind_temp<-rbind(ind_temp,j)
    }
  }
  cleaned_data2<-cleaned_data[ind_temp,]
  #ind_mls<-match(phen_all_rct_use$respondent_master_id,phen_all_use$respondent_master_id)
  cleaned_data2<-cbind(cleaned_data2,"health_dangerous"=factor(phen_all_use[ind_temp,kl])) #Phenotype)
  cleaned_data2 <- within(cleaned_data2, health_dangerous <- relevel(health_dangerous, ref = "Normal"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+health_dangerous+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[8:11,5]
  
  effect.size<-coef(mss)[8:11,1]
  pval2<-coef(mss)[8:11,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size1","Effect_size2","Effect_size3","Effect_size4","Pvalue1","Pvalue2","Pvalue3","Pvalue4")
results_ftrs$FDR1 <- p.adjust(results_ftrs$Pvalue1,method = "BH")
results_ftrs$FDR2 <- p.adjust(results_ftrs$Pvalue2,method = "BH")
results_ftrs$FDR3 <- p.adjust(results_ftrs$Pvalue3,method = "BH")
results_ftrs$FDR4 <- p.adjust(results_ftrs$Pvalue4,method = "BH")

#length(which(results_ftrs$Pvalue<0.05))
#length(which(results_ftrs$FDR<0.05))


effect_size_phen<-cbind(effect_size_phen,results_ftrs[,1:4])
pval_phen<-cbind(pval_phen,results_ftrs[,5:8])
fdr_phen<-cbind(fdr_phen,results_ftrs[,9:12])

colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-rownames(coef(mss))[8:11]

lim<-dim(effect_size_phen)[2]


##HR

esz<-array(0,dim=c(dim(data_transformed)[1],1))
pv_sp<-array(0,dim=c(dim(data_transformed)[1],1)) #p value 
kl<-301

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],8))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,7]<-as.numeric(sp[k,1])
    cleaned_data[k,8]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","sp","village")
  ind_temp<-0
  for(j in 1:dim(cleaned_data)[1]){
    na_score<-0
    for(k in 1:dim(cleaned_data)[2]){
      if(is.na(cleaned_data[j,k])==T){
        na_score<-na_score+1
      }
    }
    if(na_score==0){
      ind_temp<-rbind(ind_temp,j)
    }
  }
  cleaned_data2<-cleaned_data[ind_temp,]
  #ind_mls<-match(phen_all_rct_use$respondent_master_id,phen_all_use$respondent_master_id)
  cleaned_data2<-cbind(cleaned_data2,"health_dangerous"=factor(phen_all_use[ind_temp,kl])) #Phenotype)
  cleaned_data2 <- within(cleaned_data2, health_dangerous <- relevel(health_dangerous, ref = "Normal"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+health_dangerous+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[8,5]
  
  effect.size<-coef(mss)[8,1]
  pval2<-coef(mss)[8,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size1","Pvalue1")
results_ftrs$FDR1 <- p.adjust(results_ftrs$Pvalue1,method = "BH")

#length(which(results_ftrs$Pvalue<0.05))
#length(which(results_ftrs$FDR<0.05))


effect_size_phen<-cbind(effect_size_phen,results_ftrs[,1])
pval_phen<-cbind(pval_phen,results_ftrs[,2])
fdr_phen<-cbind(fdr_phen,results_ftrs[,3])

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("HR(>100)")

lim<-dim(effect_size_phen)[2]

##O2 sat

esz<-array(0,dim=c(dim(data_transformed)[1],1))
pv_sp<-array(0,dim=c(dim(data_transformed)[1],1)) #p value 
kl<-303

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],8))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,7]<-as.numeric(sp[k,1])
    cleaned_data[k,8]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","sp","village")
  ind_temp<-0
  for(j in 1:dim(cleaned_data)[1]){
    na_score<-0
    for(k in 1:dim(cleaned_data)[2]){
      if(is.na(cleaned_data[j,k])==T){
        na_score<-na_score+1
      }
    }
    if(na_score==0){
      ind_temp<-rbind(ind_temp,j)
    }
  }
  cleaned_data2<-cleaned_data[ind_temp,]
  #ind_mls<-match(phen_all_rct_use$respondent_master_id,phen_all_use$respondent_master_id)
  cleaned_data2<-cbind(cleaned_data2,"health_dangerous"=factor(phen_all_use[ind_temp,kl])) #Phenotype)
  cleaned_data2 <- within(cleaned_data2, health_dangerous <- relevel(health_dangerous, ref = "Normal"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+health_dangerous+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[8,5]
  
  effect.size<-coef(mss)[8,1]
  pval2<-coef(mss)[8,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size1","Pvalue1")
results_ftrs$FDR1 <- p.adjust(results_ftrs$Pvalue1,method = "BH")

#length(which(results_ftrs$Pvalue<0.05))
#length(which(results_ftrs$FDR<0.05))


effect_size_phen<-cbind(effect_size_phen,results_ftrs[,1])
pval_phen<-cbind(pval_phen,results_ftrs[,2])
fdr_phen<-cbind(fdr_phen,results_ftrs[,3])

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("O2 saturation(<97)")

lim<-dim(effect_size_phen)[2]

##Hb total

esz<-array(0,dim=c(dim(data_transformed)[1],2))
pv_sp<-array(0,dim=c(dim(data_transformed)[1],2)) #p value 
kl<-304

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],8))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,7]<-as.numeric(sp[k,1])
    cleaned_data[k,8]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","sp","village")
  ind_temp<-0
  for(j in 1:dim(cleaned_data)[1]){
    na_score<-0
    for(k in 1:dim(cleaned_data)[2]){
      if(is.na(cleaned_data[j,k])==T){
        na_score<-na_score+1
      }
    }
    if(na_score==0){
      ind_temp<-rbind(ind_temp,j)
    }
  }
  cleaned_data2<-cleaned_data[ind_temp,]
  #ind_mls<-match(phen_all_rct_use$respondent_master_id,phen_all_use$respondent_master_id)
  cleaned_data2<-cbind(cleaned_data2,"health_dangerous"=factor(phen_all_use[ind_temp,kl])) #Phenotype)
  cleaned_data2 <- within(cleaned_data2, health_dangerous <- relevel(health_dangerous, ref = "Normal"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+health_dangerous+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[8:9,5]
  
  effect.size<-coef(mss)[8:9,1]
  pval2<-coef(mss)[8:9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size1","Effect_size2","Pvalue1","Pvalue2")
results_ftrs$FDR1 <- p.adjust(results_ftrs$Pvalue1,method = "BH")
results_ftrs$FDR2 <- p.adjust(results_ftrs$Pvalue2,method = "BH")

#length(which(results_ftrs$Pvalue<0.05))
#length(which(results_ftrs$FDR<0.05))


effect_size_phen<-cbind(effect_size_phen,results_ftrs[,1:2])
pval_phen<-cbind(pval_phen,results_ftrs[,3:4])
fdr_phen<-cbind(fdr_phen,results_ftrs[,5:6])

colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-rownames(coef(mss))[8:9]



write.csv(effect_size_phen,'effect_size_phen_8_lmer_pwy.csv')
write.csv(pval_phen,'pval_phen_8_lmer_pwy.csv')
write.csv(fdr_phen,'fdr_phen_8_lmer_pwy.csv')




