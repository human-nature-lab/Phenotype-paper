
# ==================================
# By: Shivkumar Vishnempet Shridhar, Christakis and Brito group, HNL, Yale (2023)
# Honduras microbiome project, script for running associations between species and phenotypes
# ==================================



library(lmerTest)

st<-read.csv('st_123.csv',row.names = 1)
st<-as.matrix.data.frame(st)
phen_all_use<-read.csv('phen_b4.csv',row.names = 1,header=T)
data_transformed<-read.csv('mb_sp_10_clr.csv',row.names = 1,header=T)
dm<-read.csv('mb_sp.csv',row.names=1,header=T)
phen_all_rct_use<-read.csv('phen_rct_b4.csv',row.names = 1,header=T)
mb_samp_sp<-data_transformed
data_transformed<-as.matrix.data.frame(data_transformed)


##All physiological variables

kl<-32

for(kl in c(32,35,36)){
  
  esz<-array(0,dim=c(dim(dm)[1],1))
  pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
  
  for(i in 1:dim(data_transformed)[1]){#dim(data_transformed)[1]
    
    sp<-as.data.frame(data_transformed[i,])
    #cleaned_data<-cbind(cleaned_data,sp)
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num -- add
      cleaned_data[k,7]<-as.numeric(st[k])
      cleaned_data[k,8]<-as.numeric(as.character(phen_all_use[k,kl]))
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,9]<-as.numeric(sp[k,1])
      cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","a1c","sp","village")
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
    #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
    #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
    
    
    #s1 = lm(
    #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
    #  data = cleaned_data2
    #)
    
    cleaned_data2$village<-factor(cleaned_data2$village)
    #if(length(which(cleaned_data2$sp>3))){
    
    ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+a1c+(1|village),
              data = cleaned_data2)
    
    mss<-summary(ms1)
    
    coef(mss)[9,5]
    
    effect.size<-coef(mss)[9,1]
    pval2<-coef(mss)[9,5]#Pvalue of phenotype
    
    esz[i,]<-effect.size
    pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    
  }
  #}
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  
  if(kl==32){
    effect_size_phen<-results_ftrs$Effect_size
    pval_phen<-results_ftrs$Pvalue

  }else{
    effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
    pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

  }
  
}

colnames(effect_size_phen)<-c("A1c","Systolic","Diastolic")

length(which(fdr_phen[,3]<0.05))




#################################################################################################################
#Change all -- add sampling date num - column - 371

#MAP

phen_map<-as.numeric(as.character(phen_all_use$mb_d0400))+(as.numeric(as.character(phen_all_use$mb_d0300))-as.numeric(as.character(phen_all_use$mb_d0400)))/3

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    cleaned_data[k,8]<-as.numeric(as.character(phen_map[k]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","map","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+map+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")

effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[4]<-c("MAP")


##BMI


esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,8]<-as.numeric(sp[k,1])
    cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[6,5]
  
  effect.size<-coef(mss)[6,1]
  pval2<-coef(mss)[6,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")

effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[5]<-c("BMI")

##Rest physiological

for(kl in c(71:74,37)){
  
  esz<-array(0,dim=c(dim(dm)[1],1))
  pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
  
  for(i in 1:dim(data_transformed)[1]){
    
    sp<-as.data.frame(data_transformed[i,])
    #cleaned_data<-cbind(cleaned_data,sp)
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
      cleaned_data[k,7]<-as.numeric(st[k])
      cleaned_data[k,8]<-as.numeric(as.character(phen_all_use[k,kl]))
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,9]<-as.numeric(sp[k,1])
      cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","a1c","sp","village")
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
    #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
    #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
    
    
    #s1 = lm(
    #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
    #  data = cleaned_data2
    #)
    
    cleaned_data2$village<-factor(cleaned_data2$village)
    
    
    ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+a1c+(1|village),
              data = cleaned_data2)
    
    mss<-summary(ms1)
    
    coef(mss)[9,5]
    
    effect.size<-coef(mss)[9,1]
    pval2<-coef(mss)[9,5]#Pvalue of phenotype
    
    esz[i,]<-effect.size
    pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    
  }
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  
  effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
  pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  
}

colnames(effect_size_phen)[6:10]<-c("Heart rate","Perfusion index","Pulse","O2 sat","Hb total")


##Cough for 1 month

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    cleaned_data[k,8]<-ifelse(phen_all_use[k,16]=="No",0,1)
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","cough","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+cough+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")

effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Cough(1m)")

##overall health


esz<-array(0,dim=c(dim(dm)[1],4))
pv_sp<-array(0,dim=c(dim(dm)[1],4)) #p value 

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,15]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,8]<-as.numeric(sp[k,1])
    cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","sp","village")
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
  cleaned_data2<-cbind(cleaned_data2,"health"=factor(phen_all_use[ind_temp,15])) #Phenotype)
  cleaned_data2 <- within(cleaned_data2, health <- relevel(health, ref = "Good"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+health+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  effect.size<-coef(mss)[9:12,1]
  pval2<-coef(mss)[9:12,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size1","Effect_size2","Effect_size3","Effect_size4","Pvalue1","Pvalue2","Pvalue3","Pvalue4")


effect_size_phen<-cbind(effect_size_phen,results_ftrs[,1:4])
pval_phen<-cbind(pval_phen,results_ftrs[,5:8])

colnames(effect_size_phen)[12:15]<-c("Excellent","Fair","Poor","Very good") 


for(kl in 20:30){
  
  esz<-array(0,dim=c(dim(dm)[1],1))
  pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
  
  for(i in 1:dim(data_transformed)[1]){
    
    sp<-as.data.frame(data_transformed[i,])
    #cleaned_data<-cbind(cleaned_data,sp)
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
      cleaned_data[k,7]<-as.numeric(st[k])
      cleaned_data[k,8]<-ifelse(is.na(phen_all_use[k,kl])==F,1,0)
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,9]<-as.numeric(sp[k,1])
      cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","a1c","sp","village")
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
    #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
    #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
    
    
    #s1 = lm(
    #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
    #  data = cleaned_data2
    #)
    
    cleaned_data2$village<-factor(cleaned_data2$village)
    
    
    ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+a1c+(1|village),
              data = cleaned_data2)
    
    mss<-summary(ms1)
    
    coef(mss)[9,5]
    
    effect.size<-coef(mss)[9,1]
    pval2<-coef(mss)[9,5]#Pvalue of phenotype
    
    esz[i,]<-effect.size
    pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    
  }
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  
    effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
    pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

  
}

colnames(effect_size_phen)<-c("Diabetes","Allergies","Cystic fibrosis","Heart disease","Endocrine illness","Renal failure","Asthma","Stomach illness","Intestinal illness","Arthritis","MS")



##Medication

for(kl in 136:139){
  
  esz<-array(0,dim=c(dim(dm)[1],1))
  pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
  
  for(i in 1:dim(data_transformed)[1]){
    
    sp<-as.data.frame(data_transformed[i,])
    #cleaned_data<-cbind(cleaned_data,sp)
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
      cleaned_data[k,7]<-as.numeric(st[k])
      cleaned_data[k,8]<-ifelse(is.na(phen_all_use[k,kl])==F,1,0)
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,9]<-as.numeric(sp[k,1])
      cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","drug","sp","village")
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
    #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
    #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
    
    
    #s1 = lm(
    #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
    #  data = cleaned_data2
    #)
    
    cleaned_data2$village<-factor(cleaned_data2$village)
    
    
    ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+drug+(1|village),
              data = cleaned_data2)
    
    mss<-summary(ms1)
    
    coef(mss)[9,5]
    
    effect.size<-coef(mss)[9,1]
    pval2<-coef(mss)[9,5]#Pvalue of phenotype
    
    esz[i,]<-effect.size
    pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    
  }
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  
  
  effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
  pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  
  
}


#medicaiton continued

for(kl in 140:142){
  
  esz<-array(0,dim=c(dim(dm)[1],1))
  pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
  
  for(i in 1:dim(data_transformed)[1]){
    
    sp<-as.data.frame(data_transformed[i,])
    #cleaned_data<-cbind(cleaned_data,sp)
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
      cleaned_data[k,7]<-as.numeric(st[k])
      cleaned_data[k,8]<-ifelse(is.na(phen_all_use[k,kl])==F,1,0)
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,9]<-as.numeric(sp[k,1])
      cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","drug","sp","village")
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
    #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
    #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
    
    
    #s1 = lm(
    #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
    #  data = cleaned_data2
    #)
    
    cleaned_data2$village<-factor(cleaned_data2$village)
    
    
    ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+drug+(1|village),
              data = cleaned_data2)
    
    mss<-summary(ms1)
    
    coef(mss)[9,5]
    
    effect.size<-coef(mss)[9,1]
    pval2<-coef(mss)[9,5]#Pvalue of phenotype
    
    esz[i,]<-effect.size
    pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    
  }
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  
  
    effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
    pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  
  
}

colnames(effect_size_phen)<-c("Anti-fungal","Vitamins","Anti-hypertensive")

print("mental_health")
####Mental health

##Personalty types

for(kl in 143:152){
  for(j in 1:dim(phen_all_use)[1]){
    if(is.na(phen_all_use[j,kl])==F){
      if(phen_all_use[j,kl]=="Dont_Know"){
        phen_all_use[j,kl]<-NA
      }
    }
  }
}

for(kl in 143:152){
  
  esz<-array(0,dim=c(dim(dm)[1],1))
  pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
  
  for(i in 1:dim(data_transformed)[1]){
    
    sp<-as.data.frame(data_transformed[i,])
    #cleaned_data<-cbind(cleaned_data,sp)
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
      cleaned_data[k,7]<-as.numeric(st[k])
      #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
      if(is.na(phen_all_use[k,kl])==F){
        if(phen_all_use[k,kl]=="Disagree strongly"){
          cleaned_data[k,8]<-1
        }else if(phen_all_use[k,kl]=="Disagree a little"){
          cleaned_data[k,8]<-2
        }else if(phen_all_use[k,kl]=="Neither agree nor disagree"){
          cleaned_data[k,8]<-3
        }else if(phen_all_use[k,kl]=="Agree a little"){
          cleaned_data[k,8]<-4
        }else if(phen_all_use[k,kl]=="Agree strongly"){
          cleaned_data[k,8]<-5
        }
      }
      if(is.na(phen_all_use[k,kl])==T){
        cleaned_data[k,8]<-NA
      }
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,9]<-as.numeric(sp[k,1])
      cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","a1c","sp","village")
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
    #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
    #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
    
    
    #s1 = lm(
    #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
    #  data = cleaned_data2
    #)
    
    cleaned_data2$village<-factor(cleaned_data2$village)
    
    
    ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+a1c+(1|village),
              data = cleaned_data2)
    
    mss<-summary(ms1)
    
    coef(mss)[9,5]
    
    effect.size<-coef(mss)[9,1]
    pval2<-coef(mss)[9,5]#Pvalue of phenotype
    
    esz[i,]<-effect.size
    pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    
  }
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  
  
  effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
  pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  
}

colnames(effect_size_phen)[4:dim(effect_size_phen)[2]]<-c("Reserved","Trusting","Lazy","Relaxed","Not creative","Outgoing","Fault others","Thorough job","Nervous","Openess")

lim<-dim(effect_size_phen)[2]

##personality types discrete

print("personality_discrete")

for(i in 1:dim(phen_all_use)[1]){
  for(j in 143:152){
    if(is.na(phen_all_use[i,j])==F){
      if(phen_all_use[i,j]=="Dont_Know"){
        phen_all_use[i,j]<-NA
      }
    }
  }
}

for(kl in 143:152){
  
  
  esz<-array(0,dim=c(dim(dm)[1],4))
  pv_sp<-array(0,dim=c(dim(dm)[1],4)) #p value 
  
  for(i in 1:dim(data_transformed)[1]){
    
    sp<-as.data.frame(data_transformed[i,])
    #cleaned_data<-cbind(cleaned_data,sp)
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
      cleaned_data[k,7]<-as.numeric(st[k])
      cleaned_data[k,8]<-as.numeric(sp[k,1])
      cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","sp","village")
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
    cleaned_data2<-cbind(cleaned_data2,"person_cat"=factor(phen_all_use[ind_temp,kl])) #Phenotype)
    cleaned_data2 <- within(cleaned_data2, person_cat <- relevel(person_cat, ref = "Neither agree nor disagree"))
    
    
    #s1 = lm(
    #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
    #  data = cleaned_data2
    #)
    
    cleaned_data2$village<-factor(cleaned_data2$village)
    
    
    ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+person_cat+(1|village),
              data = cleaned_data2)
    
    mss<-summary(ms1)
    
    coef(mss)[9,5]
    
    effect.size<-coef(mss)[9:12,1]
    pval2<-coef(mss)[9:12,5]#Pvalue of phenotype
    
    esz[i,]<-effect.size
    pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    
  }
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size1","Effect_size2","Effect_size3","Effect_size4","Pvalue1","Pvalue2","Pvalue3","Pvalue4")
  
  
  effect_size_phen<-cbind(effect_size_phen,results_ftrs[,1:4])
  pval_phen<-cbind(pval_phen,results_ftrs[,5:8])
  colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-rownames(coef(mss))[8:11]
  lim<-dim(effect_size_phen)[2]
  
}

lim<-dim(effect_size_phen)[2]

print("cognitive")

##Cognitive score

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    cleaned_data[k,8]<-as.numeric(as.character(phen_all_use[k,153]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","sampling_date","a1c","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+sampling_date+a1c+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")

effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Cognitive score")

lim<-dim(effect_size_phen)[2]

print("cognitive_status")

## Cognitive status

esz<-array(0,dim=c(dim(dm)[1],2))
pv_sp<-array(0,dim=c(dim(dm)[1],2)) #p value 

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    cleaned_data[k,8]<-as.numeric(sp[k,1])
    cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","sp","village")
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
  cleaned_data2<-cbind(cleaned_data2,"cg_status"=factor(phen_all_use[ind_temp,154])) #Phenotype)
  cleaned_data2 <- within(cleaned_data2, cg_status <- relevel(cg_status, ref = "none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+cg_status+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9:10,1]
  pval2<-coef(mss)[9:10,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size1","Effect_size2","Pvalue1","Pvalue2")


effect_size_phen<-cbind(effect_size_phen,results_ftrs[,1:2])
pval_phen<-cbind(pval_phen,results_ftrs[,3:4])
colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-rownames(coef(mss))[9:10]


#GAD7

for(kl in c(183:189,192:200)){
  
  esz<-array(0,dim=c(dim(dm)[1],1))
  pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
  
  for(i in 1:dim(data_transformed)[1]){
    
    sp<-as.data.frame(data_transformed[i,])
    #cleaned_data<-cbind(cleaned_data,sp)
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
      cleaned_data[k,7]<-as.numeric(st[k])
      #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
      if(is.na(phen_all_use[k,kl])==F){
        if(phen_all_use[k,kl]=="Not at all"){
          cleaned_data[k,8]<-0
        }else if(phen_all_use[k,kl]=="Several days"){
          cleaned_data[k,8]<-1
        }else if(phen_all_use[k,kl]=="Over half the days"){
          cleaned_data[k,8]<-2
        }else if(phen_all_use[k,kl]=="Nearly every day"){
          cleaned_data[k,8]<-3
        }
      }
      if(is.na(phen_all_use[k,kl])==T){
        cleaned_data[k,8]<-NA
      }
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,9]<-as.numeric(sp[k,1])
      cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","sampling_date","a1c","sp","village")
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
    #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
    #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
    
    
    #s1 = lm(
    #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
    #  data = cleaned_data2
    #)
    
    cleaned_data2$village<-factor(cleaned_data2$village)
    
    
    ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+a1c+(1|village),
              data = cleaned_data2)
    
    mss<-summary(ms1)
    
    coef(mss)[9,5]
    
    effect.size<-coef(mss)[9,1]
    pval2<-coef(mss)[9,5]#Pvalue of phenotype
    
    esz[i,]<-effect.size
    pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    
  }
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  

    effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
    pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

  
}

colnames(effect_size_phen)<-c("GAD7-Q1","GAD7-Q2","GAD7-Q3","GAD7-Q4","GAD7-Q5","GAD7-Q6","GAD7-Q7","PHQ9-Q1","PHQ9-Q2","PHQ9-Q3","PHQ9-Q4","PHQ9-Q5","PHQ9-Q6","PHQ9-Q7","PHQ9-Q8","PHQ9-Q9")
lim<-dim(effect_size_phen)[2]

#factors

for(kl in c(183:189,192:200)){
  
  esz<-array(0,dim=c(dim(dm)[1],3))
  pv_sp<-array(0,dim=c(dim(dm)[1],3)) #p value 
  
  for(i in 1:dim(data_transformed)[1]){
    
    sp<-as.data.frame(data_transformed[i,])
    #cleaned_data<-cbind(cleaned_data,sp)
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
      cleaned_data[k,7]<-as.numeric(st[k])
      #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,8]<-as.numeric(sp[k,1])
      cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","sp","village")
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
    cleaned_data2<-cbind(cleaned_data2,"ment_cat"=factor(phen_all_use[ind_temp,kl])) #Phenotype)
    cleaned_data2 <- within(cleaned_data2, ment_cat <- relevel(ment_cat, ref = "Not at all"))
    
    
    #s1 = lm(
    #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
    #  data = cleaned_data2
    #)
    
    cleaned_data2$village<-factor(cleaned_data2$village)
    
    
    ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+ment_cat+(1|village),
              data = cleaned_data2)
    
    mss<-summary(ms1)
    
    coef(mss)[9:11,5]
    
    effect.size<-coef(mss)[9:11,1]
    pval2<-coef(mss)[9:11,5]#Pvalue of phenotype
    
    esz[i,]<-effect.size
    pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    
  }
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size1","Effect_size2","Effect_size3","Pvalue1","Pvalue2","Pvalue3")
  
  
  effect_size_phen<-cbind(effect_size_phen,results_ftrs[,1:3])
  pval_phen<-cbind(pval_phen,results_ftrs[,4:6])
  colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-rownames(coef(mss))[8:10]
  lim<-dim(effect_size_phen)[2]
  
}



##Gad7 score

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    cleaned_data[k,8]<-as.numeric(as.character(phen_all_use[k,190]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","a1c","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+a1c+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")

effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("GAD7 score")

lim<-dim(effect_size_phen)[2]

##GAD7 status

esz<-array(0,dim=c(dim(dm)[1],3))
pv_sp<-array(0,dim=c(dim(dm)[1],3)) #p value 

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    cleaned_data[k,8]<-as.numeric(sp[k,1])
    cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","sp","village")
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
  cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+GAD7_cat+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9:11,1]
  pval2<-coef(mss)[9:11,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size_mild","Effect_size_moderate","Effect_size_severe","Pvalue_mild","Pvalue_moderate","Pvalue_severe")
length(which(results_ftrs$Pvalue_mild<0.05))
length(which(results_ftrs$Pvalue_moderate<0.05))
length(which(results_ftrs$Pvalue_severe<0.05))

effect_size_phen<-cbind(effect_size_phen,results_ftrs[,1:3])
pval_phen<-cbind(pval_phen,results_ftrs[,4:6])
colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-rownames(coef(mss))[9:11]
#lim<-dim(effect_size_phen)[2]


#colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-c("GAD7(mild)","GAD7(moderate)","GAD7(severe)")




##PHQ9 score

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    cleaned_data[k,8]<-as.numeric(as.character(phen_all_use[k,201]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","a1c","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+a1c+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")

effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("PHQ9 score")

lim<-dim(effect_size_phen)[2]

##PHQ9 status

esz<-array(0,dim=c(dim(dm)[1],3))
pv_sp<-array(0,dim=c(dim(dm)[1],3)) #p value 

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    cleaned_data[k,8]<-as.numeric(sp[k,1])
    cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","sp","village")
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
  cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,202])) #Phenotype)
  cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+GAD7_cat+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9:11,1]
  pval2<-coef(mss)[9:11,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size_mild","Effect_size_moderate","Effect_size_severe","Pvalue_mild","Pvalue_moderate","Pvalue_severe")
length(which(results_ftrs$Pvalue_mild<0.05))
length(which(results_ftrs$Pvalue_moderate<0.05))
length(which(results_ftrs$Pvalue_severe<0.05))

effect_size_phen<-cbind(effect_size_phen,results_ftrs[,1:3])
pval_phen<-cbind(pval_phen,results_ftrs[,4:6])
colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-rownames(coef(mss))[9:11]

for(kl in 155:180){
  
  esz<-array(0,dim=c(dim(dm)[1],1))
  pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
  
  for(i in 1:dim(data_transformed)[1]){
    
    sp<-as.data.frame(data_transformed[i,])
    #cleaned_data<-cbind(cleaned_data,sp)
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
      cleaned_data[k,7]<-as.numeric(st[k])
      cleaned_data[k,8]<-ifelse(is.na(phen_all_use[k,kl])==F,1,0)
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,9]<-as.numeric(sp[k,1])
      cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","animals","sp","village")
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
    #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
    #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
    
    
    #s1 = lm(
    #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
    #  data = cleaned_data2
    #)
    
    cleaned_data2$village<-factor(cleaned_data2$village)
    
    
    ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+animals+(1|village),
              data = cleaned_data2)
    
    mss<-summary(ms1)
    
    coef(mss)[9,5]
    
    effect.size<-coef(mss)[9,1]
    pval2<-coef(mss)[9,5]#Pvalue of phenotype
    
    esz[i,]<-effect.size
    pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    
  }
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  

    effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
    pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  
}

colnames(effect_size_phen)<-c("Cat","Dog","Parakeet","Rabbit","Horse","Mice","None pet","Cow","Goat","Pig","Chicken","Duck","Turkey","Sheep","Geese","None farm","Bat","Lizard","Monkey","Snake","Bird","Possum","Rat","Squirrel","Other wild","None wild")

lim<-dim(effect_size_phen)[2]

##Food

for(kl in 75:93){
  
  esz<-array(0,dim=c(dim(dm)[1],1))
  pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
  
  for(i in 1:dim(data_transformed)[1]){
    
    sp<-as.data.frame(data_transformed[i,])
    #cleaned_data<-cbind(cleaned_data,sp)
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
      cleaned_data[k,7]<-as.numeric(st[k])
      #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
      if(is.na(phen_all_use[k,kl])==F){
        if(phen_all_use[k,kl]=="Never/rarely"){#Change to 1/50,1/10,4/7,1
          cleaned_data[k,8]<-1/50
        }else if((phen_all_use[k,kl]=="A few days per month")|(phen_all_use[k,kl]=="A few times per month")){
          cleaned_data[k,8]<-1/10
        }else if(phen_all_use[k,kl]=="A few days per week"){
          cleaned_data[k,8]<-4/7
        }else if(phen_all_use[k,kl]=="Every day"){
          cleaned_data[k,8]<-1
        }
      }
      if(is.na(phen_all_use[k,kl])==T){
        cleaned_data[k,8]<-NA
      }
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,9]<-as.numeric(sp[k,1])
      cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","food","sp","village")
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
    #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
    #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
    
    
    #s1 = lm(
    #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
    #  data = cleaned_data2
    #)
    
    cleaned_data2$village<-factor(cleaned_data2$village)
    
    
    ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+food+(1|village),
              data = cleaned_data2)
    
    mss<-summary(ms1)
    
    coef(mss)[9,5]
    
    effect.size<-coef(mss)[9,1]
    pval2<-coef(mss)[9,5]#Pvalue of phenotype
    
    esz[i,]<-effect.size
    pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    
  }
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  
  
  effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
  pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  
}
colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-c("Beans","Tortillas","Rice","Bread","Milk","Yogurt","Cream/butter","Cheese","Eggs","Vegetables","Fruits","Natural juice","Chicken","Beef/Pork","Ham/sausages/hotdog","Fish","Soda","Fruit juice","Chips")

##Diet diversity score

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
kl<-382

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    if(is.na(phen_all_use[k,kl])==F){
      cleaned_data[k,8]<-as.numeric(as.character(phen_all_use[k,kl]))
    }
    if(is.na(phen_all_use[k,kl])==T){
      cleaned_data[k,8]<-NA
    }
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","food","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+food+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-"DDS"

##GAD7 score --- 0 or 1

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    if(is.na(phen_all_use[k,190])==F){
      if(phen_all_use[k,190]<=8){
        cleaned_data[k,8]<-0
      }else{
        cleaned_data[k,8]<-1
      }
    }
    if(is.na(phen_all_use[k,190])==T){
      cleaned_data[k,8]<-NA
    }
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","a1c","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+a1c+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

print("GAD7-dict")
# colnames(effect_size_phen)<-c("GAD7-dict")

##PHQ9 score --- 0 or 1

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    if(is.na(phen_all_use[k,201])==F){
      if(phen_all_use[k,201]<=8){
        cleaned_data[k,8]<-0
      }else{
        cleaned_data[k,8]<-1
      }
    }
    if(is.na(phen_all_use[k,201])==T){
      cleaned_data[k,7]<-NA
    }
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","a1c","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+a1c+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")

effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)<-c("GAD7-dict","PHQ9-dict")

print("travel")  
##Travel

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    if(is.na(phen_all_use[k,181])==F){
      if(phen_all_use[k,181]=="Rarely/Never"){
        cleaned_data[k,8]<-1/100
      }else if(phen_all_use[k,181]=="At least once a month"){
        cleaned_data[k,8]<-1/30
      }else if(phen_all_use[k,181]=="At least once a week"){
        cleaned_data[k,8]<-1/7
      }else if(phen_all_use[k,181]=="Every day"){
        cleaned_data[k,8]<-1
      }
    }
    if(is.na(phen_all_use[k,181])==T){
      cleaned_data[k,8]<-NA
    }
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","travel","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+travel+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Travel")

##Monthly expenditure

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    cleaned_data[k,8]<-as.numeric(as.character(phen_all_use[k,182]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","a1c","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+a1c+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")

effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Monthly expenditure")

print("alcohol")  
##Alcohol 

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
kl<-7

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_rct_use)[1],10))
  for(k in 1:dim(phen_all_rct_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    if(is.na(phen_all_rct_use[k,kl])==F){
      if(phen_all_rct_use[k,kl]=="Never"){#change to 0,1/30,3/30,3/7,1
        cleaned_data[k,8]<-0
      }else if(phen_all_rct_use[k,kl]=="Monthly or less"){
        cleaned_data[k,8]<-1/30
      }else if(phen_all_rct_use[k,kl]=="Two to four times a month"){
        cleaned_data[k,8]<-4/30
      }else if(phen_all_rct_use[k,kl]=="Two to three times a week"){
        cleaned_data[k,8]<-3/7
      }else if(phen_all_rct_use[k,kl]=="Four or more times a week"){
        cleaned_data[k,8]<-6/7
      }
    }
    if(is.na(phen_all_rct_use[k,kl])==T){
      cleaned_data[k,8]<-NA
    }
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","travel","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+travel+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Alcohol freq") 


## Alcohol daily

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
kl<-8

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_rct_use)[1],10))
  for(k in 1:dim(phen_all_rct_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    if(is.na(phen_all_rct_use[k,kl])==F){
      if(phen_all_rct_use[k,kl]=="1 or 2"){
        cleaned_data[k,8]<-1.5
      }else if(phen_all_rct_use[k,kl]=="3 or 4"){
        cleaned_data[k,8]<-3.5
      }else if(phen_all_rct_use[k,kl]=="5 or 6"){
        cleaned_data[k,8]<-5.5
      }else if(phen_all_rct_use[k,kl]=="7 to 9"){
        cleaned_data[k,8]<-8
      }else if(phen_all_rct_use[k,kl]=="10 or more"){
        cleaned_data[k,8]<-12
      }
    }
    if(is.na(phen_all_rct_use[k,kl])==T){
      cleaned_data[k,8]<-0
    }
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","travel","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+travel+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Alcohol daily") 


##Alcohol 6 drinks

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
kl<-9

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_rct_use)[1],10))
  for(k in 1:dim(phen_all_rct_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    if(is.na(phen_all_rct_use[k,kl])==F){
      if(phen_all_rct_use[k,kl]=="Never"){
        cleaned_data[k,8]<-0
      }else if(phen_all_rct_use[k,kl]=="Less than monthly"){
        cleaned_data[k,8]<-1/50
      }else if(phen_all_rct_use[k,kl]=="Monthly"){
        cleaned_data[k,8]<-1/30
      }else if(phen_all_rct_use[k,kl]=="Weekly"){
        cleaned_data[k,8]<-1/7
      }else if(phen_all_rct_use[k,kl]=="Daily or almost daily"){
        cleaned_data[k,8]<-1
      }
    }
    if(is.na(phen_all_rct_use[k,kl])==T){
      cleaned_data[k,8]<-0
    }
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","travel","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+travel+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Alcohol 6 drinks") 

print("cigarette")  

##Cigarette

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
kl<-10

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_rct_use)[1],10))
  for(k in 1:dim(phen_all_rct_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    if(is.na(phen_all_rct_use[k,kl])==F){
      if(phen_all_rct_use[k,kl]=="No"){
        cleaned_data[k,8]<-0
      }
      if(phen_all_rct_use[k,kl]=="Yes"){
        cleaned_data[k,8]<-1
      }
    }
    if(is.na(phen_all_rct_use[k,kl])==T){
      cleaned_data[k,8]<-NA
    }
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","travel","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+travel+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Cigarette")

##Cigarette frequency

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
kl<-11

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_rct_use)[1],10))
  for(k in 1:dim(phen_all_rct_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    if(is.na(phen_all_rct_use[k,kl])==F){
      cleaned_data[k,8]<-as.numeric(phen_all_rct_use[k,kl])
    }
    if(is.na(phen_all_rct_use[k,kl])==T){
      cleaned_data[k,8]<-0
    }
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","travel","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+travel+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Cigarette frequency")

print("altruism")
##Altruism

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    cleaned_data[k,8]<-as.numeric(as.character(phen_all_use[k,208]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","altruism","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+altruism+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")

effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Altruism")




print("risk taking 1")   



##Risk taking

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    cleaned_data[k,8]<-as.numeric(as.character(phen_all_use[k,364]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","a1c","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+a1c+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")

effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Risk taking")

print("risk taking")  


##Education

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
kl<-5

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_rct_use)[1],10))
  for(k in 1:dim(phen_all_rct_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    if(is.na(phen_all_rct_use[k,kl])==F){
      if(phen_all_rct_use[k,kl]=="Have not completed any type of school"){
        cleaned_data[k,8]<-0
      }else if(phen_all_rct_use[k,kl]=="1st grade"){
        cleaned_data[k,8]<-1
      }else if(phen_all_rct_use[k,kl]=="2nd grade"){
        cleaned_data[k,8]<-2
      }else if(phen_all_rct_use[k,kl]=="3rd grade"){
        cleaned_data[k,8]<-3
      }else if(phen_all_rct_use[k,kl]=="4th grade"){
        cleaned_data[k,8]<-4
      }else if(phen_all_rct_use[k,kl]=="5th grade"){
        cleaned_data[k,8]<-5
      }else if(phen_all_rct_use[k,kl]=="6th grade"){
        cleaned_data[k,8]<-6
      }else if(phen_all_rct_use[k,kl]=="Some secondary"){
        cleaned_data[k,8]<-7
      }else if(phen_all_rct_use[k,kl]=="More than secondary"){
        cleaned_data[k,8]<-8
      }
    }
    if(is.na(phen_all_rct_use[k,kl])==T){
      cleaned_data[k,8]<-NA
    }
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","education","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+education+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Education") 

lim<-dim(effect_size_phen)[2]

print("education")  

##Education factors

esz<-array(0,dim=c(dim(dm)[1],9))
pv_sp<-array(0,dim=c(dim(dm)[1],9)) #p value 
kl<-5

#for(im in 1:dim(phen_all_rct_use)[1]){
#  if(is.na(phen_all_rct_use[im,41])==F){
#  if(phen_all_rct_use[im,41]==0){
#    phen_all_rct_use[im,41]<-NA
#  }
#  }
#}

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_rct_use)[1],9))
  for(k in 1:dim(phen_all_rct_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,8]<-as.numeric(sp[k,1])
    cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","sp","village")
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
  cleaned_data2<-cbind(cleaned_data2,"education"=factor(phen_all_rct_use[ind_temp,kl])) #Phenotype)
  cleaned_data2 <- within(cleaned_data2, education <- relevel(education, ref = "Have not completed any type of school"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+education+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9:17,5]
  
  effect.size<-coef(mss)[9:17,1]
  pval2<-coef(mss)[9:17,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size1","Effect_size2","Effect_size3","Effect_size4","Effect_size5","Effect_size6","Effect_size7","Effect_size8","Effect_size9","Pvalue1","Pvalue2","Pvalue3","Pvalue4","Pvalue5","Pvalue6","Pvalue7","Pvalue8","Pvalue9")


effect_size_phen<-cbind(effect_size_phen,results_ftrs[,1:9])
pval_phen<-cbind(pval_phen,results_ftrs[,10:18])

#colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-c("1st grade","2nd grade","3rd grade","4th grade","5th grade","6th grade","Some secondary","More than secondary")
colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-rownames(coef(mss))[9:17]

lim<-dim(effect_size_phen)[2]

print("social")  

##Social network factors
colnames(phen_all_use)[227:230]
colnames(phen_all_rct_use)[14:17]
phen_all_use[,227:230]<-phen_all_rct_use[,14:17]

for(kl in c(203:206,360:363,365:368,228:229)){
  
  esz<-array(0,dim=c(dim(dm)[1],1))
  pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
  
  for(i in 1:dim(data_transformed)[1]){
    
    sp<-as.data.frame(data_transformed[i,])
    #cleaned_data<-cbind(cleaned_data,sp)
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
      cleaned_data[k,7]<-as.numeric(st[k])
      cleaned_data[k,8]<-as.numeric(as.character(phen_all_use[k,kl]))
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,9]<-as.numeric(sp[k,1])
      cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","soc","sp","village")
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
    #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
    #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
    
    
    #s1 = lm(
    #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
    #  data = cleaned_data2
    #)
    
    cleaned_data2$village<-factor(cleaned_data2$village)
    
    
    ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+soc+(1|village),
              data = cleaned_data2)
    
    mss<-summary(ms1)
    
    coef(mss)[9,5]
    
    effect.size<-coef(mss)[9,1]
    pval2<-coef(mss)[9,5]#Pvalue of phenotype
    
    esz[i,]<-effect.size
    pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    
  }
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  
  effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
  pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  
}

colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-colnames(phen_all_use)[c(203:206,360:363,365:368,228:229)]

print("partner") 

for(kl in 230){
  
  esz<-array(0,dim=c(dim(dm)[1],1))
  pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
  
  for(i in 1:dim(data_transformed)[1]){
    
    sp<-as.data.frame(data_transformed[i,])
    #cleaned_data<-cbind(cleaned_data,sp)
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
      cleaned_data[k,7]<-as.numeric(st[k])
      cleaned_data[k,8]<-as.numeric(as.character(phen_all_use[k,kl]))
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,9]<-as.numeric(sp[k,1])
      cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","soc","sp","village")
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
    #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
    #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
    
    
    #s1 = lm(
    #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
    #  data = cleaned_data2
    #)
    
    cleaned_data2$village<-factor(cleaned_data2$village)
    
    if(i!=2261){#Totally wierd error -- DOwndated VTV
      ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+soc+(1|village),
                data = cleaned_data2)
      
      mss<-summary(ms1)
      
      coef(mss)[9,5]
      
      effect.size<-coef(mss)[9,1]
      pval2<-coef(mss)[9,5]#Pvalue of phenotype
      
      esz[i,]<-effect.size
      pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    }
  }
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  
  effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
  pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  
}

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-colnames(phen_all_use)[230]



print("partner last")

#Living with partner(yes or no)

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
kl<-227

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    if(is.na(phen_all_use[k,kl])==F){
      if(phen_all_use[k,kl]=="No"){
        cleaned_data[k,8]<-0
      }
      if(phen_all_use[k,kl]=="Yes"){
        cleaned_data[k,8]<-1
      }
    }
    if(is.na(phen_all_use[k,kl])==T){
      cleaned_data[k,8]<-NA
    }
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","partner","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+partner+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Live with partner")

##Washing hands

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
kl<-12

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_rct_use)[1],10))
  for(k in 1:dim(phen_all_rct_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    if(is.na(phen_all_rct_use[k,kl])==F){
      cleaned_data[k,8]<-1
    }
    if(is.na(phen_all_rct_use[k,kl])==T){
      cleaned_data[k,8]<-0
    }
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","travel","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+travel+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")



effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Washing hands")

#Distance from village center --population weighted

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    cleaned_data[k,8]<-as.numeric(as.character(phen_all_use[k,372]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","a1c","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+a1c+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")

effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Distance_center")


##HH essentials

for(kl in 310:316){
  
  esz<-array(0,dim=c(dim(dm)[1],1))
  pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
  
  for(i in 1:dim(data_transformed)[1]){
    
    sp<-as.data.frame(data_transformed[i,])
    #cleaned_data<-cbind(cleaned_data,sp)
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
      cleaned_data[k,7]<-as.numeric(st[k])
      cleaned_data[k,8]<-ifelse(is.na(phen_all_use[k,kl])==F,1,0)
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,9]<-as.numeric(sp[k,1])
      cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","hh","sp","village")
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
    #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
    #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
    
    
    #s1 = lm(
    #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
    #  data = cleaned_data2
    #)
    
    cleaned_data2$village<-factor(cleaned_data2$village)
    
    
    ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+hh+(1|village),
              data = cleaned_data2)
    
    mss<-summary(ms1)
    
    coef(mss)[9,5]
    
    effect.size<-coef(mss)[9,1]
    pval2<-coef(mss)[9,5]#Pvalue of phenotype
    
    esz[i,]<-effect.size
    pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    
  }
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  
    effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
    pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  
}

colnames(effect_size_phen)<-c("Electricity","Radio","Television","Cell phone","No cellphone","Refrigerator","None")
lim<-dim(effect_size_phen)[2]

##HH remaining

for(kl in c(306:309,317:335,337:341,343:349,351:359)){
  
  esz<-array(0,dim=c(dim(dm)[1],1))
  pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
  
  #if(sum(phen_all_use[,kl])>0){
  for(i in 1:dim(data_transformed)[1]){
    
    sp<-as.data.frame(data_transformed[i,])
    #cleaned_data<-cbind(cleaned_data,sp)
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
      cleaned_data[k,7]<-as.numeric(st[k])
      cleaned_data[k,8]<-as.numeric(as.character(phen_all_use[k,kl]))
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,9]<-as.numeric(sp[k,1])
      cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","hh","sp","village")
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
    #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
    #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
    
    
    #s1 = lm(
    #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
    #  data = cleaned_data2
    #)
    
    cleaned_data2$village<-factor(cleaned_data2$village)
    
    if(sum(cleaned_data2[,7])>0){
      
      ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+hh+(1|village),
                data = cleaned_data2)
      
      mss<-summary(ms1)
      
      coef(mss)[9,5]
      
      effect.size<-coef(mss)[9,1]
      pval2<-coef(mss)[9,5]#Pvalue of phenotype
      
      esz[i,]<-effect.size
      pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    }else{
      esz[i,]<-0
      pv_sp[i,]<-1
    }
    
  }
  #}
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  
  effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
  pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  
  
}

colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-colnames(phen_all_use)[c(306:309,317:335,337:341,343:349,351:359)]

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
    }else if(phen_all_use[i,74]=="Refused"){
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

esz<-array(0,dim=c(dim(dm)[1],3))
pv_sp<-array(0,dim=c(dim(dm)[1],3)) #p value 
kl<-298

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,8]<-as.numeric(sp[k,1])
    cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","sp","village")
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
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+health_dangerous+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9:11,5]
  
  effect.size<-coef(mss)[9:11,1]
  pval2<-coef(mss)[9:11,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size1","Effect_size2","Effect_size3","Pvalue1","Pvalue2","Pvalue3")

#length(which(results_ftrs$Pvalue<0.05))
#length(which(results_ftrs$FDR<0.05))


effect_size_phen<-results_ftrs[,1:3]
pval_phen<-results_ftrs[,4:6]

colnames(effect_size_phen)<-rownames(coef(mss))[9:11]

lim<-dim(effect_size_phen)[2]


##BP factors

esz<-array(0,dim=c(dim(dm)[1],4))
pv_sp<-array(0,dim=c(dim(dm)[1],4)) #p value 
kl<-299

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,8]<-as.numeric(sp[k,1])
    cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","sp","village")
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
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+health_dangerous+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9:12,5]
  
  effect.size<-coef(mss)[9:12,1]
  pval2<-coef(mss)[9:12,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size1","Effect_size2","Effect_size3","Effect_size4","Pvalue1","Pvalue2","Pvalue3","Pvalue4")



effect_size_phen<-cbind(effect_size_phen,results_ftrs[,1:4])
pval_phen<-cbind(pval_phen,results_ftrs[,5:8])

colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-rownames(coef(mss))[9:12]

lim<-dim(effect_size_phen)[2]

##BMI factors

esz<-array(0,dim=c(dim(dm)[1],4))
pv_sp<-array(0,dim=c(dim(dm)[1],4)) #p value 
kl<-300

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,8]<-as.numeric(sp[k,1])
    cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","sp","village")
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
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+health_dangerous+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9:12,5]
  
  effect.size<-coef(mss)[9:12,1]
  pval2<-coef(mss)[9:12,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size1","Effect_size2","Effect_size3","Effect_size4","Pvalue1","Pvalue2","Pvalue3","Pvalue4")


effect_size_phen<-cbind(effect_size_phen,results_ftrs[,1:4])
pval_phen<-cbind(pval_phen,results_ftrs[,5:8])

colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-rownames(coef(mss))[9:12]

lim<-dim(effect_size_phen)[2]


##HR

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
kl<-301

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,8]<-as.numeric(sp[k,1])
    cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","sp","village")
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
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+health_dangerous+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size1","Pvalue1")



effect_size_phen<-cbind(effect_size_phen,results_ftrs[,1])
pval_phen<-cbind(pval_phen,results_ftrs[,2])

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("HR(>100)")

lim<-dim(effect_size_phen)[2]

##O2 sat

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
kl<-303

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,8]<-as.numeric(sp[k,1])
    cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","sp","village")
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
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+health_dangerous+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size1","Pvalue1")


effect_size_phen<-cbind(effect_size_phen,results_ftrs[,1])
pval_phen<-cbind(pval_phen,results_ftrs[,2])

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("O2 saturation(<97)")

lim<-dim(effect_size_phen)[2]

##Hb total

esz<-array(0,dim=c(dim(dm)[1],2))
pv_sp<-array(0,dim=c(dim(dm)[1],2)) #p value 
kl<-304

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,8]<-as.numeric(sp[k,1])
    cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","sp","village")
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
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+health_dangerous+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9:10,5]
  
  effect.size<-coef(mss)[9:10,1]
  pval2<-coef(mss)[9:10,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size1","Effect_size2","Pvalue1","Pvalue2")

#length(which(results_ftrs$Pvalue<0.05))
#length(which(results_ftrs$FDR<0.05))


effect_size_phen<-cbind(effect_size_phen,results_ftrs[,1:2])
pval_phen<-cbind(pval_phen,results_ftrs[,3:4])

colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-rownames(coef(mss))[9:10]

for(kl in c(375:380,383)){
  
  esz<-array(0,dim=c(dim(dm)[1],1))
  pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
  
  for(i in 1:dim(data_transformed)[1]){
    
    sp<-as.data.frame(data_transformed[i,])
    #cleaned_data<-cbind(cleaned_data,sp)
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
      cleaned_data[k,7]<-as.numeric(st[k])
      cleaned_data[k,8]<-as.numeric(as.character(phen_all_use[k,kl]))
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,9]<-as.numeric(sp[k,1])
      cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","sampling_date","hh","sp","village")
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
    #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
    #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
    
    
    #s1 = lm(
    #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
    #  data = cleaned_data2
    #)
    
    cleaned_data2$village<-factor(cleaned_data2$village)
    
    if(sum(cleaned_data2[,7])>0){
      
      ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+hh+(1|village),
                data = cleaned_data2)
      
      mss<-summary(ms1)
      
      coef(mss)[9,5]
      
      effect.size<-coef(mss)[9,1]
      pval2<-coef(mss)[9,5]#Pvalue of phenotype
      
      esz[i,]<-effect.size
      pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    }else{
      esz[i,]<-0
      pv_sp[i,]<-1
    }
    
  }
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  print(colnames(phen_all_use)[kl]) 
  

    effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
    pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  
  
}

colnames(effect_size_phen)<-colnames(phen_all_use)[c(375:380,383)]


##Education
#Primary education

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
kl<-5

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_rct_use)[1],10))
  for(k in 1:dim(phen_all_rct_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    if(is.na(phen_all_rct_use[k,kl])==F){
      if(phen_all_rct_use[k,kl]=="Have not completed any type of school"){
        cleaned_data[k,8]<-0
      }else if(phen_all_rct_use[k,kl]=="1st grade"){
        cleaned_data[k,8]<-1
      }else if(phen_all_rct_use[k,kl]=="2nd grade"){
        cleaned_data[k,8]<-1
      }else if(phen_all_rct_use[k,kl]=="3rd grade"){
        cleaned_data[k,8]<-1
      }else if(phen_all_rct_use[k,kl]=="4th grade"){
        cleaned_data[k,8]<-0
      }else if(phen_all_rct_use[k,kl]=="5th grade"){
        cleaned_data[k,8]<-0
      }else if(phen_all_rct_use[k,kl]=="6th grade"){
        cleaned_data[k,8]<-0
      }else if(phen_all_rct_use[k,kl]=="Some secondary"){
        cleaned_data[k,8]<-0
      }else if(phen_all_rct_use[k,kl]=="More than secondary"){
        cleaned_data[k,8]<-0
      }
    }
    if(is.na(phen_all_rct_use[k,kl])==T){
      cleaned_data[k,8]<-NA
    }
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","education","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+education+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Education(primary)") 

print("education")  

#Middle education

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
kl<-5

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_rct_use)[1],10))
  for(k in 1:dim(phen_all_rct_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    if(is.na(phen_all_rct_use[k,kl])==F){
      if(phen_all_rct_use[k,kl]=="Have not completed any type of school"){
        cleaned_data[k,8]<-0
      }else if(phen_all_rct_use[k,kl]=="1st grade"){
        cleaned_data[k,8]<-0
      }else if(phen_all_rct_use[k,kl]=="2nd grade"){
        cleaned_data[k,8]<-0
      }else if(phen_all_rct_use[k,kl]=="3rd grade"){
        cleaned_data[k,8]<-0
      }else if(phen_all_rct_use[k,kl]=="4th grade"){
        cleaned_data[k,8]<-1
      }else if(phen_all_rct_use[k,kl]=="5th grade"){
        cleaned_data[k,8]<-1
      }else if(phen_all_rct_use[k,kl]=="6th grade"){
        cleaned_data[k,8]<-1
      }else if(phen_all_rct_use[k,kl]=="Some secondary"){
        cleaned_data[k,8]<-0
      }else if(phen_all_rct_use[k,kl]=="More than secondary"){
        cleaned_data[k,8]<-0
      }
    }
    if(is.na(phen_all_rct_use[k,kl])==T){
      cleaned_data[k,8]<-NA
    }
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","education","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+education+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Education(4th-6th grade)") 

print("middle education")  

#Secondary education

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 
kl<-5

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_rct_use)[1],10))
  for(k in 1:dim(phen_all_rct_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    if(is.na(phen_all_rct_use[k,kl])==F){
      if(phen_all_rct_use[k,kl]=="Have not completed any type of school"){
        cleaned_data[k,8]<-0
      }else if(phen_all_rct_use[k,kl]=="1st grade"){
        cleaned_data[k,8]<-0
      }else if(phen_all_rct_use[k,kl]=="2nd grade"){
        cleaned_data[k,8]<-0
      }else if(phen_all_rct_use[k,kl]=="3rd grade"){
        cleaned_data[k,8]<-0
      }else if(phen_all_rct_use[k,kl]=="4th grade"){
        cleaned_data[k,8]<-0
      }else if(phen_all_rct_use[k,kl]=="5th grade"){
        cleaned_data[k,8]<-0
      }else if(phen_all_rct_use[k,kl]=="6th grade"){
        cleaned_data[k,8]<-0
      }else if(phen_all_rct_use[k,kl]=="Some secondary"){
        cleaned_data[k,8]<-1
      }else if(phen_all_rct_use[k,kl]=="More than secondary"){
        cleaned_data[k,8]<-1
      }
    }
    if(is.na(phen_all_rct_use[k,kl])==T){
      cleaned_data[k,8]<-NA
    }
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","education","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+education+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Education(secondary)") 

print("secondary education")  

## Diarrhea

#colnames(phen_all_use)

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    cleaned_data[k,8]<-ifelse(phen_all_use[k,381]=="No",0,1)
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","cough","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+cough+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")

effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Diarrhea")


print("Diarrhea")

##Bristol stool

esz<-array(0,dim=c(dim(dm)[1],1))
pv_sp<-array(0,dim=c(dim(dm)[1],1)) #p value 

for(i in 1:dim(data_transformed)[1]){
  
  sp<-as.data.frame(data_transformed[i,])
  #cleaned_data<-cbind(cleaned_data,sp)
  cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],10))
  for(k in 1:dim(phen_all_use)[1]){
    cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
    cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
    cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
    cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
    cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
    cleaned_data[k,6]<-as.numeric(phen_all_use[k,371])#Sampling date num
    cleaned_data[k,7]<-as.numeric(st[k])
    cleaned_data[k,8]<-ifelse(phen_all_use[k,381]=="No",0,1)
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,9]<-as.numeric(sp[k,1])
    cleaned_data[k,10]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","cough","sp","village")
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
  #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+cough+bristol+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[9,5]
  effect.size<-coef(mss)[9,1]
  pval2<-coef(mss)[9,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")

effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Bristol")


print("Bristol")

write.csv(effect_size_phen,'effect_size_phen_all.csv')
write.csv(pval_phen,'pval_phen__all.csv')

