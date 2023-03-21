library(lmerTest)

st<-read.csv('bristol.csv',row.names = 1)
st<-as.matrix.data.frame(st)
#phen_all_use<-read.csv('phen_all_use.csv',row.names = 1,header=T)
phen_all_use<-read.csv('phen_all_use6_all.csv',row.names = 1,header=T)
#data_transformed<-read.csv('data_transformed.csv',row.names = 1,header=T)
data_transformed<-read.csv('pwy_all_use.csv',row.names = 1,header=T)
phen_all_rct_use<-read.csv('phen_all_rct_use22.csv',row.names = 1,header=T)
data_transformed<-as.matrix.data.frame(data_transformed)

for(kl in 375:380){
  
  esz<-array(0,dim=c(dim(data_transformed)[1],1))
  pv_sp<-array(0,dim=c(dim(data_transformed)[1],1)) #p value 
  
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
      cleaned_data[k,6]<-as.numeric(st[k])
      cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,8]<-as.numeric(sp[k,1])
      cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","hh","sp","village")
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
      
      ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+hh+(1|village),
                data = cleaned_data2)
      
      mss<-summary(ms1)
      
      coef(mss)[8,5]
      
      effect.size<-coef(mss)[8,1]
      pval2<-coef(mss)[8,5]#Pvalue of phenotype
      
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
  results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
  length(which(results_ftrs$Pvalue<0.05))
  length(which(results_ftrs$FDR<0.05))
  print(colnames(phen_all_use)[kl]) 
  
  if(kl==375){
    effect_size_phen<-results_ftrs$Effect_size
    pval_phen<-results_ftrs$Pvalue
    fdr_phen<-results_ftrs$FDR
  }else{
    effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
    pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
    fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)
    
  }
  
  
}

colnames(effect_size_phen)<-colnames(phen_all_use)[375:380]


##Education
#Primary education

esz<-array(0,dim=c(dim(data_transformed)[1],1))
pv_sp<-array(0,dim=c(dim(data_transformed)[1],1)) #p value 
kl<-5

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
    cleaned_data[k,6]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    if(is.na(phen_all_rct_use[k,kl])==F){
      if(phen_all_rct_use[k,kl]=="Have not completed any type of school"){
        cleaned_data[k,7]<-0
      }else if(phen_all_rct_use[k,kl]=="1st grade"){
        cleaned_data[k,7]<-1
      }else if(phen_all_rct_use[k,kl]=="2nd grade"){
        cleaned_data[k,7]<-1
      }else if(phen_all_rct_use[k,kl]=="3rd grade"){
        cleaned_data[k,7]<-1
      }else if(phen_all_rct_use[k,kl]=="4th grade"){
        cleaned_data[k,7]<-0
      }else if(phen_all_rct_use[k,kl]=="5th grade"){
        cleaned_data[k,7]<-0
      }else if(phen_all_rct_use[k,kl]=="6th grade"){
        cleaned_data[k,7]<-0
      }else if(phen_all_rct_use[k,kl]=="Some secondary"){
        cleaned_data[k,7]<-0
      }else if(phen_all_rct_use[k,kl]=="More than secondary"){
        cleaned_data[k,7]<-0
      }
    }
    if(is.na(phen_all_rct_use[k,kl])==T){
      cleaned_data[k,7]<-NA
    }
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,8]<-as.numeric(sp[k,1])
    cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","education","sp","village")
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
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+education+(1|village),
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
colnames(results_ftrs)<-c("Effect_size","Pvalue")
results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
length(which(results_ftrs$Pvalue<0.05))
length(which(results_ftrs$FDR<0.05))


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Education(primary)") 

print("education")  

#Middle education

esz<-array(0,dim=c(dim(data_transformed)[1],1))
pv_sp<-array(0,dim=c(dim(data_transformed)[1],1)) #p value 
kl<-5

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
    cleaned_data[k,6]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    if(is.na(phen_all_rct_use[k,kl])==F){
      if(phen_all_rct_use[k,kl]=="Have not completed any type of school"){
        cleaned_data[k,7]<-0
      }else if(phen_all_rct_use[k,kl]=="1st grade"){
        cleaned_data[k,7]<-0
      }else if(phen_all_rct_use[k,kl]=="2nd grade"){
        cleaned_data[k,7]<-0
      }else if(phen_all_rct_use[k,kl]=="3rd grade"){
        cleaned_data[k,7]<-0
      }else if(phen_all_rct_use[k,kl]=="4th grade"){
        cleaned_data[k,7]<-1
      }else if(phen_all_rct_use[k,kl]=="5th grade"){
        cleaned_data[k,7]<-1
      }else if(phen_all_rct_use[k,kl]=="6th grade"){
        cleaned_data[k,7]<-1
      }else if(phen_all_rct_use[k,kl]=="Some secondary"){
        cleaned_data[k,7]<-0
      }else if(phen_all_rct_use[k,kl]=="More than secondary"){
        cleaned_data[k,7]<-0
      }
    }
    if(is.na(phen_all_rct_use[k,kl])==T){
      cleaned_data[k,7]<-NA
    }
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,8]<-as.numeric(sp[k,1])
    cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","education","sp","village")
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
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+education+(1|village),
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
colnames(results_ftrs)<-c("Effect_size","Pvalue")
results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
length(which(results_ftrs$Pvalue<0.05))
length(which(results_ftrs$FDR<0.05))


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Education(4th-6th grade)") 

print("middle education")  

#Secondary education

esz<-array(0,dim=c(dim(data_transformed)[1],1))
pv_sp<-array(0,dim=c(dim(data_transformed)[1],1)) #p value 
kl<-5

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
    cleaned_data[k,6]<-as.numeric(st[k])
    #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
    if(is.na(phen_all_rct_use[k,kl])==F){
      if(phen_all_rct_use[k,kl]=="Have not completed any type of school"){
        cleaned_data[k,7]<-0
      }else if(phen_all_rct_use[k,kl]=="1st grade"){
        cleaned_data[k,7]<-0
      }else if(phen_all_rct_use[k,kl]=="2nd grade"){
        cleaned_data[k,7]<-0
      }else if(phen_all_rct_use[k,kl]=="3rd grade"){
        cleaned_data[k,7]<-0
      }else if(phen_all_rct_use[k,kl]=="4th grade"){
        cleaned_data[k,7]<-0
      }else if(phen_all_rct_use[k,kl]=="5th grade"){
        cleaned_data[k,7]<-0
      }else if(phen_all_rct_use[k,kl]=="6th grade"){
        cleaned_data[k,7]<-0
      }else if(phen_all_rct_use[k,kl]=="Some secondary"){
        cleaned_data[k,7]<-1
      }else if(phen_all_rct_use[k,kl]=="More than secondary"){
        cleaned_data[k,7]<-1
      }
    }
    if(is.na(phen_all_rct_use[k,kl])==T){
      cleaned_data[k,7]<-NA
    }
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,8]<-as.numeric(sp[k,1])
    cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","education","sp","village")
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
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+education+(1|village),
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
colnames(results_ftrs)<-c("Effect_size","Pvalue")
results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
length(which(results_ftrs$Pvalue<0.05))
length(which(results_ftrs$FDR<0.05))


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Education(secondary)") 

print("secondary education")  

## Diarrhea

#colnames(phen_all_use)

esz<-array(0,dim=c(dim(data_transformed)[1],1))
pv_sp<-array(0,dim=c(dim(data_transformed)[1],1)) #p value 

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
    cleaned_data[k,6]<-as.numeric(st[k])
    cleaned_data[k,7]<-ifelse(phen_all_use[k,381]=="No",0,1)
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,8]<-as.numeric(sp[k,1])
    cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","cough","sp","village")
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
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+cough+(1|village),
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
colnames(results_ftrs)<-c("Effect_size","Pvalue")
results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
length(which(results_ftrs$Pvalue<0.05))
length(which(results_ftrs$FDR<0.05))

effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Diarrhea")


print("Diarrhea")

##Bristol stool

esz<-array(0,dim=c(dim(data_transformed)[1],1))
pv_sp<-array(0,dim=c(dim(data_transformed)[1],1)) #p value 

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
    cleaned_data[k,6]<-as.numeric(st[k])
    cleaned_data[k,7]<-ifelse(phen_all_use[k,381]=="No",0,1)
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,8]<-as.numeric(sp[k,1])
    cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","cough","sp","village")
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
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+cough+bristol+(1|village),
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
colnames(results_ftrs)<-c("Effect_size","Pvalue")
results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
length(which(results_ftrs$Pvalue<0.05))
length(which(results_ftrs$FDR<0.05))

effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Bristol")


print("Bristol")

write.csv(effect_size_phen,'effect_size_phen_9_lmer_pwy.csv')
write.csv(pval_phen,'pval_phen_9_lmer_pwy.csv')
write.csv(fdr_phen,'fdr_phen_9_lmer_pwy.csv')


