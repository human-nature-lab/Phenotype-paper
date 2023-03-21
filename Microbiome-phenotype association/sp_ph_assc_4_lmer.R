library(lmerTest)

st<-read.csv('bristol.csv',row.names = 1)
st<-as.matrix.data.frame(st)
#phen_all_use<-read.csv('phen_all_use.csv',row.names = 1,header=T)
phen_all_use<-read.csv('phen_all_use6_all.csv',row.names = 1,header=T)
#data_transformed<-read.csv('data_transformed.csv',row.names = 1,header=T)
data_transformed<-read.csv('data_transformed_1187.csv',row.names = 1,header=T)
phen_all_rct_use<-read.csv('phen_all_rct_use22.csv',row.names = 1,header=T)
mb_samp_sp<-read.csv('mb_samp_sp.csv',row.names = 1)
data_transformed<-as.matrix.data.frame(data_transformed)
ind_t<-read.csv('ind_t.csv')
ind_t<-ind_t[,1]

data_transformed<-data_transformed[ind_t,]


##GAD7 and PHQ9 all sub components

#Discretized

for(kl in c(183:189,192:200)){
  
  esz<-array(0,dim=c(dim(mb_samp_sp)[1],1))
  pv_sp<-array(0,dim=c(dim(mb_samp_sp)[1],1)) #p value 
  
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
      #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
      if(is.na(phen_all_use[k,kl])==F){
        if(phen_all_use[k,kl]=="Not at all"){
          cleaned_data[k,7]<-0
        }else if(phen_all_use[k,kl]=="Several days"){
          cleaned_data[k,7]<-1
        }else if(phen_all_use[k,kl]=="Over half the days"){
          cleaned_data[k,7]<-2
        }else if(phen_all_use[k,kl]=="Nearly every day"){
          cleaned_data[k,7]<-3
        }
      }
      if(is.na(phen_all_use[k,kl])==T){
        cleaned_data[k,7]<-NA
      }
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,8]<-as.numeric(sp[k,1])
      cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","a1c","sp","village")
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
    
    
    ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c+(1|village),
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
  
  if(kl==183){
    effect_size_phen<-results_ftrs$Effect_size
    pval_phen<-results_ftrs$Pvalue
    fdr_phen<-results_ftrs$FDR
  }else{
    effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
    pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
    fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)
  }
  
}

colnames(effect_size_phen)<-c("GAD7-Q1","GAD7-Q2","GAD7-Q3","GAD7-Q4","GAD7-Q5","GAD7-Q6","GAD7-Q7","PHQ9-Q1","PHQ9-Q2","PHQ9-Q3","PHQ9-Q4","PHQ9-Q5","PHQ9-Q6","PHQ9-Q7","PHQ9-Q8","PHQ9-Q9")
lim<-dim(effect_size_phen)[2]

#factors

for(kl in c(183:189,192:200)){
  
  esz<-array(0,dim=c(dim(mb_samp_sp)[1],3))
  pv_sp<-array(0,dim=c(dim(mb_samp_sp)[1],3)) #p value 
  
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
    cleaned_data2<-cbind(cleaned_data2,"ment_cat"=factor(phen_all_use[ind_temp,kl])) #Phenotype)
    cleaned_data2 <- within(cleaned_data2, ment_cat <- relevel(ment_cat, ref = "Not at all"))
    
    
    #s1 = lm(
    #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
    #  data = cleaned_data2
    #)
    
    cleaned_data2$village<-factor(cleaned_data2$village)
    
    
    ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+ment_cat+(1|village),
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
  

    effect_size_phen<-cbind(effect_size_phen,results_ftrs[,1:3])
    pval_phen<-cbind(pval_phen,results_ftrs[,4:6])
    fdr_phen<-cbind(fdr_phen,results_ftrs[,7:9])
    colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-rownames(coef(mss))[8:10]
    lim<-dim(effect_size_phen)[2]
  
}



##Gad7 score

esz<-array(0,dim=c(dim(mb_samp_sp)[1],1))
pv_sp<-array(0,dim=c(dim(mb_samp_sp)[1],1)) #p value 

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
    cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,190]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,8]<-as.numeric(sp[k,1])
    cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","a1c","sp","village")
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
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c+(1|village),
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

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("GAD7 score")

lim<-dim(effect_size_phen)[2]

##GAD7 status

esz<-array(0,dim=c(dim(mb_samp_sp)[1],3))
pv_sp<-array(0,dim=c(dim(mb_samp_sp)[1],3)) #p value 

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
  cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
  cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+GAD7_cat+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[8,5]
  
  effect.size<-coef(mss)[8:10,1]
  pval2<-coef(mss)[8:10,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size_mild","Effect_size_moderate","Effect_size_severe","Pvalue_mild","Pvalue_moderate","Pvalue_severe")
results_ftrs$FDR_mild <- p.adjust(results_ftrs$Pvalue_mild,method = "BH")
results_ftrs$FDR_moderate <- p.adjust(results_ftrs$Pvalue_moderate,method = "BH")
results_ftrs$FDR_severe <- p.adjust(results_ftrs$Pvalue_severe,method = "BH")
length(which(results_ftrs$Pvalue_mild<0.05))
length(which(results_ftrs$Pvalue_moderate<0.05))
length(which(results_ftrs$Pvalue_severe<0.05))
length(which(results_ftrs$FDR_mild<0.05))
length(which(results_ftrs$FDR_moderate<0.05))
length(which(results_ftrs$FDR_severe<0.05))
#length(which(results_ftrs$Pvalue<0.05))
#length(which(results_ftrs$FDR<0.05))

effect_size_phen<-cbind(effect_size_phen,results_ftrs[,1:3])
pval_phen<-cbind(pval_phen,results_ftrs[,4:6])
fdr_phen<-cbind(fdr_phen,results_ftrs[,7:9])
colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-rownames(coef(mss))[8:10]
#lim<-dim(effect_size_phen)[2]


#colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-c("GAD7(mild)","GAD7(moderate)","GAD7(severe)")




##PHQ9 score

esz<-array(0,dim=c(dim(mb_samp_sp)[1],1))
pv_sp<-array(0,dim=c(dim(mb_samp_sp)[1],1)) #p value 

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
    cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,201]))
    #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
    #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
    cleaned_data[k,8]<-as.numeric(sp[k,1])
    cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
  }
  cleaned_data<-as.data.frame(cleaned_data)
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","a1c","sp","village")
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
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c+(1|village),
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

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("PHQ9 score")

lim<-dim(effect_size_phen)[2]

##PHQ9 status

esz<-array(0,dim=c(dim(mb_samp_sp)[1],3))
pv_sp<-array(0,dim=c(dim(mb_samp_sp)[1],3)) #p value 

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
  cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,202])) #Phenotype)
  cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
  
  
  #s1 = lm(
  #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
  #  data = cleaned_data2
  #)
  
  cleaned_data2$village<-factor(cleaned_data2$village)
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+GAD7_cat+(1|village),
            data = cleaned_data2)
  
  mss<-summary(ms1)
  
  coef(mss)[8,5]
  
  effect.size<-coef(mss)[8:10,1]
  pval2<-coef(mss)[8:10,5]#Pvalue of phenotype
  
  esz[i,]<-effect.size
  pv_sp[i,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  
}

results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size_mild","Effect_size_moderate","Effect_size_severe","Pvalue_mild","Pvalue_moderate","Pvalue_severe")
results_ftrs$FDR_mild <- p.adjust(results_ftrs$Pvalue_mild,method = "BH")
results_ftrs$FDR_moderate <- p.adjust(results_ftrs$Pvalue_moderate,method = "BH")
results_ftrs$FDR_severe <- p.adjust(results_ftrs$Pvalue_severe,method = "BH")
length(which(results_ftrs$Pvalue_mild<0.05))
length(which(results_ftrs$Pvalue_moderate<0.05))
length(which(results_ftrs$Pvalue_severe<0.05))
length(which(results_ftrs$FDR_mild<0.05))
length(which(results_ftrs$FDR_moderate<0.05))
length(which(results_ftrs$FDR_severe<0.05))
#length(which(results_ftrs$Pvalue<0.05))
#length(which(results_ftrs$FDR<0.05))

effect_size_phen<-cbind(effect_size_phen,results_ftrs[,1:3])
pval_phen<-cbind(pval_phen,results_ftrs[,4:6])
fdr_phen<-cbind(fdr_phen,results_ftrs[,7:9])
colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-rownames(coef(mss))[8:10]



#colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-c("PHQ9(mild)","PHQ9(moderate)","PHQ9(severe)")

write.csv(effect_size_phen,'effect_size_phen_4_lmer.csv')
write.csv(pval_phen,'pval_phen_4_lmer.csv')
write.csv(fdr_phen,'fdr_phen_4_lmer.csv')

