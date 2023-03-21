library(lmerTest)

st<-read.csv('bristol.csv',row.names = 1)
st<-as.matrix.data.frame(st)
#phen_all_use<-read.csv('phen_all_use.csv',row.names = 1,header=T)
phen_all_use<-read.csv('phen_all_use6_all.csv',row.names = 1,header=T)
#data_transformed<-read.csv('data_transformed.csv',row.names = 1,header=T)
data_transformed<-read.csv('pwy_all_use.csv',row.names = 1,header=T)
phen_all_rct_use<-read.csv('phen_all_rct_use22.csv',row.names = 1,header=T)
data_transformed<-as.matrix.data.frame(data_transformed)

##Animals

for(kl in 155:180){
  
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
      cleaned_data[k,7]<-ifelse(is.na(phen_all_use[k,kl])==F,1,0)
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,8]<-as.numeric(sp[k,1])
      cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","animals","sp","village")
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
    
    
    ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+animals+(1|village),
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
  
  if(kl==155){
    effect_size_phen<-results_ftrs$Effect_size
    pval_phen<-results_ftrs$Pvalue
    fdr_phen<-results_ftrs$FDR
  }else{
    effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
    pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
    fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)
  }
  
}

colnames(effect_size_phen)<-c("Cat","Dog","Parakeet","Rabbit","Horse","Mice","None pet","Cow","Goat","Pig","Chicken","Duck","Turkey","Sheep","Geese","None farm","Bat","Lizard","Monkey","Snake","Bird","Possum","Rat","Squirrel","Other wild","None wild")

lim<-dim(effect_size_phen)[2]

##Food

for(kl in 75:93){
  
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
      #cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
      if(is.na(phen_all_use[k,kl])==F){
        if(phen_all_use[k,kl]=="Never/rarely"){#Change to 1/50,1/10,4/7,1
          cleaned_data[k,7]<-1/50
        }else if((phen_all_use[k,kl]=="A few days per month")|(phen_all_use[k,kl]=="A few times per month")){
          cleaned_data[k,7]<-1/10
        }else if(phen_all_use[k,kl]=="A few days per week"){
          cleaned_data[k,7]<-4/7
        }else if(phen_all_use[k,kl]=="Every day"){
          cleaned_data[k,7]<-1
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
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","food","sp","village")
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
    
    
    ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+food+(1|village),
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
  
}
colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-c("Beans","Tortillas","Rice","Bread","Milk","Yogurt","Cream/butter","Cheese","Eggs","Vegetables","Fruits","Natural juice","Chicken","Beef/Pork","Ham/sausages/hotdog","Fish","Soda","Fruit juice","Chips")

##Diet diversity score

esz<-array(0,dim=c(dim(data_transformed)[1],1))
pv_sp<-array(0,dim=c(dim(data_transformed)[1],1)) #p value 
kl<-382

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
      cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
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
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","food","sp","village")
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
  
  
  ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+food+(1|village),
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

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-"DDS"



write.csv(effect_size_phen,'effect_size_phen_5_lmer_pwy.csv')
write.csv(pval_phen,'pval_phen_5_lmer_pwy.csv')
write.csv(fdr_phen,'fdr_phen_5_lmer_pwy.csv')







