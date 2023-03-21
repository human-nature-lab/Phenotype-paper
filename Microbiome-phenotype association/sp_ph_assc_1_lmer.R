
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


##All physiological variables

kl<-32

for(kl in c(32,35,36)){

esz<-array(0,dim=c(dim(mb_samp_sp)[1],1))
pv_sp<-array(0,dim=c(dim(mb_samp_sp)[1],1)) #p value 

for(i in 1:dim(data_transformed)[1]){#dim(data_transformed)[1]
  
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
length(which(results_ftrs$Pvalue[1:500]<0.05))
length(which(results_ftrs$FDR[1:500]<0.05))

if(kl==32){
  effect_size_phen<-results_ftrs$Effect_size
  pval_phen<-results_ftrs$Pvalue
  fdr_phen<-results_ftrs$FDR
}else{
  effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
  pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)
}

}

colnames(effect_size_phen)<-c("A1c","Systolic","Diastolic")


#MAP

phen_map<-as.numeric(as.character(phen_all_use$mb_d0400))+(as.numeric(as.character(phen_all_use$mb_d0300))-as.numeric(as.character(phen_all_use$mb_d0400)))/3

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
  cleaned_data[k,7]<-as.numeric(as.character(phen_map[k]))
  #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
  #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
  cleaned_data[k,8]<-as.numeric(sp[k,1])
  cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
}
cleaned_data<-as.data.frame(cleaned_data)
colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","map","sp","village")
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


ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+map+(1|village),
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

colnames(effect_size_phen)[4]<-c("MAP")


##BMI


esz<-array(0,dim=c(dim(mb_samp_sp)[1],1))
pv_sp<-array(0,dim=c(dim(mb_samp_sp)[1],1)) #p value 
  
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
    #cleaned_data2<-cbind(cleaned_data2,"GAD7_cat"=factor(phen_all_use[ind_temp,191])) #Phenotype)
    #cleaned_data2 <- within(cleaned_data2, GAD7_cat <- relevel(GAD7_cat, ref = "minimal or none"))
    
    
    #s1 = lm(
    #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
    #  data = cleaned_data2
    #)
    
    cleaned_data2$village<-factor(cleaned_data2$village)
    
    
    ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+(1|village),
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
  results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
  length(which(results_ftrs$Pvalue<0.05))
  length(which(results_ftrs$FDR<0.05))
  
    effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
    pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
    fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)
  
colnames(effect_size_phen)[5]<-c("BMI")

##Rest physiological

    for(kl in c(71:74,37)){
      
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
          cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
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
      
    }
    
    colnames(effect_size_phen)[6:10]<-c("Heart rate","Perfusion index","Pulse","O2 sat","Hb total")
    
    
##Cough for 1 month

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
        cleaned_data[k,7]<-ifelse(phen_all_use[k,16]=="No",0,1)
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
    
    colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Cough(1m)")
    
##overall health
    

    esz<-array(0,dim=c(dim(mb_samp_sp)[1],4))
    pv_sp<-array(0,dim=c(dim(mb_samp_sp)[1],4)) #p value 
    
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
        #cleaned_data[k,8]<-ifelse(phen_all_use[k,15]=="moderate",1,0)  #Phenotype
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
      cleaned_data2<-cbind(cleaned_data2,"health"=factor(phen_all_use[ind_temp,15])) #Phenotype)
      cleaned_data2 <- within(cleaned_data2, health <- relevel(health, ref = "Good"))
      
      
      #s1 = lm(
      #  sp~age+sex+batch_effect+dna_conc+BMI+bristol+a1c,
      #  data = cleaned_data2
      #)
      
      cleaned_data2$village<-factor(cleaned_data2$village)
      
      
      ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+bristol+health+(1|village),
                data = cleaned_data2)
      
      mss<-summary(ms1)
      
      coef(mss)[8,5]
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
    
colnames(effect_size_phen)[12:15]<-c("Excellent","Fair","Poor","Very good") 
    
    write.csv(effect_size_phen,'effect_size_phen_1_lmer.csv')
    write.csv(pval_phen,'pval_phen_1_lmer.csv')
    write.csv(fdr_phen,'fdr_phen_1_lmer.csv')
    
    
    
