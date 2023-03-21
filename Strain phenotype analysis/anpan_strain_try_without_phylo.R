library(phyr)
library(lmerTest)

sp_name<-read.csv('species_vs_sample_name_2285.csv',row.names=1)
phen_all_rct_use<-read.csv('phen_all_rct_use22.csv',row.names = 1,header=T)
mb_list<-read.csv('mb_list.csv',row.names = 1)
phen_all_use<-read.csv('phen_all_use.csv',row.names = 1)
sp_ind<-read.csv('sp_use_strain.csv',row.names = 1)
data_transformed<-read.csv('data_transformed.csv',row.names = 1,header=T)
data_transformed<-as.matrix.data.frame(data_transformed)
ind_t<-read.csv('ind_t.csv')
ind_t<-ind_t[,1]

data_transformed<-data_transformed[ind_t,]

ind_sp<-read.csv('strain_sp_list.csv',row.names=1)
ind_sp2<-read.csv('strain_sp_list_sp.csv',row.names=1)

data_transformed2<-data_transformed[ind_sp2[,1],]

st<-read.csv('bristol.csv',row.names = 1)
st<-as.matrix.data.frame(st)

sp_name2<-sp_name[ind_sp2[,1],1]

setwd("/WORKAREA/home/sv437/strain_assc_anpan/tree_files/")

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

for(kl in 143:152){
  for(j in 1:dim(phen_all_use)[1]){
    if(is.na(phen_all_use[j,kl])==F){
      if(phen_all_use[j,kl]=="Dont_Know"){
        phen_all_use[j,kl]<-NA
      }
    }
  }
}


phen_map<-as.numeric(as.character(phen_all_use$mb_d0400))+(as.numeric(as.character(phen_all_use$mb_d0300))-as.numeric(as.character(phen_all_use$mb_d0400)))/3

phen_all_use2<-cbind(phen_all_use,phen_map)

contt_ind<-c(32,35,36,40,71:74,37,153,190,201,182,364,203:207,360:363,365:369,228:230,306:309,317:359)
phen_all_use2<-array(NA,dim=c(dim(phen_all_use)[1],dim(phen_all_use)[2]))

setwd("/WORKAREA/home/sv437/strain_assc_anpan/tree_files/")


file_list<-list.files(path=".")

for(kl in c(32,35,36)){
  
  esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
  pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 
  
  for(io in 1:dim(ind_sp)[1]){
    
    tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
    
    cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
    if(dim(cor_mat)[1]>15){
      sp<-as.data.frame(data_transformed[ind_sp2[io,],])
      
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
      colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","a1c","species","village")
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
      
      cleaned_data2<-as.data.frame(cleaned_data2)
      
      rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
      rownames(cleaned_data)<-rownames(phen_all_use)
      
      ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
      
      ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
      ind_sub2_ind<-which(is.na(ind_sub)==F)
      
      #tr3<-tr2[ind_sub2_ind]
      cleaned_data3<-cleaned_data2[ind_sub2,]
      
      sigma_phylo = 1
      n=length(ind_sub)
      
      true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
      
      #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
      cleaned_data4<-cleaned_data3
      
      cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
      cleaned_data3$village<-as.factor(cleaned_data3$village)
      
      cleaned_data4$a1c<-cleaned_data4$a1c##+true_phylo_effects[ind_sub2_ind]
      cleaned_data4$village<-as.factor(cleaned_data4$village)
      
      #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
      #                                      data = cleaned_data3,
      #                                      family = "gaussian")#use ordinal or gaussian
      
      
      #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
      #                   data = cleaned_data4,
      #                   family = "gaussian")#use ordinal or gaussian
      
      ms1<-lmer(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
                data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
      
      mss<-summary(ms1)
      
      coef(mss)
      
      effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
      pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
      
      #summary(mod)
      
      #summary(mod2)
      
      #mod2$B[8]
      #mod2$B.pvalue[8]
      #print(io)
      
      
      #effect.size<-mod2$B[8]
      #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
      
      esz[io,]<-effect.size
      pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    }
  }
  
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  #results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
  results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
  length(which(results_ftrs$Pvalue<0.05))
  length(which(results_ftrs$FDR<0.05))
  
  rownames(results_ftrs)<-sp_name2
  
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

##########


#Continuation of anpan_strain_try.R

phen_map<-as.numeric(as.character(phen_all_use$mb_d0400))+(as.numeric(as.character(phen_all_use$mb_d0300))-as.numeric(as.character(phen_all_use$mb_d0400)))/3



esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 

for(io in 1:dim(ind_sp)[1]){
  
  tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
  
  cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
  if(dim(cor_mat)[1]>15){
    sp<-as.data.frame(data_transformed[ind_sp2[io,],])
    
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
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","a1c","species","village")
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
    
    cleaned_data2<-as.data.frame(cleaned_data2)
    
    rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
    rownames(cleaned_data)<-rownames(phen_all_use)
    
    ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
    
    ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
    ind_sub2_ind<-which(is.na(ind_sub)==F)
    
    #tr3<-tr2[ind_sub2_ind]
    cleaned_data3<-cleaned_data2[ind_sub2,]
    
    sigma_phylo = 1
    n=length(ind_sub)
    
    true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
    
    #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
    cleaned_data4<-cleaned_data3
    
    cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
    cleaned_data3$village<-as.factor(cleaned_data3$village)
    
    cleaned_data4$a1c<-cleaned_data4$a1c#+true_phylo_effects[ind_sub2_ind]
    cleaned_data4$village<-as.factor(cleaned_data4$village)
    
    #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
    #                                      data = cleaned_data3,
    #                                      family = "gaussian")#use ordinal or gaussian
    
    
    #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
    #                   data = cleaned_data4,
    #                   family = "gaussian")#use ordinal or gaussian
    
    ms1<-lmer(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
              data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
    
    mss<-summary(ms1)
    
    coef(mss)
    
    effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
    pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
    
    #summary(mod)
    
    #summary(mod2)
    
    #mod2$B[8]
    #mod2$B.pvalue[8]
    #print(io)
    
    
    #effect.size<-mod2$B[8]
    #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
    
    esz[io,]<-effect.size
    pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  }
}


results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")
#results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
length(which(results_ftrs$Pvalue<0.05))
length(which(results_ftrs$FDR<0.05))

rownames(results_ftrs)<-sp_name2


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)


colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-"MAP"

##BMI

esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 

for(io in 1:dim(ind_sp)[1]){
  
  tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
  
  cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
  if(dim(cor_mat)[1]>15){
    sp<-as.data.frame(data_transformed[ind_sp2[io,],])
    
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],8))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(st[k])
      #cleaned_data[k,7]<-as.numeric(as.character(phen_map[k]))
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,7]<-as.numeric(sp[k,1])
      cleaned_data[k,8]<-as.numeric(phen_all_use[k,11])
      
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","species","village")
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
    
    cleaned_data2<-as.data.frame(cleaned_data2)
    
    rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
    rownames(cleaned_data)<-rownames(phen_all_use)
    
    ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
    
    ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
    ind_sub2_ind<-which(is.na(ind_sub)==F)
    
    #tr3<-tr2[ind_sub2_ind]
    cleaned_data3<-cleaned_data2[ind_sub2,]
    
    sigma_phylo = 1
    n=length(ind_sub)
    
    true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
    
    #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
    cleaned_data4<-cleaned_data3
    
    cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
    cleaned_data3$village<-as.factor(cleaned_data3$village)
    
    cleaned_data4$BMI<-cleaned_data4$BMI#+true_phylo_effects[ind_sub2_ind]
    cleaned_data4$village<-as.factor(cleaned_data4$village)
    
    #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
    #                                      data = cleaned_data3,
    #                                      family = "gaussian")#use ordinal or gaussian
    
    
    #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
    #                   data = cleaned_data4,
    #                   family = "gaussian")#use ordinal or gaussian
    
    ms1<-lmer(BMI ~ age+sex+batch_effect+dna_conc+bristol+species  +(1|village), 
              data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
    
    mss<-summary(ms1)
    
    coef(mss)
    
    effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
    pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
    
    #summary(mod)
    
    #summary(mod2)
    
    #mod2$B[8]
    #mod2$B.pvalue[8]
    #print(io)
    
    
    #effect.size<-mod2$B[8]
    #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
    
    esz[io,]<-effect.size
    pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
  }
}


results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")
#results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
length(which(results_ftrs$Pvalue<0.05))
length(which(results_ftrs$FDR<0.05))

rownames(results_ftrs)<-sp_name2


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-"BMI"

lim<-dim(effect_size_phen)[2]

## Rest physiological

for(kl in c(71:74,37)){
  
  esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
  pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 
  
  for(io in 1:dim(ind_sp)[1]){
    
    tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
    
    cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
    if(dim(cor_mat)[1]>15){
      sp<-as.data.frame(data_transformed[ind_sp2[io,],])
      
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
      colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","phy","species","village")
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
      
      cleaned_data2<-as.data.frame(cleaned_data2)
      
      rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
      rownames(cleaned_data)<-rownames(phen_all_use)
      
      ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
      
      ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
      ind_sub2_ind<-which(is.na(ind_sub)==F)
      
      #tr3<-tr2[ind_sub2_ind]
      cleaned_data3<-cleaned_data2[ind_sub2,]
      
      sigma_phylo = 1
      n=length(ind_sub)
      
      true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
      
      #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
      cleaned_data4<-cleaned_data3
      
      cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
      cleaned_data3$village<-as.factor(cleaned_data3$village)
      
      cleaned_data4$phy<-cleaned_data4$phy#+true_phylo_effects[ind_sub2_ind]
      cleaned_data4$village<-as.factor(cleaned_data4$village)
      
      #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
      #                                      data = cleaned_data3,
      #                                      family = "gaussian")#use ordinal or gaussian
      
      
      #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
      #                   data = cleaned_data4,
      #                   family = "gaussian")#use ordinal or gaussian
      
      ms1<-lmer(phy ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
                data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
      
      mss<-summary(ms1)
      
      coef(mss)
      
      effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
      pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
      
      #summary(mod)
      
      #summary(mod2)
      
      #mod2$B[8]
      #mod2$B.pvalue[8]
      #print(io)
      
      
      #effect.size<-mod2$B[8]
      #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
      
      esz[io,]<-effect.size
      pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    }
  }
  
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  #results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
  results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
  length(which(results_ftrs$Pvalue<0.05))
  length(which(results_ftrs$FDR<0.05))
  
  rownames(results_ftrs)<-sp_name2
  
  
  effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
  pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)
  
  
  
}

colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-c("Heart rate","Perfusion index","Pulse","O2 sat","Hb total")



##Cough for 1 month

esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 

for(io in 1:dim(ind_sp)[1]){
  
  tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
  
  cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
  if(dim(cor_mat)[1]>15){
    sp<-as.data.frame(data_transformed[ind_sp2[io,],])
    
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
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","phy","species","village")
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
    
    cleaned_data2<-as.data.frame(cleaned_data2)
    
    rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
    rownames(cleaned_data)<-rownames(phen_all_use)
    
    ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
    
    ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
    ind_sub2_ind<-which(is.na(ind_sub)==F)
    
    #tr3<-tr2[ind_sub2_ind]
    cleaned_data3<-cleaned_data2[ind_sub2,]
    
    sigma_phylo = 1
    n=length(ind_sub)
    
    true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
    
    #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
    cleaned_data4<-cleaned_data3
    
    cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
    cleaned_data3$village<-as.factor(cleaned_data3$village)
    
    cleaned_data4$phy<-cleaned_data4$phy#+true_phylo_effects[ind_sub2_ind]
    cleaned_data4$village<-as.factor(cleaned_data4$village)
    
    #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
    #                                      data = cleaned_data3,
    #                                      family = "gaussian")#use ordinal or gaussian
    
    
    #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
    #                   data = cleaned_data4,
    #                   family = "gaussian")#use ordinal or gaussian
    if(sum(cleaned_data4$phy)>0){
    ms1<-lmer(phy ~ age+sex+batch_effect+dna_conc+BMI+bristol+species+(1|village), 
              data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
    
    mss<-summary(ms1)
    
    coef(mss)
    
    effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
    pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
    
    #summary(mod)
    
    #summary(mod2)
    
    #mod2$B[8]
    #mod2$B.pvalue[8]
    #print(io)
    
    
    #effect.size<-mod2$B[8]
    #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
    
    esz[io,]<-effect.size
    pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    }
  }
}


results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")
#results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
length(which(results_ftrs$Pvalue<0.05))
length(which(results_ftrs$FDR<0.05))

rownames(results_ftrs)<-sp_name2


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-"Cough(1m)"

lim<-dim(effect_size_phen)[2]

##Chronic conditions

for(kl in 20:29){
  
  esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
  pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 
  
  for(io in 1:dim(ind_sp)[1]){
    
    tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
    
    cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
    if(dim(cor_mat)[1]>15){
      sp<-as.data.frame(data_transformed[ind_sp2[io,],])
      
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
      colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","chronic","species","village")
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
      
      cleaned_data2<-as.data.frame(cleaned_data2)
      
      rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
      rownames(cleaned_data)<-rownames(phen_all_use)
      
      ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
      
      ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
      ind_sub2_ind<-which(is.na(ind_sub)==F)
      
      #tr3<-tr2[ind_sub2_ind]
      cleaned_data3<-cleaned_data2[ind_sub2,]
      
      sigma_phylo = 1
      n=length(ind_sub)
      
      true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
      
      #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
      cleaned_data4<-cleaned_data3
      
      cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
      cleaned_data3$village<-as.factor(cleaned_data3$village)
      
      cleaned_data4$chronic<-cleaned_data4$chronic#+true_phylo_effects[ind_sub2_ind]
      cleaned_data4$village<-as.factor(cleaned_data4$village)
      
      #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
      #                                      data = cleaned_data3,
      #                                      family = "gaussian")#use ordinal or gaussian
      
      
      #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
      #                   data = cleaned_data4,
      #                   family = "gaussian")#use ordinal or gaussian
      if((sum(cleaned_data4$chronic)>2)&(dim(cleaned_data4)[1]>30)){
      ms1<-lmer(chronic ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
                data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
      
      mss<-summary(ms1)
      
      coef(mss)
      
      effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
      pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
      
      #summary(mod)
      
      #summary(mod2)
      
      #mod2$B[8]
      #mod2$B.pvalue[8]
      #print(io)
      
      
      #effect.size<-mod2$B[8]
      #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
      
      esz[io,]<-effect.size
      pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
      }
    }
  }
  
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  #results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
  results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
  length(which(results_ftrs$Pvalue<0.05))
  length(which(results_ftrs$FDR<0.05))
  
  rownames(results_ftrs)<-sp_name2
  
  
  effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
  pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)
  
  
  
}

colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-c("Diabetes","Allergies","Cystic fibrosis","Heart disease","Endocrine illness","Renal failure","Asthma","Stomach illness","Intestinal illness","Arthritis")

lim<-dim(effect_size_phen)[2]
print("medication")
## Medication


for(kl in 136:142){
  
  esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
  pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 
  
  for(io in 1:dim(ind_sp)[1]){
    
    tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
    
    cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
    if(dim(cor_mat)[1]>15){
      sp<-as.data.frame(data_transformed[ind_sp2[io,],])
      
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
      colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","chronic","species","village")
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
      
      cleaned_data2<-as.data.frame(cleaned_data2)
      
      rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
      rownames(cleaned_data)<-rownames(phen_all_use)
      
      ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
      
      ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
      ind_sub2_ind<-which(is.na(ind_sub)==F)
      
      #tr3<-tr2[ind_sub2_ind]
      cleaned_data3<-cleaned_data2[ind_sub2,]
      
      sigma_phylo = 1
      n=length(ind_sub)
      
      true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
      
      #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
      cleaned_data4<-cleaned_data3
      
      cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
      cleaned_data3$village<-as.factor(cleaned_data3$village)
      
      cleaned_data4$chronic<-cleaned_data4$chronic#+true_phylo_effects[ind_sub2_ind]
      cleaned_data4$village<-as.factor(cleaned_data4$village)
      
      #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
      #                                      data = cleaned_data3,
      #                                      family = "gaussian")#use ordinal or gaussian
      
      
      #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
      #                   data = cleaned_data4,
      #                   family = "gaussian")#use ordinal or gaussian
      if((sum(cleaned_data4$chronic)>2)&(dim(cleaned_data4)[1]>30)){
      ms1<-lmer(chronic ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
                data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
      
      mss<-summary(ms1)
      
      coef(mss)
      
      effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
      pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
      
      #summary(mod)
      
      #summary(mod2)
      
      #mod2$B[8]
      #mod2$B.pvalue[8]
      #print(io)
      
      
      #effect.size<-mod2$B[8]
      #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
      
      esz[io,]<-effect.size
      pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
      }
    }
  }
  
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  #results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
  results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
  length(which(results_ftrs$Pvalue<0.05))
  length(which(results_ftrs$FDR<0.05))
  
  rownames(results_ftrs)<-sp_name2
  
  
  effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
  pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)
  
}
  


colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-c("Pain killers","Antibiotics","Anti-diarrheal","Anti-parasitic","Anti-fungal","Vitamins","Anti-hypertensive")

lim<-dim(effect_size_phen)[2]
print("personality types")
## Personality types


for(kl in 143:152){
  
  esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
  pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 
  
  for(io in 1:dim(ind_sp)[1]){
    
    tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
    
    cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
    if(dim(cor_mat)[1]>15){
      sp<-as.data.frame(data_transformed[ind_sp2[io,],])
      
      cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
      for(k in 1:dim(phen_all_use)[1]){
        cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
        cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
        cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
        cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
        cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
        cleaned_data[k,6]<-as.numeric(st[k])
        if(is.na(phen_all_use[k,kl])==F){
          if(phen_all_use[k,kl]=="Disagree strongly"){
            cleaned_data[k,7]<-1
          }else if(phen_all_use[k,kl]=="Disagree a little"){
            cleaned_data[k,7]<-2
          }else if(phen_all_use[k,kl]=="Neither agree nor disagree"){
            cleaned_data[k,7]<-3
          }else if(phen_all_use[k,kl]=="Agree a little"){
            cleaned_data[k,7]<-4
          }else if(phen_all_use[k,kl]=="Agree strongly"){
            cleaned_data[k,7]<-5
          }
        }
        #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
        #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
        cleaned_data[k,8]<-as.numeric(sp[k,1])
        cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
        
      }
      cleaned_data<-as.data.frame(cleaned_data)
      colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","chronic","species","village")
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
      
      cleaned_data2<-as.data.frame(cleaned_data2)
      
      rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
      rownames(cleaned_data)<-rownames(phen_all_use)
      
      ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
      
      ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
      ind_sub2_ind<-which(is.na(ind_sub)==F)
      
      #tr3<-tr2[ind_sub2_ind]
      cleaned_data3<-cleaned_data2[ind_sub2,]
      
      sigma_phylo = 1
      n=length(ind_sub)
      
      true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
      
      #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
      cleaned_data4<-cleaned_data3
      
      cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
      cleaned_data3$village<-as.factor(cleaned_data3$village)
      
      cleaned_data4$chronic<-cleaned_data4$chronic#+true_phylo_effects[ind_sub2_ind]
      cleaned_data4$village<-as.factor(cleaned_data4$village)
      
      #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
      #                                      data = cleaned_data3,
      #                                      family = "gaussian")#use ordinal or gaussian
      
      
      #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
      #                   data = cleaned_data4,
      #                   family = "gaussian")#use ordinal or gaussian
      if((sum(cleaned_data4$chronic)>1)&(dim(cleaned_data4)[1]>30))
      ms1<-lmer(chronic ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
                data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
      
      mss<-summary(ms1)
      
      coef(mss)
      
      effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
      pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
      
      #summary(mod)
      
      #summary(mod2)
      
      #mod2$B[8]
      #mod2$B.pvalue[8]
      #print(io)
      
      
      #effect.size<-mod2$B[8]
      #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
      
      esz[io,]<-effect.size
      pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    }
  }
  
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  #results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
  results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
  length(which(results_ftrs$Pvalue<0.05))
  length(which(results_ftrs$FDR<0.05))
  
  rownames(results_ftrs)<-sp_name2
  
  
  effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
  pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)
  
  
  
}

colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-c("Reserved","Trusting","Lazy","Relaxed","Not creative","Outgoing","Fault others","Thorough job","Nervous","Openess")

lim<-dim(effect_size_phen)[2]

print("cognitive score")
#Cognitive score

kl<-153

esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 

for(io in 1:dim(ind_sp)[1]){
  
  tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
  
  cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
  if(dim(cor_mat)[1]>15){
    sp<-as.data.frame(data_transformed[ind_sp2[io,],])
    
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
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","phy","species","village")
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
    
    cleaned_data2<-as.data.frame(cleaned_data2)
    
    rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
    rownames(cleaned_data)<-rownames(phen_all_use)
    
    ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
    
    ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
    ind_sub2_ind<-which(is.na(ind_sub)==F)
    
    #tr3<-tr2[ind_sub2_ind]
    cleaned_data3<-cleaned_data2[ind_sub2,]
    
    sigma_phylo = 1
    n=length(ind_sub)
    
    true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
    
    #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
    cleaned_data4<-cleaned_data3
    
    cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
    cleaned_data3$village<-as.factor(cleaned_data3$village)
    
    cleaned_data4$phy<-cleaned_data4$phy#+true_phylo_effects[ind_sub2_ind]
    cleaned_data4$village<-as.factor(cleaned_data4$village)
    
    #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
    #                                      data = cleaned_data3,
    #                                      family = "gaussian")#use ordinal or gaussian
    
    
    #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
    #                   data = cleaned_data4,
    #                   family = "gaussian")#use ordinal or gaussian
    if((sum(cleaned_data4$phy)>0)&(dim(cleaned_data4)[1]>30)){
    ms1<-lmer(phy ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
              data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
    
    mss<-summary(ms1)
    
    coef(mss)
    
    effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
    pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
    
    #summary(mod)
    
    #summary(mod2)
    
    #mod2$B[8]
    #mod2$B.pvalue[8]
    #print(io)
    
    
    #effect.size<-mod2$B[8]
    #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
    
    esz[io,]<-effect.size
    pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    }
  }
}


results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")
#results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
length(which(results_ftrs$Pvalue<0.05))
length(which(results_ftrs$FDR<0.05))

rownames(results_ftrs)<-sp_name2


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)


colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Cognitive score")



lim<-dim(effect_size_phen)[2]
print("gad7 and phq9 sub components")
#GAD7 and PHQ9 sub components

for(kl in c(183:189,192:200)){
  
  esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
  pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 
  
  for(io in 1:dim(ind_sp)[1]){
    
    tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
    
    cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
    if(dim(cor_mat)[1]>15){
      sp<-as.data.frame(data_transformed[ind_sp2[io,],])
      
      cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
      for(k in 1:dim(phen_all_use)[1]){
        cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
        cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
        cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
        cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
        cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
        cleaned_data[k,6]<-as.numeric(st[k])
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
      colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","chronic","species","village")
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
      
      cleaned_data2<-as.data.frame(cleaned_data2)
      
      rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
      rownames(cleaned_data)<-rownames(phen_all_use)
      
      ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
      
      ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
      ind_sub2_ind<-which(is.na(ind_sub)==F)
      
      #tr3<-tr2[ind_sub2_ind]
      cleaned_data3<-cleaned_data2[ind_sub2,]
      
      sigma_phylo = 1
      n=length(ind_sub)
      
      true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
      
      #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
      cleaned_data4<-cleaned_data3
      
      cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
      cleaned_data3$village<-as.factor(cleaned_data3$village)
      
      cleaned_data4$chronic<-cleaned_data4$chronic#+true_phylo_effects[ind_sub2_ind]
      cleaned_data4$village<-as.factor(cleaned_data4$village)
      
      #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
      #                                      data = cleaned_data3,
      #                                      family = "gaussian")#use ordinal or gaussian
      
      
      #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
      #                   data = cleaned_data4,
      #                   family = "gaussian")#use ordinal or gaussian
      if((sum(cleaned_data4$chronic)>0)&(dim(cleaned_data4)[1]>30)){
      ms1<-lmer(chronic ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
                data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
      
      mss<-summary(ms1)
      
      coef(mss)
      
      effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
      pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
      
      #summary(mod)
      
      #summary(mod2)
      
      #mod2$B[8]
      #mod2$B.pvalue[8]
      #print(io)
      
      
      #effect.size<-mod2$B[8]
      #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
      
      esz[io,]<-effect.size
      pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
      }
    }
  }
  
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  #results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
  results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
  length(which(results_ftrs$Pvalue<0.05))
  length(which(results_ftrs$FDR<0.05))
  
  rownames(results_ftrs)<-sp_name2
  
  
  effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
  pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)
  
  
  
}

colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-c("GAD7-Q1","GAD7-Q2","GAD7-Q3","GAD7-Q4","GAD7-Q5","GAD7-Q6","GAD7-Q7","PHQ9-Q1","PHQ9-Q2","PHQ9-Q3","PHQ9-Q4","PHQ9-Q5","PHQ9-Q6","PHQ9-Q7","PHQ9-Q8","PHQ9-Q9")

lim<-dim(effect_size_phen)[2]

print("gad7 score and phq9 score")
#GAD7 score and PHQ9 score

for(kl in c(190,201)){
  
  esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
  pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 
  
  for(io in 1:dim(ind_sp)[1]){
    
    tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
    
    cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
    if(dim(cor_mat)[1]>15){
      sp<-as.data.frame(data_transformed[ind_sp2[io,],])
      
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
      colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","phy","species","village")
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
      
      cleaned_data2<-as.data.frame(cleaned_data2)
      
      rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
      rownames(cleaned_data)<-rownames(phen_all_use)
      
      ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
      
      ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
      ind_sub2_ind<-which(is.na(ind_sub)==F)
      
      #tr3<-tr2[ind_sub2_ind]
      cleaned_data3<-cleaned_data2[ind_sub2,]
      
      sigma_phylo = 1
      n=length(ind_sub)
      
      true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
      
      #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
      cleaned_data4<-cleaned_data3
      
      cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
      cleaned_data3$village<-as.factor(cleaned_data3$village)
      
      cleaned_data4$phy<-cleaned_data4$phy#+true_phylo_effects[ind_sub2_ind]
      cleaned_data4$village<-as.factor(cleaned_data4$village)
      
      #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
      #                                      data = cleaned_data3,
      #                                      family = "gaussian")#use ordinal or gaussian
      
      
      #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
      #                   data = cleaned_data4,
      #                   family = "gaussian")#use ordinal or gaussian
      if((sum(cleaned_data4$phy)>0)&(dim(cleaned_data4)[1]>30)){
      ms1<-lmer(phy ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
                data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
      
      mss<-summary(ms1)
      
      coef(mss)
      
      effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
      pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
      
      #summary(mod)
      
      #summary(mod2)
      
      #mod2$B[8]
      #mod2$B.pvalue[8]
      #print(io)
      
      
      #effect.size<-mod2$B[8]
      #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
      
      esz[io,]<-effect.size
      pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
      }
    }
  }
  
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  #results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
  results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
  length(which(results_ftrs$Pvalue<0.05))
  length(which(results_ftrs$FDR<0.05))
  
  rownames(results_ftrs)<-sp_name2
  
  
  effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
  pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)
  
}

colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-c("GAD7 score","PHQ9 score")



lim<-dim(effect_size_phen)[2]
print("animals")
#Animals

for(kl in 155:180){
  
  esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
  pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 
  
  for(io in 1:dim(ind_sp)[1]){
    
    tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
    
    cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
    if(dim(cor_mat)[1]>15){
      sp<-as.data.frame(data_transformed[ind_sp2[io,],])
      
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
      colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","chronic","species","village")
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
      
      cleaned_data2<-as.data.frame(cleaned_data2)
      
      rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
      rownames(cleaned_data)<-rownames(phen_all_use)
      
      ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
      
      ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
      ind_sub2_ind<-which(is.na(ind_sub)==F)
      
      #tr3<-tr2[ind_sub2_ind]
      cleaned_data3<-cleaned_data2[ind_sub2,]
      
      sigma_phylo = 1
      n=length(ind_sub)
      
      true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
      
      #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
      cleaned_data4<-cleaned_data3
      
      cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
      cleaned_data3$village<-as.factor(cleaned_data3$village)
      
      cleaned_data4$chronic<-cleaned_data4$chronic#+true_phylo_effects[ind_sub2_ind]
      cleaned_data4$village<-as.factor(cleaned_data4$village)
      
      #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
      #                                      data = cleaned_data3,
      #                                      family = "gaussian")#use ordinal or gaussian
      
      
      #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
      #                   data = cleaned_data4,
      #                   family = "gaussian")#use ordinal or gaussian
      if((sum(cleaned_data4$chronic)>1)&(dim(cleaned_data4)[1]>30)&(sum(cleaned_data4$chronic)<dim(cleaned_data4)[1])){
      ms1<-lmer(chronic ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
                data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
      
      mss<-summary(ms1)
      
      coef(mss)
      
      effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
      pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
      
      #summary(mod)
      
      #summary(mod2)
      
      #mod2$B[8]
      #mod2$B.pvalue[8]
      #print(io)
      
      
      #effect.size<-mod2$B[8]
      #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
      
      esz[io,]<-effect.size
      pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
      }
    }
  }
  
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  #results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
  results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
  length(which(results_ftrs$Pvalue<0.05))
  length(which(results_ftrs$FDR<0.05))
  
  rownames(results_ftrs)<-sp_name2
  
  
  effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
  pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)
  
  
  
}

colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-c("Cat","Dog","Parakeet","Rabbit","Horse","Mice","None pet","Cow","Goat","Pig","Chicken","Duck","Turkey","Sheep","Geese","None farm","Bat","Lizard","Monkey","Snake","Bird","Possum","Rat","Squirrel","Other wild","None wild")

lim<-dim(effect_size_phen)[2]

print("food")
#Food


for(kl in 75:93){
  
  esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
  pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 
  
  for(io in 1:dim(ind_sp)[1]){
    
    tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
    
    cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
    if(dim(cor_mat)[1]>15){
      sp<-as.data.frame(data_transformed[ind_sp2[io,],])
      
      cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
      for(k in 1:dim(phen_all_use)[1]){
        cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
        cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
        cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
        cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
        cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
        cleaned_data[k,6]<-as.numeric(st[k])
        if(is.na(phen_all_use[k,kl])==F){
          if(phen_all_use[k,kl]=="Never/rarely"){#Change to 1/50,1/10,4/7,1
            cleaned_data[k,7]<-1/50
          }else if(phen_all_use[k,kl]=="A few days per month"){
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
      colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","chronic","species","village")
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
      
      cleaned_data2<-as.data.frame(cleaned_data2)
      
      rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
      rownames(cleaned_data)<-rownames(phen_all_use)
      
      ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
      
      ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
      ind_sub2_ind<-which(is.na(ind_sub)==F)
      
      #tr3<-tr2[ind_sub2_ind]
      cleaned_data3<-cleaned_data2[ind_sub2,]
      
      sigma_phylo = 1
      n=length(ind_sub)
      
      true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
      
      #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
      cleaned_data4<-cleaned_data3
      
      cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
      cleaned_data3$village<-as.factor(cleaned_data3$village)
      
      cleaned_data4$chronic<-cleaned_data4$chronic#+true_phylo_effects[ind_sub2_ind]
      cleaned_data4$village<-as.factor(cleaned_data4$village)
      
      #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
      #                                      data = cleaned_data3,
      #                                      family = "gaussian")#use ordinal or gaussian
      
      
      #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
      #                   data = cleaned_data4,
      #                   family = "gaussian")#use ordinal or gaussian
      if((sum(cleaned_data4$chronic)>0)&(dim(cleaned_data4)[1]>30)&(length(unique(cleaned_data4$chronic))>1)){
      ms1<-lmer(chronic ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
                data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
      
      mss<-summary(ms1)
      
      coef(mss)
      
      effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
      pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
      
      #summary(mod)
      
      #summary(mod2)
      
      #mod2$B[8]
      #mod2$B.pvalue[8]
      #print(io)
      
      
      #effect.size<-mod2$B[8]
      #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
      
      esz[io,]<-effect.size
      pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
      }
    }
  }
  
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  #results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
  results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
  length(which(results_ftrs$Pvalue<0.05))
  length(which(results_ftrs$FDR<0.05))
  
  rownames(results_ftrs)<-sp_name2
  
  
  effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
  pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)
  
  
  
}

colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-c("Beans","Tortillas","Rice","Bread","Milk","Yogurt","Cream/butter","Cheese","Eggs","Vegetables","Fruits","Natural juice","Chicken","Beef/Pork","Ham/sausages/hotdog","Fish","Soda","Fruit juice","Chips")

lim<-dim(effect_size_phen)[2]


# Travel

esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 

for(io in 1:dim(ind_sp)[1]){
  
  tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
  
  cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
  if(dim(cor_mat)[1]>15){
    sp<-as.data.frame(data_transformed[ind_sp2[io,],])
    
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(st[k])
      if(is.na(phen_all_use[k,181])==F){
        if(phen_all_use[k,181]=="Rarely/Never"){
          cleaned_data[k,7]<-1/100
        }else if(phen_all_use[k,181]=="At least once a month"){
          cleaned_data[k,7]<-1/30
        }else if(phen_all_use[k,181]=="At least once a week"){
          cleaned_data[k,7]<-1/7
        }else if(phen_all_use[k,181]=="Every day"){
          cleaned_data[k,7]<-1
        }
      }
      if(is.na(phen_all_use[k,181])==T){
        cleaned_data[k,7]<-NA
      }
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,8]<-as.numeric(sp[k,1])
      cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
      
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","chronic","species","village")
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
    
    cleaned_data2<-as.data.frame(cleaned_data2)
    
    rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
    rownames(cleaned_data)<-rownames(phen_all_use)
    
    ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
    
    ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
    ind_sub2_ind<-which(is.na(ind_sub)==F)
    
    #tr3<-tr2[ind_sub2_ind]
    cleaned_data3<-cleaned_data2[ind_sub2,]
    
    sigma_phylo = 1
    n=length(ind_sub)
    
    true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
    
    #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
    cleaned_data4<-cleaned_data3
    
    cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
    cleaned_data3$village<-as.factor(cleaned_data3$village)
    
    cleaned_data4$chronic<-cleaned_data4$chronic#+true_phylo_effects[ind_sub2_ind]
    cleaned_data4$village<-as.factor(cleaned_data4$village)
    
    #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
    #                                      data = cleaned_data3,
    #                                      family = "gaussian")#use ordinal or gaussian
    
    
    #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
    #                   data = cleaned_data4,
    #                   family = "gaussian")#use ordinal or gaussian
    if((sum(cleaned_data4$chronic)>0)&(dim(cleaned_data4)[1]>30)){
    ms1<-lmer(chronic ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
              data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
    
    mss<-summary(ms1)
    
    coef(mss)
    
    effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
    pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
    
    #summary(mod)
    
    #summary(mod2)
    
    #mod2$B[8]
    #mod2$B.pvalue[8]
    #print(io)
    
    
    #effect.size<-mod2$B[8]
    #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
    
    esz[io,]<-effect.size
    pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    }
  }
}


results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")
#results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
length(which(results_ftrs$Pvalue<0.05))
length(which(results_ftrs$FDR<0.05))

rownames(results_ftrs)<-sp_name2


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)



colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Travel")

lim<-dim(effect_size_phen)[2]


# Monthly expenditure

kl<-182
esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 

for(io in 1:dim(ind_sp)[1]){
  
  tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
  
  cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
  if(dim(cor_mat)[1]>15){
    sp<-as.data.frame(data_transformed[ind_sp2[io,],])
    
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
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","phy","species","village")
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
    
    cleaned_data2<-as.data.frame(cleaned_data2)
    
    rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
    rownames(cleaned_data)<-rownames(phen_all_use)
    
    ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
    
    ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
    ind_sub2_ind<-which(is.na(ind_sub)==F)
    
    #tr3<-tr2[ind_sub2_ind]
    cleaned_data3<-cleaned_data2[ind_sub2,]
    
    sigma_phylo = 1
    n=length(ind_sub)
    
    true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
    
    #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
    cleaned_data4<-cleaned_data3
    
    cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
    cleaned_data3$village<-as.factor(cleaned_data3$village)
    
    cleaned_data4$phy<-cleaned_data4$phy#+true_phylo_effects[ind_sub2_ind]
    cleaned_data4$village<-as.factor(cleaned_data4$village)
    
    #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
    #                                      data = cleaned_data3,
    #                                      family = "gaussian")#use ordinal or gaussian
    
    
    #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
    #                   data = cleaned_data4,
    #                   family = "gaussian")#use ordinal or gaussian
    if((sum(cleaned_data4$phy)>0)&(dim(cleaned_data4)[1]>30)){
    ms1<-lmer(phy ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
              data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
    
    mss<-summary(ms1)
    
    coef(mss)
    
    effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
    pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
    
    #summary(mod)
    
    #summary(mod2)
    
    #mod2$B[8]
    #mod2$B.pvalue[8]
    #print(io)
    
    
    #effect.size<-mod2$B[8]
    #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
    
    esz[io,]<-effect.size
    pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    }
  }
}


results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")
#results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
length(which(results_ftrs$Pvalue<0.05))
length(which(results_ftrs$FDR<0.05))

rownames(results_ftrs)<-sp_name2


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)

colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Monthly expenditure")

lim<-dim(effect_size_phen)[2]


## Alcohol frequency

kl<-8
esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 

for(io in 1:dim(ind_sp)[1]){
  
  tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
  
  cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
  if(dim(cor_mat)[1]>15){
    sp<-as.data.frame(data_transformed[ind_sp2[io,],])
    
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(st[k])
      if(is.na(phen_all_rct_use[k,kl])==F){
        if(phen_all_rct_use[k,kl]=="1 or 2"){
          cleaned_data[k,7]<-1.5
        }else if(phen_all_rct_use[k,kl]=="3 or 4"){
          cleaned_data[k,7]<-3.5
        }else if(phen_all_rct_use[k,kl]=="5 or 6"){
          cleaned_data[k,7]<-5.5
        }else if(phen_all_rct_use[k,kl]=="7 to 9"){
          cleaned_data[k,7]<-8
        }else if(phen_all_rct_use[k,kl]=="10 or more"){
          cleaned_data[k,7]<-12
        }
      }
      if(is.na(phen_all_rct_use[k,kl])==T){
        cleaned_data[k,7]<-0
      }
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,8]<-as.numeric(sp[k,1])
      cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
      
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","chronic","species","village")
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
    
    cleaned_data2<-as.data.frame(cleaned_data2)
    
    rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
    rownames(cleaned_data)<-rownames(phen_all_use)
    
    ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
    
    ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
    ind_sub2_ind<-which(is.na(ind_sub)==F)
    
    #tr3<-tr2[ind_sub2_ind]
    cleaned_data3<-cleaned_data2[ind_sub2,]
    
    sigma_phylo = 1
    n=length(ind_sub)
    
    true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
    
    #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
    cleaned_data4<-cleaned_data3
    
    cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
    cleaned_data3$village<-as.factor(cleaned_data3$village)
    
    cleaned_data4$chronic<-cleaned_data4$chronic#+true_phylo_effects[ind_sub2_ind]
    cleaned_data4$village<-as.factor(cleaned_data4$village)
    
    #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
    #                                      data = cleaned_data3,
    #                                      family = "gaussian")#use ordinal or gaussian
    
    
    #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
    #                   data = cleaned_data4,
    #                   family = "gaussian")#use ordinal or gaussian
    if((length(which(cleaned_data4$chronic>0))>1)&(dim(cleaned_data4)[1]>30)){
    ms1<-lmer(chronic ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
              data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
    
    mss<-summary(ms1)
    
    coef(mss)
    
    effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
    pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
    
    #summary(mod)
    
    #summary(mod2)
    
    #mod2$B[8]
    #mod2$B.pvalue[8]
    #print(io)
    
    
    #effect.size<-mod2$B[8]
    #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
    
    esz[io,]<-effect.size
    pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    }
  }
}


results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")
#results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
length(which(results_ftrs$Pvalue<0.05))
length(which(results_ftrs$FDR<0.05))

rownames(results_ftrs)<-sp_name2


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)



colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Alcohol frequency")

lim<-dim(effect_size_phen)[2]


## Alcohol daily

esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 
kl<-9
for(io in 1:dim(ind_sp)[1]){
  
  tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
  
  cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
  if(dim(cor_mat)[1]>15){
    sp<-as.data.frame(data_transformed[ind_sp2[io,],])
    
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(st[k])
      if(is.na(phen_all_rct_use[k,kl])==F){
        if(phen_all_rct_use[k,kl]=="Never"){
          cleaned_data[k,7]<-0
        }else if(phen_all_rct_use[k,kl]=="Less than monthly"){
          cleaned_data[k,7]<-1/50
        }else if(phen_all_rct_use[k,kl]=="Monthly"){
          cleaned_data[k,7]<-1/30
        }else if(phen_all_rct_use[k,kl]=="Weekly"){
          cleaned_data[k,7]<-1/7
        }else if(phen_all_rct_use[k,kl]=="Daily or almost daily"){
          cleaned_data[k,7]<-1
        }
      }
      if(is.na(phen_all_rct_use[k,kl])==T){
        cleaned_data[k,7]<-0
      }
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,8]<-as.numeric(sp[k,1])
      cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
      
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","chronic","species","village")
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
    
    cleaned_data2<-as.data.frame(cleaned_data2)
    
    rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
    rownames(cleaned_data)<-rownames(phen_all_use)
    
    ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
    
    ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
    ind_sub2_ind<-which(is.na(ind_sub)==F)
    
    #tr3<-tr2[ind_sub2_ind]
    cleaned_data3<-cleaned_data2[ind_sub2,]
    
    sigma_phylo = 1
    n=length(ind_sub)
    
    true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
    
    #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
    cleaned_data4<-cleaned_data3
    
    cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
    cleaned_data3$village<-as.factor(cleaned_data3$village)
    
    cleaned_data4$chronic<-cleaned_data4$chronic#+true_phylo_effects[ind_sub2_ind]
    cleaned_data4$village<-as.factor(cleaned_data4$village)
    
    #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
    #                                      data = cleaned_data3,
    #                                      family = "gaussian")#use ordinal or gaussian
    
    
    #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
    #                   data = cleaned_data4,
    #                   family = "gaussian")#use ordinal or gaussian
    if((length(which(cleaned_data4$chronic>0))>1)&(dim(cleaned_data4)[1]>30)){
    ms1<-lmer(chronic ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
              data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
    
    mss<-summary(ms1)
    
    coef(mss)
    
    effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
    pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
    
    #summary(mod)
    
    #summary(mod2)
    
    #mod2$B[8]
    #mod2$B.pvalue[8]
    #print(io)
    
    
    #effect.size<-mod2$B[8]
    #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
    
    esz[io,]<-effect.size
    pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    }
  }
}


results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")
#results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
length(which(results_ftrs$Pvalue<0.05))
length(which(results_ftrs$FDR<0.05))

rownames(results_ftrs)<-sp_name2


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)



colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Alcohol daily")

lim<-dim(effect_size_phen)[2]




##Cigarette use

esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 
kl<-10
for(io in 1:dim(ind_sp)[1]){
  
  tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
  
  cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
  if(dim(cor_mat)[1]>15){
    sp<-as.data.frame(data_transformed[ind_sp2[io,],])
    
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(st[k])
      if(is.na(phen_all_rct_use[k,kl])==F){
        if(phen_all_rct_use[k,kl]=="No"){
          cleaned_data[k,7]<-0
        }
        if(phen_all_rct_use[k,kl]=="Yes"){
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
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","chronic","species","village")
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
    
    cleaned_data2<-as.data.frame(cleaned_data2)
    
    rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
    rownames(cleaned_data)<-rownames(phen_all_use)
    
    ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
    
    ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
    ind_sub2_ind<-which(is.na(ind_sub)==F)
    
    #tr3<-tr2[ind_sub2_ind]
    cleaned_data3<-cleaned_data2[ind_sub2,]
    
    sigma_phylo = 1
    n=length(ind_sub)
    
    true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
    
    #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
    cleaned_data4<-cleaned_data3
    
    cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
    cleaned_data3$village<-as.factor(cleaned_data3$village)
    
    cleaned_data4$chronic<-cleaned_data4$chronic#+true_phylo_effects[ind_sub2_ind]
    cleaned_data4$village<-as.factor(cleaned_data4$village)
    
    #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
    #                                      data = cleaned_data3,
    #                                      family = "gaussian")#use ordinal or gaussian
    
    
    #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
    #                   data = cleaned_data4,
    #                   family = "gaussian")#use ordinal or gaussian
    if((length(which(cleaned_data4$chronic>0))>1)&(dim(cleaned_data4)[1]>30)){
    ms1<-lmer(chronic ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
              data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
    
    mss<-summary(ms1)
    
    coef(mss)
    
    effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
    pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
    
    #summary(mod)
    
    #summary(mod2)
    
    #mod2$B[8]
    #mod2$B.pvalue[8]
    #print(io)
    
    
    #effect.size<-mod2$B[8]
    #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
    
    esz[io,]<-effect.size
    pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    }
  }
}


results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")
#results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
length(which(results_ftrs$Pvalue<0.05))
length(which(results_ftrs$FDR<0.05))

rownames(results_ftrs)<-sp_name2


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)



colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Cigarette use")

lim<-dim(effect_size_phen)[2]


##Cigarette frequency

esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 
kl<-11
for(io in 1:dim(ind_sp)[1]){
  
  tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
  
  cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
  if(dim(cor_mat)[1]>15){
    sp<-as.data.frame(data_transformed[ind_sp2[io,],])
    
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(st[k])
      if(is.na(phen_all_rct_use[k,kl])==F){
        cleaned_data[k,7]<-as.numeric(phen_all_rct_use[k,kl])
      }
      if(is.na(phen_all_rct_use[k,kl])==T){
        cleaned_data[k,7]<-0
      }
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,8]<-as.numeric(sp[k,1])
      cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
      
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","chronic","species","village")
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
    
    cleaned_data2<-as.data.frame(cleaned_data2)
    
    rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
    rownames(cleaned_data)<-rownames(phen_all_use)
    
    ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
    
    ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
    ind_sub2_ind<-which(is.na(ind_sub)==F)
    
    #tr3<-tr2[ind_sub2_ind]
    cleaned_data3<-cleaned_data2[ind_sub2,]
    
    sigma_phylo = 1
    n=length(ind_sub)
    
    true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
    
    #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
    cleaned_data4<-cleaned_data3
    
    cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
    cleaned_data3$village<-as.factor(cleaned_data3$village)
    
    cleaned_data4$chronic<-cleaned_data4$chronic#+true_phylo_effects[ind_sub2_ind]
    cleaned_data4$village<-as.factor(cleaned_data4$village)
    
    #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
    #                                      data = cleaned_data3,
    #                                      family = "gaussian")#use ordinal or gaussian
    
    
    #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
    #                   data = cleaned_data4,
    #                   family = "gaussian")#use ordinal or gaussian
    if((length(which(cleaned_data4$chronic>0))>1)&(dim(cleaned_data4)[1]>30)){
    ms1<-lmer(chronic ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
              data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
    
    mss<-summary(ms1)
    
    coef(mss)
    
    effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
    pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
    
    #summary(mod)
    
    #summary(mod2)
    
    #mod2$B[8]
    #mod2$B.pvalue[8]
    #print(io)
    
    
    #effect.size<-mod2$B[8]
    #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
    
    esz[io,]<-effect.size
    pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    }
  }
}


results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")
#results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
length(which(results_ftrs$Pvalue<0.05))
length(which(results_ftrs$FDR<0.05))

rownames(results_ftrs)<-sp_name2


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)



colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Cigarette use")

lim<-dim(effect_size_phen)[2]



##Altruism and risk taking

for(kl in c(208,364)){
  esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
  pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 
  
  for(io in 1:dim(ind_sp)[1]){
    
    tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
    
    cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
    if(dim(cor_mat)[1]>15){
      sp<-as.data.frame(data_transformed[ind_sp2[io,],])
      
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
      colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","phy","species","village")
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
      
      cleaned_data2<-as.data.frame(cleaned_data2)
      
      rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
      rownames(cleaned_data)<-rownames(phen_all_use)
      
      ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
      
      ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
      ind_sub2_ind<-which(is.na(ind_sub)==F)
      
      #tr3<-tr2[ind_sub2_ind]
      cleaned_data3<-cleaned_data2[ind_sub2,]
      
      sigma_phylo = 1
      n=length(ind_sub)
      
      true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
      
      #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
      cleaned_data4<-cleaned_data3
      
      cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
      cleaned_data3$village<-as.factor(cleaned_data3$village)
      
      cleaned_data4$phy<-cleaned_data4$phy#+true_phylo_effects[ind_sub2_ind]
      cleaned_data4$village<-as.factor(cleaned_data4$village)
      
      #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
      #                                      data = cleaned_data3,
      #                                      family = "gaussian")#use ordinal or gaussian
      
      
      #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
      #                   data = cleaned_data4,
      #                   family = "gaussian")#use ordinal or gaussian
      if((sum(cleaned_data4$phy)>0)&(dim(cleaned_data4)[1]>30)){
      ms1<-lmer(phy ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
                data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
      
      mss<-summary(ms1)
      
      coef(mss)
      
      effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
      pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
      
      #summary(mod)
      
      #summary(mod2)
      
      #mod2$B[8]
      #mod2$B.pvalue[8]
      #print(io)
      
      
      #effect.size<-mod2$B[8]
      #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
      
      esz[io,]<-effect.size
      pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
      }
    }
  }
  
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  #results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
  results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
  length(which(results_ftrs$Pvalue<0.05))
  length(which(results_ftrs$FDR<0.05))
  
  rownames(results_ftrs)<-sp_name2
  
  
  effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
  pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)
}
colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-c("Altruism","Risk taking")

lim<-dim(effect_size_phen)[2]

##Education
esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 
kl<-5
for(io in 1:dim(ind_sp)[1]){
  
  tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
  
  cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
  if(dim(cor_mat)[1]>15){
    sp<-as.data.frame(data_transformed[ind_sp2[io,],])
    
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(st[k])
      if(is.na(phen_all_rct_use[k,kl])==F){
        if(phen_all_rct_use[k,kl]=="Have not completed any type of school"){
          cleaned_data[k,7]<-0
        }else if(phen_all_rct_use[k,kl]=="1st grade"){
          cleaned_data[k,7]<-1
        }else if(phen_all_rct_use[k,kl]=="2nd grade"){
          cleaned_data[k,7]<-2
        }else if(phen_all_rct_use[k,kl]=="3rd grade"){
          cleaned_data[k,7]<-3
        }else if(phen_all_rct_use[k,kl]=="4th grade"){
          cleaned_data[k,7]<-4
        }else if(phen_all_rct_use[k,kl]=="5th grade"){
          cleaned_data[k,7]<-5
        }else if(phen_all_rct_use[k,kl]=="6th grade"){
          cleaned_data[k,7]<-6
        }else if(phen_all_rct_use[k,kl]=="Some secondary"){
          cleaned_data[k,7]<-7
        }else if(phen_all_rct_use[k,kl]=="More than secondary"){
          cleaned_data[k,7]<-8
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
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","chronic","species","village")
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
    
    cleaned_data2<-as.data.frame(cleaned_data2)
    
    rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
    rownames(cleaned_data)<-rownames(phen_all_use)
    
    ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
    
    ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
    ind_sub2_ind<-which(is.na(ind_sub)==F)
    
    #tr3<-tr2[ind_sub2_ind]
    cleaned_data3<-cleaned_data2[ind_sub2,]
    
    sigma_phylo = 1
    n=length(ind_sub)
    
    true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
    
    #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
    cleaned_data4<-cleaned_data3
    
    cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
    cleaned_data3$village<-as.factor(cleaned_data3$village)
    
    cleaned_data4$chronic<-cleaned_data4$chronic#+true_phylo_effects[ind_sub2_ind]
    cleaned_data4$village<-as.factor(cleaned_data4$village)
    
    #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
    #                                      data = cleaned_data3,
    #                                      family = "gaussian")#use ordinal or gaussian
    
    
    #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
    #                   data = cleaned_data4,
    #                   family = "gaussian")#use ordinal or gaussian
    if((sum(cleaned_data4$chronic)>0)&(dim(cleaned_data4)[1]>30)){
    ms1<-lmer(chronic ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
              data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
    
    mss<-summary(ms1)
    
    coef(mss)
    
    effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
    pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
    
    #summary(mod)
    
    #summary(mod2)
    
    #mod2$B[8]
    #mod2$B.pvalue[8]
    #print(io)
    
    
    #effect.size<-mod2$B[8]
    #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
    
    esz[io,]<-effect.size
    pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    }
  }
}


results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")
#results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
length(which(results_ftrs$Pvalue<0.05))
length(which(results_ftrs$FDR<0.05))

rownames(results_ftrs)<-sp_name2


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)



colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Education")

lim<-dim(effect_size_phen)[2]

##Social network factors

for(kl in c(203:207,360:363,365:369)){
  esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
  pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 
  
  for(io in 1:dim(ind_sp)[1]){
    
    tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
    
    cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
    if(dim(cor_mat)[1]>15){
      sp<-as.data.frame(data_transformed[ind_sp2[io,],])
      
      cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
      for(k in 1:dim(phen_all_use)[1]){
        cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
        cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
        cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
        cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
        cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
        cleaned_data[k,6]<-as.numeric(st[k])
        if(is.na(phen_all_use[k,kl])==F){
          cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
        }else{
          cleaned_data[k,7]<-NA
        }
        #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
        #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
        cleaned_data[k,8]<-as.numeric(sp[k,1])
        cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
        
      }
      cleaned_data<-as.data.frame(cleaned_data)
      colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","phy","species","village")
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
      
      cleaned_data2<-as.data.frame(cleaned_data2)
      
      rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
      rownames(cleaned_data)<-rownames(phen_all_use)
      
      ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
      
      ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
      ind_sub2_ind<-which(is.na(ind_sub)==F)
      
      #tr3<-tr2[ind_sub2_ind]
      cleaned_data3<-cleaned_data2[ind_sub2,]
      
      sigma_phylo = 1
      n=length(ind_sub)
      
      true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
      
      #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
      cleaned_data4<-cleaned_data3
      
      cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
      cleaned_data3$village<-as.factor(cleaned_data3$village)
      
      cleaned_data4$phy<-cleaned_data4$phy#+true_phylo_effects[ind_sub2_ind]
      cleaned_data4$village<-as.factor(cleaned_data4$village)
      
      #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
      #                                      data = cleaned_data3,
      #                                      family = "gaussian")#use ordinal or gaussian
      
      
      #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
      #                   data = cleaned_data4,
      #                   family = "gaussian")#use ordinal or gaussian
      if((sum(cleaned_data4$phy)>0)&(dim(cleaned_data4)[1]>30)){
      ms1<-lmer(phy ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
                data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
      
      mss<-summary(ms1)
      
      coef(mss)
      
      effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
      pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
      
      #summary(mod)
      
      #summary(mod2)
      
      #mod2$B[8]
      #mod2$B.pvalue[8]
      #print(io)
      
      
      #effect.size<-mod2$B[8]
      #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
      
      esz[io,]<-effect.size
      pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
      }
    }
  }
  
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  #results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
  results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
  length(which(results_ftrs$Pvalue<0.05))
  length(which(results_ftrs$FDR<0.05))
  
  rownames(results_ftrs)<-sp_name2
  
  
  effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
  pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)
}
#colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-colnames(phen_all_use)[c(203:207,360:363,365:369,228:230)]
#Run 227:230 again
colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-colnames(phen_all_use)[c(203:207,360:363,365:369)]
lim<-dim(effect_size_phen)[2]





#Washing hands


esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 
kl<-12
for(io in 1:dim(ind_sp)[1]){
  
  tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
  
  cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
  if(dim(cor_mat)[1]>15){
    sp<-as.data.frame(data_transformed[ind_sp2[io,],])
    
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(st[k])
      if(is.na(phen_all_rct_use[k,kl])==F){
        cleaned_data[k,7]<-1
      }
      if(is.na(phen_all_rct_use[k,kl])==T){
        cleaned_data[k,7]<-0
      }
      #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
      #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
      cleaned_data[k,8]<-as.numeric(sp[k,1])
      cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
      
    }
    cleaned_data<-as.data.frame(cleaned_data)
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","chronic","species","village")
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
    
    cleaned_data2<-as.data.frame(cleaned_data2)
    
    rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
    rownames(cleaned_data)<-rownames(phen_all_use)
    
    ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
    
    ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
    ind_sub2_ind<-which(is.na(ind_sub)==F)
    
    #tr3<-tr2[ind_sub2_ind]
    cleaned_data3<-cleaned_data2[ind_sub2,]
    
    sigma_phylo = 1
    n=length(ind_sub)
    
    true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
    
    #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
    cleaned_data4<-cleaned_data3
    
    cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
    cleaned_data3$village<-as.factor(cleaned_data3$village)
    
    cleaned_data4$chronic<-cleaned_data4$chronic#+true_phylo_effects[ind_sub2_ind]
    cleaned_data4$village<-as.factor(cleaned_data4$village)
    
    #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
    #                                      data = cleaned_data3,
    #                                      family = "gaussian")#use ordinal or gaussian
    
    
    #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
    #                   data = cleaned_data4,
    #                   family = "gaussian")#use ordinal or gaussian
    if((length(which(cleaned_data4$chronic>0))>0)&(dim(cleaned_data4)[1]>30)&(length(which(cleaned_data4$chronic>0))<dim(cleaned_data4)[1])){
    ms1<-lmer(chronic ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
              data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
    
    mss<-summary(ms1)
    
    coef(mss)
    
    effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
    pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
    
    #summary(mod)
    
    #summary(mod2)
    
    #mod2$B[8]
    #mod2$B.pvalue[8]
    #print(io)
    
    
    #effect.size<-mod2$B[8]
    #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
    
    esz[io,]<-effect.size
    pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    }
  }
}


results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")
#results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
length(which(results_ftrs$Pvalue<0.05))
length(which(results_ftrs$FDR<0.05))

rownames(results_ftrs)<-sp_name2


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)




colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Washing hands")

lim<-dim(effect_size_phen)[2]


##HH essentials
for(kl in 310:316){
  
  esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
  pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 
  
  for(io in 1:dim(ind_sp)[1]){
    
    tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
    
    cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
    if(dim(cor_mat)[1]>15){
      sp<-as.data.frame(data_transformed[ind_sp2[io,],])
      
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
      colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","chronic","species","village")
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
      
      cleaned_data2<-as.data.frame(cleaned_data2)
      
      rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
      rownames(cleaned_data)<-rownames(phen_all_use)
      
      ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
      
      ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
      ind_sub2_ind<-which(is.na(ind_sub)==F)
      
      #tr3<-tr2[ind_sub2_ind]
      cleaned_data3<-cleaned_data2[ind_sub2,]
      
      sigma_phylo = 1
      n=length(ind_sub)
      
      true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
      
      #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
      cleaned_data4<-cleaned_data3
      
      cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
      cleaned_data3$village<-as.factor(cleaned_data3$village)
      
      cleaned_data4$chronic<-cleaned_data4$chronic#+true_phylo_effects[ind_sub2_ind]
      cleaned_data4$village<-as.factor(cleaned_data4$village)
      
      #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
      #                                      data = cleaned_data3,
      #                                      family = "gaussian")#use ordinal or gaussian
      
      
      #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
      #                   data = cleaned_data4,
      #                   family = "gaussian")#use ordinal or gaussian
      if((sum(cleaned_data4$chronic)>0)&(dim(cleaned_data4)[1]>30)){
      ms1<-lmer(chronic ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
                data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
      
      mss<-summary(ms1)
      
      coef(mss)
      
      effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
      pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
      
      #summary(mod)
      
      #summary(mod2)
      
      #mod2$B[8]
      #mod2$B.pvalue[8]
      #print(io)
      
      
      #effect.size<-mod2$B[8]
      #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
      
      esz[io,]<-effect.size
      pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
      }
    }
  }
  
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  #results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
  results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
  length(which(results_ftrs$Pvalue<0.05))
  length(which(results_ftrs$FDR<0.05))
  
  rownames(results_ftrs)<-sp_name2
  
  
  effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
  pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)
  
  
  
}

colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-c("Electricity","Radio","Television","Cell phone","No cellphone","Refrigerator","None")

lim<-dim(effect_size_phen)[2]


##HH remaining

for(kl in c(306:309,317:359)){#306:309,317:359
  
  esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
  pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 
  
  for(io in 1:dim(ind_sp)[1]){
    
    tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
    
    cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
    if(dim(cor_mat)[1]>15){
      sp<-as.data.frame(data_transformed[ind_sp2[io,],])
      
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
      colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","chronic","species","village")
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
      
      cleaned_data2<-as.data.frame(cleaned_data2)
      
      rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
      rownames(cleaned_data)<-rownames(phen_all_use)
      
      ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
      
      ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
      ind_sub2_ind<-which(is.na(ind_sub)==F)
      
      #tr3<-tr2[ind_sub2_ind]
      cleaned_data3<-cleaned_data2[ind_sub2,]
      
      sigma_phylo = 1
      n=length(ind_sub)
      
      true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
      
      #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
      cleaned_data4<-cleaned_data3
      
      cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
      cleaned_data3$village<-as.factor(cleaned_data3$village)
      
      cleaned_data4$chronic<-cleaned_data4$chronic##+true_phylo_effects[ind_sub2_ind]
      cleaned_data4$village<-as.factor(cleaned_data4$village)
      
      #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
      #                                      data = cleaned_data3,
      #                                      family = "gaussian")#use ordinal or gaussian
      
      
      #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
      #                   data = cleaned_data4,
      #                   family = "gaussian")#use ordinal or gaussian
      if((length(which(cleaned_data4$chronic>0))>2)&(dim(cleaned_data4)[1]>40)&(length(which(cleaned_data4$chronic>0))<(dim(cleaned_data4)[1]-2))){
      ms1<-lmer(chronic ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
                data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
      
      mss<-summary(ms1)
      
      coef(mss)
      
      effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
      pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
      
      #summary(mod)
      
      #summary(mod2)
      
      #mod2$B[8]
      #mod2$B.pvalue[8]
      #print(io)
      
      
      #effect.size<-mod2$B[8]
      #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
      
      esz[io,]<-effect.size
      pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
      }
    }
  }
  
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  #results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
  results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
  length(which(results_ftrs$Pvalue<0.05))
  length(which(results_ftrs$FDR<0.05))
  
  rownames(results_ftrs)<-sp_name2
  
  
  effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
  pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)
  
  
  
}

colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-colnames(phen_all_use)[c(306:309,317:359)]

lim<-dim(effect_size_phen)[2]

##Partners phenotype!

for(kl in c(228:230)){
  esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
  pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 
  
  for(io in 1:dim(ind_sp)[1]){
    
    tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
    
    cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
    if(dim(cor_mat)[1]>15){
      sp<-as.data.frame(data_transformed[ind_sp2[io,],])
      
      cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
      for(k in 1:dim(phen_all_use)[1]){
        cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
        cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
        cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
        cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
        cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
        cleaned_data[k,6]<-as.numeric(st[k])
        if(is.na(phen_all_use[k,kl])==F){
          cleaned_data[k,7]<-as.numeric(as.character(phen_all_use[k,kl]))
        }else{
          cleaned_data[k,7]<-NA
        }
        #cleaned_data[k,8]<-ifelse(phen_all_use[k,191]=="moderate",1,0)  #Phenotype
        #cleaned_data[k,9]<-ifelse(phen_all_use[k,191]=="severe",1,0)  #Phenotype
        cleaned_data[k,8]<-as.numeric(sp[k,1])
        cleaned_data[k,9]<-as.numeric(phen_all_use[k,11])
        
      }
      cleaned_data<-as.data.frame(cleaned_data)
      colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","phy","species","village")
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
      
      cleaned_data2<-as.data.frame(cleaned_data2)
      
      rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
      rownames(cleaned_data)<-rownames(phen_all_use)
      
      ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
      
      ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
      ind_sub2_ind<-which(is.na(ind_sub)==F)
      
      #tr3<-tr2[ind_sub2_ind]
      cleaned_data3<-cleaned_data2[ind_sub2,]
      
      sigma_phylo = 1
      n=length(ind_sub)
      
      true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
      
      #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
      cleaned_data4<-cleaned_data3
      
      cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
      cleaned_data3$village<-as.factor(cleaned_data3$village)
      
      cleaned_data4$phy<-cleaned_data4$phy##+true_phylo_effects[ind_sub2_ind]
      cleaned_data4$village<-as.factor(cleaned_data4$village)
      
      if((length(which(cleaned_data4$phy>0))>2)&(dim(cleaned_data4)[1]>30)){
        #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
        #                                      data = cleaned_data3,
        #                                      family = "gaussian")#use ordinal or gaussian
        
        
        #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
        #                   data = cleaned_data4,
        #                   family = "gaussian")#use ordinal or gaussian
        
        ms1<-lmer(phy ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
                  data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
        
        mss<-summary(ms1)
        
        coef(mss)
        
        effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
        pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
        
        #summary(mod)
        
        #summary(mod2)
        
        #mod2$B[8]
        #mod2$B.pvalue[8]
        #print(io)
        
        
        #effect.size<-mod2$B[8]
        #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
        
        esz[io,]<-effect.size
        pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
      }
    }
  }
  
  
  results_ftrs<-cbind(esz,pv_sp)
  results_ftrs<-as.data.frame(results_ftrs)
  colnames(results_ftrs)<-c("Effect_size","Pvalue")
  #results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
  results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
  length(which(results_ftrs$Pvalue<0.05))
  length(which(results_ftrs$FDR<0.05))
  
  rownames(results_ftrs)<-sp_name2
  
  
  effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
  pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
  fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)
}
#colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-colnames(phen_all_use)[c(203:207,360:363,365:369,228:230)]
#Run 227:230 again
colnames(effect_size_phen)[(lim+1):dim(effect_size_phen)[2]]<-colnames(phen_all_use)[c(228:230)]
lim<-dim(effect_size_phen)[2]

#Living with partners

esz<-array(NA,dim=c(dim(ind_sp2)[1],1))
pv_sp<-array(NA,dim=c(dim(ind_sp2)[1],1)) #p value 
kl<-227
for(io in 1:dim(ind_sp)[1]){
  
  tr2 <- ape::read.tree(file_list[ind_sp[io,1]])
  
  cor_mat = ape::vcv.phylo(tr2, corr = TRUE)
  if(dim(cor_mat)[1]>15){
    sp<-as.data.frame(data_transformed[ind_sp2[io,],])
    
    cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],9))
    for(k in 1:dim(phen_all_use)[1]){
      cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
      cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
      cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
      cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
      cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
      cleaned_data[k,6]<-as.numeric(st[k])
      if(is.na(phen_all_use[k,kl])==F){
        if(phen_all_use[k,kl]=="No"){
          cleaned_data[k,7]<-0
        }
        if(phen_all_use[k,kl]=="Yes"){
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
    colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","chronic","species","village")
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
    
    cleaned_data2<-as.data.frame(cleaned_data2)
    
    rownames(cleaned_data2)<-rownames(phen_all_use[ind_temp,])
    rownames(cleaned_data)<-rownames(phen_all_use)
    
    ind_sub<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data2))#was cleaned_data2
    
    ind_sub2<-ind_sub[which(is.na(ind_sub)==F)]
    ind_sub2_ind<-which(is.na(ind_sub)==F)
    
    #tr3<-tr2[ind_sub2_ind]
    cleaned_data3<-cleaned_data2[ind_sub2,]
    
    sigma_phylo = 1
    n=length(ind_sub)
    
    true_phylo_effects = sigma_phylo * MASS::mvrnorm(1, mu = rep(0, n), Sigma = cor_mat)
    
    #cleaned_data3$a1c<-cleaned_data3$a1c+true_phylo_effects
    cleaned_data4<-cleaned_data3
    
    cleaned_data3<-cbind(cleaned_data3,"phylo"=true_phylo_effects[ind_sub2_ind])
    cleaned_data3$village<-as.factor(cleaned_data3$village)
    
    cleaned_data4$chronic<-cleaned_data4$chronic##+true_phylo_effects[ind_sub2_ind]
    cleaned_data4$village<-as.factor(cleaned_data4$village)
    
    if((length(which(cleaned_data4$chronic>0))>2)&(dim(cleaned_data4)[1]>40)&(length(which(cleaned_data4$chronic>0))<(dim(cleaned_data4)[1]-2))){
      #mod <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species + (1 | phylo) +(1|village), 
      #                                      data = cleaned_data3,
      #                                      family = "gaussian")#use ordinal or gaussian
      
      
      #mod2 <- phyr::pglmm(a1c ~ age+sex+batch_effect+dna_conc+BMI+bristol+species +(1|village), 
      #                   data = cleaned_data4,
      #                   family = "gaussian")#use ordinal or gaussian
      
      ms1<-lmer(chronic ~ age+sex+batch_effect+dna_conc+BMI+bristol+species  +(1|village), 
                data = cleaned_data4)#Cannot use phylo because you have one value for every data point --- so no random effect
      
      mss<-summary(ms1)
      
      coef(mss)
      
      effect.size<-coef(mss)[which(rownames(coef(mss))=="species"),1]
      pval2<-coef(mss)[which(rownames(coef(mss))=="species"),5]#Pvalue of phenotype
      
      #summary(mod)
      
      #summary(mod2)
      
      #mod2$B[8]
      #mod2$B.pvalue[8]
      #print(io)
      
      
      #effect.size<-mod2$B[8]
      #pval2<-mod2$B.pvalue[8]#Pvalue of phenotype
      
      esz[io,]<-effect.size
      pv_sp[io,]<-pval2 #p value for anova -- colors --- FDR for +/- symbols
    }
  }
}


results_ftrs<-cbind(esz,pv_sp)
results_ftrs<-as.data.frame(results_ftrs)
colnames(results_ftrs)<-c("Effect_size","Pvalue")
#results_ftrs<-results_ftrs[which(is.na(results_ftrs$Pvalue)==F),]
results_ftrs$FDR <- p.adjust(results_ftrs$Pvalue,method = "BH")
length(which(results_ftrs$Pvalue<0.05))
length(which(results_ftrs$FDR<0.05))

rownames(results_ftrs)<-sp_name2


effect_size_phen<-cbind(effect_size_phen,results_ftrs$Effect_size)
pval_phen<-cbind(pval_phen,results_ftrs$Pvalue)
fdr_phen<-cbind(fdr_phen,results_ftrs$FDR)



colnames(effect_size_phen)[dim(effect_size_phen)[2]]<-c("Live with partner")

lim<-dim(effect_size_phen)[2]



setwd("/WORKAREA/home/sv437/strain_assc_anpan/")
write.csv(effect_size_phen,'effect_size_phen_strain_without_phylo.csv')
write.csv(pval_phen,'pval_phen_strain_without_phylo.csv')
write.csv(fdr_phen,'fdr_phen_strain_without_phylo.csv')
save.image('strain_without_phylo_all.RData')
