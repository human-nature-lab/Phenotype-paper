

library(lmerTest)
library(vegan)
library(ape)

st<-read.csv('st_123.csv',row.names = 1)
st<-as.matrix.data.frame(st)
phen_all_use<-read.csv('phen_b4.csv',row.names = 1,header=T)
phen_all_use2<-read.csv('phen_all_raw2.csv',row.names=1,header=T)
data_transformed<-read.csv('mb_sp_01_clr.csv',row.names = 1,header=T)
st_ind<-read.csv('st_ind_03.csv',row.names=1)#Matched species name list with tree file
ind_sp2<-st_ind[,1]

mb_samp_sp<-data_transformed
data_transformed<-as.matrix.data.frame(data_transformed)

xm<-unique(phen_all_use$village_code)
for(i in 1:length(xm)){
  if(i==1){
    xm_length<-length(which(phen_all_use$village_code==xm[i]))
  }else{
    xm_length<-rbind(xm_length,length(which(phen_all_use$village_code==xm[i])))
  }
}

xm2<-xm[which(xm_length>5)]#These are the 19 villages in the cohort

phen_all_use19<-phen_all_use[which(phen_all_use$village_code%in%xm2),]
phen_all_use19_2<-phen_all_use2[which(phen_all_use$village_code%in%xm2),]
mb_samp_sp19<-mb_samp_sp[,which(phen_all_use$village_code%in%xm2)]
data_transformed19<-data_transformed[,which(phen_all_use$village_code%in%xm2)]
st19<-st[which(phen_all_use$village_code%in%xm2),1]

setwd("~/tree_files")

file_list<-list.files(path=".")







################################################################################
####
## ALMER model


library(evolvability)

effect_size_phen<-array(NA,dim=c(length(ind_sp2),dim(phen_all_use19_2)[2]))
pval_phen<-array(NA,dim=c(length(ind_sp2),dim(phen_all_use19_2)[2]))
effect_size_phen_no<-array(NA,dim=c(length(ind_sp2),dim(phen_all_use19_2)[2]))
pval_phen_no<-array(NA,dim=c(length(ind_sp2),dim(phen_all_use19_2)[2]))

setwd("~/tree_files")

for(i in 1:length(ind_sp2)){#length(ind_sp2)
  if(is.na(ind_sp2[i])==F){
    lmm<-0
    tr2 <- ape::read.tree(file_list[ind_sp2[i]])
    ind_stt<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(phen_all_use19))
    tr3<-drop.tip(tr2,which(is.na(ind_stt)==T))
    
    if(length(tr3$tip.label)>50){
      
      for(kl in 1:dim(phen_all_use19_2)[2]){#1:dim(phen_all_use19_2)[2]
        if(!(kl%in%(c(3,12)))){
          
          
          sp<-as.data.frame(data_transformed19[i,])
          cleaned_data<-array(0,dim=c(dim(phen_all_use19)[1],10))
          for(k in 1:dim(phen_all_use19)[1]){
            cleaned_data[k,1]<-as.numeric(phen_all_use19[k,6])#Age
            cleaned_data[k,2]<-as.numeric(phen_all_use19[k,4])#Sex
            cleaned_data[k,3]<-as.numeric(phen_all_use19[k,224])#Batch effect
            cleaned_data[k,4]<-as.numeric(phen_all_use19[k,225])#DNA concentration
            cleaned_data[k,5]<-as.numeric(phen_all_use19[k,40])#BMI
            cleaned_data[k,6]<-as.numeric(phen_all_use19[k,371])#Sampling date
            cleaned_data[k,7]<-as.numeric(st19[k])#Bristol stool scale
            cleaned_data[k,8]<-as.numeric(as.character(phen_all_use19_2[k,kl]))#Phenotype of interest
            cleaned_data[k,9]<-as.numeric(sp[k,1])#Species
            cleaned_data[k,10]<-as.numeric(as.character(phen_all_use19[k,11]))#Village code
          }
          cleaned_data<-as.data.frame(cleaned_data)
          rownames(cleaned_data)<-rownames(phen_all_use19)
          colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","phen","sp","village")
          cleaned_data<-cleaned_data[ind_stt[which(is.na(ind_stt)==F)],]
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
          cleaned_data2$village<-factor(cleaned_data2$village)
          
          tr4<-keep.tip(tr3,ind_temp)
          cleaned_data2<-cbind(cleaned_data2,"phyl"=tr4$tip.label)
          #cor_mat = ape::vcv.phylo(tr4, corr = TRUE)
          A1 <- Matrix::Matrix(ape::vcv(tr4), sparse = TRUE)
          colnames(A1) <- rownames(A1) <-tr4$tip.label
          
          
          if((sum(cleaned_data2$phen)>0)&(length(unique(cleaned_data2$phen))>1)){
            ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+phen+(1|village),
                      data = cleaned_data2)#The results are same as lmer, except pvalue -- find how they calculated
            mss<-summary(ms1)]
            
            if(length(coef(mss)[which(rownames(coef(mss))=="phen"),5])>0){
              
              effect_size_phen_no[i,kl]<-coef(mss)[which(rownames(coef(mss))=="phen"),1]
              pval_phen_no[i,kl]<-coef(mss)[which(rownames(coef(mss))=="phen"),5]
              
              ms3<-Almer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+phen+(1|village)+(1|phyl),
                         data = cleaned_data2,A=list(phyl=A1))#Use this!!!
              mss<-summary(ms3)
              
              
              effect_size_phen[i,kl]<-coef(summary(ms3))[which(rownames(coef(mss))=="phen"),1]
              pval_phen[i,kl]<-2*pnorm(abs(coef(mss)[which(rownames(coef(mss))=="phen"),3]), lower.tail = FALSE)#Pvalue of phenotype
              
            }else{
              effect_size_phen_no[i,kl]<-0
              pval_phen_no[i,kl]<-1
              effect_size_phen[i,kl]<-0
              pval_phen[i,kl]<-1
            }
            
            
          }else{
            effect_size_phen_no[i,kl]<-0
            pval_phen_no[i,kl]<-1
            effect_size_phen[i,kl]<-0
            pval_phen[i,kl]<-1
          }
          
        }else if(kl==3){#BMI
          cleaned_data<-array(0,dim=c(dim(phen_all_use19)[1],9))
          for(k in 1:dim(phen_all_use19)[1]){
            cleaned_data[k,1]<-as.numeric(phen_all_use19[k,6])#Age
            cleaned_data[k,2]<-as.numeric(phen_all_use19[k,4])#Sex
            cleaned_data[k,3]<-as.numeric(phen_all_use19[k,224])#Batch effect
            cleaned_data[k,4]<-as.numeric(phen_all_use19[k,225])#DNA concentration
            cleaned_data[k,5]<-as.numeric(phen_all_use19[k,40])#BMI
            cleaned_data[k,6]<-as.numeric(phen_all_use19[k,371])#Sampling date
            cleaned_data[k,7]<-as.numeric(st19[k])#Bristol stool scale
            cleaned_data[k,8]<-as.numeric(sp[k,1])#Species
            cleaned_data[k,9]<-as.numeric(as.character(phen_all_use19[k,11]))#Village code
          }
          cleaned_data<-as.data.frame(cleaned_data)
          colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","sp","village")
          cleaned_data<-cleaned_data[ind_stt[which(is.na(ind_stt)==F)],]
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
          cleaned_data2$village<-factor(cleaned_data2$village)
          
          tr4<-keep.tip(tr3,ind_temp)
          cleaned_data2<-cbind(cleaned_data2,"phyl"=tr4$tip.label)
          #cor_mat = ape::vcv.phylo(tr4, corr = TRUE)
          A1 <- Matrix::Matrix(ape::vcv(tr4), sparse = TRUE)
          colnames(A1) <- rownames(A1) <-tr4$tip.label
          
          #if((sum(cleaned_data2$phen)>0)&(length(unique(cleaned_data2$phen))>1)){
          
          ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+(1|village),
                    data = cleaned_data2)
          mss<-summary(ms1)
          coef(mss)[which(rownames(coef(mss))=="BMI"),5]
          
          if(length(coef(mss)[which(rownames(coef(mss))=="BMI"),5])>0){
            
            effect_size_phen_no[i,kl]<-coef(mss)[which(rownames(coef(mss))=="BMI"),1]
            pval_phen_no[i,kl]<-coef(mss)[which(rownames(coef(mss))=="BMI"),5]#Pvalue of phenotype
            
            ms3<-Almer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+(1|village)+(1|phyl),
                       data = cleaned_data2,A=list(phyl=A1))#Use this!!!
            mss<-summary(ms3)
            effect_size_phen[i,kl]<-coef(summary(ms3))[which(rownames(coef(mss))=="BMI"),1]
            pval_phen[i,kl]<-2*pnorm(abs(coef(mss)[which(rownames(coef(mss))=="BMI"),3]), lower.tail = FALSE)#Pvalue of phenotype
            
            
          }else{
            effect_size_phen_no[i,kl]<-0
            pval_phen_no[i,kl]<-1
            effect_size_phen[i,kl]<-0
            pval_phen[i,kl]<-1
          }
          
          
          
        }else if(kl==12){#Bristol stool scale
          cleaned_data<-array(0,dim=c(dim(phen_all_use19)[1],9))
          for(k in 1:dim(phen_all_use19)[1]){
            cleaned_data[k,1]<-as.numeric(phen_all_use19[k,6])#Age
            cleaned_data[k,2]<-as.numeric(phen_all_use19[k,4])#Sex
            cleaned_data[k,3]<-as.numeric(phen_all_use19[k,224])#Batch effect
            cleaned_data[k,4]<-as.numeric(phen_all_use19[k,225])#DNA concentration
            cleaned_data[k,5]<-as.numeric(phen_all_use19[k,40])#BMI
            cleaned_data[k,6]<-as.numeric(phen_all_use19[k,371])#Sampling date
            cleaned_data[k,7]<-as.numeric(st19[k])#Bristol stool scalet
            cleaned_data[k,8]<-as.numeric(sp[k,1])#Species
            cleaned_data[k,9]<-as.numeric(as.character(phen_all_use19[k,11]))#Village code
          }
          cleaned_data<-as.data.frame(cleaned_data)
          colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","sp","village")
          cleaned_data<-cleaned_data[ind_stt[which(is.na(ind_stt)==F)],]
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
          cleaned_data2$village<-factor(cleaned_data2$village)
          
          tr4<-keep.tip(tr3,ind_temp)
          cleaned_data2<-cbind(cleaned_data2,"phyl"=tr4$tip.label)
          #cor_mat = ape::vcv.phylo(tr4, corr = TRUE)
          A1 <- Matrix::Matrix(ape::vcv(tr4), sparse = TRUE)
          colnames(A1) <- rownames(A1) <-tr4$tip.label
          
          
          #if((sum(cleaned_data2$phen)>0)&(length(unique(cleaned_data2$phen))>1)){
          ms1<-lmer(sp~age+sex+batch_effect+dna_conc+bristol+sampling_date+bristol+(1|village),
                    data = cleaned_data2)
          mss<-summary(ms1)
          coef(mss)[which(rownames(coef(mss))=="bristol"),5]
          
          if(length(coef(mss)[which(rownames(coef(mss))=="bristol"),5])>0){
            
            effect_size_phen_no[i,kl]<-coef(mss)[which(rownames(coef(mss))=="bristol"),1]
            pval_phen_no[i,kl]<-coef(mss)[which(rownames(coef(mss))=="bristol"),5]#Pvalue of phenotype
            
            ms3<-Almer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+(1|village)+(1|phyl),
                       data = cleaned_data2,A=list(phyl=A1))#Use this!!!
            mss<-summary(ms3)

            effect_size_phen[i,kl]<-coef(summary(ms3))[which(rownames(coef(mss))=="bristol"),1]
            pval_phen[i,kl]<-2*pnorm(abs(coef(mss)[which(rownames(coef(mss))=="bristol"),3]), lower.tail = FALSE)#Pvalue of phenotype
            
            
          }else{
            effect_size_phen_no[i,kl]<-0
            pval_phen_no[i,kl]<-1
            effect_size_phen[i,kl]<-0
            pval_phen[i,kl]<-1
          }
          
          
        }
        
      }
    }
  }
  print(i)
  
}


write.csv(effect_size_phen,'esz_with_all.csv')
write.csv(pval_phen,'pval_with_all.csv')
write.csv(effect_size_phen_no,'esz_without_all.csv')
write.csv(pval_phen_no,'pval_without_all.csv')


##CHecking for phylogenetic signal

phh_l<-array(NA,dim=c(length(ind_sp2),dim(phen_all_use19_2)[2]))

for(i in 1:length(ind_sp2)){
  tr2 <- ape::read.tree(file_list[ind_sp2[i]])
  ind_stt<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(phen_all_use19))
  tr3<-drop.tip(tr2,which(is.na(ind_stt)==T))
  if(length(tr3$tip.label)>50){
    if(is.na(ind_sp2[i])==F){
      for(jk in 1:dim(phen_all_use19_2)[2]){
        cleaned_data<-array(0,dim=c(dim(phen_all_use19)[1]))
        for(k in 1:dim(phen_all_use19)[1]){
          cleaned_data[k]<-as.numeric(phen_all_use19_2[k,jk])
        }
        cleaned_data<-cleaned_data[ind_stt[which(is.na(ind_stt)==F)]]
        ind_temp<-0
        for(j in 1:dim(cleaned_data)[1]){
          na_score<-0
          if(is.na(cleaned_data[j])==T){
            na_score<-na_score+1
          }
          if(na_score==0){
            ind_temp<-rbind(ind_temp,j)
          }
        }
        cleaned_data2<-cleaned_data[ind_temp]
        tr4<-keep.tip(tr3,ind_temp)
        cleaned_data2<-as.data.frame(cleaned_data2)
        rownames(cleaned_data2)<-tr4$tip.label
        
        temp<-phylosig(tr4,cleaned_data2[,1], method="lambda", test=FALSE)
        if(length(temp$lambda)>0){
          phh_l[i,jk]<-temp$lambda
        }
        
      }
    }
  }
  print(i)
}


for(i in 1:dim(phh_l)[1]){
  for(j in 1:dim(phh_l)[2]){
    if(phh_l[i,j]<0.01){#Threshold cutoff decided based on overall distribution -- only including ones with reasonable phylogenetic signal
      effect_size_phen[i,j]<-NA
      pval_phen_no[i,kl]<-NA
      effect_size_phen[i,kl]<-NA
      pval_phen[i,kl]<-NA
    }
  }
}




#########################################################################
####plotting 

##Figure 3C

#IF taking too long! Read in results :-

e1<-read.csv('esz_with_all_st03ff.csv',row.names=1)
e2<-read.csv('esz_without_all_st03ff.csv',row.names=1)
f1<-read.csv('pval_with_all_st03ff.csv',row.names=1)
f2<-read.csv('pval_without_all_st03ff.csv',row.names=1)

#Other wise ---

e1<-effect_size_phen
e2<-effect_size_phen_no
f1<-pval_phen
f2<-pval_phen_no

#


t1<-read.csv('phen_all_raw2.csv',row.names=1)
t2<-read.csv('mb_sp_01_name.csv',row.names=1)

rownames(e1)<-t2[,1]
colnames(e1)<-colnames(t1)
rownames(e2)<-t2[,1]
colnames(e2)<-colnames(t1)
rownames(f1)<-t2[,1]
colnames(f1)<-colnames(t1)
rownames(f2)<-t2[,1]
colnames(f2)<-colnames(t1)

na_species<-0#3 species have all NA rows
for(i in 1:dim(e1)[1]){
  if(length(which(is.na(e1[i,])==F))>0){
    na_species<-c(na_species,i)
  }
}

e1<-e1[na_species,]
e2<-e2[na_species,]
f1<-f1[na_species,]
f2<-f2[na_species,]



e11<-c(0,0)
e12<-c(0,0)

for(i in 1:dim(f1)[1]){
  for(j in 1:dim(f1)[2]){
    if(is.na(f1[i,j])==F){
      if(f1[i,j]<0.05){
        #i1<-rbind(i1,i)
        #j1<-rbind(j1,j)
        e11<-rbind(e11,cbind(e2[i,j],e1[i,j]))
      }else if(f1[i,j]>0){
        e12<-rbind(e12,cbind(e2[i,j],e1[i,j]))
      }
    }
  }
}

e11<-as.data.frame(e11)
colnames(e11)<-c("x1","y1")
e12<-as.data.frame(e12)
colnames(e12)<-c("x2","y2")

e11<-e11[which(is.na(e11$y1)==F),]
e12<-e12[which(is.na(e12$y2)==F),]

library(ggplot2)
ggplot(e11, aes(x=x1, y=y1)) +
  geom_point(size=2, shape=23)+xlim(-2,2)+ylim(-2,2)

ms1<-lm(e11$y1~e11$x1)
colnames(e12)<-c("x1","y1")


e13<-cbind(cbind(e12,"shp"=array(21,dim=c(dim(e12)[1],1))),"coll2"=array("black",dim=c(dim(e12)[1],1)))



op<-cor.test(e13$x1,e13$y1,method="spearman")
rho<-op$estimate
pval<-op$p.value
op2<-lm(e13$y1~e13$x1)
coef(op2)
e14<-e13[c(dim(e13)[1]:-1:1),]

pdf('figure_3C.pdf',width=6,height=5)#was width=13,heigh=6(or 4)
ggplot(e14,aes(x=x1,y=y1))+geom_point(alpha=0.1,color="#4575b4")+geom_abline(intercept=0,slope=1,linetype="dashed", size=1.3)+xlim(-10,10)+ylim(-10,10)+
  #geom_abline(intercept=coef(op2)[1],slope=0.95,col="red", size=1)
  geom_segment(x=10,y=10*coef(op2)[2]+coef(op2)[1],xend=-10,yend=-10*coef(op2)[2]+coef(op2)[1],col="red", size=1)+#was 0.93
  geom_hline(yintercept=0,linetype="dashed", size=1.3)+xlim(-1,1)+ylim(-1,1)+
  theme(legend.position = "none",panel.background = element_blank(),panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),axis.text.x=element_text(size=15,color="black"),axis.text.y=element_text(size=15,color="black"),axis.line.y = element_line(color="black"),axis.line.x = element_line(color="black"))+
  xlab("Without strain-phylogenetic effect")+ylab("With strain-phylogenetic effect")+
  theme(plot.margin=margin(c(1,0,0,0), unit = "cm"))
dev.off()


######################################

##Figure 3A

setwd("~/strains batch 123")


library(lmerTest)
library(vegan)
library(ape)
library(ggtree)
library(phyr)

st<-read.csv('st_123.csv',row.names = 1)
st<-as.matrix.data.frame(st)
phen_all_use<-read.csv('phen_b4.csv',row.names = 1,header=T)
phen_all_use2<-read.csv('phen_all_raw2.csv',row.names=1,header=T)
data_transformed<-read.csv('mb_sp_01_clr.csv',row.names = 1,header=T)
st_ind<-read.csv('st_ind_03.csv',row.names=1) 
ind_sp2<-st_ind[,1]

mb_samp_sp<-data_transformed
data_transformed<-as.matrix.data.frame(data_transformed)

xm<-unique(phen_all_use$village_code)
for(i in 1:length(xm)){
  if(i==1){
    xm_length<-length(which(phen_all_use$village_code==xm[i]))
  }else{
    xm_length<-rbind(xm_length,length(which(phen_all_use$village_code==xm[i])))
  }
}

xm2<-xm[which(xm_length>5)]#These are the 19 villages in the cohort

phen_all_use19<-phen_all_use[which(phen_all_use$village_code%in%xm2),]
phen_all_use19_2<-phen_all_use2[which(phen_all_use$village_code%in%xm2),]
mb_samp_sp19<-mb_samp_sp[,which(phen_all_use$village_code%in%xm2)]
data_transformed19<-data_transformed[,which(phen_all_use$village_code%in%xm2)]
st19<-st[which(phen_all_use$village_code%in%xm2),1]

setwd("/strains batch 123/tree_files")

file_list<-list.files(path=".")

i<-25#Example for HHW
#i<-353#Example for fruits


setwd("/strains batch 123/tree_files")

#jk<-3#BMI
jk<-110#HHW
#jk<-76#Fruits consumption

tr2 <- ape::read.tree(file_list[ind_sp2[i]])
ind_stt<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(phen_all_use19))
tr3<-drop.tip(tr2,which(is.na(ind_stt)==T))
cleaned_data<-array(0,dim=c(dim(phen_all_use19)[1]))
for(k in 1:dim(phen_all_use19)[1]){
  cleaned_data[k]<-as.numeric(phen_all_use19_2[k,jk])
}
cleaned_data<-cleaned_data[ind_stt[which(is.na(ind_stt)==F)]]
ind_temp<-0
for(j in 1:dim(cleaned_data)[1]){
  na_score<-0
  if(is.na(cleaned_data[j])==T){
    na_score<-na_score+1
  }
  if(na_score==0){
    ind_temp<-rbind(ind_temp,j)
  }
}
cleaned_data2<-cleaned_data[ind_temp]
tr4<-keep.tip(tr3,ind_temp)
cleaned_data2<-as.data.frame(cleaned_data2)
rownames(cleaned_data2)<-tr4$tip.label

cleaned_data2<-cbind(cleaned_data2,array(NA,dim=c(dim(cleaned_data2)[1],1)))
colnames(cleaned_data2)[2]<-"color2"

cleaned_data2[,2]<-ifelse(is.na(cleaned_data2[,2])==T,ifelse((cleaned_data2[,1]==1),"#440154",cleaned_data2[,2]),cleaned_data2[,2])
cleaned_data2[,2]<-ifelse(is.na(cleaned_data2[,2])==T,ifelse((cleaned_data2[,1]==2),"#3b528b",cleaned_data2[,2]),cleaned_data2[,2])
cleaned_data2[,2]<-ifelse(is.na(cleaned_data2[,2])==T,ifelse((cleaned_data2[,1]==3),"#21918c",cleaned_data2[,2]),cleaned_data2[,2])
cleaned_data2[,2]<-ifelse(is.na(cleaned_data2[,2])==T,ifelse((cleaned_data2[,1]==4),"#5ec962",cleaned_data2[,2]),cleaned_data2[,2])
cleaned_data2[,2]<-ifelse(is.na(cleaned_data2[,2])==T,ifelse((cleaned_data2[,1]==5),"#fde725",cleaned_data2[,2]),cleaned_data2[,2])


setwd("/WORKAREA/home/sv437/strains batch 123")
ggm<-ggtree(tr4)+layout_dendrogram() +
  theme_dendrogram()+geom_tippoint(shape=20,color=cleaned_data2$color2)

ggsave(file=paste0('figure_3a.pdf'),plot=ggm,width=10,height=5)#For hhw



