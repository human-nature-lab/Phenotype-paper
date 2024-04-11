
library(curatedMetagenomicData)
library(dplyr)
library(lmerTest)
library(tidyverse)
library(vegan)

sampleMetadata2<- sampleMetadata |>
  #dplyr::filter(age >= 15) |>#Splitting age into 3 levels in regression: <18,18-60 & >60
  dplyr::filter(!base::is.na(BMI))|>
  dplyr::filter(body_site=="stool")|>
  dplyr::filter(study_name%in%c("AsnicarF_2021","HMP_2019_ibdmdb","HMP_2012","QinN_2014","KaurK_2020","LokmerA_2019","Obregon-TitoAJ_2015","PasolliE_2019","RubelMA_2020"))


sample_all<-sampleMetadata2[,which(colnames(sampleMetadata2)%in%c("study_name","sample_id","age","gender","BMI"))]

#Reading in Honduran metadata
hon_phen<-read.csv('phen_b4.csv',header=T)#Contains all other variables -- on each person
colnames(hon_phen)[1]<-"sample_ID"
hon_mb<-read.csv('mb_sp.csv',row.names = 1)#Honduran species profiles

# Filtering out redundant villages
xm<-tabulate(hon_phen$village_code)
hon_phen<-hon_phen %>% filter(village_code%in%which(xm>5))#Discarding very small villages(N<=5)-- to account for surveying outsiders
hon_mb<-hon_mb[,which(hon_phen$village_code%in%(which(xm>5)))]


#Reading in Metaphlan profiles from All other cohorts
mb_all_dataset_cm<-read.csv('all_dataset_mb.csv',row.names=1)
all_sample_ids<-read.csv('all_sample_ids.csv',header=F)
mb_all_dataset_cm_c<-mb_all_dataset_cm[,which(all_sample_ids[1,]%in%sample_all$sample_id)]
rm(mb_all_dataset_cm)

#Get species & filter for species -- with atleast 1 abundance
mb_all_dataset_cm_c2<-mb_all_dataset_cm_c[which(grepl("t__", rownames(mb_all_dataset_cm_c), fixed = TRUE)==T),]#Filter for species
mb_all_dataset_cm_c2<-mb_all_dataset_cm_c2[which(rowSums(mb_all_dataset_cm_c2)>0),]#Filter for species which have at least 1 individual


#Merge all other cohorts with Honduran cohort
sp_int<-unique(c(rownames(mb_all_dataset_cm_c2),rownames(hon_mb)))

#Initializing Master dataset with all cohorts
mb_all<-array(0,dim=c(length(sp_int),(dim(hon_mb)[2]+dim(mb_all_dataset_cm_c2)[2])))
dim(mb_all)
rownames(mb_all)<-sp_int
mb_all<-as.data.frame(mb_all)
colnames(mb_all)[1:length(hon_phen$sample_ID)]<-hon_phen$sample_ID
colnames(mb_all)[(length(hon_phen$sample_ID)+1):length(colnames(mb_all))]<-sample_all$sample_id


#Honduran datasets
mb_all[match(rownames(hon_mb),rownames(mb_all)),1:length(hon_phen$sample_ID)]<-hon_mb
temp2<-which(all_sample_ids[1,]%in%sample_all$sample_id)
temp_ind<-match(sample_all$sample_id,all_sample_ids[1,temp2])
mb_all[match(rownames(mb_all_dataset_cm_c2),rownames(mb_all)),(length(hon_phen$sample_ID)+temp_ind)]<-mb_all_dataset_cm_c2


# Calculate and attach alpha diversity
div_alpha<-array(NA,dim=c(dim(mb_all)[2]))

for(i in 1:dim(mb_all)[2]){
  div_alpha[i]<-diversity(mb_all[,i])#Default is Shannon
}

#Sample MEtadata -- merged

sample_all2<-array(NA,dim=c(dim(mb_all)[2],dim(sample_all)[2]))
colnames(sample_all2)<-colnames(sample_all)
sample_all2<-as.data.frame(sample_all2)
sample_all2$sample_id[1:length(hon_phen$sample_ID)]<-hon_phen$sample_ID
sample_all2$study_name[1:length(hon_phen$sample_ID)]<-"HMB_2024"
sample_all2[1:length(hon_phen$sample_ID),match(c("age","gender","BMI"),colnames(sample_all2))]<-hon_phen[,match(c("age_at_survey","gender","BMI"),colnames(hon_phen))]
sample_all2[(length(hon_phen$sample_ID)+1):length(colnames(mb_all)),]<-sample_all
sample_all2$gender<-ifelse(sample_all2$gender=="female","1",ifelse(sample_all2$gender=="male","0",sample_all2$gender))
sample_all2$Alpha_diversity<-div_alpha
sample_all2$age_cat<-ifelse(sample_all2$age<18,"age<18",ifelse(((sample_all2$age>=18)&(sample_all2$age<60)),"age(18-60)",ifelse(sample_all2$age>=60,"age>=60",sample_all2$age)))
sample_all2$age_cat<-factor(sample_all2$age_cat,levels=c("age<18","age(18-60)","age>=60"))

#CLR transformed relative abundance

taxa<-as.numeric(unlist(mb_all))

taxa25<-mb_all+min(taxa[which(as.numeric(taxa)>0)])/2
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x), na.rm=na.rm) / length(x))
}
Gmean_core = apply(taxa25, 1, gm_mean)
data_prepared = cbind(Gmean_core,taxa25)
data_transformed = t(apply(data_prepared,1,function(x){
  log(x / x[1])[-1]
}))
data_transformed<-as.data.frame(data_transformed)

## Meta-analysis starts

#Initializing functions

R_fromlm <- function(n,t) {
  r <- t / sqrt((t^2) + (n - 1))
  return(r)
}

compute_lm <- function(exprs_row){
  if(length(unique(sample_non_west22$age_cat))>1){
    lm_res <-  lm(spc ~ as.numeric(as.character(BMI)) + age_cat + gender,data=sample_non_west22)
  }else{
    lm_res <-  lm(spc ~ as.numeric(as.character(BMI)) + as.numeric(as.character(age)) + gender,data=sample_non_west22)
  }

  ez<-R_fromlm(length(temp_jm),summary(lm_res)[["coefficients"]]["as.numeric(as.character(BMI))", "t value"])
  #ez<-summary(lm_res)[["coefficients"]]["as.numeric(as.character(BMI))", "Estimate"]
  pv<-summary(lm_res)[["coefficients"]]["as.numeric(as.character(BMI))", "Pr(>|t|)"]
  res<-c(ez,pv)
  return(res)
}

#Initializing variables

require(foreach)
require(doParallel)
cl <- makeCluster(5)
registerDoParallel(cl)

output_lm <- foreach(jm=1:length(unique(sample_all2$study_name)), .combine=cbind) %dopar% {
  temp_jm<-which(sample_all2$study_name==unique(sample_all2$study_name)[jm])
  for(im in 1:dim(data_transformed)[1]){#dim(data_transformed)[1]
    sp1<-data_transformed[im,temp_jm]
    sample_non_west22<-cbind(sample_all2[temp_jm,],t(sp1))
    colnames(sample_non_west22)[dim(sample_non_west22)[2]]<-"spc"
    res11<-compute_lm(sample_non_west22)
    if(im==1){
      ezt<-res11[1]
      pvt<-res11[2]
    }else{
      ezt<-c(ezt,res11[1])
      pvt<-c(pvt,res11[2])
    }
  }
  print(list(ezt,pvt))
}



ez_nw<-as.data.frame(output_lm[1,])
pv_nw<-as.data.frame(output_lm[2,])
colnames(ez_nw)<-unique(sample_all2$study_name)
rownames(ez_nw)<-rownames(data_transformed)
colnames(pv_nw)<-unique(sample_all2$study_name)
rownames(pv_nw)<-rownames(data_transformed)

n_non<-table(sample_all2$study_name)
n_non2<-n_non[match(rownames(n_non),unique(sample_all2$study_name))]


#
for(j in 1:dim(ez_nw)[2]){
  if(j==1){
    temp<-cbind(rownames(ez_nw),ez_nw[,j],rep(as.numeric(n_non[match(unique(sample_all2$study_name)[j],rownames(n_non))]),dim(ez_nw)[1]),rep(unique(sample_all2$study_name)[j],dim(ez_nw)[1]))
    colnames(temp)<-c("Species","Cor","Number","Study")
    ez_non<-temp
  }else{
    temp<-cbind(rownames(ez_nw),ez_nw[,j],rep(as.numeric(n_non[match(unique(sample_all2$study_name)[j],rownames(n_non))]),dim(ez_nw)[1]),rep(unique(sample_all2$study_name)[j],dim(ez_nw)[1]))
    colnames(temp)<-c("Species","Cor","Number","Study")
    ez_non<-rbind(ez_non,temp)
  }
}

for(j in 1:dim(ez_nw)[1]){
  if(j==1){
    n_non_full<-as.numeric(n_non2)
    stud_lab_full<-as.character(unique(sample_all2$study_name))
  }else{
    n_non_full<-rbind(n_non_full,as.numeric(n_non2))
    stud_lab_full<-rbind(stud_lab_full,as.character(unique(sample_all2$study_name)))
  }
}

require(foreach)
require(doParallel)
cl <- makeCluster(5)
registerDoParallel(cl)

output_cor <- foreach(im=1:dim(ez_nw)[1], .combine=rbind) %dopar% {
  a2 <- meta::metacor(cor=as.numeric(ez_nw[im,]),
                      n=n_non_full[im,],
                      studlab=stud_lab_full[im,],
                      method.tau="PM",
                      sm="ZCOR")
  
  final_vector <-c(a2$TE.random,
                   a2$seTE.random,
                   paste(a2$lower.random,a2$upper.random,sep=";"),
                   a2$zval.random,
                   a2$pval.random,
                   a2$tau2,
                   a2$I2)
  
  names(final_vector) <- c("RE_Correlation","SE_RE","CI_RE","Zscore","p-value","tau2","I^2")
  print(list(a2$cor,final_vector))
}
stopCluster(cl)

a2_cor_all<-t(as.data.frame(output_cor[,1]))#All correlations
rownames(a2_cor_all)<-rownames(ez_nw)
colnames(a2_cor_all)<-colnames(ez_nw)
final_vector_all<-t(as.data.frame(output_cor[,2]))#All random effects for every species
final_vector_all<-as.data.frame(final_vector_all)
colnames(final_vector_all)<-c("RE_Correlation","SE_RE","CI_RE","Zscore","p-value","tau2","I^2")
rownames(final_vector_all)<-rownames(ez_nw)
final_vector_all$FDR_Qvalue <- p.adjust(as.numeric(final_vector_all$`p-value`),method = "BH")

length(which(final_vector_all$FDR_Qvalue<0.05))

#Check which species have high prevalence and min abundance -- 

sp_chosen<-rownames(final_vector_all)[which(final_vector_all$FDR_Qvalue<0.05)]
tempc<-match(sp_chosen,rownames(mb_all))

##
final_vector_all_c<-final_vector_all[tempc,]
ez_nw_c<-a2_cor_all[tempc,]
ez_nw_c<-ifelse(is.nan(ez_nw_c)==T,0,ez_nw_c)
# pv_nw_c<-pv_nw[tempc,]
# pv_nw_c<-ifelse(is.nan(pv_nw_c)==T,1,pv_nw_c)


#Making better species labels

mb_samp_sp_2<-array(0,dim=c(dim(ez_nw)[1],1))
sp_sig_full<-rownames(ez_nw)
for(i in 1:dim(ez_nw)[1]){
  ind_temp2<-unlist(gregexpr("s__",as.character(sp_sig_full[i]),fixed=T))
  ind_temp3<-unlist(gregexpr("t__",as.character(sp_sig_full[i]),fixed=T))
  ind_temp4<-unlist(gregexpr("g__GGB",as.character(sp_sig_full[i]),fixed=T))
  ind_temp5<-unlist(gregexpr("f__FGB",as.character(sp_sig_full[i]),fixed=T))
  ind_temp6<-unlist(gregexpr("o__OFGB",as.character(sp_sig_full[i]),fixed=T))
  ind_temp7<-unlist(gregexpr("c__CFGB",as.character(sp_sig_full[i]),fixed=T))
  if(ind_temp3>10){
    if(ind_temp7>0){
      ind_tempb<-unlist(gregexpr("p__",as.character(sp_sig_full[i]),fixed=T))
      mb_samp_sp_2[i,1]<-paste0("{",substr(as.character(sp_sig_full[i]),ind_tempb,ind_temp7-2),"}",substr(as.character(sp_sig_full[i]),ind_temp2+3,ind_temp3-2)," (",substr(as.character(sp_sig_full[i]),ind_temp3+3,nchar(as.character(sp_sig_full[i]))),")")
      
    }
    if((ind_temp6>0)&(ind_temp7<0)){
      ind_tempb<-unlist(gregexpr("c__",as.character(sp_sig_full[i]),fixed=T))
      mb_samp_sp_2[i,1]<-paste0("{",substr(as.character(sp_sig_full[i]),ind_tempb,ind_temp6-2),"}",substr(as.character(sp_sig_full[i]),ind_temp2+3,ind_temp3-2)," (",substr(as.character(sp_sig_full[i]),ind_temp3+3,nchar(as.character(sp_sig_full[i]))),")")
    }
    if((ind_temp5>0)&((ind_temp6+ind_temp7)<0)){
      ind_tempb<-unlist(gregexpr("o__",as.character(sp_sig_full[i]),fixed=T))
      mb_samp_sp_2[i,1]<-paste0("{",substr(as.character(sp_sig_full[i]),ind_tempb,ind_temp5-2),"}",substr(as.character(sp_sig_full[i]),ind_temp2+3,ind_temp3-2)," (",substr(as.character(sp_sig_full[i]),ind_temp3+3,nchar(as.character(sp_sig_full[i]))),")")
    }
    if((ind_temp4>0)&((ind_temp6+ind_temp7+ind_temp5)<0)){
      ind_tempb<-unlist(gregexpr("f__",as.character(sp_sig_full[i]),fixed=T))
      mb_samp_sp_2[i,1]<-paste0("{",substr(as.character(sp_sig_full[i]),ind_tempb,ind_temp4-2),"}",substr(as.character(sp_sig_full[i]),ind_temp2+3,ind_temp3-2)," (",substr(as.character(sp_sig_full[i]),ind_temp3+3,nchar(as.character(sp_sig_full[i]))),")")
    }
    if((ind_temp6+ind_temp7+ind_temp5+ind_temp4)<0){
      mb_samp_sp_2[i,1]<-paste0(substr(as.character(sp_sig_full[i]),ind_temp2+3,ind_temp3-2)," (",substr(as.character(sp_sig_full[i]),ind_temp3+3,nchar(as.character(sp_sig_full[i]))),")")
    }
  }
}




for(i in 1:dim(ez_nw_c)[1]){
  if(i==1){
    y_all_sigd<-cbind(as.data.frame(c(ez_nw_c[i,],final_vector_all_c$RE_Correlation[i])),as.data.frame(c("This study (N=1871,Honduras)","AsnicarF,2021 (N=1097,GBR/USA)","HMP,2012 (N=147,USA)","HMP,2019 (N=1274,USA)","Kaurk,2020 (N=31,India)","LokmerA,2019 (N=56,Cameroon)","Obregon-TitoAJ,2015 (N=14,Peru)","PasolliE,2019 (N=112,Madagascar)","QinN,2014 (N=237,China)","RubelMA,2020 (N=162,Cameroon)","Random effects")),rep(mb_samp_sp_2[tempc[i],1],dim(ez_nw_c)[2]+1),rep(strsplit(final_vector_all_c$CI_RE[i],";",fixed=T)[[1]][1],dim(ez_nw_c)[2]+1),rep(strsplit(final_vector_all_c$CI_RE[i],";",fixed=T)[[1]][2],dim(ez_nw_c)[2]+1),as.data.frame(c(as.numeric(pv_nw_c[i,]),1)))
    colnames(y_all_sigd)<-c("Estimate","Study","Species","y_max","y_min","PV")
  }else{
    tempd<-cbind(as.data.frame(c(ez_nw_c[i,],final_vector_all_c$RE_Correlation[i])),as.data.frame(c("This study (N=1871,Honduras)","AsnicarF,2021 (N=1097,GBR/USA)","HMP,2012 (N=147,USA)","HMP,2019 (N=1274,USA)","Kaurk,2020 (N=31,India)","LokmerA,2019 (N=56,Cameroon)","Obregon-TitoAJ,2015 (N=14,Peru)","PasolliE,2019 (N=112,Madagascar)","QinN,2014 (N=237,China)","RubelMA,2020 (N=162,Cameroon)","Random effects")),rep(mb_samp_sp_2[tempc[i],1],dim(ez_nw_c)[2]+1),rep(strsplit(final_vector_all_c$CI_RE[i],";",fixed=T)[[1]][1],dim(ez_nw_c)[2]+1),rep(strsplit(final_vector_all_c$CI_RE[i],";",fixed=T)[[1]][2],dim(ez_nw_c)[2]+1),as.data.frame(c(as.numeric(pv_nw_c[i,]),1)))
    colnames(tempd)<-c("Estimate","Study","Species","y_max","y_min","PV")
    y_all_sigd<-rbind(y_all_sigd,tempd)
  }
}

y_all_sigd$PV<-ifelse(is.nan(y_all_sigd$PV)==T,1,y_all_sigd$PV)

y_all_sigd$Significance<-ifelse(y_all_sigd$PV<0.05,"Significant",ifelse(y_all_sigd$Study=="Random effects","Random effects","Not Significant"))
y_all_sigd$Study<-factor(y_all_sigd$Study,levels=c("AsnicarF,2021 (N=1097,GBR/USA)","HMP,2012 (N=147,USA)","HMP,2019 (N=1274,USA)","QinN,2014 (N=237,China)","This study (N=1871,Honduras)","Kaurk,2020 (N=31,India)","LokmerA,2019 (N=56,Cameroon)","Obregon-TitoAJ,2015 (N=14,Peru)","PasolliE,2019 (N=112,Madagascar)","RubelMA,2020 (N=162,Cameroon)","Random effects"))
y_all_sigd$Species<-factor(y_all_sigd$Species,levels=mb_samp_sp_2[tempc[match(sort(as.numeric(as.character(final_vector_all_c$RE_Correlation))),as.numeric(as.character(final_vector_all_c$RE_Correlation)))],1])

library(ggplot2)

colnames(y_all_sigd)[7]<-"Significance"
y_all_sigd$Significance<-factor(y_all_sigd$Significance,levels=c("Random effects","Significant","Not Significant"))
y_all_sigd$Cohort<-ifelse(y_all_sigd$Study=="Random effects", "Random effects","Correlation coefficient")

ggplot(y_all_sigd, aes(x=Species, y=as.numeric(as.character(Estimate)), group=Study,color=Significance,shape = Cohort)) + #, fill=Study, shape=Study
  #geom_line() +
  geom_hline(yintercept = 0,linetype=2)+
  geom_errorbar(aes(ymin=as.numeric(as.character(y_min)), ymax=as.numeric(as.character(y_max))), width=.2,
                position=position_dodge(0.05),color='black',alpha=0.2)+
  geom_point(size=3,stroke=1)+#,color=y_all_sigd$PV_color
  scale_shape_manual(values=c(16,18))+#was 21
  #scale_shape_manual(values = c(3,2,1,11,8,13,6,7,4,14,18))+ 
  #scale_fill_manual(values=c("black","#d53e4f","#1a9850"))+
  scale_color_manual(values=c("black","#d53e4f","#1a9850"))+
  ylab("Correlation coefficient")+
  ggtitle("Body Mass Index (All cohorts)") +
  theme(text = element_text(size=17),axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black"),panel.background=element_blank(),axis.line.x = element_line(color="black", linewidth = 0.5),
        axis.line.y = element_line(color="black", linewidth = 0.5),plot.title = element_text(hjust = 0.5,face='bold'))+
  coord_flip()

ggsave('bmi_meta_analysis.pdf',width=15,height=9)




