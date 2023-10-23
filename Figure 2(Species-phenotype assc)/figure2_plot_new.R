
library(lmerTest)
require(foreach)
require(doParallel)
library(ggplot2)

st<-read.csv('st_123.csv',row.names = 1)
st<-as.matrix.data.frame(st)
phen_all_use<-read.csv('phen_b4.csv',row.names = 1,header=T)

phen_all_use2<-read.csv('phen_all_raw2.csv',row.names=1,header=T)#Phenotype of interest (sample x phenotyp)
data_transformed<-read.csv('mb_sp_01_clr.csv',row.names = 1,header=T)#CLR transformed species abundance data, species x sample
#phen_all_rct_use<-read.csv('phen_rct_b4.csv',row.names = 1,header=T)
mb_samp_sp<-data_transformed
data_transformed<-as.matrix.data.frame(data_transformed)

#setwd("~/Species-ph assc batch 123/19_village_results")

###################
##Accounting for people who migrated between villages to submit to survey (N=15)


xm<-unique(phen_all_use$village_code)
for(i in 1:length(xm)){
  if(i==1){
    xm_length<-length(which(phen_all_use$village_code==xm[i]))
  }else{
    xm_length<-rbind(xm_length,length(which(phen_all_use$village_code==xm[i])))
  }
}

xm2<-xm[which(xm_length>5)]#Filtering for people who migrated between villages to take survey

###

phen_all_use19<-phen_all_use[which(phen_all_use$village_code%in%xm2),]
phen_all_use19_2<-phen_all_use2[which(phen_all_use$village_code%in%xm2),]
mb_samp_sp19<-mb_samp_sp[,which(phen_all_use$village_code%in%xm2)]
data_transformed19<-data_transformed[,which(phen_all_use$village_code%in%xm2)]
st19<-st[which(phen_all_use$village_code%in%xm2),1]


#Initializing matrix for phenotype of interest
init_mat<-function(kl){
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
  colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","phen","sp","village")
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
  return(cleaned_data2)
}

#Initializing matrix when phenotype is a control variable
init_mat_no_phen<-function(){
  sp<-as.data.frame(data_transformed19[i,])
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
  return(cleaned_data2)
}


#################################
##Species phenotype association starts

cl <- makeCluster(3)#Initializing parralle for loops
registerDoParallel(cl)

effect_size_phen<-array(NA,dim=c(dim(data_transformed19)[1],dim(phen_all_use19_2)[2]))
pval_phen<-array(NA,dim=c(dim(data_transformed19)[1],dim(phen_all_use19_2)[2]))

multiResultClass <- function(result1=NULL,result2=NULL)
{
  me <- list(
    result1 = result1,
    result2 = result2
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}

#for(kl in 1:dim(phen_all_use19_2)[2]){#1:dim(phen_all_use19_2)[2]
output4<-foreach(kl=1:dim(phen_all_use19_2)[2],.multicombine=TRUE) %dopar% {
  result <- multiResultClass()
  temp<-array(NA,dim=c(dim(data_transformed19)[1],1))
  temp2<-array(NA,dim=c(dim(data_transformed19)[1],1))
  library(lmerTest)
  if(!(kl%in%(c(3,12)))){
    for(i in 1:dim(data_transformed19)[1]){
    
    cleaned_data3<-init_mat(kl)
    if((sum(cleaned_data3$phen)>0)&(length(unique(cleaned_data3$phen))>1)){
      ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+phen+(1|village),
                data = cleaned_data3)
      mss<-summary(ms1)
      coef(mss)[which(rownames(coef(mss))=="phen"),5]
      
      effect.size<-coef(mss)[which(rownames(coef(mss))=="phen"),1]
      pval2<-coef(mss)[which(rownames(coef(mss))=="phen"),5]#Pvalue of phenotype
      
      temp[i,1]<-effect.size
      temp2[i,1]<-pval2
      
    }else{
      temp[i,1]<-effect.size
      temp2[i,1]<-pval2
    }
    }
    
  }else if(kl==3){#BMI
      for(i in 1:dim(data_transformed19)[1]){
        
        cleaned_data3<-init_mat_no_phen()
      
        ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+(1|village),
                  data = cleaned_data3)
        mss<-summary(ms1)
        coef(mss)[which(rownames(coef(mss))=="BMI"),5]
        
        effect.size<-coef(mss)[which(rownames(coef(mss))=="BMI"),1]
        pval2<-coef(mss)[which(rownames(coef(mss))=="BMI"),5]#Pvalue of phenotype
        
        temp[i,1]<-effect.size
        temp2[i,1]<-pval2
      
    }
      
    }else if(kl==12){#Bristol stool scale
      for(i in 1:dim(data_transformed19)[1]){
        
        cleaned_data3<-init_mat_no_phen()
        ms1<-lmer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+(1|village),
                  data = cleaned_data3)
        mss<-summary(ms1)
        coef(mss)[which(rownames(coef(mss))=="bristol"),5]
        
        effect.size<-coef(mss)[which(rownames(coef(mss))=="bristol"),1]
        pval2<-coef(mss)[which(rownames(coef(mss))=="bristol"),5]#Pvalue of phenotype
        
        temp[i,1]<-effect.size
        temp2[i,1]<-pval2
      
    }
    }

  result$result1<-temp
  result$result2<-temp2
  return(result)
}
    
for(kl in 1:dim(phen_all_use19_2)[2]){
effect_size_phen[,kl]<-as.numeric(output4[[kl]]$result1)
pval_phen[,kl]<-as.numeric(output4[[kl]]$result2)
}

#

stopCluster(cl)

write.csv(effect_size_phen,'esz_vil_all.csv')
write.csv(pval_phen,'pval_vil_all.csv')


##Plotting

#Fig 2A

mb_samp_sp<-read.csv('mb_sp_01.csv',row.names=1)
mb_samp_sp<-read.csv('mb_sp_01_name.csv',row.names = 1)

phen_all_names<-read.csv('phen_all_raw2.csv',row.names = 1)

mb_samp_sp_2<-array(0,dim=c(dim(mb_samp_sp)[1],1))
lk<-array("#000000",dim=c(dim(mb_samp_sp)[1],1))
i<-14
ind_t<-0
for(i in 1:dim(mb_samp_sp)[1]){
  ind_temp2<-unlist(gregexpr("s__",as.character(mb_samp_sp[i,1]),fixed=T))
  ind_temp3<-unlist(gregexpr("t__",as.character(mb_samp_sp[i,1]),fixed=T))
  ind_temp4<-unlist(gregexpr("g__GGB",as.character(mb_samp_sp[i,1]),fixed=T))
  ind_temp5<-unlist(gregexpr("f__FGB",as.character(mb_samp_sp[i,1]),fixed=T))
  ind_temp6<-unlist(gregexpr("o__OFGB",as.character(mb_samp_sp[i,1]),fixed=T))
  ind_temp7<-unlist(gregexpr("c__CFGB",as.character(mb_samp_sp[i,1]),fixed=T))
  if(ind_temp3>10){
    if(ind_temp7>0){
      ind_tempb<-unlist(gregexpr("p__",as.character(mb_samp_sp[i,1]),fixed=T))
      mb_samp_sp_2[i,1]<-paste0("{",substr(as.character(mb_samp_sp[i,1]),ind_tempb,ind_temp7-2),"}",substr(as.character(mb_samp_sp[i,1]),ind_temp2+3,ind_temp3-2)," (",substr(as.character(mb_samp_sp[i,1]),ind_temp3+3,nchar(as.character(mb_samp_sp[i,1]))),")")
      ind_t<-rbind(ind_t,i)
      lk[i,1]<-"#feb24c"
      
    }
    if((ind_temp6>0)&(ind_temp7<0)){
      ind_tempb<-unlist(gregexpr("c__",as.character(mb_samp_sp[i,1]),fixed=T))
      mb_samp_sp_2[i,1]<-paste0("{",substr(as.character(mb_samp_sp[i,1]),ind_tempb,ind_temp6-2),"}",substr(as.character(mb_samp_sp[i,1]),ind_temp2+3,ind_temp3-2)," (",substr(as.character(mb_samp_sp[i,1]),ind_temp3+3,nchar(as.character(mb_samp_sp[i,1]))),")")
      ind_t<-rbind(ind_t,i)
      lk[i,1]<-"#feb24c"
    }
    if((ind_temp5>0)&((ind_temp6+ind_temp7)<0)){
      ind_tempb<-unlist(gregexpr("o__",as.character(mb_samp_sp[i,1]),fixed=T))
      mb_samp_sp_2[i,1]<-paste0("{",substr(as.character(mb_samp_sp[i,1]),ind_tempb,ind_temp5-2),"}",substr(as.character(mb_samp_sp[i,1]),ind_temp2+3,ind_temp3-2)," (",substr(as.character(mb_samp_sp[i,1]),ind_temp3+3,nchar(as.character(mb_samp_sp[i,1]))),")")
      ind_t<-rbind(ind_t,i)
      lk[i,1]<-"#feb24c"
    }
    if((ind_temp4>0)&((ind_temp6+ind_temp7+ind_temp5)<0)){
      ind_tempb<-unlist(gregexpr("f__",as.character(mb_samp_sp[i,1]),fixed=T))
      mb_samp_sp_2[i,1]<-paste0("{",substr(as.character(mb_samp_sp[i,1]),ind_tempb,ind_temp4-2),"}",substr(as.character(mb_samp_sp[i,1]),ind_temp2+3,ind_temp3-2)," (",substr(as.character(mb_samp_sp[i,1]),ind_temp3+3,nchar(as.character(mb_samp_sp[i,1]))),")")
      ind_t<-rbind(ind_t,i)
      lk[i,1]<-"#feb24c"
    }
    if((ind_temp6+ind_temp7+ind_temp5+ind_temp4)<0){
      mb_samp_sp_2[i,1]<-paste0(substr(as.character(mb_samp_sp[i,1]),ind_temp2+3,ind_temp3-2)," (",substr(as.character(mb_samp_sp[i,1]),ind_temp3+3,nchar(as.character(mb_samp_sp[i,1]))),")")
      ind_t<-rbind(ind_t,i)
    }
  }
}

#Reading in
#effect_size_phen_lmer<-read.csv('esz_vil_all.csv',row.names=1)
#pval_phen_lmer<-read.csv('pval_vil_all.csv',row.names=1)

for(j in 1:dim(pval_phen_lmer)[2]){
  if(j==1){
    pval_phen_lmer_flat<-pval_phen_lmer[,j]
  }else{
    pval_phen_lmer_flat<-c(pval_phen_lmer_flat,pval_phen_lmer[,j])
  }
}

fdr_phen_lmer_flat<-p.adjust(pval_phen_lmer_flat,method = "BH")

fdr_phen_lmer<-array(0,dim=c(dim(pval_phen_lmer)[1],dim(pval_phen_lmer)[2]))

for(j in 1:dim(pval_phen_lmer)[2]){
  fdr_phen_lmer[,j]<-fdr_phen_lmer_flat[((j-1)*dim(pval_phen_lmer)[1]+1):(j*dim(pval_phen_lmer)[1])]
}


mb_samp_sp_name2<-mb_samp_sp_2


##NEED to change all "N="

colnames(effect_size_phen_lmer)[1:13]<-c("Hemoglobin A1c","Mean arterial pressure","Body mass index","Heart rate","Oxygen saturation","Hemoglobin total","Poor (N=223)","Fair (N=969)","Very good (N=79)","Excellent (N=132)","Cough (N=413)","Bristol stool scale","Diarrhea (N=85)")
colnames(effect_size_phen_lmer)[14:20]<-c("Diabetes (N=40)","Allergies (N=124)","Heart disease (N=59)","Asthma (N=62)","Stomach illness (N=261)","Intestinal illness (N=103)","Arthritis (N=43)")#colnames(effect_size_phen_2)[14:24]
colnames(effect_size_phen_lmer)[21:27]<-c("Painkillers (N=1184)","Antibiotics (N=242)","Anti-diarrheal (N=60)","Anti-parasitic (N=49)","Anti-fungal (N=41)","Vitamins (N=335)","Anti-hypertensive (N=82)")#colnames(effect_size_phen_2)[25:28]
colnames(effect_size_phen_lmer)[28:35]<-c("Reserved (N=1666)","Nervous (N=950)","Openness (N=1550)","Cognitive impairment (N=100)","Dementia (N=185)","Alcohol daily frequency","Cigarette usage (N=82)","Cigarette frequency")
colnames(effect_size_phen_lmer)[36:41]<-c("Mild anxiety (N=395)","Moderate anxiety (N=123)","Severe anxiety (N=59)","Mild depression (N=469)","Moderate depression (N=138)","Severe depression (N=89)")
colnames(effect_size_phen_lmer)[42:65]<-c("Cat (N=960)","Dog (N=1384)","Parakeet (N=157)","No pets (N=220)","Cow (N=535)","Goat (N=63)","Pig (N=269)","Chicken (N=1664)","Duck (N=636)","Turkey (N=413)","Sheep (N=51)","Geese (N=82)","Horse (N=468)","Rabbit (N=176)","No farm animals (N=171)","Mice (N=820)","Bat (N=546)","Lizard (N=534)","Snake (N=546)","Bird (N=660)","Possum (N=588)","Rat (N=590)","Squirrel (N=51)","No wild animals (N=591)")
colnames(effect_size_phen_lmer)[c(72,77,79:80,82)]<-c("Cream/butter","Natural juice","Beef/pork","Ham/sausages/hotdog","Diet diversity score")
colnames(effect_size_phen_lmer)[c(66:71,73:76,78,81)]<-colnames(phen_all_names)[c(66:71,73:76,78,81)]


fdr_chk_lmer<-array(0,dim=c(dim(fdr_phen_lmer)[1],dim(fdr_phen_lmer)[2]))
for(i in 1:dim(fdr_chk_lmer)[1]){
  for(j in 1:dim(fdr_chk_lmer)[2]){
    if(is.na(fdr_phen_lmer[i,j])==F){
      if(fdr_phen_lmer[i,j]<0.05){
        fdr_chk_lmer[i,j]<-1
      }
    }
  }
}



colnames(effect_size_phen_lmer)[83:92]<-c("Living with partner (N=1232)","Number of partners","Grades 1-3 (N=920)","Grades 4-6 (N=528)","Grades >6 (N=129)","Friend ties (same building)","Friend ties (different building)","Betweenness (friendship)","Transitivity (friendship)","Familial ties (same building)")
colnames(effect_size_phen_lmer)[93:102]<-c("Family ties (different building)","Betweeness (familial)","Transitivity (familial)","Degree (all ties)","Clustering coefficient (all ties)","Betweenness (all ties)","Kin percentage (to third degree)","Altruism","Risk taking","Washing hands (N=1603)")
colnames(effect_size_phen_lmer)[103:112]<-c("Distance to village center","Distance to main road","Number of churches","Altitude","Travel","Monthly expenditure","Household size","Household wealth index","TV (N=958)","No electronics (N=89)")
colnames(effect_size_phen_lmer)[113:123]<-c("Earth/Sand floor (N=594)","Ceramic floor (N=151)","Glass windows (N=86)","Unfinished windows (N=59)","Clay/mud walls (N=1256)","Cement walls (N=402)","Concrete roof (N=37)","Sleeping rooms","Spring(protected) (N=1441)","Tube well (N=91)","Dug well(unprotected) (N=32)")


effect_size_phen_lmer_all<-effect_size_phen_lmer
fdr_phen_lmer_all<-fdr_phen_lmer
length(which(fdr_phen_lmer_all<0.05))

fdr_chk_lmer_all<-rowSums(fdr_chk_lmer)
length(which(fdr_chk_lmer_all>9))#was >8

ind_c1<-which(fdr_chk_lmer_all>9)

ind_phen<-c(1,3:4,14,15,18:19,22,24,27,33,30,32,38,41,43,50,51,65:67,69:82,87,83,90,94,99,106,107,110:115,117,118,120,121)

fdr2<-fdr_chk_lmer[,ind_phen]

ind_c2<-0

for(i in 1:dim(fdr2)[2]){
  if(sum(fdr2[,i])<5){#was 15
    ind_c2<-c(ind_c2,which(fdr2[,i]==1))
  }
}

ind_c2<-ind_c2[2:length(ind_c2)]

ind_cc<-unique(c(ind_c1,ind_c2))

effect_size_phen_sig<-effect_size_phen_lmer_all[ind_cc,ind_phen]
fdr_phen_sig<-fdr_chk_lmer[ind_cc,ind_phen]

effect_size_phen_sig_plot<-effect_size_phen_sig

for(i in 1:dim(effect_size_phen_sig)[1]){
  for(j in 1:dim(effect_size_phen_sig)[2]){
    if((fdr_phen_sig[i,j]==1)){
      effect_size_phen_sig_plot[i,j]<-effect_size_phen_sig[i,j]
    }else{
      effect_size_phen_sig_plot[i,j]<-0
    }
  }
}

disp_fdr<-array("",dim=c(dim(fdr_phen_sig)[1],dim(fdr_phen_sig)[2]))
for(i in 1:dim(disp_fdr)[1]){
  for(j in 1:dim(disp_fdr)[2]){
    if((fdr_phen_sig[i,j]==1)){
      if(effect_size_phen_sig_plot[i,j]>0){
        disp_fdr[i,j]<-"+"
      }else if(effect_size_phen_sig_plot[i,j]<0){
        disp_fdr[i,j]<-"-"
      }else{
        disp_fdr[i,j]<-""
      }
    }
  }
}
#max(fdr_phen_sig)

effect_size_phen_sig_plot<-t(effect_size_phen_sig_plot)
rownames(effect_size_phen_sig_plot)<-colnames(effect_size_phen_lmer)[ind_phen]
colnames(effect_size_phen_sig_plot)<-as.character(mb_samp_sp_2[ind_cc,1])
paletteLength<-15
myColor <- colorRampPalette(c("#cb181d","white", "#33a02c"))(paletteLength)#myColor <- colorRampPalette(c("Orange", "white", "Darkblue"))(paletteLength)
myBreaks <- c(seq(min(effect_size_phen_sig_plot),0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(effect_size_phen_sig_plot)/paletteLength, max(effect_size_phen_sig_plot), length.out=floor(paletteLength/2)))

myBreaks
myBreaks[7]<-as.numeric(-0.001)
myBreaks[9]<-as.numeric(0.001)
myBreaks

library(pheatmap)

#pheatmap(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,cluster_rows = F,fontsize_number = 13,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr),number_color = "grey30")
out2<-pheatmap(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
#Save 12.00 x 14.00 for shorter plot -- still rotate and entire page for this plot

ind_ca2<-out2$tree_col[["order"]]
mb_samp_sp_name2[ind_cc[ind_ca2]]

#Choices - 09-03-23
#ind_ca3<-ind_ca2[c(1:22,81:-1:23)]#Choice 1
#ind_ca3<-ind_ca2[c(23:81,22:-1:1)]#Choice 2
phen_ord<-order(as.numeric(effect_size_phen_sig_plot[46,]))#Earth/sand floor
ind_ca3<-phen_ord#Choice 3

#ind_ca3<-phen_ord[length(phen_ord):-1:1]#Choice 4


#Bubble plot

myb2<-c(0.00000100,0.1,0.5208435,0.9,1.15,1.432,1.9361451,4.16674833,8)
disp_fdr5<-t(disp_fdr[ind_ca3,])
myColor#Color

#myBreaks6<-c(myb[c(length(myb):-1:1)]*-1,myb)
myBreaks7<-c(myb2[c(length(myb2):-1:1)]*-1,myb2)
myColor6<-colorRampPalette(c("#cb181d","white", "#33a02c"))(17)
effect_size_phen_sig_plot3<-effect_size_phen_sig_plot[,ind_ca3]

out3<-pheatmap(effect_size_phen_sig_plot[,ind_ca3],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor6,breaks=myBreaks7,display_numbers = t(disp_fdr[ind_ca3,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = F)


for(i in 1:dim(effect_size_phen_sig_plot3)[1]){
  for(j in 1:dim(effect_size_phen_sig_plot3)[2]){
    if((i==1)&(j==1)){
      if(disp_fdr5[i,j]==""){
        temp<-c(i,j,NA)
      }else if(disp_fdr5[i,j]=="-"){
        temp<-c(i,j,"-")
      }else if(disp_fdr5[i,j]=="+"){
        temp<-c(i,j,3)
      }
      #if(effect_size_phen_sig_plot3[i,j]<myBreaks6[1]){
      #  temp<-c(temp,myColor6[1])
      #}else 
      if(effect_size_phen_sig_plot3[i,j]>myBreaks7[length(myColor6)]){
        temp<-c(temp,myColor6[length(myColor6)])
      }
      for(k in 2:length(myColor6)){
        if((effect_size_phen_sig_plot3[i,j]>=myBreaks7[k-1])&(effect_size_phen_sig_plot3[i,j]<=myBreaks7[k])){
          temp<-c(temp,myColor6[k-1])
        }
      }
      temp<-c(temp,rownames(effect_size_phen_sig_plot3)[i],colnames(effect_size_phen_sig_plot3)[j])
      d<-temp
    }else{
      if(disp_fdr5[i,j]==""){
        temp<-c(i,j,NA)
      }else if(disp_fdr5[i,j]=="-"){
        temp<-c(i,j,45)
      }else if(disp_fdr5[i,j]=="+"){
        temp<-c(i,j,43)
      }
      ##if(effect_size_phen_sig_plot3[i,j]<myBreaks6[1]){
      #  temp<-c(temp,myColor6[1])
      #}else 
      if(effect_size_phen_sig_plot3[i,j]>myBreaks7[length(myColor6)]){
        temp<-c(temp,myColor6[length(myColor6)])
      }
      for(k in 2:length(myColor6)){
        if((effect_size_phen_sig_plot3[i,j]>=myBreaks7[k-1])&(effect_size_phen_sig_plot3[i,j]<=myBreaks7[k])){
          temp<-c(temp,myColor6[k-1])
        }
      }
      temp<-c(temp,rownames(effect_size_phen_sig_plot3)[i],colnames(effect_size_phen_sig_plot3)[j])
      d<-rbind(d,temp)
      #sz_t<-rbind(sz_t,length(temp))
    }
    
  }
}

##Change order of species --hclust
#out<-pheatmap(effect_size_phen_sig_plot2,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot2),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr2),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
#colnames(effect_size_phen_sig_plot2[out$tree_col["order"],])


d<-as.data.frame(d)

colnames(d)<-c("x1","y1","shp","clr","x11","y11")

library(ggplot2)
d$x1<-factor(d$x1,levels=c(dim(effect_size_phen_sig_plot3)[1]:-1:1))
d$y1<-factor(d$y1,levels=1:dim(effect_size_phen_sig_plot3)[2])
d$x11<-factor(d$x11,levels=rownames(effect_size_phen_sig_plot3)[c(length(rownames(effect_size_phen_sig_plot3)):-1:1)])
d$y11<-factor(d$y11,levels=colnames(effect_size_phen_sig_plot3))
#d$y11<-factor(d$y11,levels=colnames(effect_size_phen_sig_plot3)[out2$tree_col[["order"]]])

d2<-d

for(i in 1:dim(d)[1]){
  if(is.na(d$shp[i])==T){
    d2$clr[i]<-"#FFFFFF"
  }
}



#out2<-pheatmap(d2,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot2),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor5,breaks=myBreaks5,display_numbers = disp_fdr4,number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

clk<-d2$clr
clk<-ifelse(d2$clr=="#FFFFFF","#FFFFFF","#000000")

#write.csv(out2$tree_col[["order"]],'sp_s_sel.csv')
effect_size_phen_sig_plot12<-effect_size_phen_sig_plot3
#"#377eb8","#e5d8bd","#33a02c","#6a3d9a","#a50026","#dd3497","#dfc27d","#ff7f00","#000000","#525252","#a6bddb","#cab2d6"

pdf('figure_2A.pdf',width=15,height=11)
ggplot(d2,aes(y11,x11,fill=clr))+
  geom_point(shape=20,color=clk,stroke=0,size=5.9)+#Extra line to give it outline
  geom_point(shape=20,color=d2$clr,stroke=0,size=5.5)+
  theme(axis.text.x=element_text(angle=90,hjust=0.95,color="black",face="bold",vjust=0.2),legend.position = "none",axis.text.y=element_text(color="black",face="bold"),panel.background = element_blank())+#plot.margin = margin(c(-20,(dim(effect_size_phen_sig_plot3)[2]+4),-40,(dim(effect_size_phen_sig_plot3)[1]+10)))
  xlab("")+ylab("")+expand_limits(x=c(-0.5,(dim(effect_size_phen_sig_plot3)[2]+10)))+
  theme(plot.margin=margin(c(1,12,1,1), unit = "pt"))+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot12)[2]+1,y=dim(effect_size_phen_sig_plot12)[1]+0.5,yend=dim(effect_size_phen_sig_plot12)[1]+0.5,col="#4575b4",size=0.3)+#Physiological variables
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot12)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot12)[2]+1),ymin=(dim(effect_size_phen_sig_plot12)[1]-2.5),ymax=(dim(effect_size_phen_sig_plot12)[1]+0.5),color="#4575b4",fill="#4575b4")+
  annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+5),y=(dim(effect_size_phen_sig_plot12)[1]-0.5),label="Physiological",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+4.3),y=(dim(effect_size_phen_sig_plot12)[1]-1.7),label="variables",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot12)[1]+0.5),ymax=(dim(effect_size_phen_sig_plot12)[1]-2.5),color="#4575b4",fill="#4575b4")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot12)[2]+1,y=dim(effect_size_phen_sig_plot12)[1]-2.5,yend=dim(effect_size_phen_sig_plot12)[1]-2.5,col="#a6cee3",size=0.3)+#Chronic conditions
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot12)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot12)[2]+1),ymin=(dim(effect_size_phen_sig_plot12)[1]-6.5),ymax=(dim(effect_size_phen_sig_plot12)[1]-2.5),color="#a6cee3",fill="#a6cee3")+
  annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+3.4),y=(dim(effect_size_phen_sig_plot12)[1]-3.5),label="Chronic",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+4.3),y=(dim(effect_size_phen_sig_plot12)[1]-4.7),label="conditions",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot12)[1]-2.5),ymax=(dim(effect_size_phen_sig_plot12)[1]-6.5),color="#a6cee3",fill="#a6cee3")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot12)[2]+1,y=dim(effect_size_phen_sig_plot12)[1]-6.5,yend=dim(effect_size_phen_sig_plot12)[1]-6.5,col="#762a83",size=0.3)+#Medication
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot12)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot12)[2]+1),ymin=(dim(effect_size_phen_sig_plot12)[1]-9.5),ymax=(dim(effect_size_phen_sig_plot12)[1]-6.5),color="#762a83",fill="#762a83")+
  annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+4.3),y=(dim(effect_size_phen_sig_plot12)[1]-8),label="Medications",color="black",fontface="bold")+
  #annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+3.5),y=(dim(effect_size_phen_sig_plot12)[1]-8.2),label="(non-kin)",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot12)[1]-6.5),ymax=(dim(effect_size_phen_sig_plot12)[1]-9.5),color="#762a83",fill="#762a83")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot12)[2]+1,y=dim(effect_size_phen_sig_plot12)[1]-9.5,yend=dim(effect_size_phen_sig_plot12)[1]-9.5,col="#dfc27d",size=0.3)+#Alcohol
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot12)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot12)[2]+1),ymin=(dim(effect_size_phen_sig_plot12)[1]-10.5),ymax=(dim(effect_size_phen_sig_plot12)[1]-9.5),color="#dfc27d",fill="#dfc27d")+
  annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+3.4),y=(dim(effect_size_phen_sig_plot12)[1]-10),label="Alcohol",color="black",fontface="bold")+
  #annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+4.3),y=(dim(effect_size_phen_sig_plot12)[1]-9.7),label="variables",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot12)[1]-9.5),ymax=(dim(effect_size_phen_sig_plot12)[1]-10.5),color="#dfc27d",fill="#dfc27d")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot12)[2]+1,y=dim(effect_size_phen_sig_plot12)[1]-10.5,yend=dim(effect_size_phen_sig_plot12)[1]-10.5,col="#e7298a",size=0.3)+#Mental health
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot12)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot12)[2]+1),ymin=(dim(effect_size_phen_sig_plot12)[1]-14.5),ymax=(dim(effect_size_phen_sig_plot12)[1]-10.5),color="#e7298a",fill="#e7298a")+
  annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+3.4),y=(dim(effect_size_phen_sig_plot12)[1]-12),label="Mental",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+3.7),y=(dim(effect_size_phen_sig_plot12)[1]-13.2),label="health",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot12)[1]-10.5),ymax=(dim(effect_size_phen_sig_plot12)[1]-14.5),color="#e7298a",fill="#e7298a")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot12)[2]+1,y=dim(effect_size_phen_sig_plot12)[1]-14.5,yend=dim(effect_size_phen_sig_plot12)[1]-14.5,col="#807dba",size=0.3)+#Animals
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot12)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot12)[2]+1),ymin=(dim(effect_size_phen_sig_plot12)[1]-18.5),ymax=(dim(effect_size_phen_sig_plot12)[1]-14.5),color="#807dba",fill="#807dba")+
  annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+3.3),y=(dim(effect_size_phen_sig_plot12)[1]-16.5),label="Animals",color="black",fontface="bold")+
  #annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+4.8),y=(dim(effect_size_phen_sig_plot12)[1]-14.2),label="(kin + non-kin)",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot12)[1]-14.5),ymax=(dim(effect_size_phen_sig_plot12)[1]-18.5),color="#807dba",fill="#807dba")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot12)[2]+1,y=dim(effect_size_phen_sig_plot12)[1]-18.5,yend=dim(effect_size_phen_sig_plot12)[1]-18.5,col="#adba07",size=0.3)+#Diet
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot12)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot12)[2]+1),ymin=(dim(effect_size_phen_sig_plot12)[1]-34.5),ymax=(dim(effect_size_phen_sig_plot12)[1]-18.5),color="#adba07",fill="#adba07")+
  annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+2.7),y=(dim(effect_size_phen_sig_plot12)[1]-26.5),label="Diet",color="black",fontface="bold")+
  #annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+3.5),y=(dim(effect_size_phen_sig_plot12)[1]-18.7),label="behavior",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot12)[1]-18.5),ymax=(dim(effect_size_phen_sig_plot12)[1]-34.5),color="#adba07",fill="#adba07")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot12)[2]+1,y=dim(effect_size_phen_sig_plot12)[1]-34.5,yend=dim(effect_size_phen_sig_plot12)[1]-34.5,col="#f9ce1d",size=0.3)+#Education
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot12)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot12)[2]+1),ymin=(dim(effect_size_phen_sig_plot12)[1]-35.5),ymax=(dim(effect_size_phen_sig_plot12)[1]-34.5),color="#f9ce1d",fill="#f9ce1d")+
  annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+4.3),y=(dim(effect_size_phen_sig_plot12)[1]-35),label="Education",color="black",fontface="bold")+
  #annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+3.8),y=(dim(effect_size_phen_sig_plot12)[1]-24.2),label="factors",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot12)[1]-34.5),ymax=(dim(effect_size_phen_sig_plot12)[1]-35.5),color="#f9ce1d",fill="#f9ce1d")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot12)[2]+1,y=dim(effect_size_phen_sig_plot12)[1]-35.5,yend=dim(effect_size_phen_sig_plot12)[1]-35.5,col="#d94801",size=0.3)+#Social factors
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot12)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot12)[2]+1),ymin=(dim(effect_size_phen_sig_plot12)[1]-39.5),ymax=(dim(effect_size_phen_sig_plot12)[1]-35.5),color="#d94801",fill="#d94801")+
  annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+3.1),y=(dim(effect_size_phen_sig_plot12)[1]-36.5),label="Social",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+3.5),y=(dim(effect_size_phen_sig_plot12)[1]-37.7),label="factors",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot12)[1]-35.5),ymax=(dim(effect_size_phen_sig_plot12)[1]-39.5),color="#d94801",fill="#d94801")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot12)[2]+1,y=dim(effect_size_phen_sig_plot12)[1]-39.5,yend=dim(effect_size_phen_sig_plot12)[1]-39.5,col="#4575b4",size=0.3)+#Village factors
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot12)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot12)[2]+1),ymin=(dim(effect_size_phen_sig_plot12)[1]-40.5),ymax=(dim(effect_size_phen_sig_plot12)[1]-39.5),color="#4575b4",fill="#4575b4")+
  annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+3.4),y=(dim(effect_size_phen_sig_plot12)[1]-39.5),label="Village",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+3.7),y=(dim(effect_size_phen_sig_plot12)[1]-40.7),label="factors",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot12)[1]-39.5),ymax=(dim(effect_size_phen_sig_plot12)[1]-40.5),color="#4575b4",fill="#4575b4")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot12)[2]+1,y=dim(effect_size_phen_sig_plot12)[1]-40.5,yend=dim(effect_size_phen_sig_plot12)[1]-40.5,col="#246556",size=0.3)+#Economic factors
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot12)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot12)[2]+1),ymin=(dim(effect_size_phen_sig_plot12)[1]-50.5),ymax=(dim(effect_size_phen_sig_plot12)[1]-40.5),color="#246556",fill="#246556")+##00b50b##e7298a
  annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+4),y=(dim(effect_size_phen_sig_plot12)[1]-45.5),label="Economic",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+4.3),y=(dim(effect_size_phen_sig_plot12)[1]-46.7),label="factors",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot12)[1]-40.5),ymax=(dim(effect_size_phen_sig_plot12)[1]-50.5),color="#246556",fill="#246556")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot12)[2]+1,y=dim(effect_size_phen_sig_plot12)[1]-50.5,yend=dim(effect_size_phen_sig_plot12)[1]-50.5,col="#000000",size=0.3)+#Water sources
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot12)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot12)[2]+1),ymin=(dim(effect_size_phen_sig_plot12)[1]-51.5),ymax=(dim(effect_size_phen_sig_plot12)[1]-50.5),color="#000000",fill="#000000")+
  annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+5),y=(dim(effect_size_phen_sig_plot12)[1]-51),label="Water source",color="black",fontface="bold")+
  #annotate("text",x=(dim(effect_size_phen_sig_plot12)[2]+3.3),y=(dim(effect_size_phen_sig_plot12)[1]-52.7),label="Source",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot12)[1]-50.5),ymax=(dim(effect_size_phen_sig_plot12)[1]-51.5),color="#000000",fill="#000000")

dev.off()

#Fig 2B,C,D

### Figure 2B

div_alpha_ph2<-read.csv('div_alpha_ph2_ch.csv',row.names=1)
div_alpha_ph2$chronic <- factor(div_alpha_ph2$chronic, levels=c("Healthy","Diabetes","Allergies","Heart disease","Asthma","Stomach illness","Intestinal illness","Arthritis"),ordered=TRUE)


#Only significant comparisons
p_ch<-array(NA,dim=c(3,3))

x1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$chronic=="Healthy")]
y1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$chronic=="Allergies")]
x1<-x1[which(is.na(x1)==F)]
y1<-y1[which(is.na(y1)==F)]

p_ch[1,3]<-wilcox.test(as.numeric(as.character(x1)),as.numeric(as.character(y1)))$p.value
p_ch[1,1:2]<-c("Healthy","Allergies")

x1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$chronic=="Healthy")]
y1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$chronic=="Stomach illness")]
x1<-x1[which(is.na(x1)==F)]
y1<-y1[which(is.na(y1)==F)]

p_ch[2,3]<-wilcox.test(as.numeric(as.character(x1)),as.numeric(as.character(y1)))$p.value
p_ch[2,1:2]<-c("Healthy","Stomach illness")

x1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$chronic=="Healthy")]
y1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$chronic=="Intestinal illness")]
x1<-x1[which(is.na(x1)==F)]
y1<-y1[which(is.na(y1)==F)]

p_ch[3,3]<-wilcox.test(as.numeric(as.character(x1)),as.numeric(as.character(y1)))$p.value
p_ch[3,1:2]<-c("Healthy","Intestinal illness")

p_ch<-as.data.frame(p_ch)
colnames(p_ch)<-c("group1","group2","pval")

p_ch[,3]<-format(as.numeric(as.character(p_ch[,3])),scientific=T,digits = 3)

my_comparisons <- list( c("Healthy", "Allergies"), c("Healthy", "Stomach illness"), c("Healthy", "Intestinal illness") )

div_alpha_ph22<-div_alpha_ph2[which(is.na(div_alpha_ph2$alpha_div)==F),]
div_mm<-mean(div_alpha_ph22$alpha_div[which(div_alpha_ph22$chronic=="Healthy")])

library(ggpubr)
pdf('figure_2B.pdf',width=4,height=6)
ggplot(div_alpha_ph2, aes(x=chronic, y=as.numeric(as.character(alpha_div)),group=chronic)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#1f78b4","#b2df8a","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6"),outlier.alpha=0)+labs(x="",y="Shannon diversity")+
  theme(text = element_text(size=20),axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black",angle=60,hjust=1),panel.background=element_blank(),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))+ scale_x_discrete(breaks=seq(1,8,1))+
  stat_pvalue_manual(data=p_ch,label = "pval", y.position = c(5.2, 5.5, 5.8),size=5)+ylim(1.5,6.5)+#+
  geom_segment(x=0,xend=9,y=(div_mm+0.05),yend=(div_mm+0.05),col="#878787",linetype="longdash",size=1)+#+
  theme(plot.margin=margin(c(0.5,0.5,0.5,1), unit = "cm"))+scale_x_discrete(labels=c("Healthy","Diabetes","Allergies","Heart disease","Asthma","Stomach illness","Intestinal illness","Arthritis"))
dev.off()

###Figure 2C

div_alpha_ph2<-read.csv('div_alpha_ph2_anti.csv',row.names=1)

div_alpha_ph2$medication <- factor(div_alpha_ph2$medication, levels=c("No medication","Antibiotic","Anti-diarrheal","Anti-parasitic","Anti-fungal","Vitamins","Anti-hypertensive"),ordered=TRUE)

p_ch<-array(NA,dim=c(3,3))

x1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$medication=="No medication")]
y1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$medication=="Antibiotic")]
x1<-x1[which(is.na(x1)==F)]
y1<-y1[which(is.na(y1)==F)]

p_ch[1,3]<-wilcox.test(as.numeric(as.character(x1)),as.numeric(as.character(y1)))$p.value
p_ch[1,1:2]<-c("No medication","Antibiotic")

x1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$medication=="No medication")]
y1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$medication=="Anti-diarrheal")]
x1<-x1[which(is.na(x1)==F)]
y1<-y1[which(is.na(y1)==F)]

p_ch[2,3]<-wilcox.test(as.numeric(as.character(x1)),as.numeric(as.character(y1)))$p.value
p_ch[2,1:2]<-c("No medication","Anti-diarrheal")

x1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$medication=="No medication")]
y1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$medication=="Anti-parasitic")]
x1<-x1[which(is.na(x1)==F)]
y1<-y1[which(is.na(y1)==F)]

p_ch[3,3]<-wilcox.test(as.numeric(as.character(x1)),as.numeric(as.character(y1)))$p.value
p_ch[3,1:2]<-c("No medication","Anti-parasitic")

p_ch<-as.data.frame(p_ch)
colnames(p_ch)<-c("group1","group2","pval")

p_ch[,3]<-format(as.numeric(as.character(p_ch[,3])),scientific=T,digits = 3)



div_alpha_ph3<-div_alpha_ph2[which(is.na(div_alpha_ph2[,1])==F),]
length(which(div_alpha_ph3[,2]=="Antibiotic"))


div_alpha_ph22<-div_alpha_ph2[which(is.na(div_alpha_ph2$alpha_div)==F),]
div_mm<-mean(div_alpha_ph22$alpha_div[which(div_alpha_ph22$medication=="No medication")])

pdf('figure_2C.pdf',width=4,height=6)
ggplot(div_alpha_ph2, aes(x=medication, y=as.numeric(as.character(alpha_div)),group=medication)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f"),outlier.alpha=0)+labs(x="",y="Shannon diversity")+
  theme(text = element_text(size=20),axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black",angle=60,hjust=1),panel.background=element_blank(),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))+ scale_x_discrete(breaks=seq(1,7,1))+
  stat_pvalue_manual(data=p_ch,label = "pval", y.position = c(5.2, 5.5, 5.8),size=5)+ylim(1.8,6.3)+#+
  geom_segment(x=0,xend=8,y=(div_mm+0.05),yend=(div_mm+0.05),col="#878787",linetype="longdash",size=1)+#+
  theme(plot.margin=margin(c(0.5,0.5,0.5,1), unit = "cm"))+scale_x_discrete(labels=c("No medication","Antibiotic","Anti-diarrheal","Anti-parasitic","Anti-fungal","Vitamins","Anti-hypertensive"))
dev.off()

##Figure 2D - wealth distribution

div_alpha_ph2<-read.csv('div_alpha_ph2_hhw.csv',row.names=1)

div_alpha_ph2$wealth <- factor(div_alpha_ph2$wealth, levels=c("Wealth index = 1","Wealth index = 2","Wealth index = 3","Wealth index = 4","Wealth index = 5"),ordered=TRUE)


p_ch<-array(NA,dim=c(4,3))

x1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$wealth=="Wealth index = 1")]
y1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$wealth=="Wealth index = 2")]
x1<-x1[which(is.na(x1)==F)]
y1<-y1[which(is.na(y1)==F)]

p_ch[1,3]<-wilcox.test(as.numeric(as.character(x1)),as.numeric(as.character(y1)))$p.value
p_ch[1,1:2]<-c("Wealth index = 1","Wealth index = 2")

x1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$wealth=="Wealth index = 1")]
y1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$wealth=="Wealth index = 3")]
x1<-x1[which(is.na(x1)==F)]
y1<-y1[which(is.na(y1)==F)]

p_ch[2,3]<-wilcox.test(as.numeric(as.character(x1)),as.numeric(as.character(y1)))$p.value
p_ch[2,1:2]<-c("Wealth index = 1","Wealth index = 3")

x1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$wealth=="Wealth index = 1")]
y1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$wealth=="Wealth index = 4")]
x1<-x1[which(is.na(x1)==F)]
y1<-y1[which(is.na(y1)==F)]

p_ch[3,3]<-wilcox.test(as.numeric(as.character(x1)),as.numeric(as.character(y1)))$p.value
p_ch[3,1:2]<-c("Wealth index = 1","Wealth index = 3")

x1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$wealth=="Wealth index = 1")]
y1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$wealth=="Wealth index = 5")]
x1<-x1[which(is.na(x1)==F)]
y1<-y1[which(is.na(y1)==F)]

p_ch[4,3]<-wilcox.test(as.numeric(as.character(x1)),as.numeric(as.character(y1)))$p.value
p_ch[4,1:2]<-c("Wealth index = 1","Wealth index = 5")

p_ch<-as.data.frame(p_ch)
colnames(p_ch)<-c("group1","group2","pval")

p_ch[,3]<-format(as.numeric(as.character(p_ch[,3])),scientific=T,digits = 3)



my_comparisons <- list( c("Wealth index = 1", "Wealth index = 2"), c("Wealth index = 1", "Wealth index = 3"), c("Wealth index = 1", "Wealth index = 4"),c("Wealth index = 1", "Wealth index = 5") )

div_alpha_ph22<-div_alpha_ph2[which(is.na(div_alpha_ph2$alpha_div)==F),]

pdf('figure_2D.pdf',width=4,height=6)
ggplot(div_alpha_ph2, aes(x=wealth, y=as.numeric(as.character(alpha_div)),group=wealth)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#b2df8a","#33a02c","#e31a1c","#fdbf6f"),outlier.alpha=0)+labs(x="",y="Shannon diversity")+
  theme(text = element_text(size=20),axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black",angle=60,hjust=1),panel.background=element_blank(),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))+ scale_x_discrete(breaks=seq(1,5,1))+
  stat_pvalue_manual(data=p_ch,label = "pval", y.position = c(5.2,5.5,5.8,6.1),size=5)+ylim(1.8,6.7)+#+
  theme(plot.margin=margin(c(0.5,0.5,0.5,1), unit = "cm"))+scale_x_discrete(labels=c("Wealth index = 1","Wealth index = 2","Wealth index = 3","Wealth index = 4","Wealth index = 5"))
dev.off()









