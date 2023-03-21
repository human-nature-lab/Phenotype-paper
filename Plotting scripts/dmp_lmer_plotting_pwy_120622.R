load('dmp_lmer_plotting_all_100922.RData')

#######################################################


##Normal plotting
##health

phen_all_use<-read.csv('phen_all_use6_all.csv',row.names = 1)

effect_size_phen_lmer<-read.csv('effect_size_phen_1_lmer_pwy.csv',row.names = 1)
pval_phen_lmer<-read.csv('pval_phen_1_lmer_pwy.csv',row.names = 1)
fdr_phen_lmer<-read.csv('fdr_phen_1_lmer_pwy.csv',row.names = 1)

effect_size_phen_lmer<-cbind(effect_size_phen_lmer,read.csv('effect_size_phen_2_lmer_pwy.csv',row.names = 1))
pval_phen_lmer<-cbind(pval_phen_lmer,read.csv('pval_phen_2_lmer_pwy.csv',row.names = 1))
fdr_phen_lmer<-cbind(fdr_phen_lmer,read.csv('fdr_phen_2_lmer_pwy.csv',row.names = 1))

effect_size_phen_lmer<-cbind(effect_size_phen_lmer,read.csv('effect_size_phen_3_lmer_pwy.csv',row.names = 1))
pval_phen_lmer<-cbind(pval_phen_lmer,read.csv('pval_phen_3_lmer_pwy.csv',row.names = 1))
fdr_phen_lmer<-cbind(fdr_phen_lmer,read.csv('fdr_phen_3_lmer_pwy.csv',row.names = 1))

effect_size_phen_lmer<-cbind(effect_size_phen_lmer,read.csv('effect_size_phen_4_lmer_pwy.csv',row.names = 1))
pval_phen_lmer<-cbind(pval_phen_lmer,read.csv('pval_phen_4_lmer_pwy.csv',row.names = 1))
fdr_phen_lmer<-cbind(fdr_phen_lmer,read.csv('fdr_phen_4_lmer_pwy.csv',row.names = 1))

effect_size_phen_lmer<-cbind(effect_size_phen_lmer,read.csv('effect_size_phen_6_lmer_pwy.csv',row.names = 1))
pval_phen_lmer<-cbind(pval_phen_lmer,read.csv('pval_phen_6_lmer_pwy.csv',row.names = 1))
fdr_phen_lmer<-cbind(fdr_phen_lmer,read.csv('fdr_phen_6_lmer_pwy.csv',row.names = 1))

effect_size_phen_lmer<-cbind(effect_size_phen_lmer,read.csv('effect_size_phen_9_lmer_pwy.csv',row.names = 1))
pval_phen_lmer<-cbind(pval_phen_lmer,read.csv('pval_phen_9_lmer_pwy.csv',row.names = 1))
fdr_phen_lmer<-cbind(fdr_phen_lmer,read.csv('fdr_phen_9_lmer_pwy.csv',row.names = 1))


effect_size_phen_2<-read.csv('effect_size_phen_2.csv',row.names=1)
pwy<-read.csv('pwy_all_use_un_norm.csv',row.names=1)#Get names after :
pwy_name<-rownames(pwy)

pwy_name_1<-array(NA,dim=length(pwy_name))
pwy_name_2<-array(NA,dim=length(pwy_name))
for(i in 1:length(pwy_name)){
  lim<-nchar(pwy_name[i])
  ind<-unlist(gregexpr(":",as.character(pwy_name[i]),fixed=T))
  if(ind>0){
    pwy_name_2[i]<-substr(as.character(pwy_name[i]),ind+2,lim)
    pwy_name_1[i]<-substr(as.character(pwy_name[i]),1,ind-1)
  }else{
    pwy_name_2[i]<-pwy_name[i]
    pwy_name_1[i]<-pwy_name[i]
  }
}


#grepl(":",pwy_name[3],fixed=T)
#unlist(gregexpr(":",as.character(pwy_name[4]),fixed=T))

#length(which(grepl("{",sml[,1],fixed=T)==T))

#effect_size_phen_lmer<-effect_size_phen_lmer[1:2285,]
#pval_phen_lmer<-pval_phen_lmer[1:2285,]
#fdr_phen_lmer<-fdr_phen_lmer[1:2285,]

#effect_size_phen_lmer<-effect_size_phen_lmer[ind_thresh_use,]#Apply it for atleast 100 ppl and above 0.0001
#pval_phen_lmer<-pval_phen_lmer[ind_thresh_use,]
#fdr_phen_lmer<-fdr_phen_lmer[ind_thresh_use,]
#mb_samp_sp_name2<-mb_samp_sp_name[ind_thresh_use]


colnames(effect_size_phen_lmer)[16:26]<-c("Diabetes (N=24)","Allergies (N=104)","Cystic fibrosis (N=18)","Heart disease (N=47)","Endocrine illness (N=16)","Renal failure (N=9)","Asthma (N=47)","Stomach illness (N=175)","Intestinal illness (N=70)","Arthritis (N=39)","MS (N=2)")#colnames(effect_size_phen_2)[14:24]
colnames(effect_size_phen_lmer)[27:30]<-c("Painkillers (N=714)","Antibiotics (N=137)","Anti-diarrheal (N=49)","Anti-parasitic (N=34)")#colnames(effect_size_phen_2)[25:28]
colnames(effect_size_phen_lmer)[85:86]<-c("Dementia","Cognitive impairment")
colnames(effect_size_phen_lmer)[152:154]<-c("Mild (GAD7)","Moderate (GAD7)","Severe (GAD7)")
colnames(effect_size_phen_lmer)[156:158]<-c("Mild (PHQ9)","Moderate (PHQ9)","Severe (PHQ9)")
#colnames(effect_size_phen_lmer)[30:33]<-c("Anti-parasitic(34)","Anti-fungal(82)","Vitamins(121)","Anti-hypertensive(72)")
#colnames(effect_size_phen_lmer)[9:11]<-c("Oxygen saturation","Hb total","Cough (1 month)")
colnames(effect_size_phen_lmer)[30:33]<-c("Anti-parasitic (N=34)","Anti-fungal (N=82)","Vitamins (N=121)","Anti-hypertensive (N=72)")
colnames(effect_size_phen_lmer)[9:11]<-c("Oxygen saturation","Hb total","Cough (N=257)")
colnames(effect_size_phen_lmer)[1]<-"Hb A1c"
colnames(effect_size_phen_lmer)[87:102]<-colnames(effect_size_phen_2)[c(99:105,108:116)]
#colnames(effect_size_phen_lmer)[164:167]<-c("Alcohol daily frequency","Alcohol 6 drinks","Cigarette usage(46)","Cigarette frequency")
colnames(effect_size_phen_lmer)[164:167]<-c("Alcohol daily frequency","Alcohol 6 drinks","Cigarette usage (N=46)","Cigarette frequency")
colnames(effect_size_phen_lmer)[209:210]<-c("Diarrhea (N=73)","Bristol stool scale")
colnames(effect_size_phen_lmer)[12:15]<-c("Excellent (N=75)","Fair (N=612)","Poor (N=153)","Very good (N=56)")
colnames(effect_size_phen_lmer)[85:86]<-c("Dementia (N=89)","Cognitive impairment (N=56)")

pval_chk_lmer<-array(0,dim=c(dim(pval_phen_lmer)[1],dim(pval_phen_lmer)[2]))
for(i in 1:dim(pval_chk_lmer)[1]){
  for(j in 1:dim(pval_chk_lmer)[2]){
    if(is.na(pval_phen_lmer[i,j])==F){
      if(pval_phen_lmer[i,j]<0.05){
        pval_chk_lmer[i,j]<-1
      }
    }
  }
}
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
#colnames(effect_size_phen)
pval_chk_lmer_all<-rowSums(pval_chk_lmer)
length(which(pval_chk_lmer_all>3))

#colnames(effect_size_phen_lmer)[34:43]
colnames(effect_size_phen_lmer)[44:83]<-c("Reserved(Agree a little)","Reserved(Agree strongly)","Reserved(Disagree a little)","Reserved(Disagree strongly)","Trusting(Agree a little)","Trusting(Agree strongly)","Trusting(Disagree a little)","Trusting(Disagree strongly)","Lazy(Agree a little)","Lazy(Agree strongly)","Lazy(Disagree a little)","Lazy(Disagree strongly)","Relaxed(Agree a little)","Relaxed(Agree strongly)","Relaxed(Disagree a little)","Relaxed(Disagree strongly)","Not creative(Agree a little)","Not creative(Agree strongly)","Not creative(Disagree a little)","Not creative(Disagree strongly)","Outgoing(Agree a little)","Outgoing(Agree strongly)","Outgoing(Disagree a little)","Outgoing(Disagree strongly)","Fault others(Agree a little)","Fault others(Agree strongly)","Fault others(Disagree a little)","Fault others(Disagree strongly)","Thorough job(Agree a little)","Thorough job(Agree strongly)","Thorough job(Disagree a little)","Thorough job(Disagree strongly)","Nervous(Agree a little)","Nervous(Agree strongly)","Nervous(Disagree a little)","Nervous(Disagree strongly)","Openess(Agree a little)","Openess(Agree strongly)","Openess(Disagree a little)","Openess(Disagree strongly)")

fdr_chk_lmer_all<-rowSums(fdr_chk_lmer[,c(1,4:6,9:10,14,13,15,12,11,209,210,16:17,19,22:25,27:34,42:43,86,85,164,166:167,152:154,156:158)])#c(1,4:11,16:17,19:20,22:25)#44:83(personality types discrete)#
length(which(fdr_chk_lmer_all>0))

effect_size_phen_sig<-effect_size_phen_lmer[which(fdr_chk_lmer_all>0),c(1,4:6,9:10,14,13,15,12,11,209,210,16:17,19,22:25,27:34,42:43,86,85,164,166:167,152:154,156:158)]
pval_phen_sig<-pval_chk_lmer[which(fdr_chk_lmer_all>0),c(1,4:6,9:10,14,13,15,12,11,209,210,16:17,19,22:25,27:34,42:43,86,85,164,166:167,152:154,156:158)]
fdr_phen_sig<-fdr_chk_lmer[which(fdr_chk_lmer_all>0),c(1,4:6,9:10,14,13,15,12,11,209,210,16:17,19,22:25,27:34,42:43,86,85,164,166:167,152:154,156:158)]

effect_size_phen_sig_plot<-effect_size_phen_sig

for(i in 1:dim(effect_size_phen_sig)[1]){
  for(j in 1:dim(effect_size_phen_sig)[2]){
    if((pval_phen_sig[i,j]==1)){
      effect_size_phen_sig_plot[i,j]<-effect_size_phen_sig[i,j]
    }else{
      effect_size_phen_sig_plot[i,j]<-0
    }
  }
}

disp_fdr<-array("",dim=c(dim(fdr_phen_sig)[1],dim(fdr_phen_sig)[2]))
for(i in 1:dim(disp_fdr)[1]){
  for(j in 1:dim(disp_fdr)[2]){
    if((fdr_phen_sig[i,j]==1)&(pval_phen_sig[i,j]==1)){
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
rownames(effect_size_phen_sig_plot)<-colnames(effect_size_phen_lmer)[c(1,4:6,9:10,14,13,15,12,11,209,210,16:17,19,22:25,27:34,42:43,86,85,164,166:167,152:154,156:158)]
#colnames(effect_size_phen_sig_plot)<-as.character(mb_samp_sp_name2[which(fdr_chk_lmer_all>0)])
colnames(effect_size_phen_sig_plot)<-as.character(pwy_name_2[which(fdr_chk_lmer_all>0)])

paletteLength<-13
myColor <- colorRampPalette(c("#cb181d","white", "#33a02c"))(paletteLength)#myColor <- colorRampPalette(c("Orange", "white", "Darkblue"))(paletteLength)
myBreaks <- c(seq(min(effect_size_phen_sig_plot),0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(effect_size_phen_sig_plot)/paletteLength, max(effect_size_phen_sig_plot), length.out=floor(paletteLength/2)))

myBreaks
myBreaks[7]<-as.numeric(-0.001)
myBreaks[9]<-as.numeric(0.001)
myBreaks

library(pheatmap)

#pheatmap(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,cluster_rows = F,fontsize_number = 13,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr),number_color = "grey30")
pheatmap(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

#
#effect_size_phen_lmer2<-effect_size_phen_lmer[,c(1,4:6,9:10,14,13,15,12,11,209,210,16:17,19,22:25,27:34,42:43,86,85,164,166:167,152:154,156:158)]
#pval_phen_lmer2<-pval_phen_lmer[,c(1,4:6,9:10,14,13,15,12,11,209,210,16:17,19,22:25,27:34,42:43,86,85,164,166:167,152:154,156:158)]
#fdr_phen_lmer2<-fdr_phen_lmer[,c(1,4:6,9:10,14,13,15,12,11,209,210,16:17,19,22:25,27:34,42:43,86,85,164,166:167,152:154,156:158)]

#Pathway name
colnames(effect_size_phen_sig_plot)<-as.character(pwy_name_1[which(fdr_chk_lmer_all>0)])
pheatmap(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

pheatmap(effect_size_phen_sig_plot[,1:84],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[1:84,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
pheatmap(effect_size_phen_sig_plot[,85:167],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[85:167,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

colnames(effect_size_phen_sig_plot)<-as.character(pwy_name[which(fdr_chk_lmer_all>0)])
write.csv(effect_size_phen_sig_plot,'effect_size_phen_sig_plot_pwy_health.csv')
write.csv(disp_fdr,'disp_fdr_pwy_health.csv')

for(i in 1:dim(disp_fdr)[1]){
  if(i==1){
    sp_c<-length(which(disp_fdr[i,]!=""))
  }else{
    sp_c<-rbind(sp_c,length(which(disp_fdr[i,]!="")))
  }
}
colnames(effect_size_phen_sig_plot)[which(sp_c==max(sp_c))]


##########
##Food and animals

effect_size_phen_lmer<-read.csv('effect_size_phen_5_lmer_pwy.csv',row.names = 1)
pval_phen_lmer<-read.csv('pval_phen_5_lmer_pwy.csv',row.names = 1)
fdr_phen_lmer<-read.csv('fdr_phen_5_lmer_pwy.csv',row.names = 1)

pwy<-read.csv('pwy_all_use_un_norm.csv',row.names=1)#Get names after :
pwy_name<-rownames(pwy)

pwy_name_1<-array(NA,dim=length(pwy_name))
pwy_name_2<-array(NA,dim=length(pwy_name))
for(i in 1:length(pwy_name)){
  lim<-nchar(pwy_name[i])
  ind<-unlist(gregexpr(":",as.character(pwy_name[i]),fixed=T))
  if(ind>0){
    pwy_name_2[i]<-substr(as.character(pwy_name[i]),ind+2,lim)
    pwy_name_1[i]<-substr(as.character(pwy_name[i]),1,ind-1)
  }else{
    pwy_name_2[i]<-pwy_name[i]
    pwy_name_1[i]<-pwy_name[i]
  }
}

#effect_size_phen_lmer<-effect_size_phen_lmer[1:2285,]
#pval_phen_lmer<-pval_phen_lmer[1:2285,]
#fdr_phen_lmer<-fdr_phen_lmer[1:2285,]

#effect_size_phen_lmer<-effect_size_phen_lmer[ind_thresh_use,]#Apply it for atleast 100 ppl and above 0.0001
#pval_phen_lmer<-pval_phen_lmer[ind_thresh_use,]
#fdr_phen_lmer<-fdr_phen_lmer[ind_thresh_use,]

#colnames(effect_size_phen_lmer)[1:26]<-colnames(effect_size_phen_2)[45:70]
colnames(effect_size_phen_lmer)[1:26]<-c("Cat (N=645)","Dog (N=944)","Parakeet (N=102)","Rabbit (N=120)","Horse (N=314)","Mice (N=585)","None pet (N=88)","Cow (N=397)","Goat (N=40)","Pig (N=208)","Chicken (N=1094)","Duck (N=478)","Turkey (N=322)","Sheep (N=33)","Geese (N=68)","None farm (N=69)","Bat (N=410)","Lizard (N=423)","Monkey (N=12)","Snake (N=437)","Bird (N=794)","Possum (N=467)","Rat (N=457)","Squirrel (N=480)","Other Wild (N=1)","None wild (N=243)")
colnames(effect_size_phen_lmer)[38:41]<-c("Natural juice","Chicken","Beef/pork","Ham/sausages/hotdog")
colnames(effect_size_phen_lmer)[33]<-"Cream/butter"
colnames(effect_size_phen_lmer)[46]<-"Diet diversity score"

pval_chk_lmer<-array(0,dim=c(dim(fdr_phen_lmer)[1],dim(pval_phen_lmer)[2]))
for(i in 1:dim(pval_chk_lmer)[1]){
  for(j in 1:dim(pval_chk_lmer)[2]){
    if(is.na(pval_phen_lmer[i,j])==F){
      if(pval_phen_lmer[i,j]<0.05){
        pval_chk_lmer[i,j]<-1
      }
    }
  }
}
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
#colnames(effect_size_phen)
pval_chk_lmer_all<-rowSums(pval_chk_lmer)
length(which(pval_chk_lmer_all>3))

fdr_chk_lmer_all<-rowSums(fdr_chk_lmer[,c(1:3,7:15,5,4,16,6,17:18,20:24,26:42,46)])#c(1,4:11,16:17,19:20,22:25)#44:83(personality types discrete)#
length(which(fdr_chk_lmer_all>0))#>5

effect_size_phen_sig<-effect_size_phen_lmer[which(fdr_chk_lmer_all>0),c(1:3,7:15,5,4,16,6,17:18,20:24,26:42,46)]
pval_phen_sig<-pval_chk_lmer[which(fdr_chk_lmer_all>0),c(1:3,7:15,5,4,16,6,17:18,20:24,26:42,46)]
fdr_phen_sig<-fdr_chk_lmer[which(fdr_chk_lmer_all>0),c(1:3,7:15,5,4,16,6,17:18,20:24,26:42,46)]

effect_size_phen_sig_plot<-effect_size_phen_sig

for(i in 1:dim(effect_size_phen_sig)[1]){
  for(j in 1:dim(effect_size_phen_sig)[2]){
    if((pval_phen_sig[i,j]==1)){
      effect_size_phen_sig_plot[i,j]<-effect_size_phen_sig[i,j]
    }else{
      effect_size_phen_sig_plot[i,j]<-0
    }
  }
}

disp_fdr<-array("",dim=c(dim(fdr_phen_sig)[1],dim(fdr_phen_sig)[2]))
for(i in 1:dim(disp_fdr)[1]){
  for(j in 1:dim(disp_fdr)[2]){
    if((fdr_phen_sig[i,j]==1)&(pval_phen_sig[i,j]==1)){
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
rownames(effect_size_phen_sig_plot)<-colnames(effect_size_phen_lmer)[c(1:3,7:15,5,4,16,6,17:18,20:24,26:42,46)]
colnames(effect_size_phen_sig_plot)<-as.character(pwy_name_2[which(fdr_chk_lmer_all>0)])
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
pheatmap(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

colnames(effect_size_phen_sig_plot)<-as.character(pwy_name_1[which(fdr_chk_lmer_all>0)])

pheatmap(effect_size_phen_sig_plot[,1:81],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[1:81,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
pheatmap(effect_size_phen_sig_plot[,82:162],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[82:162,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
pheatmap(effect_size_phen_sig_plot[,163:243],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[163:243,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

#effect_size_phen_lmer2<-cbind(effect_size_phen_lmer2,effect_size_phen_lmer[,c(1:3,7:15,5,4,16,6,17:18,20:24,26:42,46)])
#pval_phen_lmer2<-cbind(pval_phen_lmer2,pval_phen_lmer[,c(1:3,7:15,5,4,16,6,17:18,20:24,26:42,46)])
#fdr_phen_lmer2<-cbind(fdr_phen_lmer2,fdr_phen_lmer[,c(1:3,7:15,5,4,16,6,17:18,20:24,26:42,46)])

colnames(effect_size_phen_sig_plot)<-as.character(pwy_name[which(fdr_chk_lmer_all>0)])
write.csv(effect_size_phen_sig_plot,'effect_size_phen_sig_plot_pwy_food_animal.csv')
write.csv(disp_fdr,'disp_fdr_pwy_food_animal.csv')

for(i in 1:dim(disp_fdr)[1]){
  if(i==1){
    sp_c<-length(which(disp_fdr[i,]!=""))
  }else{
    sp_c<-rbind(sp_c,length(which(disp_fdr[i,]!="")))
  }
}
colnames(effect_size_phen_sig_plot)[which(sp_c==max(sp_c))]

colnames(effect_size_phen_sig_plot)[which(sp_c==4)]


#######
##Socio-economic factors

effect_size_phen_lmer<-read.csv('effect_size_phen_6_lmer_pwy.csv',row.names = 1)
pval_phen_lmer<-read.csv('pval_phen_6_lmer_pwy.csv',row.names = 1)
fdr_phen_lmer<-read.csv('fdr_phen_6_lmer_pwy.csv',row.names = 1)

effect_size_phen_lmer<-cbind(effect_size_phen_lmer,read.csv('effect_size_phen_7_lmer_pwy.csv',row.names = 1))
pval_phen_lmer<-cbind(pval_phen_lmer,read.csv('pval_phen_7_lmer_pwy.csv',row.names = 1))
fdr_phen_lmer<-cbind(fdr_phen_lmer,read.csv('fdr_phen_7_lmer_pwy.csv',row.names = 1))

effect_size_phen_lmer<-cbind(effect_size_phen_lmer,read.csv('effect_size_phen_9_lmer_pwy.csv',row.names = 1))
pval_phen_lmer<-cbind(pval_phen_lmer,read.csv('pval_phen_9_lmer_pwy.csv',row.names = 1))
fdr_phen_lmer<-cbind(fdr_phen_lmer,read.csv('fdr_phen_9_lmer_pwy.csv',row.names = 1))

colnames(effect_size_phen_lmer)[1:2]<-c("GAD7(binary)","PHQ9(binary)")
colnames(effect_size_phen_lmer)[3:9]<-c("Travel","Monthly expenditure","Alcohol frequency","Alcohol daily","Alcohol(6 drinks)","Cigarette use","Cigarette frequency")
colnames(effect_size_phen_lmer)[11]<-"Risk taking"
colnames(effect_size_phen_lmer)[22:29]<-c("Friend ties(same building)","Friend ties(different building)","Betweeness(friendship)","Transitivity(friendship)","Risky interaction","Family ties(same building)","Family ties(different building)","Transitivity(familial)")
#colnames(effect_size_phen_lmer)[27:30]<-colnames(effect_size_phen_2)[25:28]
colnames(effect_size_phen_lmer)[33:34]<-c("Clustering coefficient(all ties)","Betweeness(all ties)")
colnames(effect_size_phen_lmer)[36:39]<-c("Partner live duration","Number of partners","Partner live age","Living with partner (N=375)")
#colnames(effect_size_phen_lmer)[42:95]<-colnames(effect_size_phen_2)[c(138:144,134:137,145:187)]
colnames(effect_size_phen_lmer)[42:95]<-c("Electricity (N=1069)","Radio (N=519)","TV (N=565)","Cellphone (N=927)","No cellphone (N=12)","Refrigerator (N=397)","None electronics (N=38)","Chimney (N=959)","No chimney (N=151)","Stove (N=61)","No stove (N=1)","Wood (N=1109)","Gas(fuel) [N=46]","Electricity fuel (N=14)","Kerosene (N=2)","None fuel (N=1)","Separate kitchen (N=1106)","Cement floor (N=552)","Earth/Sand floor (N=507)","Ceramic floor (N=88)","Tiles floor (N=19)","Mud bricks floor (N=1)","Wood floor (N=0)","Other floor (N=5)","Wooden windows (N=1035)","Glass windows (N=51)","Metal windows (N=6)","Unfinished windows (N=35)","No windows (N=45)","Clay/mud walls(N=693)","Clay brick walls (N=3)","Cement walls (N=452)","Cane/palm/trunks walls (N=0)","Wood unpolished walls (N=21)","Wood polished walls (N=0)","Discarded materials walls (N=1)","No walls (N=0)","Other walls (N=2)","Plastic roof (N=7)","Metal roof (N=982)","Clay roof (N=120)","Thatch/palm roof (N=10)","Concrete roof (N=27)","Wood roof (N=22)","Other roof(N=4)","Sleeping rooms","Spring(protected) [N=917]","Spring(unprotected) [N=19]","Tube well (N=72)","Dug well(protected) [N=78]","Dug well(unprotected) [N=26]","Surface water (N=10)","Bottled water (N=15)","Other water (N=6)")
colnames(effect_size_phen_lmer)[13:21]<-c("Education(1st grade)","Education(2nd grade)","Education(3rd grade)","Education(4th grade)","Education(5th grade)","Education(6th grade)","Education(More than secondary)","Education(Secondary)","Education(Some secondary)")
colnames(effect_size_phen_lmer)[30:32]<-c("Betweeness(familial)","Degree(all ties)","%kin(upto 3rd degree)")
colnames(effect_size_phen_lmer)[35]<-"k-cycle centrality"
colnames(effect_size_phen_lmer)[40:41]<-c("Washing hands (N=1130)","Distance to village center")
colnames(effect_size_phen_lmer)[102:104]<-c("Primary (N=373)","Middle (N=311)","Secondary (N=113)")
colnames(effect_size_phen_lmer)[96:101]<-c("Household size","Household wealth index","Distance to main road","Distance to health center","Number of churches","Deforestation (%)")
colnames(effect_size_phen_lmer)[32]<-"Kin percentage (upto third degree)"

pval_chk_lmer<-array(0,dim=c(dim(fdr_phen_lmer)[1],dim(pval_phen_lmer)[2]))
for(i in 1:dim(pval_chk_lmer)[1]){
  for(j in 1:dim(pval_chk_lmer)[2]){
    if(is.na(pval_phen_lmer[i,j])==F){
      if(pval_phen_lmer[i,j]<0.05){
        pval_chk_lmer[i,j]<-1
      }
    }
  }
}
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
colnames(effect_size_phen_lmer)
pval_chk_lmer_all<-rowSums(pval_chk_lmer)
length(which(pval_chk_lmer_all>3))


fdr_chk_lmer_all<-rowSums(fdr_chk_lmer[,c(36:39,102:104,22:25,27:31,33:34,32,10:11,40:41,98:99,101,3:4,96:97,47,48,51,54,58,60,67,71,73,81,87:88,91,92,94)])
length(which(fdr_chk_lmer_all>0))

effect_size_phen_sig<-effect_size_phen_lmer[which(fdr_chk_lmer_all>0),c(36:39,102:104,22:25,27:31,33:34,32,10:11,40:41,98:99,101,3:4,96:97,47,48,51,54,58,60,67,71,73,81,87:88,91,92,94)]
pval_phen_sig<-pval_chk_lmer[which(fdr_chk_lmer_all>0),c(36:39,102:104,22:25,27:31,33:34,32,10:11,40:41,98:99,101,3:4,96:97,47,48,51,54,58,60,67,71,73,81,87:88,91,92,94)]
fdr_phen_sig<-fdr_chk_lmer[which(fdr_chk_lmer_all>0),c(36:39,102:104,22:25,27:31,33:34,32,10:11,40:41,98:99,101,3:4,96:97,47,48,51,54,58,60,67,71,73,81,87:88,91,92,94)]

effect_size_phen_sig_plot<-effect_size_phen_sig

for(i in 1:dim(effect_size_phen_sig)[1]){
  for(j in 1:dim(effect_size_phen_sig)[2]){
    if((pval_phen_sig[i,j]==1)){
      effect_size_phen_sig_plot[i,j]<-effect_size_phen_sig[i,j]
    }else{
      effect_size_phen_sig_plot[i,j]<-0
    }
  }
}

disp_fdr<-array("",dim=c(dim(fdr_phen_sig)[1],dim(fdr_phen_sig)[2]))
for(i in 1:dim(disp_fdr)[1]){
  for(j in 1:dim(disp_fdr)[2]){
    if((fdr_phen_sig[i,j]==1)&(pval_phen_sig[i,j]==1)){
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
rownames(effect_size_phen_sig_plot)<-colnames(effect_size_phen_lmer)[c(36:39,102:104,22:25,27:31,33:34,32,10:11,40:41,98:99,101,3:4,96:97,47,48,51,54,58,60,67,71,73,81,87:88,91,92,94)]
colnames(effect_size_phen_sig_plot)<-as.character(pwy_name_2[which(fdr_chk_lmer_all>0)])
paletteLength<-13
myColor <- colorRampPalette(c("#cb181d","white", "#33a02c"))(paletteLength)#myColor <- colorRampPalette(c("Orange", "white", "Darkblue"))(paletteLength)
myBreaks <- c(seq(min(effect_size_phen_sig_plot),0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(effect_size_phen_sig_plot)/paletteLength, max(effect_size_phen_sig_plot), length.out=floor(paletteLength/2)))

myBreaks
myBreaks[7]<-as.numeric(-0.001)
myBreaks[9]<-as.numeric(0.001)
myBreaks

library(pheatmap)

#pheatmap(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,cluster_rows = F,fontsize_number = 13,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr),number_color = "grey30")
pheatmap(effect_size_phen_sig_plot[,1:80],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[1:80,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
pheatmap(effect_size_phen_sig_plot[,81:160],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[81:160,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
pheatmap(effect_size_phen_sig_plot[,161:227],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[161:227,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

#effect_size_phen_lmer2<-cbind(effect_size_phen_lmer2,effect_size_phen_lmer[,c(36:39,102:104,22:25,27:31,33:34,32,10:11,40:41,98:99,101,3:4,96:97,47,48,51,54,58,60,67,71,73,81,87:88,91,92,94)])
#pval_phen_lmer2<-cbind(pval_phen_lmer2,pval_phen_lmer[,c(36:39,102:104,22:25,27:31,33:34,32,10:11,40:41,98:99,101,3:4,96:97,47,48,51,54,58,60,67,71,73,81,87:88,91,92,94)])
#fdr_phen_lmer2<-cbind(fdr_phen_lmer2,fdr_phen_lmer[,c(36:39,102:104,22:25,27:31,33:34,32,10:11,40:41,98:99,101,3:4,96:97,47,48,51,54,58,60,67,71,73,81,87:88,91,92,94)])

#rownames(effect_size_phen_lmer2)<-pwy_name#full name

#write.csv(effect_size_phen_lmer2,'esz_pwy.csv')
#write.csv(pval_phen_lmer2,'pval_pwy.csv')
#write.csv(fdr_phen_lmer2,'fdr_pwy.csv')


#Pathway name
colnames(effect_size_phen_sig_plot)<-as.character(pwy_name_1[which(fdr_chk_lmer_all>0)])
pheatmap(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)


colnames(effect_size_phen_sig_plot)<-as.character(pwy_name[which(fdr_chk_lmer_all>0)])
write.csv(effect_size_phen_sig_plot,'effect_size_phen_sig_plot_pwy_socioeco.csv')
write.csv(disp_fdr,'disp_fdr_pwy_socioeco.csv')

for(i in 1:dim(disp_fdr)[1]){
  if(i==1){
    sp_c<-length(which(disp_fdr[i,]!=""))
  }else{
    sp_c<-rbind(sp_c,length(which(disp_fdr[i,]!="")))
  }
}
colnames(effect_size_phen_sig_plot)[which(sp_c==max(sp_c))]


