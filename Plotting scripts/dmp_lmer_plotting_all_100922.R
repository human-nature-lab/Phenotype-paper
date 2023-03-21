
##Plotting for NEW LMER model

##First load dmp_analysis_all_meta4_08-15-22.RData

#RV -- ph on left, sp on right, species explains phenotype

mb_samp_sp_3<-read.csv('MB_1187_2285.csv',row.names=1)

#mb_samp_sp_3<-mb_samp_sp2[ind_t,]
ind_thresh_use<-0
for(i in 1:dim(mb_samp_sp_3)[1]){
  t100<-length(which(mb_samp_sp_3[i,]>0))
  if((t100>150)&((sum(mb_samp_sp_3[i,])/1187)>0.0001)){#was 0.0001 for 910 species
    ind_thresh_use<-rbind(ind_thresh_use,i)
  }
}


#write.csv(ind_thresh_use,'ind_thresh_use2_00001_150.csv')
#write.csv(ind_thresh_use,'ind_thresh_use3_0001.csv')
#write.csv(ind_thresh_use,'ind_thresh_use4_001.csv')
ind_thresh_use2<-read.csv('ind_thresh_use2_00001_150.csv')
ind_thresh_use3<-read.csv('ind_thresh_use3_0001.csv')
ind_thresh_use4<-read.csv('ind_thresh_use4_001.csv')

mb_samp_sp_name_use2<-mb_samp_sp_name[ind_thresh_use2[,2]]
mb_samp_sp_name_use3<-mb_samp_sp_name[ind_thresh_use3[,2]]
mb_samp_sp_name_use4<-mb_samp_sp_name[ind_thresh_use4[,2]]

#Use ind_thresh_use Get which species are present in 976 strain list and use it in anpan
#write.csv(mb_samp_sp_name2,'species_sig_name.csv')

#file_list<-list.files(path=".")
#ik<-1
#ik2<-1
#sccm<-0

#for(ik in 1:length(file_list)){
#  for(ik2 in 1:length(mb_samp_sp_name2)){
#  if(grepl("t__",as.character(file_list[ik]),fixed=T)){
#    ind_s1<-unlist(gregexpr("t__",as.character(file_list[ik]),fixed=T))
#    ind_s2<-unlist(gregexpr(".StrainPhlAn4",as.character(file_list[ik]),fixed=T))
#    sp_comp<-substr(file_list[ik],ind_s1+3,ind_s2-1)
#    
#    ind_s3<-unlist(gregexpr("(",as.character(mb_samp_sp_name2[ik2]),fixed=T))
#    ind_s4<-unlist(gregexpr(")",as.character(mb_samp_sp_name2[ik2]),fixed=T))
#    sp_comp2<-substr(mb_samp_sp_name2[ik2],ind_s3+1,ind_s4-1)
#    
#    if(sp_comp==sp_comp2){
#      sccm<-sccm+1
#      if(sccm==1){
#        strain_sp_list<-ik
#        strain_sp_list_sp<-ind_thresh_use[ik2+1]
#      }else{
#        strain_sp_list<-rbind(strain_sp_list,ik)
#        strain_sp_list_sp<-rbind(strain_sp_list_sp,ind_thresh_use[ik2+1])
#      }
#    }
#    
#  }
#  }
#}

#write.csv(strain_sp_list[,1],'strain_sp_list.csv')
#write.csv(strain_sp_list_sp[,1],'strain_sp_list_sp.csv')



effect_size_phen_lmer_rv<-read.csv('effect_size_phen_1_lmer_rv.csv',row.names = 1)
pval_phen_lmer_rv<-read.csv('pval_phen_1_lmer_rv.csv',row.names = 1)
fdr_phen_lmer_rv<-read.csv('fdr_phen_1_lmer_rv.csv',row.names = 1)

effect_size_phen_lmer_rv<-effect_size_phen_lmer_rv[1:2285,]
pval_phen_lmer_rv<-pval_phen_lmer_rv[1:2285,]
fdr_phen_lmer_rv<-fdr_phen_lmer_rv[1:2285,]

effect_size_phen_lmer_rv<-effect_size_phen_lmer_rv[ind_thresh_use,]#Apply it for atleast 100 ppl and above 0.0001
pval_phen_lmer_rv<-pval_phen_lmer_rv[ind_thresh_use,]
fdr_phen_lmer_rv<-fdr_phen_lmer_rv[ind_thresh_use,]


colnames(effect_size_phen_lmer_rv)[12:22]<-colnames(effect_size_phen_2)[14:24]

pval_chk_lmer_rv<-array(0,dim=c(dim(pval_phen_lmer_rv)[1],dim(pval_phen_lmer_rv)[2]))
for(i in 1:dim(pval_chk_lmer_rv)[1]){
  for(j in 1:dim(pval_chk_lmer_rv)[2]){
    if(is.na(pval_phen_lmer_rv[i,j])==F){
      if(pval_phen_lmer_rv[i,j]<0.05){
        pval_chk_lmer_rv[i,j]<-1
      }
    }
  }
}
fdr_chk_lmer_rv<-array(0,dim=c(dim(fdr_phen_lmer_rv)[1],dim(fdr_phen_lmer_rv)[2]))
for(i in 1:dim(fdr_chk_lmer_rv)[1]){
  for(j in 1:dim(fdr_chk_lmer_rv)[2]){
    if(is.na(fdr_phen_lmer_rv[i,j])==F){
      if(fdr_phen_lmer_rv[i,j]<0.05){
        fdr_chk_lmer_rv[i,j]<-1
      }
    }
  }
}
colnames(effect_size_phen)
pval_chk_lmer_rv_all<-rowSums(pval_chk_lmer_rv)
length(which(pval_chk_lmer_rv_all>3))


fdr_chk_lmer_rv_all<-rowSums(fdr_chk_lmer_rv[,c(1,4:13,15:16,18:21)])
length(which(fdr_chk_lmer_rv_all>3))

effect_size_phen_sig<-effect_size_phen_lmer_rv[which(fdr_chk_lmer_rv_all>3),c(1,4:13,15:16,18:21)]
pval_phen_sig<-pval_chk_lmer_rv[which(fdr_chk_lmer_rv_all>3),c(1,4:13,15:16,18:21)]
fdr_phen_sig<-fdr_chk_lmer_rv[which(fdr_chk_lmer_rv_all>3),c(1,4:13,15:16,18:21)]

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
rownames(effect_size_phen_sig_plot)<-colnames(effect_size_phen_lmer_rv)[c(1,4:13,15:16,18:21)]
colnames(effect_size_phen_sig_plot)<-as.character(mb_samp_sp_name[which(fdr_chk_lmer_rv_all>3)])
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

##################
###Same model as DMP but with Random effect (village)

effect_size_phen_lmer<-read.csv('effect_size_phen_1_lmer.csv',row.names = 1)
pval_phen_lmer<-read.csv('pval_phen_1_lmer.csv',row.names = 1)
fdr_phen_lmer<-read.csv('fdr_phen_1_lmer.csv',row.names = 1)

effect_size_phen_lmer<-cbind(effect_size_phen_lmer,read.csv('effect_size_phen_2_lmer.csv',row.names = 1))
pval_phen_lmer<-cbind(pval_phen_lmer,read.csv('pval_phen_2_lmer.csv',row.names = 1))
fdr_phen_lmer<-cbind(fdr_phen_lmer,read.csv('fdr_phen_2_lmer.csv',row.names = 1))

effect_size_phen_lmer<-effect_size_phen_lmer[1:2285,]
pval_phen_lmer<-pval_phen_lmer[1:2285,]
fdr_phen_lmer<-fdr_phen_lmer[1:2285,]

effect_size_phen_lmer<-effect_size_phen_lmer[ind_thresh_use,]#Apply it for atleast 100 ppl and above 0.0001
pval_phen_lmer<-pval_phen_lmer[ind_thresh_use,]
fdr_phen_lmer<-fdr_phen_lmer[ind_thresh_use,]


colnames(effect_size_phen_lmer)[16:26]<-colnames(effect_size_phen_2)[14:24]
colnames(effect_size_phen_lmer)[27:30]<-colnames(effect_size_phen_2)[25:28]


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
colnames(effect_size_phen)
pval_chk_lmer_all<-rowSums(pval_chk_lmer)
length(which(pval_chk_lmer_all>3))


fdr_chk_lmer_all<-rowSums(fdr_chk_lmer[,c(1,4:11,16:17,19:20,22:25)])
length(which(fdr_chk_lmer_all>1))

effect_size_phen_sig<-effect_size_phen_lmer[which(fdr_chk_lmer_all>1),c(1,4:11,16:17,19:20,22:25)]
pval_phen_sig<-pval_chk_lmer[which(fdr_chk_lmer_all>1),c(1,4:11,16:17,19:20,22:25)]
fdr_phen_sig<-fdr_chk_lmer[which(fdr_chk_lmer_all>1),c(1,4:11,16:17,19:20,22:25)]

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
rownames(effect_size_phen_sig_plot)<-colnames(effect_size_phen_lmer)[c(1,4:11,16:17,19:20,22:25)]
colnames(effect_size_phen_sig_plot)<-as.character(mb_samp_sp_name[which(fdr_chk_lmer_all>1)])
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

##Mod 2 --- sampling date as Random effect, and removing batch effect

effect_size_phen_lmer_mod2<-read.csv('effect_size_phen_1_lmer_mod2.csv',row.names = 1)
pval_phen_lmer_mod2<-read.csv('pval_phen_1_lmer_mod2.csv',row.names = 1)
fdr_phen_lmer_mod2<-read.csv('fdr_phen_1_lmer_mod2.csv',row.names = 1)

effect_size_phen_lmer_mod2<-effect_size_phen_lmer_mod2[1:2285,]
pval_phen_lmer_mod2<-pval_phen_lmer_mod2[1:2285,]
fdr_phen_lmer_mod2<-fdr_phen_lmer_mod2[1:2285,]

effect_size_phen_lmer_mod2<-cbind(effect_size_phen_lmer_mod2,read.csv('effect_size_phen_2_lmer_mod2.csv',row.names = 1))
pval_phen_lmer_mod2<-cbind(pval_phen_lmer_mod2,read.csv('pval_phen_2_lmer_mod2.csv',row.names = 1))
fdr_phen_lmer_mod2<-cbind(fdr_phen_lmer_mod2,read.csv('fdr_phen_2_lmer_mod2.csv',row.names = 1))



effect_size_phen_lmer_mod2<-effect_size_phen_lmer_mod2[ind_thresh_use,]#Apply it for atleast 100 ppl and above 0.0001
pval_phen_lmer_mod2<-pval_phen_lmer_mod2[ind_thresh_use,]
fdr_phen_lmer_mod2<-fdr_phen_lmer_mod2[ind_thresh_use,]


colnames(effect_size_phen_lmer_mod2)[16:26]<-colnames(effect_size_phen_2)[14:24]
colnames(effect_size_phen_lmer_mod2)[27:30]<-colnames(effect_size_phen_2)[25:28]


pval_chk_lmer_mod2<-array(0,dim=c(dim(pval_phen_lmer_mod2)[1],dim(pval_phen_lmer_mod2)[2]))
for(i in 1:dim(pval_chk_lmer_mod2)[1]){
  for(j in 1:dim(pval_chk_lmer_mod2)[2]){
    if(is.na(pval_phen_lmer_mod2[i,j])==F){
      if(pval_phen_lmer_mod2[i,j]<0.05){
        pval_chk_lmer_mod2[i,j]<-1
      }
    }
  }
}
fdr_chk_lmer_mod2<-array(0,dim=c(dim(fdr_phen_lmer_mod2)[1],dim(fdr_phen_lmer_mod2)[2]))
for(i in 1:dim(fdr_chk_lmer_mod2)[1]){
  for(j in 1:dim(fdr_chk_lmer_mod2)[2]){
    if(is.na(fdr_phen_lmer_mod2[i,j])==F){
      if(fdr_phen_lmer_mod2[i,j]<0.05){
        fdr_chk_lmer_mod2[i,j]<-1
      }
    }
  }
}
#colnames(effect_size_phen)
pval_chk_lmer_all<-rowSums(pval_chk_lmer_mod2)
length(which(pval_chk_lmer_all>3))


fdr_chk_lmer_all<-rowSums(fdr_chk_lmer_mod2[,c(1,4:11,16:17,19:20,22:25)])
length(which(fdr_chk_lmer_all>2))

effect_size_phen_sig<-effect_size_phen_lmer_mod2[which(fdr_chk_lmer_all>2),c(1,4:11,16:17,19:20,22:25)]
pval_phen_sig<-pval_chk_lmer_mod2[which(fdr_chk_lmer_all>2),c(1,4:11,16:17,19:20,22:25)]
fdr_phen_sig<-fdr_chk_lmer_mod2[which(fdr_chk_lmer_all>2),c(1,4:11,16:17,19:20,22:25)]

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
rownames(effect_size_phen_sig_plot)<-colnames(effect_size_phen_lmer_mod2)[c(1,4:11,16:17,19:20,22:25)]
colnames(effect_size_phen_sig_plot)<-as.character(mb_samp_sp_name[which(fdr_chk_lmer_all>2)])
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




#######################################################


##Normal plotting
##health

effect_size_phen_lmer<-read.csv('effect_size_phen_1_lmer.csv',row.names = 1)
pval_phen_lmer<-read.csv('pval_phen_1_lmer.csv',row.names = 1)
fdr_phen_lmer<-read.csv('fdr_phen_1_lmer.csv',row.names = 1)

effect_size_phen_lmer<-cbind(effect_size_phen_lmer,read.csv('effect_size_phen_2_lmer.csv',row.names = 1))
pval_phen_lmer<-cbind(pval_phen_lmer,read.csv('pval_phen_2_lmer.csv',row.names = 1))
fdr_phen_lmer<-cbind(fdr_phen_lmer,read.csv('fdr_phen_2_lmer.csv',row.names = 1))

effect_size_phen_lmer<-cbind(effect_size_phen_lmer,read.csv('effect_size_phen_3_lmer.csv',row.names = 1))
pval_phen_lmer<-cbind(pval_phen_lmer,read.csv('pval_phen_3_lmer.csv',row.names = 1))
fdr_phen_lmer<-cbind(fdr_phen_lmer,read.csv('fdr_phen_3_lmer.csv',row.names = 1))

effect_size_phen_lmer<-cbind(effect_size_phen_lmer,read.csv('effect_size_phen_4_lmer.csv',row.names = 1))
pval_phen_lmer<-cbind(pval_phen_lmer,read.csv('pval_phen_4_lmer.csv',row.names = 1))
fdr_phen_lmer<-cbind(fdr_phen_lmer,read.csv('fdr_phen_4_lmer.csv',row.names = 1))

effect_size_phen_lmer<-cbind(effect_size_phen_lmer,read.csv('effect_size_phen_6_lmer.csv',row.names = 1))
pval_phen_lmer<-cbind(pval_phen_lmer,read.csv('pval_phen_6_lmer.csv',row.names = 1))
fdr_phen_lmer<-cbind(fdr_phen_lmer,read.csv('fdr_phen_6_lmer.csv',row.names = 1))

effect_size_phen_lmer<-cbind(effect_size_phen_lmer,read.csv('effect_size_phen_9_lmer.csv',row.names = 1))
pval_phen_lmer<-cbind(pval_phen_lmer,read.csv('pval_phen_9_lmer.csv',row.names = 1))
fdr_phen_lmer<-cbind(fdr_phen_lmer,read.csv('fdr_phen_9_lmer.csv',row.names = 1))

effect_size_phen_lmer<-effect_size_phen_lmer[1:2285,]
pval_phen_lmer<-pval_phen_lmer[1:2285,]
fdr_phen_lmer<-fdr_phen_lmer[1:2285,]

effect_size_phen_lmer<-effect_size_phen_lmer[ind_thresh_use2[,2],]#[ind_thresh_use,]#Apply it for atleast 100 ppl and above 0.0001
pval_phen_lmer<-pval_phen_lmer[ind_thresh_use2[,2],]#[ind_thresh_use,]
fdr_phen_lmer<-fdr_phen_lmer[ind_thresh_use2[,2],]#[ind_thresh_use,]
mb_samp_sp_name2<-mb_samp_sp_name[ind_thresh_use2[,2]]#[ind_thresh_use]


#colnames(effect_size_phen_lmer)[16:26]<-colnames(effect_size_phen_2)[14:24]
#colnames(effect_size_phen_lmer)[27:30]<-colnames(effect_size_phen_2)[25:28]
colnames(effect_size_phen_lmer)[16:26]<-c("Diabetes (N=20)","Allergies (N=104)","Cystic fibrosis (N=18)","Heart disease (N=47)","Endocrine illness (N=16)","Renal failure (N=9)","Asthma (N=47)","Stomach illness (N=175)","Intestinal illness (N=70)","Arthritis (N=39)","MS (N=2)")#colnames(effect_size_phen_2)[14:24]
colnames(effect_size_phen_lmer)[27:30]<-c("Painkillers (N=714)","Antibiotics (N=137)","Anti-diarrheal (N=49)","Anti-parasitic (N=34)")#colnames(effect_size_phen_2)[25:28]
colnames(effect_size_phen_lmer)[85:86]<-c("Dementia","Cognitive impairment")
colnames(effect_size_phen_lmer)[152:154]<-c("Mild","Moderate","Severe")#Add GAD7 in category label
colnames(effect_size_phen_lmer)[156:158]<-c("Mild","Moderate","Severe")#Add PHQ9 in category label
#colnames(effect_size_phen_lmer)[30:33]<-c("Anti-parasitic(34)","Anti-fungal(82)","Vitamins(121)","Anti-hypertensive(72)")
#colnames(effect_size_phen_lmer)[9:11]<-c("Oxygen saturation","Hb total","Cough (1 month)")
colnames(effect_size_phen_lmer)[30:33]<-c("Anti-parasitic (N=34)","Anti-fungal (N=82)","Vitamins (N=122)","Anti-hypertensive (N=43)")
colnames(effect_size_phen_lmer)[9:11]<-c("Oxygen saturation","Hb total","Cough (N=257)")
colnames(effect_size_phen_lmer)[1]<-"Hb A1c"
colnames(effect_size_phen_lmer)[87:102]<-colnames(effect_size_phen_2)[c(99:105,108:116)]
#colnames(effect_size_phen_lmer)[164:167]<-c("Alcohol daily frequency","Alcohol 6 drinks","Cigarette usage(46)","Cigarette frequency")
colnames(effect_size_phen_lmer)[164:167]<-c("Alcohol daily frequency","Alcohol 6 drinks","Cigarette usage (N=53)","Cigarette frequency")
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
colnames(effect_size_phen)
pval_chk_lmer_all<-rowSums(pval_chk_lmer)
length(which(pval_chk_lmer_all>3))

#colnames(effect_size_phen_lmer)[34:43]
colnames(effect_size_phen_lmer)[44:83]<-c("Reserved(Agree a little)","Reserved(Agree strongly)","Reserved(Disagree a little)","Reserved(Disagree strongly)","Trusting(Agree a little)","Trusting(Agree strongly)","Trusting(Disagree a little)","Trusting(Disagree strongly)","Lazy(Agree a little)","Lazy(Agree strongly)","Lazy(Disagree a little)","Lazy(Disagree strongly)","Relaxed(Agree a little)","Relaxed(Agree strongly)","Relaxed(Disagree a little)","Relaxed(Disagree strongly)","Not creative(Agree a little)","Not creative(Agree strongly)","Not creative(Disagree a little)","Not creative(Disagree strongly)","Outgoing(Agree a little)","Outgoing(Agree strongly)","Outgoing(Disagree a little)","Outgoing(Disagree strongly)","Fault others(Agree a little)","Fault others(Agree strongly)","Fault others(Disagree a little)","Fault others(Disagree strongly)","Thorough job(Agree a little)","Thorough job(Agree strongly)","Thorough job(Disagree a little)","Thorough job(Disagree strongly)","Nervous(Agree a little)","Nervous(Agree strongly)","Nervous(Disagree a little)","Nervous(Disagree strongly)","Openess(Agree a little)","Openess(Agree strongly)","Openess(Disagree a little)","Openess(Disagree strongly)")

fdr_chk_lmer_all<-rowSums(fdr_chk_lmer[,c(1,4:6,9:10,14,13,15,12,11,209,210,16:17,19,22:25,27:34,42:43,86,85,164,166:167,152:154,156:158)])#c(1,4:11,16:17,19:20,22:25)#44:83(personality types discrete)#
length(which(fdr_chk_lmer_all>5))#>5

effect_size_phen_sig<-effect_size_phen_lmer[which(fdr_chk_lmer_all>5),c(1,4:6,9:10,14,13,15,12,11,209,210,16:17,19,22:25,27:34,42:43,86,85,164,166:167,152:154,156:158)]
pval_phen_sig<-pval_chk_lmer[which(fdr_chk_lmer_all>5),c(1,4:6,9:10,14,13,15,12,11,209,210,16:17,19,22:25,27:34,42:43,86,85,164,166:167,152:154,156:158)]
fdr_phen_sig<-fdr_chk_lmer[which(fdr_chk_lmer_all>5),c(1,4:6,9:10,14,13,15,12,11,209,210,16:17,19,22:25,27:34,42:43,86,85,164,166:167,152:154,156:158)]

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
colnames(effect_size_phen_sig_plot)<-as.character(mb_samp_sp_name2[which(fdr_chk_lmer_all>5)])
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

for(i in 1:dim(disp_fdr)[1]){
  if(i==1){
    sp_c<-length(which(disp_fdr[i,]!=""))
  }else{
    sp_c<-rbind(sp_c,length(which(disp_fdr[i,]!="")))
  }
}
colnames(effect_size_phen_sig_plot)[which(sp_c==max(sp_c))]
length(which(grepl("{",colnames(effect_size_phen_sig_plot),fixed=T)==T))

xmo<-match(colnames(effect_size_phen_sig_plot),mb_samp_sp_name_use2)#USE
#xmo2<-match(colnames(effect_size_phen_sig_plot),mb_samp_sp_name_use3)
#xmo3<-match(colnames(effect_size_phen_sig_plot),mb_samp_sp_name_use4)

xt<-(which(is.na(xmo)==F))
#xt2<-(which(is.na(xmo2)==F))

effect_size_phen_sig_plot2<-effect_size_phen_sig_plot[,xt]
disp_fdr2<-disp_fdr[xt,]

#pheatmap(effect_size_phen_sig_plot2[,jkl33[xt2]],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks2,display_numbers = t(disp_fdr[jkl33[xt2],]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
pheatmap(effect_size_phen_sig_plot2,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr2),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

ol<-c(0.001,0.1,10^(1.12-(0.243*5)),10^(1.12-(0.243*4)),10^(1.12-(0.243*3)),10^(1.12-(0.243*2)),10^(1.12-(0.243*1)),10^1.12)
myBreaks11<-c(-ol[length(ol):(-1):1],ol)

pheatmap(effect_size_phen_sig_plot2,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks11,display_numbers = t(disp_fdr2),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)


#effect_size_phen_lmer22<-effect_size_phen_lmer[,c(1,4:6,9:10,14,13,15,12,11,209,210,16:17,19,22:25,27:34,42:43,86,85,164,166:167,152:154,156:158)]
#pval_phen_lmer22<-pval_phen_lmer[,c(1,4:6,9:10,14,13,15,12,11,209,210,16:17,19,22:25,27:34,42:43,86,85,164,166:167,152:154,156:158)]
#fdr_phen_lmer22<-fdr_phen_lmer[,c(1,4:6,9:10,14,13,15,12,11,209,210,16:17,19,22:25,27:34,42:43,86,85,164,166:167,152:154,156:158)]

#Splitting

#pheatmap(effect_size_phen_sig_plot2[,1:103],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr2[1:103,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

#pheatmap(effect_size_phen_sig_plot2[,104:206],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr2[104:206,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
#pheatmap(effect_size_phen_sig_plot2[,207:309],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr2[207:309,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
#pheatmap(effect_size_phen_sig_plot2[,310:412],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr2[310:412,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
#pheatmap(effect_size_phen_sig_plot2[,413:515],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr2[413:515,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
#pheatmap(effect_size_phen_sig_plot2[,516:618],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr2[516:618,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
#pheatmap(effect_size_phen_sig_plot2[,619:722],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr2[619:722,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

#effect_size_phen_lmer2<-effect_size_phen_lmer[,c(1,4:6,9:10,14,13,15,12,11,209,210,16:17,19,22:25,27:34,42:43,86,85,164,166:167,152:154,156:158)]
#rownames(effect_size_phen_lmer2)<-mb_samp_sp_name2
#write.csv(effect_size_phen_lmer2,'health_eff.csv')


pheatmap(effect_size_phen_sig_plot2[c(31:32,36:41),which(disp_fdr2[c(31:32,36:41),]!="")],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot)[c(31:32,36:41)],legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr2[which(disp_fdr2[c(31:32,36:41),]!=""),c(31:32,36:41)]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)


##Phylum analysis
phylum_ind<-which((grepl("c__",mb_samp[,1],fixed=T)==F)&(grepl("p__",mb_samp[,1],fixed=T)==T))
mb_samp_phylum_temp<-as.data.frame(mb_samp[phylum_ind,])
lp<-0
for(i in 1:length(phylum_ind)){
  lp<-rbind(lp,length(which(mb_samp_phylum_temp[i,3:1190]>0)))
}
lp<-lp[2:length(lp)]

mb_samp_phylum<-as.data.frame(mb_samp_phylum_temp[which(lp>200),])

sp_ph_names<-colnames(effect_size_phen_sig_plot2)
ind_sp_ph<-array(0,dim=c(length(sp_ph_names)))

for(i in 1:length(sp_ph_names)){
  for(j in 1:dim(mb_samp)[1]){
    ind_temp1<-unlist(gregexpr("(",as.character(sp_ph_names[i]),fixed=T))
    ind_temp2<-unlist(gregexpr(")",as.character(sp_ph_names[i]),fixed=T))
    if(grepl(substr(sp_ph_names[i],ind_temp1+1,ind_temp2-1),as.character(mb_samp[j,1]),fixed=T)==T){
      ind_sp_ph[i]<-j
    }
  }
}

sp_ph_full_names<-as.character(mb_samp[ind_sp_ph,1])

phyl_ph<-array(0,dim=c(length(sp_ph_full_names)))

for(i in 1:dim(mb_samp_phylum)[1]){
  for(j in 1:length(sp_ph_full_names)){
    if((grepl(as.character(mb_samp_phylum[i,1]),sp_ph_full_names[j],fixed=T))==T){
      #phyl_ph[j]<-as.character(mb_samp_phylum[i,1])
      phyl_ph[j]<-substr(as.character(mb_samp_phylum[i,1]),unlist(gregexpr("p__",as.character(mb_samp_phylum[i,1]),fixed=T))+3,nchar(as.character(mb_samp_phylum[i,1])))
      
      #phyl_ph[j]<-i
    }
  }
}
#effect_size_phen_sig_plot2<-effect_size_phen_sig_plot
phyl_ph<-as.data.frame(phyl_ph)
rownames(phyl_ph)<-colnames(effect_size_phen_sig_plot2)

colnames(phyl_ph)<-"Phylum"

pheatmap(effect_size_phen_sig_plot2,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr2),number_color = "black",treeheight_row = 0,cluster_rows = F,cluster_cols = T,treeheight_col = 0,annotation_col=phyl_ph)

#Corrplot Figure Supp.




#########
#splitting

split_num<-c(718:822)
phyl_ph<-array(0,dim=c(length(sp_ph_full_names[split_num])))

for(i in 1:dim(mb_samp_phylum)[1]){
  for(j in 1:length(sp_ph_full_names[split_num])){
    if((grepl(as.character(mb_samp_phylum[i,1]),sp_ph_full_names[j],fixed=T))==T){
      phyl_ph[j]<-as.character(mb_samp_phylum[i,1])
      #phyl_ph[j]<-i
    }
  }
}

phyl_ph<-as.data.frame(phyl_ph)
rownames(phyl_ph)<-colnames(effect_size_phen_sig_plot[,split_num])

colnames(phyl_ph)<-"Phylum"

pheatmap(effect_size_phen_sig_plot[,split_num],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[split_num,]),number_color = "black",treeheight_row = 0,cluster_rows = F,cluster_cols = T,treeheight_col = 0)

pheatmap(effect_size_phen_sig_plot[,split_num],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[split_num,]),number_color = "black",treeheight_row = 0,cluster_rows = F,cluster_cols = T,treeheight_col = 0,annotation_col=phyl_ph)



##Correlating effect sizes --- maybe wihtout thresholds for species? -- removed cigarette frequency, very good and fair,cognitive impairment, dementia
colnames(effect_size_phen_lmer)[c(12,14)]<-c("Excellent(health)[75]","Poor(health)[152]")
cp1<-effect_size_phen_lmer[which(fdr_chk_lmer_all>0),c(1,4,5,6,9:10,14,12,11,16:17,19,22:25,27:34,42:43,164,166,152:154,156:158)]#,86,85-- cognitive impairment, dementia
colnames(cp1)<-colnames(effect_size_phen_lmer[c(1,4,5,6,9:10,14,12,11,16:17,19,22:25,27:34,42:43,164,166,152:154,156:158)])
#cp_ind<-match(cp1_names,mb_samp_sp_name)
#cp_cor<-cor(cp1[cp_ind,],cp1[cp_ind,],method="pearson")
cp_cor<-cor(cp1,cp1,method="pearson")
pheatmap(cp_cor,treeheight_col = 0)

library(corrplot)

corrplot(as.matrix(cp_cor),method = "square",tl.cex = 0.85,cl.cex = 0.85,tl.col = "black",order= "hclust")

##Co-occurrence of species -- plot

effect_size_phen_sig<-effect_size_phen_lmer[which(fdr_chk_lmer_all>0),c(1,4,5,6,9:10,14,12,11,16:17,19,22:25,27:34,42:43,164,166,152:154,156:158)]
pval_phen_sig<-pval_chk_lmer[which(fdr_chk_lmer_all>0),c(1,4,5,6,9:10,14,12,11,16:17,19,22:25,27:34,42:43,164,166,152:154,156:158)]
fdr_phen_sig<-fdr_chk_lmer[which(fdr_chk_lmer_all>0),c(1,4,5,6,9:10,14,12,11,16:17,19,22:25,27:34,42:43,164,166,152:154,156:158)]

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

cp2<-array(0,dim=c(dim(disp_fdr)[2],dim(disp_fdr)[2]))

for(i in 1:dim(disp_fdr)[2]){
  for(j in 1:dim(disp_fdr)[2]){
   if(i==j){
     cp2[i,j]<-0
   }
    if(i!=j){
      for(k in 1:dim(disp_fdr)[1]){
        if((disp_fdr[k,i]=="+")&(disp_fdr[k,j]=="+")){
          cp2[i,j]<-cp2[i,j]+1
        }else if((disp_fdr[k,i]=="-")&(disp_fdr[k,j]=="-")){
          cp2[i,j]<-cp2[i,j]+1
        }else if((disp_fdr[k,i]=="+")&(disp_fdr[k,j]=="-")){
          cp2[i,j]<-cp2[i,j]-1
        }else if((disp_fdr[k,i]=="-")&(disp_fdr[k,j]=="+")){
          cp2[i,j]<-cp2[i,j]-1
        }
      }
    }
  }
}

cp4<-cp2
for(i in 1:dim(cp2)[1]){
  for(j in 1:dim(cp2)[2]){
    if(cp2[i,j]>0){
      cp4[i,j]<-log(cp2[i,j])/log(dim(disp_fdr)[1])
    }
    if(cp2[i,j]<0){
      cp4[i,j]<-(log(abs(cp2[i,j]))*-1)/log(dim(disp_fdr)[1])
    }
    if(i==j){
      cp4[i,j]<-1
    }
  }
}

colnames(cp4)<-colnames(cp1)
rownames(cp4)<-colnames(cp1)
cp_cor2<-cp4

corrplot(as.matrix(cp_cor2),method = "square",tl.cex = 0.85,cl.cex = 0.85,tl.col = "black",order= "hclust")


#cp3<-cp2
#for(i in 1:dim(cp3)[1]){
#  for(j in 1:dim(cp3)[2]){
#    if(cp2[i,j]>0){
#      cp3[i,j]<-cp2[i,j]/max(cp2)
#    }
#    if(cp2[i,j]<0){
#      cp3[i,j]<-cp2[i,j]/abs(min(cp2))
#    }
#  }
#}

#colnames(cp3)<-colnames(cp1)
#rownames(cp3)<-colnames(cp1)
#cp_cor2<-cp3

#corrplot(as.matrix(cp_cor2),method = "square",tl.cex = 0.85,cl.cex = 0.85,tl.col = "black",order= "hclust")




#breaksList = seq(-max(abs(cp_cor)), +max(abs(cp_cor)),by=2 * max(abs(cp_cor))/55)
#cols <- colorRampPalette(c("Darkblue","White","Darkred"))(length(breaksList))
#cols[seq(length(cols)/2-1,length(cols)/2+1)] <- "White"
#pheatmap(cp_cor,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(cp_cor),legend = T,
#         show_rownames = T, angle_col = 90,
#         fontsize_number = 20,border_color = "#EEEEEE",na_col = "white",
#         color = cols,treeheight_row = 0, treeheight_col = 0,legend_labels = "log(p-value)",
#         fontsize_col = 14,fontsize_row = 12,breaks = breaksList)

#p<-pheatmap(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

#ml_name<-p$tree_col$labels[p$tree_col$order]
#ml_name<-colnames(effect_size_phen_sig_plot)
#View(ml_name)
#ml_nn<-c(1:44,46,48:76)

#pheatmap(effect_size_phen_sig_plot[,ml_nn],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[ml_nn,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

##Known phylum vs unknown

sp_ph_names<-colnames(effect_size_phen_sig_plot)
sp_ph_effect<-colSums(effect_size_phen_sig_plot)

circ_names<-array(0,dim=c(length(sp_ph_names)))

for(i in 1:length(sp_ph_names)[1]){
  circ_names[i]<-ifelse(grepl("{",sp_ph_names[i],fixed=T),1,0)
}
#View(circ_names)
#View(sp_ph_names)
library(ggplot2)
library(dplyr)

slices <- c(length(which(circ_names==1)),length(which(circ_names==0)))
lbls <- c("Unknown species", "Known species")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="Pie Chart of Countries")

mycols2 <- c("#0073C2FF", "#CD534CFF")

count.data <- data.frame(
  class = lbls,
  n = slices,
  prop = pct,
  label2 = c(paste0("Unknown"," (",pct[1],"%)"),paste0("Known"," (",pct[2],"%)"))
)
count.data <- count.data %>%
  arrange(desc(class)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
#count.data

ggplot(count.data, aes(x = "", y = prop, fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = label2), color = "white")+labs(title="Signficicant species",subtitle = "Health phenotypes")+
  scale_fill_manual(values = mycols2) + theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),panel.background = element_blank(),panel.grid.major = element_blank(),axis.text.x=element_blank(), #remove x axis labels
                                              axis.ticks.x=element_blank())

##Effect size contribution

slices <- c(length(which(circ_names==1)),length(which(circ_names==0)))
lbls <- c("Unknown species", "Known species")
#pct<-slices/sum(slices)*100
pct <- round(c(abs(mean(sp_ph_effect[which(circ_names==1)])),abs(mean(sp_ph_effect[which(circ_names==0)])))/sum(abs(mean(sp_ph_effect[which(circ_names==0)])),abs(mean(sp_ph_effect[which(circ_names==1)])))*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="Pie Chart of Countries")

mycols2 <- c("#0073C2FF", "#CD534CFF")

count.data <- data.frame(
  class = lbls,
  n = slices,
  prop = pct,
  label2 = c(paste0("Unknown"," (",pct[1],"%)"),paste0("Known"," (",pct[2],"%)"))
)
count.data <- count.data %>%
  arrange(desc(class)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
#count.data

ggplot(count.data, aes(x = "", y = prop, fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = label2), color = "white")+labs(title="Effect size contribution",subtitle = "Health phenotypes")+
  scale_fill_manual(values = mycols2) + theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),panel.background = element_blank(),panel.grid.major = element_blank(),axis.text.x=element_blank(), #remove x axis labels
                                              axis.ticks.x=element_blank())


#Phylum analysis
#phylum_ind<-which((grepl("c__",mb_samp[,1],fixed=T)==F)&(grepl("p__",mb_samp[,1],fixed=T)==T))

#mb_samp_phylum<-as.data.frame(mb_samp[phylum_ind,1])




#phsm<-cbind(as.character(mb_samp_phylum[,1]),phyl_ph)
#phsm<-cbind(phsm,array("Health factors",dim=c(length(phyl_ph),1)))#array("health",dim=c(length(phyl_ph),1))#array("Socio-economic factors",dim=c(length(phyl_ph),1))#array("Health factors",dim=c(length(phyl_ph),1))
#phsm<-as.data.frame(phsm)
#colnames(phsm)<-c("Phylum","count","Category")

#ggplot(data = phsm, aes(x = as.character(Category), y = as.numeric(as.character(count)), fill = Phylum)) + 
#  geom_bar(stat='identity')+ xlab("")+ylab("Species count") #+ coord_flip()

##Mental health only



fdr_chk_lmer_all<-rowSums(fdr_chk_lmer[,159:160])#c(1,4:11,16:17,19:20,22:25)#44:83(personality types discrete)#
length(which(fdr_chk_lmer_all>0))#c(87:93,151:154,94:102,155:158) -- for all GAD7 and PHQ9

effect_size_phen_sig<-effect_size_phen_lmer[which(fdr_chk_lmer_all>0),159:160]
pval_phen_sig<-pval_chk_lmer[which(fdr_chk_lmer_all>0),159:160]
fdr_phen_sig<-fdr_chk_lmer[which(fdr_chk_lmer_all>0),159:160]

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
rownames(effect_size_phen_sig_plot)<-colnames(effect_size_phen_lmer)[159:160]
colnames(effect_size_phen_sig_plot)<-as.character(mb_samp_sp_name2[which(fdr_chk_lmer_all>0)])
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

#Printing siginifcant dict species

ind_ment_dict<-which(fdr_chk_lmer_all>0)

mb_samp_sp_name_ment_full<-rownames(mb_samp_1188_3)[ind_thresh_use2[,2]]
mb_samp_sp_name_ment_full<-mb_samp_sp_name_ment_full[ind_ment_dict]


mb_samp_1188_ment<-mb_samp_1188_3[ind_thresh_use2[ind_ment_dict,2],]

write.csv(mb_samp_1188_ment,'mental_mb_samp_sp_dict.csv')
write.csv(mb_samp_sp_name_ment_full,'mental_mb_samp_sp_dict_names.csv')


##Correlating effect sizes --mental health plot
colnames(effect_size_phen_lmer)[c(12,14)]<-c("Excellent(health)[75]","Poor(health)[152]")
cp1<-effect_size_phen_lmer[which(fdr_chk_lmer_all>0),c(1,4,5,6,9:10,14,12,11,16:17,19,22:25,27:34,42:43,164,166,152:154,156:158)]#,86,85-- cognitive impairment, dementia
colnames(cp1)<-colnames(effect_size_phen_lmer[c(1,4,5,6,9:10,14,12,11,16:17,19,22:25,27:34,42:43,164,166,152:154,156:158)])
#cp_ind<-match(cp1_names,mb_samp_sp_name)
#cp_cor<-cor(cp1[cp_ind,],cp1[cp_ind,],method="pearson")
cp_cor<-cor(cp1,cp1,method="pearson")
pheatmap(cp_cor,treeheight_col = 0)

library(corrplot)

corrplot(as.matrix(cp_cor),method = "square",tl.cex = 0.85,cl.cex = 0.85,tl.col = "black",order= "hclust")

##########
##Food and animals

effect_size_phen_lmer<-read.csv('effect_size_phen_5_lmer.csv',row.names = 1)
pval_phen_lmer<-read.csv('pval_phen_5_lmer.csv',row.names = 1)
fdr_phen_lmer<-read.csv('fdr_phen_5_lmer.csv',row.names = 1)

effect_size_phen_lmer<-effect_size_phen_lmer[1:2285,]
pval_phen_lmer<-pval_phen_lmer[1:2285,]
fdr_phen_lmer<-fdr_phen_lmer[1:2285,]

effect_size_phen_lmer<-effect_size_phen_lmer[ind_thresh_use,]#Apply it for atleast 100 ppl and above 0.0001
pval_phen_lmer<-pval_phen_lmer[ind_thresh_use,]
fdr_phen_lmer<-fdr_phen_lmer[ind_thresh_use,]

#colnames(effect_size_phen_lmer)[1:26]<-colnames(effect_size_phen_2)[45:70]
colnames(effect_size_phen_lmer)[1:26]<-c("Cat (N=645)","Dog (N=944)","Parakeet (N=102)","Rabbit (N=120)","Horse (N=314)","Mice (N=585)","No pets (N=88)","Cow (N=397)","Goat (N=40)","Pig (N=208)","Chicken (N=1094)","Duck (N=478)","Turkey (N=322)","Sheep (N=33)","Geese (N=68)","No farm animals (N=69)","Bat (N=410)","Lizard (N=423)","Monkey (N=12)","Snake (N=437)","Bird (N=794)","Possum (N=467)","Rat (N=457)","Squirrel (N=480)","Other Wild (N=1)","No wild animals (N=243)")
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
length(which(fdr_chk_lmer_all>6))#>6 -- 74 species

effect_size_phen_sig<-effect_size_phen_lmer[which(fdr_chk_lmer_all>6),c(1:3,7:15,5,4,16,6,17:18,20:24,26:42,46)]
pval_phen_sig<-pval_chk_lmer[which(fdr_chk_lmer_all>6),c(1:3,7:15,5,4,16,6,17:18,20:24,26:42,46)]
fdr_phen_sig<-fdr_chk_lmer[which(fdr_chk_lmer_all>6),c(1:3,7:15,5,4,16,6,17:18,20:24,26:42,46)]

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
colnames(effect_size_phen_sig_plot)<-as.character(mb_samp_sp_name2[which(fdr_chk_lmer_all>6)])
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


pheatmap(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#f5f5f5",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

#library(corrplot)
#corrplot(effect_size_phen_sig_plot,mdethod="number",is.corr=F)

#effect_size_phen_lmer23<-cbind(effect_size_phen_lmer22,effect_size_phen_lmer[ind_thresh_use2[,2],c(1:3,7:15,5,4,16,6,17:18,20:24,26:42,46)])#run upto 948
#pval_phen_lmer23<-cbind(pval_phen_lmer22,pval_phen_lmer[ind_thresh_use2[,2],c(1:3,7:15,5,4,16,6,17:18,20:24,26:42,46)])
#fdr_phen_lmer23<-cbind(fdr_phen_lmer22,fdr_phen_lmer[ind_thresh_use2[,2],c(1:3,7:15,5,4,16,6,17:18,20:24,26:42,46)])


for(i in 1:dim(disp_fdr)[1]){
  if(i==1){
    sp_c<-length(which(disp_fdr[i,]!=""))
  }else{
    sp_c<-rbind(sp_c,length(which(disp_fdr[i,]!="")))
  }
}
colnames(effect_size_phen_sig_plot)[which(sp_c==max(sp_c))]

xmo<-match(colnames(effect_size_phen_sig_plot),mb_samp_sp_name_use2)#USE
#xmo2<-match(colnames(effect_size_phen_sig_plot),mb_samp_sp_name_use3)
#xmo3<-match(colnames(effect_size_phen_sig_plot),mb_samp_sp_name_use4)

xt<-(which(is.na(xmo)==F))
#xt2<-(which(is.na(xmo2)==F))

effect_size_phen_sig_plot2<-effect_size_phen_sig_plot[,xt]
disp_fdr2<-disp_fdr[xt,]

#pheatmap(effect_size_phen_sig_plot2[,jkl33[xt2]],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks2,display_numbers = t(disp_fdr[jkl33[xt2],]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
pheatmap(effect_size_phen_sig_plot[,xt],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[xt,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

ol<-c(0.001,0.1,10^(1.12-(0.243*5)),10^(1.12-(0.243*4)),10^(1.12-(0.243*3)),10^(1.12-(0.243*2)),10^(1.12-(0.243*1)),10^1.12)
myBreaks11<-c(-ol[length(ol):(-1):1],ol)

pheatmap(effect_size_phen_sig_plot[,xt],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks11,display_numbers = t(disp_fdr[xt,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

out2<-pheatmap(effect_size_phen_sig_plot[,xt],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[xt,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
#Splitting
#pheatmap(effect_size_phen_sig_plot[,xt[1:105]],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[xt[1:105],]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
#pheatmap(effect_size_phen_sig_plot[,xt[106:210]],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[xt[106:210],]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
#pheatmap(effect_size_phen_sig_plot[,xt[211:315]],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[xt[211:315],]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
#pheatmap(effect_size_phen_sig_plot[,xt[316:420]],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[xt[316:420],]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
#pheatmap(effect_size_phen_sig_plot[,xt[421:525]],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[xt[421:525],]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
#pheatmap(effect_size_phen_sig_plot[,xt[526:631]],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[xt[526:631],]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

#effect_size_phen_lmer2<-effect_size_phen_lmer[,c(1:3,7:15,5,4,16,6,17:18,20:24,26:42)]
#rownames(effect_size_phen_lmer2)<-mb_samp_sp_name2
#write.csv(effect_size_phen_lmer2,'food_an_eff.csv')#Need to change to ind_thesh_use2[,2] in lines 945,946,947

#pheatmap(effect_size_phen_sig_plot2[,jkl33[xt]],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[jkl33[xt],]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

##Stats
ymo<-match(colnames(effect_size_phen_sig_plot2),mb_samp_sp_name_use2)
effect_size_phen_sig_plot3<-effect_size_phen_sig_plot2[,which(is.na(ymo)==F)]
disp_fdr3<-disp_fdr[which(is.na(ymo)==F),]
for(i in 1:dim(disp_fdr3)[1]){
  if(i==1){
    sp_c<-length(which(disp_fdr3[i,]!=""))
  }else{
    sp_c<-rbind(sp_c,length(which(disp_fdr3[i,]!="")))
  }
}
colnames(effect_size_phen_sig_plot3)[which(sp_c==max(sp_c))]
length(which(disp_fdr3=="-"))
length(which(grepl("{",colnames(effect_size_phen_sig_plot3),fixed=T)==T))


##Phylum analysis
phylum_ind<-which((grepl("c__",mb_samp[,1],fixed=T)==F)&(grepl("p__",mb_samp[,1],fixed=T)==T))
mb_samp_phylum_temp<-as.data.frame(mb_samp[phylum_ind,])
lp<-0
for(i in 1:length(phylum_ind)){
  lp<-rbind(lp,length(which(mb_samp_phylum_temp[i,3:1190]>0)))
}
lp<-lp[2:length(lp)]

mb_samp_phylum<-as.data.frame(mb_samp_phylum_temp[which(lp>200),])

sp_ph_names<-colnames(effect_size_phen_sig_plot2)
ind_sp_ph<-array(0,dim=c(length(sp_ph_names)))

for(i in 1:length(sp_ph_names)){
  for(j in 1:dim(mb_samp)[1]){
    ind_temp1<-unlist(gregexpr("(",as.character(sp_ph_names[i]),fixed=T))
    ind_temp2<-unlist(gregexpr(")",as.character(sp_ph_names[i]),fixed=T))
    if(grepl(substr(sp_ph_names[i],ind_temp1+1,ind_temp2-1),as.character(mb_samp[j,1]),fixed=T)==T){
      ind_sp_ph[i]<-j
    }
  }
}

sp_ph_full_names<-as.character(mb_samp[ind_sp_ph,1])

phyl_ph<-array(0,dim=c(length(sp_ph_full_names)))

for(i in 1:dim(mb_samp_phylum)[1]){
  for(j in 1:length(sp_ph_full_names)){
    if((grepl(as.character(mb_samp_phylum[i,1]),sp_ph_full_names[j],fixed=T))==T){
      #phyl_ph[j]<-as.character(mb_samp_phylum[i,1])
      phyl_ph[j]<-substr(as.character(mb_samp_phylum[i,1]),unlist(gregexpr("p__",as.character(mb_samp_phylum[i,1]),fixed=T))+3,nchar(as.character(mb_samp_phylum[i,1])))
      #phyl_ph[j]<-i
    }
  }
}
#effect_size_phen_sig_plot2<-effect_size_phen_sig_plot
phyl_ph<-as.data.frame(phyl_ph)
rownames(phyl_ph)<-colnames(effect_size_phen_sig_plot2)

colnames(phyl_ph)<-"Phylum"

pheatmap(effect_size_phen_sig_plot2,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[xt,]),number_color = "black",treeheight_row = 0,cluster_rows = F,cluster_cols = T,treeheight_col = 0,annotation_col=phyl_ph)

######
#splitting

split_num<-c(613:714)
phyl_ph<-array(0,dim=c(length(sp_ph_full_names[split_num])))

for(i in 1:dim(mb_samp_phylum)[1]){
  for(j in 1:length(sp_ph_full_names[split_num])){
    if((grepl(as.character(mb_samp_phylum[i,1]),sp_ph_full_names[j],fixed=T))==T){
      phyl_ph[j]<-as.character(mb_samp_phylum[i,1])
      #phyl_ph[j]<-i
    }
  }
}

phyl_ph<-as.data.frame(phyl_ph)
rownames(phyl_ph)<-colnames(effect_size_phen_sig_plot[,split_num])

colnames(phyl_ph)<-"Phylum"

pheatmap(effect_size_phen_sig_plot[,split_num],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[split_num,]),number_color = "black",treeheight_row = 0,cluster_rows = F,cluster_cols = T,treeheight_col = 0)

pheatmap(effect_size_phen_sig_plot[,split_num],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[split_num,]),number_color = "black",treeheight_row = 0,cluster_rows = F,cluster_cols = T,treeheight_col = 0,annotation_col=phyl_ph)



##Correlating effect sizes --- maybe wihtout thresholds for species? -- removed cigarette frequency, very good and fair,cognitive impairment, dementia

cp1<-effect_size_phen_lmer[which(fdr_chk_lmer_all>0),c(1:3,7:15,5,4,16,6,17:24,26:42)]#,86,85-- cognitive impairment, dementia
colnames(cp1)<-colnames(effect_size_phen_lmer[c(1:3,7:15,5,4,16,6,17:24,26:42)])
#cp_ind<-match(cp1_names,mb_samp_sp_name)
#cp_cor<-cor(cp1[cp_ind,],cp1[cp_ind,],method="pearson")
cp_cor<-cor(cp1,cp1,method="pearson")
pheatmap(cp_cor,treeheight_col = 0)

library(corrplot)

corrplot(as.matrix(cp_cor),method = "square",tl.cex = 0.85,cl.cex = 0.85,tl.col = "black",order= "hclust")

##Co-occurrence of species -- plot

effect_size_phen_sig<-effect_size_phen_lmer[which(fdr_chk_lmer_all>0),c(1:3,7:15,5,4,16,6,17:24,26:42)]
pval_phen_sig<-pval_chk_lmer[which(fdr_chk_lmer_all>0),c(1:3,7:15,5,4,16,6,17:24,26:42)]
fdr_phen_sig<-fdr_chk_lmer[which(fdr_chk_lmer_all>0),c(1:3,7:15,5,4,16,6,17:24,26:42)]

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

cp2<-array(0,dim=c(dim(disp_fdr)[2],dim(disp_fdr)[2]))

for(i in 1:dim(disp_fdr)[2]){
  for(j in 1:dim(disp_fdr)[2]){
    if(i==j){
      cp2[i,j]<-0
    }
    if(i!=j){
      for(k in 1:dim(disp_fdr)[1]){
        if((disp_fdr[k,i]=="+")&(disp_fdr[k,j]=="+")){
          cp2[i,j]<-cp2[i,j]+1
        }else if((disp_fdr[k,i]=="-")&(disp_fdr[k,j]=="-")){
          cp2[i,j]<-cp2[i,j]+1
        }else if((disp_fdr[k,i]=="+")&(disp_fdr[k,j]=="-")){
          cp2[i,j]<-cp2[i,j]-1
        }else if((disp_fdr[k,i]=="-")&(disp_fdr[k,j]=="+")){
          cp2[i,j]<-cp2[i,j]-1
        }
      }
    }
  }
}

cp4<-cp2
for(i in 1:dim(cp2)[1]){
  for(j in 1:dim(cp2)[2]){
    if(cp2[i,j]>0){
      cp4[i,j]<-log(cp2[i,j])/log(dim(disp_fdr)[1])
    }
    if(cp2[i,j]<0){
      cp4[i,j]<-(log(abs(cp2[i,j]))*-1)/log(dim(disp_fdr)[1])
    }
    if(i==j){
      cp4[i,j]<-1
    }
  }
}

colnames(cp4)<-colnames(cp1)
rownames(cp4)<-colnames(cp1)
cp_cor2<-cp4

corrplot(as.matrix(cp_cor2),method = "square",tl.cex = 0.85,cl.cex = 0.85,tl.col = "black",order= "hclust")




##Known phylum vs unknown

sp_ph_names<-colnames(effect_size_phen_sig_plot)
sp_ph_effect<-colSums(effect_size_phen_sig_plot)

circ_names<-array(0,dim=c(length(sp_ph_names)))

for(i in 1:length(sp_ph_names)[1]){
  circ_names[i]<-ifelse(grepl("{",sp_ph_names[i],fixed=T),1,0)
}
#View(circ_names)
#View(sp_ph_names)
library(ggplot2)
library(dplyr)

slices <- c(length(which(circ_names==1)),length(which(circ_names==0)))
lbls <- c("Unknown species", "Known species")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="Pie Chart of Countries")

mycols2 <- c("#0073C2FF", "#CD534CFF")

count.data <- data.frame(
  class = lbls,
  n = slices,
  prop = pct,
  label2 = c(paste0("Unknown"," (",pct[1],"%)"),paste0("Known"," (",pct[2],"%)"))
)
count.data <- count.data %>%
  arrange(desc(class)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
#count.data

ggplot(count.data, aes(x = "", y = prop, fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = label2), color = "white")+labs(title="Signficicant species",subtitle = "Animal and food phenotypes")+
  scale_fill_manual(values = mycols2) + theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),panel.background = element_blank(),panel.grid.major = element_blank(),axis.text.x=element_blank(), #remove x axis labels
                                              axis.ticks.x=element_blank())

##Effect size contribution

slices <- c(length(which(circ_names==1)),length(which(circ_names==0)))
lbls <- c("Unknown species", "Known species")
#pct<-slices/sum(slices)*100
pct <- round(c(mean(abs(sp_ph_effect[which(circ_names==1)])),mean(abs(sp_ph_effect[which(circ_names==0)])))/(mean(abs(sp_ph_effect[which(circ_names==1)]))+mean(abs(sp_ph_effect[which(circ_names==0)])))*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="Pie Chart of Countries")

mycols2 <- c("#0073C2FF", "#CD534CFF")

count.data <- data.frame(
  class = lbls,
  n = slices,
  prop = pct,
  label2 = c(paste0("Unknown"," (",pct[1],"%)"),paste0("Known"," (",pct[2],"%)"))
)
count.data <- count.data %>%
  arrange(desc(class)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
#count.data

ggplot(count.data, aes(x = "", y = prop, fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = label2), color = "white")+labs(title="Effect size contribution",subtitle = "Animal and food phenotypes")+
  scale_fill_manual(values = mycols2) + theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),panel.background = element_blank(),panel.grid.major = element_blank(),axis.text.x=element_blank(), #remove x axis labels
                                              axis.ticks.x=element_blank())

##Phylum analysis
ind_sp_ph<-array(0,dim=c(length(sp_ph_names)))

for(i in 1:length(sp_ph_names)){
  for(j in 1:dim(mb_samp)[1]){
    ind_temp1<-unlist(gregexpr("(",as.character(sp_ph_names[i]),fixed=T))
    ind_temp2<-unlist(gregexpr(")",as.character(sp_ph_names[i]),fixed=T))
    if(grepl(substr(sp_ph_names[i],ind_temp1+1,ind_temp2-1),as.character(mb_samp[j,1]),fixed=T)==T){
      ind_sp_ph[i]<-j
    }
  }
}

sp_ph_full_names<-as.character(mb_samp[ind_sp_ph,1])

phyl_ph<-array(0,dim=c(dim(mb_samp_phylum)[1]))

for(i in 1:length(phyl_ph)){
  for(j in 1:length(sp_ph_full_names)){
    if((grepl(as.character(mb_samp_phylum[i,1]),sp_ph_full_names[j],fixed=T))==T){
      
      
      phyl_ph[i]<-phyl_ph[i]+1
    }
  }
}




phsm<-cbind(as.character(mb_samp_phylum[,1]),phyl_ph)
phsm<-cbind(phsm,array("Food & Animal factors",dim=c(length(phyl_ph),1)))#array("health",dim=c(length(phyl_ph),1))#array("Socio-economic factors",dim=c(length(phyl_ph),1))#array("Health factors",dim=c(length(phyl_ph),1))
phsm<-as.data.frame(phsm)
colnames(phsm)<-c("Phylum","count","Category")

ggplot(data = phsm, aes(x = as.character(Category), y = as.numeric(as.character(count)), fill = Phylum)) + 
  geom_bar(stat='identity')+ xlab("")+ylab("Species count") #+ coord_flip()



##Socioeconomic factors -- rerun altruism??




######################################################################################################################
#Complex heatmap

library("ComplexHeatmap")

effect_size_phen_lmer<-read.csv('effect_size_phen_5_lmer.csv',row.names = 1)
pval_phen_lmer<-read.csv('pval_phen_5_lmer.csv',row.names = 1)
fdr_phen_lmer<-read.csv('fdr_phen_5_lmer.csv',row.names = 1)

effect_size_phen_lmer<-effect_size_phen_lmer[1:2285,]
pval_phen_lmer<-pval_phen_lmer[1:2285,]
fdr_phen_lmer<-fdr_phen_lmer[1:2285,]

effect_size_phen_lmer<-effect_size_phen_lmer[ind_thresh_use,]#Apply it for atleast 100 ppl and above 0.0001
pval_phen_lmer<-pval_phen_lmer[ind_thresh_use,]
fdr_phen_lmer<-fdr_phen_lmer[ind_thresh_use,]

colnames(effect_size_phen_lmer)[1:26]<-colnames(effect_size_phen_2)[45:70]
colnames(effect_size_phen_lmer)[38:41]<-c("Natural juice","Chicken","Beef/pork","Ham/sausages/hotdog")
colnames(effect_size_phen_lmer)[33]<-"Cream/butter"

pval_chk_lmer<-array(0,dim=c(dim(pval_phen_lmer_rv)[1],dim(pval_phen_lmer)[2]))
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

fdr_chk_lmer_all<-rowSums(fdr_chk_lmer[,c(1:3,7:15,5,4,16,6,17:24,26:42)])#c(1,4:11,16:17,19:20,22:25)#44:83(personality types discrete)#
length(which(fdr_chk_lmer_all>5))

effect_size_phen_sig<-effect_size_phen_lmer[which(fdr_chk_lmer_all>5),c(1:3,7:15,5,4,16,6,17:24,26:42)]
pval_phen_sig<-pval_chk_lmer[which(fdr_chk_lmer_all>5),c(1:3,7:15,5,4,16,6,17:24,26:42)]
fdr_phen_sig<-fdr_chk_lmer[which(fdr_chk_lmer_all>5),c(1:3,7:15,5,4,16,6,17:24,26:42)]

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
rownames(effect_size_phen_sig_plot)<-colnames(effect_size_phen_lmer)[c(1:3,7:15,5,4,16,6,17:24,26:42)]
colnames(effect_size_phen_sig_plot)<-as.character(mb_samp_sp_name2[which(fdr_chk_lmer_all>5)])
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


pheatmap(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

############################################################################
#######
##Socio-economic factors

effect_size_phen_lmer<-read.csv('effect_size_phen_6_lmer.csv',row.names = 1)
pval_phen_lmer<-read.csv('pval_phen_6_lmer.csv',row.names = 1)
fdr_phen_lmer<-read.csv('fdr_phen_6_lmer.csv',row.names = 1)

effect_size_phen_lmer<-cbind(effect_size_phen_lmer,read.csv('effect_size_phen_7_lmer.csv',row.names = 1))
pval_phen_lmer<-cbind(pval_phen_lmer,read.csv('pval_phen_7_lmer.csv',row.names = 1))
fdr_phen_lmer<-cbind(fdr_phen_lmer,read.csv('fdr_phen_7_lmer.csv',row.names = 1))

effect_size_phen_lmer<-cbind(effect_size_phen_lmer,read.csv('effect_size_phen_9_lmer.csv',row.names = 1))
pval_phen_lmer<-cbind(pval_phen_lmer,read.csv('pval_phen_9_lmer.csv',row.names = 1))
fdr_phen_lmer<-cbind(fdr_phen_lmer,read.csv('fdr_phen_9_lmer.csv',row.names = 1))

effect_size_phen_lmer<-effect_size_phen_lmer[1:2285,]
pval_phen_lmer<-pval_phen_lmer[1:2285,]
fdr_phen_lmer<-fdr_phen_lmer[1:2285,]

#effect_size_phen_lmer24<-cbind(effect_size_phen_lmer23,effect_size_phen_lmer[ind_thresh_use2[,2],c(36:39,102:104,22:25,27:31,33:34,32,10:11,40:41,98:99,101,3:4,96:97,47,48,51,54,58,60,67,71,73,81,87:88,91,92,94)])#run upto 948
#pval_phen_lmer24<-cbind(pval_phen_lmer23,pval_phen_lmer[ind_thresh_use2[,2],c(36:39,102:104,22:25,27:31,33:34,32,10:11,40:41,98:99,101,3:4,96:97,47,48,51,54,58,60,67,71,73,81,87:88,91,92,94)])
#fdr_phen_lmer24<-cbind(fdr_phen_lmer23,fdr_phen_lmer[ind_thresh_use2[,2],c(36:39,102:104,22:25,27:31,33:34,32,10:11,40:41,98:99,101,3:4,96:97,47,48,51,54,58,60,67,71,73,81,87:88,91,92,94)])
#rownames(effect_size_phen_lmer24)<-mb_samp_sp_name_use2
#temp<-1:2285
mb_samp_sp_name_full2<-rownames(data_transformed)[ind_111111[ind_thresh_use2[,2],1]]
#write.csv(mb_samp_sp_name_full2,'mb_name_full_788_2.csv')
#write.csv(mb_samp_sp_name_full,'mb_name_full_788.csv')
#write.csv(effect_size_phen_lmer24,'esz_ph_f.csv')
#write.csv(pval_phen_lmer24,'pval_ph_f.csv')
#write.csv(fdr_phen_lmer24,'fdr_ph_f.csv')

effect_size_phen_lmer<-effect_size_phen_lmer[ind_thresh_use,]#Apply it for atleast 100 ppl and above 0.0001
pval_phen_lmer<-pval_phen_lmer[ind_thresh_use,]
fdr_phen_lmer<-fdr_phen_lmer[ind_thresh_use,]


colnames(effect_size_phen_lmer)[1:2]<-c("GAD7(binary)","PHQ9(binary)")
colnames(effect_size_phen_lmer)[3:9]<-c("Travel","Monthly expenditure","Alcohol frequency","Alcohol daily","Alcohol(6 drinks)","Cigarette use","Cigarette frequency")
colnames(effect_size_phen_lmer)[11]<-"Risk taking"
colnames(effect_size_phen_lmer)[22:29]<-c("Friend ties (same building)","Friend ties (different building)","Betweeness (friendship)","Transitivity (friendship)","Risky interaction","Family ties (same building)","Family ties (different building)","Transitivity (familial)")
#colnames(effect_size_phen_lmer)[27:30]<-colnames(effect_size_phen_2)[25:28]
colnames(effect_size_phen_lmer)[33:34]<-c("Clustering coefficient(all ties)","Betweeness(all ties)")
colnames(effect_size_phen_lmer)[36:39]<-c("Partner live duration","Number of partners","Partner live age","Living with partner (N=375)")
#colnames(effect_size_phen_lmer)[42:95]<-colnames(effect_size_phen_2)[c(138:144,134:137,145:187)]
colnames(effect_size_phen_lmer)[42:95]<-c("Electricity (N=1069)","Radio (N=519)","TV (N=565)","Cellphone (N=927)","No cellphone (N=12)","Refrigerator (N=397)","No electronics (N=38)","Chimney (N=959)","No chimney (N=151)","Stove (N=61)","No stove (N=1)","Wood (N=1109)","Gas(fuel) [N=46]","Electricity fuel (N=14)","Kerosene (N=2)","None fuel (N=1)","Separate kitchen (N=1106)","Cement floor (N=552)","Earth/Sand floor (N=507)","Ceramic floor (N=88)","Tiles floor (N=19)","Mud bricks floor (N=1)","Wood floor (N=0)","Other floor (N=5)","Wooden windows (N=1035)","Glass windows (N=51)","Metal windows (N=6)","Unfinished windows (N=35)","No windows (N=45)","Clay/mud walls(N=693)","Clay brick walls (N=3)","Cement walls (N=452)","Cane/palm/trunks walls (N=0)","Wood unpolished walls (N=21)","Wood polished walls (N=0)","Discarded materials walls (N=1)","No walls (N=0)","Other walls (N=2)","Plastic roof (N=7)","Metal roof (N=982)","Clay roof (N=120)","Thatch/palm roof (N=10)","Concrete roof (N=27)","Wood roof (N=22)","Other roof(N=4)","Sleeping rooms","Spring(protected) [N=917]","Spring(unprotected) [N=19]","Tube well (N=72)","Dug well(protected) [N=78]","Dug well(unprotected) [N=26]","Surface water (N=10)","Bottled water (N=15)","Other water (N=6)")
colnames(effect_size_phen_lmer)[13:21]<-c("Education(1st grade)","Education(2nd grade)","Education(3rd grade)","Education(4th grade)","Education(5th grade)","Education(6th grade)","Education(More than secondary)","Education(Secondary)","Education(Some secondary)")
colnames(effect_size_phen_lmer)[30:32]<-c("Betweeness (familial)","Degree (all ties)","%kin (to 3rd degree)")
colnames(effect_size_phen_lmer)[35]<-"k-cycle centrality"
colnames(effect_size_phen_lmer)[40:41]<-c("Washing hands (N=1130)","Distance to village center")
colnames(effect_size_phen_lmer)[102:104]<-c("Grades 1-3 (N=409)","Grades 4-6 (N=334)","Grades >6 (N=167)")#colnames(effect_size_phen_lmer)[102:104]<-c("Primary (N=373)","Middle (N=311)","Secondary (N=113)")
colnames(effect_size_phen_lmer)[96:101]<-c("Household size","Household wealth index","Distance to main road","Distance to health center","Number of churches","Deforestation (%)")
colnames(effect_size_phen_lmer)[32]<-"Kin percentage (to third degree)"

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


fdr_chk_lmer_all<-rowSums(fdr_chk_lmer[,c(36:39,102:104,22:25,27:31,33:34,32,10:11,40:41,98:99,101,3:4,96:97,47,48,51,54,58,60,67,71,73,81,87:88,91,92)])
length(which(fdr_chk_lmer_all>0))

effect_size_phen_sig<-effect_size_phen_lmer[which(fdr_chk_lmer_all>0),c(36:39,102:104,22:25,27:31,33:34,32,10:11,40:41,98:99,101,3:4,96:97,47,48,51,54,58,60,67,71,73,81,87:88,91,92)]
pval_phen_sig<-pval_chk_lmer[which(fdr_chk_lmer_all>0),c(36:39,102:104,22:25,27:31,33:34,32,10:11,40:41,98:99,101,3:4,96:97,47,48,51,54,58,60,67,71,73,81,87:88,91,92)]
fdr_phen_sig<-fdr_chk_lmer[which(fdr_chk_lmer_all>0),c(36:39,102:104,22:25,27:31,33:34,32,10:11,40:41,98:99,101,3:4,96:97,47,48,51,54,58,60,67,71,73,81,87:88,91,92)]

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
rownames(effect_size_phen_sig_plot)<-colnames(effect_size_phen_lmer)[c(36:39,102:104,22:25,27:31,33:34,32,10:11,40:41,98:99,101,3:4,96:97,47,48,51,54,58,60,67,71,73,81,87:88,91,92)]
colnames(effect_size_phen_sig_plot)<-as.character(mb_samp_sp_name2[which(fdr_chk_lmer_all>0)])
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

##

xmo<-match(colnames(effect_size_phen_sig_plot),mb_samp_sp_name_use2)
xt<-(which(is.na(xmo)==F))
effect_size_phen_sig_plot2<-effect_size_phen_sig_plot[,xt]
disp_fdr2<-disp_fdr[xt,]

#effect_size_phen_lmer2<-effect_size_phen_lmer[,c(36:39,102:104,22:25,27:31,33:34,32,10:11,40:41,98:99,101,3:4,96:97,47,48,51,54,58,60,67,71,73,81,87:88,91,92,94)]
#rownames(effect_size_phen_lmer2)<-mb_samp_sp_name2
#write.csv(effect_size_phen_lmer2,'socio_eco_eff.csv')#Need to switch to ind_thresh_use2[,2] in lines 1487,1488,1489 to use this
#Splitting

#for(i in 1:dim(effect_size_phen_sig_plot)[1]){
#  for(j in 1:length(jkl23)){
#    if(disp_fdr[jkl23[j],i]=="-"){
#      if(effect_size_phen_sig_plot2[i,jkl23[j]]>(0.796*-1)){
#        effect_size_phen_sig_plot2[i,jkl23[j]]<-(0.8*-1)
#      }
#    }
#  }
#}

#myColor3 <- colorRampPalette(c("#cb181d","white", "#33a02c"))(21)#myColor <- colorRampPalette(c("Orange", "white", "Darkblue"))(paletteLength)
#myBreaks3 <- c(seq(min(effect_size_phen_sig_plot2),0, length.out=ceiling(21/2) + 1), 
#              seq(max(effect_size_phen_sig_plot2)/21, max(effect_size_phen_sig_plot2), length.out=floor(21/2)))


#pheatmap(effect_size_phen_sig_plot2[,1:102],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor3,breaks=myBreaks3,display_numbers = t(disp_fdr2[1:102,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
#pheatmap(effect_size_phen_sig_plot2[,103:204],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor3,breaks=myBreaks3,display_numbers = t(disp_fdr2[103:204,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
#pheatmap(effect_size_phen_sig_plot2[,205:307],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor3,breaks=myBreaks3,display_numbers = t(disp_fdr2[205:307,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
#pheatmap(effect_size_phen_sig_plot2[,308:410],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor3,breaks=myBreaks3,display_numbers = t(disp_fdr2[308:410,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
#pheatmap(effect_size_phen_sig_plot2[,411:513],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor3,breaks=myBreaks3,display_numbers = t(disp_fdr2[411:513,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
#pheatmap(effect_size_phen_sig_plot2[,514:615],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor3,breaks=myBreaks3,display_numbers = t(disp_fdr2[514:615,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
#pheatmap(effect_size_phen_sig_plot2[,616:718],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor3,breaks=myBreaks3,display_numbers = t(disp_fdr2[616:718,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)






#socio_eco_sp_com<-colnames(effect_size_phen_sig_plot)
#write.csv(colnames(effect_size_phen_sig_plot),'socio_eco_sp_com.csv')
#socio_eco_sp_com<-read.csv('socio_eco_sp_com.csv',row.names = 1)
##Picking species
#ind1<-which(disp_fdr[,26]!="")
#ind2<-which(disp_fdr[,27]!="")
#ind12<-unique(c(ind1,ind2))

#ind_c<-array(0,dim=c(dim(disp_fdr)[2],dim(disp_fdr)[1]))
#for(i in 1:dim(ind_c)[1]){
#  temp<-which(disp_fdr[,i]!="")
#  ind_c[i,1:length(temp)]<-temp
#}
#write.csv(ind_c,'ind_c.csv')
#ind_c<-read.csv('ind_c.csv',row.names=1)

##Choosing species from ind_c

#ind_c1<-c(81,241,389,612,765,246,428,498,729,660,159,164,1,21,62,101,)
#ph_ind<-c(1,2,4,8,11,12,13,14,23:27,39,42)
#jkl<-unique(sort(ind_c[ph_ind,]))
#for(i in 1:dim(ind_c)[1]){
#  if(i%in%ph_ind){
#  if(i==1){
#    jkl2<-ind_c[i,which(ind_c[i,]>0)]
#  }else{
#    jkl2<-c(jkl2,ind_c[i,which(ind_c[i,]>0)])
#  }
#  }
#}

#pheatmap(effect_size_phen_sig_plot[,jkl2[1:102]],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[jkl2[1:102],]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = F)
#write.csv(colnames(effect_size_phen_sig_plot[,jkl2[1:102]]),'t11.csv')

#pheatmap(effect_size_phen_sig_plot[,jkl2[103:204]],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[jkl2[103:204],]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = F)
#write.csv(colnames(effect_size_phen_sig_plot[,jkl2[103:204]]),'t12.csv')

#pheatmap(effect_size_phen_sig_plot[,jkl2[205:306]],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[jkl2[205:306],]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = F)
#write.csv(colnames(effect_size_phen_sig_plot[,jkl2[205:306]]),'t13.csv')

#jkl2_ind<-c(1,4,7,17,20,23,27,42,46,48,52,57,66,70,72,79,85,102,111,128,132,157,160,178,194,197,204,206,224,236,246,247,262,269,286,279,299,300,302,306)
#okl<-colnames(effect_size_phen_sig_plot[,jkl2[jkl2_ind]])
#write.csv(okl,'ph_okl.csv')

okl<-read.csv('ph_okl.csv',row.names=1)
#Run from here

#rownames(effect_size_phen_sig_plot)
jkl21<-match(as.character(okl[,1]),colnames(effect_size_phen_sig_plot))
jkl22<-match(socio_eco_sp_com[,1],colnames(effect_size_phen_sig_plot))
jkl23<-c(jkl21,jkl22)
jkl23<-unique(jkl23)

myColor <- colorRampPalette(c("#cb181d","white", "#33a02c"))(paletteLength)#myColor <- colorRampPalette(c("Orange", "white", "Darkblue"))(paletteLength)
myBreaks <- c(seq(min(effect_size_phen_sig_plot),0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(effect_size_phen_sig_plot)/paletteLength, max(effect_size_phen_sig_plot), length.out=floor(paletteLength/2)))


pheatmap(effect_size_phen_sig_plot[,jkl23],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[jkl23,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

paletteLength<-13
myColor <- colorRampPalette(c("#cb181d","white", "#33a02c"))(paletteLength)#myColor <- colorRampPalette(c("Orange", "white", "Darkblue"))(paletteLength)
myBreaks2 <- c(seq(min(effect_size_phen_sig_plot[,jkl23]),0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(effect_size_phen_sig_plot[,jkl23])/paletteLength, max(effect_size_phen_sig_plot[,jkl23]), length.out=floor(paletteLength/2)))

#Color adjust -- for 
length(which(disp_fdr[jkl23,]=="-"))
effect_size_phen_sig_plot2<-effect_size_phen_sig_plot
for(i in 1:dim(effect_size_phen_sig_plot)[1]){
  for(j in 1:length(jkl23)){
    if(disp_fdr[jkl23[j],i]=="-"){
      if(effect_size_phen_sig_plot2[i,jkl23[j]]>(0.796*-1)){
        effect_size_phen_sig_plot2[i,jkl23[j]]<-(0.8*-1)
      }
    }
  }
}

pheatmap(effect_size_phen_sig_plot2[,jkl23],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks2,display_numbers = t(disp_fdr[jkl23,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

######
#Adding deforestation species, and removing some species

jkl24<-which(disp_fdr[,26]!="") #Adding deforestation#was 28 before
jkl25<-which(disp_fdr[,23]!="") #Adding distance from village center#was 28 before
jkl26<-which(disp_fdr[,21]!="") #Adding distance from village center#was 28 before

nm<-as.data.frame(colnames(effect_size_phen_sig_plot))
#write.csv(nm,'nm.csv')
jkl23_nc<-c(60,79,3,41,47,26,65,10,21,39,80,23,17)

jkl33<-c(jkl24,jkl23[which((c(1:81)%in%jkl23_nc)==F)],jkl25[c(1,2,4,5)],jkl26[1:3])#added jkl25 for dist from village center

paletteLength<-13
myColor <- colorRampPalette(c("#cb181d","white", "#33a02c"))(paletteLength)#myColor <- colorRampPalette(c("Orange", "white", "Darkblue"))(paletteLength)
myBreaks2 <- c(seq(min(effect_size_phen_sig_plot[,jkl33]),0, length.out=ceiling(paletteLength/2) + 1), 
               seq(max(effect_size_phen_sig_plot[,jkl33])/paletteLength, max(effect_size_phen_sig_plot[,jkl33]), length.out=floor(paletteLength/2)))

#Color adjust -- for 
length(which(disp_fdr[jkl33,]=="-"))
effect_size_phen_sig_plot2<-effect_size_phen_sig_plot
for(i in 1:dim(effect_size_phen_sig_plot)[1]){
  for(j in 1:length(jkl33)){
    if(disp_fdr[jkl33[j],i]=="-"){
      if(effect_size_phen_sig_plot2[i,jkl33[j]]>(0.9*-1)){
        effect_size_phen_sig_plot2[i,jkl33[j]]<-(0.9*-1)
      }
    }
  }
}

pheatmap(effect_size_phen_sig_plot2[,jkl33],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks2,display_numbers = t(disp_fdr[jkl33,]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

xmo<-match(colnames(effect_size_phen_sig_plot2[,jkl33]),mb_samp_sp_name_use2)#Use

xt<-(which(is.na(xmo)==F))

#pheatmap(effect_size_phen_sig_plot2[,jkl33[xt2]],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks2,display_numbers = t(disp_fdr[jkl33[xt2],]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
pheatmap(effect_size_phen_sig_plot2[,jkl33[xt]],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks2,display_numbers = t(disp_fdr[jkl33[xt],]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
out2<-pheatmap(effect_size_phen_sig_plot2[,jkl33[xt]],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks2,display_numbers = t(disp_fdr[jkl33[xt],]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

paletteLength<-15
myColor <- colorRampPalette(c("#cb181d","white", "#33a02c"))(paletteLength)#myColor <- colorRampPalette(c("Orange", "white", "Darkblue"))(paletteLength)


out2<-pheatmap(effect_size_phen_sig_plot2[,jkl33[xt]],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks11,display_numbers = t(disp_fdr[jkl33[xt],]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)


##Stats
ymo<-match(colnames(effect_size_phen_sig_plot2),mb_samp_sp_name_use2)
effect_size_phen_sig_plot3<-effect_size_phen_sig_plot2[,which(is.na(ymo)==F)]
disp_fdr3<-disp_fdr[which(is.na(ymo)==F),]
for(i in 1:dim(disp_fdr3)[1]){
  if(i==1){
    sp_c<-length(which(disp_fdr3[i,]!=""))
  }else{
    sp_c<-rbind(sp_c,length(which(disp_fdr3[i,]!="")))
  }
}
colnames(effect_size_phen_sig_plot3)[which(sp_c==max(sp_c))]
length(which(disp_fdr3=="-"))
length(which(grepl("{",colnames(effect_size_phen_sig_plot3),fixed=T)==T))

######

##Phylum analysis

phylum_ind<-which((grepl("c__",mb_samp[,1],fixed=T)==F)&(grepl("p__",mb_samp[,1],fixed=T)==T))
mb_samp_phylum_temp<-as.data.frame(mb_samp[phylum_ind,])
lp<-0
for(i in 1:length(phylum_ind)){
  lp<-rbind(lp,length(which(mb_samp_phylum_temp[i,3:1190]>0)))
}
lp<-lp[2:length(lp)]

mb_samp_phylum<-as.data.frame(mb_samp_phylum_temp[which(lp>200),])

sp_ph_names<-colnames(effect_size_phen_sig_plot)
ind_sp_ph<-array(0,dim=c(length(sp_ph_names)))

for(i in 1:length(sp_ph_names)){
  for(j in 1:dim(mb_samp)[1]){
    ind_temp1<-unlist(gregexpr("(",as.character(sp_ph_names[i]),fixed=T))
    ind_temp2<-unlist(gregexpr(")",as.character(sp_ph_names[i]),fixed=T))
    if(grepl(substr(sp_ph_names[i],ind_temp1+1,ind_temp2-1),as.character(mb_samp[j,1]),fixed=T)==T){
      ind_sp_ph[i]<-j
    }
  }
}

sp_ph_full_names<-as.character(mb_samp[ind_sp_ph,1])

phyl_ph<-array(0,dim=c(length(sp_ph_full_names[jkl33])))
phyl_ph<-array(0,dim=c(length(sp_ph_full_names[jkl33[xt]])))
#phyl_ph<-array(0,dim=c(length(sp_ph_full_names[which(is.na(ymo)==F)])))
sp_ph_full_names2<-sp_ph_full_names[jkl33[xt]]#[which(is.na(ymo)==F)]#[jkl33[xt]]

for(i in 1:dim(mb_samp_phylum)[1]){
  for(j in 1:length(sp_ph_full_names2)){
    if((grepl(as.character(mb_samp_phylum[i,1]),sp_ph_full_names2[j],fixed=T))==T){
      #phyl_ph[j]<-as.character(mb_samp_phylum[i,1])
      phyl_ph[j]<-substr(as.character(mb_samp_phylum[i,1]),unlist(gregexpr("p__",as.character(mb_samp_phylum[i,1]),fixed=T))+3,nchar(as.character(mb_samp_phylum[i,1])))
      #phyl_ph[j]<-i
    }
  }
}
#effect_size_phen_sig_plot2<-effect_size_phen_sig_plot
phyl_ph<-as.data.frame(phyl_ph)
rownames(phyl_ph)<-colnames(effect_size_phen_sig_plot2[,jkl33[xt]])

colnames(phyl_ph)<-"Phylum"

pheatmap(effect_size_phen_sig_plot2[,jkl33[xt]],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[jkl33[xt],]),number_color = "black",treeheight_row = 0,cluster_rows = F,cluster_cols = T,treeheight_col = 0,annotation_col=phyl_ph)




#jkl3<-jkl#do common ones separately

#pheatmap(effect_size_phen_sig_plot[,jkl3[1:102]],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[jkl3[1:102],]),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

#write.csv(colnames(effect_size_phen_sig_plot[,jkl3[1:102]]),'t1.csv')

ind_c2<-array(0,dim=c(dim(disp_fdr)[2],dim(disp_fdr)[1]))
for(i in 1:dim(ind_c2)[1]){
  temp<-which(disp_fdr[jkl3,i]!="")
  ind_c2[i,1:length(temp)]<-temp
}
write.csv(ind_c2,'ind_c2.csv')

ind_cc<-c(76,)




######
#splitting

split_num<-c(701:790)
phyl_ph<-array(0,dim=c(length(sp_ph_full_names[split_num])))

for(i in 1:dim(mb_samp_phylum)[1]){
  for(j in 1:length(sp_ph_full_names[split_num])){
    if((grepl(as.character(mb_samp_phylum[i,1]),sp_ph_full_names[j],fixed=T))==T){
      phyl_ph[j]<-as.character(mb_samp_phylum[i,1])
      #phyl_ph[j]<-i
    }
  }
}

phyl_ph<-as.data.frame(phyl_ph)
rownames(phyl_ph)<-colnames(effect_size_phen_sig_plot[,split_num])

colnames(phyl_ph)<-"Phylum"

pheatmap(effect_size_phen_sig_plot[,split_num],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[split_num,]),number_color = "black",treeheight_row = 0,cluster_rows = F,cluster_cols = T,treeheight_col = 0)

pheatmap(effect_size_phen_sig_plot[,split_num],annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr[split_num,]),number_color = "black",treeheight_row = 0,cluster_rows = F,cluster_cols = T,treeheight_col = 0,annotation_col=phyl_ph)



##Phylum analysis
phylum_ind<-which((grepl("c__",mb_samp[,1],fixed=T)==F)&(grepl("p__",mb_samp[,1],fixed=T)==T))
mb_samp_phylum_temp<-as.data.frame(mb_samp[phylum_ind,])
lp<-0
for(i in 1:length(phylum_ind)){
  lp<-rbind(lp,length(which(mb_samp_phylum_temp[i,3:1190]>0)))
}
lp<-lp[2:length(lp)]

mb_samp_phylum<-as.data.frame(mb_samp_phylum_temp[which(lp>200),])

sp_ph_names<-colnames(effect_size_phen_sig_plot)
ind_sp_ph<-array(0,dim=c(length(sp_ph_names)))

for(i in 1:length(sp_ph_names)){
  for(j in 1:dim(mb_samp)[1]){
    ind_temp1<-unlist(gregexpr("(",as.character(sp_ph_names[i]),fixed=T))
    ind_temp2<-unlist(gregexpr(")",as.character(sp_ph_names[i]),fixed=T))
    if(grepl(substr(sp_ph_names[i],ind_temp1+1,ind_temp2-1),as.character(mb_samp[j,1]),fixed=T)==T){
      ind_sp_ph[i]<-j
    }
  }
}

sp_ph_full_names<-as.character(mb_samp[ind_sp_ph,1])

phyl_ph<-array(0,dim=c(length(sp_ph_full_names)))

for(i in 1:dim(mb_samp_phylum)[1]){
  for(j in 1:length(sp_ph_full_names)){
    if((grepl(as.character(mb_samp_phylum[i,1]),sp_ph_full_names[j],fixed=T))==T){
      #phyl_ph[j]<-as.character(mb_samp_phylum[i,1])
      phyl_ph[j]<-substr(as.character(mb_samp_phylum[i,1]),unlist(gregexpr("p__",as.character(mb_samp_phylum[i,1]),fixed=T))+3,nchar(as.character(mb_samp_phylum[i,1])))
      #phyl_ph[j]<-i
    }
  }
}
#effect_size_phen_sig_plot2<-effect_size_phen_sig_plot
phyl_ph<-as.data.frame(phyl_ph)
rownames(phyl_ph)<-colnames(effect_size_phen_sig_plot)

colnames(phyl_ph)<-"Phylum"

pheatmap(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr),number_color = "black",treeheight_row = 0,cluster_rows = F,cluster_cols = T,treeheight_col = 0,annotation_col=phyl_ph)



##Correlation matrix
cp1<-effect_size_phen_lmer[which(fdr_chk_lmer_all>0),c(36:39,12,22:25,27:30,31,33:34,32,35,3:4,10:11,40,43,44,46,48,50,53,57,59,66:67,72,80:81,86:87,90:93)]#,86,85-- cognitive impairment, dementia
colnames(cp1)<-colnames(effect_size_phen_lmer[c(36:39,12,22:25,27:30,31,33:34,32,35,3:4,10:11,40,43,44,46,48,50,53,57,59,66:67,72,80:81,86:87,90:93)])
#cp_ind<-match(cp1_names,mb_samp_sp_name)
#cp_cor<-cor(cp1[cp_ind,],cp1[cp_ind,],method="pearson")
cp_cor<-cor(cp1,cp1,method="pearson")
pheatmap(cp_cor,treeheight_col = 0)

library(corrplot)

corrplot(as.matrix(cp_cor),method = "square",tl.cex = 0.75,cl.cex = 0.85,tl.col = "black",order= "hclust")#was tl.cex=0.8,cl.cex=0.85; 0.7,0.2



##Co-occurrence of species -- plot

effect_size_phen_sig<-effect_size_phen_lmer[which(fdr_chk_lmer_all>0),c(36:39,12,22:25,27:30,31,33:34,32,35,3:4,10:11,40,43,44,46,48,50,53,57,59,66:67,72,80:81,86:87,90:93)]
pval_phen_sig<-pval_chk_lmer[which(fdr_chk_lmer_all>0),c(36:39,12,22:25,27:30,31,33:34,32,35,3:4,10:11,40,43,44,46,48,50,53,57,59,66:67,72,80:81,86:87,90:93)]
fdr_phen_sig<-fdr_chk_lmer[which(fdr_chk_lmer_all>0),c(36:39,12,22:25,27:30,31,33:34,32,35,3:4,10:11,40,43,44,46,48,50,53,57,59,66:67,72,80:81,86:87,90:93)]

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

cp2<-array(0,dim=c(dim(disp_fdr)[2],dim(disp_fdr)[2]))

for(i in 1:dim(disp_fdr)[2]){
  for(j in 1:dim(disp_fdr)[2]){
    if(i==j){
      cp2[i,j]<-0
    }
    if(i!=j){
      for(k in 1:dim(disp_fdr)[1]){
        if((disp_fdr[k,i]=="+")&(disp_fdr[k,j]=="+")){
          cp2[i,j]<-cp2[i,j]+1
        }else if((disp_fdr[k,i]=="-")&(disp_fdr[k,j]=="-")){
          cp2[i,j]<-cp2[i,j]+1
        }else if((disp_fdr[k,i]=="+")&(disp_fdr[k,j]=="-")){
          cp2[i,j]<-cp2[i,j]-1
        }else if((disp_fdr[k,i]=="-")&(disp_fdr[k,j]=="+")){
          cp2[i,j]<-cp2[i,j]-1
        }
      }
    }
  }
}

cp4<-cp2
for(i in 1:dim(cp2)[1]){
  for(j in 1:dim(cp2)[2]){
    if(cp2[i,j]>0){
      cp4[i,j]<-log(cp2[i,j])/log(dim(disp_fdr)[1])
    }
    if(cp2[i,j]<0){
      cp4[i,j]<-(log(abs(cp2[i,j]))*-1)/log(dim(disp_fdr)[1])
    }
    if(i==j){
      cp4[i,j]<-1
    }
  }
}

colnames(cp4)<-colnames(cp1)
rownames(cp4)<-colnames(cp1)
cp_cor2<-cp4

corrplot(as.matrix(cp_cor2),method = "square",tl.cex = 0.75,cl.cex = 0.85,tl.col = "black",order= "hclust")

##Known phylum vs unknown

sp_ph_names<-colnames(effect_size_phen_sig_plot)
sp_ph_effect<-colSums(effect_size_phen_sig_plot)

circ_names<-array(0,dim=c(length(sp_ph_names)))

for(i in 1:length(sp_ph_names)[1]){
  circ_names[i]<-ifelse(grepl("{",sp_ph_names[i],fixed=T),1,0)
}
#View(circ_names)
#View(sp_ph_names)
library(ggplot2)
library(dplyr)

slices <- c(length(which(circ_names==1)),length(which(circ_names==0)))
lbls <- c("Unknown species", "Known species")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="Pie Chart of Countries")

mycols2 <- c("#0073C2FF", "#CD534CFF")

count.data <- data.frame(
  class = lbls,
  n = slices,
  prop = pct,
  label2 = c(paste0("Unknown"," (",pct[1],"%)"),paste0("Known"," (",pct[2],"%)"))
)
count.data <- count.data %>%
  arrange(desc(class)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
#count.data

ggplot(count.data, aes(x = "", y = prop, fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = label2), color = "white")+labs(title="Signficicant species",subtitle = "Socio-economic phenotypes")+
  scale_fill_manual(values = mycols2) + theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),panel.background = element_blank(),panel.grid.major = element_blank(),axis.text.x=element_blank(), #remove x axis labels
                                              axis.ticks.x=element_blank())

##Effect size contribution

slices <- c(length(which(circ_names==1)),length(which(circ_names==0)))
lbls <- c("Unknown species", "Known species")
#pct<-slices/sum(slices)*100
pct <- round(c(abs(mean(sp_ph_effect[which(circ_names==1)])),abs(mean(sp_ph_effect[which(circ_names==0)])))/sum(abs(mean(sp_ph_effect[which(circ_names==0)])),abs(mean(sp_ph_effect[which(circ_names==1)])))*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="Pie Chart of Countries")

mycols2 <- c("#0073C2FF", "#CD534CFF")

count.data <- data.frame(
  class = lbls,
  n = slices,
  prop = pct,
  label2 = c(paste0("Unknown"," (",pct[1],"%)"),paste0("Known"," (",pct[2],"%)"))
)
count.data <- count.data %>%
  arrange(desc(class)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
#count.data

ggplot(count.data, aes(x = "", y = prop, fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = label2), color = "white")+labs(title="Effect size contribution",subtitle = "Socio-economic phenotypes")+
  scale_fill_manual(values = mycols2) + theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),panel.background = element_blank(),panel.grid.major = element_blank(),axis.text.x=element_blank(), #remove x axis labels
                                              axis.ticks.x=element_blank())

##Phylum analysis
ind_sp_ph<-array(0,dim=c(length(sp_ph_names)))

for(i in 1:length(sp_ph_names)){
  for(j in 1:dim(mb_samp)[1]){
    ind_temp1<-unlist(gregexpr("(",as.character(sp_ph_names[i]),fixed=T))
    ind_temp2<-unlist(gregexpr(")",as.character(sp_ph_names[i]),fixed=T))
    if(grepl(substr(sp_ph_names[i],ind_temp1+1,ind_temp2-1),as.character(mb_samp[j,1]),fixed=T)==T){
      ind_sp_ph[i]<-j
    }
  }
}

sp_ph_full_names<-as.character(mb_samp[ind_sp_ph,1])

phyl_ph<-array(0,dim=c(dim(mb_samp_phylum)[1]))

for(i in 1:length(phyl_ph)){
  for(j in 1:length(sp_ph_full_names)){
    if((grepl(as.character(mb_samp_phylum[i,1]),sp_ph_full_names[j],fixed=T))==T){
      
      
      phyl_ph[i]<-phyl_ph[i]+1
    }
  }
}




phsm<-cbind(as.character(mb_samp_phylum[,1]),phyl_ph)
phsm<-cbind(phsm,array("Socio-economic factors",dim=c(length(phyl_ph),1)))#array("health",dim=c(length(phyl_ph),1))#array("Socio-economic factors",dim=c(length(phyl_ph),1))#array("Health factors",dim=c(length(phyl_ph),1))
phsm<-as.data.frame(phsm)
colnames(phsm)<-c("Phylum","count","Category")

ggplot(data = phsm, aes(x = as.character(Category), y = as.numeric(as.character(count)), fill = Phylum)) + 
  geom_bar(stat='identity')+ xlab("")+ylab("Species count") #+ coord_flip()





#####Dangerous health

effect_size_phen_lmer<-read.csv('effect_size_phen_8_lmer.csv',row.names = 1)
pval_phen_lmer<-read.csv('pval_phen_8_lmer.csv',row.names = 1)
fdr_phen_lmer<-read.csv('fdr_phen_8_lmer.csv',row.names = 1)

effect_size_phen_lmer<-effect_size_phen_lmer[1:2285,]
pval_phen_lmer<-pval_phen_lmer[1:2285,]
fdr_phen_lmer<-fdr_phen_lmer[1:2285,]

effect_size_phen_lmer<-effect_size_phen_lmer[ind_thresh_use2[,2],]#Apply it for atleast 100 ppl and above 0.0001
pval_phen_lmer<-pval_phen_lmer[ind_thresh_use2[,2],]
fdr_phen_lmer<-fdr_phen_lmer[ind_thresh_use2[,2],]

colnames(effect_size_phen_2)

#colnames(effect_size_phen_lmer)<-c("Hb A1c(>7)","Hb A1c(6.5-7)","Hb A1c(5.7-6.5)","Diastolic(>89)","Systolic(>119),Diastolic(>79)","Systolic(>119),Diastolic(<=79)","Systolic(>129)","BMI(<18)","BMI(>35)","BMI(25-30)","BMI(30-35)","Heart rate(>100)","Oxygen saturation(<97)","Hb total(elevated)","Hb total(low)")

colnames(effect_size_phen_lmer)<-c("Hb A1c(>7)","Hb A1c(6.5-7)","Hb A1c(5.7-6.5)","Diastolic(\u2265 90)","Systolic(\u2265 120),Diastolic(\u2265 80)","Systolic(\u2265 120),Diastolic(\u2264 80)","Systolic(\u2265 130)","BMI(<18)","BMI(>35)","BMI(25-30)","BMI(30-35)","Heart rate(>100)","Oxygen saturation(<97)","Hb total(elevated)","Hb total(low)")


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
colnames(effect_size_phen_lmer)
pval_chk_lmer_all<-rowSums(pval_chk_lmer)
length(which(pval_chk_lmer_all>3))


fdr_chk_lmer_all<-rowSums(fdr_chk_lmer[,c(3,2,1,6,5,7,4,8,10,11,9,12:13,15)])
length(which(fdr_chk_lmer_all>1))

effect_size_phen_sig<-effect_size_phen_lmer[which(fdr_chk_lmer_all>1),c(3,2,1,6,5,7,4,8,10,11,9,12:13,15)]
pval_phen_sig<-pval_chk_lmer[which(fdr_chk_lmer_all>1),c(3,2,1,6,5,7,4,8,10,11,9,12:13,15)]
fdr_phen_sig<-fdr_chk_lmer[which(fdr_chk_lmer_all>1),c(3,2,1,6,5,7,4,8,10,11,9,12:13,15)]

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
rownames(effect_size_phen_sig_plot)<-colnames(effect_size_phen_lmer)[c(3,2,1,6,5,7,4,8,10,11,9,12:13,15)]
colnames(effect_size_phen_sig_plot)<-as.character(mb_samp_sp_name2[which(fdr_chk_lmer_all>1)])
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

sp_ph_names<-colnames(effect_size_phen_sig_plot)
ind_sp_ph<-array(0,dim=c(length(sp_ph_names)))

for(i in 1:length(sp_ph_names)){
  for(j in 1:dim(mb_samp)[1]){
    ind_temp1<-unlist(gregexpr("(",as.character(sp_ph_names[i]),fixed=T))
    ind_temp2<-unlist(gregexpr(")",as.character(sp_ph_names[i]),fixed=T))
    if(grepl(substr(sp_ph_names[i],ind_temp1+1,ind_temp2-1),as.character(mb_samp[j,1]),fixed=T)==T){
      ind_sp_ph[i]<-j
    }
  }
}

sp_ph_full_names<-as.character(mb_samp[ind_sp_ph,1])

phyl_ph<-array(0,dim=c(length(sp_ph_full_names)))
#phyl_ph<-array(0,dim=c(length(sp_ph_full_names[which(is.na(ymo)==F)])))
sp_ph_full_names2<-sp_ph_full_names#[which(is.na(ymo)==F)]#[jkl33[xt]]

for(i in 1:dim(mb_samp_phylum)[1]){
  for(j in 1:length(sp_ph_full_names2)){
    if((grepl(as.character(mb_samp_phylum[i,1]),sp_ph_full_names2[j],fixed=T))==T){
      #phyl_ph[j]<-as.character(mb_samp_phylum[i,1])
      phyl_ph[j]<-substr(as.character(mb_samp_phylum[i,1]),unlist(gregexpr("p__",as.character(mb_samp_phylum[i,1]),fixed=T))+3,nchar(as.character(mb_samp_phylum[i,1])))
      #phyl_ph[j]<-i
    }
  }
}
#effect_size_phen_sig_plot2<-effect_size_phen_sig_plot
phyl_ph<-as.data.frame(phyl_ph)
rownames(phyl_ph)<-colnames(effect_size_phen_sig_plot)

colnames(phyl_ph)<-"Phylum"

pheatmap(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr),number_color = "black",treeheight_row = 0,cluster_rows = F,cluster_cols = T,treeheight_col = 0,annotation_col=phyl_ph)



length(which(disp_fdr!=""))

for(i in 1:dim(disp_fdr)[1]){
  if(i==1){
    sp_c<-length(which(disp_fdr[i,]!=""))
  }else{
    sp_c<-rbind(sp_c,length(which(disp_fdr[i,]!="")))
  }
}

colnames(effect_size_phen_sig_plot)[which(sp_c==max(sp_c))]




##################################
##Filler plots LMER model

##Alpha diversity
library(vegan)
mb_samp5<-read.csv('mb_samp_1188_3.csv',row.names = 1)
div_alpha<-array(0,dim=c(dim(mb_samp5)[2]))
#temp<-array(0,dim=c(dim(mb_samp_sp4)[2]))
#rm(diversity)

for(i in 1:dim(mb_samp5)[2]){
  div_alpha[i]<-diversity(mb_samp5[,i])#Default is Shannon!-- try thresholding? was [,i]-- miniscule difference
  #temp[i]<-diversity(mb_samp_sp4[,i],index="simpson")
}

#dds_div<-cbind(div_alpha,dds)
#colnames(dds_div)<-c("Diversity","DDS")
#write.csv(dds_div,'dds_diversity.csv')
#View(div_alpha)
#View(temp)

#phen_chronic<-array(NA,dim=c(dim(phen_all_use)[1],6)) #-- pet vs non pet didnt work
#phen_chronic[,2]<-phen_all_use[,161]
#phen_chronic[,4]<-phen_all_use[,170]
#phen_chronic[,6]<-phen_all_use[,180]

#for(i in 1:dim(phen_chronic)[1]){
#  if(is.na(phen_chronic[i,2])==T){
#    phen_chronic[i,1]<-"Pets"
#  }
#  if(is.na(phen_chronic[i,4])==T){
#    phen_chronic[i,3]<-"Farm"
#  }
 # if(is.na(phen_chronic[i,6])==T){
#    phen_chronic[i,5]<-"Wild"
#  }
#  
#}

phen_chronic<-phen_all_use[,c(155,156,157,161,165,166,170,160,175,180)]



#phen_chronic<-phen_all_use[,c(161,170,180)]-- lets try to do both? -- some individual animals

for(i in 1:dim(phen_chronic)[2]){
  if(i==1){
    ind_ch<-which(is.na(phen_chronic[,i])==F)
  }else{
    ind_ch<-c(ind_ch,which(is.na(phen_chronic[,i])==F))  
  }
  
}

ind_ch2<-unique(ind_ch)
ind_all<-c(1:1186)
ind_healthy<-ind_all[!(ind_all %in% ind_ch2)]

#Alpha diversity

div_alpha_ph<-array(NA,c(length(div_alpha),dim(phen_chronic)[2]))

i<-3
for(i in 1:dim(div_alpha_ph)[2]){
  #if(i==1){
  #  div_alpha_ph[ind_healthy,i]<-div_alpha[ind_healthy]
  #}else{
    ind_tp<-which(is.na(phen_chronic[,i])==F)
    div_alpha_ph[ind_tp,i]<-div_alpha[ind_tp]
  #}
}

div_alpha_ph<-as.data.frame(div_alpha_ph)


div_alpha_ph2<-rbind(cbind(div_alpha_ph[,1],array("Cats",dim=c(dim(div_alpha_ph)[1]))),cbind(div_alpha_ph[,2],array("Dogs",dim=c(dim(div_alpha_ph)[1]))))
div_alpha_ph2<-rbind(div_alpha_ph2,cbind(div_alpha_ph[,3],array("Parakeets",dim=c(dim(div_alpha_ph)[1]))))
div_alpha_ph2<-rbind(div_alpha_ph2,cbind(div_alpha_ph[,4],array("None pet",dim=c(dim(div_alpha_ph)[1]))))
div_alpha_ph2<-rbind(div_alpha_ph2,cbind(div_alpha_ph[,5],array("Chicken",dim=c(dim(div_alpha_ph)[1]))))
div_alpha_ph2<-rbind(div_alpha_ph2,cbind(div_alpha_ph[,6],array("Duck",dim=c(dim(div_alpha_ph)[1]))))
div_alpha_ph2<-rbind(div_alpha_ph2,cbind(div_alpha_ph[,7],array("None farm",dim=c(dim(div_alpha_ph)[1]))))
div_alpha_ph2<-rbind(div_alpha_ph2,cbind(div_alpha_ph[,8],array("Mice",dim=c(dim(div_alpha_ph)[1]))))
div_alpha_ph2<-rbind(div_alpha_ph2,cbind(div_alpha_ph[,9],array("Bird",dim=c(dim(div_alpha_ph)[1]))))
div_alpha_ph2<-rbind(div_alpha_ph2,cbind(div_alpha_ph[,10],array("None wild",dim=c(dim(div_alpha_ph)[1]))))


colnames(div_alpha_ph2)<-c("alpha_div","chronic")
div_alpha_ph2<-as.data.frame(div_alpha_ph2)

div_alpha_ph2$chronic <- factor(div_alpha_ph2$chronic, levels=c("Cats","Dogs","Parakeets","None pet","Chicken","Duck","None farm","Mice","Bird","None wild"),ordered=TRUE)
library(ggplot2)
ggplot(div_alpha_ph2, aes(x=chronic, y=as.numeric(as.character(alpha_div)),group=chronic)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a"))+labs(x="",y="Alpha diversity")+
  theme(text = element_text(size=15))+ scale_x_discrete(breaks=seq(1,10,1))#+ylim(0,0.35)#,"#ff7f00","#cab2d6","#6a3d9a","#ffff99"#x="Relationship",

fp1<-div_alpha_ph[which(is.na(div_alpha_ph[,7])==F),7]
fp2<-div_alpha_ph[which(is.na(div_alpha_ph[,6])==F),6]
fp3<-div_alpha_ph[which(is.na(div_alpha_ph[,5])==F),5]

t.test(fp1,fp2)##Not at all significant
t.test(fp1,fp3)

##Beta diversity

for(i in 1:3){
  if(i==1){
    len_ind_tp<-which(is.na(phen_chronic[,i])==F)
  }else{
    len_ind_tp<-c(len_ind_tp,which(is.na(phen_chronic[,i])==F))
  }
}

len_ind_tp<-unique(len_ind_tp)



ind_healthy<-which(is.na(div_alpha_ph[,4])==F)
bray_chronic<-array(NA,dim=c(length(ind_healthy)*(7000),2))

tpc<-0
ind_tp<-which(is.na(phen_chronic[,1])==F)
for(i in 1:length(ind_healthy)){
  for(j in 1:length(ind_tp)){
    bray_chronic[(i-1)*length(ind_tp)+j,]<-cbind(vegdist(t(cbind(mb_samp_sp4[,ind_healthy[i]],mb_samp_sp4[,ind_tp[j]]))),"Cats")
    tpc<-tpc+1  
  }
}

tpc_all<-tpc
ind_tp<-which(is.na(phen_chronic[,2])==F)
for(i in 1:length(ind_healthy)){
  for(j in 1:length(ind_tp)){
    bray_chronic[(i-1)*length(ind_tp)+j+tpc_all,]<-cbind(vegdist(t(cbind(mb_samp_sp4[,ind_healthy[i]],mb_samp_sp4[,ind_tp[j]]))),"Dogs")
    tpc<-tpc+1 
  }
}
tpc_all<-tpc
ind_tp<-which(is.na(phen_chronic[,3])==F)
for(i in 1:length(ind_healthy)){
  for(j in 1:length(ind_tp)){
    bray_chronic[(i-1)*length(ind_tp)+j+tpc_all,]<-cbind(vegdist(t(cbind(mb_samp_sp4[,ind_healthy[i]],mb_samp_sp4[,ind_tp[j]]))),"Parakeets")
    tpc<-tpc+1 
  }
}


ind_healthy<-which(is.na(div_alpha_ph[,7])==F)

tpc_all<-tpc
ind_tp<-which(is.na(phen_chronic[,5])==F)
for(i in 1:length(ind_healthy)){
  for(j in 1:length(ind_tp)){
    bray_chronic[(i-1)*length(ind_tp)+j+tpc_all,]<-cbind(vegdist(t(cbind(mb_samp_sp4[,ind_healthy[i]],mb_samp_sp4[,ind_tp[j]]))),"Chicken")
    tpc<-tpc+1 
  }
}
tpc_all<-tpc
ind_tp<-which(is.na(phen_chronic[,6])==F)
for(i in 1:length(ind_healthy)){
  for(j in 1:length(ind_tp)){
    bray_chronic[(i-1)*length(ind_tp)+j+tpc_all,]<-cbind(vegdist(t(cbind(mb_samp_sp4[,ind_healthy[i]],mb_samp_sp4[,ind_tp[j]]))),"Ducks")
    tpc<-tpc+1 
  }
}


ind_healthy<-which(is.na(div_alpha_ph[,10])==F)
tpc_all<-tpc
ind_tp<-which(is.na(phen_chronic[,8])==F)
for(i in 1:length(ind_healthy)){
  for(j in 1:length(ind_tp)){
    bray_chronic[(i-1)*length(ind_tp)+j+tpc_all,]<-cbind(vegdist(t(cbind(mb_samp_sp4[,ind_healthy[i]],mb_samp_sp4[,ind_tp[j]]))),"Mice")
    tpc<-tpc+1 
  }
}
tpc_all<-tpc
ind_tp<-which(is.na(phen_chronic[,9])==F)
for(i in 1:length(ind_healthy)){
  for(j in 1:length(ind_tp)){
    bray_chronic[(i-1)*length(ind_tp)+j+tpc_all,]<-cbind(vegdist(t(cbind(mb_samp_sp4[,ind_healthy[i]],mb_samp_sp4[,ind_tp[j]]))),"Birds")
    tpc<-tpc+1 
  }
}


colnames(bray_chronic)<-c("bray","chronic")

bray_chronic<-bray_chronic[which(is.na(bray_chronic[,2])==F),]
bray_chronic<-as.data.frame(bray_chronic)

bray_chronic$chronic <- factor(bray_chronic$chronic, levels=c("Cats","Dogs","Parakeets","Chicken","Ducks","Mice","Birds"),ordered=TRUE)

ggplot(bray_chronic, aes(x=chronic, y=as.numeric(as.character(bray)),group=chronic)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#1f78b4","#b2df8a","#fb9a99","#e31a1c","#ff7f00","#cab2d6"))+labs(x="",y="Bray-curtis(vs no animals)")+
  theme(text = element_text(size=15))+ scale_x_discrete(breaks=seq(1,7,1))#+ylim(0,0.35)#,"#ff7f00","#cab2d6","#6a3d9a","#ffff99"#x="Relationship",


###############################
#Balloon plot of FFQ
library(ggpubr)

#colnames(phen_all_use)[75:93]

phen_food<-array(0,dim=c(dim(phen_all_use[1:4,75:90])))

for(i in 1:dim(phen_all_use)[1]){
  for(j in 75:90){
    if(is.na(phen_all_use[i,j])==F){
      if(phen_all_use[i,j]=="Never/rarely"){#Change to 1/50,1/10,4/7,1
        phen_food[1,(j-74)]<-phen_food[1,(j-74)]+1
      }else if(phen_all_use[i,j]=="A few days per month"){
        phen_food[2,(j-74)]<-phen_food[2,(j-74)]+1
      }else if(phen_all_use[i,j]=="A few days per week"){
        phen_food[3,(j-74)]<-phen_food[3,(j-74)]+1
      }else if(phen_all_use[i,j]=="Every day"){
        phen_food[4,(j-74)]<-phen_food[4,(j-74)]+1
      }
    }
  }
}

phen_food<-as.data.frame(phen_food)
colnames(phen_food)<-c("Beans","Tortillas","Rice","Bread","Milk","Yogurt","Cream/butter","Cheese","Eggs","Vegetables","Fruits","Natural juice","Chicken","Beef/Pork","Ham/sausages/hotdog","Fish")#,"Soda","Fruit juice","Chips"
rownames(phen_food)<-c("Never/rarely","A few days per month","A few days per week","Every day")

ggballoonplot(phen_food[4:-1:1,])
ggballoonplot(phen_food[4:-1:1,],color = "#000000", fill = "#0073C2FF",ggtheme = theme(panel.grid=element_blank(),panel.background = element_blank(),axis.text=element_text(size=15,color="black")))

#
#library(gplots)

#heatmap.2(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

#split<-array(0,dim=c(dim(effect_size_phen_sig_plot)[1]))
#split[1:4]<-"#1f78b4"
#split[5:15]<-"#e31a1c"
#split[16:25]<-"#ff7f00"
#split[26:41]<-"#000000"


#heatmap.2(effect_size_phen_sig_plot,hclustfun = hclust,col=myColor,breaks=myBreaks,trace="none",dendrogram="none",Rowv=FALSE,RowSideColors=split,density.info = "none")



#length(rownames(effect_size_phen_sig_plot))


#Complex heatmap -- FAIL

#split = rep(1:5, each = 10)

#colnames(effect_size_phen_lmer)
#split<-array(0,dim=c(dim(effect_size_phen_lmer)[2]))

#split[c(1:3,7)]<-1
#split[c(8:15,5,4,16)]<-2
#split[c(6,17:24,26)]<-3
#split[27:42]<-4

#ha = HeatmapAnnotation(
#  empty = anno_empty(border = FALSE),
#  foo = anno_block(gp = gpar(fill = 2:6), labels = LETTERS[1:5])
#)
#Heatmap(effect_size_phen_sig_plot, name = "Effect size", row_split = split, right_annotation = ha, 
#        column_title = NULL)


#text_list = list(
#  text1 = "Pets",
#  text2 = "Farm animals",
#  text3 = "Wild animals",
#  text4 = "Food"
#)
# note how we set the width of this empty annotation
#ha = rowAnnotation(foo = anno_empty(border = FALSE, 
  #                                  width = max_text_width(unlist(text_list)) + unit(4, "mm")))
#Heatmap(effect_size_phen_sig_plot,cluster_row_slices  = F,cluster_rows = F, name = "Effect size", row_km = 4, right_annotation = ha)
#for(i in 1:4) {
#  decorate_annotation("foo", slice = i, {
#    grid.rect(x = 0, width = unit(2, "mm"), gp = gpar(fill = i, col = NA), just = "left")
#    grid.text(paste(text_list[[i]], collapse = "\n"), x = unit(4, "mm"), just = "left")
# })
#}


##################
##Phylum and known vs unknown

colnames(effect_size_phen_lmer)[c(45:47,51:59,49,60,50,61:68,70,73:88)]#[c(123:126,97,120,188:195,218:222,71:72,122,217,121,119,134,136,140,143,146,150,152,159,165,174,179:180,182:186)]#[c(45:47,51:59,49,60,50,61:68,70,73:88)]#[c(1,4,5,9,12,13,6,14:18,20:23,25:31,93:96,33:43)]
#1-34 --- measurements+chronic+medication
pval_chk_all<-rowSums(pval_chk[,c(45:47,51:59,49,60,50,61:68,70,73:88)])#,131:137
length(which(pval_chk_all>0))


fdr_chk_all<-rowSums(fdr_chk[,c(45:47,51:59,49,60,50,61:68,70,73:88)])
length(which(fdr_chk_all>0))#was 2

effect_size_phen_sig<-effect_size_phen_2[which(fdr_chk_all>0),c(45:47,51:59,49,60,50,61:68,70,73:88)]
pval_phen_sig<-pval_chk[which(fdr_chk_all>0),c(45:47,51:59,49,60,50,61:68,70,73:88)]
fdr_phen_sig<-fdr_chk[which(fdr_chk_all>0),c(45:47,51:59,49,60,50,61:68,70,73:88)]

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


effect_size_phen_sig_plot<-t(effect_size_phen_sig_plot)
rownames(effect_size_phen_sig_plot)<-colnames(effect_size_phen_2)[c(45:47,51:59,49,60,50,61:68,70,73:88)]
colnames(effect_size_phen_sig_plot)<-as.character(mb_samp_sp_name[which(fdr_chk_all>0)])# was as.character(mb_samp_sp_name[which(pval_chk_all>5)])

sp_ph_names<-colnames(effect_size_phen_sig_plot)

circ_names<-array(0,dim=c(length(sp_ph_names)))

for(i in 1:length(sp_ph_names)[1]){
  circ_names[i]<-ifelse(grepl("{",sp_ph_names[i],fixed=T),1,0)
}
#View(circ_names)
#View(sp_ph_names)
library(ggplot2)
library(dplyr)

slices <- c(length(which(circ_names==1)),length(which(circ_names==0)))
lbls <- c("Unknown species", "Known species")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="Pie Chart of Countries")

mycols2 <- c("#0073C2FF", "#CD534CFF")

count.data <- data.frame(
  class = lbls,
  n = slices,
  prop = pct,
  label2 = c(paste0("Unknown"," (",pct[1],"%)"),paste0("Known"," (",pct[2],"%)"))
)
count.data <- count.data %>%
  arrange(desc(class)) %>%
  mutate(lab.ypos = cumsum(prop) - 0.5*prop)
#count.data

ggplot(count.data, aes(x = "", y = prop, fill = class)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+
  geom_text(aes(y = lab.ypos, label = label2), color = "white")+
  scale_fill_manual(values = mycols2) +
  theme_void()

## Phylum

phylum_ind<-which((grepl("c__",mb_samp[,1],fixed=T)==F)&(grepl("p__",mb_samp[,1],fixed=T)==T))

mb_samp_phylum<-as.data.frame(mb_samp[phylum_ind,1])
#write.csv(mb_samp_phylum,'Phylum_list.csv')


ind_sp_ph<-array(0,dim=c(length(sp_ph_names)))

for(i in 1:length(sp_ph_names)){
  for(j in 1:dim(mb_samp)[1]){
    ind_temp1<-unlist(gregexpr("(",as.character(sp_ph_names[i]),fixed=T))
    ind_temp2<-unlist(gregexpr(")",as.character(sp_ph_names[i]),fixed=T))
    if(grepl(substr(sp_ph_names[i],ind_temp1+1,ind_temp2-1),as.character(mb_samp[j,1]),fixed=T)==T){
      ind_sp_ph[i]<-j
    }
  }
}

sp_ph_full_names<-as.character(mb_samp[ind_sp_ph,1])

phyl_ph<-array(0,dim=c(dim(mb_samp_phylum)[1]))

for(i in 1:length(phyl_ph)){
  for(j in 1:length(sp_ph_full_names)){
    if((grepl(as.character(mb_samp_phylum[i,1]),sp_ph_full_names[j],fixed=T))==T){
      
      
      phyl_ph[i]<-phyl_ph[i]+1
    }
  }
}




phsm<-cbind(as.character(mb_samp_phylum[,1]),phyl_ph)
phsm<-cbind(phsm,array("Food & animal factors",dim=c(length(phyl_ph),1)))#array("health",dim=c(length(phyl_ph),1))#array("Socio-economic factors",dim=c(length(phyl_ph),1))#array("Health factors",dim=c(length(phyl_ph),1))
phsm<-as.data.frame(phsm)
colnames(phsm)<-c("Phylum","count","Category")

ggplot(data = phsm, aes(x = as.character(Category), y = as.numeric(as.character(count)), fill = as.character(Phylum))) + 
  geom_bar(stat='identity') + coord_flip()

##################

#phsm_all<-phsm
#phsm_all<-rbind(phsm_all,phsm)
#phsm_all<-rbind(phsm_all,phsm)

#phsm_all$Category <- factor(phsm_all$Category, levels=c("Health factors","Food & animal factors","Socio-economic factors"),ordered=TRUE)
phsm_all2<-rbind(phsm_all[26:dim(phsm_all)[1],],phsm_all[1:25,])

#phsm_all2$Category <- factor(phsm_all2$Category, levels=c("Health factors","Food & animal factors","Socio-economic factors"),ordered=TRUE)
ggplot(data = phsm_all2, aes(x = as.character(Category), y = as.numeric(as.character(count)), fill = Phylum)) + 
  geom_bar(stat='identity') + coord_flip()



length(which((grepl("g__",as.character(mb_samp[,1]),fixed=T)==T)&!(grepl("s__",as.character(mb_samp[,1]),fixed=T)==T)))




########################################
###Finding phylum names

mb_samp_phylum

mb_samp_phylum<-as.data.frame(mb_samp[phylum_ind,1])
mb_samp_phylum_abun<-as.data.frame(mb_samp[phylum_ind,3:dim(mb_samp)[2]])

mb_samp_phylum_abun<-as.data.frame(mb_samp_phylum_abun)
rownames(mb_samp_phylum_abun)<-mb_samp_phylum[,1]
mb_samp_phylum_row<-as.data.frame(rowSums(mb_samp_phylum_abun))



#############################################
##Comparison with Dutch paper

effect_size_phen_lmer<-read.csv('effect_size_phen_1_lmer.csv',row.names = 1)
pval_phen_lmer<-read.csv('pval_phen_1_lmer.csv',row.names = 1)
fdr_phen_lmer<-read.csv('fdr_phen_1_lmer.csv',row.names = 1)

effect_size_phen_lmer<-cbind(effect_size_phen_lmer,read.csv('effect_size_phen_2_lmer.csv',row.names = 1))
pval_phen_lmer<-cbind(pval_phen_lmer,read.csv('pval_phen_2_lmer.csv',row.names = 1))
fdr_phen_lmer<-cbind(fdr_phen_lmer,read.csv('fdr_phen_2_lmer.csv',row.names = 1))

effect_size_phen_lmer<-cbind(effect_size_phen_lmer,read.csv('effect_size_phen_3_lmer.csv',row.names = 1))
pval_phen_lmer<-cbind(pval_phen_lmer,read.csv('pval_phen_3_lmer.csv',row.names = 1))
fdr_phen_lmer<-cbind(fdr_phen_lmer,read.csv('fdr_phen_3_lmer.csv',row.names = 1))

effect_size_phen_lmer<-cbind(effect_size_phen_lmer,read.csv('effect_size_phen_4_lmer.csv',row.names = 1))
pval_phen_lmer<-cbind(pval_phen_lmer,read.csv('pval_phen_4_lmer.csv',row.names = 1))
fdr_phen_lmer<-cbind(fdr_phen_lmer,read.csv('fdr_phen_4_lmer.csv',row.names = 1))

effect_size_phen_lmer<-cbind(effect_size_phen_lmer,read.csv('effect_size_phen_6_lmer.csv',row.names = 1))
pval_phen_lmer<-cbind(pval_phen_lmer,read.csv('pval_phen_6_lmer.csv',row.names = 1))
fdr_phen_lmer<-cbind(fdr_phen_lmer,read.csv('fdr_phen_6_lmer.csv',row.names = 1))

effect_size_phen_lmer<-cbind(effect_size_phen_lmer,read.csv('effect_size_phen_5_lmer.csv',row.names = 1))
pval_phen_lmer<-cbind(pval_phen_lmer,read.csv('pval_phen_5_lmer.csv',row.names = 1))
fdr_phen_lmer<-cbind(fdr_phen_lmer,read.csv('fdr_phen_5_lmer.csv',row.names = 1))

effect_size_phen_lmer<-effect_size_phen_lmer[1:2285,]
pval_phen_lmer<-pval_phen_lmer[1:2285,]
fdr_phen_lmer<-fdr_phen_lmer[1:2285,]

#effect_size_phen_lmer<-effect_size_phen_lmer[ind_thresh_use,]#Apply it for atleast 100 ppl and above 0.0001
#pval_phen_lmer<-pval_phen_lmer[ind_thresh_use,]
#fdr_phen_lmer<-fdr_phen_lmer[ind_thresh_use,]
#mb_samp_sp_name2<-mb_samp_sp_name[ind_thresh_use]


colnames(effect_size_phen_lmer)[16:26]<-colnames(effect_size_phen_2)[14:24]
colnames(effect_size_phen_lmer)[27:30]<-colnames(effect_size_phen_2)[25:28]
colnames(effect_size_phen_lmer)[85:86]<-c("Dementia","Cognitive impairment")
colnames(effect_size_phen_lmer)[152:154]<-c("Mild (GAD7)","Moderate (GAD7)","Severe (GAD7)")
colnames(effect_size_phen_lmer)[156:158]<-c("Mild (PHQ9)","Moderate (PHQ9)","Severe (PHQ9)")
colnames(effect_size_phen_lmer)[30:33]<-c("Anti-parasitic(34)","Anti-fungal(82)","Vitamins(121)","Anti-hypertensive(72)")
colnames(effect_size_phen_lmer)[9:11]<-c("Oxygen saturation","Hb total","Cough (1 month)")
colnames(effect_size_phen_lmer)[1]<-"Hb A1c"
colnames(effect_size_phen_lmer)[87:102]<-colnames(effect_size_phen_2)[c(99:105,108:116)]
colnames(effect_size_phen_lmer)[164:167]<-c("Alcohol daily frequency","Alcohol 6 drinks","Cigarette usage(46)","Cigarette frequency")

colnames(effect_size_phen_lmer)[200:225]<-colnames(effect_size_phen_2)[45:70]#Change
colnames(effect_size_phen_lmer)[237:240]<-c("Natural juice","Chicken","Beef/pork","Ham/sausages/hotdog")#Change
colnames(effect_size_phen_lmer)[232]<-"Cream/butter"#Change

rownames(effect_size_phen_lmer)<-mb_samp_sp_name
write.csv(effect_size_phen_lmer,'effect_size_phen_lmer_dmp_comp.csv')

#1603-1606
length(which(fdr_phen_lmer[1606,]<0.05))
rownames(effect_size_phen_lmer)[1782]

ind_dmp_sp<-c(1606,1782,1580,1577,1636,1652,747,44,39,41,45,1934,1589,2069,1960,2226,674,755,1405,1644,2117,337,1924,1417,1578,179,185,178,196,187,2110,1428,1429,134,1652,1786,1782,1886,743,1311,1953,44,1958,136,1481,1418,39,378,324,264,379)

ind_dmp_sp<-unique(ind_dmp_sp)


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


fdr_chk_lmer_all<-rowSums(fdr_chk_lmer[ind_dmp_sp,c(227,230:233,164,166,161:162,12:15,200:202,1,4,16,19,22:24,94,95,97,156:158)])
length(which(fdr_chk_lmer_all>0))

effect_size_phen_sig<-effect_size_phen_lmer[which(fdr_chk_lmer_all>0),c(227,230:233,164,166,161:162,12:15,200:202,1,4,16,19,22:24,94,95,97,156:158)]
pval_phen_sig<-pval_chk_lmer[which(fdr_chk_lmer_all>0),c(227,230:233,164,166,161:162,12:15,200:202,1,4,16,19,22:24,94,95,97,156:158)]
fdr_phen_sig<-fdr_chk_lmer[which(fdr_chk_lmer_all>0),c(227,230:233,164,166,161:162,12:15,200:202,1,4,16,19,22:24,94,95,97,156:158)]

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
rownames(effect_size_phen_sig_plot)<-colnames(effect_size_phen_lmer)[c(227,230:233,164,166,161:162,12:15,200:202,1,4,16,19,22:24,94,95,97,156:158)]
mb_samp_sp_name_temp<-mb_samp_sp_name[ind_dmp_sp]
colnames(effect_size_phen_sig_plot)<-as.character(mb_samp_sp_name_temp[which(fdr_chk_lmer_all>0)])
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

View(effect_size_phen_sig_plot)

##Combining phenotypes to be consistent with DMP

effect_size_phen_sig_plot_dmp<-array(0,dim=c(16,dim(effect_size_phen_sig_plot)[2]))

rownames(effect_size_phen_sig_plot_dmp)<-c("Carb intake","Fat intake","Alcohol intake","Cigarette use","Income","General health","Pets",rownames(effect_size_phen_sig_plot)[17:23],"Burnout","Depression")
colnames(effect_size_phen_sig_plot_dmp)<-colnames(effect_size_phen_sig_plot)

rownames(effect_size_phen_sig_plot_dmp)
colnames(effect_size_phen_sig_plot_dmp)

#effect_size_phen_sig_plot_dmp[2:3,17]<-1#+1 for same side, -1 for opposite
effect_size_phen_sig_plot_dmp[3,22]<-1#
effect_size_phen_sig_plot_dmp[5,1]<-1#
effect_size_phen_sig_plot_dmp[11,37]<-1*-1#
effect_size_phen_sig_plot_dmp[16,33]<-1#
effect_size_phen_sig_plot_dmp[1:3,17]<-1#
effect_size_phen_sig_plot_dmp[6,18]<-1#
effect_size_phen_sig_plot_dmp[5,13]<-1*-1
effect_size_phen_sig_plot_dmp[5,4]<-1*-1
effect_size_phen_sig_plot_dmp[13,8]<-1*-1
effect_size_phen_sig_plot_dmp[7,26]<-1*-1

library(corrplot)

rownames(effect_size_phen_sig_plot_dmp)[7]<-"Pets(1099)"
#effect_size_phen_sig_plot_dmp_plot<-effect_size_phen_sig_plot_dmp[c(1:3,5,6,7,11,13,16),c(1,4,8,13,17,18,22,26,33,37)]

#corrplot(as.matrix(effect_size_phen_sig_plot_dmp_plot),method="square",tl.cex = 0.6,cl.cex = 0,tl.col = "black",order= "hclust")#was tl.cex=0.8,cl.cex=0.85; 0.7,0.2
pheatmap(effect_size_phen_sig_plot_dmp[c(1:3,5,6,7,11,13,16),c(1,4,8,10,13,17,18,22,26,33,37)])

corrplot(as.matrix(effect_size_phen_sig_plot_dmp[c(1:3,5,6,7,11,13,16),c(1,4,8,10,13,17,18,22,26,33,37)]),method="square",tl.cex = 0.6,cl.cex = 0,tl.col = "black")

##################################
##Variance explained

library(ggplot2)

pll<-read.csv('pl_perm2_all.csv',row.names=1)
pll_pwy<-read.csv('pl_perm2_all_pwy.csv',row.names=1)
rownames(pll)
tech_factors_ind<-c(1,2,3,170,171)
phys_ind<-c(14:22)
acute_ind<-c(3,156)
chronic_ind<-c(4:13)
med_ind<-c(42:48)
person_ind<-c(49:58)
cog_ind<-c(59)
anx_ind<-c(164:166)
dep_ind<-c(167:169)
unfav_ind<-c(157:159)

health_var<-c(sum(pll$R2[phys_ind]),sum(pll$R2[acute_ind]),sum(pll$R2[chronic_ind]),sum(pll$R2[med_ind]),sum(pll$R2[person_ind]),sum(pll$R2[cog_ind]),sum(pll$R2[unfav_ind]),sum(pll$R2[anx_ind]),sum(pll$R2[dep_ind]))
#sum(health_var)
food_an_var<-sum(pll$R2[c(23:41,60:85)])
soc_eco_var<-sum(pll$R2[c(86:155,160:163)])
tech_var<-sum(pll$R2[tech_factors_ind])

health_var<-health_var*100
food_an_var<-food_an_var*100
soc_eco_var<-soc_eco_var*100
tech_var<-tech_var*100

#Pathway

health_var_pwy<-c(sum(pll_pwy$R2[phys_ind]),sum(pll_pwy$R2[acute_ind]),sum(pll_pwy$R2[chronic_ind]),sum(pll_pwy$R2[med_ind]),sum(pll_pwy$R2[person_ind]),sum(pll_pwy$R2[cog_ind]),sum(pll_pwy$R2[unfav_ind]),sum(pll_pwy$R2[anx_ind]),sum(pll_pwy$R2[dep_ind]))

food_an_var_pwy<-sum(pll_pwy$R2[c(23:41,60:85)])
soc_eco_var_pwy<-sum(pll_pwy$R2[c(86:155,160:163)])
tech_var_pwy<-sum(pll_pwy$R2[tech_factors_ind])

health_var_pwy<-health_var_pwy*100
food_an_var_pwy<-food_an_var_pwy*100
soc_eco_var_pwy<-soc_eco_var_pwy*100
tech_var_pwy<-tech_var_pwy*100

#inDFss<-cbind(c("Physiological","Acute","Chronic","Medication","Personality","Cognitive","Unfavorable","Anxiety","Depression","Food & animals","Socio-economic","Technical factors"),c(health_var,food_an_var,soc_eco_var,tech_var))
inDFss<-cbind(c("Physiological measurements","Acute conditions","Chronic conditions","Medications","Personalities","Cognitive","Unfavorable habits","Anxiety","Depression","Food & animals factors","Socio-economic factors","Technical factors"),c(health_var,food_an_var,soc_eco_var,tech_var))
inDFss<-as.data.frame(inDFss)
colnames(inDFss)<-c("Data","R2")
inDFss<-cbind(inDFss,"Type"=array("Species",dim=c(dim(inDFss)[1],1)))

inDFss_pwy<-cbind(c("Physiological measurements","Acute conditions","Chronic conditions","Medications","Personalities","Cognitive","Unfavorable habits","Anxiety","Depression","Food & animals factors","Socio-economic factors","Technical factors"),c(health_var_pwy,food_an_var_pwy,soc_eco_var_pwy,tech_var_pwy))
inDFss_pwy<-as.data.frame(inDFss_pwy)
colnames(inDFss_pwy)<-c("Data","R2")
inDFss_pwy<-cbind(inDFss_pwy,"Type"=array("Pathways",dim=c(dim(inDFss_pwy)[1],1)))

inDFss_all<-rbind(inDFss,inDFss_pwy)

smm<-array(NA,dim=c(dim(inDFss_all)[1]))
smm2<-array(NA,dim=c(dim(inDFss_all)[1]))
for(i in 1:dim(inDFss_all)[1]){
  smm[i]<-round(as.numeric(as.character(inDFss_all$R2[i])),2)
  smm2[i]<-paste0(round(as.numeric(as.character(inDFss_all$R2[i])),2),"%")
}

inDFss_all$R2<-smm

levels(inDFss_all$Data)

#g <- ggplot(inDFss_all,aes(x = Type, y=R2,fill=Data))  + geom_bar(position="stack",stat = "identity") + 
#  geom_text(aes(label = smm2), color = "black") + theme_classic()#health_var
#g

##Ordered

inDFss_all2<-inDFss_all

inDFss_all2$Data <- factor(inDFss_all2$Data, levels=c("Physiological measurements","Acute conditions","Chronic conditions","Medications","Personalities","Cognitive","Unfavorable habits","Anxiety","Depression","Food & animals factors","Socio-economic factors","Technical factors"),ordered=TRUE)

y_pos1<-array(NA,dim=c(dim(inDFss)[1]))
for(i in 1:dim(inDFss)[1]){
  #y_pos1[i]<-sum(as.numeric(as.character(inDFss$R2[1:i])))
  if(i==1){
    y_pos1[i]<-sum(as.numeric(as.character(inDFss$R2)))-as.numeric(as.character(inDFss$R2[i]))/2
  }else{
    y_pos1[i]<-sum(as.numeric(as.character(inDFss$R2)))-sum(as.numeric(as.character(inDFss$R2[1:(i-1)])))-as.numeric(as.character(inDFss$R2[i]))/2
  }
}
y_pos2<-array(NA,dim=c(dim(inDFss_pwy)[1]))
for(i in 1:dim(inDFss_pwy)[1]){
  #y_pos2[i]<-sum(as.numeric(as.character(inDFss_pwy$R2)))-sum(as.numeric(as.character(inDFss_pwy$R2[1:i])))
  if(i==1){
    y_pos2[i]<-sum(as.numeric(as.character(inDFss_pwy$R2)))-as.numeric(as.character(inDFss_pwy$R2[i]))/2
  }else{
    y_pos2[i]<-sum(as.numeric(as.character(inDFss_pwy$R2)))-sum(as.numeric(as.character(inDFss_pwy$R2[1:(i-1)])))-as.numeric(as.character(inDFss_pwy$R2[i]))/2
  }
}
#y_pos12<-c(y_pos1[length(y_pos1):-1:1],y_pos2[length(y_pos2):-1:1])
y_pos12<-c(y_pos1,y_pos2)
x_pos1<-array(0.8,dim=c(dim(y_pos1)[1]))
x_pos2<-array(1.8,dim=c(dim(y_pos2)[1]))
x_pos12<-c(x_pos1,x_pos2)

#g <- ggplot(inDFss_all2,aes(x = Type, y=R2,fill=Data))  + geom_bar(position="stack",stat = "identity",width=0.2) + 
#  geom_text(aes(label = smm2,y=y_pos12+1,x=x_pos12), color = "black")#+scale_fill_manual() #+ theme_classic()#health_var
#g

#ind_c<-c(1:9,12,13:21,24)
#y_pos123<-y_pos12[ind_c]
#g <- ggplot(inDFss_all2[ind_c,],aes(x = Type, y=R2,fill=Data))  + geom_bar(position="stack",stat = "identity",width=0.2) + 
#  geom_text(aes(label = smm2[ind_c],y=y_pos123+1), color = "black") + theme_classic()#health_var
#g
#y_pos123<-c(25.773537,22.856851,21.933466,21.216413,20.029349,19.756831,19.442583,19.175680,18.839253,12.268941,2.633531,0.000000,32.292183,30.198619,27.012991,25.854255,25.277434,24.465598,24.213648,24.035125,23.660068,14.968332,4.570267,0.000000)
y_pos123<-y_pos12
x_pos123<-x_pos12
y_pos123[6:9]<-c(20,19,18,17)
x_pos123[6:9]<-0.66
x_pos123[c(1:5,10:12)]<-0.82

y_pos123[18:21]<-c(25,24,23,22)
x_pos123[18:21]<-1.66
x_pos123[c(13:17,22:24)]<-1.82


g <- ggplot(inDFss_all2,aes(x = Type, y=R2,fill=Data))  + geom_bar(position="stack",stat = "identity",width=0.2) + 
  geom_text(aes(label = smm2,y=y_pos123,x=x_pos123),size=5, color = "black")+theme(text = element_text(size=20),panel.grid.major = element_blank(),panel.background = element_blank())+scale_fill_manual(values=c("#377eb8","#e5d8bd","#33a02c","#6a3d9a","#a50026","#dd3497","#dfc27d","#ff7f00","#000000","#525252","#a6bddb","#cab2d6"))+
  ylab("Variance explained (R2)")#+ theme_classic()#health_var
g

#############################
###Food & animals

tech_factors_ind<-c(1,2,3,170,171)
pet_ind<-c(60:62,66)
farm_ind<-c(64,67:75)
wild_ind<-c(63,65,76:85)
food_ind<-c(23:41)

health_var<-sum(pll$R2[c(3:22,42:59,157:159,164:169)])
#sum(health_var)
food_an_var<-c(sum(pll$R2[pet_ind]),sum(pll$R2[farm_ind]),sum(pll$R2[wild_ind]),sum(pll$R2[food_ind]))
soc_eco_var<-sum(pll$R2[c(86:155,160:163)])
tech_var<-sum(pll$R2[tech_factors_ind])

health_var<-health_var*100
food_an_var<-food_an_var*100
soc_eco_var<-soc_eco_var*100
tech_var<-tech_var*100

#Pathway

health_var_pwy<-sum(pll_pwy$R2[c(3:22,42:59,157:159,164:169)])
food_an_var_pwy<-c(sum(pll_pwy$R2[pet_ind]),sum(pll_pwy$R2[farm_ind]),sum(pll_pwy$R2[wild_ind]),sum(pll_pwy$R2[food_ind]))
#food_an_var_pwy<-sum(pll_pwy$R2[c(23:41,60:85)])
soc_eco_var_pwy<-sum(pll_pwy$R2[c(86:155,160:163)])
tech_var_pwy<-sum(pll_pwy$R2[tech_factors_ind])

health_var_pwy<-health_var_pwy*100
food_an_var_pwy<-food_an_var_pwy*100
soc_eco_var_pwy<-soc_eco_var_pwy*100
tech_var_pwy<-tech_var_pwy*100

#inDFss<-cbind(c("Physiological","Acute","Chronic","Medication","Personality","Cognitive","Unfavorable","Anxiety","Depression","Food & animals","Socio-economic","Technical factors"),c(health_var,food_an_var,soc_eco_var,tech_var))
inDFss<-cbind(c("Health factors","Pets","Farm animals","Wild animals","Food","Socio-economic factors","Technical factors"),c(health_var,food_an_var,soc_eco_var,tech_var))
inDFss<-as.data.frame(inDFss)
colnames(inDFss)<-c("Data","R2")
inDFss<-cbind(inDFss,"Type"=array("Species",dim=c(dim(inDFss)[1],1)))

inDFss_pwy<-cbind(c("Health factors","Pets","Farm animals","Wild animals","Food","Socio-economic factors","Technical factors"),c(health_var_pwy,food_an_var_pwy,soc_eco_var_pwy,tech_var_pwy))
inDFss_pwy<-as.data.frame(inDFss_pwy)
colnames(inDFss_pwy)<-c("Data","R2")
inDFss_pwy<-cbind(inDFss_pwy,"Type"=array("Pathways",dim=c(dim(inDFss_pwy)[1],1)))

inDFss_all<-rbind(inDFss,inDFss_pwy)

smm<-array(NA,dim=c(dim(inDFss_all)[1]))
smm2<-array(NA,dim=c(dim(inDFss_all)[1]))
for(i in 1:dim(inDFss_all)[1]){
  smm[i]<-round(as.numeric(as.character(inDFss_all$R2[i])),2)
  smm2[i]<-paste0(round(as.numeric(as.character(inDFss_all$R2[i])),2),"%")
}

inDFss_all$R2<-smm

levels(inDFss_all$Data)

#g <- ggplot(inDFss_all,aes(x = Type, y=R2,fill=Data))  + geom_bar(position="stack",stat = "identity") + 
#  geom_text(aes(label = smm2), color = "black") + theme_classic()#health_var
#g

##Ordered

inDFss_all2<-inDFss_all

inDFss_all2$Data <- factor(inDFss_all2$Data, levels=c("Health factors","Pets","Farm animals","Wild animals","Food","Socio-economic factors","Technical factors"),ordered=TRUE)

y_pos1<-array(NA,dim=c(dim(inDFss)[1]))
for(i in 1:dim(inDFss)[1]){
  #y_pos1[i]<-sum(as.numeric(as.character(inDFss$R2[1:i])))
  if(i==1){
    y_pos1[i]<-sum(as.numeric(as.character(inDFss$R2)))-as.numeric(as.character(inDFss$R2[i]))/2
  }else{
    y_pos1[i]<-sum(as.numeric(as.character(inDFss$R2)))-sum(as.numeric(as.character(inDFss$R2[1:(i-1)])))-as.numeric(as.character(inDFss$R2[i]))/2
  }
}
y_pos2<-array(NA,dim=c(dim(inDFss_pwy)[1]))
for(i in 1:dim(inDFss_pwy)[1]){
  #y_pos2[i]<-sum(as.numeric(as.character(inDFss_pwy$R2)))-sum(as.numeric(as.character(inDFss_pwy$R2[1:i])))
  if(i==1){
    y_pos2[i]<-sum(as.numeric(as.character(inDFss_pwy$R2)))-as.numeric(as.character(inDFss_pwy$R2[i]))/2
  }else{
    y_pos2[i]<-sum(as.numeric(as.character(inDFss_pwy$R2)))-sum(as.numeric(as.character(inDFss_pwy$R2[1:(i-1)])))-as.numeric(as.character(inDFss_pwy$R2[i]))/2
  }
}
#y_pos12<-c(y_pos1[length(y_pos1):-1:1],y_pos2[length(y_pos2):-1:1])
y_pos12<-c(y_pos1,y_pos2)
x_pos1<-array(0.8,dim=c(dim(y_pos1)[1]))
x_pos2<-array(1.8,dim=c(dim(y_pos2)[1]))
x_pos12<-c(x_pos1,x_pos2)

g <- ggplot(inDFss_all2,aes(x = Type, y=R2,fill=Data))  + geom_bar(position="stack",stat = "identity",width=0.2) + 
  geom_text(aes(label = smm2,y=y_pos12,x=x_pos12),size=5, color = "black")+theme(text = element_text(size=20),panel.grid.major = element_blank(),panel.background = element_blank())+scale_fill_manual(values=c("#525252","#377eb8","#d73027","#f46d43","#000000","#a6bddb","#cab2d6"))+
  ylab("Variance explained (R2)")#+ theme_classic()#health_var
g

###########################
###Socio economic factors

tech_factors_ind<-c(1,2,3,170,171)
partner_ind<-c(94:96)
edu_ind<-c(160:162)
friend_ind<-c(88:91)
family_ind<-c(140:143)
all_ind<-c(145:149)
income_ind<-c(86:87)
risky_ind<-c(144,93,163)
village_ind<-c(152:154)
hh_ind<-c(151,97:131)
water_ind<-c(132:139)

health_var<-sum(pll$R2[c(3:22,42:59,157:159,164:169)])
#sum(health_var)
#food_an_var<-c(sum(pll$R2[pet_ind]),sum(pll$R2[farm_ind]),sum(pll$R2[wild_ind]),sum(pll$R2[food_ind]))
food_an_var<-sum(pll$R2[c(23:41,60:85)])
soc_eco_var<-c(sum(pll$R2[partner_ind]),sum(pll$R2[edu_ind]),sum(pll$R2[friend_ind]),sum(pll$R2[family_ind]),sum(pll$R2[all_ind]),sum(pll$R2[risky_ind]),sum(pll$R2[village_ind]),sum(pll$R2[income_ind]),sum(pll$R2[hh_ind]),sum(pll$R2[water_ind]))
tech_var<-sum(pll$R2[tech_factors_ind])

health_var<-health_var*100
food_an_var<-food_an_var*100
soc_eco_var<-soc_eco_var*100
tech_var<-tech_var*100

#Pathway

health_var_pwy<-sum(pll_pwy$R2[c(3:22,42:59,157:159,164:169)])
#food_an_var_pwy<-c(sum(pll_pwy$R2[pet_ind]),sum(pll_pwy$R2[farm_ind]),sum(pll_pwy$R2[wild_ind]),sum(pll_pwy$R2[food_ind]))
food_an_var_pwy<-sum(pll_pwy$R2[c(23:41,60:85)])
soc_eco_var_pwy<-c(sum(pll_pwy$R2[partner_ind]),sum(pll_pwy$R2[edu_ind]),sum(pll_pwy$R2[friend_ind]),sum(pll_pwy$R2[family_ind]),sum(pll_pwy$R2[all_ind]),sum(pll_pwy$R2[risky_ind]),sum(pll_pwy$R2[village_ind]),sum(pll_pwy$R2[income_ind]),sum(pll_pwy$R2[hh_ind]),sum(pll_pwy$R2[water_ind]))
tech_var_pwy<-sum(pll_pwy$R2[tech_factors_ind])

health_var_pwy<-health_var_pwy*100
food_an_var_pwy<-food_an_var_pwy*100
soc_eco_var_pwy<-soc_eco_var_pwy*100
tech_var_pwy<-tech_var_pwy*100

#inDFss<-cbind(c("Physiological","Acute","Chronic","Medication","Personality","Cognitive","Unfavorable","Anxiety","Depression","Food & animals","Socio-economic","Technical factors"),c(health_var,food_an_var,soc_eco_var,tech_var))
inDFss<-cbind(c("Health factors","Food & animal factors","Co-habiting partners","Education","Friendship","Family","All relationships","Risky behavior","Village factors","Income","Household essentials","Water sources","Technical factors"),c(health_var,food_an_var,soc_eco_var,tech_var))
inDFss<-as.data.frame(inDFss)
colnames(inDFss)<-c("Data","R2")
inDFss<-cbind(inDFss,"Type"=array("Species",dim=c(dim(inDFss)[1],1)))

inDFss_pwy<-cbind(c("Health factors","Food & animal factors","Co-habiting partners","Education","Friendship","Family","All relationships","Risky behavior","Village factors","Income","Household essentials","Water sources","Technical factors"),c(health_var_pwy,food_an_var_pwy,soc_eco_var_pwy,tech_var_pwy))
inDFss_pwy<-as.data.frame(inDFss_pwy)
colnames(inDFss_pwy)<-c("Data","R2")
inDFss_pwy<-cbind(inDFss_pwy,"Type"=array("Pathways",dim=c(dim(inDFss_pwy)[1],1)))

inDFss_all<-rbind(inDFss,inDFss_pwy)

smm<-array(NA,dim=c(dim(inDFss_all)[1]))
smm2<-array(NA,dim=c(dim(inDFss_all)[1]))
for(i in 1:dim(inDFss_all)[1]){
  smm[i]<-round(as.numeric(as.character(inDFss_all$R2[i])),2)
  smm2[i]<-paste0(round(as.numeric(as.character(inDFss_all$R2[i])),2),"%")
}

inDFss_all$R2<-smm

levels(inDFss_all$Data)

#g <- ggplot(inDFss_all,aes(x = Type, y=R2,fill=Data))  + geom_bar(position="stack",stat = "identity") + 
#  geom_text(aes(label = smm2), color = "black") + theme_classic()#health_var
#g

##Ordered

inDFss_all2<-inDFss_all

inDFss_all2$Data <- factor(inDFss_all2$Data, levels=c("Health factors","Food & animal factors","Co-habiting partners","Education","Friendship","Family","All relationships","Risky behavior","Village factors","Income","Household essentials","Water sources","Technical factors"),ordered=TRUE)

y_pos1<-array(NA,dim=c(dim(inDFss)[1]))
for(i in 1:dim(inDFss)[1]){
  #y_pos1[i]<-sum(as.numeric(as.character(inDFss$R2[1:i])))
  if(i==1){
    y_pos1[i]<-sum(as.numeric(as.character(inDFss$R2)))-as.numeric(as.character(inDFss$R2[i]))/2
  }else{
    y_pos1[i]<-sum(as.numeric(as.character(inDFss$R2)))-sum(as.numeric(as.character(inDFss$R2[1:(i-1)])))-as.numeric(as.character(inDFss$R2[i]))/2
  }
}
y_pos2<-array(NA,dim=c(dim(inDFss_pwy)[1]))
for(i in 1:dim(inDFss_pwy)[1]){
  #y_pos2[i]<-sum(as.numeric(as.character(inDFss_pwy$R2)))-sum(as.numeric(as.character(inDFss_pwy$R2[1:i])))
  if(i==1){
    y_pos2[i]<-sum(as.numeric(as.character(inDFss_pwy$R2)))-as.numeric(as.character(inDFss_pwy$R2[i]))/2
  }else{
    y_pos2[i]<-sum(as.numeric(as.character(inDFss_pwy$R2)))-sum(as.numeric(as.character(inDFss_pwy$R2[1:(i-1)])))-as.numeric(as.character(inDFss_pwy$R2[i]))/2
  }
}
#y_pos12<-c(y_pos1[length(y_pos1):-1:1],y_pos2[length(y_pos2):-1:1])
y_pos12<-c(y_pos1,y_pos2)
x_pos1<-array(0.8,dim=c(dim(y_pos1)[1]))
x_pos2<-array(1.8,dim=c(dim(y_pos2)[1]))
x_pos12<-c(x_pos1,x_pos2)

y_pos123<-y_pos12
x_pos123<-x_pos12
y_pos123[3:10]<-c(14,13,12,11,10,9,8,7)
x_pos123[3:10]<-0.66
x_pos123[c(1:2,11:13)]<-0.82

y_pos123[16:23]<-c(16,15,14,13,12,11,10,9)
x_pos123[16:23]<-1.66
x_pos123[c(14:15,24:26)]<-1.82

g <- ggplot(inDFss_all2,aes(x = Type, y=R2,fill=Data))  + geom_bar(position="stack",stat = "identity",width=0.2) + 
  geom_text(aes(label = smm2,y=y_pos123,x=x_pos123),size=5, color = "black")+theme(text = element_text(size=20),panel.grid.major = element_blank(),panel.background = element_blank())+scale_fill_manual(values=c("#525252","#bdbdbd","#377eb8","#d73027","#f46d43","#ae017e","#006837","#49006a","#fcc5c0","#a6cee3","#fed976","#000000","#cab2d6"))+
  ylab("Variance explained (R2)")#+ theme_classic()#health_var
g


#################
##FIgure 1 phylum plotting



phylum_ind<-which((grepl("c__",mb_samp[,1],fixed=T)==F)&(grepl("p__",mb_samp[,1],fixed=T)==T))
mb_samp_phylum_temp<-as.data.frame(mb_samp[phylum_ind,])
lp<-0
for(i in 1:length(phylum_ind)){
  lp<-rbind(lp,length(which(mb_samp_phylum_temp[i,3:1190]>0)))
}
lp<-lp[2:length(lp)]

mb_samp_phylum<-as.data.frame(mb_samp_phylum_temp[which(lp>200),])

mb_samp_phylum22<-as.matrix(mb_samp_phylum[,3:dim(mb_samp_phylum)[2]])
rownames(mb_samp_phylum22)<-c("Euryarchaeota (Archaea)","Actinobacteria (Bacteria)","Bacteroidetes (Bacteria)","Candidatus Melainabacteria (Bacteria)","Candidatus Saccharibacteria (Bacteria)","Elusimicrobia (Bacteria)","Firmicutes (Bacteria)","Fusobacteria (Bacteria)","Kiritimatiellaeota (Bacteria)","Lentisphaerae (Bacteria)","Planctomycetes (Bacteria)","Proteobacteria (Bacteria)","Spirochaetes (Bacteria)","Tenericutes (Bacteria)","Verrucomicrobia (Bacteria)","Eukaryota unclassified (Eukaryota)")

mb_samp_phylum22$plot_bar(facet = "Group", xtext_keep = FALSE, legend_text_italic = FALSE)

ggplot(mb_samp_phylum22,aes(x = Type, y=R2,fill=Data))  + geom_bar(position="stack",stat = "identity",width=0.2) + 
  geom_text(aes(label = smm2,y=y_pos123,x=x_pos123),size=5, color = "black")+theme(text = element_text(size=20),panel.grid.major = element_blank(),panel.background = element_blank())+scale_fill_manual(values=c("#525252","#bdbdbd","#377eb8","#d73027","#f46d43","#ae017e","#006837","#49006a","#fcc5c0","#a6cee3","#fed976","#000000","#cab2d6"))+
  ylab("Variance explained (R2)")

#

for(j in 1:dim(mb_samp_phylum22)[2]){
    if(j==1){
      mb_samp_phylum23<-cbind(cbind(mb_samp_phylum22[,j],array(colnames(mb_samp_phylum22)[j],dim=c(dim(mb_samp_phylum22)[1],1))),rownames(mb_samp_phylum22))
    }else{
      mb_samp_phylum23<-rbind(mb_samp_phylum23,cbind(cbind(mb_samp_phylum22[,j],array(colnames(mb_samp_phylum22)[j],dim=c(dim(mb_samp_phylum22)[1],1))),rownames(mb_samp_phylum22)))
    }
}

colnames(mb_samp_phylum23)<-c("abundance","sample","phylum")
mb_samp_phylum24<-as.data.frame(mb_samp_phylum23)

ggplot(mb_samp_phylum24, aes(x=sample,y=as.numeric(as.character(abundance)), fill = as.factor(phylum))) + 
  geom_bar(position="stack",stat = "identity")   # Stacked 100% barplot
  #scale_fill_manual("", values=colours) +
  #theme(axis.text.x=element_text(angle=90, vjust=0.5)) +  # Vertical x-axis tick labels
  #scale_y_continuous(labels = scales::percent_format()) +
  #labs(y="Relative abundance")

ind_j<-0
mtt<-sort(mb_samp_phylum22[7,])

for(i in 1:dim(mb_samp_phylum22)[2]){
  ind_j<-rbind(ind_j,which(mb_samp_phylum22[7,]==mtt[i]))
}

for(j in 1:dim(mb_samp_phylum22)[2]){
  if(j==1){
    mb_samp_phylum26<-cbind(cbind(mb_samp_phylum22[,ind_j[j+1]],array(j,dim=c(dim(mb_samp_phylum22)[1],1))),rownames(mb_samp_phylum22))
  }else{
    mb_samp_phylum26<-rbind(mb_samp_phylum26,cbind(cbind(mb_samp_phylum22[,ind_j[j+1]],array(j,dim=c(dim(mb_samp_phylum22)[1],1))),rownames(mb_samp_phylum22)))
  }
}
mb_samp_phylum27<-as.data.frame(mb_samp_phylum26)
colnames(mb_samp_phylum27)<-c("abundance","sample","phylum")
mb_samp_phylum28<-mb_samp_phylum27[order(mb_samp_phylum27[,1]),]

ggplot(mb_samp_phylum28, aes(x=sample,y=as.numeric(as.character(abundance)), fill = as.factor(phylum))) + 
  geom_bar(position="stack",stat = "identity")   # Stacked 100% barplot


mb_samp_phylum27$phylum<-factor(mb_samp_phylum27$phylum,levels=c("Euryarchaeota (Archaea)","Actinobacteria (Bacteria)","Bacteroidetes (Bacteria)","Candidatus Melainabacteria (Bacteria)","Candidatus Saccharibacteria (Bacteria)","Elusimicrobia (Bacteria)","Firmicutes (Bacteria)","Fusobacteria (Bacteria)","Kiritimatiellaeota (Bacteria)","Lentisphaerae (Bacteria)","Planctomycetes (Bacteria)","Proteobacteria (Bacteria)","Spirochaetes (Bacteria)","Tenericutes (Bacteria)","Verrucomicrobia (Bacteria)","Eukaryota unclassified (Eukaryota)"),ordered=TRUE)

ggplot(mb_samp_phylum27,aes(x=as.array(as.character(sample)), y=as.numeric(as.character(abundance)),group=phylum, fill=phylum)) + 
  #geom_line()+
  geom_area() + 
  #geom_point() +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  xlab('Samples') + ylab("Relative abundance") +
  scale_fill_manual(values =colorRampPalette(RColorBrewer::brewer.pal(16, "Dark2"))(16)) + 
  #scale_fill_brewer(palette="Dark2") +
  ggtitle('Phylum-level profile of Honduras cohort')


tk3<-mb_samp_phylum28[order(-as.numeric(factor(mb_samp_phylum28$abundance))),]

ggplot(tk3, aes(x=sample,y=as.numeric(as.character(abundance)),group=as.factor(phylum), fill = as.factor(phylum))) + 
  geom_area()   # Stacked 100% barplot
###################################

a3<-mb_samp_phylum
a4<-as.numeric(unlist(a3[,3:dim(a3)[2]]))

tk<-array(0,dim=c(length(a4),3))

for(i in 1:(dim(a3)[2]-2)){
  for(j in 1:dim(a3)[1]){
    tk[i*16-(16-j),1]<-as.character(a3[j,1])
    #tk[i*16-(16-j),3]<-as.character(colnames(a3)[i+2])
    tk[i*16-(16-j),3]<-i
  }
}

tk<-array(0,dim=c(length(a4),3))

for(j in 1:dim(a3)[1]){
  for(i in 1:(dim(a3)[2]-2)){
    tk[(j-1)*(dim(a3)[2]-2)+i,1]<-as.character(a3[j,1])
    #tk[i*16-(16-j),3]<-as.character(colnames(a3)[i+2])
    tk[(j-1)*(dim(a3)[2]-2)+i,2]<-as.numeric(a3[j,i+2])
    tk[(j-1)*(dim(a3)[2]-2)+i,3]<-i
  }
}

#as.numeric(a3[3,5])


#tk[,2]<-ceiling(a4)/100#Need to do this

colnames(tk)<-c("tax","value","sample")
tk2<-as.data.frame(tk)

ggplot(tk2,aes(x=sample, y=as.numeric(as.character(value)),group=tax, fill=tax)) + 
  #geom_line()+
  geom_area() + 
  #geom_point() +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  xlab('Samples') + ylab("Relative abundance") +
  #scale_fill_manual(values =colorRampPalette(RColorBrewer::brewer.pal(16, "Dark2"))(16)) + 
  #scale_fill_brewer(palette="Dark2") +
  ggtitle('Phylum-level profile of Honduras cohort')

a5<-a3[3,3:dim(a3)[2]]
colnames(a5)<-c(1:length(a5))
a6<-sort(a5,decreasing = T)

b1<-a3
b1[3:dim(b1)[2]]<-a3[as.numeric(colnames(a6))+2]
colnames(b1)[3:dim(b1)[2]]<-colnames(a3)[as.numeric(colnames(a6))+2]

b2<-as.numeric(unlist(b1[,3:dim(a3)[2]]))

tk<-array(0,dim=c(length(b2),3))

for(j in 1:dim(b1)[1]){
  for(i in 1:(dim(b1)[2]-2)){
    tk[(j-1)*(dim(b1)[2]-2)+i,1]<-as.character(b1[j,1])
    tk[(j-1)*(dim(b1)[2]-2)+i,2]<-as.numeric(b1[j,i+2])
    #tk[(j-1)*(dim(b1)[2]-2)+i,3]<-colnames(b1[i+2])
    tk[(j-1)*(dim(b1)[2]-2)+i,3]<-i
  }
}
colnames(tk)<-c("tax","value","sample")
tk2<-as.data.frame(tk)

set.seed(1)
tk3<-tk2[order(-as.numeric(factor(tk2$sample,levels=c(1:(dim(b1)[2]-2))))),]

unique(tk3$tax)

p<-ggplot(tk3,aes(x=as.array(as.character(sample)), y=as.numeric(as.character(value)),group=tax, fill=tax)) + 
  #geom_line()+
  geom_area() + 
  #geom_point() +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  xlab('Samples') + ylab("Relative abundance") +
  scale_fill_manual(values =colorRampPalette(RColorBrewer::brewer.pal(16, "Dark2"))(16)) + 
  #scale_fill_brewer(palette="Dark2") +
  ggtitle('Phylum-level profile of Honduras cohort')

p

##

ggplot(tk3[c(3020:dim(tk3)[1]),],aes(x=as.array(as.character(sample)), y=as.numeric(as.character(value)),group=tax, fill=tax)) + 
  #geom_line()+
  geom_area() + 
  #geom_point() +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),text=element_text(size=20)) + 
  xlab('Samples(1187)') + ylab("Relative abundance") +
  scale_fill_manual(values =colorRampPalette(RColorBrewer::brewer.pal(16, "Dark2"))(16)) + 
  #scale_fill_brewer(palette="Dark2") +
  ggtitle('Phylum-level profile of Honduras cohort')

tk4<-tk3
unique(tk4$tax)


ggplot(tk3,aes(x=as.array(as.character(sample)), y=as.numeric(as.character(value)),group=tax, fill=tax)) + 
  #geom_line()+
  geom_area() + 
  #geom_point() +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),text=element_text(size=20)) + 
  xlab('Samples(1187)') + ylab("Relative abundance") +
  scale_fill_manual(values =colorRampPalette(RColorBrewer::brewer.pal(16, "Dark2"))(16)) + 
  #scale_fill_brewer(palette="Dark2") +
  ggtitle('Phylum-level profile of Honduras cohort')




###
##

effect_size_phen_lmer<-read.csv('effect_size_phen_healthy_lmer.csv',row.names=1)
pval_phen_lmer<-read.csv('pval_phen_healthy_lmer.csv',row.names=1)
fdr_phen_lmer<-read.csv('fdr_phen_healthy_lmer.csv',row.names=1)

effect_size_phen_lmer<-effect_size_phen_lmer[1:2285,]
pval_phen_lmer<-pval_phen_lmer[1:2285,]
fdr_phen_lmer<-fdr_phen_lmer[1:2285,]


#healthy vs unhealthy plotting

effect_size_phen_lmer2<-effect_size_phen_lmer[ind_thresh_use2[,2]]
pval_phen_lmer2<-pval_phen_lmer[ind_thresh_use2[,2]]
fdr_phen_lmer2<-fdr_phen_lmer[ind_thresh_use2[,2]]


length(which(fdr_phen_lmer2<0.05))
effect_size_phen_lmer2[which(fdr_phen_lmer2<0.05)]
fdr_phen_lmer2[which(fdr_phen_lmer2<0.05)]
mb_samp_sp_name_use2[which(fdr_phen_lmer2<0.05)]



efz<-effect_size_phen_lmer2[which(fdr_phen_lmer2<0.05)]
efz2<-efz[order(efz)]

fdr<-fdr_phen_lmer2[which(fdr_phen_lmer2<0.05)]
fdr2<-fdr[order(efz)]

efz_names<-mb_samp_sp_name_use2[which(fdr_phen_lmer2<0.05)]
efz_names2<-efz_names[order(efz)]

library(ggplot2)

health_bar<-cbind(cbind(cbind(efz2[c(1:5,32:36)],fdr2[c(1:5,32:36)]),efz_names2[c(1:5,32:36)]),c(array(1,dim=c(5,1)),array(2,dim=c(5,1))))
colnames(health_bar)<-c("effect","fdr","name","groupz")
health_bar<-as.data.frame(health_bar)

#levels(health_bar$name)<-factor(c("Prevotella_stercorea (SGB1679)","Coprococcus_sp_OM04_5BH (SGB5115)","Prevotella_pectinovora (SGB1662)","Prevotella_sp_Marseille_P4119 (SGB1676)","Lachnospiraceae_unclassified_SGB4924 (SGB4924_group)","Bacteroides_finegoldii (SGB1862)","Butyricimonas_virosa (SGB1784)","Barnesiella_intestinihominis (SGB1965)","Bacteroides_stercoris (SGB1830_group)","{p__Firmicutes}GGB3256_SGB4303 (SGB4303)"), levels=as.character(health_bar$name),ordered=TRUE)

#ggplot(health_bar,aes(x=c(1:10),y=as.numeric(as.character(effect)),fill=groupz))+
#  geom_bar(stat="identity",color="black")+scale_fill_manual(values = c("#4575b4","#d73027"))+#+stat_bin(geom="text",aes(),vjust=-1)
#  theme(text = element_text(size=20),panel.grid.major = element_blank(),panel.background = element_blank())+
#  coord_flip()+ylab("Effect size")#+ scale_x_discrete()

stp<-c("***","***","**","*","*","*","**","**","**","***")

ggplot(health_bar,aes(x=reorder(as.character(name),as.numeric(as.character(effect))),y=as.numeric(as.character(effect)),fill=groupz))+
  geom_bar(stat="identity",color="black")+scale_fill_manual(values = c("#4575b4","#d73027"))+#+stat_bin(geom="text",aes(),vjust=-1)
  theme(text = element_text(size=20),panel.grid.major = element_blank(),panel.background = element_blank())+ annotate("text", x=reorder(as.character(health_bar$name),as.numeric(as.character(health_bar$effect))),y=as.numeric(as.character(health_bar$effect)), label = stp,size=10,hjust=c(0,0,0,0,0,1,1,1,1,1)) +
  ylab("Effect size (Healthy vs Unhealthy)")+coord_flip()+xlab("")
  #stat_summary(fun = sum, geom="bar", label=stp,vjust = 0)#+stat_bin(geom="text",aes(y=stp),vjust=-1)

ggplot(health_bar,aes(x=reorder(as.character(name),as.numeric(as.character(effect))),y=as.numeric(as.character(effect)),fill=groupz))+
  geom_bar(stat="identity",color="black")+scale_fill_manual(values = c("#4575b4","#d73027"))+#+stat_bin(geom="text",aes(),vjust=-1)
  theme(axis.text.x = element_text(angle = 90),text = element_text(size=20),panel.grid.major = element_blank(),panel.background = element_blank())+ annotate("text", x=reorder(as.character(health_bar$name),as.numeric(as.character(health_bar$effect))),y=as.numeric(as.character(health_bar$effect)), label = stp,size=10,vjust=c(1,1,1,1,1,0.5,0.5,0.5,0.5,0.5)) +
  ylab("Effect size (Healthy vs Unhealthy)")+xlab("")


###Ran again just on 788 species--correctly -- 12-25-22
#Need to re run species comparison -- standardized -- 134 species significant

mb_samp_sp_name3<-mb_samp_sp_name[ind_thresh_use2[,2]]

tt1<-read.csv('results_maaslin_healthy_vs_unhealthy_stand.csv',row.names = 1)
spl<-array(NA,dim=c(dim(tt1)[1]))
for(i in 1:dim(tt1)[1]){
  for(j in 1:length(mb_samp_sp_name3)){
    ind_s1<-unlist(gregexpr("t__",as.character(tt1[i,1]),fixed=T))
    kll<-substr(tt1[i,1],ind_s1+3,nchar(as.character(tt1[i,1])))
    ind_s2<-unlist(gregexpr("(",as.character(mb_samp_sp_name3[j]),fixed=T))
    ind_s3<-unlist(gregexpr(")",as.character(mb_samp_sp_name3[j]),fixed=T))
    kll2<-substr(as.character(mb_samp_sp_name3[j]),ind_s2+1,ind_s3-1)
    if(kll==kll2){
      spl[i]<-j
    }
  }
}

tt2<-cbind(tt1,"sp_name"=as.character(mb_samp_sp_name3[spl[1:dim(tt1)[1]]]))
tt2

efz<-tt2$effect_size[which(tt2$FDR<0.05)]
efz_names<-as.character(tt2$sp_name[which(tt2$FDR<0.05)])
efz2<-efz[order(efz)]
efz_names2<-efz_names[order(efz)]

fdr<-tt2$FDR[which(tt2$FDR<0.05)]
fdr2<-fdr[order(efz)]

library(ggplot2)

#health_bar<-cbind(cbind(cbind(efz2[c(1:10,(length(efz)-9):length(efz))],fdr2[c(1:10,(length(efz)-9):length(efz))]),efz_names2[c(1:10,(length(efz)-9):length(efz))]),c(array(1,dim=c(10,1)),array(2,dim=c(10,1))))
health_bar<-cbind(cbind(cbind(efz2[c(1:10,(length(efz)-9):length(efz))],fdr2[c(1:10,(length(efz)-9):length(efz))]),as.character(c("Bacteroides ovatus (SGB1871)","Bacteroides caccae (SGB1877)","Parabacteroides distasonis (SGB1934)","Alistipes putredinis (SGB2318)","Bacteroides xylanisolvens (SGB1867)","Paraprevotella clara (SGB1798)","Phocaeicola vulgatus (SGB1814)","Bacteroides uniformis (SGB1836_group)","Bacteroides thetaiotaomicron (SGB1861)","Clostridium sp. AM22 11AC (SGB4749)","Prevotella sp. 885 (SGB1677)","{f__Bacilli unclassified}GGB4672 (SGB6461)","Bacilli unclassified (SGB6428)","{f__Bacteroidaceae}GGB1380 (SGB1883 group)","{f__Prevotellaceae}GGB1243 (SGB1663)","{f__Rikenellaceae}GGB1632 (SGB2240)","Coprococcus sp. (SGB5115)","{f__Rikenellaceae}GGB1631 (SGB2239)","{f__Prevotellaceae}GGB1247 (SGB1668)","{f__Rikenellaceae}GGB1630 (SGB2238)"))),c(array(1,dim=c(10,1)),array(2,dim=c(10,1))))
colnames(health_bar)<-c("effect","fdr","name","groupz")
health_bar<-as.data.frame(health_bar)

stp<-c("***","***","**","**","***","**","**","**","**","**","*","***","***","**","**","**","***","***","***","***")

ggplot(health_bar,aes(x=reorder(as.character(name),as.numeric(as.character(effect))),y=as.numeric(as.character(effect)),fill=groupz))+
  geom_bar(stat="identity",color="black",size=1)+scale_fill_manual(values = c("#4575b4","#d73027"))+#+stat_bin(geom="text",aes(),vjust=-1)
  theme(text = element_text(size=20),panel.grid.major = element_blank(),panel.background = element_blank()) +
  ylab("Effect size (Healthy vs disease associated)")+coord_flip()+xlab("")

#health_bar[1:20,3]<-as.character(c("Bacteroides ovatus (SGB1871)","Bacteroides caccae (SGB1877)","Parabacteroides distasonis (SGB1934)","Alistipes putredinis (SGB2318)","Bacteroides xylanisolvens (SGB1867)","Paraprevotella clara (SGB1798)","Phocaeicola vulgatus (SGB1814)","Bacteroides uniformis (SGB1836_group)","Bacteroides thetaiotaomicron (SGB1861)","Clostridium sp. AM22 11AC (SGB4749)","Prevotella sp. 885 (SGB1677)","{f__Bacilli unclassified}GGB4672 (SGB6461)","Bacilli unclassified (SGB6428)","{f__Bacteroidaceae}GGB1380 (SGB1883 group)","{f__Prevotellaceae}GGB1243 (SGB1663)","{f__Rikenellaceae}GGB1632 (SGB2240)","Coprococcus sp. (SGB5115)","{f__Rikenellaceae}GGB1631 (SGB2239)","{f__Prevotellaceae}GGB1247 (SGB1668)","{f__Rikenellaceae}GGB1630 (SGB2238)"))

#ggplot(health_bar,aes(x=reorder(as.character(name),as.numeric(as.character(effect))),y=as.numeric(as.character(effect)),fill=groupz))+
#  geom_bar(stat="identity",color="black")+scale_fill_manual(values = c("#4575b4","#d73027"))+#+stat_bin(geom="text",aes(),vjust=-1)
#  theme(text = element_text(size=20),panel.grid.major = element_blank(),panel.background = element_blank())+ annotate("text", x=reorder(as.character(health_bar$name),as.numeric(as.character(health_bar$effect))),y=as.numeric(as.character(health_bar$effect)), label = stp,size=10,hjust=c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1)) +
#  ylab("Effect size (Healthy vs Unhealthy)")+coord_flip()+xlab("")


ggplot(health_bar,aes(x=reorder(as.character(name),as.numeric(as.character(effect))),y=as.numeric(as.character(effect)),fill=groupz))+
  geom_bar(stat="identity",color="black")+scale_fill_manual(values = c("#4575b4","#d73027"))+#+stat_bin(geom="text",aes(),vjust=-1)
  theme(axis.text.x = element_text(angle = 90),text = element_text(size=20),panel.grid.major = element_blank(),panel.background = element_blank())+ annotate("text", x=reorder(as.character(health_bar$name),as.numeric(as.character(health_bar$effect))),y=as.numeric(as.character(health_bar$effect)), label = stp,size=10,vjust=c(1,1,1,1,1,0.5,0.5,0.5,0.5,0.5)) +
  ylab("Effect size (Healthy vs Unhealthy)")+xlab("")


##Adding Diet diversity score (DDS)

phen_all_use<-read.csv('phen_all_use6_all.csv',row.names = 1)
colnames(phen_all_use)
ind1<-c(76,77,78)#	Cereals
#colnames(phen_all_use)[ind1]
ind3<-c(84)#	Vegetables
#colnames(phen_all_use)[ind3]
ind4<-c(85,86)#	Fruits
ind5<-c(87:88)#	Meat, poultry, offal
#colnames(phen_all_use)[ind9]
ind6<-c(83)#	Eggs
ind7<-c(90)#	Fish and seafood
ind8<-c(75)#Pulses,legumes,nuts
ind9<-c(79:82)#Milk and milk products
ind10<-c(89,93)#Oil/fats
ind11<-c(91,92)#	Sugar/honey

dds<-array(NA,dim=c(dim(phen_all_use)[1],1))
for(i in 1:dim(phen_all_use)[1]){
  sc<-0
  temp<-0
  for(j in 1:length(ind1)){
    if(is.na(phen_all_use[i,ind1[j]])==F){
    if(phen_all_use[i,ind1[j]]=="Every day"){
      temp<-1
    }
    }
  }
  sc<-sc+temp
  temp<-0
  for(j in 1:length(ind3)){
    if(is.na(phen_all_use[i,ind3[j]])==F){
    if(phen_all_use[i,ind3[j]]=="Every day"){
      temp<-1
    }
    }
  }
  sc<-sc+temp
  temp<-0
  for(j in 1:length(ind4)){
    if(is.na(phen_all_use[i,ind4[j]])==F){
    if(phen_all_use[i,ind4[j]]=="Every day"){
      temp<-1
    }
    }
  }
  sc<-sc+temp
  temp<-0
  for(j in 1:length(ind5)){
    if(is.na(phen_all_use[i,ind5[j]])==F){
    if(phen_all_use[i,ind5[j]]=="Every day"){
      temp<-1
    }
    }
  }
  sc<-sc+temp
  temp<-0
  for(j in 1:length(ind6)){
    if(is.na(phen_all_use[i,ind6[j]])==F){
    if(phen_all_use[i,ind6[j]]=="Every day"){
      temp<-1
    }
    }
  }
  sc<-sc+temp
  temp<-0
  for(j in 1:length(ind7)){
    if(is.na(phen_all_use[i,ind7[j]])==F){
    if(phen_all_use[i,ind7[j]]=="Every day"){
      temp<-1
    }
    }
  }
  sc<-sc+temp
  temp<-0
  for(j in 1:length(ind8)){
    if(is.na(phen_all_use[i,ind8[j]])==F){
    if(phen_all_use[i,ind8[j]]=="Every day"){
      temp<-1
    }
    }
  }
  sc<-sc+temp
  temp<-0
  for(j in 1:length(ind9)){
    if(is.na(phen_all_use[i,ind9[j]])==F){
    if(phen_all_use[i,ind9[j]]=="Every day"){
      temp<-1
    }
    }
  }
  sc<-sc+temp
  temp<-0
  for(j in 1:length(ind10)){
    if(is.na(phen_all_use[i,ind10[j]])==F){
    if(phen_all_use[i,ind10[j]]=="Every day"){
      temp<-1
    }
    }
  }
  sc<-sc+temp
  temp<-0
  for(j in 1:length(ind11)){
    if(is.na(phen_all_use[i,ind11[j]])==F){
    if(phen_all_use[i,ind11[j]]=="Every day"){
      temp<-1
    }
    }
  }
  sc<-sc+temp
  temp<-0
  dds[i,1]<-sc
}


#phen_all_use<-cbind(phen_all_use,"DDS"=dds)
#write.csv(phen_all_use,'phen_all_use6_all.csv')








