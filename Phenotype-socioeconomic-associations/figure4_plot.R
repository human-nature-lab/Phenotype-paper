
# ==================================
# By: Shivkumar Vishnempet Shridhar, Christakis and Brito group, Yale (2023)
# Honduras microbiome project, plotting for
# socioeconomic-phenotype association
# ==================================

#
library(ggplot2)
library(ggpubr)
library(pheatmap)
#

effect_size_phen_lmer<-read.csv('species_phenotype_assc.csv',row.names=1)
fdr_phen_lmer<-read.csv('species_phenotype_fdr.csv',row.names=1)


#Choosing only food & animal phenotypes
effect_size_phen_lmer<-effect_size_phen_lmer[,83:123]
fdr_phen_lmer<-fdr_phen_lmer[,83:123]

#Creating concise species names

mb_samp_sp<-rownames(effect_size_phen_lmer)
mb_samp_sp<-as.data.frame(mb_samp_sp)

mb_samp_sp_2<-array(0,dim=c(dim(mb_samp_sp)[1],1))
#lk<-array("#000000",dim=c(dim(mb_samp_sp)[1],1))
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
      #lk[i,1]<-"#feb24c"
      
    }
    if((ind_temp6>0)&(ind_temp7<0)){
      ind_tempb<-unlist(gregexpr("c__",as.character(mb_samp_sp[i,1]),fixed=T))
      mb_samp_sp_2[i,1]<-paste0("{",substr(as.character(mb_samp_sp[i,1]),ind_tempb,ind_temp6-2),"}",substr(as.character(mb_samp_sp[i,1]),ind_temp2+3,ind_temp3-2)," (",substr(as.character(mb_samp_sp[i,1]),ind_temp3+3,nchar(as.character(mb_samp_sp[i,1]))),")")
      ind_t<-rbind(ind_t,i)
      #lk[i,1]<-"#feb24c"
    }
    if((ind_temp5>0)&((ind_temp6+ind_temp7)<0)){
      ind_tempb<-unlist(gregexpr("o__",as.character(mb_samp_sp[i,1]),fixed=T))
      mb_samp_sp_2[i,1]<-paste0("{",substr(as.character(mb_samp_sp[i,1]),ind_tempb,ind_temp5-2),"}",substr(as.character(mb_samp_sp[i,1]),ind_temp2+3,ind_temp3-2)," (",substr(as.character(mb_samp_sp[i,1]),ind_temp3+3,nchar(as.character(mb_samp_sp[i,1]))),")")
      ind_t<-rbind(ind_t,i)
      #lk[i,1]<-"#feb24c"
    }
    if((ind_temp4>0)&((ind_temp6+ind_temp7+ind_temp5)<0)){
      ind_tempb<-unlist(gregexpr("f__",as.character(mb_samp_sp[i,1]),fixed=T))
      mb_samp_sp_2[i,1]<-paste0("{",substr(as.character(mb_samp_sp[i,1]),ind_tempb,ind_temp4-2),"}",substr(as.character(mb_samp_sp[i,1]),ind_temp2+3,ind_temp3-2)," (",substr(as.character(mb_samp_sp[i,1]),ind_temp3+3,nchar(as.character(mb_samp_sp[i,1]))),")")
      ind_t<-rbind(ind_t,i)
      #lk[i,1]<-"#feb24c"
    }
    if((ind_temp6+ind_temp7+ind_temp5+ind_temp4)<0){
      mb_samp_sp_2[i,1]<-paste0(substr(as.character(mb_samp_sp[i,1]),ind_temp2+3,ind_temp3-2)," (",substr(as.character(mb_samp_sp[i,1]),ind_temp3+3,nchar(as.character(mb_samp_sp[i,1]))),")")
      ind_t<-rbind(ind_t,i)
    }
  }
}


##Cleaning labels

colnames(effect_size_phen_lmer)[1:17]<-c("Number of partners","Living with partner (N=1232)","Grades 1-3 (N=920)","Grades 4-6 (N=528)","Grades >6 (N=129)","Friend ties (same building)","Friend ties (different building)","Betweenness (friendship)","Transitivity (friendship)","Familial ties (same building)","Familial ties (different building)","Betweenness (familial)","Transitivity (familial)","Degree (all ties)","Clustering coefficient (all ties)","Betweenness (all ties)","Kin percentage (to third degree)")
colnames(effect_size_phen_lmer)[19:23]<-c("Risk taking","Washing hands (N=1603)","Distance to village center","Distance to main road","Number of churches")
colnames(effect_size_phen_lmer)[26:41]<-c("Monthly expenditure","Household size","Household wealth index","TV (N=958)","No electronics (N=89)","Eath/Sand floor (N=594)","Ceramic floor (N=151)","Glass windows (N=86)","Unfinished windows (N=59)","Clay/mud walls (N=1256)","Cement walls (N=402)","Concrete roof (N=37)","Sleeping rooms","Spring (protected) (N=1441)","Tube well (N=91)","Dug well (protected) (N=93)")

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


ind_phen<-1:41

fdr_chk_lmer_all<-rowSums(fdr_chk_lmer[,ind_phen])#c(1,4:11,16:17,19:20,22:25)#44:83(personality types discrete)#
length(which(fdr_chk_lmer_all>11))#>5


ind_c<-which(fdr_chk_lmer_all>11)#was ind_c1

effect_size_phen_sig<-effect_size_phen_lmer[ind_c,ind_phen]
fdr_phen_sig<-fdr_chk_lmer[ind_c,ind_phen]

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
    if(fdr_phen_sig[i,j]==1){
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
colnames(effect_size_phen_sig_plot)<-as.character(mb_samp_sp_2[ind_c,1])
paletteLength<-15
myColor <- colorRampPalette(c("#cb181d","white", "#33a02c"))(paletteLength)#myColor <- colorRampPalette(c("Orange", "white", "Darkblue"))(paletteLength)
myBreaks <- c(seq(min(effect_size_phen_sig_plot),0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(effect_size_phen_sig_plot)/paletteLength, max(effect_size_phen_sig_plot), length.out=floor(paletteLength/2)))

myBreaks
myBreaks[7]<-as.numeric(-0.001)
myBreaks[9]<-as.numeric(0.001)
myBreaks

library(pheatmap)

kl<-read.csv('sp_s_sel.csv',row.names=1)#Ordering species for better visualization

effect_size_phen_sig_plot<-effect_size_phen_sig_plot[,kl[,1]]
disp_fdr<-disp_fdr[kl[,1],]

#pheatmap(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,cluster_rows = F,fontsize_number = 13,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr),number_color = "grey30")
pheatmap(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = F)

ol<-c(10^(-6),10^(-3),10^(0.9-(0.5*5)),10^(0.9-(0.5*4)),10^(0.9-(0.5*3)),10^(0.9-(0.5*2)),10^(0.9-(0.5*1)),10^0.9)
myBreaks6<-c(-ol[length(ol):(-1):1],ol)


disp_fdr5<-t(disp_fdr)


paletteLength<-15
myColor6 <- colorRampPalette(c("#cb181d","white", "#33a02c"))(paletteLength)

effect_size_phen_sig_plot3<-effect_size_phen_sig_plot


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
      if(effect_size_phen_sig_plot3[i,j]>myBreaks6[length(myColor6)]){
        temp<-c(temp,myColor6[length(myColor6)])
      }
      for(k in 2:length(myColor6)){
        if((effect_size_phen_sig_plot3[i,j]>=myBreaks6[k-1])&(effect_size_phen_sig_plot3[i,j]<=myBreaks6[k])){
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
      if(effect_size_phen_sig_plot3[i,j]>myBreaks6[length(myColor6)]){
        temp<-c(temp,myColor6[length(myColor6)])
      }
      for(k in 2:length(myColor6)){
        if((effect_size_phen_sig_plot3[i,j]>=myBreaks6[k-1])&(effect_size_phen_sig_plot3[i,j]<=myBreaks6[k])){
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

out2<-pheatmap(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = F)


d<-as.data.frame(d)

colnames(d)<-c("x1","y1","shp","clr","x11","y11")

library(ggplot2)
d$x1<-factor(d$x1,levels=c(dim(effect_size_phen_sig_plot3)[1]:-1:1))
d$y1<-factor(d$y1,levels=1:dim(effect_size_phen_sig_plot3)[2])
d$x11<-factor(d$x11,levels=rownames(effect_size_phen_sig_plot3)[c(length(rownames(effect_size_phen_sig_plot3)):-1:1)])
#d$y11<-factor(d$y11,levels=colnames(effect_size_phen_sig_plot2))
d$y11<-factor(d$y11,levels=colnames(effect_size_phen_sig_plot3))

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

#"#377eb8","#e5d8bd","#33a02c","#6a3d9a","#a50026","#dd3497","#dfc27d","#ff7f00","#000000","#525252","#a6bddb","#cab2d6"

pdf('figure_4A.pdf',width=13,height=10)
ggplot(d2,aes(y11,x11,fill=clr))+
  geom_point(shape=20,color=clk,stroke=0,size=5.9)+#Extra line to give it outline
  geom_point(shape=20,color=d2$clr,stroke=0,size=5.5)+
  theme(axis.text.x=element_text(angle=90,hjust=0.95,color="black",face="bold",vjust=0.2),legend.position = "none",axis.text.y=element_text(color="black",face="bold"),panel.background = element_blank())+#plot.margin = margin(c(-20,(dim(effect_size_phen_sig_plot3)[2]+4),-40,(dim(effect_size_phen_sig_plot3)[1]+10)))
  xlab("")+ylab("")+expand_limits(x=c(-0.5,(dim(effect_size_phen_sig_plot3)[2]+10)))+#scale_x_discrete(1:(dim(effect_size_phen_sig_plot3)[2]+20))+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot3)[2]+1,y=dim(effect_size_phen_sig_plot3)[1]+0.5,yend=dim(effect_size_phen_sig_plot3)[1]+0.5,col="#4575b4",size=0.3)+#Partner
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot3)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot3)[2]+1),ymin=(dim(effect_size_phen_sig_plot3)[1]-1.5),ymax=(dim(effect_size_phen_sig_plot3)[1]+0.5),color="#4575b4",fill="#4575b4")+#+
  theme(plot.margin=margin(c(1,12,1,1), unit = "pt"))+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+3.7),y=(dim(effect_size_phen_sig_plot3)[1]-0.5),label="Co-habiting",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+3.5),y=(dim(effect_size_phen_sig_plot3)[1]-1.7),label="partners",color="black",fontface="bold",hjust=0.6)+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot3)[1]+0.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-1.5),color="#4575b4",fill="#4575b4")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot3)[2]+1,y=dim(effect_size_phen_sig_plot3)[1]-1.5,yend=dim(effect_size_phen_sig_plot3)[1]-1.5,col="#525252",size=0.3)+#Education
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot3)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot3)[2]+1),ymin=(dim(effect_size_phen_sig_plot3)[1]-4.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-1.5),color="#525252",fill="#525252")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+3.4),y=(dim(effect_size_phen_sig_plot3)[1]-4),label="Education",color="black",fontface="bold")+
  #annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+4),y=(dim(effect_size_phen_sig_plot3)[1]-9.7),label="animals",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot3)[1]-4.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-1.5),color="#525252",fill="#525252")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot3)[2]+1,y=dim(effect_size_phen_sig_plot3)[1]-4.5,yend=dim(effect_size_phen_sig_plot3)[1]-4.5,col="#d94801",size=0.3)+#Friendship
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot3)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot3)[2]+1),ymin=(dim(effect_size_phen_sig_plot3)[1]-8.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-4.5),color="#d94801",fill="#d94801")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+3.5),y=(dim(effect_size_phen_sig_plot3)[1]-7),label="Friendship",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+3.5),y=(dim(effect_size_phen_sig_plot3)[1]-8.2),label="(non-kin)",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot3)[1]-8.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-4.5),color="#d94801",fill="#d94801")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot3)[2]+1,y=dim(effect_size_phen_sig_plot3)[1]-8.5,yend=dim(effect_size_phen_sig_plot3)[1]-8.5,col="#e7298a",size=0.3)+#Family
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot3)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot3)[2]+1),ymin=(dim(effect_size_phen_sig_plot3)[1]-12.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-8.5),color="#e7298a",fill="#e7298a")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+3.7),y=(dim(effect_size_phen_sig_plot3)[1]-10.5),label="Family (kin)",color="black",fontface="bold")+
  #annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+4),y=(dim(effect_size_phen_sig_plot3)[1]-19.7),label="animals",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot3)[1]-12.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-8.5),color="#e7298a",fill="#e7298a")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot3)[2]+1,y=dim(effect_size_phen_sig_plot3)[1]-12.5,yend=dim(effect_size_phen_sig_plot3)[1]-12.5,col="#807dba",size=0.3)+#All
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot3)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot3)[2]+1),ymin=(dim(effect_size_phen_sig_plot3)[1]-16.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-12.5),color="#807dba",fill="#807dba")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+3),y=(dim(effect_size_phen_sig_plot3)[1]-13),label="All ties",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+4.8),y=(dim(effect_size_phen_sig_plot3)[1]-14.2),label="(kin + non-kin)",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot3)[1]-16.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-12.5),color="#807dba",fill="#807dba")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot3)[2]+1,y=dim(effect_size_phen_sig_plot3)[1]-16.5,yend=dim(effect_size_phen_sig_plot3)[1]-16.5,col="#fcadad",size=0.3)+#Risky behavior
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot3)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot3)[2]+1),ymin=(dim(effect_size_phen_sig_plot3)[1]-19.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-16.5),color="#fcadad",fill="#fcadad")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+2.7),y=(dim(effect_size_phen_sig_plot3)[1]-17.5),label="Risky",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+3.5),y=(dim(effect_size_phen_sig_plot3)[1]-18.7),label="behavior",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot3)[1]-19.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-16.5),color="#fcadad",fill="#fcadad")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot3)[2]+1,y=dim(effect_size_phen_sig_plot3)[1]-19.5,yend=dim(effect_size_phen_sig_plot3)[1]-19.5,col="#00b0f0",size=0.3)+#Village factors
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot3)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot3)[2]+1),ymin=(dim(effect_size_phen_sig_plot3)[1]-23.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-19.5),color="#00b0f0",fill="#00b0f0")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+3),y=(dim(effect_size_phen_sig_plot3)[1]-21),label="Village",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+3.3),y=(dim(effect_size_phen_sig_plot3)[1]-22.2),label="factors",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot3)[1]-23.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-19.5),color="#00b0f0",fill="#00b0f0")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot3)[2]+1,y=dim(effect_size_phen_sig_plot3)[1]-23.5,yend=dim(effect_size_phen_sig_plot3)[1]-23.5,col="#7030a0",size=0.3)+#Income
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot3)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot3)[2]+1),ymin=(dim(effect_size_phen_sig_plot3)[1]-28.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-23.5),color="#7030a0",fill="#7030a0")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+3),y=(dim(effect_size_phen_sig_plot3)[1]-24.5),label="Income",color="black",fontface="bold")+
  #annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+3.8),y=(dim(effect_size_phen_sig_plot3)[1]-24.2),label="factors",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot3)[1]-25.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-23.5),color="#7030a0",fill="#7030a0")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot3)[2]+1,y=dim(effect_size_phen_sig_plot3)[1]-25.5,yend=dim(effect_size_phen_sig_plot3)[1]-25.5,col="#ffc000",size=0.3)+#HH essentials
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot3)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot3)[2]+1),ymin=(dim(effect_size_phen_sig_plot3)[1]-40.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-25.5),color="#ffc000",fill="#ffc000")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+3.7),y=(dim(effect_size_phen_sig_plot3)[1]-30.5),label="Household",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+4.1),y=(dim(effect_size_phen_sig_plot3)[1]-31.7),label="essentials",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot3)[1]-37.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-25.5),color="#ffc000",fill="#ffc000")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot3)[2]+1,y=dim(effect_size_phen_sig_plot3)[1]-37.5,yend=dim(effect_size_phen_sig_plot3)[1]-37.5,col="#000000",size=0.3)+#Water sources
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot3)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot3)[2]+1),ymin=(dim(effect_size_phen_sig_plot3)[1]-41.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-37.5),color="#000000",fill="#000000")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+2.9),y=(dim(effect_size_phen_sig_plot3)[1]-38.5),label="Water",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+3.5),y=(dim(effect_size_phen_sig_plot3)[1]-39.7),label="sources",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot3)[1]-40.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-37.5),color="#000000",fill="#000000")

dev.off()

##Figure 4B


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

pdf('figure_4B.pdf',width=4,height=6)
ggplot(div_alpha_ph2, aes(x=wealth, y=as.numeric(as.character(alpha_div)),group=wealth)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#b2df8a","#33a02c","#e31a1c","#fdbf6f"),outlier.alpha=0)+labs(x="",y="Shannon diversity")+
  theme(text = element_text(size=20),axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black",angle=60,hjust=1),panel.background=element_blank(),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))+ scale_x_discrete(breaks=seq(1,5,1))+
  stat_pvalue_manual(data=p_ch,label = "pval", y.position = c(5.2,5.5,5.8,6.1),size=5)+ylim(1.8,6.7)+#+
  theme(plot.margin=margin(c(0.5,0.5,0.5,1), unit = "cm"))+scale_x_discrete(labels=c("Wealth index = 1","Wealth index = 2","Wealth index = 3","Wealth index = 4","Wealth index = 5"))
dev.off()



