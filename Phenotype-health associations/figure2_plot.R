# ==================================
# By: Shivkumar Vishnempet Shridhar, Christakis and Brito group, HNL, Yale (2023)
# Honduras microbiome project, plotting for
# health-phenotype association
# ==================================

#
library(ggplot2)
library(ggpubr)
library(pheatmap)
#

effect_size_phen_lmer<-read.csv('species_phenotype_assc.csv',row.names=1)
fdr_phen_lmer<-read.csv('species_phenotype_fdr.csv',row.names=1)


#Choosing only health phenotypes
effect_size_phen_lmer<-effect_size_phen_lmer[,1:41]
fdr_phen_lmer<-fdr_phen_lmer[,1:41]


#Creating concise species names

mb_samp_sp<-rownames(effect_size_phen_lmer)
mb_samp_sp<-as.data.frame(mb_samp_sp)

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





##Cleaning labels

colnames(effect_size_phen_lmer)[31:32]<-c("Dementia (N=100)","Cognitive impairment (N=185)")
colnames(effect_size_phen_lmer)[36:38]<-c("Mild anxiety (N=395)","Moderate anxiety (N=123)","Severe anxiety (N=59)")
colnames(effect_size_phen_lmer)[39:41]<-c("Mild depression (N=469)","Moderate depression (N=138)","Severe depression (N=89)")
colnames(effect_size_phen_lmer)[21:27]<-c("Painkillers (N=1184)","Antibiotics (N=242)","Anti-diarrheal (N=60)","Anti-parasitic (N=49)","Anti-fungal (N=41)","Vitamins (N=335)","Anti-hypertensive (N=82)")
colnames(effect_size_phen_lmer)[1:13]<-c("Hemoglobin A1c","Mean arterial presure","Body Mass Index","Heart rate","Oxygen saturation","Hemoglobin total","Poor (N=223)","Fair (N=969)","Very good (N=79)","Excellent (N=132)","Cough (N=413)","Diarrhea (N=85)","Bristol stool scale")
colnames(effect_size_phen_lmer)[14:20]<-c("Diabetes (N=40)","Allergies (N=124)","Heart disease (N=59)","Asthma (N=62)","Stomach illness (N=261)","Intestinall illness (N=103)","Arthritis (N=43)")
colnames(effect_size_phen_lmer)[28:30]<-c("Reserved (N=1666)","Nervous (N=950)","Openness (N=1550)")
colnames(effect_size_phen_lmer)[33:35]<-c("Alcohol daily frequency","Cigarette usage (N=82)","Cigarette frequency")


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
length(which(fdr_chk_lmer_all>7))#>5


ind_c<-which(fdr_chk_lmer_all>7)#was ind_c1

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


kl<-read.csv('sp_h_sel.csv',row.names=1)#Ordering species for better visualization

effect_size_phen_sig_plot<-effect_size_phen_sig_plot[,kl[,1]]
disp_fdr<-disp_fdr[kl[,1],]

#pheatmap(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,cluster_rows = F,fontsize_number = 13,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr),number_color = "grey30")
pheatmap(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = F)





for(i in 1:dim(disp_fdr)[1]){
  if(i==1){
    sp_c<-length(which(disp_fdr[i,]!=""))
  }else{
    sp_c<-rbind(sp_c,length(which(disp_fdr[i,]!="")))
  }
}
colnames(effect_size_phen_sig_plot)[which(sp_c==max(sp_c))]
length(which(grepl("{",colnames(effect_size_phen_sig_plot),fixed=T)==T))

length(which(grepl("{",mb_samp_sp_2,fixed=T)==T))

min(effect_size_phen_sig_plot[which(effect_size_phen_sig_plot>0)])

effect_size_phen_sig_plot_h<-effect_size_phen_sig_plot

max(effect_size_phen_sig_plot_h)

#ol<-c(10^(-6),10^(-3),10^(0.72-(0.4*5)),10^(0.72-(0.4*4)),10^(0.72-(0.4*3)),10^(0.72-(0.4*2)),10^(0.72-(0.4*1)),10^0.72)
#Need to change because of new highs around 8
ol<-c(10^(-6),10^(-3),10^(0.9-(0.5*5)),10^(0.9-(0.5*4)),10^(0.9-(0.5*3)),10^(0.9-(0.5*2)),10^(0.9-(0.5*1)),10^0.9)
myBreaks11<-c(-ol[length(ol):(-1):1],ol)

length(which(fdr_phen_lmer[,ind_phen]<0.05))

##### Bubble plot

disp_fdr3<-t(disp_fdr)
myColor#Color
myBreaks#breaks
myBreaks4<-myBreaks11
myColor4<-myColor


for(i in 1:dim(effect_size_phen_sig_plot)[1]){
  for(j in 1:dim(effect_size_phen_sig_plot)[2]){
    if((i==1)&(j==1)){
      if(disp_fdr3[i,j]==""){
        temp<-c(i,j,NA)
      }else if(disp_fdr3[i,j]=="-"){
        temp<-c(i,j,"-")
      }else if(disp_fdr3[i,j]=="+"){
        temp<-c(i,j,3)
      }
      if(effect_size_phen_sig_plot[i,j]>myBreaks4[length(myColor4)]){
        temp<-c(temp,myColor4[length(myColor4)])
      }
      for(k in 2:length(myColor4)){
        if((effect_size_phen_sig_plot[i,j]>=myBreaks4[k-1])&(effect_size_phen_sig_plot[i,j]<=myBreaks4[k])){
          temp<-c(temp,myColor4[k-1])
        }
      }
      temp<-c(temp,rownames(effect_size_phen_sig_plot)[i],colnames(effect_size_phen_sig_plot)[j])
      d<-temp
    }else{
      if(disp_fdr3[i,j]==""){
        temp<-c(i,j,NA)
      }else if(disp_fdr3[i,j]=="-"){
        temp<-c(i,j,45)
      }else if(disp_fdr3[i,j]=="+"){
        temp<-c(i,j,43)
      }
      #if(effect_size_phen_sig_plot[i,j]<myBreaks4[1]){
      #  temp<-c(temp,myColor4[1])
      #}else 
      if(effect_size_phen_sig_plot[i,j]>myBreaks4[length(myColor4)]){
        temp<-c(temp,myColor4[length(myColor4)])
      }
      for(k in 2:length(myColor4)){
        if((effect_size_phen_sig_plot[i,j]>=myBreaks4[k-1])&(effect_size_phen_sig_plot[i,j]<=myBreaks4[k])){
          temp<-c(temp,myColor4[k-1])
          #tpp<-rbind(tpp,k)
        }
      }
      temp<-c(temp,rownames(effect_size_phen_sig_plot)[i],colnames(effect_size_phen_sig_plot)[j])
      d<-rbind(d,temp)
      #sz_t<-rbind(sz_t,length(temp))
    }
    
  }
}

##Change order of species --hclust
out<-pheatmap(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks11,display_numbers = t(disp_fdr),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

d<-as.data.frame(d)

colnames(d)<-c("x1","y1","shp","clr","x11","y11")


d$x1<-factor(d$x1,levels=c(dim(effect_size_phen_sig_plot)[1]:-1:1))
d$y1<-factor(d$y1,levels=1:dim(effect_size_phen_sig_plot)[2])
d$x11<-factor(d$x11,levels=c(rownames(effect_size_phen_sig_plot)[c(length(rownames(effect_size_phen_sig_plot)):-1:1)]))
d$y11<-factor(d$y11,levels=colnames(effect_size_phen_sig_plot))

d2<-d
for(i in 1:dim(d)[1]){
  if(is.na(d$shp[i])==T){
    d2$clr[i]<-"#FFFFFF"
  }
}

clk<-d2$clr
clk<-ifelse(d2$clr=="#FFFFFF","#FFFFFF","#000000")

#"#377eb8","#e5d8bd","#33a02c","#6a3d9a","#a50026","#dd3497","#dfc27d","#ff7f00","#000000","#525252","#a6bddb","#cab2d6"


#Saving plot as PDF

pdf('figure_2A.pdf',width=13,height=10)
ggplot(d2,aes(y11,x11,fill=clr))+
  geom_point(shape=20,color=clk,stroke=0,size=5.9)+#Extra line to give it outline
  geom_point(shape=20,color=d2$clr,stroke=0,size=5.5)+
  #geom_point(shape=as.numeric(as.character(d$shp)),fill="black",size=5)+
  theme(axis.text.x=element_text(angle=90,hjust=0.95,color="black",face="bold",vjust=0.2),legend.position = "none",axis.text.y=element_text(color="black",face="bold"),panel.background = element_blank())+#plot.margin = margin(c(-20,(dim(effect_size_phen_sig_plot)[2]+4),-40,(dim(effect_size_phen_sig_plot)[1]+10)))
  xlab("")+ylab("")+expand_limits(x=c(-0.5,(dim(effect_size_phen_sig_plot)[2]+10)))+#scale_x_discrete(1:(dim(effect_size_phen_sig_plot)[2]+20))+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot)[2]+1,y=dim(effect_size_phen_sig_plot)[1]+0.5,yend=dim(effect_size_phen_sig_plot)[1]+0.5,col="#4575b4",size=0.2)+
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot)[2]+1),ymin=(dim(effect_size_phen_sig_plot)[1]-5.5),ymax=(dim(effect_size_phen_sig_plot)[1]+0.5),color="#4575b4",fill="#4575b4")+#+
  theme(plot.margin=margin(c(1,10,1,1), unit = "pt"))+
  annotate("text",x=(dim(effect_size_phen_sig_plot)[2]+5.2),y=(dim(effect_size_phen_sig_plot)[1]-2),label="Physiological",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot)[2]+5.2),y=(dim(effect_size_phen_sig_plot)[1]-3.2),label="variables",color="black",hjust=0.7,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot)[1]-5.5),ymax=(dim(effect_size_phen_sig_plot)[1]+0.5),color="#4575b4",fill="#4575b4")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot)[2]+1,y=dim(effect_size_phen_sig_plot)[1]-5.5,yend=dim(effect_size_phen_sig_plot)[1]-5.5,col="#969696",size=0.2)+#Overall health
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot)[2]+1),ymin=(dim(effect_size_phen_sig_plot)[1]-9.5),ymax=(dim(effect_size_phen_sig_plot)[1]-5.5),color="#969696",fill="#969696")+
  annotate("text",x=(dim(effect_size_phen_sig_plot)[2]+5.2),y=(dim(effect_size_phen_sig_plot)[1]-7.5),label="Overall health",color="black",fontface="bold",hjust=0.48)+
  #annotate("text",x=(dim(effect_size_phen_sig_plot)[2]+5),y=(dim(effect_size_phen_sig_plot)[1]-8.2),label="health",color="black",hjust=0.7,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot)[1]-9.5),ymax=(dim(effect_size_phen_sig_plot)[1]-5.5),color="#969696",fill="#969696")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot)[2]+1,y=dim(effect_size_phen_sig_plot)[1]-9.5,yend=dim(effect_size_phen_sig_plot)[1]-9.5,col="#dfc27d",size=0.2)+#Acute conditions
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot)[2]+1),ymin=(dim(effect_size_phen_sig_plot)[1]-12.5),ymax=(dim(effect_size_phen_sig_plot)[1]-9.5),color="#dfc27d",fill="#dfc27d")+
  annotate("text",x=(dim(effect_size_phen_sig_plot)[2]+3.8),y=(dim(effect_size_phen_sig_plot)[1]-10.5),label="Acute",color="black",fontface="bold",hjust=0.75)+
  annotate("text",x=(dim(effect_size_phen_sig_plot)[2]+5.2),y=(dim(effect_size_phen_sig_plot)[1]-11.7),label="conditions",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot)[1]-12.5),ymax=(dim(effect_size_phen_sig_plot)[1]-9.5),color="#dfc27d",fill="#dfc27d")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot)[2]+1,y=dim(effect_size_phen_sig_plot)[1]-12.5,yend=dim(effect_size_phen_sig_plot)[1]-12.5,col="#a6cee3",size=0.2)+#Chronic conditions
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot)[2]+1),ymin=(dim(effect_size_phen_sig_plot)[1]-19.5),ymax=(dim(effect_size_phen_sig_plot)[1]-12.5),color="#a6cee3",fill="#a6cee3")+
  annotate("text",x=(dim(effect_size_phen_sig_plot)[2]+5),y=(dim(effect_size_phen_sig_plot)[1]-15.5),label="Chronic",color="black",fontface="bold",hjust=0.75)+
  annotate("text",x=(dim(effect_size_phen_sig_plot)[2]+5.2),y=(dim(effect_size_phen_sig_plot)[1]-16.7),label="conditions",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot)[1]-19.5),ymax=(dim(effect_size_phen_sig_plot)[1]-12.5),color="#a6cee3",fill="#a6cee3")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot)[2]+1,y=dim(effect_size_phen_sig_plot)[1]-19.5,yend=dim(effect_size_phen_sig_plot)[1]-19.5,col="#762a83",size=0.2)+#Medications
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot)[2]+1),ymin=(dim(effect_size_phen_sig_plot)[1]-26.5),ymax=(dim(effect_size_phen_sig_plot)[1]-19.5),color="#762a83",fill="#762a83")+
  annotate("text",x=(dim(effect_size_phen_sig_plot)[2]+5),y=(dim(effect_size_phen_sig_plot)[1]-23),label="Medications",color="black",fontface="bold")+
  #annotate("text",x=(dim(effect_size_phen_sig_plot)[2]+5.2),y=(dim(effect_size_phen_sig_plot)[1]-22.7),label="conditions",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot)[1]-26.5),ymax=(dim(effect_size_phen_sig_plot)[1]-19.5),color="#762a83",fill="#762a83")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot)[2]+1,y=dim(effect_size_phen_sig_plot)[1]-26.5,yend=dim(effect_size_phen_sig_plot)[1]-26.5,col="#980043",size=0.2)+#Personalities
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot)[2]+1),ymin=(dim(effect_size_phen_sig_plot)[1]-29.5),ymax=(dim(effect_size_phen_sig_plot)[1]-26.5),color="#980043",fill="#980043")+
  annotate("text",x=(dim(effect_size_phen_sig_plot)[2]+5.5),y=(dim(effect_size_phen_sig_plot)[1]-27.5),label="Personalities",color="black",fontface="bold")+
  #annotate("text",x=(dim(effect_size_phen_sig_plot)[2]+5.2),y=(dim(effect_size_phen_sig_plot)[1]-22.7),label="conditions",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot)[1]-29.5),ymax=(dim(effect_size_phen_sig_plot)[1]-26.5),color="#980043",fill="#980043")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot)[2]+1,y=dim(effect_size_phen_sig_plot)[1]-29.5,yend=dim(effect_size_phen_sig_plot)[1]-29.5,col="#e7298a",size=0.2)+#Cognitive
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot)[2]+1),ymin=(dim(effect_size_phen_sig_plot)[1]-31.5),ymax=(dim(effect_size_phen_sig_plot)[1]-29.5),color="#e7298a",fill="#e7298a")+
  annotate("text",x=(dim(effect_size_phen_sig_plot)[2]+4.2),y=(dim(effect_size_phen_sig_plot)[1]-30.5),label="Cognitive",color="black",fontface="bold")+
  #annotate("text",x=(dim(effect_size_phen_sig_plot)[2]+5.2),y=(dim(effect_size_phen_sig_plot)[1]-22.7),label="conditions",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot)[1]-31.5),ymax=(dim(effect_size_phen_sig_plot)[1]-29.5),color="#e7298a",fill="#e7298a")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot)[2]+1,y=dim(effect_size_phen_sig_plot)[1]-31.5,yend=dim(effect_size_phen_sig_plot)[1]-31.5,col="#adba07",size=0.2)+#Unfavorable habits
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot)[2]+1),ymin=(dim(effect_size_phen_sig_plot)[1]-34.5),ymax=(dim(effect_size_phen_sig_plot)[1]-31.5),color="#adba07",fill="#adba07")+
  annotate("text",x=(dim(effect_size_phen_sig_plot)[2]+5.2),y=(dim(effect_size_phen_sig_plot)[1]-32.5),label="Unfavorable",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot)[2]+4.5),y=(dim(effect_size_phen_sig_plot)[1]-33.7),label="habits",color="black",fontface="bold",hjust=0.6)+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot)[1]-34.5),ymax=(dim(effect_size_phen_sig_plot)[1]-31.5),color="#adba07",fill="#adba07")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot)[2]+1,y=dim(effect_size_phen_sig_plot)[1]-34.5,yend=dim(effect_size_phen_sig_plot)[1]-34.5,col="#d94801",size=0.2)+#Anxiety
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot)[2]+1),ymin=(dim(effect_size_phen_sig_plot)[1]-37.5),ymax=(dim(effect_size_phen_sig_plot)[1]-34.5),color="#d94801",fill="#d94801")+
  annotate("text",x=(dim(effect_size_phen_sig_plot)[2]+5.2),y=(dim(effect_size_phen_sig_plot)[1]-35.5),label="Anxiety",color="black",fontface="bold",hjust=0.85)+
  annotate("text",x=(dim(effect_size_phen_sig_plot)[2]+4.5),y=(dim(effect_size_phen_sig_plot)[1]-36.7),label="(GAD7)",color="black",fontface="bold",hjust=0.6)+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot)[1]-37.5),ymax=(dim(effect_size_phen_sig_plot)[1]-34.5),color="#d94801",fill="#d94801")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot)[2]+1,y=dim(effect_size_phen_sig_plot)[1]-37.5,yend=dim(effect_size_phen_sig_plot)[1]-37.5,col="#000000",size=0.2)+#Depression
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot)[2]+1),ymin=(dim(effect_size_phen_sig_plot)[1]-40.5),ymax=(dim(effect_size_phen_sig_plot)[1]-37.5),color="#000000",fill="#000000")+
  annotate("text",x=(dim(effect_size_phen_sig_plot)[2]+5.2),y=(dim(effect_size_phen_sig_plot)[1]-38.5),label="Depression",color="black",fontface="bold",hjust=0.5)+
  annotate("text",x=(dim(effect_size_phen_sig_plot)[2]+4.5),y=(dim(effect_size_phen_sig_plot)[1]-39.7),label="(PHQ9)",color="black",fontface="bold",hjust=0.6)+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot)[1]-40.5),ymax=(dim(effect_size_phen_sig_plot)[1]-37.5),color="#000000",fill="#000000")

dev.off()

### Figure 2B

div_alpha_ph2<-read.csv('div_alpha_ph2_ch.csv',row.names=1)

div_alpha_ph2$chronic <- factor(div_alpha_ph2$chronic, levels=c("Healthy","Diabetes","Allergies","Heart disease","Asthma","Stomach illness","Intestinal illness","Arthritis"),ordered=TRUE)


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



###Figure 2D


tt1<-read.csv('results_maaslin_healthy_vs_unhealthy.csv',row.names=1)
mbl_ind<-match(as.character(tt1[,1]),mb_samp_sp[,1])
tt2<-cbind(tt1,"sp_name"=mb_samp_sp_2[mbl_ind,1])


efz<-tt2$effect_size[which(tt2$FDR<0.05)]
efz_names<-as.character(tt2$sp_name[which(tt2$FDR<0.05)])
efz2<-efz[order(efz)]
efz_names2<-efz_names[order(efz)]

fdr<-tt2$FDR[which(tt2$FDR<0.05)]
fdr2<-fdr[order(efz)]


health_bar<-cbind(cbind(cbind(efz2[c(1:10,(length(efz)-9):length(efz))],fdr2[c(1:10,(length(efz)-9):length(efz))]),as.character(c(efz_names2[1:10],efz_names2[(length(efz_names2)-9):length(efz_names2)]))),c(array("Disease associated",dim=c(10,1)),array("Health associated",dim=c(10,1))))


colnames(health_bar)<-c("effect","fdr","name","category")
health_bar<-as.data.frame(health_bar)

stp<-c("***","***","**","**","***","**","**","**","**","**","*","***","***","**","**","**","***","***","***","***")

pdf('figure_2D.pdf',width=13,height=6)
ggplot(health_bar,aes(x=reorder(as.character(name),as.numeric(as.character(effect))),y=as.numeric(as.character(effect)),fill=category))+
  geom_bar(stat="identity",color="black",size=1)+scale_fill_manual(values = c("#4575b4","#d73027"))+#+stat_bin(geom="text",aes(),vjust=-1)
  theme(text = element_text(size=20),panel.grid.major = element_blank(),panel.background = element_blank(),axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black"),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),legend.position = "top",legend.text = element_text(size=20)) +
  ylab("Effect size")+coord_flip()+xlab("")+#+mdthemes::md_theme_classic()
  theme(plot.margin=margin(c(0,2,0,0), unit = "cm"))
dev.off()


