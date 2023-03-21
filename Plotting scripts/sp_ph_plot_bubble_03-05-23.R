##########################################
#Run this after running Normal health plotting in dmp_lmer_plotting_all

disp_fdr3<-t(disp_fdr2)
myColor#Color
myBreaks#breaks
#myBreaks3<-myBreaks
#myBreaks3[8]<-0.00036*(-1)
myBreaks4<-myBreaks11
myColor4<-myColor
rownames(effect_size_phen_sig_plot2)[39:41]<-c(" Mild"," Moderate"," Severe")

#myBreaks4<-c(myBreaks3[1:6],-0.36,-0.001,myBreaks3[9:length(myBreaks3)])
#myColor4<-c(myColor[1:6],"#F7DDDE","#F7DDDE",myColor[8:length(myColor)])#was "#F0BCBE","#F7DDDE"

for(i in 1:dim(effect_size_phen_sig_plot2)[1]){
  for(j in 1:dim(effect_size_phen_sig_plot2)[2]){
    if((i==1)&(j==1)){
      if(disp_fdr3[i,j]==""){
        temp<-c(i,j,NA)
      }else if(disp_fdr3[i,j]=="-"){
        temp<-c(i,j,"-")
      }else if(disp_fdr3[i,j]=="+"){
        temp<-c(i,j,3)
      }
      if(effect_size_phen_sig_plot2[i,j]<myBreaks4[1]){
        temp<-c(temp,myColor4[1])
      }else if(effect_size_phen_sig_plot2[i,j]>myBreaks4[length(myColor4)]){
        temp<-c(temp,myColor4[length(myColor4)])
      }
      for(k in 2:length(myColor4)){
        if((effect_size_phen_sig_plot2[i,j]>=myBreaks4[k-1])&(effect_size_phen_sig_plot2[i,j]<=myBreaks4[k])){
          temp<-c(temp,myColor4[k])
        }
      }
      temp<-c(temp,rownames(effect_size_phen_sig_plot2)[i],colnames(effect_size_phen_sig_plot2)[j])
      d<-temp
    }else{
      if(disp_fdr3[i,j]==""){
        temp<-c(i,j,NA)
      }else if(disp_fdr3[i,j]=="-"){
        temp<-c(i,j,45)
      }else if(disp_fdr3[i,j]=="+"){
        temp<-c(i,j,43)
      }
      if(effect_size_phen_sig_plot2[i,j]<myBreaks4[1]){
        temp<-c(temp,myColor4[1])
      }else if(effect_size_phen_sig_plot2[i,j]>myBreaks4[length(myColor4)]){
        temp<-c(temp,myColor4[length(myColor4)])
      }
      for(k in 2:length(myColor4)){
        if((effect_size_phen_sig_plot2[i,j]>=myBreaks4[k-1])&(effect_size_phen_sig_plot2[i,j]<=myBreaks4[k])){
          temp<-c(temp,myColor4[k])
          #tpp<-rbind(tpp,k)
        }
      }
      temp<-c(temp,rownames(effect_size_phen_sig_plot2)[i],colnames(effect_size_phen_sig_plot2)[j])
      d<-rbind(d,temp)
      #sz_t<-rbind(sz_t,length(temp))
    }
    
  }
}

##Change order of species --hclust
out<-pheatmap(effect_size_phen_sig_plot2,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks11,display_numbers = t(disp_fdr2),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)
#colnames(effect_size_phen_sig_plot2[out$tree_col["order"],])


d<-as.data.frame(d)

colnames(d)<-c("x1","y1","shp","clr","x11","y11")

library(ggplot2)
d$x1<-factor(d$x1,levels=c(dim(effect_size_phen_sig_plot2)[1]:-1:1))
d$y1<-factor(d$y1,levels=1:dim(effect_size_phen_sig_plot2)[2])
d$x11<-factor(d$x11,levels=c(rownames(effect_size_phen_sig_plot2)[c(length(rownames(effect_size_phen_sig_plot2)):-1:1)]))
#d$y11<-factor(d$y11,levels=colnames(effect_size_phen_sig_plot2))
d$y11<-factor(d$y11,levels=colnames(effect_size_phen_sig_plot2)[out$tree_col[["order"]]])

#ggplot(d,aes(y1,x1,fill=clr,size=8))+
#  geom_point(shape=21,fill=d$clr)+
#  geom_point(shape=as.numeric(as.character(d$shp)),fill="black",size=5)


#Save 10.0 x 14.0
ggplot(d,aes(y11,x11,fill=clr))+
  geom_point(shape=20,color=d$clr,stroke=0,size=6)+
  #geom_point(shape=as.numeric(as.character(d$shp)),fill="black",size=5)+
  theme(axis.text.x=element_text(angle=90,hjust=0.95,color="black",face="bold",vjust=0.2),legend.position = "none",axis.text.y=element_text(color="black",face="bold"),panel.background = element_blank())+#plot.margin = margin(c(-20,(dim(effect_size_phen_sig_plot2)[2]+4),-40,(dim(effect_size_phen_sig_plot2)[1]+10)))
  xlab("")+ylab("")+expand_limits(x=c(-0.5,(dim(effect_size_phen_sig_plot2)[2]+10)))+#scale_x_discrete(1:(dim(effect_size_phen_sig_plot2)[2]+20))+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot2)[2]+1,y=dim(effect_size_phen_sig_plot2)[1]+0.5,yend=dim(effect_size_phen_sig_plot2)[1]+0.5,col="#4575b4")+
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot2)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot2)[2]+1),ymin=(dim(effect_size_phen_sig_plot2)[1]-5.5),ymax=(dim(effect_size_phen_sig_plot2)[1]+0.5),color="#4575b4",fill="#4575b4")+#+
  theme(plot.margin=margin(c(1,10,1,1), unit = "pt"))+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+5),y=(dim(effect_size_phen_sig_plot2)[1]-2),label="Physiological",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+5),y=(dim(effect_size_phen_sig_plot2)[1]-3.2),label="variables",color="black",hjust=0.7,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot2)[1]-5.5),ymax=(dim(effect_size_phen_sig_plot2)[1]+0.5),color="#4575b4",fill="#4575b4")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot2)[2]+1,y=dim(effect_size_phen_sig_plot2)[1]-5.5,yend=dim(effect_size_phen_sig_plot2)[1]-5.5,col="#d73027",size=0.08)+
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot2)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot2)[2]+1),ymin=(dim(effect_size_phen_sig_plot2)[1]-9.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-5.5),color="#d73027",fill="#d73027")+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+5),y=(dim(effect_size_phen_sig_plot2)[1]-7.5),label="Overall health",color="black",fontface="bold",hjust=0.48)+
  #annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+5),y=(dim(effect_size_phen_sig_plot2)[1]-8.2),label="health",color="black",hjust=0.7,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot2)[1]-9.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-5.5),color="#d73027",fill="#d73027")
  #+xlim(0,(dim(effect_size_phen_sig_plot2)[2]+5))
  #opts(plot.margin = unit(c(-10,dim(effect_size_phen_sig_plot2)[2]+2 , -20, dim(effect_size_phen_sig_plot2)[1]+4), "cm"))


d2<-d
for(i in 1:dim(d)[1]){
  if(is.na(d$shp[i])==T){
    d2$clr[i]<-"#FFFFFF"
  }
}

#"#377eb8","#e5d8bd","#33a02c","#6a3d9a","#a50026","#dd3497","#dfc27d","#ff7f00","#000000","#525252","#a6bddb","#cab2d6"

ggplot(d2,aes(y11,x11,fill=clr))+
  geom_point(shape=20,color=d2$clr,stroke=0,size=6)+
  #geom_point(shape=as.numeric(as.character(d$shp)),fill="black",size=5)+
  theme(axis.text.x=element_text(angle=90,hjust=0.95,color="black",face="bold",vjust=0.2),legend.position = "none",axis.text.y=element_text(color="black",face="bold"),panel.background = element_blank())+#plot.margin = margin(c(-20,(dim(effect_size_phen_sig_plot2)[2]+4),-40,(dim(effect_size_phen_sig_plot2)[1]+10)))
  xlab("")+ylab("")+expand_limits(x=c(-0.5,(dim(effect_size_phen_sig_plot2)[2]+10)))+#scale_x_discrete(1:(dim(effect_size_phen_sig_plot2)[2]+20))+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot2)[2]+1,y=dim(effect_size_phen_sig_plot2)[1]+0.5,yend=dim(effect_size_phen_sig_plot2)[1]+0.5,col="#4575b4",size=0.2)+
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot2)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot2)[2]+1),ymin=(dim(effect_size_phen_sig_plot2)[1]-5.5),ymax=(dim(effect_size_phen_sig_plot2)[1]+0.5),color="#4575b4",fill="#4575b4")+#+
  theme(plot.margin=margin(c(1,10,1,1), unit = "pt"))+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+5.2),y=(dim(effect_size_phen_sig_plot2)[1]-2),label="Physiological",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+5.2),y=(dim(effect_size_phen_sig_plot2)[1]-3.2),label="variables",color="black",hjust=0.7,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot2)[1]-5.5),ymax=(dim(effect_size_phen_sig_plot2)[1]+0.5),color="#4575b4",fill="#4575b4")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot2)[2]+1,y=dim(effect_size_phen_sig_plot2)[1]-5.5,yend=dim(effect_size_phen_sig_plot2)[1]-5.5,col="#969696",size=0.2)+#Overall health
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot2)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot2)[2]+1),ymin=(dim(effect_size_phen_sig_plot2)[1]-9.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-5.5),color="#969696",fill="#969696")+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+5.2),y=(dim(effect_size_phen_sig_plot2)[1]-7.5),label="Overall health",color="black",fontface="bold",hjust=0.48)+
  #annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+5),y=(dim(effect_size_phen_sig_plot2)[1]-8.2),label="health",color="black",hjust=0.7,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot2)[1]-9.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-5.5),color="#969696",fill="#969696")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot2)[2]+1,y=dim(effect_size_phen_sig_plot2)[1]-9.5,yend=dim(effect_size_phen_sig_plot2)[1]-9.5,col="#dfc27d",size=0.2)+#Acute conditions
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot2)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot2)[2]+1),ymin=(dim(effect_size_phen_sig_plot2)[1]-12.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-9.5),color="#dfc27d",fill="#dfc27d")+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+3.8),y=(dim(effect_size_phen_sig_plot2)[1]-10.5),label="Acute",color="black",fontface="bold",hjust=0.75)+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+5.2),y=(dim(effect_size_phen_sig_plot2)[1]-11.7),label="conditions",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot2)[1]-12.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-9.5),color="#dfc27d",fill="#dfc27d")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot2)[2]+1,y=dim(effect_size_phen_sig_plot2)[1]-12.5,yend=dim(effect_size_phen_sig_plot2)[1]-12.5,col="#a6cee3",size=0.2)+#Chronic conditions
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot2)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot2)[2]+1),ymin=(dim(effect_size_phen_sig_plot2)[1]-19.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-12.5),color="#a6cee3",fill="#a6cee3")+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+5),y=(dim(effect_size_phen_sig_plot2)[1]-15.5),label="Chronic",color="black",fontface="bold",hjust=0.75)+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+5.2),y=(dim(effect_size_phen_sig_plot2)[1]-16.7),label="conditions",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot2)[1]-19.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-12.5),color="#a6cee3",fill="#a6cee3")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot2)[2]+1,y=dim(effect_size_phen_sig_plot2)[1]-19.5,yend=dim(effect_size_phen_sig_plot2)[1]-19.5,col="#762a83",size=0.2)+#Medications
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot2)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot2)[2]+1),ymin=(dim(effect_size_phen_sig_plot2)[1]-26.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-19.5),color="#762a83",fill="#762a83")+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+5),y=(dim(effect_size_phen_sig_plot2)[1]-23),label="Medications",color="black",fontface="bold")+
  #annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+5.2),y=(dim(effect_size_phen_sig_plot2)[1]-22.7),label="conditions",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot2)[1]-26.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-19.5),color="#762a83",fill="#762a83")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot2)[2]+1,y=dim(effect_size_phen_sig_plot2)[1]-26.5,yend=dim(effect_size_phen_sig_plot2)[1]-26.5,col="#980043",size=0.2)+#Personalities
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot2)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot2)[2]+1),ymin=(dim(effect_size_phen_sig_plot2)[1]-29.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-26.5),color="#980043",fill="#980043")+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+5.5),y=(dim(effect_size_phen_sig_plot2)[1]-27.5),label="Personalities",color="black",fontface="bold")+
  #annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+5.2),y=(dim(effect_size_phen_sig_plot2)[1]-22.7),label="conditions",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot2)[1]-29.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-26.5),color="#980043",fill="#980043")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot2)[2]+1,y=dim(effect_size_phen_sig_plot2)[1]-29.5,yend=dim(effect_size_phen_sig_plot2)[1]-29.5,col="#e7298a",size=0.2)+#Cognitive
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot2)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot2)[2]+1),ymin=(dim(effect_size_phen_sig_plot2)[1]-31.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-29.5),color="#e7298a",fill="#e7298a")+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+4.2),y=(dim(effect_size_phen_sig_plot2)[1]-30.5),label="Cognitive",color="black",fontface="bold")+
  #annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+5.2),y=(dim(effect_size_phen_sig_plot2)[1]-22.7),label="conditions",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot2)[1]-31.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-29.5),color="#e7298a",fill="#e7298a")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot2)[2]+1,y=dim(effect_size_phen_sig_plot2)[1]-31.5,yend=dim(effect_size_phen_sig_plot2)[1]-31.5,col="#adba07",size=0.2)+#Unfavorable habits
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot2)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot2)[2]+1),ymin=(dim(effect_size_phen_sig_plot2)[1]-34.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-31.5),color="#adba07",fill="#adba07")+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+5.2),y=(dim(effect_size_phen_sig_plot2)[1]-32.5),label="Unfavorable",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+4.5),y=(dim(effect_size_phen_sig_plot2)[1]-33.7),label="habits",color="black",fontface="bold",hjust=0.6)+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot2)[1]-34.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-31.5),color="#adba07",fill="#adba07")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot2)[2]+1,y=dim(effect_size_phen_sig_plot2)[1]-34.5,yend=dim(effect_size_phen_sig_plot2)[1]-34.5,col="#d94801",size=0.2)+#Anxiety
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot2)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot2)[2]+1),ymin=(dim(effect_size_phen_sig_plot2)[1]-37.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-34.5),color="#d94801",fill="#d94801")+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+5.2),y=(dim(effect_size_phen_sig_plot2)[1]-35.5),label="Anxiety",color="black",fontface="bold",hjust=0.85)+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+4.5),y=(dim(effect_size_phen_sig_plot2)[1]-36.7),label="(GAD7)",color="black",fontface="bold",hjust=0.6)+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot2)[1]-37.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-34.5),color="#d94801",fill="#d94801")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot2)[2]+1,y=dim(effect_size_phen_sig_plot2)[1]-37.5,yend=dim(effect_size_phen_sig_plot2)[1]-37.5,col="#000000",size=0.2)+#Depression
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot2)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot2)[2]+1),ymin=(dim(effect_size_phen_sig_plot2)[1]-40.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-37.5),color="#000000",fill="#000000")+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+5.2),y=(dim(effect_size_phen_sig_plot2)[1]-38.5),label="Depression",color="black",fontface="bold",hjust=0.5)+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+4.5),y=(dim(effect_size_phen_sig_plot2)[1]-39.7),label="(PHQ9)",color="black",fontface="bold",hjust=0.6)+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot2)[1]-40.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-37.5),color="#000000",fill="#000000")

#Save 9.5 x 13.5
#Save 10.00 x 13.00 then move the text around -- changed y axis text to 9.5, and the spacing to 1.1, 
#2. then for species, select all text, center it then right it; 3. add zscore legend 4. convert pdf to svg



##Food and animals

##RUN AFTER running food and animals in dmp_lmer
disp_fdr4<-t(disp_fdr[xt,])
myColor#Color
myBreaks11#breaks
#myBreaks5<-myBreaks[c(1:6,8,7,9:length(myBreaks))]#c(myBreaks[1:6],-0.36,-0.001,myBreaks[9:length(myBreaks)])
#myBreaks5<-myBreaks[c((length(myBreaks):-1:9),9:length(myBreaks))]
#myBreaks5[1:8]<-myBreaks5[1:8]*(-1)
myBreaks5<-myBreaks11
myColor5<-myColor
#myColor <- colorRampPalette(c("#cb181d","white", "#33a02c"))(length(myBreaks5))
#myColor5<-myColor[c(1:6,8,7,9:length(myColor))]#c(myColor[1:6],"#F7DDDE","#F7DDDE",myColor[8:length(myColor)])#was "#F0BCBE","#F7DDDE"

effect_size_phen_sig_plot2<-effect_size_phen_sig_plot[,xt]
# geom_segment(x=1,y=1*0.95+coef(op2)[1],xend=-1,yend=-1*0.95+coef(op2)[1],col="red", size=1)+

for(i in 1:dim(effect_size_phen_sig_plot2)[1]){
  for(j in 1:dim(effect_size_phen_sig_plot2)[2]){
    if((i==1)&(j==1)){
      if(disp_fdr4[i,j]==""){
        temp<-c(i,j,NA)
      }else if(disp_fdr4[i,j]=="-"){
        temp<-c(i,j,"-")
      }else if(disp_fdr4[i,j]=="+"){
        temp<-c(i,j,3)
      }
      if(effect_size_phen_sig_plot2[i,j]<myBreaks5[1]){
        temp<-c(temp,myColor5[1])
      }else if(effect_size_phen_sig_plot2[i,j]>myBreaks5[length(myColor5)]){
        temp<-c(temp,myColor5[length(myColor5)])
      }
      for(k in 2:length(myColor5)){
        if((effect_size_phen_sig_plot2[i,j]>=myBreaks5[k-1])&(effect_size_phen_sig_plot2[i,j]<=myBreaks5[k])){
          temp<-c(temp,myColor5[k])
        }
      }
      temp<-c(temp,rownames(effect_size_phen_sig_plot2)[i],colnames(effect_size_phen_sig_plot2)[j])
      d<-temp
    }else{
      if(disp_fdr4[i,j]==""){
        temp<-c(i,j,NA)
      }else if(disp_fdr4[i,j]=="-"){
        temp<-c(i,j,45)
      }else if(disp_fdr4[i,j]=="+"){
        temp<-c(i,j,43)
      }
      if(effect_size_phen_sig_plot2[i,j]<myBreaks5[1]){
        temp<-c(temp,myColor5[1])
      }else if(effect_size_phen_sig_plot2[i,j]>myBreaks5[length(myColor5)]){
        temp<-c(temp,myColor5[length(myColor5)])
      }
      for(k in 2:length(myColor5)){
        if((effect_size_phen_sig_plot2[i,j]>=myBreaks5[k-1])&(effect_size_phen_sig_plot2[i,j]<=myBreaks5[k])){
          temp<-c(temp,myColor5[k])
          #tpp<-rbind(tpp,k)
        }
      }
      temp<-c(temp,rownames(effect_size_phen_sig_plot2)[i],colnames(effect_size_phen_sig_plot2)[j])
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
d$x1<-factor(d$x1,levels=c(dim(effect_size_phen_sig_plot2)[1]:-1:1))
d$y1<-factor(d$y1,levels=1:dim(effect_size_phen_sig_plot2)[2])
d$x11<-factor(d$x11,levels=rownames(effect_size_phen_sig_plot2)[c(length(rownames(effect_size_phen_sig_plot2)):-1:1)])
#d$y11<-factor(d$y11,levels=colnames(effect_size_phen_sig_plot2))
d$y11<-factor(d$y11,levels=colnames(effect_size_phen_sig_plot2)[out2$tree_col[["order"]]])

d2<-d
for(i in 1:dim(d)[1]){
  if(is.na(d$shp[i])==T){
    d2$clr[i]<-"#FFFFFF"
  }
}
#out2<-pheatmap(d2,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot2),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor5,breaks=myBreaks5,display_numbers = disp_fdr4,number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)


#"#377eb8","#e5d8bd","#33a02c","#6a3d9a","#a50026","#dd3497","#dfc27d","#ff7f00","#000000","#525252","#a6bddb","#cab2d6"

ggplot(d2,aes(y11,x11,fill=clr))+
  geom_point(shape=20,color=d2$clr,stroke=0,size=6)+
  #geom_point(shape=as.numeric(as.character(d$shp)),fill="black",size=5)+
  theme(axis.text.x=element_text(angle=90,hjust=0.95,color="black",face="bold",vjust=0.2),legend.position = "none",axis.text.y=element_text(color="black",face="bold"),panel.background = element_blank())+#plot.margin = margin(c(-20,(dim(effect_size_phen_sig_plot2)[2]+4),-40,(dim(effect_size_phen_sig_plot2)[1]+10)))
  xlab("")+ylab("")+expand_limits(x=c(-0.5,(dim(effect_size_phen_sig_plot2)[2]+10)))+#scale_x_discrete(1:(dim(effect_size_phen_sig_plot2)[2]+20))+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot2)[2]+1,y=dim(effect_size_phen_sig_plot2)[1]+0.5,yend=dim(effect_size_phen_sig_plot2)[1]+0.5,col="#4575b4",size=0.2)+
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot2)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot2)[2]+1),ymin=(dim(effect_size_phen_sig_plot2)[1]-3.5),ymax=(dim(effect_size_phen_sig_plot2)[1]+0.5),color="#4575b4",fill="#4575b4")+#+
  theme(plot.margin=margin(c(1,10,1,1), unit = "pt"))+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+4.3),y=(dim(effect_size_phen_sig_plot2)[1]-1.5),label="Pets",color="black",fontface="bold",hjust=1.2)+
  #annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+4.5),y=(dim(effect_size_phen_sig_plot2)[1]-33.7),label="habits",color="black",fontface="bold",hjust=0.6)+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot2)[1]+0.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-3.5),color="#4575b4",fill="#4575b4")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot2)[2]+1,y=dim(effect_size_phen_sig_plot2)[1]-3.5,yend=dim(effect_size_phen_sig_plot2)[1]-3.5,col="#ae017e",size=0.2)+#Farm animals
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot2)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot2)[2]+1),ymin=(dim(effect_size_phen_sig_plot2)[1]-14.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-3.5),color="#ae017e",fill="#ae017e")+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+3.5),y=(dim(effect_size_phen_sig_plot2)[1]-8.5),label="Farm",color="black",fontface="bold",hjust=0.75)+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+4),y=(dim(effect_size_phen_sig_plot2)[1]-9.7),label="animals",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot2)[1]-14.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-3.5),color="#ae017e",fill="#ae017e")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot2)[2]+1,y=dim(effect_size_phen_sig_plot2)[1]-14.5,yend=dim(effect_size_phen_sig_plot2)[1]-14.5,col="#d94801",size=0.2)+#Wild animals
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot2)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot2)[2]+1),ymin=(dim(effect_size_phen_sig_plot2)[1]-23.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-14.5),color="#d94801",fill="#d94801")+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+3.5),y=(dim(effect_size_phen_sig_plot2)[1]-18.5),label="Wild",color="black",fontface="bold",hjust=0.75)+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+4),y=(dim(effect_size_phen_sig_plot2)[1]-19.7),label="animals",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot2)[1]-23.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-14.5),color="#d94801",fill="#d94801")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot2)[2]+1,y=dim(effect_size_phen_sig_plot2)[1]-23.5,yend=dim(effect_size_phen_sig_plot2)[1]-23.5,col="#000000",size=0.2)+#Food
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot2)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot2)[2]+1),ymin=(dim(effect_size_phen_sig_plot2)[1]-40.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-23.5),color="#000000",fill="#000000")+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+3.5),y=(dim(effect_size_phen_sig_plot2)[1]-34),label="Food",color="black",fontface="bold",hjust=0.75)+
  #annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+4),y=(dim(effect_size_phen_sig_plot2)[1]-19.7),label="animals",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot2)[1]-40.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-23.5),color="#000000",fill="#000000")

#Save 10.00 x 12.00, y axis font =10, spacing =1.12 ---- instead try 9(or 9.5) x 12
#Rerun food and animals for equal axis

#RUnning wihtout no pets, no farm, no wild

d3<-d2[which((d2$x11!="No pets (N=88)")&(d2$x11!="No farm animals (N=69)")&(d2$x11!="No wild animals (N=243)")),]

ggplot(d3,aes(y11,x11,fill=clr))+
  geom_point(shape=20,color=d3$clr,stroke=0,size=6)+
  #geom_point(shape=as.numeric(as.character(d$shp)),fill="black",size=5)+
  theme(axis.text.x=element_text(angle=90,hjust=0.95,color="black",face="bold",vjust=0.2),legend.position = "none",axis.text.y=element_text(color="black",face="bold"),panel.background = element_blank())+#plot.margin = margin(c(-20,(dim(effect_size_phen_sig_plot2)[2]+4),-40,(dim(effect_size_phen_sig_plot2)[1]+10)))
  xlab("")+ylab("")+expand_limits(x=c(-0.5,(dim(effect_size_phen_sig_plot2)[2]+10)))+#scale_x_discrete(1:(dim(effect_size_phen_sig_plot2)[2]+20))+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot2)[2]+1,y=dim(effect_size_phen_sig_plot2)[1]-2.5,yend=dim(effect_size_phen_sig_plot2)[1]-2.5,col="#4575b4",size=0.2)+
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot2)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot2)[2]+1),ymin=(dim(effect_size_phen_sig_plot2)[1]-5.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-2.5),color="#4575b4",fill="#4575b4")+#+
  theme(plot.margin=margin(c(1,10,1,1), unit = "pt"))+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+4.3),y=(dim(effect_size_phen_sig_plot2)[1]-4),label="Pets",color="black",fontface="bold",hjust=1.2)+
  #annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+4.5),y=(dim(effect_size_phen_sig_plot2)[1]-33.7),label="habits",color="black",fontface="bold",hjust=0.6)+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot2)[1]+0.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-5.5),color="#4575b4",fill="#4575b4")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot2)[2]+1,y=dim(effect_size_phen_sig_plot2)[1]-5.5,yend=dim(effect_size_phen_sig_plot2)[1]-5.5,col="#ae017e",size=0.2)+#Farm animals
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot2)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot2)[2]+1),ymin=(dim(effect_size_phen_sig_plot2)[1]-15.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-5.5),color="#ae017e",fill="#ae017e")+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+3.5),y=(dim(effect_size_phen_sig_plot2)[1]-10),label="Farm",color="black",fontface="bold",hjust=0.75)+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+4),y=(dim(effect_size_phen_sig_plot2)[1]-11.2),label="animals",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot2)[1]-15.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-5.5),color="#ae017e",fill="#ae017e")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot2)[2]+1,y=dim(effect_size_phen_sig_plot2)[1]-15.5,yend=dim(effect_size_phen_sig_plot2)[1]-15.5,col="#d94801",size=0.2)+#Wild animals
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot2)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot2)[2]+1),ymin=(dim(effect_size_phen_sig_plot2)[1]-23.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-15.5),color="#d94801",fill="#d94801")+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+3.2),y=(dim(effect_size_phen_sig_plot2)[1]-18.5),label="Wild",color="black",fontface="bold",hjust=0.75)+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+4),y=(dim(effect_size_phen_sig_plot2)[1]-19.7),label="animals",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot2)[1]-23.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-15.5),color="#d94801",fill="#d94801")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot2)[2]+1,y=dim(effect_size_phen_sig_plot2)[1]-23.5,yend=dim(effect_size_phen_sig_plot2)[1]-23.5,col="#000000",size=0.2)+#Food
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot2)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot2)[2]+1),ymin=(dim(effect_size_phen_sig_plot2)[1]-40.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-23.5),color="#000000",fill="#000000")+
  annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+3.5),y=(dim(effect_size_phen_sig_plot2)[1]-34),label="Food",color="black",fontface="bold",hjust=0.75)+
  #annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+4),y=(dim(effect_size_phen_sig_plot2)[1]-19.7),label="animals",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot2)[1]-40.5),ymax=(dim(effect_size_phen_sig_plot2)[1]-23.5),color="#000000",fill="#000000")

#Save 9.00 x 12.00



##Socio-economic phenotypes

##RUN AFTER running socio-economic in dmp_lmer
disp_fdr5<-t(disp_fdr[jkl33[xt],])
myColor#Color
myBreaks11#breaks
#myBreaks6<-myBreaks2[c(1:7,9:length(myBreaks2))]#c(myBreaks[1:6],-0.36,-0.001,myBreaks[9:length(myBreaks)])
#myColor6<-myColor#c(myColor[1:6],"#F7DDDE","#F7DDDE",myColor[8:length(myColor)])#was "#F0BCBE","#F7DDDE"
myBreaks6<-myBreaks11
myColor6<-myColor


effect_size_phen_sig_plot3<-effect_size_phen_sig_plot2[,jkl33[xt]]
# geom_segment(x=1,y=1*0.95+coef(op2)[1],xend=-1,yend=-1*0.95+coef(op2)[1],col="red", size=1)+

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
      if(effect_size_phen_sig_plot3[i,j]<myBreaks6[1]){
        temp<-c(temp,myColor6[1])
      }else if(effect_size_phen_sig_plot3[i,j]>myBreaks6[length(myColor6)]){
        temp<-c(temp,myColor6[length(myColor6)])
      }
      for(k in 2:length(myColor6)){
        if((effect_size_phen_sig_plot3[i,j]>=myBreaks6[k-1])&(effect_size_phen_sig_plot3[i,j]<=myBreaks6[k])){
          temp<-c(temp,myColor6[k])
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
      if(effect_size_phen_sig_plot3[i,j]<myBreaks6[1]){
        temp<-c(temp,myColor6[1])
      }else if(effect_size_phen_sig_plot3[i,j]>myBreaks6[length(myColor6)]){
        temp<-c(temp,myColor6[length(myColor6)])
      }
      for(k in 2:length(myColor6)){
        if((effect_size_phen_sig_plot3[i,j]>=myBreaks6[k-1])&(effect_size_phen_sig_plot3[i,j]<=myBreaks6[k])){
          temp<-c(temp,myColor6[k])
          #tpp<-rbind(tpp,k)
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
#d$y11<-factor(d$y11,levels=colnames(effect_size_phen_sig_plot2))
d$y11<-factor(d$y11,levels=colnames(effect_size_phen_sig_plot3)[out2$tree_col[["order"]]])

d2<-d
for(i in 1:dim(d)[1]){
  if(is.na(d$shp[i])==T){
    d2$clr[i]<-"#FFFFFF"
  }
}
#out2<-pheatmap(d2,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot2),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor5,breaks=myBreaks5,display_numbers = disp_fdr4,number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)


#"#377eb8","#e5d8bd","#33a02c","#6a3d9a","#a50026","#dd3497","#dfc27d","#ff7f00","#000000","#525252","#a6bddb","#cab2d6"

ggplot(d2,aes(y11,x11,fill=clr))+
  geom_point(shape=20,color=d2$clr,stroke=0,size=6)+
  #geom_point(shape=as.numeric(as.character(d$shp)),fill="black",size=5)+
  theme(axis.text.x=element_text(angle=90,hjust=0.95,color="black",face="bold",vjust=0.2),legend.position = "none",axis.text.y=element_text(color="black",face="bold"),panel.background = element_blank())+#plot.margin = margin(c(-20,(dim(effect_size_phen_sig_plot3)[2]+4),-40,(dim(effect_size_phen_sig_plot3)[1]+10)))
  xlab("")+ylab("")+expand_limits(x=c(-0.5,(dim(effect_size_phen_sig_plot3)[2]+10)))+#scale_x_discrete(1:(dim(effect_size_phen_sig_plot3)[2]+20))+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot3)[2]+1,y=dim(effect_size_phen_sig_plot3)[1]+0.5,yend=dim(effect_size_phen_sig_plot3)[1]+0.5,col="#4575b4",size=0.3)+#Partner
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot3)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot3)[2]+1),ymin=(dim(effect_size_phen_sig_plot3)[1]-3.5),ymax=(dim(effect_size_phen_sig_plot3)[1]+0.5),color="#4575b4",fill="#4575b4")+#+
  theme(plot.margin=margin(c(1,12,1,1), unit = "pt"))+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+4.7),y=(dim(effect_size_phen_sig_plot3)[1]-1),label="Co-habiting",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+4.5),y=(dim(effect_size_phen_sig_plot3)[1]-2.2),label="partners",color="black",fontface="bold",hjust=0.6)+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot3)[1]+0.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-3.5),color="#4575b4",fill="#4575b4")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot3)[2]+1,y=dim(effect_size_phen_sig_plot3)[1]-3.5,yend=dim(effect_size_phen_sig_plot3)[1]-3.5,col="#525252",size=0.3)+#Education
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot3)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot3)[2]+1),ymin=(dim(effect_size_phen_sig_plot3)[1]-6.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-3.5),color="#525252",fill="#525252")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+4.3),y=(dim(effect_size_phen_sig_plot3)[1]-5),label="Education",color="black",fontface="bold")+
  #annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+4),y=(dim(effect_size_phen_sig_plot3)[1]-9.7),label="animals",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot3)[1]-6.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-3.5),color="#525252",fill="#525252")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot3)[2]+1,y=dim(effect_size_phen_sig_plot3)[1]-6.5,yend=dim(effect_size_phen_sig_plot3)[1]-6.5,col="#d94801",size=0.3)+#Friendship
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot3)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot3)[2]+1),ymin=(dim(effect_size_phen_sig_plot3)[1]-10.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-6.5),color="#d94801",fill="#d94801")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+4.5),y=(dim(effect_size_phen_sig_plot3)[1]-8),label="Friendship",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+4.5),y=(dim(effect_size_phen_sig_plot3)[1]-9.2),label="(non-kin)",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot3)[1]-10.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-6.5),color="#d94801",fill="#d94801")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot3)[2]+1,y=dim(effect_size_phen_sig_plot3)[1]-10.5,yend=dim(effect_size_phen_sig_plot3)[1]-10.5,col="#e7298a",size=0.3)+#Family
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot3)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot3)[2]+1),ymin=(dim(effect_size_phen_sig_plot3)[1]-14.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-10.5),color="#e7298a",fill="#e7298a")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+4.7),y=(dim(effect_size_phen_sig_plot3)[1]-12.5),label="Family (kin)",color="black",fontface="bold")+
  #annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+4),y=(dim(effect_size_phen_sig_plot3)[1]-19.7),label="animals",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot3)[1]-14.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-10.5),color="#e7298a",fill="#e7298a")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot3)[2]+1,y=dim(effect_size_phen_sig_plot3)[1]-14.5,yend=dim(effect_size_phen_sig_plot3)[1]-14.5,col="#807dba",size=0.3)+#All
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot3)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot3)[2]+1),ymin=(dim(effect_size_phen_sig_plot3)[1]-18.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-14.5),color="#807dba",fill="#807dba")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+3.8),y=(dim(effect_size_phen_sig_plot3)[1]-15.5),label="All ties",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+6.5),y=(dim(effect_size_phen_sig_plot3)[1]-16.7),label="(kin + non-kin)",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot3)[1]-18.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-14.5),color="#807dba",fill="#807dba")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot3)[2]+1,y=dim(effect_size_phen_sig_plot3)[1]-18.5,yend=dim(effect_size_phen_sig_plot3)[1]-18.5,col="#fcadad",size=0.3)+#Risky behavior
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot3)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot3)[2]+1),ymin=(dim(effect_size_phen_sig_plot3)[1]-21.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-18.5),color="#fcadad",fill="#fcadad")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+3.0),y=(dim(effect_size_phen_sig_plot3)[1]-19.5),label="Risky",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+4.5),y=(dim(effect_size_phen_sig_plot3)[1]-20.7),label="behavior",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot3)[1]-21.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-18.5),color="#fcadad",fill="#fcadad")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot3)[2]+1,y=dim(effect_size_phen_sig_plot3)[1]-21.5,yend=dim(effect_size_phen_sig_plot3)[1]-21.5,col="#00b0f0",size=0.3)+#Village factors
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot3)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot3)[2]+1),ymin=(dim(effect_size_phen_sig_plot3)[1]-25.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-21.5),color="#00b0f0",fill="#00b0f0")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+3.5),y=(dim(effect_size_phen_sig_plot3)[1]-23),label="Village",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+3.8),y=(dim(effect_size_phen_sig_plot3)[1]-24.2),label="factors",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot3)[1]-25.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-21.5),color="#00b0f0",fill="#00b0f0")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot3)[2]+1,y=dim(effect_size_phen_sig_plot3)[1]-25.5,yend=dim(effect_size_phen_sig_plot3)[1]-25.5,col="#7030a0",size=0.3)+#Income
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot3)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot3)[2]+1),ymin=(dim(effect_size_phen_sig_plot3)[1]-27.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-25.5),color="#7030a0",fill="#7030a0")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+3.7),y=(dim(effect_size_phen_sig_plot3)[1]-26.5),label="Income",color="black",fontface="bold")+
  #annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+3.8),y=(dim(effect_size_phen_sig_plot3)[1]-24.2),label="factors",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot3)[1]-27.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-25.5),color="#7030a0",fill="#7030a0")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot3)[2]+1,y=dim(effect_size_phen_sig_plot3)[1]-27.5,yend=dim(effect_size_phen_sig_plot3)[1]-27.5,col="#ffc000",size=0.3)+#HH essentials
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot3)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot3)[2]+1),ymin=(dim(effect_size_phen_sig_plot3)[1]-40.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-27.5),color="#ffc000",fill="#ffc000")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+4.8),y=(dim(effect_size_phen_sig_plot3)[1]-33.5),label="Household",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+5.1),y=(dim(effect_size_phen_sig_plot3)[1]-34.7),label="essentials",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot3)[1]-40.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-27.5),color="#ffc000",fill="#ffc000")+
  geom_segment(x=0,xend=dim(effect_size_phen_sig_plot3)[2]+1,y=dim(effect_size_phen_sig_plot3)[1]-40.5,yend=dim(effect_size_phen_sig_plot3)[1]-40.5,col="#000000",size=0.3)+#Water sources
  geom_rect(xmin=c(dim(effect_size_phen_sig_plot3)[2]+0.5),xmax=c(dim(effect_size_phen_sig_plot3)[2]+1),ymin=(dim(effect_size_phen_sig_plot3)[1]-44.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-40.5),color="#000000",fill="#000000")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+4.8),y=(dim(effect_size_phen_sig_plot3)[1]-41.5),label="Water",color="black",fontface="bold")+
  annotate("text",x=(dim(effect_size_phen_sig_plot3)[2]+5.1),y=(dim(effect_size_phen_sig_plot3)[1]-42.7),label="sources",color="black",hjust=0.6,fontface="bold")+
  geom_rect(xmin=c(-0.5),xmax=0,ymin=(dim(effect_size_phen_sig_plot3)[1]-44.5),ymax=(dim(effect_size_phen_sig_plot3)[1]-40.5),color="#000000",fill="#000000")

#Save 9.5 x 13.00










