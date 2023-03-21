library(pheatmap)

##Figure 1 

##Google maps plot
library(ggmap)
register_google(key="###")
#remotes::install_github("dyerlab/popgraph")##NEEDS too many requirements, need new laptop?
library('tidygraph')
library(ggplot2)
library('igraph')
library(vegan)

#register_google(key="###")


##

g_lon<-read.csv('g_lon.csv',row.names = 1)
g_lat<-read.csv('g_lat.csv',row.names = 1)
g_lon2<-read.csv('g_lon2.csv',row.names = 1)
g_lat2<-read.csv('g_lat2.csv',row.names = 1)
g_lon<-g_lon[,1]
g_lon2<-g_lon2[,1]
g_lat<-g_lat[,1]
g_lat2<-g_lat2[,1]
bray_col<-read.csv('bray_col.csv',row.names = 1)
edges<-read.csv('ggmap_edges.csv')
edges<-edges[,2:dim(edges)[2]]
bn_shape<-read.csv('bn_shape.csv',row.names = 1)
bray_col<-bray_col[,1]
bn_shape<-as.character(bn_shape[,1])
vil_lon<-read.csv('vil_lon.csv',row.names = 1)
vil_lat<-read.csv('vil_lat.csv',row.names = 1)
vil_cood<-read.csv('vil_cood.csv',row.names = 1)
vil_lon<-vil_lon[,1]
vil_lat<-vil_lat[,1]

houses <- get_map(location = c(mean(g_lon), mean(g_lat)), source="google", zoom=17,maptype="hybrid")#map_type- terrain-background , satellite (USE) or hybrid,terrain#Gives with no labels but yellow color--- maybe the one we want?
p<-ggmap(houses)

p
plot_vector<- as.data.frame(cbind(g_lon,g_lat))
p + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = 0.9, colour="white",alpha=0.8) + geom_point(aes(g_lon2,g_lat2), data=plot_vector)+ geom_point(aes(x=g_lon2,y=g_lat2),plot_vector,size=4,alpha=1,color=adjustcolor(bray_col, alpha.f = 1),shape=bn_shape)+coord_fixed(ylim=c(14.9448,14.9475), ratio=1/cos(pi*14.9448/180))+xlab("Longitude")+ylab("Latitude")+theme(text=element_text(size=20))

#p + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = 0.7, colour="white",alpha=0.8) + geom_point(aes(g_lon2,g_lat2), data=plot_vector)+ geom_point(aes(x=g_lon2,y=g_lat2),plot_vector,size=3,alpha=1,color=adjustcolor(bray_col, alpha.f = 1),shape=bn_shape)+coord_fixed(ylim=c(14.9438,14.9485), ratio=1/cos(pi*14.9448/180))+xlab("Longitude")+ylab("Latitude")+theme(text=element_text(size=20))

#Village plotting
houses <- get_map(location = c( mean(vil_lon),mean(vil_lat)), source="google", zoom=8,maptype="hybrid",style = c(feature="country",element="labels",visibility='off'))#map_type- terrain-background , satellite (USE) or hybrid,terrain#Gives with no labels but yellow color--- maybe the one we want?
p<-ggmap(houses)
p

p + geom_point(aes(x=vil_lon,y=vil_lat),vil_cood,size=2,alpha=1,color=adjustcolor(vil_cood$col_vil, alpha.f = 1))+coord_fixed(ylim=c(13.7,15.8),xlim=c(-90.8,-88.1), ratio=1/cos(pi*14.9448/180))+xlab("Longitude")+ylab("Latitude")+theme(text=element_text(size=15))#color=adjustcolor(age_col, alpha.f = 1),

#p + geom_point(aes(x=vil_lon,y=vil_lat),vil_cood,size=7,alpha=1,color="#74add1",shape="*")+coord_fixed(ylim=c(13.7,15.8),xlim=c(-90.8,-88.1), ratio=1/cos(pi*14.9448/180))+xlab("Longitude")+ylab("Latitude")+theme(text=element_text(size=15))#color=adjustcolor(age_col, alpha.f = 1),

#p + geom_point(aes(x=vil_lon,y=vil_lat),vil_cood,size=7,alpha=1,color="#80cdc1",shape="*")+coord_fixed(ylim=c(13.7,15.8),xlim=c(-90.8,-88.1), ratio=1/cos(pi*14.9448/180))+xlab("Longitude")+ylab("Latitude")+theme(text=element_text(size=15))#color=adjustcolor(age_col, alpha.f = 1),

p + geom_point(aes(x=vil_lon,y=vil_lat),vil_cood,size=2,alpha=1,color="black")+coord_fixed(ylim=c(13.7,15.8),xlim=c(-90.8,-88.1), ratio=1/cos(pi*14.9448/180))+xlab("Longitude")+ylab("Latitude")+theme(text=element_text(size=15))+ geom_point(aes(x=vil_lon,y=vil_lat),vil_cood,size=1.2,alpha=1,color="#fdae61")#color=adjustcolor(age_col, alpha.f = 1),

houses <- get_map(location = c( mean(vil_lon),mean(vil_lat)), source="google", zoom=9,maptype="hybrid",style = c(feature="country",element="labels",visibility='off'))#map_type- terrain-background , satellite (USE) or hybrid,terrain#Gives with no labels but yellow color--- maybe the one we want?
p<-ggmap(houses)
p
p + geom_point(aes(x=vil_lon,y=vil_lat),vil_cood,size=2,alpha=1,color="black")+xlab("Longitude")+ylab("Latitude")+theme(text=element_text(size=15))+ geom_point(aes(x=vil_lon,y=vil_lat),vil_cood,size=1.2,alpha=1,color="#fdae61")+geom_label(label="El Salvador",x=-89.38,y=14.35,size=3.1,fontface="bold",family="sans",label.padding = unit(0.15,"lines"))+#color=adjustcolor(age_col, alpha.f = 1),
  geom_label(label="Guatemala",x=-89.6,y=14.5,size=3.1,fontface="bold",family="sans",label.padding = unit(0.15,"lines"))+
  geom_label(label="Honduras",x=-89,y=14.61,size=3.1,fontface="bold",family="sans",label.padding = unit(0.15,"lines"))

##Other figure 1 plots -- regression line

###Figure 1C

library(vegan)
phen_all_use<-read.csv('phen_all_use6_all.csv',row.names=1)
vc<-unique(phen_all_use$village_code)
mb_samp_sp<-read.csv('mb_samp_1188_3.csv',row.names=1)

br<-array(NA,dim=c(dim(phen_all_use)[1],1))

#ind_vil<-which(phen_all_use$village_code==vc[5])
#mb_mean_vil<-rowSums(mb_samp_sp[,ind_vil])/length(ind_vil)

#kl_gen<-array(NA,dim=c(dim(br)[1],1))

j<-1
for(j in 1:dim(phen_all_use)[1]){
  tp<-which(vc==phen_all_use$village_code[j])
  ind_vil<-which(phen_all_use$village_code==vc[tp])
  ind_vil<-ind_vil[ind_vil!=j]
  mb_mean_vil<-rowSums(mb_samp_sp[,ind_vil])/length(ind_vil)
  br[j,1]<-as.numeric(vegdist(t(cbind(mb_mean_vil,mb_samp_sp[,j])), "bray"))
  
}

#phen_all_use$bray_vil_mean<-br
#phen_all_use<-cbind(phen_all_use,"bray_vil_mean"=br[,1])

#write.csv(phen_all_use,'phen_all_use3_all.csv')

##Plots w.r.t Popluation weighted center
library(ggplot2)
unique(phen_all_use$village_name[which(phen_all_use$village_code==vc[6])])
ind_vil<-which(phen_all_use$village_code==vc[6])
phen_all_use2<-phen_all_use[ind_vil,]

ggplot(phen_all_use2, aes(x=distance_center, y=bray_vil_mean)) + geom_point()+xlab("Distance from village center")+ylab("Bray-Curtis dissimilarity")#+ylim(0,1)

#FIgure 1B
#All 1187 samples
library(RColorBrewer)
ggplot(phen_all_use, aes(x=distance_center, y=bray_vil_mean)) + geom_point(aes(color=as.character(village_code)))+xlim(0,2)+xlab("Distance from village center")+ylab("Bray-curtis dissimilarity")+ scale_color_brewer(palette="Paired")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())# +scale_x_continuous(trans='log2') +

ggplot(phen_all_use, aes(x=distance_center_build, y=bray_vil_mean)) + geom_point(aes(color=as.character(village_code)),size=1)+xlim(0,2)+xlab("Distance from unweighted village center(km)")+ylab("Bray-curtis dissimilarity")+ scale_color_brewer(palette="Paired")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),legend.key = element_rect(colour = "white", fill = NA))+# +scale_x_continuous(trans='log2') +
  geom_smooth(method='lm',formula=y~x,color="black")

ggplot(phen_all_use, aes(x=distance_center, y=bray_vil_mean)) + geom_point(aes(color=as.character(village_code)),size=1)+xlim(0,2)+xlab("Distance from weighted village center(km)")+ylab("Bray-Curtis dissimilarity")+ scale_color_brewer(palette="Paired")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),text = element_text(size = 20),legend.position="none",axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black"))+# +scale_x_continuous(trans='log2') +#legend.key = element_rect(colour = "white", fill = NA)
  geom_smooth(method='lm',formula=y~x,color="black")+guides(fill=guide_legend(title="Village code"))+ labs(fill = "Village code")




##########################################################################################################################
##Figure 2

##Read in all the files

##Association

paletteLength<-15
myColor <- colorRampPalette(c("#cb181d","white", "#33a02c"))(paletteLength)#myColor <- colorRampPalette(c("Orange", "white", "Darkblue"))(paletteLength)
myBreaks <- c(seq(min(effect_size_phen_sig_plot),0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(effect_size_phen_sig_plot)/paletteLength, max(effect_size_phen_sig_plot), length.out=floor(paletteLength/2)))

myBreaks
myBreaks[7]<-as.numeric(-0.001)
myBreaks[9]<-as.numeric(0.001)
myBreaks

#Size for saving - 10.0 x 14.00

pheatmap(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 13,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

#With phylum annotation

pheatmap(effect_size_phen_sig_plot2,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr2),number_color = "black",treeheight_row = 0,cluster_rows = F,cluster_cols = T,treeheight_col = 0,annotation_col=phyl_ph)

##Shannon diversity

#Read in
div_alpha_ph2<-read.csv('div_alpha_ph2_chronic.csv',row.names=1)

#Chronic conditions

div_alpha_ph2$chronic <- factor(div_alpha_ph2$chronic, levels=c("Healthy(848)","Diabetes(24)","Allergies(104)","Heart disease(47)","Asthma(47)","Stomach illness(175)","Intestinal illness(70)","Arthritis(39)"),ordered=TRUE)
library(ggplot2)
ggplot(div_alpha_ph2, aes(x=chronic, y=as.numeric(as.character(alpha_div)),group=chronic)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#1f78b4","#b2df8a","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6"),outlier.alpha=0)+labs(x="",y="Alpha diversity")+
  theme(text = element_text(size=15))+ scale_x_discrete(breaks=seq(1,8,1))+ylim(1.4,5.2)#,"#ff7f00","#cab2d6","#6a3d9a","#ffff99"#x="Relationship",

x1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$chronic=="Healthy(848)")]
y1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$chronic=="Allergies(104)")]
x1<-x1[which(is.na(x1)==F)]
y1<-y1[which(is.na(y1)==F)]

wilcox.test(as.numeric(as.character(x1)),as.numeric(as.character(y1)))

library(ggpubr)
my_comparisons <- list( c("Healthy(848)", "Allergies(104)"), c("Healthy(848)", "Stomach illness(175)"), c("Healthy(848)", "Intestinal illness(70)") )

#Save 4x4 pdf for Illana plot
div_alpha_ph22<-div_alpha_ph2[which(is.na(div_alpha_ph2$alpha_div)==F),]
div_mm<-mean(div_alpha_ph22$alpha_div[which(div_alpha_ph22$chronic=="Healthy(848)")])
#Save pdf - 5.00 x 4.00
ggplot(div_alpha_ph2, aes(x=chronic, y=as.numeric(as.character(alpha_div)),group=chronic)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#1f78b4","#b2df8a","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6"),outlier.alpha=0)+labs(x="",y="Shannon diversity")+
  theme(text = element_text(size=20),axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black"),panel.background=element_blank(),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))+ scale_x_discrete(breaks=seq(1,8,1))+
  stat_compare_means(comparisons = my_comparisons, label.y = c(5.2, 5.5, 5.8),size=5)+ylim(1.5,6.5)+#+
  geom_segment(x=0,xend=9,y=(div_mm+0.05),yend=(div_mm+0.05),col="#878787",linetype="longdash",size=1)+#+
  theme(plot.margin=margin(c(1,0.5,2,1), unit = "cm"))
#  annotate("text",x=1,y=(1),label="Healthy",color="black",fontface="bold")

#Save pdf - 6.00 x 4.00 (for below)
ggplot(div_alpha_ph2, aes(x=chronic, y=as.numeric(as.character(alpha_div)),group=chronic)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#1f78b4","#b2df8a","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6"),outlier.alpha=0)+labs(x="",y="Shannon diversity")+
  theme(text = element_text(size=20),axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black"),panel.background=element_blank(),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))+ scale_x_discrete(breaks=seq(1,8,1))+
  stat_compare_means(comparisons = my_comparisons, label.y = c(5.2, 5.5, 5.8),size=5)+ylim(1.5,6.5)+#+
  geom_segment(x=0,xend=9,y=(div_mm+0.05),yend=(div_mm+0.05),col="#878787",linetype="longdash",size=1)+#+
  theme(plot.margin=margin(c(2,0.5,4,1), unit = "cm"))

  
library(svglite)
svglite('Figure 2B 4.svg', width=4, height=4)
ggplot(div_alpha_ph2, aes(x=chronic, y=as.numeric(as.character(alpha_div)),group=chronic)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#1f78b4","#b2df8a","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6"),outlier.alpha=0)+labs(x="",y="Shannon diversity")+
  theme(text = element_text(size=20),axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black"),panel.background=element_blank(),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))+ scale_x_discrete(breaks=seq(1,8,1))+
  stat_compare_means(comparisons = my_comparisons, label.y = c(5.2, 5.5, 5.8),size=5)+ylim(1.5,6.5)
dev.off()

#Medications
div_alpha_ph3<-read.csv('div_alpha_ph2_medication.csv',row.names = 1)

div_alpha_ph3$medication <- factor(div_alpha_ph3$medication, levels=c("No medication","Antibiotic","Anti-diarrheal","Anti-parasitic","Anti-fungal","Vitamins","Anti-hypertensive"),ordered=TRUE)

ggplot(div_alpha_ph3, aes(x=medication, y=as.numeric(as.character(alpha_div)),group=medication)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f"),outlier.alpha = 0)+labs(x="",y="Alpha diversity")+
  theme(text = element_text(size=15))+ scale_x_discrete(breaks=seq(1,7,1))+ylim(1.5,5.2)#,"#ff7f00","#cab2d6","#6a3d9a","#ffff99"#x="Relationship",

x1<-div_alpha_ph3$alpha_div[which(div_alpha_ph3$medication=="No medication")]
y1<-div_alpha_ph3$alpha_div[which(div_alpha_ph3$medication=="Anti-hypertensive")]
x1<-x1[which(is.na(x1)==F)]
y1<-y1[which(is.na(y1)==F)]

wilcox.test(as.numeric(as.character(x1)),as.numeric(as.character(y1)))

library(ggpubr)
my_comparisons <- list( c("No medication", "Antibiotic"),c("No medication", "Anti-diarrheal"),c("No medication", "Anti-parasitic"))
#bray_chronic2<-bray_chronic[which(is.na(bray_chronic[,1])==F),]
ggplot(div_alpha_ph3, aes(x=medication, y=as.numeric(as.character(alpha_div)),group=medication)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f"),outlier.alpha=0)+labs(x="",y="Shannon diversity")+
  theme(text = element_text(size=20),axis.text.y = element_text(color="black",size=15),axis.text.x = element_text(color="black"),panel.background=element_blank(),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))+ scale_x_discrete(breaks=seq(1,7,1))+
  stat_compare_means(comparisons = my_comparisons, label.y = c(5.2,5.5,5.8),size=5)+ylim(1.8,6.1)#+

ggplot(div_alpha_ph3, aes(x=medication, y=as.numeric(as.character(alpha_div)),group=medication)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f"),outlier.alpha=0)+labs(x="",y="Shannon diversity")+
  theme(text = element_text(size=20),axis.text.y = element_text(color="black",size=15),axis.text.x = element_text(color="black"),panel.background=element_blank(),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))+ scale_x_discrete(breaks=seq(1,7,1))+
  stat_compare_means(comparisons = my_comparisons, label.y = c(5.2,5.5,5.8),size=5)+ylim(1.8,6.3)+#+
  theme(plot.margin=margin(c(2,0.5,4,1), unit = "cm"))

library(svglite)
svglite('Figure 2C 4.svg', width=4, height=4)
ggplot(div_alpha_ph3, aes(x=medication, y=as.numeric(as.character(alpha_div)),group=medication)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f"),outlier.alpha=0)+labs(x="",y="Shannon diversity")+
  theme(text = element_text(size=20),axis.text.y = element_text(color="black",size=15),axis.text.x = element_text(color="black"),panel.background=element_blank(),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))+ scale_x_discrete(breaks=seq(1,7,1))+
  stat_compare_means(comparisons = my_comparisons, label.y = c(5.2,5.5,5.8),size=5)+ylim(1.8,6.1)#+
dev.off()


##Healthy vs unhealthy

tt2<-read.csv('health_vs_unhealth_tt2.csv',row.names=1)

efz<-tt2$effect_size[which(tt2$FDR<0.05)]
efz_names<-as.character(tt2$sp_name[which(tt2$FDR<0.05)])
efz2<-efz[order(efz)]
efz_names2<-efz_names[order(efz)]

fdr<-tt2$FDR[which(tt2$FDR<0.05)]
fdr2<-fdr[order(efz)]

library(ggplot2)
#library(ggtext)
library(mdthemes)

#health_bar<-cbind(cbind(cbind(efz2[c(1:10,(length(efz)-9):length(efz))],fdr2[c(1:10,(length(efz)-9):length(efz))]),efz_names2[c(1:10,(length(efz)-9):length(efz))]),c(array(1,dim=c(10,1)),array(2,dim=c(10,1))))
health_bar<-cbind(cbind(cbind(efz2[c(1:10,(length(efz)-9):length(efz))],fdr2[c(1:10,(length(efz)-9):length(efz))]),as.character(c("Bacteroides ovatus (SGB1871)","Bacteroides caccae (SGB1877)","Parabacteroides distasonis (SGB1934)","Alistipes putredinis (SGB2318)","Bacteroides xylanisolvens (SGB1867)","Paraprevotella clara (SGB1798)","Phocaeicola vulgatus (SGB1814)","Bacteroides uniformis (SGB1836_group)","Bacteroides thetaiotaomicron (SGB1861)","Clostridium sp. AM22 11AC (SGB4749)","Prevotella sp. 885 (SGB1677)","{f__Bacilli unclassified}GGB4672 (SGB6461)","Bacilli unclassified (SGB6428)","{f__Bacteroidaceae}GGB1380 (SGB1883 group)","{f__Prevotellaceae}GGB1243 (SGB1663)","{f__Rikenellaceae}GGB1632 (SGB2240)","Coprococcus sp. (SGB5115)","{f__Rikenellaceae}GGB1631 (SGB2239)","{f__Prevotellaceae}GGB1247 (SGB1668)","{f__Rikenellaceae}GGB1630 (SGB2238)"))),c(array("Disease associated",dim=c(10,1)),array("Healthy",dim=c(10,1))))
colnames(health_bar)<-c("effect","fdr","name","category")
health_bar<-as.data.frame(health_bar)

stp<-c("***","***","**","**","***","**","**","**","**","**","*","***","***","**","**","**","***","***","***","***")
#7,14.0
ggplot(health_bar,aes(x=reorder(as.character(name),as.numeric(as.character(effect))),y=as.numeric(as.character(effect)),fill=category))+
  geom_bar(stat="identity",color="black",size=1)+scale_fill_manual(values = c("#4575b4","#d73027"))+#+stat_bin(geom="text",aes(),vjust=-1)
  theme(text = element_text(size=20),panel.grid.major = element_blank(),panel.background = element_blank(),axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black"),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),legend.position = "top",legend.text = element_text(size=20)) +
  ylab("Effect size")+coord_flip()+xlab("")#+mdthemes::md_theme_classic()

h<-ggplot(health_bar,aes(x=reorder(as.character(name),as.numeric(as.character(effect))),y=as.numeric(as.character(effect)),fill=category))+
  geom_bar(stat="identity",color="black",size=1)+scale_fill_manual(values = c("#4575b4","#d73027"))+#+stat_bin(geom="text",aes(),vjust=-1)
  theme(text = element_text(size=20),panel.grid.major = element_blank(),panel.background = element_blank(),axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black"),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),legend.position = "top",legend.text = element_text(size=20)) +
  ylab("Effect size")+coord_flip()+xlab("")

svglite('Figure 2E3.svg', width=12, height=6)
plot(h)
dev.off()

health_bar<-cbind(cbind(cbind(efz2[c(1:10,(length(efz)-9):length(efz))],fdr2[c(1:10,(length(efz)-9):length(efz))]),as.character(c("*Bacteroides ovatus* (SGB1871)","*Bacteroides caccae* (SGB1877)","*Parabacteroides distasonis* (SGB1934)","*Alistipes putredinis* (SGB2318)","*Bacteroides xylanisolvens* (SGB1867)","*Paraprevotella clara* (SGB1798)","*Phocaeicola vulgatus* (SGB1814)","*Bacteroides uniformis* (SGB1836_group)","*Bacteroides thetaiotaomicron* (SGB1861)","*Clostridium sp. AM22 11AC* (SGB4749)","*Prevotella sp. 885* (SGB1677)","{f__Bacilli unclassified}GGB4672 (SGB6461)","*Bacilli unclassified* (SGB6428)","{f__Bacteroidaceae}GGB1380 (SGB1883 group)","{f__Prevotellaceae}GGB1243 (SGB1663)","{f__Rikenellaceae}GGB1632 (SGB2240)","*Coprococcus sp.* (SGB5115)","{f__Rikenellaceae}GGB1631 (SGB2239)","{f__Prevotellaceae}GGB1247 (SGB1668)","{f__Rikenellaceae}GGB1630 (SGB2238)"))),c(array("Disease associated",dim=c(10,1)),array("Healthy",dim=c(10,1))))
colnames(health_bar)<-c("effect","fdr","name","category")
health_bar<-as.data.frame(health_bar)

ggplot(health_bar,aes(x=reorder(as.character(name),as.numeric(as.character(effect))),y=as.numeric(as.character(effect)),fill=category))+
  geom_bar(stat="identity",color="black",size=1)+scale_fill_manual(values = c("#4575b4","#d73027"))+#+stat_bin(geom="text",aes(),vjust=-1)
  theme(text = element_text(size=20),panel.grid.major = element_blank(),panel.background = element_blank(),axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black"),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),legend.position = "top",legend.text = element_text(size=20)) +
  ylab("Effect size")+coord_flip()+xlab("")+mdthemes::md_theme_classic(axis.text.y = element_text(color="black",size=rel(20)),axis.text.x = element_text(color="black",size=20))

library(svglite)

h<-ggplot(health_bar,aes(x=reorder(as.character(name),as.numeric(as.character(effect))),y=as.numeric(as.character(effect)),fill=category))+
  geom_bar(stat="identity",color="black",size=1)+scale_fill_manual(values = c("#4575b4","#d73027"))+#+stat_bin(geom="text",aes(),vjust=-1)
  theme(text = element_text(size=20),panel.grid.major = element_blank(),panel.background = element_blank(),axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black"),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),legend.position = "top",legend.text = element_text(size=20)) +
  ylab("Effect size")+coord_flip()+xlab("")#+mdthemes::md_theme_classic(axis.text.y = element_text(color="black",size=rel(20)),axis.text.x = element_text(color="black",size=20))

library(svglite)
ggsave(file="Figure 2E 2.svg", plot=h, width=6, height=4)

svglite('Figure 2E3.svg', width=12, height=6)
plot(h)
dev.off()

##########################################################################################################################
##Figure 3

##Read in food and animal files
effect_size_phen_sig_plot<-read.csv('effect_sig_plot_food_an.csv',row.names=1)#Need to extract sp names separately
disp_fdr2<-read.csv('disp_fdr_food_an.csv',row.names=1)
phyl_ph<-read.csv('phyl_ph_food_an.csv',row.names = 1)

for(i in 1:dim(disp_fdr2)[1]){
  for(j in 1:dim(disp_fdr2)[2]){
    if(is.na(disp_fdr2[i,j])==T){
      disp_fdr2[i,j]<-""
    }
  }
}

paletteLength<-15
myColor <- colorRampPalette(c("#cb181d","white", "#33a02c"))(paletteLength)#myColor <- colorRampPalette(c("Orange", "white", "Darkblue"))(paletteLength)
myBreaks <- c(seq(min(effect_size_phen_sig_plot),0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(effect_size_phen_sig_plot)/paletteLength, max(effect_size_phen_sig_plot), length.out=floor(paletteLength/2)))

myBreaks
myBreaks[7]<-as.numeric(-0.001)
myBreaks[9]<-as.numeric(0.001)
myBreaks

pheatmap(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr2),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

#With phylum annotation

pheatmap(effect_size_phen_sig_plot,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr2),number_color = "black",treeheight_row = 0,cluster_rows = F,cluster_cols = T,treeheight_col = 0,annotation_col=phyl_ph2)



##Diet diversity

#library(vegan)
#mb_samp5<-read.csv('mb_samp_1188_3.csv',row.names = 1)
#mb_samp5<-read.csv('mb_samp_sp_788_rel.csv',row.names=1)
#mb_samp5<-t(mb_samp5)
#div_alpha<-array(0,dim=c(dim(mb_samp5)[2]))
#temp<-array(0,dim=c(dim(mb_samp_sp4)[2]))
#rm(diversity)

#for(i in 1:dim(mb_samp5)[2]){
#  div_alpha[i]<-diversity(mb_samp5[,i])#Default is Shannon!-- try thresholding? was [,i]-- miniscule difference
  #temp[i]<-diversity(mb_samp_sp4[,i],index="simpson")
#}
#div_alpha_ph5[,1]<-div_alpha

div_alpha_ph5<-read.csv('dds_diversity.csv',row.names=1)
div_alpha_ph5<-as.data.frame(div_alpha_ph5)
colnames(div_alpha_ph5)<-c("alpha_div","DDS")

#Chronic conditions

div_alpha_ph5$DDS <- factor(div_alpha_ph5$DDS)
library(ggplot2)
ggplot(div_alpha_ph5, aes(x=DDS, y=as.numeric(as.character(alpha_div)),group=DDS)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a50026","#d73027","#f46d43","#fdae61","#fee08b","#d9ef8b","#a6d96a","#66bd63","#1a9850","#006837"),outlier.alpha=0)+labs(x="",y="Alpha diversity")+
  theme(text = element_text(size=15))+ scale_x_discrete(breaks=seq(1,10,1))#+ylim(1.4,5.2)#,"#ff7f00","#cab2d6","#6a3d9a","#ffff99"#x="Relationship",

div_alpha_ph6<-div_alpha_ph5
for(i in 1:dim(div_alpha_ph5)[1]){
  if(as.numeric(as.character(div_alpha_ph5$DDS[i]))<=3){
    div_alpha_ph6$DDS[i]<-1
  }else if((as.numeric(as.character(div_alpha_ph5$DDS[i]))>3)&(as.numeric(as.character(div_alpha_ph5$DDS[i]))<=6)){
    div_alpha_ph6$DDS[i]<-2
  }else if(as.numeric(as.character(div_alpha_ph5$DDS[i]))>6){
    div_alpha_ph6$DDS[i]<-3
  }
}

#x1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$chronic=="Healthy(848)")]
#y1<-div_alpha_ph2$alpha_div[which(div_alpha_ph2$chronic=="Allergies(104)")]
#x1<-x1[which(is.na(x1)==F)]
#y1<-y1[which(is.na(y1)==F)]

wilcox.test(as.numeric(as.character(x1)),as.numeric(as.character(y1)))

library(ggpubr)
my_comparisons <- list( c("1", "4"), c("3", "4"))

#Save 4x5 pdf for Illana plot
ggplot(div_alpha_ph6, aes(x=DDS, y=as.numeric(as.character(alpha_div)),group=DDS)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a50026","#fee08b","#006837"),outlier.alpha=0)+labs(x="Diet diversity score",y="Shannon diversity")+
  theme(text = element_text(size=20),axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black"),panel.background=element_blank(),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))+ scale_x_discrete(breaks=seq(1,3,1))#+
 # stat_compare_means(comparisons = my_comparisons, label.y = c(5.2, 5.5),size=5)+ylim(2,6.7)#+

library(svglite)
svglite('Figure S dds2.svg', width=5, height=4)
ggplot(div_alpha_ph5, aes(x=DDS, y=as.numeric(as.character(alpha_div)),group=DDS)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a50026","#d73027","#f46d43","#fdae61","#fee08b","#d9ef8b","#a6d96a","#66bd63","#1a9850","#006837"),outlier.alpha=0)+labs(x="Diet diversity score",y="Shannon diversity")+
  theme(text = element_text(size=20),axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black"),panel.background=element_blank(),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))+ scale_x_discrete(breaks=seq(1,10,1))+
  stat_compare_means(comparisons = my_comparisons, label.y = c(5.2, 5.5),size=5)+ylim(2,6.7)#+xlim(-1,10)#+
dev.off()




##########################################################################################################################
## Figure 4

#Read in the socio-economic figures

#
paletteLength<-15
myColor <- colorRampPalette(c("#cb181d","white", "#33a02c"))(paletteLength)#myColor <- colorRampPalette(c("Orange", "white", "Darkblue"))(paletteLength)
myBreaks2 <- c(seq(min(effect_size_phen_sig_plot3),0, length.out=ceiling(paletteLength/2) + 1), 
               seq(max(effect_size_phen_sig_plot3)/paletteLength, max(effect_size_phen_sig_plot3), length.out=floor(paletteLength/2)))

pheatmap(effect_size_phen_sig_plot3,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot3),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks2,display_numbers = t(disp_fdr3),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

#With phylum
pheatmap(effect_size_phen_sig_plot3,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot3),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr3),number_color = "black",treeheight_row = 0,cluster_rows = F,cluster_cols = T,treeheight_col = 0,annotation_col=phyl_ph3)

#Household wealth index

div_alpha_ph4<-read.csv('div_alpha_ph2_wealth.csv',row.names=1)
library(ggplot2)
ggplot(div_alpha_ph4, aes(x=wealth, y=as.numeric(as.character(alpha_div)),group=wealth)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#b2df8a","#33a02c","#e31a1c","#fdbf6f"),outlier.alpha = 0)+labs(x="",y="Shannon diversity")+
  theme(text = element_text(size=15))+ scale_x_discrete(breaks=seq(1,5,1))+ylim(1.5,5.2)#,"#ff7f00","#cab2d6","#6a3d9a","#ffff99"#x="Relationship",

library(ggpubr)
my_comparisons <- list(c("Wealth index = 1", "Wealth index = 5"),c("Wealth index = 2", "Wealth index = 5"),c("Wealth index = 3", "Wealth index = 5"))
#bray_chronic2<-bray_chronic[which(is.na(bray_chronic[,1])==F),]

ggplot(div_alpha_ph4, aes(wealth, y=as.numeric(as.character(alpha_div)),group=wealth)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#b2df8a","#33a02c","#e31a1c","#fdbf6f"),outlier.alpha=0)+labs(x="Household wealth index",y="Shannon diversity")+
  theme(text = element_text(size=20),axis.text.y = element_text(color="black",size=20),axis.text.x = element_text(color="black",size=20),panel.background=element_blank(),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))+ scale_x_discrete(breaks=seq(1,5,1))+
  stat_compare_means(comparisons = my_comparisons, label.y = c(5.2,5.5,5.8),size=5)+ylim(1.8,6.2)#+xlab("Household welath index",vjust=0.5)#+scale_x_continuous(limits=c("1","2","3","4","5"))#+

g<-ggplot(div_alpha_ph4, aes(wealth, y=as.numeric(as.character(alpha_div)),group=wealth)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#b2df8a","#33a02c","#e31a1c","#fdbf6f"),outlier.alpha=0)+labs(x="Household wealth index",y="Shannon diversity")+
  theme(text = element_text(size=20),axis.text.y = element_text(color="black",size=20),axis.text.x = element_text(color="black",size=20),panel.background=element_blank(),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))+ scale_x_discrete(breaks=seq(1,5,1))+
  scale_x_discrete(labels=c(1:5))+
  stat_compare_means(comparisons = my_comparisons, label.y = c(5.2,5.5,5.8),size=5)+ylim(1.8,6.2)#+xlab("Household welath index",vjust=0.5)#+scale_x_continuous(limits=c("1","2","3","4","5"))#+

g

library(svglite)

ggsave(file="Figure 4B3.pdf", plot=g, width=5, height=4)




###################################################################################################################
## Figure 5

#5A


#5B -- spearman regression
e1<-read.csv('esz_strain_with_phylo_4.csv',row.names=1)
e2<-read.csv('esz_strain_without_phylo.csv',row.names=1)
f1<-read.csv('fdr_with_phylo_4.csv',row.names=1)
f2<-read.csv('fdr_without_phylo.csv',row.names=1)

e11<-c(0,0)
e12<-c(0,0)

for(i in 1:dim(f1)[1]){
  for(j in 1:dim(f1)[2]){
    if(is.na(f1[i,j])==F){
      if(f1[i,j]<0.05){
        #i1<-rbind(i1,i)
        #j1<-rbind(j1,j)
        e11<-rbind(e11,cbind(e2[i,j],e1[i,j]))
      }else{
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

library(ggpubr)
ggscatter(e11, x = "x1", y = "y1", 
          add = "reg.line", conf.int = FALSE, 
          cor.coef = TRUE, cor.method = "spearman",alpha=0.3,
          xlab = "Without phylogenic effect", ylab = "With phylogenic effect")+xlim(-1,1)+ylim(-1,1)#+xlim(-1.8,1.8)+ylim(-1.8,1.8)

ggscatter(e11, x = "x1", y = "y1", 
          add = "reg.line", conf.int = FALSE, 
          cor.coef = TRUE, cor.method = "spearman",alpha=0.3,
          xlab = "Without phylogenic effect", ylab = "With phylogenic effect")+geom_abline(slope=1,linetype=3)+xlim(-1,1)+ylim(-1,1)#+xlim(-1.8,1.8)+ylim(-1.8,1.8)

colnames(e12)<-c("x1","y1")
e13<-rbind(cbind(cbind(e11,"shp"=array(16,dim=c(dim(e11)[1],1))),"coll2"=array("red",dim=c(dim(e11)[1],1))),cbind(cbind(e12,"shp"=array(21,dim=c(dim(e12)[1],1))),"coll2"=array("black",dim=c(dim(e12)[1],1))))
#e14<-rbind(e11,e12)

op<-cor.test(e13$x1,e13$y1,method="spearman")
rho<-op$estimate
pval<-op$p.value
op2<-lm(e13$y1~e13$x1)
coef(op2)
e14<-e13[c(dim(e13)[1]:-1:1),]
#5,6 -- 5.5 x 6.00 now 
ggplot(e14,aes(x=x1,y=y1))+geom_point(alpha=0.1,color="#4575b4")+geom_abline(intercept=0,slope=1,linetype="dashed", size=1.3)+xlim(-1,1)+ylim(-1,1)+
  #geom_abline(intercept=coef(op2)[1],slope=0.95,col="red", size=1)
  geom_segment(x=1,y=1*0.95+coef(op2)[1],xend=-1,yend=-1*0.95+coef(op2)[1],col="red", size=1)+
  theme(legend.position = "none",panel.background = element_blank(),panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),axis.text.x=element_text(size=15,color="black"),axis.text.y=element_text(size=15,color="black"),axis.line.y = element_line(color="black"),axis.line.x = element_line(color="black"))+
  xlab("Without strain-phylogenetic effect")+ylab("With strain-phylogenetic effect")+
  theme(plot.margin=margin(c(1,0,0,0), unit = "cm"))


ggplot(e14,aes(x=x1,y=y1,color=factor(coll2)))+geom_point(aes(alpha=ifelse(shp==16, 0.6, 0.1)))+geom_abline(intercept=0,slope=1,linetype="dashed", size=1.3)+xlim(-1,1)+ylim(-1,1)+
  #geom_abline(intercept=coef(op2)[1],slope=0.95,col="red", size=1)
  geom_segment(x=1,y=1*0.95+coef(op2)[1],xend=-1,yend=-1*0.95+coef(op2)[1],col="red", size=1)+
  theme(legend.position = "none",panel.background = element_blank(),panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),axis.title.x=element_text(size=15),axis.title.y=element_text(size=15))+
  xlab("Without phylogenetic effect")+ylab("With phylogenetic effect")
  #+geom_point()


###########################################
# Strain vs no strain most distinguished species and phenotype

ml<-read.csv('esz_strain_with_phylo_4.csv',row.names=1)
mll<-read.csv('esz_strain_without_phylo.csv',row.names=1)
kl<-read.csv('fdr_with_phylo_4.csv',row.names=1)#Need to get the 126 phenotypes (only phenotypes we want), also include DDS
kll<-read.csv('fdr_without_phylo.csv',row.names=1)

st_name_full<-read.csv('strain_names_776.csv',row.names=1)


mlp<-ml
mllp<-mll
for(i in 1:dim(mll)[1]){
  for(j in 1:dim(mll)[2]){
    if((is.na(ml[i,j])==F)&(is.na(mll[i,j])==F)){
      sm<-i
    }else{
      mlp[i,j]<-NA
      mllp[i,j]<-NA
    }
  }
}

opp<-which(is.na(mlp[1,])==F)
lu<-wilcox.test(as.numeric(as.character(mlp[1,opp])),as.numeric(as.character(mllp[1,opp])))

##Wilcoxon test
opp_t<-array(NA,dim=c(dim(mlp)[1],1))

for(i in 1:dim(mllp)[1]){
  opp<-which(is.na(mllp[i,])==F)
  if(length(opp)>0){
    lu<-wilcox.test(as.numeric(as.character(mlp[i,opp])),as.numeric(as.character(mllp[i,opp])))
    opp_t[i,1]<-lu$p.value
  }
}

opp_t2<-opp_t[which(is.na(opp_t)==F),1]
opp_t2<-as.data.frame(opp_t2)
rownames(opp_t2)<-st_name_full[which(is.na(opp_t)==F),1]
opp_t2<-cbind(opp_t2,1:dim(opp_t2)[1])
opp_t3<-opp_t2[order(opp_t2$"opp_t2"),]

barplot(-log10(opp_t3[1:20,1]))

#
gh<-cbind(cbind(-log10(opp_t3[1:20,1]),rownames(opp_t3)[1:20]),c(1,1,1,2,2,2,array(3,dim=c(14))))
gh<-as.data.frame(gh)
colnames(gh)<-c("pval","species","cat1")

library(ggplot2)
ggplot(gh,aes(x=reorder(as.character(species),-1*as.numeric(as.character(pval))),y=as.numeric(as.character(pval)),fill=cat1))+
  geom_bar(stat="identity",color="black",size=1)+scale_fill_manual(values = c("#d73027","#2166ac","#f5f5f5"))+geom_hline(yintercept=1.3,linetype="dashed",size=1,color="blue")+geom_hline(yintercept=2,linetype="dashed",size=1,color="red")+#+stat_bin(geom="text",aes(),vjust=-1)
  theme(text = element_text(size=20),panel.grid.major = element_blank(),panel.background = element_blank(),axis.text.y = element_text(color='black'),axis.text.x = element_text(color='black')) +
  ylab("-log10(pvalue)")+xlab("")+coord_flip()#,angle=90,vjust=0.2
#save 5.00 x 10.00
#write.csv(opp_t3,'strain_with_without_phylo.csv')

##Deviation


#Do from phenotype angle

opp_p<-array(NA,dim=c(dim(mlp)[2],1))

for(i in 1:dim(mllp)[2]){
  opp<-which(is.na(mllp[,i])==F)
  if(length(opp)>0){
    lu<-wilcox.test(as.numeric(as.character(mlp[opp,i])),as.numeric(as.character(mllp[opp,i])))
    opp_p[i,1]<-lu$p.value
  }
}

opp_p2<-opp_p[which(is.na(opp_p)==F),1]
opp_p2<-as.data.frame(opp_p2)
rownames(opp_p2)<-colnames(mllp)[which(is.na(opp_p)==F)]
opp_p2<-cbind(opp_p2,1:dim(opp_p2)[1])
opp_p3<-opp_p2[order(opp_p2$"opp_p2"),]

barplot(-log10(opp_p3[1:20,1]))

#
gh<-cbind(cbind(-log10(opp_p3[1:20,1]),as.character(rownames(opp_p3)[1:20])),c(array(1,dim=c(2)),array(2,dim=c(9)),array(3,dim=c(9))))
gh<-as.data.frame(gh)
colnames(gh)<-c("pval","ph","cat1")
gh$ph<-as.character(gh$ph)
#gh$ph<-as.character(c("No stove","Vegetables","No cellphone","Unfinished windows","Cream/butter","Sheep","Eggs","Separate kitchen","Stove","Stomach illness","Anti-parasitic","Parakeet","Goat","Duck","None wild","Anti-hypertensive","Lizard","Pain killers","Anti-diarrheal","Beans"))
gh$ph<-as.character(c("Clustering coefficient","Transitivity (friendship)","Fish","Plastic roof","No cellphone","Anti-fungal","Anti-parasitic","Intestinal illness","Antibiotics","Allergies","Fruit juice","Rabbit","Ham/sausages/hotdog","Anti-diarrheal","Vegetables","Cough (1month)","Asthma","Transitivity (familial)","Chips","Stomach illness"))

library(ggplot2)
ggplot(gh,aes(x=reorder(as.character(ph),-1*as.numeric(as.character(pval))),y=as.numeric(as.character(pval)),fill=cat1))+
  geom_bar(stat="identity",color="black",size=1)+scale_fill_manual(values = c("#b2abd2","#4daf4a","#e41a1c"))+geom_hline(yintercept=1.3,linetype="dashed",size=1,color="blue")+geom_hline(yintercept=2,linetype="dashed",size=1,color="red")+geom_hline(yintercept=3,linetype="dashed",size=1,color="green")+geom_hline(yintercept=4,linetype="dashed",size=1,color="purple")+#+stat_bin(geom="text",aes(),vjust=-1)
  theme(text = element_text(size=20),panel.grid.major = element_blank(),panel.background = element_blank(),axis.text.y = element_text(color='black'),axis.text.x = element_text(color='black')) +
  ylab("-log10(pvalue)")+xlab("")+coord_flip()#,angle=90,vjust=0.2






#Refer figure5_strain_ph_assc.R

##############################################################################################################################################
##Supplementary figures

#Dangerous health associations

#Input dangerous health tables

#

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
pheatmap(effect_size_phen_sig_plot4,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot4),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr4),number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = T)

#Annotation

pheatmap(effect_size_phen_sig_plot4,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(effect_size_phen_sig_plot4),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "white",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,display_numbers = t(disp_fdr4),number_color = "black",treeheight_row = 0,cluster_rows = F,cluster_cols = T,treeheight_col = 0,annotation_col=phyl_ph4)


##Bray curtis with phenotypes

#CHronic conditions
bray_chronic<-read.csv('bray_chronic.csv',row.names = 1)
bray_chronic$chronic <- factor(bray_chronic$chronic, levels=c("Healthy","Diabetes","Allergies","Cystic fibrosis","Heart disease","Endocrine illness","Asthma","Stomach illness","Intestinal illness","Arthritis"),ordered=TRUE)


my_comparisons <- list( c("Healthy", "Diabetes"),c("Healthy", "Allergies"),c("Healthy", "Heart disease"),c("Healthy", "Asthma"), c("Healthy", "Stomach illness"), c("Healthy", "Intestinal illness"),c("Healthy", "Arthritis") )
#bray_chronic2<-bray_chronic[which(is.na(bray_chronic[,1])==F),]
ggplot(bray_chronic, aes(x=chronic, y=as.numeric(as.character(bray)),group=chronic)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#1f78b4","#b2df8a","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6"),outlier.alpha=0)+labs(x="",y="Bray-curtis dissimilarity(vs healthy)")+
  theme(text = element_text(size=15))+ scale_x_discrete(breaks=seq(1,8,1))+
  stat_compare_means(comparisons = my_comparisons, label.y = c(1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8))+ylim(0.3,1.9)#+

#Medication
bray_medicaiton<-read.csv('bray_medicaiton.csv',row.names = 1)

bray_medication$medication <- factor(bray_medication$medication, levels=c("No medication","Antibiotic","Anti-diarrheal","Anti-parasitic","Anti-fungal","Vitamins","Anti-hypertensive"),ordered=TRUE)

my_comparisons <- list( c("No medication", "Antibiotic"),c("No medication", "Anti-diarrheal"),c("No medication", "Anti-parasitic"),c("No medication", "Anti-fungal"),c("No medication", "Vitamins"),c("No medication", "Anti-hypertensive"))
#bray_chronic2<-bray_chronic[which(is.na(bray_chronic[,1])==F),]
ggplot(bray_medication, aes(x=medication, y=as.numeric(as.character(bray)),group=medication)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f"),outlier.alpha=0)+labs(x="",y="Bray-curtis(vs no medication)")+
  theme(text = element_text(size=15))+ scale_x_discrete(breaks=seq(1,7,1))+
  stat_compare_means(comparisons = my_comparisons, label.y = c(1.1,1.2,1.3,1.4,1.5,1.6))+ylim(0.3,1.7)#+


##PCOA plot

##PCOA plot
#load("~/Chris network papers/Honduras/Ecological niches project 01-18-22/community_detection_02-06-22/General physical all villages/DMP project/dmp_lmer_plotting_all_100922.R.RData")
st<-read.csv('bristol.csv',row.names = 1)
phen_all_use3<-read.csv('phen_all_use6_all.csv',row.names=1)
mb_samp_1188_3<-read.csv('mb_samp_1188_3.csv',row.names=1)

colnames(phen_all_use3)[15]
unique(phen_all_use3[,15])
unique(phen_all_use3[,15])
#write.csv(mb_samp_1188_3,'MB_1187_2285.csv')
#which(mb_samp_sp_name=="Prevotella_copri")

#klm<-rowSums(mb_samp_1188_3)
#ffm<-which(klm>1000)

which((grepl("SGB1626",rownames(mb_samp_1188_3),fixed=T))==T)
#SGB79883
rownames(mb_samp_1188_3)[337]

#for(i in 1:dim(data_transformed)[1]){
#cleaned_data<-phen_all_use[,c(6,4,224,225,40)]
#cleaned_data<-cbind(cleaned_data,"bristol"=st)
#cleaned_data<-cbind(cleaned_data,"a1c"=phen_all_use[,32])
sp<-as.data.frame(mb_samp_1188_3[which((grepl("SGB1626",rownames(mb_samp_1188_3),fixed=T))==T),])#SGB79883,SGB1626-prev clade A,1613,
#ffm[18]
#sp<-as.data.frame(mb_samp_1188_3[ffm[18],])
st<-st[,1]
#sp2<-sp[ind_u]
#cleaned_data<-cbind(cleaned_data,sp)
cleaned_data<-array(0,dim=c(dim(phen_all_use3)[1],8))
for(k in 1:dim(phen_all_use3)[1]){
  cleaned_data[k,1]<-as.numeric(phen_all_use3[k,6])
  cleaned_data[k,2]<-as.numeric(phen_all_use3[k,4])
  cleaned_data[k,3]<-as.numeric(phen_all_use3[k,224])
  cleaned_data[k,4]<-as.numeric(phen_all_use3[k,225])
  cleaned_data[k,5]<-as.numeric(phen_all_use3[k,40])
  cleaned_data[k,6]<-as.numeric(st[k])
  if(is.na(phen_all_use3[k,15])==F){
    if(phen_all_use3[k,15]=="Poor"){
      cleaned_data[k,7]<-0
    }else if(phen_all_use3[k,15]=="Fair"){
      cleaned_data[k,7]<-1
    }else if(phen_all_use3[k,15]=="Good"){
      cleaned_data[k,7]<-2
    }else if(phen_all_use3[k,15]=="Very good"){
      cleaned_data[k,7]<-3
    }else if(phen_all_use3[k,15]=="Excellent"){
      cleaned_data[k,7]<-4
    }
  }
  if(is.na(phen_all_use3[k,15])==T){
    cleaned_data[k,7]<-NA
  }
  cleaned_data[k,8]<-as.numeric(sp[k])
}
cleaned_data<-as.data.frame(cleaned_data)
colnames(cleaned_data)<-c("age","sex","batch_effect","dna_conc","BMI","bristol","Overall_health","sp")
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


#}
# calculate bray-curtis dissimilarity
library(vegan)
dist_sp <- vegdist(t(mb_samp_1188_3[,ind_temp]),  method = "bray")
# do PCOA
#myPCOA <- pcoa(dist_sp)
myPCOAv2= cmdscale(dist_sp)
#biplot(myPCOA)
#dev.off()
# fit phenotype vectors on the PCOA data
en2=envfit(myPCOAv2,cleaned_data2[,c(1:2,5:7)], na.rm=T, permutations=9999)
en_coord_cont = as.data.frame(scores(en2, "vectors")) * ordiArrowMul(en2)
en_coord_sp = as.data.frame(scores(en2,display="vectors"))#, display = "species"
coord2=en_coord_cont/2
row.names(coord2)=c("AGE","SEX","BMI","BRISTOL","HEALTH STATUS")

# plot PCOA with phenotype vectors
library(viridis)
#rownames(mb_samp_1188_3)[427]
#`Prevotella Copri(Clade D)`<-t(cleaned_data2[,8])
#ffm[18]
`Prevotella copri clade A(SGB1626)`<-t(cleaned_data2[,8])
library(ggplot2)
ggplot(data = as.data.frame(myPCOAv2), aes(x = V1, y = V2)) + 
  geom_point(data = as.data.frame(myPCOAv2), size = 2, alpha = 0.6, aes(colour=`Prevotella copri clade A(SGB1626)`))  +
  scale_color_viridis(option = "D",direction=-1)+
  geom_segment(aes(x = 0, y = 0, xend = Dim1, yend = Dim2), data = coord2, size =1, colour = "black") +
  geom_text(data = coord2, aes(x = Dim1, y = Dim2), colour = "black", fontface = "bold", label = row.names(coord2)) + 
  theme_bw( ) +xlab ("PCoA1") + ylab ("PCoA2")  #+ xlim(-0.4,0.6) +ylim(-0.4,0.4)


##Do for every sub question?
PHQ9_score<-phen_all_use3$phq9_score[ind_temp]

ggplot(data = as.data.frame(myPCOAv2), aes(x = V1, y = V2)) + 
  geom_point(data = as.data.frame(myPCOAv2), size = 2, alpha = 0.6, aes(colour=PHQ9_score))  +
  scale_color_viridis(option = "D",direction=-1)+
  #geom_segment(aes(x = 0, y = 0, xend = Dim1, yend = Dim2), data = coord2, size =1, colour = "black") +
  #geom_text(data = coord2, aes(x = Dim1, y = Dim2), colour = "black", fontface = "bold", label = row.names(coord2)) + 
  theme_bw( ) +xlab ("PCoA1") + ylab ("PCoA2")  #+ xlim(-0.4,0.6) +ylim(-0.4,0.4)

######################################################################################
## Figure plotting supplementary
#
#Correlation plots

library(corrplot)

phen_all_use<-read.csv('phen_all_use6_all.csv',row.names=1)
ind_list<-c(40,32,35,36,71:74,37)

for(i in 1:dim(phen_all_use)[1]){
  for(j in 1:length(ind_list)){
    if(is.na(phen_all_use[i,ind_list[j]])==F){
      if(phen_all_use[i,ind_list[j]]=="Dont_Know"){
        phen_all_use[i,ind_list[j]]<-NA
      }
    }
  }
}

phen_map<-as.numeric(as.character(phen_all_use$mb_d0400))+(as.numeric(as.character(phen_all_use$mb_d0300))-as.numeric(as.character(phen_all_use$mb_d0400)))/3
phen_all_use<-cbind(phen_all_use,"MAP"=phen_map)

phys<-array(NA,dim=c(dim(phen_all_use)[1],6))
ind_list2<-c(32,383,40,71,74,37)

for(i in 1:dim(phys)[1]){
  for(j in 1:length(ind_list2)){
    if(is.na(phen_all_use[i,ind_list2[j]])==F){
      phys[i,j]<-as.numeric(as.character(phen_all_use[i,ind_list2[j]]))
    }
  }
}

colnames(phys)<-c("Hb A1c","MAP","BMI","Heart rate","Oxygen saturation","Hb total")

phys2<-array(NA,dim=c(dim(phen_all_use)[1],14))
ind_list3<-c(20,21,23,26:29,136:142)
for(i in 1:dim(phys2)[1]){
  for(j in 1:length(ind_list3)){
    if(is.na(phen_all_use[i,ind_list3[j]])==F){
      phys2[i,j]<-1
    }else if(is.na(phen_all_use[i,ind_list3[j]])==T){
      phys2[i,j]<-0
    }
  }
}

colnames(phys2)<-c("Diabetes (N=24)","Allergies (N=104)","Heart disease (N=47)","Asthma (N=47)","Stomach illness (N=175)","Intestinal illness (N=70)","Arthritis (N=39)","Pain killers (N=714)","Antibiotics (N=137)","Anti-diarrheal (N=49)","Anti-parasitic (N=34)","Anti-fungal (N=82)","Vitamins (N=121)","Anti-hypertensive (N=72)")

phys3<-cbind(phys,phys2)
phys3<-as.data.frame(phys3)

ind_na<-0
for(i in 1:dim(phys3)[1]){
if(length(which(is.na(phys3[i,])==T))==0){
  ind_na<-rbind(ind_na,i)
}
}

phys4<-phys3[ind_na,]

##Correlating effect sizes --mental health plot
#colnames(effect_size_phen_lmer)[c(12,14)]<-c("Excellent(health)[75]","Poor(health)[152]")
#cp1<-effect_size_phen_lmer[which(fdr_chk_lmer_all>0),c(1,4,5,6,9:10,14,12,11,16:17,19,22:25,27:34,42:43,164,166,152:154,156:158)]#,86,85-- cognitive impairment, dementia
#colnames(cp1)<-colnames(effect_size_phen_lmer[c(1,4,5,6,9:10,14,12,11,16:17,19,22:25,27:34,42:43,164,166,152:154,156:158)])
#cp_ind<-match(cp1_names,mb_samp_sp_name)
#cp_cor<-cor(cp1[cp_ind,],cp1[cp_ind,],method="pearson")
cp_cor<-cor(phys4,phys4,method="pearson")
pheatmap(cp_cor,treeheight_col = 0)

library(corrplot)

corrplot(as.matrix(cp_cor),method = "square",tl.cex = 0.85,cl.cex = 0.85,tl.col = "black",order= "hclust")

cp_cor2<-ifelse(cp_cor>0.999,NA,cp_cor)

paletteLength<-15
myColor <- colorRampPalette(c("#cb181d","white", "#33a02c"))(paletteLength)#myColor <- colorRampPalette(c("Orange", "white", "Darkblue"))(paletteLength)
myBreaks <- c(seq(min(cp_cor),0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(cp_cor2[which(is.na(cp_cor2)==F)])/paletteLength, max(cp_cor2[which(is.na(cp_cor2)==F)]), length.out=floor(paletteLength/2)))

myBreaks
myBreaks[7]<-as.numeric(-0.001)
myBreaks[9]<-as.numeric(0.001)
myBreaks

#myBreaks<-c(myBreaks,1)
#myColor<-c(myColor,"#525252")


pheatmap(cp_cor2,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(cp_cor),label_col=colnames(cp_cor),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "grey50",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = F)

#Animals



ind_list4<-c(1:3,7:15,5,4,16,6,17:18,20:24,26)
il4<-c("Cat (N=645)","Dog (N=944)","Parakeet (N=102)","Rabbit (N=120)","Horse (N=314)","Mice (N=585)","None pet (N=88)","Cow (N=397)","Goat (N=40)","Pig (N=208)","Chicken (N=1094)","Duck (N=478)","Turkey (N=322)","Sheep (N=33)","Geese (N=68)","None farm (N=69)","Bat (N=410)","Lizard (N=423)","Monkey (N=12)","Snake (N=437)","Bird (N=794)","Possum (N=467)","Rat (N=457)","Squirrel (N=480)","Other Wild (N=1)","None wild (N=243)")
an2<-array(NA,dim=c(dim(phen_all_use)[1],length(ind_list4)))
for(i in 1:dim(an2)[1]){
  for(j in 1:length(ind_list4)){
    if(is.na(phen_all_use[i,154+ind_list4[j]])==F){
      an2[i,j]<-1
    }else if(is.na(phen_all_use[i,154+ind_list4[j]])==T){
      an2[i,j]<-0
    }
  }
}

an2<-as.data.frame(an2)

ind_na<-0
for(i in 1:dim(an2)[1]){
  if(length(which(is.na(an2[i,])==T))==0){
    ind_na<-rbind(ind_na,i)
  }
}

an3<-an2[ind_na,]
colnames(an3)<-il4[ind_list4]

cp_cor<-cor(an3,an3,method="pearson")
pheatmap(cp_cor,treeheight_col = 0)

library(corrplot)

corrplot(as.matrix(cp_cor),method = "square",tl.cex = 0.85,cl.cex = 0.85,tl.col = "black",order= "hclust")

cp_cor2<-ifelse(cp_cor>0.999,NA,cp_cor)

paletteLength<-15
myColor <- colorRampPalette(c("#cb181d","white", "#33a02c"))(paletteLength)#myColor <- colorRampPalette(c("Orange", "white", "Darkblue"))(paletteLength)
myBreaks <- c(seq(min(cp_cor),0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(cp_cor2[which(is.na(cp_cor2)==F)])/paletteLength, max(cp_cor2[which(is.na(cp_cor2)==F)]), length.out=floor(paletteLength/2)))

myBreaks
myBreaks[7]<-as.numeric(-0.001)
myBreaks[9]<-as.numeric(0.001)
myBreaks

#myBreaks<-c(myBreaks,1)
#myColor<-c(myColor,"#525252")


pheatmap(cp_cor2,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(cp_cor),label_col=colnames(cp_cor),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "grey50",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = F)

# Food
ind_list6<-c(75:90,382)
il4<-c("Beans","Tortillas","Rice","Bread","Milk","Yogurt","Cream/butter","Cheese","Eggs","Vegetables","Fruits","Natural juice","Chicken","Beef/Pork","Ham/sausages/hotdog","Fish","Diet diversity score")
an2<-array(NA,dim=c(dim(phen_all_use)[1],length(ind_list6)))

for(i in 1:dim(an2)[1]){
  for(j in 1:(length(ind_list6)-1)){
    if(is.na(phen_all_use[i,ind_list6[j]])==F){
      if(phen_all_use[i,ind_list6[j]]=="Never/rarely"){#Change to 1/50,1/10,4/7,1
        an2[i,j]<-1/50
      }else if((phen_all_use[i,ind_list6[j]]=="A few days per month")|(phen_all_use[i,ind_list6[j]]=="A few times per month")){
        an2[i,j]<-1/10
      }else if(phen_all_use[i,ind_list6[j]]=="A few days per week"){
        an2[i,j]<-4/7
      }else if(phen_all_use[i,ind_list6[j]]=="Every day"){
        an2[i,j]<-1
      }
    }
  }
}
an2[,17]<-phen_all_use[,382]

an2<-as.data.frame(an2)

ind_na<-0
for(i in 1:dim(an2)[1]){
  if(length(which(is.na(an2[i,])==T))==0){
    ind_na<-rbind(ind_na,i)
  }
}

an3<-an2[ind_na,]
colnames(an3)<-il4
an3[,17]<-an3[,17]/9

cp_cor<-cor(an3,an3,method="pearson")
pheatmap(cp_cor,treeheight_col = 0)

library(corrplot)

corrplot(as.matrix(cp_cor),method = "square",tl.cex = 0.85,cl.cex = 0.85,tl.col = "black",order= "hclust")

cp_cor2<-ifelse(cp_cor>0.999,NA,cp_cor)

paletteLength<-15
myColor <- colorRampPalette(c("#cb181d","white", "#33a02c"))(paletteLength)#myColor <- colorRampPalette(c("Orange", "white", "Darkblue"))(paletteLength)
myBreaks <- c(seq(min(cp_cor),0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(cp_cor2[which(is.na(cp_cor2)==F)])/paletteLength, max(cp_cor2[which(is.na(cp_cor2)==F)]), length.out=floor(paletteLength/2)))

myBreaks
myBreaks[7]<-as.numeric(-0.001)
myBreaks[9]<-as.numeric(0.001)
myBreaks

#myBreaks<-c(myBreaks,1)
#myColor<-c(myColor,"#525252")


pheatmap(cp_cor2,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(cp_cor),label_col=colnames(cp_cor),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "grey50",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = F)




#Economic factors



ind_list5<-c(310:316,306:309,317:359,375,376)
il5<-c("Electricity (N=1069)","Radio (N=519)","TV (N=565)","Cellphone (N=927)","No cellphone (N=12)","Refrigerator (N=397)","None electronics (N=38)","Chimney (N=959)","No chimney (N=151)","Stove (N=61)","No stove (N=1)","Wood (N=1109)","Gas(fuel) [N=46]","Electricity fuel (N=14)","Kerosene (N=2)","None fuel (N=1)","Separate kitchen (N=1106)","Cement floor (N=552)","Earth/Sand floor (N=507)","Ceramic floor (N=88)","Tiles floor (N=19)","Mud bricks floor (N=1)","Wood floor (N=0)","Other floor (N=5)","Wooden windows (N=1035)","Glass windows (N=51)","Metal windows (N=6)","Unfinished windows (N=35)","No windows (N=45)","Clay/mud walls(N=693)","Clay brick walls (N=3)","Cement walls (N=452)","Cane/palm/trunks walls (N=0)","Wood unpolished walls (N=21)","Wood polished walls (N=0)","Discarded materials walls (N=1)","No walls (N=0)","Other walls (N=2)","Plastic roof (N=7)","Metal roof (N=982)","Clay roof (N=120)","Thatch/palm roof (N=10)","Concrete roof (N=27)","Wood roof (N=22)","Other roof(N=4)","Sleeping rooms","Spring(protected) [N=917]","Spring(unprotected) [N=19]","Tube well (N=72)","Dug well(protected) [N=78]","Dug well(unprotected) [N=26]","Surface water (N=10)","Bottled water (N=15)","Other water (N=6)","Household size","Household wealth index")
ind_list5_c<-c(96,97,47,48,51,54,58,60,67,71,73,81,87:88,91,92,94)
ind_list5_c<-ind_list5_c-41
ec1<-array(NA,dim=c(dim(phen_all_use)[1],length(ind_list5)))
for(i in 1:dim(ec1)[1]){
  for(j in 1:7){
    if(is.na(phen_all_use[i,ind_list5[j]])==F){
      ec1[i,j]<-1
    }else if(is.na(phen_all_use[i,ind_list5[j]])==T){
      ec1[i,j]<-0
    }
  }
}

for(i in 1:dim(ec1)[1]){
  for(j in 8:length(ind_list5)){
    if(is.na(phen_all_use[i,ind_list5[j]])==F){
      ec1[i,j]<-as.numeric(as.character(phen_all_use[i,ind_list5[j]]))
    }
  }
}

ec1<-as.data.frame(ec1[,ind_list5_c])

ind_na<-0
for(i in 1:dim(ec1)[1]){
  if(length(which(is.na(ec1[i,])==T))==0){
    ind_na<-rbind(ind_na,i)
  }
}

ec2<-ec1[ind_na,]
colnames(ec2)<-il5[ind_list5_c]

cp_cor<-cor(ec2,ec2,method="pearson")
for(i in 1:dim(cp_cor)[1]){
  for(j in 1:dim(cp_cor)[2]){
    if(i==j){
      cp_cor[i,j]<-NA
    }
  }
}

pheatmap(cp_cor,treeheight_col = 0)

library(corrplot)
#Save 5.00 x 6.00
corrplot(as.matrix(cp_cor),method = "square",tl.cex = 0.85,cl.cex = 0.85,tl.col = "black",order= "hclust")


cp_cor2<-cp_cor

paletteLength<-15
myColor <- colorRampPalette(c("#cb181d","white", "#33a02c"))(paletteLength)#myColor <- colorRampPalette(c("Orange", "white", "Darkblue"))(paletteLength)
myBreaks <- c(seq(min(cp_cor2[which(is.na(cp_cor2)==F)]),0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(cp_cor2[which(is.na(cp_cor2)==F)])/paletteLength, max(cp_cor2[which(is.na(cp_cor2)==F)]), length.out=floor(paletteLength/2)))

myBreaks
myBreaks[7]<-as.numeric(-0.001)
myBreaks[9]<-as.numeric(0.001)
myBreaks

#myBreaks<-c(myBreaks,1)
#myColor<-c(myColor,"#525252")


pheatmap(cp_cor2,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(cp_cor),label_col=colnames(cp_cor),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "grey50",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = F)


########
## Hierarchichal clustering of Effect sizes of all phenotypes

jh<-read.csv('health_eff.csv',row.names = 1)
colnames(jh)<-c("Hb A1c","MAP","BMI","Heart rate","Oxygen saturation","Hb total","Poor (N=153)","Fair (N=612)","Very good (N=56)","Excellent (N=75)","Cough (N=257)","Diarrhea (N=73)","Bristol stool scale","Diabetes (N=24)","Allergies (N=104)","Heart disease (N=47)","Asthma (N=47)","Stomach illness (N=175)","Intestinal illness (N=70)","Arthritis (N=39)","Painkillers (N=714)","Antibiotics (N=137)","Anti-diarrheal (N=49)","Anti-parasitic (N=34)","Anti-fungal (N=82)","Vitamins (N=121)","Anti-hypertensive (N=72)","Reserved","Nervous","Openess","Cognitive impairment (N=56)","Dementia (N=89)","Alcohol daily frequency","Cigarette usage (N=46)","Cigarette frequency","Mild (GAD7)","Moderate (GAD7)","Severe (GAD7)","Mild (PHQ9)","Moderate (PHQ9)","Severe (PHQ9)")

jh2<-read.csv('food_an_eff.csv',row.names=1)
colnames(jh2)<-c("Cat (N=645)","Dog (N=944)","Parakeet (N=102)","None pet (N=88)","Cow (N=397)","Goat (N=40)","Pig (N=208)","Chicken (N=1094)","Duck (N=478)","Turkey (N=322)","Sheep (N=33)","Geese (N=68)","Horse (N=314)","Rabbit (N=120)","None farm (N=69)","Mice (N=585)","Bat (N=410)","Lizard (N=423)","Snake (N=437)","Bird (N=794)","Possum (N=467)","Rat (N=457)","Squirrel (N=480)","None wild (N=243)","Beans","Tortillas","Rice","Bread","Milk","Yogurt","Cream/butter","Cheese","Eggs","Vegetables","Fruits","Natural juice","Chicken","Beef/pork","Ham/sausages/hotdog","Fish")

jh3<-read.csv('socio_eco_eff.csv',row.names=1)
colnames(jh3)<-c("Partner live duration","Number of partners","Partner live age","Living with partner (N=375)","Grades 1-3 (N=373)","Grades 4-6 (N=311)","Grades >6 (N=113)","Friend ties(same building)","Friend ties(different building)","Betweeness(friendship)","Transitivity(friendship)","Family ties(same building)","Family ties(different building)","Transitivity(familial)","Betweeness(familial)","Degree(all ties)","Clustering coefficient(all ties)","Betweeness(all ties)","Kin percentage (upto third degree)","Altruism","Risk taking","Washing hands (N=1130)","Distance to village center","Distance to main road","Distance to health center","Deforestation (%)","Travel","Monthly expenditure","Household size","Household wealth index","Refrigerator (N=397)","None electronics (N=38)","Stove (N=61)","Gas(fuel) [N=46]","Separate kitchen (N=1106)","Earth/Sand floor (N=507)","Glass windows (N=51)","Clay/mud walls(N=693)","Cement walls (N=452)","Metal roof (N=982)","Sleeping rooms","Spring(protected) [N=917]","Dug well(protected) [N=78]","Dug well(unprotected) [N=26]","Bottled water (N=15)")

jh4<-cbind(cbind(jh,jh2),jh3)
write.csv(jh4,'species-phenotype-assc.csv')

jk4<-dist(t(jh4),method="euclidean")
jk5<-hclust(jk4)
plot(jk5)

#nodePar2 <- list(lab.cex = 0.6, pch = c(NA, 19), 
#                cex = 0.7, col = "blue")#Doesnt work
plot(jk5)

library(ggdendro)
ggdendrogram(jk5,rotate=T)
#plot(jk5,set("branches_k_color", k=3))


plot(as.phylo(jk5),type='fan')
plot(as.phylo(jk5),type='cladogram')
plot(as.phylo(jk5),type='unrooted')



library(ape)
#temp<-rtree(126)
temp<-as.dendrogram(jk5)
plot(temp)
#temp2<-read.tree(jk4)

##Actual plot --- evol_tree_plot in Andromache


######################
##Circular plot (Unhealthy ranges of physiological measurements)

phy_var3<-read.csv('phy_var3.csv',row.names=1)

phys<-c("Body mass index (BMI)","Heart rate","Hemoglobin A1c","Hemoglobin total (female)","Hemoglobin total (male)","Mean arterial pressure (MAP)","Oxygen saturation")

library(ggplot2)
ggplot(phy_var3, aes(x=phen_use, y=as.numeric(as.character(Values)),group=phen_use)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f"))+labs(y="Anthropometrics (normalized)",x="")+
  theme(text = element_text(size=15))+ scale_x_discrete(breaks=seq(1,7,1))+
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.background = element_blank(),axis.text.y = element_text(color="black",size=15),
        legend.position="left")+
  geom_hline(yintercept=1, linetype='dashed', col = 'red')+
  guides(fill = guide_legend(title = "Physiological"))+
  #theme(legend.position = c(.2,.85))
  geom_rect(data=phy_var3, 
            mapping=aes(xmin=x_norm, xmax=x_norm2, ymin=y_norm, ymax=y_norm2), fill = alpha("grey60", 0.008),color=alpha("red", .5))+
  #theme(plot.margin=margin(c(0.5,0.5,3,0.5), unit = "cm"))+
  scale_x_discrete(labels=phys)+
  theme(axis.text.x = element_text(angle = 60,color="black",vjust=(1),hjust=0.95))+
  theme(plot.margin=margin(c(0.5,0.5,0,0.5), unit = "cm"))

#Save 6.00 x 6.00






########################################################################################################################
##Creating tables

phen_all_use<-read.csv('phen_all_use6_all.csv',row.names=1)
ind_list<-c(40,32,35,36,71:74,37)

for(i in 1:dim(phen_all_use)[1]){
  for(j in 1:length(ind_list)){
    if(is.na(phen_all_use[i,ind_list[j]])==F){
      if(phen_all_use[i,ind_list[j]]=="Dont_Know"){
        phen_all_use[i,ind_list[j]]<-NA
      }
    }
  }
}

phen_map<-as.numeric(as.character(phen_all_use$mb_d0400))+(as.numeric(as.character(phen_all_use$mb_d0300))-as.numeric(as.character(phen_all_use$mb_d0400)))/3
phen_all_use<-cbind(phen_all_use,"MAP"=phen_map)

phys<-array(NA,dim=c(dim(phen_all_use)[1],6))
ind_list2<-c(32,383,40,71,74,37)

for(i in 1:dim(phys)[1]){
  for(j in 1:length(ind_list2)){
    if(is.na(phen_all_use[i,ind_list2[j]])==F){
      phys[i,j]<-as.numeric(as.character(phen_all_use[i,ind_list2[j]]))
    }
  }
}

colnames(phys)<-c("Hb A1c","MAP","BMI","Heart rate","Oxygen saturation","Hb total")

library(tableone)

tab1<-CreateTableOne(data = as.data.frame(phys))
print(tab1)
write.csv(print(tab1), file = "Supp_table1.csv")


phys2<-cbind(phys,phen_all_use$gender)
colnames(phys2)[7]<-"Gender"
phys2[,7]<-ifelse(phys2[,7]==1,"Female","Male")

phys3<-as.data.frame(phys2)
phys3$Gender<-factor(phys3$Gender,levels=c("Male","Female"))

phys4<-phys[which(phys3$Gender=="Male"),]
tab2<-CreateTableOne(data = as.data.frame(phys4))
#write.csv(print(tab2), file = "Supp_table1m.csv")

phys4<-phys[which(phys3$Gender=="Female"),]
tab2<-CreateTableOne(data = as.data.frame(phys4))
#write.csv(print(tab2), file = "Supp_table1f.csv")


health1<-array(NA,dim=c(dim(phen_all_use)[1],14))
ind_list3<-c(20,21,23,26:29,136:142)
for(i in 1:dim(health1)[1]){
  for(j in 1:length(ind_list3)){
    if(is.na(phen_all_use[i,ind_list3[j]])==F){
      health1[i,j]<-1
    }else if(is.na(phen_all_use[i,ind_list3[j]])==T){
      health1[i,j]<-0
    }
  }
}

colnames(health1)<-c("Diabetes (N=24)","Allergies (N=104)","Heart disease (N=47)","Asthma (N=47)","Stomach illness (N=175)","Intestinal illness (N=70)","Arthritis (N=39)","Pain killers (N=714)","Antibiotics (N=137)","Anti-diarrheal (N=49)","Anti-parasitic (N=34)","Anti-fungal (N=82)","Vitamins (N=121)","Anti-hypertensive (N=72)")
p2<-CreateTableOne(data = as.data.frame(health1),factorVars = colnames(health1))
#write.csv(print(p2),file="Supp_table_med_ch.csv")

CreateTableOne(data = as.data.frame(health1[which(phys3$Gender=="Male"),]),factorVars = colnames(health1))

CreateTableOne(data = as.data.frame(health1[which(phys3$Gender=="Female"),]),factorVars = colnames(health1))

health2<-phen_all_use[,15]

length(which(health2[which(phys3$Gender=="Male")]=="Excellent"))
length(which(health2[which(phys3$Gender=="Female")]=="Excellent"))

length(which(health2[which(phys3$Gender=="Male")]=="Poor"))
length(which(health2[which(phys3$Gender=="Female")]=="Poor"))

##Acute

st<-read.csv('bristol.csv',row.names = 1)
st<-as.matrix.data.frame(st)

acute<-array(NA,dim=c(dim(phen_all_use)[1],2))
ind_list2<-c(16,381)

for(i in 1:dim(acute)[1]){
  for(j in 1:length(ind_list2)){
    if(is.na(phen_all_use[i,ind_list2[j]])==F){
      if(phen_all_use[i,ind_list2[j]]=="Yes"){
      acute[i,j]<-1
    }else if(phen_all_use[i,ind_list2[j]]=="No"){
      acute[i,j]<-0
    }
    }
  }
}

colnames(acute)<-c("Cough","Diarrhea")

acute<-cbind(acute,"Bristol"=st)
colnames(acute)[3]<-"Bristol"

CreateTableOne(data = as.data.frame(acute),factorVars = colnames(acute))
CreateTableOne(data = as.data.frame(acute[which(phys3$Gender=="Male"),]),factorVars = colnames(acute))
CreateTableOne(data = as.data.frame(acute[which(phys3$Gender=="Female"),]),factorVars = colnames(acute))

##Personality

#person<-array(NA,dim=c(dim(phen_all_use)[1],2))
person<-phen_all_use[,c(143,151,152)]


p1<-CreateTableOne(data = as.data.frame(person),factorVars = colnames(person))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(person[which(phys3$Gender=="Male"),]),factorVars = colnames(person))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(person[which(phys3$Gender=="Female"),]),factorVars = colnames(person))
write.csv(print(p1),'p1.csv')

##Anxiety, depression

#person<-array(NA,dim=c(dim(phen_all_use)[1],2))
anx<-phen_all_use[,c(191,202)]


p1<-CreateTableOne(data = as.data.frame(anx),factorVars = colnames(anx))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(anx[which(phys3$Gender=="Male"),]),factorVars = colnames(anx))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(anx[which(phys3$Gender=="Female"),]),factorVars = colnames(anx))
write.csv(print(p1),'p1.csv')


##Alcohol and cigarettes
phen_all_rct_use<-read.csv('phen_all_rct_use22.csv',row.names=1)
alc<-array(NA,dim=c(1187,1))
kl<-8
for(k in 1:1187){
if(is.na(phen_all_rct_use[k,kl])==F){
  if(phen_all_rct_use[k,kl]=="1 or 2"){
    alc[k,1]<-1.5
  }else if(phen_all_rct_use[k,kl]=="3 or 4"){
    alc[k,1]<-3.5
  }else if(phen_all_rct_use[k,kl]=="5 or 6"){
    alc[k,1]<-5.5
  }else if(phen_all_rct_use[k,kl]=="7 to 9"){
    alc[k,1]<-8
  }else if(phen_all_rct_use[k,kl]=="10 or more"){
    alc[k,1]<-12
  }
}
}


p1<-CreateTableOne(data = as.data.frame(alc))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(alc[which(phys3$Gender=="Male"),]))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(alc[which(phys3$Gender=="Female"),]))
write.csv(print(p1),'p1.csv')

#Cigarette use

#person<-array(NA,dim=c(dim(phen_all_use)[1],2))
cig<-phen_all_rct_use[,10]


p1<-CreateTableOne(data = as.data.frame(cig),factorVars = colnames(cig))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(cig[which(phys3$Gender=="Male")]),factorVars = colnames(cig))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(cig[which(phys3$Gender=="Female")]),factorVars = colnames(cig))
write.csv(print(p1),'p1.csv')

##Cigarette frequency
cig<-phen_all_rct_use[,11]


p1<-CreateTableOne(data = as.data.frame(cig))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(cig[which(phys3$Gender=="Male")]))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(cig[which(phys3$Gender=="Female")]))
write.csv(print(p1),'p1.csv')

##Food and pets

food<-phen_all_use[,75:90]
#kl<-155:180

colnames(food)<-c("Beans","Tortillas","Rice","Bread","Milk","Yogurt","Cream/butter","Cheese","Eggs","Vegetables","Fruits","Natural juice","Chicken","Beef/pork","Ham/sausages/hotdog","Fish")

p1<-CreateTableOne(data = as.data.frame(food),factorVars = colnames(food))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(food[which(phys3$Gender=="Male"),]),factorVars = colnames(food))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(food[which(phys3$Gender=="Female"),]),factorVars = colnames(food))
write.csv(print(p1),'p1.csv')

#Animals
ind_list4<-c(1:3,7:15,5,4,16,6,17:18,20:24,26)
ann<-array(NA,dim=c(1187,length(c(155:180))))
#kl<-155:180
  for(k in 1:1187){
     for(kl in 155:180){
       if(is.na(phen_all_use[k,kl])==F){
         ann[k,(kl-154)]<-1
         }else if(is.na(phen_all_use[k,kl])==T){
           ann[k,(kl-154)]<-0
         }
         }
}
colnames(ann)<-c("Cat (N=645)","Dog (N=944)","Parakeet (N=102)","Rabbit (N=120)","Horse (N=314)","Mice (N=585)","None pet (N=88)","Cow (N=397)","Goat (N=40)","Pig (N=208)","Chicken (N=1094)","Duck (N=478)","Turkey (N=322)","Sheep (N=33)","Geese (N=68)","None farm (N=69)","Bat (N=410)","Lizard (N=423)","Monkey (N=12)","Snake (N=437)","Bird (N=794)","Possum (N=467)","Rat (N=457)","Squirrel (N=480)","Other Wild (N=1)","None wild (N=243)")
ann2<-ann[,ind_list4]
p1<-CreateTableOne(data = as.data.frame(ann2),factorVars = colnames(ann2))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(ann2[which(phys3$Gender=="Male"),]),factorVars = colnames(ann2))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(ann2[which(phys3$Gender=="Female"),]),factorVars = colnames(ann2))
write.csv(print(p1),'p1.csv')

#Partner
part<-array(NA,dim=c(1187,length(228:230)))
for(i in 1:1187){
  for(j in 1:dim(part)[2]){
    if((is.na(phen_all_use[i,(j+227)])==F)){
       if((phen_all_use[i,(j+227)]!="Dont_Know")|(phen_all_use[i,(j+227)]!="Removed")){
      part[i,j]<-as.numeric(as.character(phen_all_use[i,(j+227)]))
    }else{
      part[i,j]<-NA
    }
    }
  }
}


p1<-CreateTableOne(data = as.data.frame(part))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part[which(phys3$Gender=="Male"),]))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part[which(phys3$Gender=="Female"),]))
write.csv(print(p1),'p1.csv')

part<-phen_all_use[,227]
p1<-CreateTableOne(data = as.data.frame(part))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part[which(phys3$Gender=="Male")]))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part[which(phys3$Gender=="Female")]))
write.csv(print(p1),'p1.csv')

#Education

ed<-phen_all_rct_use[,5]
p1<-CreateTableOne(data = as.data.frame(ed))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(ed[which(phys3$Gender=="Male")]),factorVars = colnames(ed))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(ed[which(phys3$Gender=="Female")]),factorVars = colnames(ed))
write.csv(print(p1),'p1.csv')

#Social network factors
#203:207,360:363,365:369,228:229
part<-phen_all_use[,c(203:207,360:363,365:369)]
p1<-CreateTableOne(data = as.data.frame(part))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part[which(phys3$Gender=="Male"),]))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part[which(phys3$Gender=="Female"),]))
write.csv(print(p1),'p1.csv')

part<-phen_all_use[,208]#Altruism
p1<-CreateTableOne(data = as.data.frame(part))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part[which(phys3$Gender=="Male")]))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part[which(phys3$Gender=="Female")]))
write.csv(print(p1),'p1.csv')

part<-phen_all_rct_use[,12]#Washing hands
part<-ifelse(is.na(part)==F,"Yes","No")
p1<-CreateTableOne(data = as.data.frame(part),factorVars = colnames(part))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part[which(phys3$Gender=="Male")]),factorVars = colnames(part))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part[which(phys3$Gender=="Female")]),factorVars = colnames(part))
write.csv(print(p1),'p1.csv')

part<-phen_all_use[,c(372,377,378,379)]#Altruism
colnames(part)<-c("village center","main road","health center","deforestation")
p1<-CreateTableOne(data = as.data.frame(part),factorVars = colnames(part)[4])
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part[which(phys3$Gender=="Male"),]),factorVars = colnames(part)[4])
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part[which(phys3$Gender=="Female"),]),factorVars = colnames(part)[4])
write.csv(print(p1),'p1.csv')

#Income

part<-array(NA,dim=c(1187,1))
for(i in 1:1187){
    if((is.na(phen_all_use[i,182])==F)){
      if((phen_all_use[i,182]!="Dont_Know")|(phen_all_use[i,182]!="Refused")){
        part[i,1]<-as.numeric(as.character(phen_all_use[i,182]))
      }else{
        part[i,1]<-NA
      }
    }
}
p1<-CreateTableOne(data = as.data.frame(part))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part[which(phys3$Gender=="Male")]))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part[which(phys3$Gender=="Female")]))
write.csv(print(p1),'p1.csv')

part<-phen_all_use[,181]#Travel
p1<-CreateTableOne(data = as.data.frame(part))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part[which(phys3$Gender=="Male")]))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part[which(phys3$Gender=="Female")]))
write.csv(print(p1),'p1.csv')

#Household and water

ind_list5<-c(310:316,306:309,317:359,375,376)
il5<-c("Electricity (N=1069)","Radio (N=519)","TV (N=565)","Cellphone (N=927)","No cellphone (N=12)","Refrigerator (N=397)","None electronics (N=38)","Chimney (N=959)","No chimney (N=151)","Stove (N=61)","No stove (N=1)","Wood (N=1109)","Gas(fuel) [N=46]","Electricity fuel (N=14)","Kerosene (N=2)","None fuel (N=1)","Separate kitchen (N=1106)","Cement floor (N=552)","Earth/Sand floor (N=507)","Ceramic floor (N=88)","Tiles floor (N=19)","Mud bricks floor (N=1)","Wood floor (N=0)","Other floor (N=5)","Wooden windows (N=1035)","Glass windows (N=51)","Metal windows (N=6)","Unfinished windows (N=35)","No windows (N=45)","Clay/mud walls(N=693)","Clay brick walls (N=3)","Cement walls (N=452)","Cane/palm/trunks walls (N=0)","Wood unpolished walls (N=21)","Wood polished walls (N=0)","Discarded materials walls (N=1)","No walls (N=0)","Other walls (N=2)","Plastic roof (N=7)","Metal roof (N=982)","Clay roof (N=120)","Thatch/palm roof (N=10)","Concrete roof (N=27)","Wood roof (N=22)","Other roof(N=4)","Sleeping rooms","Spring(protected) [N=917]","Spring(unprotected) [N=19]","Tube well (N=72)","Dug well(protected) [N=78]","Dug well(unprotected) [N=26]","Surface water (N=10)","Bottled water (N=15)","Other water (N=6)","Household size","Household wealth index")
ind_list5_c<-c(96,97,47,48,51,54,58,60,67,71,73,81,87:88,91,92,94)
ind_list5_c<-ind_list5_c-41

part<-phen_all_use[,310:316]#Travel
colnames(part)<-c("Electricity (N=1069)","Radio (N=519)","TV (N=565)","Cellphone (N=927)","No cellphone (N=12)","Refrigerator (N=397)","None electronics (N=38)")
part<-ifelse(is.na(phen_all_use[,310:316])==F,"Yes","No")
p1<-CreateTableOne(data = as.data.frame(part),factorVars = colnames(part))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part[which(phys3$Gender=="Male"),]),factorVars = colnames(part))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part[which(phys3$Gender=="Female"),]),factorVars = colnames(part))
write.csv(print(p1),'p1.csv')

part<-phen_all_use[,375:376]
colnames(part)<-c("hh_size","HW")
p1<-CreateTableOne(data = as.data.frame(part),factorVars = colnames(part)[2])
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part[which(phys3$Gender=="Male"),]),factorVars = colnames(part)[2])
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part[which(phys3$Gender=="Female"),]),factorVars = colnames(part)[2])
write.csv(print(p1),'p1.csv')

part<-phen_all_use[,ind_list5[ind_list5_c]]
colnames(part)<-il5[ind_list5_c]
part<-part[,5:dim(part)[2]]

part2<-part[,1:8]
part2<-ifelse(part[,1:8]==0,"No","Yes")
p1<-CreateTableOne(data = as.data.frame(part2),factorVars = colnames(part2))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part2[which(phys3$Gender=="Male"),]),factorVars = colnames(part2))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part2[which(phys3$Gender=="Female"),]),factorVars = colnames(part2))
write.csv(print(p1),'p1.csv')

part<-part[,9]
#colnames(part)<-c("hh_size","HW")
p1<-CreateTableOne(data = as.data.frame(part))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part[which(phys3$Gender=="Male")]))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part[which(phys3$Gender=="Female")]))
write.csv(print(p1),'p1.csv')

part2<-part[,10:13]
part2<-ifelse(part[,10:13]==0,"No","Yes")
p1<-CreateTableOne(data = as.data.frame(part2),factorVars = colnames(part2))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part2[which(phys3$Gender=="Male"),]),factorVars = colnames(part2))
write.csv(print(p1),'p1.csv')
p1<-CreateTableOne(data = as.data.frame(part2[which(phys3$Gender=="Female"),]),factorVars = colnames(part2))
write.csv(print(p1),'p1.csv')

##Variance explained supplementary


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
farm_ind<-c(64,67:75)
wild_ind<-c(63,65,76:85)
food_ind<-c(23:41)
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
pet_ind<-c(60:62,66)
farm_ind<-c(64,67:75)
wild_ind<-c(63,65,76:85)
food_ind<-c(23:41)

health_var<-c(sum(pll$R2[phys_ind]),sum(pll$R2[acute_ind]),sum(pll$R2[chronic_ind]),sum(pll$R2[med_ind]),sum(pll$R2[person_ind]),sum(pll$R2[cog_ind]),sum(pll$R2[unfav_ind]),sum(pll$R2[anx_ind]),sum(pll$R2[dep_ind]))
#sum(health_var)
food_an_var<-c(sum(pll$R2[pet_ind]),sum(pll$R2[farm_ind]),sum(pll$R2[wild_ind]),sum(pll$R2[food_ind]))
soc_eco_var<-c(sum(pll$R2[partner_ind]),sum(pll$R2[edu_ind]),sum(pll$R2[friend_ind]),sum(pll$R2[family_ind]),sum(pll$R2[all_ind]),sum(pll$R2[risky_ind]),sum(pll$R2[village_ind]),sum(pll$R2[income_ind]),sum(pll$R2[hh_ind]),sum(pll$R2[water_ind]))
tech_var<-sum(pll$R2[tech_factors_ind])

health_var<-health_var*100
food_an_var<-food_an_var*100
soc_eco_var<-soc_eco_var*100
tech_var<-tech_var*100

#Pathway

health_var_pwy<-c(sum(pll_pwy$R2[phys_ind]),sum(pll_pwy$R2[acute_ind]),sum(pll_pwy$R2[chronic_ind]),sum(pll_pwy$R2[med_ind]),sum(pll_pwy$R2[person_ind]),sum(pll_pwy$R2[cog_ind]),sum(pll_pwy$R2[unfav_ind]),sum(pll_pwy$R2[anx_ind]),sum(pll_pwy$R2[dep_ind]))

food_an_var_pwy<-c(sum(pll_pwy$R2[pet_ind]),sum(pll_pwy$R2[farm_ind]),sum(pll_pwy$R2[wild_ind]),sum(pll_pwy$R2[food_ind]))
soc_eco_var_pwy<-c(sum(pll_pwy$R2[partner_ind]),sum(pll_pwy$R2[edu_ind]),sum(pll_pwy$R2[friend_ind]),sum(pll_pwy$R2[family_ind]),sum(pll_pwy$R2[all_ind]),sum(pll_pwy$R2[risky_ind]),sum(pll_pwy$R2[village_ind]),sum(pll_pwy$R2[income_ind]),sum(pll_pwy$R2[hh_ind]),sum(pll_pwy$R2[water_ind]))
tech_var_pwy<-sum(pll_pwy$R2[tech_factors_ind])

health_var_pwy<-health_var_pwy*100
food_an_var_pwy<-food_an_var_pwy*100
soc_eco_var_pwy<-soc_eco_var_pwy*100
tech_var_pwy<-tech_var_pwy*100

#inDFss<-cbind(c("Physiological","Acute","Chronic","Medication","Personality","Cognitive","Unfavorable","Anxiety","Depression","Food & animals","Socio-economic","Technical factors"),c(health_var,food_an_var,soc_eco_var,tech_var))
inDFss<-cbind(c("Physiological measurements","Acute conditions","Chronic conditions","Medications","Personalities","Cognitive","Unfavorable habits","Anxiety","Depression","Pets","Farm animals","Wild animals","Food","Co-habiting partners","Education","Friendship","Family","All relationships","Risky behavior","Village factors","Income","Household essentials","Water sources","Technical factors"),c(health_var,food_an_var,soc_eco_var,tech_var))
inDFss<-as.data.frame(inDFss)
colnames(inDFss)<-c("Data","R2")
inDFss<-cbind(inDFss,"Type"=array("Species",dim=c(dim(inDFss)[1],1)))

inDFss_pwy<-cbind(c("Physiological measurements","Acute conditions","Chronic conditions","Medications","Personalities","Cognitive","Unfavorable habits","Anxiety","Depression","Pets","Farm animals","Wild animals","Food","Co-habiting partners","Education","Friendship","Family","All relationships","Risky behavior","Village factors","Income","Household essentials","Water sources","Technical factors"),c(health_var_pwy,food_an_var_pwy,soc_eco_var_pwy,tech_var_pwy))
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

inDFss_all2$Data <- factor(inDFss_all2$Data, levels=c("Physiological measurements","Acute conditions","Chronic conditions","Medications","Personalities","Cognitive","Unfavorable habits","Anxiety","Depression","Pets","Farm animals","Wild animals","Food","Co-habiting partners","Education","Friendship","Family","All relationships","Risky behavior","Village factors","Income","Household essentials","Water sources","Technical factors"),ordered=TRUE)

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
  geom_text(aes(label = smm2,y=y_pos123,x=x_pos123),size=5, color = "black")+theme(text = element_text(size=20),panel.grid.major = element_blank(),panel.background = element_blank())+scale_fill_manual(values=c("#377eb8","#0de0c9","#33a02c","#6a3d9a","#a50026","#dd3497","#dfc27d","#ff7f00","#000000","#a6cee3","#d73027","#525252","#3690c0","#41ab5d","#dd3497","#ae017e","#006837","#49006a","#fcc5c0","#fc9272","#fed976","#737373","661807","#cab2d6"))+
  ylab("Variance explained (R2)")#+ theme_classic()#health_var
g

ggplot(inDFss_all2,aes(x = Type, y=R2,fill=Data))  + geom_bar(position="stack",stat = "identity",width=0.5) + 
  theme(text = element_text(size=20),panel.grid.major = element_blank(),panel.background = element_blank(),axis.title.y=element_text(size=20),axis.text.x=element_text(size=15,color="black"),axis.text.y=element_text(size=15,color="black"),axis.line.y = element_line(color="black"),axis.line.x = element_line(color="black"))+scale_fill_manual(values=c("#377eb8","#0de0c9","#33a02c","#6a3d9a","#a50026","#dd3497","#dfc27d","#ff7f00","#000000","#a6cee3","#d73027","#525252","#3690c0","#41ab5d","#dd3497","#ae017e","#006837","#49006a","#fcc5c0","#fc9272","#fed976","#737373","661807","#cab2d6"))+
  ylab("Variance explained (R2)")+guides(fill=guide_legend(title="Sub-category"))#+#+ theme_classic()#health_var



##Phenotype-phenotype correlation

jh<-read.csv('health_eff.csv',row.names = 1)
colnames(jh)<-c("Hb A1c","MAP","BMI","Heart rate","Oxygen saturation","Hb total","Poor (N=153)","Fair (N=612)","Very good (N=56)","Excellent (N=75)","Cough (N=257)","Diarrhea (N=73)","Bristol stool scale","Diabetes (N=24)","Allergies (N=104)","Heart disease (N=47)","Asthma (N=47)","Stomach illness (N=175)","Intestinal illness (N=70)","Arthritis (N=39)","Painkillers (N=714)","Antibiotics (N=137)","Anti-diarrheal (N=49)","Anti-parasitic (N=34)","Anti-fungal (N=82)","Vitamins (N=121)","Anti-hypertensive (N=72)","Reserved","Nervous","Openess","Cognitive impairment (N=56)","Dementia (N=89)","Alcohol daily frequency","Cigarette usage (N=46)","Cigarette frequency","Mild (GAD7)","Moderate (GAD7)","Severe (GAD7)","Mild (PHQ9)","Moderate (PHQ9)","Severe (PHQ9)")

jh2<-read.csv('food_an_eff.csv',row.names=1)
colnames(jh2)<-c("Cat (N=645)","Dog (N=944)","Parakeet (N=102)","None pet (N=88)","Cow (N=397)","Goat (N=40)","Pig (N=208)","Chicken (N=1094)","Duck (N=478)","Turkey (N=322)","Sheep (N=33)","Geese (N=68)","Horse (N=314)","Rabbit (N=120)","None farm (N=69)","Mice (N=585)","Bat (N=410)","Lizard (N=423)","Snake (N=437)","Bird (N=794)","Possum (N=467)","Rat (N=457)","Squirrel (N=480)","None wild (N=243)","Beans","Tortillas","Rice","Bread","Milk","Yogurt","Cream/butter","Cheese","Eggs","Vegetables","Fruits","Natural juice","Chicken","Beef/pork","Ham/sausages/hotdog","Fish")

jh3<-read.csv('socio_eco_eff.csv',row.names=1)
colnames(jh3)<-c("Partner live duration","Number of partners","Partner live age","Living with partner (N=375)","Grades 1-3 (N=373)","Grades 4-6 (N=311)","Grades >6 (N=113)","Friend ties(same building)","Friend ties(different building)","Betweeness(friendship)","Transitivity(friendship)","Family ties(same building)","Family ties(different building)","Transitivity(familial)","Betweeness(familial)","Degree(all ties)","Clustering coefficient(all ties)","Betweeness(all ties)","Kin percentage (upto third degree)","Altruism","Risk taking","Washing hands (N=1130)","Distance to village center","Distance to main road","Distance to health center","Deforestation (%)","Travel","Monthly expenditure","Household size","Household wealth index","Refrigerator (N=397)","None electronics (N=38)","Stove (N=61)","Gas(fuel) [N=46]","Separate kitchen (N=1106)","Earth/Sand floor (N=507)","Glass windows (N=51)","Clay/mud walls(N=693)","Cement walls (N=452)","Metal roof (N=982)","Sleeping rooms","Spring(protected) [N=917]","Dug well(protected) [N=78]","Dug well(unprotected) [N=26]","Bottled water (N=15)")

jh4<-cbind(cbind(jh,jh2),jh3)
jh4<-jh4[,1:125]#Removing bottled water
colnames(jh4)[c(14,26,34,27)]<-c("Diabetes (N=20)","Vitamins (N=122)","Cigarette usage (N=56)","Anti-hypertensive (N=43)")

jh5<-read.csv('tpp.csv',row.names=1)
colnames(jh5)[1:81]<-colnames(jh4)[1:81]
colnames(jh5)[82]<-"Diet diversity score"
colnames(jh5)[83:126]<-colnames(jh4)[82:125]


cp_cor<-cor(jh5,jh5,method="pearson")

pheatmap(cp_cor,treeheight_col = 0)

library(corrplot)
#Save 5.00 x 6.00
corrplot(as.matrix(cp_cor),method = "square",tl.cex = 0.85,cl.cex = 0.85,tl.col = "black",order= "hclust")


cp_cor2<-cp_cor

paletteLength<-15
myColor <- colorRampPalette(c("#cb181d","white", "#33a02c"))(paletteLength)#myColor <- colorRampPalette(c("Orange", "white", "Darkblue"))(paletteLength)
myBreaks <- c(seq(-1,0, length.out=ceiling(paletteLength/2) + 1), 
              seq(1/paletteLength,1, length.out=floor(paletteLength/2)))

myBreaks
myBreaks[7]<-as.numeric(-0.001)
myBreaks[9]<-as.numeric(0.001)
myBreaks

#myBreaks<-c(myBreaks,1)
#myColor<-c(myColor,"#525252")


pheatmap(cp_cor,annotation_row = NULL,annotation_names_row = T,labels_row = rownames(cp_cor),legend = T,fontsize_number = 15,border_color = "#EEEEEE",na_col = "grey50",fontsize_col = 10,angle_col = 90,fontsize_row = 10,fontface="bold",color=myColor,breaks=myBreaks,number_color = "black",treeheight_col = 0,treeheight_row = 0,cluster_rows = F,cluster_cols = F,show_colnames = F)
#Save 16.00 x 14.00




