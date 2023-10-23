#######################################################
###Plotting


file_list<-list.files(path=".")
io<-682

kl<-read.table(file_list[io],header = T)

ind_match<-match(kl$sample,rownames(phen_all_rct_use))
ind_keep<-which(is.na(ind_match)==F)

tlk<-cbind(kl$percentage_of_polymorphic_sites[ind_keep],as.character(phen_all_use[ind_match[ind_keep],83]))
tlk[,2]<-ifelse(tlk[,2]=="Every day",4,tlk[,2])
tlk[,2]<-ifelse(tlk[,2]=="A few days per week",3,tlk[,2])
tlk[,2]<-ifelse((tlk[,2]=="A few days per month")|(tlk[,2]=="A few times per month"),2,tlk[,2])
tlk[,2]<-ifelse(tlk[,2]=="Never/rarely",1,tlk[,2])

div_alpha_ph2<-as.data.frame(tlk[which(is.na(tlk[,2])==F),])
colnames(div_alpha_ph2)<-c("poly","wealth")
div_alpha_ph2<-div_alpha_ph2[which(div_alpha_ph2[,2]!=5),]

ggplot(div_alpha_ph2, aes(x=wealth, y=as.numeric(as.character(poly)),group=wealth)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#b2df8a","#33a02c","#e31a1c"))+labs(x="Food consumption (frequency)",y="% of polymorphic sites")+
  theme(text = element_text(size=15),panel.background=element_blank(),axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y=element_text(size=15),axis.text.y = element_text(color="black"))+ scale_x_discrete(breaks=seq(1,5,1))+#+ylim(1.5,5.2)#,"#ff7f00","#cab2d6","#6a3d9a","#ffff99"#x="Relationship",
  ggtitle(hk_short[io]) +
  theme(plot.title = element_text(hjust = 0.5))#+ylim(0,1)

x1<-div_alpha_ph2$poly[which(div_alpha_ph2$wealth=="1")]
y1<-div_alpha_ph2$poly[which(div_alpha_ph2$wealth=="2")]
x1<-x1[which(is.na(x1)==F)]
y1<-y1[which(is.na(y1)==F)]

wilcox.test(as.numeric(as.character(x1)),as.numeric(as.character(y1)))


library(ggpubr)
my_comparisons <- list(c("1", "3"),c("1", "4"),c("2", "4"))
#bray_chronic2<-bray_chronic[which(is.na(bray_chronic[,1])==F),]
pdf(file='figure_4a.pdf',width=5,height=6)
ggplot(div_alpha_ph2, aes(x=wealth, y=as.numeric(as.character(poly)),group=wealth)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#b2df8a","#33a02c","#e31a1c"))+labs(x="",y="% of polymorphic site")+
  theme(text = element_text(size=18),axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black"),panel.background=element_blank(),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))+ scale_x_discrete(breaks=seq(1,4,1))+
  stat_compare_means(comparisons = my_comparisons, label.y = c(7.1,7.6,8.1))+ylim(0,9)+#++
  theme(plot.margin=margin(c(0.5,0.5,2,0.5), unit = "cm"))
dev.off()
#Save 6.00 x 5.00


#Figure 4B

io<-4

kl<-read.table(file_list[io],header = T)

ind_match<-match(kl$sample,rownames(phen_all_use))
ind_keep<-which(is.na(ind_match)==F)

tlk<-cbind(kl$percentage_of_polymorphic_sites[ind_keep],phen_all_use$household_wealth_index[ind_match[ind_keep]])

div_alpha_ph2<-as.data.frame(tlk)
colnames(div_alpha_ph2)<-c("poly","wealth")

#if((length(which(tlk[,2]==1))>0)&(length(which(tlk[,2]==2))>0)&(length(which(tlk[,2]==3))>0)&(length(which(tlk[,2]==4))>0)&(length(which(tlk[,2]==5))>0)){

ggplot(div_alpha_ph2, aes(x=wealth, y=as.numeric(as.character(poly)),group=wealth)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#b2df8a","#33a02c","#e31a1c","#fdbf6f"))+labs(x="",y="% of polymorphic sites")+
  theme(text = element_text(size=15),panel.background=element_blank(),axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y=element_text(size=15),axis.text.y = element_text(color="black"))+ scale_x_discrete(breaks=seq(1,5,1))+#+ylim(1.5,5.2)#,"#ff7f00","#cab2d6","#6a3d9a","#ffff99"#x="Relationship",
  ggtitle(hk_short[io]) +
  theme(plot.title = element_text(hjust = 0.5))#+ylim(0,3.5)


x1<-div_alpha_ph2$poly[which(div_alpha_ph2$wealth=="3")]
y1<-div_alpha_ph2$poly[which(div_alpha_ph2$wealth=="4")]
x1<-x1[which(is.na(x1)==F)]
y1<-y1[which(is.na(y1)==F)]

wilcox.test(as.numeric(as.character(x1)),as.numeric(as.character(y1)))


library(ggpubr)
my_comparisons <- list(c("1", "4"),c("1", "5"),c("2", "5"),c("3", "5"))
#bray_chronic2<-bray_chronic[which(is.na(bray_chronic[,1])==F),]#[which(is.na(div_alpha_ph2[,2])==F),]
pdf(file='figure_4b.pdf',width=5,height=6)
ggplot(div_alpha_ph2[which(is.na(div_alpha_ph2[,2])==F),], aes(x=as.character(wealth), y=as.numeric(as.character(poly)),group=wealth)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=c("#a6cee3","#b2df8a","#33a02c","#e31a1c","#fdbf6f"))+labs(x="Household wealth index",y="% of polymorphic sites")+
  theme(text = element_text(size=18),axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black"),panel.background=element_blank(),axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5))+ scale_x_discrete(breaks=seq(1,5,1))+
  stat_compare_means(comparisons = my_comparisons, label.y = c(4.3,4.6,4.9,5.2))+ylim(0,6)+#++
  theme(plot.margin=margin(c(0.5,0.5,0.5,0.5), unit = "cm"))
dev.off()

#Save 6.00 x 5.00





