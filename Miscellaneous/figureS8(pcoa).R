
# ==================================
# By: Shivkumar Vishnempet Shridhar, Christakis and Brito group, HNL, Yale (2023)
# Honduras microbiome project, script for creating PCoA (Principal coordinate analysis plot) for overall microbiome
# ==================================



#Pcoa


st<-read.csv('st_123.csv',row.names = 1)#Bristol stool scale
phen_all_use<-read.csv('phen_b4.csv',row.names=1)#sample x phenotypes
mb_samp<-read.csv('mb_sp_10.csv',row.names=1)#sample x species

which((grepl("SGB1626",rownames(mb_samp),fixed=T))==T)

sp<-as.data.frame(mb_samp[which((grepl("SGB1626",rownames(mb_samp),fixed=T))==T),])

st<-st[,1]

cleaned_data<-array(0,dim=c(dim(phen_all_use)[1],8))
for(k in 1:dim(phen_all_use)[1]){
  cleaned_data[k,1]<-as.numeric(phen_all_use[k,6])
  cleaned_data[k,2]<-as.numeric(phen_all_use[k,4])
  cleaned_data[k,3]<-as.numeric(phen_all_use[k,224])
  cleaned_data[k,4]<-as.numeric(phen_all_use[k,225])
  cleaned_data[k,5]<-as.numeric(phen_all_use[k,40])
  cleaned_data[k,6]<-as.numeric(st[k])
  if(is.na(phen_all_use[k,15])==F){
    if(phen_all_use[k,15]=="Poor"){
      cleaned_data[k,7]<-0
    }else if(phen_all_use[k,15]=="Fair"){
      cleaned_data[k,7]<-1
    }else if(phen_all_use[k,15]=="Good"){
      cleaned_data[k,7]<-2
    }else if(phen_all_use[k,15]=="Very good"){
      cleaned_data[k,7]<-3
    }else if(phen_all_use[k,15]=="Excellent"){
      cleaned_data[k,7]<-4
    }
  }
  if(is.na(phen_all_use[k,15])==T){
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
dist_sp <- vegdist(t(mb_samp[,ind_temp]),  method = "bray")
# do PCOA
#myPCOA <- pcoa(dist_sp)
myPCOAv2= cmdscale(dist_sp)
myPCOAv2[,2]<-myPCOAv2[,2]*(-1)#Just to re orient axes
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
  theme_bw( ) +xlab ("PCoA1") + ylab ("PCoA2")+  #+ xlim(-0.4,0.6) +ylim(-0.4,0.4)
  theme(axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black"))

#Save 4.00 x 7.00


