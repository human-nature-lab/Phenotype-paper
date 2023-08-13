
# ==================================
# By: Shivkumar Vishnempet Shridhar, Christakis and Brito group, HNL, Yale (2023)
# Honduras microbiome project, plotting for
# food&animal-phenotype association
# ==================================


##CLR transforming species abundances:

kl<-read.csv('mb_sp_10.csv',row.names=1)
mb_samp2<-kl
taxa<-as.numeric(unlist(mb_samp2))
taxa25<-mb_samp2+min(taxa[which(as.numeric(taxa)>0)])/2
dim(taxa)

min(taxa[which(as.numeric(taxa)>0)])/2

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x), na.rm=na.rm) / length(x))
 }
 Gmean_core = apply(taxa25, 1, gm_mean)
 data_prepared = cbind(Gmean_core,taxa25)
 data_transformed = t(apply(data_prepared,1,function(x){
   log(x / x[1])[-1]
 }))
 dim(data_transformed)

##For pathways
pathways_use<-read.csv('pathways_use.csv',row.names=1)

mb_samp2<-pathways_use
taxa<-as.numeric(unlist(mb_samp2))
taxa25<-mb_samp2+min(taxa[which(as.numeric(taxa)>0)])/2
dim(taxa)

min(taxa[which(as.numeric(taxa)>0)])/2
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x), na.rm=na.rm) / length(x))
 }
 Gmean_core = apply(taxa25, 1, gm_mean)
 data_prepared = cbind(Gmean_core,taxa25)
 data_transformed = t(apply(data_prepared,1,function(x){
   log(x / x[1])[-1]
 }))
 dim(data_transformed)



######################################################################################################
### Figure 1

phen_all_use<-read.csv('phen_b4.csv',row.names=1)#sample x phenotypes
phen_all_use<-phen_all_use[which(!(phen_all_use$village_code%in%c(3,18,45,61,81,94,128,134,157))),]
mb_samp_sp<-read.csv('mb_sp.csv',row.names=1)#sample x species

#
library(ggmap)
register_google(key="....") #put in your key
library(ggplot2)
library('igraph')
library(lmerTest)
#

#Need to recompute Bray-distance within mariposal -- and re-assign colors



g_lon<-read.csv('g_lon.csv',row.names = 1)
g_lat<-read.csv('g_lat.csv',row.names = 1)
g_lon2<-read.csv('g_lon2.csv',row.names = 1)
g_lat2<-read.csv('g_lat2.csv',row.names = 1)
g_lon<-g_lon[,1]
g_lon2<-g_lon2[,1]
g_lat<-g_lat[,1]
g_lat2<-g_lat2[,1]
edges<-read.csv('ggmap_edges.csv')
edges<-edges[,2:dim(edges)[2]]
bn_shape<-read.csv('bn_shape.csv',row.names = 1)
zm_vc_v<-read.csv('zm_vc_v.csv',row.names=1)

ind_vil<-match(zm_vc_v[,1],phen_all_use$respondent_master_id)

mb_mean_vil<-rowSums(mb_samp_sp[,ind_vil])/length(ind_vil)

br<-array(NA,dim=c(dim(zm_vc_v)[1],3))
br[,1]<-as.character(zm_vc_v[,1])

for(i in 1:dim(br)[1]){
  br[i,2]<-as.numeric(vegdist(t(cbind(mb_mean_vil,mb_samp_sp[,ind_vil[i]])), "bray"))
}

bray_col<-array(NA,dim=c(dim(zm_vc_v)[1]))
for(i in 1:length(bray_col)){
  if(br[i,2]==0){
    bray_col[i]<-"#000000"
  }else if((br[i,2]>0)&(br[i,2]<0.36)){
    bray_col[i]<-"#006837"
  }else if((br[i,2]>=0.36)&(br[i,2]<0.42)){
    bray_col[i]<-"#1a9850"
  }else if((br[i,2]>=0.42)&(br[i,2]<0.48)){
    bray_col[i]<-"#66bd63"
  }else if((br[i,2]>=0.48)&(br[i,2]<0.54)){
    bray_col[i]<-"#a6d96a"
  }else if((br[i,2]>=0.54)&(br[i,2]<0.6)){
    bray_col[i]<-"#d9ef8b"
  }else if((br[i,2]>=0.6)&(br[i,2]<0.66)){
    bray_col[i]<-"#fee08b"
  }else if((br[i,2]>=0.66)&(br[i,2]<0.72)){
    bray_col[i]<-"#fdae61"
  }else if((br[i,2]>=0.72)&(br[i,2]<0.78)){
    bray_col[i]<-"#f46d43"
  }else if((br[i,2]>=0.78)&(br[i,2]<0.84)){
    bray_col[i]<-"#d73027"
  }else if((br[i,2]>=0.84)&(br[i,2]<0.9)){
    bray_col[i]<-"#a50026"
  }else if(br[i]>=0.9){
    bray_col[i]<-"#a50026"
  }
}

br[,3]<-bray_col

houses <- get_map(location = c(mean(g_lon), mean(g_lat)), source="google", zoom=17,maptype="hybrid")#map_type- terrain-background , satellite (USE) or hybrid,terrain#Gives with no labels but yellow color--- maybe the one we want?
p<-ggmap(houses)

p
plot_vector<- as.data.frame(cbind(g_lon,g_lat))
p + geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2), data=edges, size = 0.9, colour="white",alpha=0.8) + geom_point(aes(g_lon2,g_lat2), data=plot_vector)+ geom_point(aes(x=g_lon2,y=g_lat2),plot_vector,size=4,alpha=1,color=adjustcolor(bray_col, alpha.f = 1),shape=as.character(bn_shape[,1]))+coord_fixed(ylim=c(14.9448,14.9475), ratio=1/cos(pi*14.9448/180))+xlab("Longitude")+ylab("Latitude")+theme(text=element_text(size=20))


#Bigger map plot --19 villages


##Fig 1A2

vc<-unique(phen_all_use$village_code)
vc<-vc[which(is.na(vc)==F)]

for(i in 1:length(vc)){
  if(i==1){
    vil_lat<-mean(phen_all_use$building_latitude[which(phen_all_use$village_code==vc[i])])
    vil_lon<-mean(phen_all_use$building_longitude[which(phen_all_use$village_code==vc[i])])
    vil_size<-length(which(phen_all_use$village_code==vc[i]))
    vil_name<-as.character(unique(phen_all_use$village_name[which(phen_all_use$village_code==vc[i])]))#was village name
  }else{
    vil_lat<-rbind(vil_lat,mean(phen_all_use$building_latitude[which(phen_all_use$village_code==vc[i])]))
    vil_lon<-rbind(vil_lon,mean(phen_all_use$building_longitude[which(phen_all_use$village_code==vc[i])]))
    vil_size<-rbind(vil_size,length(which(phen_all_use$village_code==vc[i])))
    vil_name<-rbind(vil_name,as.character(unique(phen_all_use$village_name[which(phen_all_use$village_code==vc[i])])))
  }
}

vil_cood<-cbind(vil_lon,vil_lat)
vil_cood<-cbind(vil_cood,vil_size)
vil_cood<-as.data.frame(vil_cood)
colnames(vil_cood)<-c("vil_lon","vil_lat","vil_size")

#Village plotting
houses <- get_map(location = c( mean(vil_lon),mean(vil_lat)), source="google", zoom=8,maptype="hybrid",style = c(feature="country",element="labels",visibility='off'))#map_type- terrain-background , satellite (USE) or hybrid,terrain#Gives with no labels but yellow color--- maybe the one we want?
p<-ggmap(houses)
p

p + geom_point(aes(x=vil_lon,y=vil_lat),vil_cood,size=2,alpha=1,color="black")+coord_fixed(ylim=c(13.7,15.8),xlim=c(-90.8,-88.1), ratio=1/cos(pi*14.9448/180))+xlab("Longitude")+ylab("Latitude")+theme(text=element_text(size=18),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"))+ geom_point(aes(x=vil_lon,y=vil_lat),vil_cood,size=1.2,alpha=1,color="#fdae61")#color=adjustcolor(age_col, alpha.f = 1),

#Save 7.00 x 9.00



##Fig 1C

#vc<-unique(phen_all_use3$village_code)
i<-1

avg_coord<-as.matrix(vil_cood[,1:2])

tsp<-mapdist("14.94739+-88.91059","14.94443+-88.91376",mode="walking")

paste0(as.character(phen_all_use$building_latitude[261]),"+",as.character(phen_all_use$building_longitude[261]))

for(j in 1:dim(phen_all_use)[1]){
  for(i in 1:length(vc)){
    if(j!=636){
      if(phen_all_use$village_code[j]==vc[i]){
        tsp<-mapdist(paste0(as.character(phen_all_use$building_latitude[j]),"+",as.character(phen_all_use$building_longitude[j])),paste0(as.character(avg_coord[i,2]),"+",as.character(avg_coord[i,1])),mode="walking")
        phen_all_use$distance_center_build[j]<-tsp$km 
      }
    } 
  }
}

library(ggplot2)
phen_all_use2<-phen_all_use[which(is.na(phen_all_use$village_code)==F),]

ggplot(phen_all_use2, aes(x=distance_center, y=bray_vil_mean)) + geom_point(aes(color=as.character(village_code)),size=1)+xlim(0,2)+xlab("Distance from weighted village center(km)")+ylab("Bray-Curtis dissimilarity")+ scale_color_brewer(palette="Paired")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),text = element_text(size = 20),legend.position="none",axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black"))+# +scale_x_continuous(trans='log2') +#legend.key = element_rect(colour = "white", fill = NA)
  geom_smooth(method='lm',formula=y~x,color="black")+guides(fill=guide_legend(title="Village code"))+ labs(fill = "Village code")+
  theme(axis.line.x=element_line(colour="black"),axis.line.y=element_line(colour="black"))

#Save 5.00 x 6.00

mss<-lmer(bray_vil_mean~distance_center,data=phen_all_use2)
summary(mss)

