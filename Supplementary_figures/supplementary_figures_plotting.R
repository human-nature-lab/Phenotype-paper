
# ==================================
# By: Shivkumar Vishnempet Shridhar, Christakis and Brito group, Yale (2023)
# Honduras microbiome project, plotting for
# Supplemenrary figures
# ==================================

#
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(ggtree)
#

effect_size_phen_lmer<-read.csv('species_phenotype_assc.csv',row.names=1)
fdr_phen_lmer<-read.csv('species_phenotype_fdr.csv',row.names=1)

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

colnames(effect_size_phen_lmer)[(1:24)+41]<-c("Cat (N=960)","Dog (N=1384)","Parakeet (N=157)","No pets (N=220)","Cow (N=535)","Goat (N=63)","Pig (N=269)","Chicken (N=1664)","Duck (N=636)","Turkey (N=413)","Sheep (N=51)","Geese (N=82)","Horse (N=468)","Rabbit (N=176)","No farm animals (N=171)","Mice (N=820)","Bat (N=546)","Lizard (N=534)","Snake (N=546)","Bird (N=660)","Possum (N=588)","Rat (N=590)","Squirrel (N=590)","No wild animals (N=591)")
colnames(effect_size_phen_lmer)[(c(31,36,38:39,41))+41]<-c("Cream/butter","Natural juice","Beef/pork","Ham/sausages/hotdog","Diet diversity score")

colnames(effect_size_phen_lmer)[1:17+82]<-c("Number of partners","Living with partner (N=1232)","Grades 1-3 (N=920)","Grades 4-6 (N=528)","Grades >6 (N=129)","Friend ties (same building)","Friend ties (different building)","Betweenness (friendship)","Transitivity (friendship)","Familial ties (same building)","Familial ties (different building)","Betweenness (familial)","Transitivity (familial)","Degree (all ties)","Clustering coefficient (all ties)","Betweenness (all ties)","Kin percentage (to third degree)")
colnames(effect_size_phen_lmer)[19:23+82]<-c("Risk taking","Washing hands (N=1603)","Distance to village center","Distance to main road","Number of churches")
colnames(effect_size_phen_lmer)[26:41+82]<-c("Monthly expenditure","Household size","Household wealth index","TV (N=958)","No electronics (N=89)","Eath/Sand floor (N=594)","Ceramic floor (N=151)","Glass windows (N=86)","Unfinished windows (N=59)","Clay/mud walls (N=1256)","Cement walls (N=402)","Concrete roof (N=37)","Sleeping rooms","Spring (protected) (N=1441)","Tube well (N=91)","Dug well (protected) (N=93)")


##Figure S2

ph4<-dist(t(effect_size_phen_lmer),method="euclidean")
ph5<-hclust(ph4)

ph_cat<-array(NA,dim=c(dim(ph)[2],1))

ph_cat[1:41,1]<-"#4472c4"
ph_cat[42:82,1]<-"#ff0000"
ph_cat[83:123,1]<-"#548235"


pdf('figure_S2.pdf', width=8, height=15)
ggtree(ph5)+geom_tiplab(hjust=-0.06)+
  geom_tippoint(color=ph_cat[,1], size=3, alpha=.75)+xlim(-110,50)+#+ #+
  theme_tree2(legend.position='right')
dev.off()  


##Figure S4

pll<-read.csv('pl_perm2_all_123.csv',row.names=1)
pll_pwy<-read.csv('pl_perm2_all_pwy_123.csv',row.names=1)

rownames(pll)
tech_factors_ind<-c(154:159)
phys_ind<-c(1:10)
acute_ind<-c(11:12)
chronic_ind<-c(13:19)
med_ind<-c(20:26)
person_ind<-c(27:29)
cog_ind<-c(30:31)
anx_ind<-c(35:38)
dep_ind<-c(39:42)
unfav_ind<-c(32:34)
#farm_ind<-c(64,67:75)
#wild_ind<-c(63,65,76:85)
#food_ind<-c(23:41)
partner_ind<-c(84:85)
edu_ind<-c(86:88)
friend_ind<-c(89:92)
family_ind<-c(93:96)
all_ind<-c(97:100)
income_ind<-c(108:109)
risky_ind<-c(101:103)
village_ind<-c(104:107)
hh_ind<-c(110:146)
water_ind<-c(147:153)
pet_ind<-c(43:46)
farm_ind<-c(47:57)
wild_ind<-c(58:66)
food_ind<-c(67:83)

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

##Ordered

inDFss_all2<-inDFss_all

inDFss_all2$Data <- factor(inDFss_all2$Data, levels=c("Physiological measurements","Acute conditions","Chronic conditions","Medications","Personalities","Cognitive","Unfavorable habits","Anxiety","Depression","Pets","Farm animals","Wild animals","Food","Co-habiting partners","Education","Friendship","Family","All relationships","Risky behavior","Village factors","Income","Household essentials","Water sources","Technical factors"),ordered=TRUE)


library(ggplot2)

pdf('figure_S4.pdf',width=12,height=8)
ggplot(inDFss_all2,aes(x = Type, y=R2,fill=Data))  + geom_bar(position="stack",stat = "identity",width=0.5) + 
  theme(text = element_text(size=25),panel.grid.major = element_blank(),panel.background = element_blank(),axis.title.y=element_text(size=25),axis.text.x=element_text(size=20,color="black"),axis.text.y=element_text(size=20,color="black"),axis.line.y = element_line(color="black"),axis.line.x = element_line(color="black"))+scale_fill_manual(values=c("#377eb8","#0de0c9","#33a02c","#6a3d9a","#a50026","#dd3497","#dfc27d","#ff7f00","#000000","#a6cee3","#d73027","#525252","#3690c0","#41ab5d","#dd3497","#ae017e","#006837","#49006a","#fcc5c0","#fc9272","#fed976","#737373","661807","#cab2d6"))+
  ylab("Variance explained (R2)")+guides(fill=guide_legend(title="Sub-category"))#+#+ theme_classic()#health_var
dev.off()


##Figure S5

div_animal<-read.csv('div_alpha_ph2_animal.csv',row.names=1)

div_animal[,2]<-factor(div_animal[,2],levels=as.character(unique(div_animal[,2])),ordered = T)
color_animal<-c("#a6cee3","#1f78b4","#b2df8a","#878787","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#878787","#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#878787")

pdf('figure_S5.pdf',width=7,height=6)
ggplot(div_animal, aes(x=animals, y=as.numeric(as.character(diversity)),group=animals)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE,fill=color_animal)+labs(x="",y="Shannon diversity")+
  theme(text = element_text(size=15),panel.background=element_blank(),legend.position = "none",axis.line.x = element_line(color="black", size = 0.5),axis.line.y = element_line(color="black", size = 0.5),axis.title.y=element_text(size=15),axis.text.y = element_text(color="black"),axis.text.x = element_text(color="black",angle=60,hjust=1))+ 
  scale_x_discrete(breaks=seq(1,24,1))+scale_x_discrete(labels=as.character(unique(div_animal[,2])))+#++
  theme(plot.margin=margin(c(0.5,0.5,0.5,0.5), unit = "cm"))
dev.off()








