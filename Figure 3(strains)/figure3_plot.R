
# ==================================
# By: Shivkumar Vishnempet Shridhar, Francesco Beghini, Christakis and Brito group, Yale (2023)
# Honduras microbiome project, plotting for
# Figure 3
# ==================================


library(lmerTest)
require(foreach)
require(doParallel)
library(evolvability)
library(parallel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(phytools)
library(TreeTools)

st<-read.csv('st_123.csv',row.names = 1)#Bristol stool scale
st<-as.matrix.data.frame(st)
phen_all_use<-read.csv('phen_b4.csv',header=T)#Contains all other variables -- on each person
colnames(phen_all_use)[1]<-"sample_ID"
st_ind<-read.csv('st_ind_03.csv',row.names=1)#Index matching for species - to strainphlan files
ind_sp2<-st_ind[,1]

phen_all_raw<-read.csv('phen_all_raw2.csv',header=T)#Phenotype of interest (sample x phenotyp)
phen_all_raw[,1]<-phen_all_use$sample_ID
colnames(phen_all_raw)[1]<-"sample_ID"
colnames(phen_all_raw)[13]<-"bristol"

data_transformed<-read.csv('mb_sp_01_clr.csv',row.names = 1,header=T)#CLR transformed species abundance data, species x sample
#mb_samp_sp<-data_transformed
#data_transformed<-as.matrix.data.frame(data_transformed)

###################
#Merge age, sex, and other cov including village factors to phen_all_raw then filter

colnames(phen_all_use)[c(6,4,224,225,371,11)]
xm<-tabulate(phen_all_use$village_code)
phen_all<-phen_all_raw %>% inner_join(phen_all_use%>%select(c("age_at_survey","gender","batch_effect","dna_conc","sampling_date_num","village_code","sample_ID")))
phen_all<-phen_all %>% filter(village_code%in%which(xm>5))#Discarding very small villages(N<=5)-- to account for surveying outsiders

##Accounting for people who migrated between villages to submit to survey (N=15)

data_transformed<-data_transformed[,which(phen_all_use$village_code%in%which(xm>5))]#Discarding very small villages(N<=5)-- to account for surveying outsiders

init_mat<-function(kl){
  phen<-colnames(phen_all)[kl]
  raw_data<-phen_all%>%select(phen)
  rownames(raw_data)<-phen_all$sample_ID
  colnames(raw_data)<-"phen"
  return(raw_data)
}

kl<-which(colnames(phen_all_raw)=="Household.wealth.index")
i<-25

#Set working directory to where the Strainphlan tree files are

file_list<-list.files(path=".")
tr2 <- ape::read.tree(file_list[ind_sp2[i]])
ind_stt<-match(tr2$tip.label[1:length(tr2$tip.label)],phen_all$sample_ID)
tr3<-drop.tip(tr2,which(is.na(ind_stt)==T))

cleaned_data<-init_mat(kl)
cleaned_data<-cleaned_data[ind_stt[which(is.na(ind_stt)==F)],]
ind_temp<-which(is.na(cleaned_data)==F)
tr4<-keep.tip(tr3,ind_temp)
cleaned_data<-as.data.frame(cleaned_data[ind_temp])
colnames(cleaned_data)<-"HHW"
ind_m<-ind_stt[which(is.na(ind_stt)==F)]
rownames(cleaned_data)<-phen_all$sample_ID[ind_m[ind_temp]]
#rownames(cleaned_data)<-tr4$tip.label

cleaned_data<-cbind(cleaned_data,array(NA,dim=c(dim(cleaned_data)[1],2)))
colnames(cleaned_data)[2]<-"color2"

cleaned_data[,2]<-ifelse(is.na(cleaned_data[,2])==T,ifelse((cleaned_data[,1]==1),"#440154",cleaned_data[,2]),cleaned_data[,2])
cleaned_data[,2]<-ifelse(is.na(cleaned_data[,2])==T,ifelse((cleaned_data[,1]==2),"#3b528b",cleaned_data[,2]),cleaned_data[,2])
cleaned_data[,2]<-ifelse(is.na(cleaned_data[,2])==T,ifelse((cleaned_data[,1]==3),"#21918c",cleaned_data[,2]),cleaned_data[,2])
cleaned_data[,2]<-ifelse(is.na(cleaned_data[,2])==T,ifelse((cleaned_data[,1]==4),"#5ec962",cleaned_data[,2]),cleaned_data[,2])
cleaned_data[,2]<-ifelse(is.na(cleaned_data[,2])==T,ifelse((cleaned_data[,1]==5),"#fde725",cleaned_data[,2]),cleaned_data[,2])

#Mid-point rooting 

tr5<-midpoint_root(tr4)
tr6<-Preorder(tr5)
tip_sub<-Subtree(tr6, 1647)#Identifying the branching which determines the strain/sub-tree
tip_sub2<-tip_sub$tip.label

colnames(cleaned_data)[3]<-"color_strain"#Coloring the samples belonging to a different strain by a different color
cleaned_data[,3]<-"#beaed4"
cleaned_data[match(tip_sub2,rownames(cleaned_data)),3]<-"#e41a1c"
cleaned_data<-as.data.frame(cleaned_data)

#Set directory where figures can be saved

ggm7<-ggtree(tr5,layout="fan",size=0.1,open.angle = 20)+#+layout_dendrogram() +
  theme_dendrogram()+
  geom_tippoint(shape=20,color=cleaned_data3$color_an)+
  geom_hilight(node=1647, alpha=1, fill=NA, color="black",
               size=0.05,extendto=0.69)#+
ggm7

gheatmap(rotate_tree(ggm7,10),data=cleaned_data5,offset=0, width=.1,#was 110 for rotate and 0 for col angle
         colnames_angle=0, colnames_offset_y = 2,legend_title = "Household wealth index")+
  scale_fill_manual(values=c("#440154","#3b528b","#21918c","#5ec962","#fde725"))

ggsave(file='figure_3a.pdf',width=8,height=7)


#####################################################################################

### Species-phenotype association with and without strains


################################
#Initializing matrix for phenotype of interest

init_mat<-function(kl,i){
  sp<-t(data_transformed%>%slice(i))
  phen<-colnames(phen_all)[kl]
  raw_data<-cbind(phen_all%>%select(c("age_at_survey","gender","batch_effect","dna_conc","BMI","sampling_date_num","bristol",phen,"village_code")),sp)
  colnames(raw_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","phen","village","species")
  rownames(raw_data)<-phen_all$sample_ID
  raw_data$phen<-ifelse(raw_data$phen%in%c("Removed","Dont_Know","Refused"),NA,raw_data$phen)
  cleaned_data<-raw_data%>%filter(complete.cases(raw_data))
  cleaned_data$phen<-as.numeric(as.character(cleaned_data$phen))
  cleaned_data$village<-factor(cleaned_data$village)
  return(cleaned_data)
}

init_mat_no_phen<-function(i){
  sp<-t(data_transformed%>%slice(i))
  raw_data<-cbind(phen_all%>%select(c("age_at_survey","gender","batch_effect","dna_conc","BMI","sampling_date_num","bristol","village_code")),sp)
  colnames(raw_data)<-c("age","sex","batch_effect","dna_conc","BMI","sampling_date","bristol","village","species")
  rownames(raw_data)<-phen_all$sample_ID
  cleaned_data<-raw_data%>%filter(complete.cases(raw_data))
  cleaned_data$village<-factor(cleaned_data$village)
  return(cleaned_data)
}

cl <- makeCluster(3)#Initializing parrallel for loops
registerDoParallel(cl)

#See if the following lines can be deleted
effect_size_phen<-array(NA,dim=c(dim(data_transformed)[1],(dim(phen_all_raw)[2]-1)))
pval_phen<-array(NA,dim=c(dim(data_transformed)[1],(dim(phen_all_raw)[2]-1)))



output <- foreach(kl=2:dim(phen_all_raw)[2], .combine=rbind) %dopar% {
  
  library(evolvability)
  library(lmerTest)
  library(dplyr)
  effect_size_phen<-array(NA,dim=length(ind_sp2))
  pval_phen<-array(NA,dim=length(ind_sp2))
  effect_size_phen_no<-array(NA,length(ind_sp2))
  pval_phen_no<-array(NA,dim=length(ind_sp2))
  phyl2<-array(NA,dim=length(ind_sp2))
  
  
  for(i in 1:length(ind_sp2)){
  if(is.na(ind_sp2[i])==F){
    
    lmm<-0
    tr2 <- ape::read.tree(file_list[ind_sp2[i]])
    
    if(length(tr3$tip.label)>50){
      

  if(!(colnames(phen_all)[kl]%in%(c("BMI","bristol")))){
    for(i in 1:dim(data_transformed)[1]){
      
      cleaned_data<-init_mat(kl,i)
      if((sum(as.numeric(as.character(cleaned_data3$phen)>0)))&(length(unique(as.numeric(as.character(cleaned_data3$phen))))>1)){
        
        ind_stt<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data))
        tr3<-drop.tip(tr2,which(is.na(ind_stt)==T))
        cleaned_data2<-cleaned_data[ind_stt[which(is.na(ind_stt)==F)],]
        cleaned_data2<-cbind(cleaned_data2,"phyl"=tr3$tip.label)
        A1 <- Matrix::Matrix(ape::vcv(tr3), sparse = TRUE)
        colnames(A1) <- rownames(A1) <-tr3$tip.label
        
        
        ms1<-lmer(species~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+phen+(1|village),
                  data = cleaned_data2)
        mss<-summary(ms1)
        #coef(mss)[which(rownames(coef(mss))=="phen"),5]
        ms3<-Almer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+phen+(1|village)+(1|phyl),
                   data = cleaned_data2,A=list(phyl=A1))
        mss2<-summary(ms3)
        
        effect_size_phen_no[i]<-coef(mss)[which(rownames(coef(mss))=="phen"),1]
        pval_phen_no[i]<-coef(mss)[which(rownames(coef(mss))=="phen"),5]
        effect_size_phen[i]<-coef(summary(ms3))[which(rownames(coef(ms3))=="phen"),1]
        pval_phen[i]<-2*pnorm(abs(coef(ms3)[which(rownames(coef(ms3))=="phen"),3]), lower.tail = FALSE)#Pvalue of phenotype
        phyl2[i]<-phylosig(tr3,cleaned_data2$phen, method="lambda", test=FALSE)#Checking for phylogenetic signal
        
      }else{
        effect_size_phen_no[i]<-0
        pval_phen_no[i]<-1
        effect_size_phen[i]<-0
        pval_phen[i]<-1
      }
    }
    
  }else{
    for(i in 1:dim(data_transformed)[1]){
      
      cleaned_data<-init_mat_no_phen(i)
      
      ind_stt<-match(tr2$tip.label[1:length(tr2$tip.label)],rownames(cleaned_data))
      tr3<-drop.tip(tr2,which(is.na(ind_stt)==T))
      cleaned_data2<-cleaned_data[ind_stt[which(is.na(ind_stt)==F)],]
      cleaned_data2<-cbind(cleaned_data2,"phyl"=tr3$tip.label)
      A1 <- Matrix::Matrix(ape::vcv(tr3), sparse = TRUE)
      colnames(A1) <- rownames(A1) <-tr3$tip.label
      
      ms1<-lmer(species~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+(1|village),
                data = cleaned_data2)
      mss<-summary(ms1)
      ms3<-Almer(sp~age+sex+batch_effect+dna_conc+BMI+sampling_date+bristol+phen+(1|village)+(1|phyl),
                 data = cleaned_data2,A=list(phyl=A1))
      mss2<-summary(ms3)
      #coef(mss)[which(rownames(coef(mss))==colnames(phen_all_raw)[kl]),5]
      
      effect_size_phen_no[i]<-coef(mss)[which(rownames(coef(mss))=="phen"),1]
      pval_phen_no[i]<-coef(mss)[which(rownames(coef(mss))=="phen"),5]
      effect_size_phen[i]<-coef(summary(ms3))[which(rownames(coef(ms3))=="phen"),1]
      pval_phen[i]<-2*pnorm(abs(coef(ms3)[which(rownames(coef(ms3))=="phen"),3]), lower.tail = FALSE)#Pvalue of phenotype
      phyl2[i]<-phylosig(tr3,cleaned_data2$phen, method="lambda", test=FALSE)#Checking for phylogenetic signal
      
    }
    
  }
    }
  }
  }
  
  print(cbind(effect_size_phen_no,pval_phen_no,effect_size_phen,pval_phen,phyl2))
  
}


###Plotting figure 3b


e12<-as.data.frame(output[which((output[,2]<0.05)&output[,5]>0),])#Taking only the significant values in normal-species-phenotype model
colnames(e12)<-c("x1","x1p","y1","y1p")

op<-cor.test(e12$x1,e12$y1,method="spearman")
rho<-op$estimate
pval<-op$p.value
op2<-lm(e12$y1~e12$x1)
coef(op2)
e14<-e12[c(dim(e12)[1]:-1:1),]


pdf('figure_3b.pdf',width=6,height=5)#was width=13,heigh=6(or 4)
ggplot(e14,aes(x=x1,y=y1))+geom_point(alpha=0.5,color="#4575b4",size=3)+geom_abline(intercept=0,slope=1,linetype="dashed", size=1.3)+xlim(-10,10)+ylim(-10,10)+
  geom_segment(x=10,y=10*coef(op2)[2]+coef(op2)[1],xend=-10,yend=-10*coef(op2)[2]+coef(op2)[1],col="red", size=1)+
  theme(legend.position = "none",panel.background = element_blank(),panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),axis.text.x=element_text(size=15,color="black"),axis.text.y=element_text(size=15,color="black"),axis.line.y = element_line(color="black"),axis.line.x = element_line(color="black"))+
  xlab("Without strain-phylogenetic effect")+ylab("With strain-phylogenetic effect")+
  theme(plot.margin=margin(c(1,0,0,0), unit = "cm"))
dev.off()

