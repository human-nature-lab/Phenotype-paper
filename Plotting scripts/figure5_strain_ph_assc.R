
health_ind<-c(1,4,11,6,9,10,12:21,22:28,39,56,57,105:108) #-- 34 signficant species
food_an_ind<-c(58:102) #-- 58 signficant species
soc_eco_ind<-c(109:124,130:184,103:104)

effect_size_phen<-read.csv('esz_strain_with_phylo_4.csv',row.names=1)#effect_size_phen<-read.csv('effect_size_phen_strain_with_phylo_4.csv',row.names=1)
fdr_phen<-read.csv('fdr_with_phylo_4.csv',row.names=1)#fdr_phen<-read.csv('fdr_phen_strain_with_phylo_4.csv',row.names=1)

colnames(effect_size_phen)
colnames(effect_size_phen)[health_ind]
colnames(effect_size_phen)[soc_eco_ind]
colnames(effect_size_phen)[food_an_ind]

colnames(effect_size_phen)[c(39,56,57)]

##Create a long matrix of significant effect sizes for each phenotype (FDR significant)
##Colors can be length of FDR sig speceis in each category like physiological variables

count<-0
count_pos<-0
count_neg<-0

for(i in 1:dim(effect_size_phen)[1]){
  for(j in c(1:10)){
    if(is.na(fdr_phen[i,j])==F){
      if(fdr_phen[i,j]<0.05){
        count<-count+1
        #if(count==1){
        #  count_neg<-count_neg+1
        #  strain_ph_plot_neg<-cbind(effect_size_phen[i,j],"Physiological variables")
        #}else{
        if(effect_size_phen[i,j]>0){
          count_pos<-count_pos+1
          if(count_pos==1){
            strain_ph_plot_pos<-cbind(cbind(effect_size_phen[i,j],"Physiological variables"),fdr_phen[i,j])
          }else{
            strain_ph_plot_pos<-rbind(strain_ph_plot_pos,cbind(cbind(effect_size_phen[i,j],"Physiological variables"),fdr_phen[i,j]))
          }
        }else if(effect_size_phen[i,j]<0){
          count_neg<-count_neg+1
          
          if(count_neg==1){
            strain_ph_plot_neg<-cbind(cbind(effect_size_phen[i,j],"Physiological variables"),fdr_phen[i,j])
          }else{
            strain_ph_plot_neg<-rbind(strain_ph_plot_neg,cbind(cbind(effect_size_phen[i,j],"Physiological variables"),fdr_phen[i,j]))
          }
        }
        #strain_ph_plot<-rbind(strain_ph_plot,cbind(effect_size_phen[i,j],"Physiological variables"))
        #}
        
      }
    }
  }
}

count_all<-count
count_pos_all<-count_pos
count_neg_all<-count_neg

#Chronic condition

count<-0
count_pos<-0
count_neg<-0

for(i in 1:dim(effect_size_phen)[1]){
  for(j in c(12:21)){
    if(is.na(fdr_phen[i,j])==F){
      if(fdr_phen[i,j]<0.05){
        count<-count+1
        #strain_ph_plot<-rbind(strain_ph_plot,cbind(effect_size_phen[i,j],"Chronic conditions"))
        if(effect_size_phen[i,j]>0){
          count_pos<-count_pos+1
          strain_ph_plot_pos<-rbind(strain_ph_plot_pos,cbind(cbind(effect_size_phen[i,j],"Chronic conditions"),fdr_phen[i,j]))
        }else if(effect_size_phen[i,j]<0){
          count_neg<-count_neg+1
          strain_ph_plot_neg<-rbind(strain_ph_plot_neg,cbind(cbind(effect_size_phen[i,j],"Chronic conditions"),fdr_phen[i,j]))
        }
      }
    }
  }
}

count_all<-rbind(count_all,count)
count_pos_all<-rbind(count_pos_all,count_pos)
count_neg_all<-rbind(count_neg_all,count_neg)

#Medication

count<-0
count_pos<-0
count_neg<-0

for(i in 1:dim(effect_size_phen)[1]){
  for(j in c(22:28)){
    if(is.na(fdr_phen[i,j])==F){
      if(fdr_phen[i,j]<0.05){
        count<-count+1
        #strain_ph_plot<-rbind(strain_ph_plot,cbind(effect_size_phen[i,j],"Medication"))
        if(effect_size_phen[i,j]>0){
          count_pos<-count_pos+1
          strain_ph_plot_pos<-rbind(strain_ph_plot_pos,cbind(cbind(effect_size_phen[i,j],"Medication"),fdr_phen[i,j]))
        }else if(effect_size_phen[i,j]<0){
          count_neg<-count_neg+1
          strain_ph_plot_neg<-rbind(strain_ph_plot_neg,cbind(cbind(effect_size_phen[i,j],"Medication"),fdr_phen[i,j]))
        }
      }
    }
  }
}

count_all<-rbind(count_all,count)
count_pos_all<-rbind(count_pos_all,count_pos)
count_neg_all<-rbind(count_neg_all,count_neg)


#Mental health (personality+cognitive+gad7+phq9)

count<-0
count_pos<-0
count_neg<-0

for(i in 1:dim(effect_size_phen)[1]){
  for(j in c(29:57)){
    if(is.na(fdr_phen[i,j])==F){
      if(fdr_phen[i,j]<0.05){
        count<-count+1
        #strain_ph_plot<-rbind(strain_ph_plot,cbind(effect_size_phen[i,j],"Mental health"))
        if(effect_size_phen[i,j]>0){
          count_pos<-count_pos+1
          strain_ph_plot_pos<-rbind(strain_ph_plot_pos,cbind(cbind(effect_size_phen[i,j],"Mental health"),fdr_phen[i,j]))
        }else if(effect_size_phen[i,j]<0){
          count_neg<-count_neg+1
          strain_ph_plot_neg<-rbind(strain_ph_plot_neg,cbind(cbind(effect_size_phen[i,j],"Mental health"),fdr_phen[i,j]))
        }
      }
    }
  }
}

count_all<-rbind(count_all,count)
count_pos_all<-rbind(count_pos_all,count_pos)
count_neg_all<-rbind(count_neg_all,count_neg)


#Unfavorable habits --alcohol cigarette

count<-0
count_pos<-0
count_neg<-0

for(i in 1:dim(effect_size_phen)[1]){
  for(j in c(105:108)){
    if(is.na(fdr_phen[i,j])==F){
      if(fdr_phen[i,j]<0.05){
        count<-count+1
        #strain_ph_plot<-rbind(strain_ph_plot,cbind(effect_size_phen[i,j],"Unfavorable habits"))
        if(effect_size_phen[i,j]>0){
          count_pos<-count_pos+1
          strain_ph_plot_pos<-rbind(strain_ph_plot_pos,cbind(cbind(effect_size_phen[i,j],"Unfavorable habits"),fdr_phen[i,j]))
        }else if(effect_size_phen[i,j]<0){
          count_neg<-count_neg+1
          strain_ph_plot_neg<-rbind(strain_ph_plot_neg,cbind(cbind(effect_size_phen[i,j],"Unfavorable habits"),fdr_phen[i,j]))
        }
      }
    }
  }
}

count_all<-rbind(count_all,count)
count_pos_all<-rbind(count_pos_all,count_pos)
count_neg_all<-rbind(count_neg_all,count_neg)


#Animals

count<-0
count_pos<-0
count_neg<-0

for(i in 1:dim(effect_size_phen)[1]){
  for(j in c(58:83)){
    if(is.na(fdr_phen[i,j])==F){
      if(fdr_phen[i,j]<0.05){
        count<-count+1
        #strain_ph_plot<-rbind(strain_ph_plot,cbind(effect_size_phen[i,j],"Animal exposure"))
        if(effect_size_phen[i,j]>0){
          count_pos<-count_pos+1
          strain_ph_plot_pos<-rbind(strain_ph_plot_pos,cbind(cbind(effect_size_phen[i,j],"Animal exposure"),fdr_phen[i,j]))
        }else if(effect_size_phen[i,j]<0){
          count_neg<-count_neg+1
          strain_ph_plot_neg<-rbind(strain_ph_plot_neg,cbind(cbind(effect_size_phen[i,j],"Animal exposure"),fdr_phen[i,j]))
        }
      }
    }
  }
}

count_all<-rbind(count_all,count)
count_pos_all<-rbind(count_pos_all,count_pos)
count_neg_all<-rbind(count_neg_all,count_neg)


#Food

count<-0
count_pos<-0
count_neg<-0

for(i in 1:dim(effect_size_phen)[1]){
  for(j in c(84:102)){
    if(is.na(fdr_phen[i,j])==F){
      if(fdr_phen[i,j]<0.05){
        count<-count+1
        #strain_ph_plot<-rbind(strain_ph_plot,cbind(effect_size_phen[i,j],"Food"))
        if(effect_size_phen[i,j]>0){
          count_pos<-count_pos+1
          strain_ph_plot_pos<-rbind(strain_ph_plot_pos,cbind(cbind(effect_size_phen[i,j],"Food"),fdr_phen[i,j]))
        }else if(effect_size_phen[i,j]<0){
          count_neg<-count_neg+1
          strain_ph_plot_neg<-rbind(strain_ph_plot_neg,cbind(cbind(effect_size_phen[i,j],"Food"),fdr_phen[i,j]))
        }
      }
    }
  }
}

count_all<-rbind(count_all,count)
count_pos_all<-rbind(count_pos_all,count_pos)
count_neg_all<-rbind(count_neg_all,count_neg)


#Partner

count<-0
count_pos<-0
count_neg<-0

for(i in 1:dim(effect_size_phen)[1]){
  for(j in c(181:184)){
    if(is.na(fdr_phen[i,j])==F){
      if(fdr_phen[i,j]<0.05){
        count<-count+1
        #strain_ph_plot<-rbind(strain_ph_plot,cbind(effect_size_phen[i,j],"Partner"))
        if(effect_size_phen[i,j]>0){
          count_pos<-count_pos+1
          strain_ph_plot_pos<-rbind(strain_ph_plot_pos,cbind(cbind(effect_size_phen[i,j],"Partner"),fdr_phen[i,j]))
        }else if(effect_size_phen[i,j]<0){
          count_neg<-count_neg+1
          strain_ph_plot_neg<-rbind(strain_ph_plot_neg,cbind(cbind(effect_size_phen[i,j],"Partner"),fdr_phen[i,j]))
        }
      }
    }
  }
}

count_all<-rbind(count_all,count)
count_pos_all<-rbind(count_pos_all,count_pos)
count_neg_all<-rbind(count_neg_all,count_neg)


#Friendship

count<-0
count_pos<-0
count_neg<-0

for(i in 1:dim(effect_size_phen)[1]){
  for(j in c(112:115)){
    if(is.na(fdr_phen[i,j])==F){
      if(fdr_phen[i,j]<0.05){
        count<-count+1
        #strain_ph_plot<-rbind(strain_ph_plot,cbind(effect_size_phen[i,j],"Friendship network"))
        if(effect_size_phen[i,j]>0){
          count_pos<-count_pos+1
          strain_ph_plot_pos<-rbind(strain_ph_plot_pos,cbind(cbind(effect_size_phen[i,j],"Friendship network"),fdr_phen[i,j]))
        }else if(effect_size_phen[i,j]<0){
          count_neg<-count_neg+1
          strain_ph_plot_neg<-rbind(strain_ph_plot_neg,cbind(cbind(effect_size_phen[i,j],"Friendship network"),fdr_phen[i,j]))
        }
      }
    }
  }
}

count_all<-rbind(count_all,count)
count_pos_all<-rbind(count_pos_all,count_pos)
count_neg_all<-rbind(count_neg_all,count_neg)

#Family

count<-0
count_pos<-0
count_neg<-0

for(i in 1:dim(effect_size_phen)[1]){
  for(j in c(117:120)){
    if(is.na(fdr_phen[i,j])==F){
      if(fdr_phen[i,j]<0.05){
        count<-count+1
        #strain_ph_plot<-rbind(strain_ph_plot,cbind(effect_size_phen[i,j],"Familial network"))
        if(effect_size_phen[i,j]>0){
          count_pos<-count_pos+1
          strain_ph_plot_pos<-rbind(strain_ph_plot_pos,cbind(cbind(effect_size_phen[i,j],"Familial network"),fdr_phen[i,j]))
        }else if(effect_size_phen[i,j]<0){
          count_neg<-count_neg+1
          strain_ph_plot_neg<-rbind(strain_ph_plot_neg,cbind(cbind(effect_size_phen[i,j],"Familial network"),fdr_phen[i,j]))
        }
      }
    }
  }
}

count_all<-rbind(count_all,count)
count_pos_all<-rbind(count_pos_all,count_pos)
count_neg_all<-rbind(count_neg_all,count_neg)


#All ties

count<-0
count_pos<-0
count_neg<-0

for(i in 1:dim(effect_size_phen)[1]){
  for(j in c(121:125)){
    if(is.na(fdr_phen[i,j])==F){
      if(fdr_phen[i,j]<0.05){
        count<-count+1
        #strain_ph_plot<-rbind(strain_ph_plot,cbind(effect_size_phen[i,j],"Full network"))
        if(effect_size_phen[i,j]>0){
          count_pos<-count_pos+1
          strain_ph_plot_pos<-rbind(strain_ph_plot_pos,cbind(cbind(effect_size_phen[i,j],"Full network"),fdr_phen[i,j]))
        }else if(effect_size_phen[i,j]<0){
          count_neg<-count_neg+1
          strain_ph_plot_neg<-rbind(strain_ph_plot_neg,cbind(cbind(effect_size_phen[i,j],"Full network"),fdr_phen[i,j]))
        }
      }
    }
  }
}

count_all<-rbind(count_all,count)
count_pos_all<-rbind(count_pos_all,count_pos)
count_neg_all<-rbind(count_neg_all,count_neg)


#Economic factors

count<-0
count_pos<-0
count_neg<-0

for(i in 1:dim(effect_size_phen)[1]){
  for(j in c(103:104,127:172)){
    if(is.na(fdr_phen[i,j])==F){
      if(fdr_phen[i,j]<0.05){
        count<-count+1
        #strain_ph_plot<-rbind(strain_ph_plot,cbind(effect_size_phen[i,j],"Economic factors"))
        if(effect_size_phen[i,j]>0){
          count_pos<-count_pos+1
          strain_ph_plot_pos<-rbind(strain_ph_plot_pos,cbind(cbind(effect_size_phen[i,j],"Economic factors"),fdr_phen[i,j]))
        }else if(effect_size_phen[i,j]<0){
          count_neg<-count_neg+1
          strain_ph_plot_neg<-rbind(strain_ph_plot_neg,cbind(cbind(effect_size_phen[i,j],"Economic factors"),fdr_phen[i,j]))
        }
      }
    }
  }
}

count_all<-rbind(count_all,count)
count_pos_all<-rbind(count_pos_all,count_pos)
count_neg_all<-rbind(count_neg_all,count_neg)


#Water sources

count<-0
count_pos<-0
count_neg<-0

for(i in 1:dim(effect_size_phen)[1]){
  for(j in c(173:180)){
    if(is.na(fdr_phen[i,j])==F){
      if(fdr_phen[i,j]<0.05){
        count<-count+1
        #strain_ph_plot<-rbind(strain_ph_plot,cbind(effect_size_phen[i,j],"Water source"))
        if(effect_size_phen[i,j]>0){
          count_pos<-count_pos+1
          strain_ph_plot_pos<-rbind(strain_ph_plot_pos,cbind(cbind(effect_size_phen[i,j],"Water source"),fdr_phen[i,j]))
        }else if(effect_size_phen[i,j]<0){
          count_neg<-count_neg+1
          strain_ph_plot_neg<-rbind(strain_ph_plot_neg,cbind(cbind(effect_size_phen[i,j],"Water source"),fdr_phen[i,j]))
        }
      }
    }
  }
}

count_all<-rbind(count_all,count)
count_pos_all<-rbind(count_pos_all,count_pos)
count_neg_all<-rbind(count_neg_all,count_neg)

strain_ph_plot_neg<-as.data.frame(strain_ph_plot_neg)
strain_ph_plot_pos<-as.data.frame(strain_ph_plot_pos)
colnames(strain_ph_plot_pos)<-c("Effect size","Phenotype")
colnames(strain_ph_plot_neg)<-c("Effect size","Phenotype")

ord<-unique(strain_ph_plot_neg[,2])

strain_ph_plot_neg[,2] <- factor(strain_ph_plot_neg[,2], levels = ord[length(ord):-1:1])
library(ggplot2)
ggplot(strain_ph_plot_neg, aes(x=Phenotype, y=as.numeric(as.character(`Effect size`)),group=Phenotype)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE)+labs(y="Strain-phenotype effect sizes")+
  theme(text = element_text(size=15))+ scale_x_discrete(breaks=seq(1,13,1))#+ylim(0,0.35)#,"#ff7f00","#cab2d6","#6a3d9a","#ffff99"#x="Relationship",



ggplot(strain_ph_plot_neg, aes(x=Phenotype, y=as.numeric(as.character(`Effect size`)),group=Phenotype)) + #do only for with mbiome and without separately
  geom_boxplot(notch=FALSE)+coord_flip()

#strain_ph_plot_neg$`Effect size`<-as.numeric(as.character(strain_ph_plot_neg$`Effect size`))

strain_ph_plot_neg2<-strain_ph_plot_neg[which((strain_ph_plot_neg$Phenotype!="Physiological variables")&(strain_ph_plot_neg$Phenotype!="Unfavorable habits")&(strain_ph_plot_neg$Phenotype!="Partner")),]

ggplot(strain_ph_plot_neg2, aes(x=Phenotype, y=as.numeric(as.character(`Effect size`)),group=Phenotype))+ylab("Negative effect size") +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(shape = 19,
              color = "steelblue",
              position = position_jitter(width = 0.21)) +
  theme(text = element_text(size=20))+coord_flip()#+scale_x_discrete(limits = as.character(round(seq("-2", "0", by = 1),1)))

#temp<-cbind(NA,"Partner")
#colnames(temp)<-c("Effect size","Phenotype")
#strain_ph_plot_pos<-rbind(strain_ph_plot_pos,temp)
strain_ph_plot_pos[,2] <- factor(strain_ph_plot_pos[,2], levels = ord[length(ord):-1:1])

strain_ph_plot_pos2<-strain_ph_plot_pos[which((strain_ph_plot_pos$Phenotype!="Physiological variables")&(strain_ph_plot_pos$Phenotype!="Unfavorable habits")&(strain_ph_plot_pos$Phenotype!="Partner")),]

ggplot(strain_ph_plot_pos2, aes(x=Phenotype, y=as.numeric(as.character(`Effect size`)),group=Phenotype))+ylab("Effect size") +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(shape = 19,
              color = "steelblue",
              position = position_jitter(width = 0.21)) +
  theme(text = element_text(size=20))+coord_flip()+scale_x_discrete()


strain_ph_plot_neg3<-cbind(strain_ph_plot_neg2,array("Negative",dim=c(dim(strain_ph_plot_neg2)[1],1)))
colnames(strain_ph_plot_neg3)[4]<-"type"
strain_ph_plot_pos3<-cbind(strain_ph_plot_pos2,array("Positive",dim=c(dim(strain_ph_plot_pos2)[1],1)))
colnames(strain_ph_plot_pos3)[4]<-"type"

strain_ph_plot_all3<-rbind(strain_ph_plot_neg3,strain_ph_plot_pos3)
strain_ph_plot_all3<-strain_ph_plot_all3[which(as.numeric(as.character(strain_ph_plot_all3[,1]))<2.5),]

ggplot(strain_ph_plot_all3, aes(x=Phenotype, y=as.numeric(as.character(`Effect size`)),group=Phenotype))+ylab("Effect size") +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(shape = 19,
              color = "steelblue",
              position = position_jitter(width = 0.21)) +
  theme(text = element_text(size=20),panel.spacing.x =unit(-0.75,"lines"))+coord_flip()+ facet_wrap(~type,scales="free_x")

#
strain_ph_plot_neg4<-cbind(cbind(cbind(as.numeric(as.character(strain_ph_plot_neg[,1])),as.character(strain_ph_plot_neg[,2])),as.character(strain_ph_plot_neg[,3])),array("Negative",dim=c(dim(strain_ph_plot_neg)[1],1)))
colnames(strain_ph_plot_neg4)[4]<-"type"
strain_ph_plot_pos4<-cbind(cbind(cbind(as.numeric(as.character(strain_ph_plot_pos[,1])),as.character(strain_ph_plot_pos[,2])),as.character(strain_ph_plot_pos[,3])),array("Positive",dim=c(dim(strain_ph_plot_pos)[1],1)))
colnames(strain_ph_plot_pos4)[4]<-"type"

strain_ph_plot_all4<-rbind(strain_ph_plot_neg4,strain_ph_plot_pos4)
strain_ph_plot_all4<-strain_ph_plot_all4[which((as.numeric(as.character(strain_ph_plot_all4[,1]))<1.5)&(as.numeric(as.character(strain_ph_plot_all4[,1]))>c(-1.5))),]

strain_ph_plot_all4[279,1]<-1.4085
strain_ph_plot_all4<-as.data.frame(strain_ph_plot_all4)
colnames(strain_ph_plot_all4)<-c("Effect size","Phenotype","fdr","type")
strain_ph_plot_all4[,2] <- factor(strain_ph_plot_all4[,2], levels = as.character(unique(strain_ph_plot_all4[,2])))
ord<-unique(strain_ph_plot_all4[,2])
strain_ph_plot_all4[,2] <- factor(strain_ph_plot_all4[,2], levels = ord[length(ord):-1:1])

ggplot(strain_ph_plot_all4, aes(x=Phenotype, y=as.numeric(as.character(`Effect size`)),group=Phenotype))+ylab("Effect size") +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(shape = 19,
              color = "steelblue",
              position = position_jitter(width = 0.21)) +
  theme(text = element_text(size=20,color="black"),axis.text.x=element_text(color="black"),axis.text.y=element_text(color="black"),panel.spacing.x =unit(-0.5,"lines"), strip.text = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+
  coord_flip()+ facet_wrap(~as.character(type),scales="free_x")+#+scale_y_discrete(labels=c("-1","-0.5","0","0.5","1"))
  xlab("")

ggplot(strain_ph_plot_all4, aes(x=Phenotype, y=as.numeric(as.character(`Effect size`)),group=Phenotype))+ylab("Effect size") +
  #geom_boxplot(outlier.alpha = 0) +
  geom_jitter(shape = 19,
              color = "steelblue",
              position = position_jitter(width = 0.21)) +
  theme(text = element_text(size=20,color="black"),axis.text.x=element_text(color="black"),axis.text.y=element_text(color="black"),axis.line.x = element_line(color="black"),axis.line.y=element_line(color="black"),panel.spacing.x =unit(-0.5,"lines"), strip.text = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())+
  coord_flip()+ facet_wrap(~as.character(type),scales="free_x")+#+scale_y_discrete(labels=c("-1","-0.5","0","0.5","1"))
  xlab("")+
  theme(plot.margin=margin(c(0,3,0,0), unit = "cm"))

#+
#Save 4.00 x 8.00

##Finding how many significant associations

effect_size_phen_sig<-effect_size_phen

for(i in 1:dim(effect_size_phen)[1]){
  for(j in 1:dim(effect_size_phen)[2]){
    if(is.na(fdr_phen[i,j])==F){
      if(fdr_phen[i,j]>0.05){
        effect_size_phen_sig[i,j]<-NA
      }
    }
  }
}

effect_size_phen_sig2<-effect_size_phen_sig[,103:183]#Food&a - 58:102, health 1:57,socioeco- 103:183
for(i in 1:dim(effect_size_phen_sig2)[1]){
  if(i==1){
    na_score<-length(which(is.na(effect_size_phen_sig2[i,])==F))
  }else{
    na_score<-rbind(na_score,length(which(is.na(effect_size_phen_sig2[i,])==F)))
  }
}
length(which(na_score>0))


  #geom_segment(x=0,xend=1.5,y=0,yend=0,col="#a6cee3",size=0.2)+#Health
  #geom_rect(xmin=c(-1.5),xmax=1.5,ymin=0,ymax=(-5.5),color="#a6cee3",fill="#a6cee3")+
  #annotate("text",x=(2),y=(-3),label="Health",color="black",fontface="bold",hjust=0.75)+
  #annotate("text",x=(dim(effect_size_phen_sig_plot2)[2]+5.2),y=(dim(effect_size_phen_sig_plot2)[1]-16.7),label="conditions",color="black",hjust=0.6,fontface="bold")+
  #geom_rect(xmin=c(-1.5),xmax=1.5,ymin=0,ymax=(-5.5),color="#a6cee3",fill="#a6cee3")


##Without friendship

strain_ph_plot_all5<-strain_ph_plot_all4[which(strain_ph_plot_all4$Phenotype!="Friendship network"),]

cmm<-array("black",dim=c(dim(strain_ph_plot_all5)[1]))
for(i in 1:length(cmm)){
  if((as.numeric(as.character(strain_ph_plot_all5$fdr[i]))<0.05)&(as.numeric(as.character(strain_ph_plot_all5$fdr[i]))>=0.01)){
    cmm[i]<-"steelblue"
  }else if((as.numeric(as.character(strain_ph_plot_all5$fdr[i]))<0.01)&(as.numeric(as.character(strain_ph_plot_all5$fdr[i]))>=0.001)){
    cmm[i]<-"red"
  }else if((as.numeric(as.character(strain_ph_plot_all5$fdr[i]))<0.001)&(as.numeric(as.character(strain_ph_plot_all5$fdr[i]))>=0.0001)){
    cmm[i]<-"green"
  }else if((as.numeric(as.character(strain_ph_plot_all5$fdr[i]))<0.0001)){
    cmm[i]<-"purple"
  }
}

#4.00 x 8.00
ggplot(strain_ph_plot_all5, aes(x=Phenotype, y=as.numeric(as.character(`Effect size`)),group=Phenotype))+ylab("Effect size") +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(shape = 19,
              color="steelblue",#color = cmm,
              position = position_jitter(width = 0.21)) +
  theme(text = element_text(size=20,color="black"),axis.text.x=element_text(color="black"),axis.text.y=element_text(color="black"),panel.spacing.x =unit(-0.5,"lines"), strip.text = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line.x = element_line(color="black"))+
  coord_flip()+ facet_wrap(~as.character(type),scales="free_x")+#+scale_y_discrete(labels=c("-1","-0.5","0","0.5","1"))
  xlab("")


kmn<-as.character(unique(strain_ph_plot_all5[,2]))

ln<-strain_ph_plot_all5[which(as.character(strain_ph_plot_all5[,2])%in%kmn[1:5]),]

mean(as.numeric(as.character(ln[which(as.numeric(as.character(ln[,1]))>0),1])))
mean(as.numeric(as.character(ln[which(as.numeric(as.character(ln[,1]))<0),1])))
length(ln[which(as.numeric(as.character(ln[,3]))<0.0001),1])

ln<-strain_ph_plot_all5[which(as.character(strain_ph_plot_all5[,2])%in%kmn[6:7]),]
mean(as.numeric(as.character(ln[which(as.numeric(as.character(ln[,1]))>0),1])))
mean(as.numeric(as.character(ln[which(as.numeric(as.character(ln[,1]))<0),1])))
length(ln[which(as.numeric(as.character(ln[,3]))<0.0001),1])

ln<-strain_ph_plot_all5[which(as.character(strain_ph_plot_all5[,2])%in%kmn[8:12]),]
mean(as.numeric(as.character(ln[which(as.numeric(as.character(ln[,1]))>0),1])))
mean(as.numeric(as.character(ln[which(as.numeric(as.character(ln[,1]))<0),1])))
length(ln[which(as.numeric(as.character(ln[,3]))<0.0001),1])


############################################################################################################################
##It is fine -- all points are there

###OLD PLOTTING

ggplot(strain_ph_plot_neg, aes(x=Phenotype, y=as.numeric(as.character(`Effect size`)),group=Phenotype))+ylab("Negative effect size") +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(shape = 19,
              color = "steelblue",
              position = position_jitter(width = 0.21)) +
  theme(text = element_text(size=20))+coord_flip()

temp<-cbind(NA,"Partner")
colnames(temp)<-c("Effect size","Phenotype")
strain_ph_plot_pos<-rbind(strain_ph_plot_pos,temp)
strain_ph_plot_pos[,2] <- factor(strain_ph_plot_pos[,2], levels = ord[length(ord):-1:1])

strain_ph_plot_pos[34,1]<-as.character(as.numeric(strain_ph_plot_pos[34,1])/3.5)
strain_ph_plot_pos[1,1]<-as.character(as.numeric(strain_ph_plot_pos[1,1])/4)

ggplot(strain_ph_plot_pos, aes(x=Phenotype, y=as.numeric(as.character(`Effect size`)),group=Phenotype))+ylab("Positive effect size") +
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(shape = 19,
              color = "steelblue",
              position = position_jitter(width = 0.21)) +
  theme(text = element_text(size=20))+coord_flip()


#ggplot(strain_ph_plot_neg, aes(x=Phenotype, y=as.numeric(as.character(`Effect size`)),group=Phenotype)) +
#library(ggpubr)
#ggboxplot(strain_ph_plot_neg, x = "Phenotype", y = as.numeric(as.character("Effect size")),
#          color = "Phenotype", palette = "jco",
#          add = "jitter")

###############################

##Just to check number of significant species in all health

pval_chk<-array(0,dim=c(dim(pval_phen)[1],dim(pval_phen)[2]))
for(i in 1:dim(pval_chk)[1]){
  for(j in 1:dim(pval_chk)[2]){
    if(is.na(pval_phen[i,j])==F){
      if(pval_phen[i,j]<0.05){
        pval_chk[i,j]<-1
      }
    }
  }
}
fdr_chk<-array(0,dim=c(dim(fdr_phen)[1],dim(fdr_phen)[2]))
for(i in 1:dim(fdr_chk)[1]){
  for(j in 1:dim(fdr_chk)[2]){
    if(is.na(fdr_phen[i,j])==F){
      if(fdr_phen[i,j]<0.05){
        fdr_chk[i,j]<-1
      }
    }
  }
}
colnames(effect_size_phen)
pval_chk_lmer_all<-rowSums(pval_chk)
length(which(pval_chk_lmer_all>3))

fdr_chk_lmer_all<-rowSums(fdr_chk[,food_an_ind])#c(1,4:11,16:17,19:20,22:25)#44:83(personality types discrete)#
length(which(fdr_chk_lmer_all>0))

effect_size_phen_sig<-effect_size_phen[which(fdr_chk_lmer_all>5),food_an_ind]
pval_phen_sig<-pval_chk[which(fdr_chk_lmer_all>5),food_an_ind]
fdr_phen_sig<-fdr_chk[which(fdr_chk_lmer_all>5),food_an_ind]

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
rownames(effect_size_phen_sig_plot)<-colnames(effect_size_phen)[health_ind]
