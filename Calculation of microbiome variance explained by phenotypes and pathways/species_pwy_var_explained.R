library(vegan)

phen_all_use<-read.csv('phen_all_use6_all.csv',row.names = 1,header=T)#phen_all_use4_all.csv
data_transformed<-read.csv('mb_samp_1188_3.csv',row.names=1,header=T)#read.csv('data_transformed_1187.csv',row.names = 1,header=T)
data_transformed<-as.matrix.data.frame(data_transformed)
#ind_t<-read.csv('ind_t.csv')
#ind_t<-ind_t[,1]
#data_transformed<-data_transformed[ind_t,]
st<-read.csv('bristol.csv',row.names = 1)
st<-as.matrix.data.frame(st)
phen_all_rct_use<-read.csv('phen_all_rct_use22.csv',row.names = 1,header=T)

for(i in 1:dim(phen_all_use)[1]){
  if(is.na(phen_all_use[i,35])==F){
    if(phen_all_use[i,35]=="Dont_Know"){
      phen_all_use[i,35]<-NA
    }
  }
  if(is.na(phen_all_use[i,36])==F){
    if(phen_all_use[i,36]=="Dont_Know"){
      phen_all_use[i,36]<-NA
    }
  }
  if(is.na(phen_all_use[i,37])==F){
    if(phen_all_use[i,37]=="Dont_Know"){
      phen_all_use[i,37]<-NA
    }
  }
  if(is.na(phen_all_use[i,71])==F){
    if(phen_all_use[i,71]=="Dont_Know"){
      phen_all_use[i,71]<-NA
    }
  }
  if(is.na(phen_all_use[i,72])==F){
    if(phen_all_use[i,72]=="Dont_Know"){
      phen_all_use[i,72]<-NA
    }
  }
  if(is.na(phen_all_use[i,73])==F){
    if(phen_all_use[i,73]=="Dont_Know"){
      phen_all_use[i,73]<-NA
    }
  }
  if(is.na(phen_all_use[i,74])==F){
    if(phen_all_use[i,74]=="Dont_Know"){
      phen_all_use[i,74]<-NA
    }
  }
}

phen_all_use<-cbind(phen_all_use,st)
colnames(phen_all_use)[dim(phen_all_use)[2]]<-"Bristol"
phen_map<-as.numeric(as.character(phen_all_use$mb_d0400))+(as.numeric(as.character(phen_all_use$mb_d0300))-as.numeric(as.character(phen_all_use$mb_d0400)))/3
phen_all_use<-cbind(phen_all_use,"MAP"=phen_map)
phen_all_k<-array(NA,dim=c(dim(phen_all_use)[1],dim(phen_all_use)[2]))
colnames(phen_all_k)<-colnames(phen_all_use)

##Continuous/Discrete

for(kl in c(32,35,36,381,6,40,4,382,71:74,37,153,182,208,364,203:207,360:363,365:369,228:230,372,306:309,317:359,375:380,383)){
  for(k in 1:dim(phen_all_use)[1]){
    phen_all_k[k,kl]<-as.numeric(as.character(phen_all_use[k,kl]))
  }
}




##Discrete
#Skipping overall health -- since chronic conditions included
for(kl in c(143:152)){
  for(k in 1:dim(phen_all_use)[1]){
    if(is.na(phen_all_use[k,kl])==F){
      if(phen_all_use[k,kl]=="Disagree strongly"){
        phen_all_k[k,kl]<-1
      }else if(phen_all_use[k,kl]=="Disagree a little"){
        phen_all_k[k,kl]<-2
      }else if(phen_all_use[k,kl]=="Neither agree nor disagree"){
        phen_all_k[k,kl]<-3
      }else if(phen_all_use[k,kl]=="Agree a little"){
        phen_all_k[k,kl]<-4
      }else if(phen_all_use[k,kl]=="Agree strongly"){
        phen_all_k[k,kl]<-5
      }
    }
    if(is.na(phen_all_use[k,kl])==T){
      phen_all_k[k,kl]<-NA
    }
  }
}

for(kl in c(75:93)){
  for(k in 1:dim(phen_all_use)[1]){
    if(is.na(phen_all_use[k,kl])==F){
      if(phen_all_use[k,kl]=="Never/rarely"){
        phen_all_k[k,kl]<-1/50
      }else if(phen_all_use[k,kl]=="A few days per month"){
        phen_all_k[k,kl]<-1/10
      }else if(phen_all_use[k,kl]=="A few days per week"){
        phen_all_k[k,kl]<-4/7
      }else if(phen_all_use[k,kl]=="Every day"){
        phen_all_k[k,kl]<-1
      }
    }
    if(is.na(phen_all_use[k,kl])==T){
      phen_all_k[k,kl]<-NA
    }
  }
}

for(kl in 181){
  for(k in 1:dim(phen_all_use)[1]){
    if(is.na(phen_all_use[k,kl])==F){
      if(phen_all_use[k,kl]=="Never/rarely"){
        phen_all_k[k,kl]<-1/100
      }else if(phen_all_use[k,kl]=="At least once a month"){
        phen_all_k[k,kl]<-1/30
      }else if(phen_all_use[k,kl]=="At least once a week"){
        phen_all_k[k,kl]<-1/7
      }else if(phen_all_use[k,kl]=="Every day"){
        phen_all_k[k,kl]<-1
      }
    }
    if(is.na(phen_all_use[k,kl])==T){
      phen_all_k[k,kl]<-NA
    }
  }
}



ind_rct<-match(rownames(phen_all_use),rownames(phen_all_rct_use))
phen_all_k<-cbind(phen_all_k,"Alcohol(daily)"=array(NA,dim=c(dim(phen_all_k)[1],1)))
for(kl in 8){
  for(k in 1:dim(phen_all_use)[1]){
    if(is.na(ind_rct[k])==F){
    if(is.na(phen_all_rct_use[ind_rct[k],kl])==F){
      if(phen_all_rct_use[ind_rct[k],kl]=="1 or 2"){
        phen_all_k[k,dim(phen_all_k)[2]]<-1.5
      }else if(phen_all_rct_use[ind_rct[k],kl]=="3 or 4"){
        phen_all_k[k,dim(phen_all_k)[2]]<-3.5
      }else if(phen_all_rct_use[ind_rct[k],kl]=="5 or 6"){
        phen_all_k[k,dim(phen_all_k)[2]]<-5.5
      }else if(phen_all_rct_use[ind_rct[k],kl]=="7 to 9"){
        phen_all_k[k,dim(phen_all_k)[2]]<-8
      }else if(phen_all_rct_use[ind_rct[k],kl]=="10 or more"){
        phen_all_k[k,dim(phen_all_k)[2]]<-12
      }
    }
    if(is.na(phen_all_rct_use[ind_rct[k],kl])==T){
      phen_all_k[k,dim(phen_all_k)[2]]<-0
    }
  }
  }
}
colnames(phen_all_k)[dim(phen_all_k)[2]]<-"Alcohol(daily)"

phen_all_k<-cbind(phen_all_k,"Cigarette use"=array(NA,dim=c(dim(phen_all_k)[1],1)))
for(kl in 10){
  for(k in 1:dim(phen_all_use)[1]){
    if(is.na(ind_rct[k])==F){
      if(is.na(phen_all_rct_use[ind_rct[k],kl])==F){
        if(phen_all_rct_use[ind_rct[k],kl]=="No"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }else if(phen_all_rct_use[ind_rct[k],kl]=="Yes"){
          phen_all_k[k,dim(phen_all_k)[2]]<-1
        }
      }
      if(is.na(phen_all_rct_use[ind_rct[k],kl])==T){
        phen_all_k[k,dim(phen_all_k)[2]]<-NA
      }
    }
  }
}
colnames(phen_all_k)[dim(phen_all_k)[2]]<-"Cigarette use"

phen_all_k<-cbind(phen_all_k,"Cigarette frequency"=array(NA,dim=c(dim(phen_all_k)[1],1)))
for(kl in 11){
  for(k in 1:dim(phen_all_use)[1]){
    if(is.na(ind_rct[k])==F){
      if(is.na(phen_all_rct_use[ind_rct[k],kl])==F){
          phen_all_k[k,dim(phen_all_k)[2]]<-as.numeric(phen_all_rct_use[ind_rct[k],kl])
      }
      if(is.na(phen_all_rct_use[ind_rct[k],kl])==T){
        phen_all_k[k,dim(phen_all_k)[2]]<-NA
      }
    }
  }
}
colnames(phen_all_k)[dim(phen_all_k)[2]]<-"Cigarette frequency"

phen_all_k<-cbind(phen_all_k,"Education(primary)"=array(NA,dim=c(dim(phen_all_k)[1],1)))
for(kl in 5){
  for(k in 1:dim(phen_all_use)[1]){
    if(is.na(ind_rct[k])==F){
      if(is.na(phen_all_rct_use[ind_rct[k],kl])==F){
        if(phen_all_rct_use[ind_rct[k],kl]=="1st grade"){
          phen_all_k[k,dim(phen_all_k)[2]]<-1
        }else if(phen_all_rct_use[ind_rct[k],kl]=="2nd grade"){
          phen_all_k[k,dim(phen_all_k)[2]]<-1
        }else if(phen_all_rct_use[ind_rct[k],kl]=="3rd grade"){
          phen_all_k[k,dim(phen_all_k)[2]]<-1
        }else if(phen_all_rct_use[ind_rct[k],kl]=="4th grade"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }else if(phen_all_rct_use[ind_rct[k],kl]=="5th grade"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }else if(phen_all_rct_use[ind_rct[k],kl]=="6th grade"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }else if(phen_all_rct_use[ind_rct[k],kl]=="Some secondary"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }else if(phen_all_rct_use[ind_rct[k],kl]=="Secondary"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }else if(phen_all_rct_use[ind_rct[k],kl]=="More than secondary"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }else if(phen_all_rct_use[ind_rct[k],kl]=="Have not completed any type of school"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }
      }
      if(is.na(phen_all_rct_use[ind_rct[k],kl])==T){
        phen_all_k[k,dim(phen_all_k)[2]]<-NA
      }
    }
  }
}
colnames(phen_all_k)[dim(phen_all_k)[2]]<-"Education(primary)"

phen_all_k<-cbind(phen_all_k,"Education(Middle)"=array(NA,dim=c(dim(phen_all_k)[1],1)))
for(kl in 5){
  for(k in 1:dim(phen_all_use)[1]){
    if(is.na(ind_rct[k])==F){
      if(is.na(phen_all_rct_use[ind_rct[k],kl])==F){
        if(phen_all_rct_use[ind_rct[k],kl]=="1st grade"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }else if(phen_all_rct_use[ind_rct[k],kl]=="2nd grade"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }else if(phen_all_rct_use[ind_rct[k],kl]=="3rd grade"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }else if(phen_all_rct_use[ind_rct[k],kl]=="4th grade"){
          phen_all_k[k,dim(phen_all_k)[2]]<-1
        }else if(phen_all_rct_use[ind_rct[k],kl]=="5th grade"){
          phen_all_k[k,dim(phen_all_k)[2]]<-1
        }else if(phen_all_rct_use[ind_rct[k],kl]=="6th grade"){
          phen_all_k[k,dim(phen_all_k)[2]]<-1
        }else if(phen_all_rct_use[ind_rct[k],kl]=="Some secondary"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }else if(phen_all_rct_use[ind_rct[k],kl]=="Secondary"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }else if(phen_all_rct_use[ind_rct[k],kl]=="More than secondary"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }else if(phen_all_rct_use[ind_rct[k],kl]=="Have not completed any type of school"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }
      }
      if(is.na(phen_all_rct_use[ind_rct[k],kl])==T){
        phen_all_k[k,dim(phen_all_k)[2]]<-NA
      }
    }
  }
}
colnames(phen_all_k)[dim(phen_all_k)[2]]<-"Education(middle)"

phen_all_k<-cbind(phen_all_k,"Education(Secondary)"=array(NA,dim=c(dim(phen_all_k)[1],1)))
for(kl in 5){
  for(k in 1:dim(phen_all_use)[1]){
    if(is.na(ind_rct[k])==F){
      if(is.na(phen_all_rct_use[ind_rct[k],kl])==F){
        if(phen_all_rct_use[ind_rct[k],kl]=="1st grade"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }else if(phen_all_rct_use[ind_rct[k],kl]=="2nd grade"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }else if(phen_all_rct_use[ind_rct[k],kl]=="3rd grade"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }else if(phen_all_rct_use[ind_rct[k],kl]=="4th grade"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }else if(phen_all_rct_use[ind_rct[k],kl]=="5th grade"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }else if(phen_all_rct_use[ind_rct[k],kl]=="6th grade"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }else if(phen_all_rct_use[ind_rct[k],kl]=="Some secondary"){
          phen_all_k[k,dim(phen_all_k)[2]]<-1
        }else if(phen_all_rct_use[ind_rct[k],kl]=="Secondary"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }else if(phen_all_rct_use[ind_rct[k],kl]=="More than secondary"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }else if(phen_all_rct_use[ind_rct[k],kl]=="Have not completed any type of school"){
          phen_all_k[k,dim(phen_all_k)[2]]<-0
        }
      }
      if(is.na(phen_all_rct_use[ind_rct[k],kl])==T){
        phen_all_k[k,dim(phen_all_k)[2]]<-NA
      }
    }
  }
}
colnames(phen_all_k)[dim(phen_all_k)[2]]<-"Education(secondary)"

phen_all_k<-cbind(phen_all_k,"Washing hands"=array(NA,dim=c(dim(phen_all_k)[1],1)))
for(kl in 12){
  for(k in 1:dim(phen_all_use)[1]){
    if(is.na(ind_rct[k])==F){
      if(is.na(phen_all_rct_use[ind_rct[k],kl])==F){
          phen_all_k[k,dim(phen_all_k)[2]]<-1
      }
      if(is.na(phen_all_rct_use[ind_rct[k],kl])==T){
        phen_all_k[k,dim(phen_all_k)[2]]<-0
      }
    }
  }
}
colnames(phen_all_k)[dim(phen_all_k)[2]]<-"Washing hands"

phen_all_k<-cbind(phen_all_k,"GAD7(mild)"=array(NA,dim=c(dim(phen_all_k)[1],1)))
for(k in 1:dim(phen_all_use)[1]){
    phen_all_k[k,dim(phen_all_k)[2]]<-ifelse(phen_all_use[k,191]=="mild",1,0)
}
colnames(phen_all_k)[dim(phen_all_k)[2]]<-"GAD7(mild)"
phen_all_k<-cbind(phen_all_k,"GAD7(moderate)"=array(NA,dim=c(dim(phen_all_k)[1],1)))
for(k in 1:dim(phen_all_use)[1]){
  phen_all_k[k,dim(phen_all_k)[2]]<-ifelse(phen_all_use[k,191]=="moderate",1,0)
}
colnames(phen_all_k)[dim(phen_all_k)[2]]<-"GAD7(moderate)"
phen_all_k<-cbind(phen_all_k,"GAD7(severe)"=array(NA,dim=c(dim(phen_all_k)[1],1)))
for(k in 1:dim(phen_all_use)[1]){
  phen_all_k[k,dim(phen_all_k)[2]]<-ifelse(phen_all_use[k,191]=="severe",1,0)
}
colnames(phen_all_k)[dim(phen_all_k)[2]]<-"GAD7(severe)"
phen_all_k<-cbind(phen_all_k,"PHQ9(mild)"=array(NA,dim=c(dim(phen_all_k)[1],1)))
for(k in 1:dim(phen_all_use)[1]){
  phen_all_k[k,dim(phen_all_k)[2]]<-ifelse(phen_all_use[k,202]=="mild",1,0)
}
colnames(phen_all_k)[dim(phen_all_k)[2]]<-"PHQ9(mild)"
phen_all_k<-cbind(phen_all_k,"PHQ9(moderate)"=array(NA,dim=c(dim(phen_all_k)[1],1)))
for(k in 1:dim(phen_all_use)[1]){
  phen_all_k[k,dim(phen_all_k)[2]]<-ifelse(phen_all_use[k,202]=="moderate",1,0)
}
colnames(phen_all_k)[dim(phen_all_k)[2]]<-"PHQ9(moderate)"
phen_all_k<-cbind(phen_all_k,"PHQ9(severe)"=array(NA,dim=c(dim(phen_all_k)[1],1)))
for(k in 1:dim(phen_all_use)[1]){
  phen_all_k[k,dim(phen_all_k)[2]]<-ifelse(phen_all_use[k,202]=="severe",1,0)
}
colnames(phen_all_k)[dim(phen_all_k)[2]]<-"PHQ9(severe)"

##Dichotomous

for(kl in c(16,227)){
  for(k in 1:dim(phen_all_use)[1]){
    phen_all_k[k,kl]<-ifelse(phen_all_use[k,kl]=="No",0,1)
  }
}

for(kl in c(20:30,136:142,155:180,310:316)){
  for(k in 1:dim(phen_all_use)[1]){
    phen_all_k[k,kl]<-ifelse(is.na(phen_all_use[k,kl])==F,1,0)
  }
}

##Variance explained -- sp~ph1+ph2+..

ind_l<-0
for(i in 1:dim(phen_all_k)[2]){
  if(length(which(is.na(phen_all_k[,i])==F))>0){
    ind_l<-rbind(ind_l,i)
  }
}
phen_all_k2<-phen_all_k[,ind_l]
colnames(phen_all_k2)
colnames(phen_all_k2)[3]<-"Cough"
colnames(phen_all_k2)[4:14]<-c("Diabetes (24)","Allergies (104)","Cystic fibrosis (18)","Heart disease (47)","Endocrine illness (16)","Renal failure (9)","Asthma (47)","Stomach illness (175)","Intestinal illness (70)","Arthritis (39)","MS (2)")
colnames(phen_all_k2)[15:18]<-c("HbA1c","Systolic","Diastolic","Hbtot")
colnames(phen_all_k2)[20:23]<-c("Heart rate","Perf.index","Pulse","O2sat")
colnames(phen_all_k2)[24:42]<-c("Beans","Tortillas","Rice","Bread","Milk","Yogurt","Cream/butter","Cheese","Eggs","Vegetables","Fruits","Natural juice","Chicken","Beef/Pork","Ham/sausages/hotdog","Fish","Soda","Fruit juice","Chips")
colnames(phen_all_k2)[43:49]<-c("Painkillers (714)","Antibiotics (137)","Anti-diarrheal (49)","Anti-parasite (34)","Anti-fungal (82)","Vitamins (121)","Anti-hypertensive (72)")
colnames(phen_all_k2)[50:59]<-c("Reserved","Trusting","Lazy","Relaxed","Not creative","Outgoing","Fault others","Thorough job","Nervous","Openess")
colnames(phen_all_k2)[61:86]<-c("Cat(645)","Dog(944)","Parakeet(102)","Rabbit(120)","Horse(314)","Mice(585)","None pet(88)","Cow(397)","Goat(40)","Pig(208)","Chicken(1094)","Duck(478)","Turkey(322)","Sheep(33)","Geese(68)","None farm(69)","Bat(410)","Lizard(423)","Monkey(12)","Snake(437)","Bird(794)","Possum(467)","Rat(457)","Squirrel(480)","Other Wild(1)","None wild(243)")
colnames(phen_all_k2)[87:88]<-c("Travel","Monthly expenditure")
colnames(phen_all_k2)[95:98]<-c("Living with partner","Partner live duration","Partner #","Partner live age")
op<-colnames(phen_all_k2)
#write.csv(op,"op.csv")
colnames(phen_all_k2)

data_use<-as.data.frame(t(data_transformed))
library(vegan)
#pl_perm<-adonis2(data_transformed~`gender`+`age_at_survey`+`Cough`+`Diabetes (24)`+`Allergies (104)`+`Cystic fibrosis (18)`+`Heart disease (47)`+`Endocrine illness (16)`+`Renal failure (9)`+`Asthma (47)`+`Stomach illness (175)`+`Intestinal illness (70)`+`Arthritis (39)`+`HbA1c`+`Systolic`+`Diastolic`+`Hbtot`+`BMI`+`Heart rate`+`Perf.index`+`Pulse`+`O2sat`+`Beans`+`Tortillas`+`Rice`+`Bread`+`Milk`+`Yogurt`+`Cream/butter`+`Cheese`+`Eggs`+`Vegetables`+`Fruits`+`Natural juice`+`Chicken`+`Beef/Pork`+`Ham/sausages/hotdog`+`Fish`+`Soda`+`Fruit juice`+`Chips`+`Painkillers (714)`+`Antibiotics (137)`+`Anti-diarrheal (49)`+`Anti-parasite (34)`+`Anti-fungal (82)`+`Vitamins (121)`+`Anti-hypertensive (72)`+`Reserved`+`Trusting`+`Lazy`+`Relaxed`+`Not creative`+`Outgoing`+`Fault others`+`Thorough job`+`Nervous`+`Openess`+`cognitive_score`+`Cat(645)`+`Dog(944)`+`Parakeet(102)`+`Rabbit(120)`+`Horse(314)`+`Mice(585)`+`None pet(88)`+`Cow(397)`+`Goat(40)`+`Pig(208)`+`Chicken(1094)`+`Duck(478)`+`Turkey(322)`+`Sheep(33)`+`Geese(68)`+`None farm(69)`+`Bat(410)`+`Lizard(423)`+`Monkey(12)`+`Snake(437)`+`Bird(794)`+`Possum(467)`+`Rat(457)`+`Squirrel(480)`+`Other Wild(1)`+`None wild(243)`+`Travel`+`Monthly expenditure`+`Friend.ties.same.building.`+`Friend.ties.different.building.`+`Betweeness.friendship.`+`Transitivity.friendship.`+`risky_int`+`altruism`+`Living with partner`+`Partner live duration`+`Partner #`+`Partner live age`+`food.cook.chimney`+`food.cook.no.chimney`+`stove`+`none.stove`+`Electricity`+`Radio`+`Television`+`Cell.mobile.phone`+`Non.mobile.phone`+`Refrigerator`+`None`+`Wood`+`Gas`+`fuel.Electricity`+`Kerosene`+`None.fuel`+`Separate.kitchen`+`Cement.floor`+`Earth.Sand.floor`+`Ceramic.floor`+`Tiles.floor`+`Mud.bricks.floor`+`Wood.floor`+`Other.floor`+`Wooden.windows`+`Glass.windows`+`Metal.windows`+`Unfinished.windows`+`No.windows`+`Clay.mud.walls`+`Clay.brick.walls`+`Cement.walls`+ `Cane.palm.trunks.walls`+`Wood.unpolished.walls`+`Wood.polished.walls`+`Discarded.materials.walls`+`No.walls`+`Other.walls`+`Plastic.roof`+`Metal.roof`+`Clay.roof`+`Thatch.palm.roof`+`Concrete.roof`+`Wood.roof`+`Other.roof`+`Sleeping.rooms`+`Spring.water.protected`+`Spring.water.unprotected`+`tube.well`+`Dug.well.protected`+`Dug.well.unprotected`+`Surface.water`+`bottle.water`+ `Other.water`+`Familial.ties.same.building.`+`Familial.ties.different.building.`+`Transitivity.familial.`+`Betweeness.familial.`+`risk_taking`+`Degree`+`Percent.kin`+`Clustering.coefficient`+`Betweeness`+`kcycle`+`distance_center`+`hh_size`+`household_wealth_index`+`time_to_main_road`+`time_to_health_ctr`+`prop_defor_500`+`Bristol`+`MAP`+`Alcohol(daily)`+`Cigarette use`+`Cigarette frequency`+`Education(primary)`+`Education(middle)``Education(secondary)`+`Washing hands`+`GAD7(mild)`+`GAD7(moderate)`+`GAD7(severe)`+`PHQ9(mild)`+`PHQ9(moderate)`+`PHQ9(severe)`, data=phen_all_k2, permutations = 999)

#pl_perm<-adonis2(data_transformed~`gender`+`age_at_survey`+`Cough`+`Diabetes (24)`+`Allergies (104)`+`Cystic fibrosis (18)`+`Heart disease (47)`+`Endocrine illness (16)`+`Renal failure (9)`+`Asthma (47)`+`Stomach illness (175)`+`Intestinal illness (70)`+`Arthritis (39)`+`HbA1c`+`Systolic`+`Diastolic`+`Hbtot`+`BMI`+`Heart rate`+`Bristol`+`MAP`+`Alcohol(daily)`+`Cigarette use`+`Cigarette frequency`+`Education(primary)`+`Education(middle)``Education(secondary)`+`Washing hands`+`GAD7(mild)`+`GAD7(moderate)`+`GAD7(severe)`+`PHQ9(mild)`+`PHQ9(moderate)`+`PHQ9(severe)`, data=phen_all_k2, permutations = 999)
#pl_perm<-adonis2(data_use~`gender`+`age_at_survey`+`Cough`, data=as.data.frame(phen_all_k2), permutations = 999,na.action=na.omit)

#pl_perm2<-adonis2(data_use~`gender`+`age_at_survey`+`Cough`+`Diabetes (24)`+`Allergies (104)`+`Cystic fibrosis (18)`+`Heart disease (47)`+`Endocrine illness (16)`+`Renal failure (9)`+`Asthma (47)`+`Stomach illness (175)`+`Intestinal illness (70)`+`Arthritis (39)`+`HbA1c`+`Systolic`+`Diastolic`+`Hbtot`+`BMI`+`Heart rate`+`Bristol`+`MAP`+`Alcohol(daily)`+`Cigarette use`+`Cigarette frequency`+`GAD7(mild)`+`GAD7(moderate)`+`GAD7(severe)`+`PHQ9(mild)`+`PHQ9(moderate)`+`PHQ9(severe)`, data=as.data.frame(phen_all_k2), permutations = 999,na.action=na.omit)

#length(which(is.na(phen_all_k2[,3])==T))

#lk<-pl_perm2

phen_all_k2[,174]<-ifelse(is.na(phen_all_k2[,174])==T,0,phen_all_k2[,174])
phen_all_k2[,c(24:42)]<-ifelse(is.na(phen_all_k2[,c(24:42)])==T,0,phen_all_k2[,c(24:42)])
phen_all_k2[,c(61:87)]<-ifelse(is.na(phen_all_k2[,c(61:87)])==T,0,phen_all_k2[,c(61:87)])
phen_all_k2[,c(95:98)]<-ifelse(is.na(phen_all_k2[,c(95:98)])==T,0,phen_all_k2[,c(95:98)])
phen_all_k2[,c(162,175:177,158:161,99:152)]<-ifelse(is.na(phen_all_k2[,c(162,175:177,158:161,99:152)])==T,0,phen_all_k2[,c(162,175:177,158:161,99:152)])

phen_all_k34<-na.omit(phen_all_k2)
dim(phen_all_k34)

colnames(phen_all_k2)

phen_all_use22<-read.csv('phen_all_use6_all.csv',row.names = 1,header=T)#phen_all_use4_all.csv

phen_all_k2<-cbind(phen_all_k2,"Diarrhea"=array(NA,dim=c(dim(phen_all_k2)[1],1)))

for(kl in 381){
  for(k in 1:dim(phen_all_use22)[1]){
    phen_all_k2[k,dim(phen_all_k2)[2]]<-ifelse(phen_all_use22[k,kl]=="No",0,1)
  }
}
colnames(phen_all_k2)[dim(phen_all_k2)[2]]<-"Diarrhea(73)"

phen_all_k2<-cbind(phen_all_k2,"batch_effect"=array(NA,dim=c(dim(phen_all_k2)[1],1)))
for(kl in 224){
  for(k in 1:dim(phen_all_use)[1]){
    phen_all_k2[k,dim(phen_all_k2)[2]]<-as.numeric(as.character(phen_all_use[k,kl]))
  }
}
colnames(phen_all_k2)[(dim(phen_all_k2)[2])]<-c("batch_effect")

phen_all_k2<-cbind(phen_all_k2,"dna_conc"=array(NA,dim=c(dim(phen_all_k2)[1],1)))
for(kl in 225){
  for(k in 1:dim(phen_all_use)[1]){
    phen_all_k2[k,dim(phen_all_k2)[2]]<-as.numeric(as.character(phen_all_use[k,kl]))
  }
}
colnames(phen_all_k2)[(dim(phen_all_k2)[2])]<-"dna_conc"

phen_all_k34<-na.omit(phen_all_k2)
dim(phen_all_k34)


##na.omit -- omitting too many values
phen_all_k3<-phen_all_k2
#phen_all_k3[,3:14]<-ifelse(is.na(phen_all_k3[,3:14])==T,0,phen_all_k3[,3:14])
#phen_all_k3[,24:87]<-ifelse(is.na(phen_all_k3[,24:87])==T,0,phen_all_k3[,24:87])
#phen_all_k3[,174]<-ifelse(is.na(phen_all_k3[,174])==T,0,phen_all_k3[,174])

#Include diarrhea after data request
pl_perm2_health<-adonis2(data_use~`gender`+`age_at_survey`+`Cough`+`Diabetes (24)`+`Allergies (104)`+`Cystic fibrosis (18)`+`Heart disease (47)`+`Endocrine illness (16)`+`Renal failure (9)`+`Asthma (47)`+`Stomach illness (175)`+`Intestinal illness (70)`+`Arthritis (39)`+`HbA1c`+`Systolic`+`Diastolic`+`Hbtot`+`BMI`+`Heart rate`+`Painkillers (714)`+`Antibiotics (137)`+`Anti-diarrheal (49)`+`Anti-parasite (34)`+`Anti-fungal (82)`+`Vitamins (121)`+`Anti-hypertensive (72)`+`Reserved`+`Trusting`+`Lazy`+`Relaxed`+`Not creative`+`Outgoing`+`Fault others`+`Thorough job`+`Nervous`+`Openess`+`cognitive_score`+`Bristol`+`MAP`+`Alcohol(daily)`+`Cigarette use`+`Cigarette frequency`+`GAD7(mild)`+`GAD7(moderate)`+`GAD7(severe)`+`PHQ9(mild)`+`PHQ9(moderate)`+`PHQ9(severe)`, data=as.data.frame(phen_all_k3), permutations = 999,na.action=na.omit)

#Food and animals
pl_perm2_animals<-adonis2(data_use~`Beans`+`Tortillas`+`Rice`+`Bread`+`Milk`+`Yogurt`+`Cream/butter`+`Cheese`+`Eggs`+`Vegetables`+`Fruits`+`Natural juice`+`Chicken`+`Beef/Pork`+`Ham/sausages/hotdog`+`Fish`+`Soda`+`Fruit juice`+`Chips`+`Cat(645)`+`Dog(944)`+`Parakeet(102)`+`Rabbit(120)`+`Horse(314)`+`Mice(585)`+`None pet(88)`+`Cow(397)`+`Goat(40)`+`Pig(208)`+`Chicken(1094)`+`Duck(478)`+`Turkey(322)`+`Sheep(33)`+`Geese(68)`+`None farm(69)`+`Bat(410)`+`Lizard(423)`+`Monkey(12)`+`Snake(437)`+`Bird(794)`+`Possum(467)`+`Rat(457)`+`Squirrel(480)`+`Other Wild(1)`+`None wild(243)`, data=as.data.frame(phen_all_k3), permutations = 999,na.action=na.omit)

write.csv(pl_perm2_health,'pl_perm2_health.csv')
write.csv(pl_perm2_animals,'pl_perm2_animals.csv')

#Socio-economic
pl_perm2_socio_eco<-adonis2(data_use~`Travel`+`Monthly expenditure`+`Friend.ties.same.building.`+`Friend.ties.different.building.`+`Betweeness.friendship.`+`Transitivity.friendship.`+`risky_int`+`altruism`+`Living with partner`+`Partner live duration`+`Partner #`+`food.cook.chimney`+`food.cook.no.chimney`+`stove`+`none.stove`+`Electricity`+`Radio`+`Television`+`Cell.mobile.phone`+`Non.mobile.phone`+`Refrigerator`+`None`+`Wood`+`Gas`+`fuel.Electricity`+`Kerosene`+`None.fuel`+`Separate.kitchen`+`Cement.floor`+`Earth.Sand.floor`+`Ceramic.floor`+`Tiles.floor`+`Mud.bricks.floor`+`Wood.floor`+`Other.floor`+`Wooden.windows`+`Glass.windows`+`Metal.windows`+`Unfinished.windows`+`No.windows`+`Clay.mud.walls`+`Clay.brick.walls`+`Cement.walls`+ `Cane.palm.trunks.walls`+`Wood.unpolished.walls`+`Wood.polished.walls`+`Discarded.materials.walls`+`No.walls`+`Other.walls`+`Plastic.roof`+`Metal.roof`+`Clay.roof`+`Thatch.palm.roof`+`Concrete.roof`+`Wood.roof`+`Other.roof`+`Sleeping.rooms`+`Spring.water.protected`+`Spring.water.unprotected`+`tube.well`+`Dug.well.protected`+`Dug.well.unprotected`+`Surface.water`+`bottle.water`+ `Other.water`+`Familial.ties.same.building.`+`Familial.ties.different.building.`+`Transitivity.familial.`+`Betweeness.familial.`+`risk_taking`+`Degree`+`Percent.kin`+`Clustering.coefficient`+`Betweeness`+`kcycle`+`distance_center`+`hh_size`+`household_wealth_index`+`time_to_main_road`+`time_to_health_ctr`+`prop_defor_500`+`Education(primary)`+`Education(middle)`+`Education(secondary)`+`Washing hands`, data=as.data.frame(phen_all_k3), permutations = 999,na.action=na.omit)

write.csv(pl_perm2_socio_eco,'pl_perm2_socio_eco.csv')
#All

pl_perm2_all<-adonis2(data_use~`gender`+`age_at_survey`+`Cough`+`Diabetes (24)`+`Allergies (104)`+`Cystic fibrosis (18)`+`Heart disease (47)`+`Endocrine illness (16)`+`Renal failure (9)`+`Asthma (47)`+`Stomach illness (175)`+`Intestinal illness (70)`+`Arthritis (39)`+`HbA1c`+`Systolic`+`Diastolic`+`Hbtot`+`BMI`+`Heart rate`+`Perf.index`+`Pulse`+`O2sat`+`Beans`+`Tortillas`+`Rice`+`Bread`+`Milk`+`Yogurt`+`Cream/butter`+`Cheese`+`Eggs`+`Vegetables`+`Fruits`+`Natural juice`+`Chicken`+`Beef/Pork`+`Ham/sausages/hotdog`+`Fish`+`Soda`+`Fruit juice`+`Chips`+`Painkillers (714)`+`Antibiotics (137)`+`Anti-diarrheal (49)`+`Anti-parasite (34)`+`Anti-fungal (82)`+`Vitamins (121)`+`Anti-hypertensive (72)`+`Reserved`+`Trusting`+`Lazy`+`Relaxed`+`Not creative`+`Outgoing`+`Fault others`+`Thorough job`+`Nervous`+`Openess`+`cognitive_score`+`Cat(645)`+`Dog(944)`+`Parakeet(102)`+`Rabbit(120)`+`Horse(314)`+`Mice(585)`+`None pet(88)`+`Cow(397)`+`Goat(40)`+`Pig(208)`+`Chicken(1094)`+`Duck(478)`+`Turkey(322)`+`Sheep(33)`+`Geese(68)`+`None farm(69)`+`Bat(410)`+`Lizard(423)`+`Monkey(12)`+`Snake(437)`+`Bird(794)`+`Possum(467)`+`Rat(457)`+`Squirrel(480)`+`Other Wild(1)`+`None wild(243)`+`Travel`+`Monthly expenditure`+`Friend.ties.same.building.`+`Friend.ties.different.building.`+`Betweeness.friendship.`+`Transitivity.friendship.`+`risky_int`+`altruism`+`Living with partner`+`Partner live duration`+`Partner #`+`food.cook.chimney`+`food.cook.no.chimney`+`stove`+`none.stove`+`Electricity`+`Radio`+`Television`+`Cell.mobile.phone`+`Non.mobile.phone`+`Refrigerator`+`None`+`Wood`+`Gas`+`fuel.Electricity`+`Kerosene`+`None.fuel`+`Separate.kitchen`+`Cement.floor`+`Earth.Sand.floor`+`Ceramic.floor`+`Tiles.floor`+`Mud.bricks.floor`+`Wood.floor`+`Other.floor`+`Wooden.windows`+`Glass.windows`+`Metal.windows`+`Unfinished.windows`+`No.windows`+`Clay.mud.walls`+`Clay.brick.walls`+`Cement.walls`+ `Cane.palm.trunks.walls`+`Wood.unpolished.walls`+`Wood.polished.walls`+`Discarded.materials.walls`+`No.walls`+`Other.walls`+`Plastic.roof`+`Metal.roof`+`Clay.roof`+`Thatch.palm.roof`+`Concrete.roof`+`Wood.roof`+`Other.roof`+`Sleeping.rooms`+`Spring.water.protected`+`Spring.water.unprotected`+`tube.well`+`Dug.well.protected`+`Dug.well.unprotected`+`Surface.water`+`bottle.water`+ `Other.water`+`Familial.ties.same.building.`+`Familial.ties.different.building.`+`Transitivity.familial.`+`Betweeness.familial.`+`risk_taking`+`Degree`+`Percent.kin`+`Clustering.coefficient`+`Betweeness`+`kcycle`+`distance_center`+`hh_size`+`household_wealth_index`+`time_to_main_road`+`time_to_health_ctr`+`prop_defor_500`+`Bristol`+`MAP`+`Alcohol(daily)`+`Cigarette use`+`Cigarette frequency`+`Education(primary)`+`Education(middle)`+`Education(secondary)`+`Washing hands`+`GAD7(mild)`+`GAD7(moderate)`+`GAD7(severe)`+`PHQ9(mild)`+`PHQ9(moderate)`+`PHQ9(severe)`+`batch_effect`+`dna_conc`, data=as.data.frame(phen_all_k3), permutations = 999,na.action=na.omit)

write.csv(pl_perm2_all,'pl_perm2_all.csv')


######
##Pathways

data_transformed<-read.csv('pwy_all_use_un_norm.csv',row.names = 1,header=T)
data_transformed<-as.matrix.data.frame(data_transformed)
data_use_pwy<-as.data.frame(t(data_transformed))

#Include diarrhea after data request
pl_perm2_health_pwy<-adonis2(data_use_pwy~`gender`+`age_at_survey`+`Cough`+`Diabetes (24)`+`Allergies (104)`+`Cystic fibrosis (18)`+`Heart disease (47)`+`Endocrine illness (16)`+`Renal failure (9)`+`Asthma (47)`+`Stomach illness (175)`+`Intestinal illness (70)`+`Arthritis (39)`+`HbA1c`+`Systolic`+`Diastolic`+`Hbtot`+`BMI`+`Heart rate`+`Painkillers (714)`+`Antibiotics (137)`+`Anti-diarrheal (49)`+`Anti-parasite (34)`+`Anti-fungal (82)`+`Vitamins (121)`+`Anti-hypertensive (72)`+`Reserved`+`Trusting`+`Lazy`+`Relaxed`+`Not creative`+`Outgoing`+`Fault others`+`Thorough job`+`Nervous`+`Openess`+`cognitive_score`+`Bristol`+`MAP`+`Alcohol(daily)`+`Cigarette use`+`Cigarette frequency`+`GAD7(mild)`+`GAD7(moderate)`+`GAD7(severe)`+`PHQ9(mild)`+`PHQ9(moderate)`+`PHQ9(severe)`, data=as.data.frame(phen_all_k3), permutations = 999,na.action=na.omit)

#Food and animals
pl_perm2_animals_pwy<-adonis2(data_use_pwy~`Beans`+`Tortillas`+`Rice`+`Bread`+`Milk`+`Yogurt`+`Cream/butter`+`Cheese`+`Eggs`+`Vegetables`+`Fruits`+`Natural juice`+`Chicken`+`Beef/Pork`+`Ham/sausages/hotdog`+`Fish`+`Soda`+`Fruit juice`+`Chips`+`Cat(645)`+`Dog(944)`+`Parakeet(102)`+`Rabbit(120)`+`Horse(314)`+`Mice(585)`+`None pet(88)`+`Cow(397)`+`Goat(40)`+`Pig(208)`+`Chicken(1094)`+`Duck(478)`+`Turkey(322)`+`Sheep(33)`+`Geese(68)`+`None farm(69)`+`Bat(410)`+`Lizard(423)`+`Monkey(12)`+`Snake(437)`+`Bird(794)`+`Possum(467)`+`Rat(457)`+`Squirrel(480)`+`Other Wild(1)`+`None wild(243)`, data=as.data.frame(phen_all_k3), permutations = 999,na.action=na.omit)

write.csv(pl_perm2_health_pwy,'pl_perm2_health_pwy.csv')
write.csv(pl_perm2_animals_pwy,'pl_perm2_animals_pwy.csv')

#Socio-economic
pl_perm2_socio_eco_pwy<-adonis2(data_use_pwy~`Travel`+`Monthly expenditure`+`Friend.ties.same.building.`+`Friend.ties.different.building.`+`Betweeness.friendship.`+`Transitivity.friendship.`+`risky_int`+`altruism`+`Living with partner`+`Partner live duration`+`Partner #`+`food.cook.chimney`+`food.cook.no.chimney`+`stove`+`none.stove`+`Electricity`+`Radio`+`Television`+`Cell.mobile.phone`+`Non.mobile.phone`+`Refrigerator`+`None`+`Wood`+`Gas`+`fuel.Electricity`+`Kerosene`+`None.fuel`+`Separate.kitchen`+`Cement.floor`+`Earth.Sand.floor`+`Ceramic.floor`+`Tiles.floor`+`Mud.bricks.floor`+`Wood.floor`+`Other.floor`+`Wooden.windows`+`Glass.windows`+`Metal.windows`+`Unfinished.windows`+`No.windows`+`Clay.mud.walls`+`Clay.brick.walls`+`Cement.walls`+ `Cane.palm.trunks.walls`+`Wood.unpolished.walls`+`Wood.polished.walls`+`Discarded.materials.walls`+`No.walls`+`Other.walls`+`Plastic.roof`+`Metal.roof`+`Clay.roof`+`Thatch.palm.roof`+`Concrete.roof`+`Wood.roof`+`Other.roof`+`Sleeping.rooms`+`Spring.water.protected`+`Spring.water.unprotected`+`tube.well`+`Dug.well.protected`+`Dug.well.unprotected`+`Surface.water`+`bottle.water`+ `Other.water`+`Familial.ties.same.building.`+`Familial.ties.different.building.`+`Transitivity.familial.`+`Betweeness.familial.`+`risk_taking`+`Degree`+`Percent.kin`+`Clustering.coefficient`+`Betweeness`+`kcycle`+`distance_center`+`hh_size`+`household_wealth_index`+`time_to_main_road`+`time_to_health_ctr`+`prop_defor_500`+`Education(primary)`+`Education(middle)`+`Education(secondary)`+`Washing hands`, data=as.data.frame(phen_all_k3), permutations = 999,na.action=na.omit)

write.csv(pl_perm2_socio_eco_pwy,'pl_perm2_socio_eco_pwy.csv')
#All

pl_perm2_all_pwy<-adonis2(data_use_pwy~`gender`+`age_at_survey`+`Cough`+`Diabetes (24)`+`Allergies (104)`+`Cystic fibrosis (18)`+`Heart disease (47)`+`Endocrine illness (16)`+`Renal failure (9)`+`Asthma (47)`+`Stomach illness (175)`+`Intestinal illness (70)`+`Arthritis (39)`+`HbA1c`+`Systolic`+`Diastolic`+`Hbtot`+`BMI`+`Heart rate`+`Perf.index`+`Pulse`+`O2sat`+`Beans`+`Tortillas`+`Rice`+`Bread`+`Milk`+`Yogurt`+`Cream/butter`+`Cheese`+`Eggs`+`Vegetables`+`Fruits`+`Natural juice`+`Chicken`+`Beef/Pork`+`Ham/sausages/hotdog`+`Fish`+`Soda`+`Fruit juice`+`Chips`+`Painkillers (714)`+`Antibiotics (137)`+`Anti-diarrheal (49)`+`Anti-parasite (34)`+`Anti-fungal (82)`+`Vitamins (121)`+`Anti-hypertensive (72)`+`Reserved`+`Trusting`+`Lazy`+`Relaxed`+`Not creative`+`Outgoing`+`Fault others`+`Thorough job`+`Nervous`+`Openess`+`cognitive_score`+`Cat(645)`+`Dog(944)`+`Parakeet(102)`+`Rabbit(120)`+`Horse(314)`+`Mice(585)`+`None pet(88)`+`Cow(397)`+`Goat(40)`+`Pig(208)`+`Chicken(1094)`+`Duck(478)`+`Turkey(322)`+`Sheep(33)`+`Geese(68)`+`None farm(69)`+`Bat(410)`+`Lizard(423)`+`Monkey(12)`+`Snake(437)`+`Bird(794)`+`Possum(467)`+`Rat(457)`+`Squirrel(480)`+`Other Wild(1)`+`None wild(243)`+`Travel`+`Monthly expenditure`+`Friend.ties.same.building.`+`Friend.ties.different.building.`+`Betweeness.friendship.`+`Transitivity.friendship.`+`risky_int`+`altruism`+`Living with partner`+`Partner live duration`+`Partner #`+`food.cook.chimney`+`food.cook.no.chimney`+`stove`+`none.stove`+`Electricity`+`Radio`+`Television`+`Cell.mobile.phone`+`Non.mobile.phone`+`Refrigerator`+`None`+`Wood`+`Gas`+`fuel.Electricity`+`Kerosene`+`None.fuel`+`Separate.kitchen`+`Cement.floor`+`Earth.Sand.floor`+`Ceramic.floor`+`Tiles.floor`+`Mud.bricks.floor`+`Wood.floor`+`Other.floor`+`Wooden.windows`+`Glass.windows`+`Metal.windows`+`Unfinished.windows`+`No.windows`+`Clay.mud.walls`+`Clay.brick.walls`+`Cement.walls`+ `Cane.palm.trunks.walls`+`Wood.unpolished.walls`+`Wood.polished.walls`+`Discarded.materials.walls`+`No.walls`+`Other.walls`+`Plastic.roof`+`Metal.roof`+`Clay.roof`+`Thatch.palm.roof`+`Concrete.roof`+`Wood.roof`+`Other.roof`+`Sleeping.rooms`+`Spring.water.protected`+`Spring.water.unprotected`+`tube.well`+`Dug.well.protected`+`Dug.well.unprotected`+`Surface.water`+`bottle.water`+ `Other.water`+`Familial.ties.same.building.`+`Familial.ties.different.building.`+`Transitivity.familial.`+`Betweeness.familial.`+`risk_taking`+`Degree`+`Percent.kin`+`Clustering.coefficient`+`Betweeness`+`kcycle`+`distance_center`+`hh_size`+`household_wealth_index`+`time_to_main_road`+`time_to_health_ctr`+`prop_defor_500`+`Bristol`+`MAP`+`Alcohol(daily)`+`Cigarette use`+`Cigarette frequency`+`Education(primary)`+`Education(middle)`+`Education(secondary)`+`Washing hands`+`GAD7(mild)`+`GAD7(moderate)`+`GAD7(severe)`+`PHQ9(mild)`+`PHQ9(moderate)`+`PHQ9(severe)`+`batch_effect`+`dna_conc`, data=as.data.frame(phen_all_k3), permutations = 999,na.action=na.omit)

write.csv(pl_perm2_all_pwy,'pl_perm2_all_pwy.csv')

save.image(file="adonis_all_vars_sp_pwy.RData")

##na.omit -- omitting too many values
#phen_all_k3<-phen_all_k2
#phen_all_k3[,3:14]<-ifelse(is.na(phen_all_k3[,3:14])==T,0,phen_all_k3[,3:14])
#phen_all_k3[,24:87]<-ifelse(is.na(phen_all_k3[,24:87])==T,0,phen_all_k3[,24:87])
#phen_all_k3[,174]<-ifelse(is.na(phen_all_k3[,174])==T,0,phen_all_k3[,174])

#Include diarrhea after data request
#pl_perm2_health<-adonis2(data_use~`gender`+`age_at_survey`+`Cough`+`Diabetes (24)`+`Allergies (104)`+`Cystic fibrosis (18)`+`Heart disease (47)`+`Endocrine illness (16)`+`Renal failure (9)`+`Asthma (47)`+`Stomach illness (175)`+`Intestinal illness (70)`+`Arthritis (39)`+`HbA1c`+`Systolic`+`Diastolic`+`Hbtot`+`BMI`+`Heart rate`+`Painkillers (714)`+`Antibiotics (137)`+`Anti-diarrheal (49)`+`Anti-parasite (34)`+`Anti-fungal (82)`+`Vitamins (121)`+`Anti-hypertensive (72)`+`Reserved`+`Trusting`+`Lazy`+`Relaxed`+`Not creative`+`Outgoing`+`Fault others`+`Thorough job`+`Nervous`+`Openess`+`cognitive_score`+`Bristol`+`MAP`+`Alcohol(daily)`+`Cigarette use`+`Cigarette frequency`+`GAD7(mild)`+`GAD7(moderate)`+`GAD7(severe)`+`PHQ9(mild)`+`PHQ9(moderate)`+`PHQ9(severe)`, data=as.data.frame(phen_all_k3), permutations = 999,na.action=na.omit)

#Food and animals
#pl_perm2_animals<-adonis2(data_use~`Beans`+`Tortillas`+`Rice`+`Bread`+`Milk`+`Yogurt`+`Cream/butter`+`Cheese`+`Eggs`+`Vegetables`+`Fruits`+`Natural juice`+`Chicken`+`Beef/Pork`+`Ham/sausages/hotdog`+`Fish`+`Soda`+`Fruit juice`+`Chips`+`Cat(645)`+`Dog(944)`+`Parakeet(102)`+`Rabbit(120)`+`Horse(314)`+`Mice(585)`+`None pet(88)`+`Cow(397)`+`Goat(40)`+`Pig(208)`+`Chicken(1094)`+`Duck(478)`+`Turkey(322)`+`Sheep(33)`+`Geese(68)`+`None farm(69)`+`Bat(410)`+`Lizard(423)`+`Monkey(12)`+`Snake(437)`+`Bird(794)`+`Possum(467)`+`Rat(457)`+`Squirrel(480)`+`Other Wild(1)`+`None wild(243)`, data=as.data.frame(phen_all_k3), permutations = 999,na.action=na.omit)

#Socio-economic
#pl_perm2_socio_eco<-adonis2(data_use~`Travel`+`Monthly expenditure`+`Friend.ties.same.building.`+`Friend.ties.different.building.`+`Betweeness.friendship.`+`Transitivity.friendship.`+`risky_int`+`altruism`+`Living with partner`+`Partner live duration`+`Partner #`+`food.cook.chimney`+`food.cook.no.chimney`+`stove`+`none.stove`+`Electricity`+`Radio`+`Television`+`Cell.mobile.phone`+`Non.mobile.phone`+`Refrigerator`+`None`+`Wood`+`Gas`+`fuel.Electricity`+`Kerosene`+`None.fuel`+`Separate.kitchen`+`Cement.floor`+`Earth.Sand.floor`+`Ceramic.floor`+`Tiles.floor`+`Mud.bricks.floor`+`Wood.floor`+`Other.floor`+`Wooden.windows`+`Glass.windows`+`Metal.windows`+`Unfinished.windows`+`No.windows`+`Clay.mud.walls`+`Clay.brick.walls`+`Cement.walls`+ `Cane.palm.trunks.walls`+`Wood.unpolished.walls`+`Wood.polished.walls`+`Discarded.materials.walls`+`No.walls`+`Other.walls`+`Plastic.roof`+`Metal.roof`+`Clay.roof`+`Thatch.palm.roof`+`Concrete.roof`+`Wood.roof`+`Other.roof`+`Sleeping.rooms`+`Spring.water.protected`+`Spring.water.unprotected`+`tube.well`+`Dug.well.protected`+`Dug.well.unprotected`+`Surface.water`+`bottle.water`+ `Other.water`+`Familial.ties.same.building.`+`Familial.ties.different.building.`+`Transitivity.familial.`+`Betweeness.familial.`+`risk_taking`+`Degree`+`Percent.kin`+`Clustering.coefficient`+`Betweeness`+`kcycle`+`distance_center`+`hh_size`+`household_wealth_index`+`time_to_main_road`+`time_to_health_ctr`+`prop_defor_500`+`Education(primary)`+`Education(middle)``Education(secondary)`+`Washing hands`, data=as.data.frame(phen_all_k3), permutations = 999,na.action=na.omit)


#All

#pl_perm2_all<-adonis2(data_use~`gender`+`age_at_survey`+`Cough`+`Diabetes (24)`+`Allergies (104)`+`Cystic fibrosis (18)`+`Heart disease (47)`+`Endocrine illness (16)`+`Renal failure (9)`+`Asthma (47)`+`Stomach illness (175)`+`Intestinal illness (70)`+`Arthritis (39)`+`HbA1c`+`Systolic`+`Diastolic`+`Hbtot`+`BMI`+`Heart rate`+`Perf.index`+`Pulse`+`O2sat`+`Beans`+`Tortillas`+`Rice`+`Bread`+`Milk`+`Yogurt`+`Cream/butter`+`Cheese`+`Eggs`+`Vegetables`+`Fruits`+`Natural juice`+`Chicken`+`Beef/Pork`+`Ham/sausages/hotdog`+`Fish`+`Soda`+`Fruit juice`+`Chips`+`Painkillers (714)`+`Antibiotics (137)`+`Anti-diarrheal (49)`+`Anti-parasite (34)`+`Anti-fungal (82)`+`Vitamins (121)`+`Anti-hypertensive (72)`+`Reserved`+`Trusting`+`Lazy`+`Relaxed`+`Not creative`+`Outgoing`+`Fault others`+`Thorough job`+`Nervous`+`Openess`+`cognitive_score`+`Cat(645)`+`Dog(944)`+`Parakeet(102)`+`Rabbit(120)`+`Horse(314)`+`Mice(585)`+`None pet(88)`+`Cow(397)`+`Goat(40)`+`Pig(208)`+`Chicken(1094)`+`Duck(478)`+`Turkey(322)`+`Sheep(33)`+`Geese(68)`+`None farm(69)`+`Bat(410)`+`Lizard(423)`+`Monkey(12)`+`Snake(437)`+`Bird(794)`+`Possum(467)`+`Rat(457)`+`Squirrel(480)`+`Other Wild(1)`+`None wild(243)`+`Travel`+`Monthly expenditure`+`Friend.ties.same.building.`+`Friend.ties.different.building.`+`Betweeness.friendship.`+`Transitivity.friendship.`+`risky_int`+`altruism`+`Living with partner`+`Partner live duration`+`Partner #`+`food.cook.chimney`+`food.cook.no.chimney`+`stove`+`none.stove`+`Electricity`+`Radio`+`Television`+`Cell.mobile.phone`+`Non.mobile.phone`+`Refrigerator`+`None`+`Wood`+`Gas`+`fuel.Electricity`+`Kerosene`+`None.fuel`+`Separate.kitchen`+`Cement.floor`+`Earth.Sand.floor`+`Ceramic.floor`+`Tiles.floor`+`Mud.bricks.floor`+`Wood.floor`+`Other.floor`+`Wooden.windows`+`Glass.windows`+`Metal.windows`+`Unfinished.windows`+`No.windows`+`Clay.mud.walls`+`Clay.brick.walls`+`Cement.walls`+ `Cane.palm.trunks.walls`+`Wood.unpolished.walls`+`Wood.polished.walls`+`Discarded.materials.walls`+`No.walls`+`Other.walls`+`Plastic.roof`+`Metal.roof`+`Clay.roof`+`Thatch.palm.roof`+`Concrete.roof`+`Wood.roof`+`Other.roof`+`Sleeping.rooms`+`Spring.water.protected`+`Spring.water.unprotected`+`tube.well`+`Dug.well.protected`+`Dug.well.unprotected`+`Surface.water`+`bottle.water`+ `Other.water`+`Familial.ties.same.building.`+`Familial.ties.different.building.`+`Transitivity.familial.`+`Betweeness.familial.`+`risk_taking`+`Degree`+`Percent.kin`+`Clustering.coefficient`+`Betweeness`+`kcycle`+`distance_center`+`hh_size`+`household_wealth_index`+`time_to_main_road`+`time_to_health_ctr`+`prop_defor_500`+`Bristol`+`MAP`+`Alcohol(daily)`+`Cigarette use`+`Cigarette frequency`+`Education(primary)`+`Education(middle)``Education(secondary)`+`Washing hands`+`GAD7(mild)`+`GAD7(moderate)`+`GAD7(severe)`+`PHQ9(mild)`+`PHQ9(moderate)`+`PHQ9(severe)`, data=as.data.frame(phen_all_k3), permutations = 999,na.action=na.omit)

#phen_all_k2[,174]<-ifelse(is.na(phen_all_k2[,174])==T,0,phen_all_k2[,174])
#phen_all_k2[,c(24:42)]<-ifelse(is.na(phen_all_k2[,c(24:42)])==T,0,phen_all_k2[,c(24:42)])
#phen_all_k2[,c(61:87)]<-ifelse(is.na(phen_all_k2[,c(61:87)])==T,0,phen_all_k2[,c(61:87)])
#phen_all_k2[,c(95:98)]<-ifelse(is.na(phen_all_k2[,c(95:98)])==T,0,phen_all_k2[,c(95:98)])
#phen_all_k2[,c(162,175:177,158:161,99:152)]<-ifelse(is.na(phen_all_k2[,c(162,175:177,158:161,99:152)])==T,0,phen_all_k2[,c(162,175:177,158:161,99:152)])
#phen_all_k34<-na.omit(phen_all_k2)
#dim(phen_all_k34)

#colnames(phen_all_k2)





