
# ==================================
# By: Shivkumar Vishnempet Shridhar, Christakis and Brito group, HNL, Yale (2023)
# Honduras microbiome project, script of variance explained
# ==================================


library(vegan)

################################################

#Species

phen_all_use<-read.csv('phen_b4.csv',row.names=1)#sample x phenotype
phen_all_rct_use<-read.csv('phen_rct_b4.csv',row.names=1)#sample x phenotype

ind_list<-c(40,32,35,36,71:74,37)
colnames(phen_all_use)[ind_list]

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
phen_all_use<-cbind(phen_all_use,"MAP"=phen_map)#Creating mean arterial pressure (MAP) from systolic and diastolic

phys<-array(NA,dim=c(dim(phen_all_use)[1],6))
ind_list2<-c(32,384,40,71,74,37)
colnames(phen_all_use)[ind_list2]

for(i in 1:dim(phys)[1]){
  for(j in 1:length(ind_list2)){
    if(is.na(phen_all_use[i,ind_list2[j]])==F){
      phys[i,j]<-as.numeric(as.character(phen_all_use[i,ind_list2[j]]))
    }
  }
}

colnames(phys)<-c("Hb A1c","MAP","BMI","Heart rate","Oxygen saturation","Hb total")

#Overall health

for(i in 1:4){
  if(i==1){
    ph2<-as.character(phen_all_use[,15])
  }else{
    ph2<-cbind(ph2,as.character(phen_all_use[,15]))
  }
}

ph2[,1]<-ifelse(is.na(ph2[,1])==F,ifelse(ph2[,1]=="Poor",1,0),NA)
ph2[,2]<-ifelse(is.na(ph2[,2])==F,ifelse(ph2[,2]=="Fair",1,0),NA)
ph2[,3]<-ifelse(is.na(ph2[,3])==F,ifelse(ph2[,3]=="Very good",1,0),NA)
ph2[,4]<-ifelse(is.na(ph2[,4])==F,ifelse(ph2[,4]=="Excellent",1,0),NA)

colnames(ph2)<-c("Poor","Fair","Very good","Excellent")

#Acute conditions cough, diarrhea, bristol

#ind_list21<-c(16,381)
ph2<-cbind(ph2,"Cough (N=413)"=phen_all_use[,16])

st<-read.csv('st_123.csv',row.names=1)


ph2<-cbind(ph2,as.character(phen_all_use[,381]))
ph2[,5]<-ifelse(is.na(ph2[,5])==F,ifelse(ph2[,5]=="Yes",1,0),NA)
ph2[,6]<-ifelse(is.na(ph2[,6])==F,ifelse(ph2[,6]=="Yes",1,0),NA)

colnames(ph2)[6]<-c("Diarrhea (N=85)")

#Chronic conditions


phys2<-array(NA,dim=c(dim(phen_all_use)[1],14))
ind_list3<-c(20,21,23,26:29,136:142)
colnames(phen_all_use)[ind_list3]
for(i in 1:dim(phys2)[1]){
  for(j in 1:length(ind_list3)){
    if(is.na(phen_all_use[i,ind_list3[j]])==F){
      phys2[i,j]<-1
    }else if(is.na(phen_all_use[i,ind_list3[j]])==T){
      phys2[i,j]<-0
    }
  }
}


colnames(phys2)<-c("Diabetes","Allergies","Heart disease","Asthma","Stomach illness","Intestinal illness","Arthritis","Pain killers","Antibiotics","Anti-diarrheal","Anti-parasitic","Anti-fungal","Vitamins","Anti-hypertensive")

phys3<-cbind(cbind(phys,ph2),phys2)
phys3<-as.data.frame(phys3)

phen_all_raw<-phys3

##Personalities

ind_list7<-c(143,151:152)
length(which(phen_all_use[,143]%in%c("Agree a little","Agree strongly")))
length(which(phen_all_use[,151]%in%c("Agree a little","Agree strongly")))
length(which(phen_all_use[,152]%in%c("Agree a little","Agree strongly")))
il4<-c("Reserved (N=1666)","Nervous (N=950)","Openess (N=1550)")
an2<-array(NA,dim=c(dim(phen_all_use)[1],length(ind_list7)))

for(i in 1:dim(an2)[1]){
  for(j in 1:(length(ind_list7))){
    if(is.na(phen_all_use[i,ind_list7[j]])==F){
      if(phen_all_use[i,ind_list7[j]]=="Disagree strongly"){
        an2[i,j]<-(-2)
      }else if(phen_all_use[i,ind_list7[j]]=="Disagree a little"){
        an2[i,j]<-(-1)
      }else if(phen_all_use[i,ind_list7[j]]=="Neither agree nor disagree"){
        an2[i,j]<-0
      }else if(phen_all_use[i,ind_list7[j]]=="Agree a little"){
        an2[i,j]<-1
      }else if(phen_all_use[i,ind_list7[j]]=="Agree strongly"){
        an2[i,j]<-2
      }
      
    }
  }
}

an2<-as.data.frame(an2)

colnames(an2)<-il4

phen_all_raw<-cbind(phen_all_raw,an2)


#Cognitive, alcohol, cigarettes, anxiety,depresssion

ph3<-cbind(as.character(phen_all_use[,154]),as.character(phen_all_use[,154]))
ph3[,1]<-ifelse(is.na(ph3[,1])==F,ifelse(ph3[,1]=="impairment",1,0),NA)
ph3[,2]<-ifelse(is.na(ph3[,2])==F,ifelse(ph3[,2]=="dementia",1,0),NA)
ph3<-cbind(ph3,as.character(phen_all_rct_use[,8]))
for(i in 1:dim(ph3)[1]){
  if(is.na(phen_all_rct_use[i,8])==F){
    if(phen_all_rct_use[i,8]=="1 or 2"){
      ph3[i,3]<-1.5
    }else if(phen_all_rct_use[i,8]=="3 or 4"){
      ph3[i,3]<-3.5
    }else if(phen_all_rct_use[i,8]=="5 or 6"){
      ph3[i,3]<-5.5
    }else if(phen_all_rct_use[i,8]=="7 to 9"){
      ph3[i,3]<-8
    }else if(phen_all_rct_use[i,8]=="10 or more"){
      ph3[i,3]<-12
    }
  }else{ph3[i,3]<-0}
}

ph3<-cbind(ph3,as.character(phen_all_rct_use[,10]))
ph3[,4]<-ifelse(is.na(ph3[,4])==F,ifelse(ph3[,4]=="Yes",1,0),0)
ph3<-cbind(ph3,as.character(phen_all_rct_use[,11]))
ph3[,5]<-ifelse(is.na(ph3[,5])==T,0,ph3[,5])

#Anxiety, depression

for(i in 1:3){
  ph3<-cbind(ph3,as.character(phen_all_use$gad7_cat))
}
ph3[,6]<-ifelse(is.na(ph3[,6])==F,ifelse((ph3[,6]=="mild"),1,0),0)
ph3[,7]<-ifelse(is.na(ph3[,7])==F,ifelse(ph3[,7]=="moderate",1,0),0)
ph3[,8]<-ifelse(is.na(ph3[,8])==F,ifelse(ph3[,8]=="severe",1,0),0)

ph3<-cbind(ph3,phen_all_use$gad7_score)

for(i in 1:3){
  ph3<-cbind(ph3,as.character(phen_all_use$phq9_cat))
}
ph3[,10]<-ifelse(is.na(ph3[,10])==F,ifelse((ph3[,10]=="mild"),1,0),0)
ph3[,11]<-ifelse(is.na(ph3[,11])==F,ifelse(ph3[,11]=="moderate",1,0),0)
ph3[,12]<-ifelse(is.na(ph3[,12])==F,ifelse((ph3[,12]=="severe")|(ph3[,12]=="moderately severe"),1,0),0)

ph3<-cbind(ph3,phen_all_use$phq9_score)

colnames(ph3)<-c("Cognitive impairment","Dementia","Alcohol daily frequency","Cigarette usage","Cigarette frequency","mild anxiety","Moderate anxiety","Severe anxiety","gad7_score","Mild depression","Moderate depression","Severe depression","phq9_score")

ph3<-as.data.frame(ph3)

phen_all_raw<-cbind(phen_all_raw,ph3)

#Animals



ind_list4<-c(1:3,7:15,5,4,16,6,17:18,20:24,26)

il4<-c("Cat","Dog","Parakeet","Rabbit","Horse","Mice","No pets","Cow","Goat","Pig","Chicken(animal)","Duck","Turkey","Sheep","Geese","No farm animals","Bat","Lizard","Monkey","Snake","Bird","Possum","Rat","Squirrel","Other Wild","No wild animals")
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

colnames(an2)<-il4[ind_list4]

phen_all_raw<-cbind(phen_all_raw,an2)

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

colnames(an2)<-il4

phen_all_raw<-cbind(phen_all_raw,an2)

#Social factors

#Living with partners, number of partners

ph4<-cbind(as.character(phen_all_rct_use$e0300),as.character(phen_all_rct_use$e0500))
ph4[,1]<-ifelse(is.na(ph4[,1])==F,ifelse(ph4[,1]=="Yes",1,0),0)
ph4[,2]<-ifelse(is.na(ph4[,2])==T,0,ph4[,2])
colnames(ph4)<-c("Living with partner","Number of partners")

#Education

for(i in 1:3){
  ph4<-cbind(ph4,as.character(phen_all_rct_use[,5]))
}

ph4[,3]<-ifelse(is.na(ph4[,3])==F,ifelse(ph4[,3]%in%c("1st grade","2nd grade","3rd grade"),1,0),0)
ph4[,4]<-ifelse(is.na(ph4[,4])==F,ifelse(ph4[,4]%in%c("4th grade","5th grade","6th grade"),1,0),0)
ph4[,5]<-ifelse(is.na(ph4[,5])==F,ifelse(ph4[,5]%in%c("Some secondary","Secondary","More than secondary"),1,0),0)

colnames(ph4)[3:5]<-c("Grades 1-3","Grades 4-6","Grades >6")

#Rest social

ph4<-cbind(ph4,phen_all_use[,c(203:206,360:363,365:368)])
for(i in 6:17){
  ph4[,i]<-ifelse(is.na(ph4[,i])==T,0,ph4[,i])
}

colnames(ph4)[6:17]<-c("Friend ties (same building)","Friend ties (different building)","Betweeness (friendship)","Transitivity (friendship)","Familial ties (same building)","Familial ties (different building)","Betweeness (familial)","Transitivity (familial)","Degree (all ties)","Clustering coefficient (all ties)","Betweeness (all ties)","Kin percentage (to third degree)")

phen_all_raw<-cbind(phen_all_raw,ph4)

#Risk factors

ph5<-cbind(cbind(as.character(phen_all_use$altruism),as.character(phen_all_use$risk_taking)),as.character(phen_all_rct_use$i3000b))
ph5[,3]<-ifelse(is.na(ph5[,3])==F,1,0)
colnames(ph5)<-c("Altruism","Risk taking","Washing hands")

#Village factors

ph5<-cbind(ph5,cbind(cbind(cbind(phen_all_use$distance_center,phen_all_use$time_to_main_road),phen_all_use$total_churches),phen_all_use$elevation))

colnames(ph5)[4:7]<-c("Distance to village center","Distance to main road","Number of churches","Altitude")

#Economic factors

ph5<-cbind(ph5,cbind(as.character(phen_all_use[,181]),as.character(phen_all_use[,182])))
for(k in 1:dim(ph5)[1]){
  if(is.na(phen_all_use[k,181])==F){
    if(phen_all_use[k,181]=="Rarely/Never"){
      ph5[k,8]<-1/100
    }else if(phen_all_use[k,181]=="At least once a month"){
      ph5[k,8]<-1/30
    }else if(phen_all_use[k,181]=="At least once a week"){
      ph5[k,8]<-1/7
    }else if(phen_all_use[k,181]=="Every day"){
      ph5[k,8]<-1
    }
  }
}

colnames(ph5)[8:9]<-c("Travel","Monthly expenditure")

# Include all household items
#ind_list8<-c(375,376,312,316,324:325,331,333,335,337,348,351,352,354,355)
ind_list8<-c(375:376,310:359)
#il6<-c("Household size","Household wealth index","TV (N=958)","No electronics (N=89)","Earth/sand floor (N=594)","Ceramic floor (N=151)","Glass windows (N=86)","Unfinished windows (N=59)","Clay/mud walls (N=1256)","Cement walls (N=402)","Concrete roof (N=37)","Sleeping rooms","Spring protected (N=1441)","Tube well (N=91)","Dug well protected (N=93)")

ph5<-cbind(ph5,phen_all_use[,ind_list8])
ph5[,12]<-ifelse(is.na(ph5[,12])==F,1,0)
ph5[,13]<-ifelse(is.na(ph5[,13])==F,1,0)
ph5[,14]<-ifelse(is.na(ph5[,14])==F,1,0)
ph5[,15]<-ifelse(is.na(ph5[,15])==F,1,0)
ph5[,16]<-ifelse(is.na(ph5[,16])==F,1,0)
ph5[,17]<-ifelse(is.na(ph5[,17])==F,1,0)
ph5[,18]<-ifelse(is.na(ph5[,18])==F,1,0)

#colnames(ph5)[10:24]<-il6

phen_all_raw<-cbind(phen_all_raw,ph5)

phen_all_raw<-cbind(phen_all_raw,phen_all_use[,c(6,4,224,225,371)])
colnames(phen_all_raw)[162:166]<-c("age","sex","batch_effect","dna_conc","sampling_date")
phen_all_raw<-cbind(phen_all_raw,"Bristol stool scale"=st)
colnames(phen_all_raw)[167]<-"Bristol stool scale"

#Separate kitchen
unique(phen_all_raw[,124])
phen_all_raw[,124]<-ifelse(is.na(phen_all_raw[,124])==F,ifelse(phen_all_raw[,124]=="Yes",1,0),NA)

data_use<-read.csv('mb_sp.csv',row.names=1)

phen_all_raw2<-array(NA,dim=c(dim(phen_all_raw)[1],dim(phen_all_raw)[2]))

for(i in 1:dim(phen_all_raw)[1]){
  for(j in 1:dim(phen_all_raw)[2]){
    if(is.na(phen_all_raw[i,j])==F){
      phen_all_raw2[i,j]<-as.numeric(as.character(phen_all_raw[i,j]))
    }
  }
}

colnames(phen_all_raw2)<-colnames(phen_all_raw)

temp<-phen_all_raw2[,151:167]

#temp2<-phen_all_use[,322]


#All

pl_perm2_all<-adonis2(t(data_use)~`Hb A1c`+`MAP`+`BMI`+`Heart rate`+`Oxygen saturation`+`Hb total`+`Poor`+`Fair`+`Very good`+`Excellent`+`Cough (N=413)`+`Diarrhea (N=85)`+`Diabetes`+`Allergies`+`Heart disease`+`Asthma`+`Stomach illness`+`Intestinal illness`+`Arthritis`+`Pain killers`+`Antibiotics`+`Anti-diarrheal`+`Anti-parasitic`+`Anti-fungal`+`Vitamins`+`Anti-hypertensive`+`Reserved (N=1666)`+`Nervous (N=950)`+`Openess (N=1550)`+`Cognitive impairment`+`Dementia`+`Alcohol daily frequency`+`Cigarette usage`+`Cigarette frequency`+`mild anxiety`+`Moderate anxiety`+`Severe anxiety`+`gad7_score`+`Mild depression`+`Moderate depression`+`Severe depression`+`phq9_score`+`Cat`+`Dog`+`Parakeet`+`No pets`+`Cow`+`Goat`+`Pig`+`Chicken(animal)`+`Duck`+`Turkey`+`Sheep`+`Geese`+`Horse`+`Rabbit`+`No farm animals`+`Mice`+`Bat`+`Lizard`+`Snake`+`Bird`+`Possum`+`Rat`+`Squirrel`+`No wild animals`+`Beans`+`Tortillas`+`Rice`+`Bread`+`Milk`+`Yogurt`+`Cream/butter`+`Cheese`+`Eggs`+`Vegetables`+`Fruits`+`Natural juice`+`Chicken`+`Beef/Pork`+`Ham/sausages/hotdog`+`Fish`+`Diet diversity score`+`Living with partner`+`Number of partners`+`Grades 1-3`+`Grades 4-6`+`Grades >6`+`Friend ties (same building)`+`Friend ties (different building)`+`Betweeness (friendship)`+`Transitivity (friendship)`+`Familial ties (same building)`+`Familial ties (different building)`+`Betweeness (familial)`+`Transitivity (familial)`+`Degree (all ties)`+`Clustering coefficient (all ties)`+`Betweeness (all ties)`+`Kin percentage (to third degree)`+`Altruism`+`Risk taking`+`Washing hands`+`Distance to village center`+`Distance to main road`+`Number of churches`+`Altitude`+`Travel`+`Monthly expenditure`+`hh_size`+`household_wealth_index`+`Electricity`+`Radio`+`Television`+`Cell.mobile.phone`+`Non.mobile.phone`+`Refrigerator`+`None`+`Wood`+`Gas`+`fuel.Electricity`+`Kerosene`+`None.fuel`+`Separate.kitchen`+`Cement.floor`+`Earth.Sand.floor`+`Ceramic.floor`+`Tiles.floor`+`Mud.bricks.floor`+`Wood.floor`+`Other.floor`+`Wooden.windows`+`Glass.windows`+`Metal.windows`+`Unfinished.windows`+`No.windows`+`Clay.mud.walls`+`Clay.brick.walls`+`Cement.walls`+`Cane.palm.trunks.walls`+`Wood.unpolished.walls`+`Wood.polished.walls`+`Discarded.materials.walls`+`No.walls`+`Other.walls`+`Plastic.roof`+`Metal.roof`+`Clay.roof`+`Thatch.palm.roof`+`Concrete.roof`+`Wood.roof`+`Other.roof`+`Sleeping.rooms`+`Spring.water.protected`+`Spring.water.unprotected`+`tube.well`+`Dug.well.protected`+`Dug.well.unprotected`+`Surface.water`+`bottle.water`+`Other.water`+`age`+`sex`+`batch_effect`+`dna_conc`+`sampling_date`+`Bristol stool scale`, data=as.data.frame(phen_all_raw2), permutations = 999,na.action=na.omit)



write.csv(pl_perm2_all,'pl_perm2_all.csv')





##########################################################################################################
#####

#Pathways

data_transformed<-read.csv('pathways_use.csv',row.names = 1,header=T)
data_transformed<-as.matrix.data.frame(data_transformed)
data_use_pwy<-as.data.frame(t(data_transformed))

#All

pl_perm2_all_pwy<-adonis2(data_use_pwy~`Hb A1c`+`MAP`+`BMI`+`Heart rate`+`Oxygen saturation`+`Hb total`+`Poor`+`Fair`+`Very good`+`Excellent`+`Cough (N=413)`+`Diarrhea (N=85)`+`Diabetes`+`Allergies`+`Heart disease`+`Asthma`+`Stomach illness`+`Intestinal illness`+`Arthritis`+`Pain killers`+`Antibiotics`+`Anti-diarrheal`+`Anti-parasitic`+`Anti-fungal`+`Vitamins`+`Anti-hypertensive`+`Reserved (N=1666)`+`Nervous (N=950)`+`Openess (N=1550)`+`Cognitive impairment`+`Dementia`+`Alcohol daily frequency`+`Cigarette usage`+`Cigarette frequency`+`mild anxiety`+`Moderate anxiety`+`Severe anxiety`+`gad7_score`+`Mild depression`+`Moderate depression`+`Severe depression`+`phq9_score`+`Cat`+`Dog`+`Parakeet`+`No pets`+`Cow`+`Goat`+`Pig`+`Chicken(animal)`+`Duck`+`Turkey`+`Sheep`+`Geese`+`Horse`+`Rabbit`+`No farm animals`+`Mice`+`Bat`+`Lizard`+`Snake`+`Bird`+`Possum`+`Rat`+`Squirrel`+`No wild animals`+`Beans`+`Tortillas`+`Rice`+`Bread`+`Milk`+`Yogurt`+`Cream/butter`+`Cheese`+`Eggs`+`Vegetables`+`Fruits`+`Natural juice`+`Chicken`+`Beef/Pork`+`Ham/sausages/hotdog`+`Fish`+`Diet diversity score`+`Living with partner`+`Number of partners`+`Grades 1-3`+`Grades 4-6`+`Grades >6`+`Friend ties (same building)`+`Friend ties (different building)`+`Betweeness (friendship)`+`Transitivity (friendship)`+`Familial ties (same building)`+`Familial ties (different building)`+`Betweeness (familial)`+`Transitivity (familial)`+`Degree (all ties)`+`Clustering coefficient (all ties)`+`Betweeness (all ties)`+`Kin percentage (to third degree)`+`Altruism`+`Risk taking`+`Washing hands`+`Distance to village center`+`Distance to main road`+`Number of churches`+`Altitude`+`Travel`+`Monthly expenditure`+`hh_size`+`household_wealth_index`+`Electricity`+`Radio`+`Television`+`Cell.mobile.phone`+`Non.mobile.phone`+`Refrigerator`+`None`+`Wood`+`Gas`+`fuel.Electricity`+`Kerosene`+`None.fuel`+`Separate.kitchen`+`Cement.floor`+`Earth.Sand.floor`+`Ceramic.floor`+`Tiles.floor`+`Mud.bricks.floor`+`Wood.floor`+`Other.floor`+`Wooden.windows`+`Glass.windows`+`Metal.windows`+`Unfinished.windows`+`No.windows`+`Clay.mud.walls`+`Clay.brick.walls`+`Cement.walls`+`Cane.palm.trunks.walls`+`Wood.unpolished.walls`+`Wood.polished.walls`+`Discarded.materials.walls`+`No.walls`+`Other.walls`+`Plastic.roof`+`Metal.roof`+`Clay.roof`+`Thatch.palm.roof`+`Concrete.roof`+`Wood.roof`+`Other.roof`+`Sleeping.rooms`+`Spring.water.protected`+`Spring.water.unprotected`+`tube.well`+`Dug.well.protected`+`Dug.well.unprotected`+`Surface.water`+`bottle.water`+`Other.water`+`age`+`sex`+`batch_effect`+`dna_conc`+`sampling_date`+`Bristol stool scale`, data=as.data.frame(phen_all_raw2), permutations = 999,na.action=na.omit)

write.csv(pl_perm2_all_pwy,'pl_perm2_all_pwy.csv')

save.image(file="all_vars_sp_pwy_123.RData")




