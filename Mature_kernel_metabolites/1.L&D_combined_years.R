setwd("~/Desktop/L&D/pheno_results_mature_kernel")

library(gdata)
library(lme4)
library(car)
library(asreml)
library(ggplot2)
library(ggpubr)
library(emmeans)
library(broom)
library(reshape2)
library(lmerTest)

# 2018 --------------------------------------------------------------------

LD_2018=read.xls('Light and Dark_Tocochromanols[4].xlsx',sheet = 1, header = TRUE,skip=2)[,1:19]

LD_2018$Treatment= factor(LD_2018$Treatment, levels=c('Dark','Control','Light'))
LD_2018$Pedigree=factor(LD_2018$Pedigree,levels=c('B73','B97','M37W','NC358','Ki11','MS71','OH7B'))
LD_2018$Rep=factor(LD_2018$Rep)
LD_2018$ID=paste(LD_2018$Pedigree,'_',LD_2018$Treatment,'_',LD_2018$Rep,sep='')
LD_2018$ID2=paste(LD_2018$Pedigree,'_',LD_2018$Treatment,sep='')

names(LD_2018)[6]='Genotype'

LD_2018_mean=LD_2018[!duplicated(LD_2018$ID),]
LD_2018_clean=LD_2018

traits=c("a.T","d.T","g.T","a.T3","d.T3","g.T3","Total.Tocopherols","Total.Tocotrienols","Total.Tocochromanols")

all_BLUE=data.frame(matrix(nrow=21,ncol=0))
all_BLUE[,1]=c(rep('Dark',7),rep('Control',7),rep('Light',7))
all_BLUE[,2]=rep(c('B73','B97','M37W','NC358','Ki11','MS71','OH7B'),3)
names(all_BLUE)=c('Treatment','Genotype')
#loop through all traits
for(trait in traits){
  print(trait)
   
  # outlier removal ---------------------------------------------------------
  
  ##MAD
  mad_thres=3
  curr=LD_2018[,c('Genotype','Treatment','Rep',paste(trait))]
  for(i in unique(curr$Treatment)){
    for (j in unique(curr$Genotype)){
      for(x in unique(curr$Rep)){
        ind=which(curr$Treatment==i&curr$Genotype==j & curr$Rep==x)
        thres=mad(curr[ind,paste(trait)], constant=1)
        mad=abs(curr[ind,paste(trait)] - median(curr[ind,paste(trait)])) / thres
        if(length(which(mad>mad_thres))>0){
          LD_2018_clean[ind,][which(mad>mad_thres),paste(trait)]=NA
        }
      }
    }
  }
  print(length(which(is.na(LD_2018_clean[,trait]))))

  # Sub-plot average --------------------------------------------------------
  LD_2018_mean[,trait]=unique(ave(LD_2018_clean[,trait],LD_2018_clean$ID,FUN=function(x) mean(x, na.rm=T)))
}


write.csv(LD_2018_clean,'1.L&D_OR_2018.csv',quote=F,row.names = F)
write.csv(LD_2018_mean,'1.L&D_averaged_2018.csv',quote=F,row.names = F)


# 2019 --------------------------------------------------------------------

LD_2019=read.xls('2019_L&D Tocos and Carots.xlsx',sheet = 1, header = TRUE,skip=2)[1:75,c(1:23,25)]
names(LD_2019)[8]='Genotype'
names(LD_2019)[24]='Total.Carots'


LD_2019$Treatment= factor(LD_2019$Treatment, levels=c('dark','control','light'))
LD_2019$Genotype=factor(LD_2019$Genotype,levels=c('B73','B97','OH7B'))
LD_2019$Rep=factor(LD_2019$Rep)
LD_2019$ID=paste(LD_2019$Genotype,'_',LD_2019$Treatment,'_',LD_2019$Rep,sep='')
LD_2019$ID2=paste(LD_2019$Genotype,'_',LD_2019$Treatment,sep='')

LD_2019_mean=LD_2019[!duplicated(LD_2019$ID),]
LD_2019_clean=LD_2019

traits=c("a.T","d.T","g.T","a.T3","d.T3","g.T3","Total.Tocopherols","Total.Tocotrienols","Total.Tocochromanols",'Total.Carots')

#loop through all traits
for(trait in traits){
  print(trait)
  
  # outlier removal ---------------------------------------------------------
  
  ##MAD
  mad_thres=3
  curr=LD_2019[,c('Genotype','Treatment','Rep',paste(trait))]
  for(i in unique(curr$Treatment)){
    for (j in unique(curr$Genotype)){
      for(x in unique(curr$Rep)){
        ind=which(curr$Treatment==i&curr$Genotype==j & curr$Rep==x)
        thres=mad(curr[ind,paste(trait)], constant=1)
        mad=abs(curr[ind,paste(trait)] - median(curr[ind,paste(trait)])) / thres
        if(length(which(mad>mad_thres))>0){
          LD_2019_clean[ind,][which(mad>mad_thres),paste(trait)]=NA
        }
      }
    }
  }

  # Sub-plot average --------------------------------------------------------
  LD_2019_mean[,trait]=unique(ave(LD_2019_clean[,trait],LD_2019_clean$ID,FUN=function(x) mean(x, na.rm=T)))
}

write.csv(LD_2019_clean,'1.L&D_OR_2019.csv',quote=F,row.names = F)
write.csv(LD_2019_mean,'1.L&D_averaged_2019.csv',quote=F,row.names = F)
