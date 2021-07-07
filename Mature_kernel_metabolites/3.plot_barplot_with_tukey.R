setwd("~/Desktop/L&D/pheno_results_mature_kernel/")

library(ggplot2)
library(stringr)


blue=read.csv('L&D_mature_BLUE.csv')
blue$BLUE=as.numeric(blue$BLUE)
blue$SE=as.numeric(blue$SE)

traits=unique(blue$Trait)

tukey_res=read.csv('L&D_mature_kernel_Tukey_result_input_for_plot.csv')[1:210,]

# 2018 --------------------------------------------------------------------
genotype=c("B73","B97","Ki11","M37W","MS71","NC358","OH7B")#same order as shown on barplot
blue_2018=blue[which(blue$Year==2018),]

tukey_2018=data.frame(matrix(nrow=189,ncol=4))
names(tukey_2018)=c("Trait","Genotype","Treatment","letter.assignment")
tukey_2018$Trait=traits
tukey_2018=tukey_2018[order(tukey_2018$Trait),]
tukey_2018$Genotype=genotype
tukey_2018=tukey_2018[order(tukey_2018$Trait,tukey_2018$Genotype),]
tukey_2018$Treatment=c('dark','control','light')

# assign letters from Tukey -----------------------------------------------
tukey_res_2018=tukey_res[,c(1:4,5)]
for (t in traits){
  for (i in genotype){
    curr=tukey_res_2018[which(tukey_res_2018$Trait==t & tukey_res_2018$Geno==i),]
    
    dc=curr$X2018[which(curr$Treatment1=='Control'&curr$Treatment2=='Dark')]
    dl=curr$X2018[which(curr$Treatment1=='Dark'&curr$Treatment2=='Light')]
    lc=curr$X2018[which(curr$Treatment1=='Control'&curr$Treatment2=='Light')]
    
    tukey_2018$letter.assignment[which(tukey_2018$Trait==t &tukey_2018$Genotype==i & tukey_2018$Treatment=='dark')]='a' #fix dark to a
    
    if (dc<0.05){
      tukey_2018$letter.assignment[which(tukey_2018$Trait==t &tukey_2018$Genotype==i & tukey_2018$Treatment=='control')]='b' 
      if (lc < 0.05 & dl < 0.05){
        tukey_2018$letter.assignment[which(tukey_2018$Trait==t &tukey_2018$Genotype==i & tukey_2018$Treatment=='light')]='c' 
      }else if(lc < 0.05 & dl >= 0.05){
        tukey_2018$letter.assignment[which(tukey_2018$Trait==t &tukey_2018$Genotype==i & tukey_2018$Treatment=='light')]='a' 
      }else if(lc >= 0.05 & dl < 0.05){
        tukey_2018$letter.assignment[which(tukey_2018$Trait==t &tukey_2018$Genotype==i & tukey_2018$Treatment=='light')]='b'
      }else{
        tukey_2018$letter.assignment[which(tukey_2018$Trait==t &tukey_2018$Genotype==i & tukey_2018$Treatment=='light')]='ab' 
      }
    }
    if (dc>=0.05){
      tukey_2018$letter.assignment[which(tukey_2018$Trait==t &tukey_2018$Genotype==i & tukey_2018$Treatment=='control')]='a' 
      if (lc < 0.05 & dl < 0.05){
        tukey_2018$letter.assignment[which(tukey_2018$Trait==t &tukey_2018$Genotype==i & tukey_2018$Treatment=='light')]='b' 
      }else if(lc < 0.05 & dl >= 0.05){
        tukey_2018$letter.assignment[which(tukey_2018$Trait==t &tukey_2018$Genotype==i & tukey_2018$Treatment=='light')]='b'
        tukey_2018$letter.assignment[which(tukey_2018$Trait==t &tukey_2018$Genotype==i & tukey_2018$Treatment=='dark')]='ab'
      }else if(lc >= 0.05 & dl < 0.05){
        tukey_2018$letter.assignment[which(tukey_2018$Trait==t &tukey_2018$Genotype==i & tukey_2018$Treatment=='light')]='b'
        tukey_2018$letter.assignment[which(tukey_2018$Trait==t &tukey_2018$Genotype==i & tukey_2018$Treatment=='control')]='ab'
      }else{
        tukey_2018$letter.assignment[which(tukey_2018$Trait==t &tukey_2018$Genotype==i & tukey_2018$Treatment=='light')]='a' 
      }
    }
  }
}

tukey_2018$ID=toupper(paste(tukey_2018$Genotype,tukey_2018$Treatment,sep='_'))


# plot barplots with letter assignment ------------------------------------

for(t in traits){
  tiff(paste('mature_kernel_BLUE_bar_tukey_2018_',t,'.png',sep=''),width=5,height=7,family='serif',units="in",res=300)
  curr=blue_2018[which(blue_2018$Trait==t),]
  curr$ID=toupper(paste(curr$Genotype,curr$Treatment,sep='_'))
  tukey_curr=tukey_2018[which(tukey_2018$Trait==t),]
  all=merge(curr,tukey_curr[,c('ID','letter.assignment')],by='ID')
  all$Treatment=str_to_title(all$Treatment)
  all$Treatment=factor(all$Treatment,levels=c('Dark','Control','Light'))
  
  p=ggplot(all, aes(x=Treatment, y=BLUE, fill=Treatment)) + 
    geom_bar(stat="identity", color="black",position=position_dodge()) +  
    facet_grid(~Genotype,switch = 'x')+
    geom_errorbar(aes(ymin=(BLUE-SE), ymax=(BLUE+SE)), width=.2,position=position_dodge(.9))+
    theme_bw()+  
    labs(y=paste(t))+
    scale_fill_manual(values=c('dodgerblue','gray',"cyan3"))+
    theme(plot.background = element_blank(),     
          panel.spacing.x=unit(0,"lines"),
          strip.text.x = element_text(size = 8,face="bold",color="black",family='serif'),
          strip.background.x=element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(family='serif'),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(face="bold",family='serif'),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = 'none')
  
  # assign significance level from Tukey ------------------------------------
  all$yval=all$BLUE+all$SE
  count=1
  for (i in genotype){
      values1=data.frame(x=0.75,y=all$yval[which(all$Genotype==i & all$Treatment=='Dark')], text=all$letter.assignment[which(all$Genotype==i & all$Treatment=='Dark')],Genotype=i)
      values2=data.frame(x=1.75,y=all$yval[which(all$Genotype==i & all$Treatment=='Control')], text=all$letter.assignment[which(all$Genotype==i & all$Treatment=='Control')],Genotype=i)
      values3=data.frame(x=2.75,y=all$yval[which(all$Genotype==i & all$Treatment=='Light')], text=all$letter.assignment[which(all$Genotype==i & all$Treatment=='Light')],Genotype=i)
      if(all$letter.assignment[all$Genotype==i & all$Treatment=='Control']=='ab'){values2$x=1.5}#slight adjustment of the position of ab
      assign(paste('sigvals',3*count-2,sep=''),values1)
      assign(paste('sigvals',3*count-1,sep=''),values2)
      assign(paste('sigvals',3*count,sep=''),values3)
      
      count=count+1
  }
  
  p2=p+geom_text(data=sigvals1, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals2, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals3, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals4, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals5, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals6, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals7, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals8, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals9, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals10, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals11, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals12, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals13, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals14, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals15, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals16, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals17, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals18, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals19, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals20, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals21, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')
    
  p2=p2+scale_fill_manual(values=c('gray','dodgerblue',"cyan3"))
  
  print(p2)
  dev.off()
}



# 2019 --------------------------------------------------------------------
genotype=c("B73","B97","OH7B")#same order as shown on barplot
blue_2019=blue[which(blue$Year==2019),]

tukey_2019=data.frame(matrix(nrow=81,ncol=4))
names(tukey_2019)=c("Trait","Genotype","Treatment","letter.assignment")
tukey_2019$Trait=traits
tukey_2019=tukey_2019[order(tukey_2019$Trait),]
tukey_2019$Genotype=genotype
tukey_2019=tukey_2019[order(tukey_2019$Trait,tukey_2019$Genotype),]
tukey_2019$Treatment=c('dark','control','light')

# assign letters from Tukey -----------------------------------------------
tukey_res_2019=tukey_res[,c(1:4,6)]
for (t in traits){
  for (i in genotype){
    curr=tukey_res_2019[which(tukey_res_2019$Trait==t & tukey_res_2019$Geno==i),]
    
    dc=curr$X2019[which(curr$Treatment1=='Control'&curr$Treatment2=='Dark')]
    dl=curr$X2019[which(curr$Treatment1=='Dark'&curr$Treatment2=='Light')]
    lc=curr$X2019[which(curr$Treatment1=='Control'&curr$Treatment2=='Light')]
    
    tukey_2019$letter.assignment[which(tukey_2019$Trait==t &tukey_2019$Genotype==i & tukey_2019$Treatment=='dark')]='a' #fix dark to a
    
    if (dc<0.05){
      tukey_2019$letter.assignment[which(tukey_2019$Trait==t &tukey_2019$Genotype==i & tukey_2019$Treatment=='control')]='b' 
      if (lc < 0.05 & dl < 0.05){
        tukey_2019$letter.assignment[which(tukey_2019$Trait==t &tukey_2019$Genotype==i & tukey_2019$Treatment=='light')]='c' 
      }else if(lc < 0.05 & dl >= 0.05){
        tukey_2019$letter.assignment[which(tukey_2019$Trait==t &tukey_2019$Genotype==i & tukey_2019$Treatment=='light')]='a' 
      }else if(lc >= 0.05 & dl < 0.05){
        tukey_2019$letter.assignment[which(tukey_2019$Trait==t &tukey_2019$Genotype==i & tukey_2019$Treatment=='light')]='b'
      }else{
        tukey_2019$letter.assignment[which(tukey_2019$Trait==t &tukey_2019$Genotype==i & tukey_2019$Treatment=='light')]='ab' 
      }
    }
    if (dc>=0.05){
      tukey_2019$letter.assignment[which(tukey_2019$Trait==t &tukey_2019$Genotype==i & tukey_2019$Treatment=='control')]='a' 
      if (lc < 0.05 & dl < 0.05){
        tukey_2019$letter.assignment[which(tukey_2019$Trait==t &tukey_2019$Genotype==i & tukey_2019$Treatment=='light')]='b' 
      }else if(lc < 0.05 & dl >= 0.05){
        tukey_2019$letter.assignment[which(tukey_2019$Trait==t &tukey_2019$Genotype==i & tukey_2019$Treatment=='light')]='b'
        tukey_2019$letter.assignment[which(tukey_2019$Trait==t &tukey_2019$Genotype==i & tukey_2019$Treatment=='dark')]='ab'
      }else if(lc >= 0.05 & dl < 0.05){
        tukey_2019$letter.assignment[which(tukey_2019$Trait==t &tukey_2019$Genotype==i & tukey_2019$Treatment=='light')]='b'
        tukey_2019$letter.assignment[which(tukey_2019$Trait==t &tukey_2019$Genotype==i & tukey_2019$Treatment=='control')]='ab'
      }else{
        tukey_2019$letter.assignment[which(tukey_2019$Trait==t &tukey_2019$Genotype==i & tukey_2019$Treatment=='light')]='a' 
      }
    }
  }
}

tukey_2019$ID=toupper(paste(tukey_2019$Genotype,tukey_2019$Treatment,sep='_'))


# plot barplots with letter assignment ------------------------------------

for(t in traits){
  tiff(paste('mature_kernel_BLUE_bar_tukey_2019_',t,'.png',sep=''),width=2,height=7,family='serif',units="in",res=300)
  curr=blue_2019[which(blue_2019$Trait==t),]
  curr$ID=toupper(paste(curr$Genotype,curr$Treatment,sep='_'))
  tukey_curr=tukey_2019[which(tukey_2019$Trait==t),]
  all=merge(curr,tukey_curr[,c('ID','letter.assignment')],by='ID')
  all$Treatment=str_to_title(all$Treatment)
  all$Treatment=factor(all$Treatment,levels=c('Dark','Control','Light'))
  all$Genotype=factor(all$Genotype,levels=c("B73","B97","OH7B"))
  
  p=ggplot(all, aes(x=Treatment, y=BLUE, fill=Treatment)) + 
    geom_bar(stat="identity", color="black",position=position_dodge()) +  
    facet_grid(~Genotype,switch = 'x')+
    geom_errorbar(aes(ymin=(BLUE-SE), ymax=(BLUE+SE)), width=.2,position=position_dodge(.9))+
    theme_bw()+  
    labs(y=paste(t))+
    scale_fill_manual(values=c('dodgerblue','gray',"cyan3"))+
    theme(plot.background = element_blank(),     
          panel.spacing.x=unit(0,"lines"),
          strip.text.x = element_text(size = 8,face="bold",color="black", angle =360,family='serif'),
          strip.background.x=element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(family='serif'),
          axis.ticks.x = element_blank(),
          axis.title = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = 'none')
  
  # assign significance level from Tukey ------------------------------------
  all$yval=all$BLUE+all$SE
  count=1
  for (i in genotype){
    values1=data.frame(x=0.8,y=all$yval[which(all$Genotype==i & all$Treatment=='Dark')], text=all$letter.assignment[which(all$Genotype==i & all$Treatment=='Dark')],Genotype=i)
    values2=data.frame(x=1.8,y=all$yval[which(all$Genotype==i & all$Treatment=='Control')], text=all$letter.assignment[which(all$Genotype==i & all$Treatment=='Control')],Genotype=i)
    values3=data.frame(x=2.8,y=all$yval[which(all$Genotype==i & all$Treatment=='Light')], text=all$letter.assignment[which(all$Genotype==i & all$Treatment=='Light')],Genotype=i)
    if(all$letter.assignment[all$Genotype==i & all$Treatment=='Control']=='ab'){values2$x=1.5}#slight adjustment of the position of ab
    assign(paste('sigvals',3*count-2,sep=''),values1)
    assign(paste('sigvals',3*count-1,sep=''),values2)
    assign(paste('sigvals',3*count,sep=''),values3)
    
    count=count+1
  }
  
  p2=p+geom_text(data=sigvals1, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals2, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals3, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals4, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals5, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals6, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals7, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals8, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')+
    geom_text(data=sigvals9, aes(x=x,y=y, label=text, fill=NA), hjust=0,vjust=-1, family ='serif')
  
  p2=p2+scale_fill_manual(values=c('gray','dodgerblue',"cyan3"))
  
  print(p2)
  dev.off()
  
}
