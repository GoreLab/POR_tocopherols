setwd("~/Desktop/L&D/pheno_embryo/")

library(gdata)
library(ggplot2)
library(stringr)


blue=read.xls('L&D_embryo.xlsx',sheet=4)
blue=blue[1:210,1:5]#carotenoid and empty cells removed
blue$BLUE=as.numeric(blue$BLUE)
blue$SE=as.numeric(blue$SE)

traits=unique(blue$Trait)
genotype=c("B73","B97","Ki11","M37W","MS71","NC358","OH7B")#same order as shown on barplot

tukey=data.frame(matrix(nrow=210,ncol=4))
names(tukey)=c("Trait","Genotype","Treatment","letter.assignment")
tukey$Trait=rep(traits,21)
tukey=tukey[order(tukey$Trait),]
tukey$Genotype=rep(genotype,30)
tukey=tukey[order(tukey$Trait,tukey$Genotype),]
tukey$Treatment=c('dark','control','light')

# assign letters from Tukey -----------------------------------------------
tukey_res=read.csv('L&D_embryo_Tukey_result_input_for_plot.csv')

for (t in traits){
  for (i in genotype){
    curr=tukey_res[which(tukey_res$Trait==t & tukey_res$Geno==i),]
    
    dc=curr$P[which(curr$Treatment1=='control'&curr$Treatment2=='dark')]
    dl=curr$P[which(curr$Treatment1=='dark'&curr$Treatment2=='light')]
    lc=curr$P[which(curr$Treatment1=='control'&curr$Treatment2=='light')]
    
    tukey$letter.assignment[which(tukey$Trait==t &tukey$Genotype==i & tukey$Treatment=='dark')]='a' #fix dark to a
    
    if (dc<0.05){
      tukey$letter.assignment[which(tukey$Trait==t &tukey$Genotype==i & tukey$Treatment=='control')]='b' 
      if (lc < 0.05 & dl < 0.05){
        tukey$letter.assignment[which(tukey$Trait==t &tukey$Genotype==i & tukey$Treatment=='light')]='c' 
      }else if(lc < 0.05 & dl >= 0.05){
        tukey$letter.assignment[which(tukey$Trait==t &tukey$Genotype==i & tukey$Treatment=='light')]='a' 
      }else if(lc >= 0.05 & dl < 0.05){
        tukey$letter.assignment[which(tukey$Trait==t &tukey$Genotype==i & tukey$Treatment=='light')]='b'
      }else{
        tukey$letter.assignment[which(tukey$Trait==t &tukey$Genotype==i & tukey$Treatment=='light')]='ab' 
      }
    }
    if (dc>=0.05){
      tukey$letter.assignment[which(tukey$Trait==t &tukey$Genotype==i & tukey$Treatment=='control')]='a' 
      if (lc < 0.05 & dl < 0.05){
        tukey$letter.assignment[which(tukey$Trait==t &tukey$Genotype==i & tukey$Treatment=='light')]='b' 
      }else if(lc < 0.05 & dl >= 0.05){
        tukey$letter.assignment[which(tukey$Trait==t &tukey$Genotype==i & tukey$Treatment=='light')]='b'
        tukey$letter.assignment[which(tukey$Trait==t &tukey$Genotype==i & tukey$Treatment=='dark')]='ab'
      }else if(lc >= 0.05 & dl < 0.05){
        tukey$letter.assignment[which(tukey$Trait==t &tukey$Genotype==i & tukey$Treatment=='light')]='b'
        tukey$letter.assignment[which(tukey$Trait==t &tukey$Genotype==i & tukey$Treatment=='control')]='ab'
      }else{
        tukey$letter.assignment[which(tukey$Trait==t &tukey$Genotype==i & tukey$Treatment=='light')]='a' 
      }
    }
  }
}

tukey$ID=paste(tukey$Genotype,tukey$Treatment,sep='_')


# plot barplots with letter assignment ------------------------------------

for(t in traits){
  tiff(paste('embryo_BLUE_bar_tukey_2018_',t,'.png',sep=''),width=5,height=7,family='serif',units="in",res=300)
  curr=blue[which(blue$Trait==t),]
  curr$ID=paste(curr$Genotype,curr$Treatment,sep='_')
  tukey_curr=tukey[which(tukey$Trait==t),]
  all=merge(curr,tukey_curr[,c('ID','letter.assignment')],by='ID')
  all$Treatment=str_to_title(all$Treatment)
  all$Treatment=factor(all$Treatment,levels=c('Dark','Control','Light'))
  
  p=ggplot(all, aes(x=Treatment, y=BLUE, fill=Treatment)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +  
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
          panel.grid.minor = element_blank())
  
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
