#Visualization of 24 DAP embryo from Light/Dark experiment
library(readxl)
library(ggpubr)
library(tidyverse)

#read in values from SAS analysis
blues <- read_xlsx("../data/Figure2_BLUES.xlsx")
tukey <- read_xlsx("../data/Figure2_Tukey.xlsx")
tukey <- tukey[-1,]

#ensure correct formatting of values
blues$BLUE=as.numeric(blues$BLUE)
blues$SE=as.numeric(blues$SE)

#Function to input greek letters for tocochromonal naming
replace_with_greek <- function(input_string) {
  # Replace 'a' with 'α' (alpha)
  input_string <- gsub("^a", "α-", input_string)
  
  # Replace 'd' with 'δ' (delta)
  input_string <- gsub("^d", "δ-", input_string)
  
  # Replace 'g' with 'γ' (gamma)
  input_string <- gsub("^g", "γ-", input_string)
  input_string <- gsub("^B", "β-", input_string)
  
  
  input_string <- gsub("TT3", " Tocochromanol", input_string)
  input_string <- gsub("lT3", "l Tocotrienols", input_string)
  
  input_string <- gsub("T3", "Tocotrienol", input_string)
  
  input_string <- gsub("lT", "l Tocopherols", input_string)
  input_string <- gsub("T$", "Tocopherol", input_string)
  
  # Replace "T" with "Tocopherol"
  
  
  return(input_string)
}
blues$Trait <- replace_with_greek(blues$Trait)

#Formatting Tukey's table for visualization purposes
tukey <- tukey %>% separate(`Simple Effect Level`, into = c("useless", "Genotype"), sep = " ")
tukey <- tukey[,c(2,3,4,10,11)]
colnames(tukey) <- c("Genotype", "Treatment1", "Treatment2","p", "trait")
tukey$p <- as.numeric(tukey$p)
traits=unique(tukey$trait)
genotype=unique(tukey$Genotype)#same order as shown on barplot
tukey_2018=data.frame(matrix(nrow=357,ncol=4))
names(tukey_2018)=c("Trait","Genotype","Treatment","letter.assignment")
tukey_2018$Trait=traits
tukey_2018=tukey_2018[order(tukey_2018$Trait),]
tukey_2018$Genotype=genotype
tukey_2018=tukey_2018[order(tukey_2018$Trait,tukey_2018$Genotype),]
tukey_2018$Treatment=c('dark','control','light')

#Get letters denoting Tukey's p-value for visualization
for (t in traits){
  for (i in genotype){
    curr=tukey[which(tukey$trait==t & tukey$Genotype==i),]
    dc=curr$p[which(curr$Treatment1=='control'&curr$Treatment2=='dark')]
    dl=curr$p[which(curr$Treatment1=='dark'&curr$Treatment2=='light')]
    lc=curr$p[which(curr$Treatment1=='control'&curr$Treatment2=='light')]
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

tukey_2018$Trait <- replace_with_greek(tukey_2018$Trait)
all <- inner_join(blues, tukey_2018, by = c("Trait", "Genotype", "Treatment"))
traits <- unique(all$Trait)
#forgot to capitalize the "C" in carotene
all$Trait <- ifelse(all$Trait == "β-carotene", "β-Carotene", all$Trait)
traits[15] = "β-Carotene"
p_list <-list()
iter =1

#for loop to make a plot for every trait
for ( t in traits){
  curr=all[which(all$Trait==t),]
  t = paste0(t, " (ng/mg dry seed)")
  curr$Treatment=stringr::str_to_title(curr$Treatment)
  curr$Treatment=factor(curr$Treatment,levels=c('Dark','Control','Light'))
  p=ggplot(curr, aes(x=Treatment, y=BLUE, fill=Treatment)) + 
    geom_bar(stat="identity", color="black", position=position_dodge()) +  
    facet_grid(~Genotype)+
    labs(y=paste(ifelse(t == "Chlorophyll a (ng/mg dry seed)", expression("Chlorophyll ",italic(a)," (ng/mg fresh tissue)"), t)), x = NULL)+
    geom_errorbar(aes(ymin = BLUE,ymax=(BLUE+SE)), width=.2,position=position_dodge(.9))+
    geom_text(aes(label = letter.assignment, y =BLUE +SE), vjust = -0.5,size=6) +
    theme_bw()+  
    theme(plot.background = element_blank(),     
          panel.spacing.x=unit(0.5,"lines"),
          strip.text.x = element_text(size = 18,face="bold",color="black", angle =360,family='Helvetica'),
          axis.title.y = element_text(family='Helvetica', size = 18),
          axis.text.x = element_blank(),
          axis.text.y = element_text(family='Helvetica', size = 18),
          axis.ticks.x = element_blank(),
    )
  p2=p+scale_fill_manual(values=c('black','dodgerblue',"#E69F00"))
  p2 = p2+
    theme(axis.title = element_text(size =25)) 
  #Important for correct scales (ensures all letters are visible)
  if(min(curr$BLUE) > 30){
    p2 <- p2+ ylim(c(0 ,max(curr$BLUE) + 0.4 *max(curr$BLUE)))
  }else if( unique(curr$Trait) == "β-Carotene"){
    p2 <- p2+ ylim(c(0 ,max(curr$BLUE) + 0.4 *max(curr$BLUE) +1))
  }else{
    p2 <- p2+ ylim(c(0,max(curr$BLUE) + 0.4 *max(curr$BLUE) +0.2))
  }
  p_list[[iter]] = p2
  iter = iter + 1
}

#Combine separate plots with ggarrange and save figures
ggarrange(p_list[[8]],p_list[[10]],p_list[[7]],p_list[[17]],  nrow = 2, ncol = 2, common.legend = T, legend = "right")
ggsave("Figure_2.png", width = 24, height = 12, units = "in", bg ="white")

ggarrange(p_list[[1]],p_list[[2]],p_list[[3]],p_list[[7]],p_list[[4]],p_list[[5]],p_list[[6]],p_list[[8]],  nrow = 2, ncol = 4, common.legend = T, legend = "right")
ggsave("Figure_S6.png", width = 30, height = 12, units = "in", bg ="white")

ggarrange(p_list[[11]],p_list[[12]],p_list[[13]],p_list[[14]],p_list[[15]],p_list[[16]],  nrow = 2, ncol = 3, common.legend = T, legend = "right")
ggsave("Figure_S8.png", width = 26, height = 12, units = "in", bg ="white")
