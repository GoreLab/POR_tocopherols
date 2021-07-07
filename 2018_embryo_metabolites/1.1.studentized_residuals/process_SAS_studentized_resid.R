setwd("~/Desktop/L&D/pheno_embryo/studentized_residuals")

traits=c('aT','dT','gT','aT3','dT3','gT3','TotalT','TotalT3','TotalTT3','chlorophyll',
'Neo', 'Vio', 'Lutein', 'Zeaxanthin', 'BCarotene', 'OtherCarotenes', 'TotalCarotenoids')


# studentized residual bonferroni threshold -------------------------------
nbind = 42 #non-missing observations
fixed=3
random=2
threshold = qt(p= 1-0.05/(2*nbind), df=(nbind-(fixed+random+1)-1))

for(t in traits){
  print(t)
  curr=read.csv(paste(t,'.csv',sep=''))
  curr=curr[c(1,2,3,21)]
  curr$ID=paste(curr$Treatment,curr$Rep,curr$Genotype,sep='_')
  names(curr)[4]=t
  print(length(which(curr[,4]>threshold)))
  if(t=='aT'){all=curr}else{all=merge(all,curr[,4:5],by='ID')}
}

#no observation will be removed as outlier
write.csv(all[,-1],'consolidated_studentized_resid.csv',quote=F,row.names = F)

