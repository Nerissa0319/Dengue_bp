library(survival)
library(kableExtra)
library(tidyverse)
library(survminer)
library(ipw)
library(lme4)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(grid)
library(gridtext)
library(stargazer)
library(broom)
library(RColorBrewer)
library(powerSurvEpi)
library(nnet)

# write cox model summary to csv files
writesummary = function(modelcox,filename){
  coef_90 = confint(modelcox,level=0.9)
  hr_90 = data.frame(exp(coef_90))
  summ = summary(modelcox)
  call1 = summ$call
  n = summ$n
  nevent = summ$nevent
  newline = data.frame('n',n,'number of events',nevent)
  newline = rbind(newline,rep('',dim(newline)[2]))
  write.table(newline,file=filename,append=FALSE,sep=',',row.names=FALSE,col.names = FALSE)
  coef1 = summ$coefficients
  newdf = data.frame('V'=rownames(coef1),coef1)
  newdf = rbind(newdf,rep('',dim(newdf)[2]))
  newdf = cbind(newdf,'Vn'=rep('',dim(newdf)[1]))
  write.table(newdf,file=filename,append=TRUE,sep=',',row.names=F,col.names = TRUE)
  confint1 = summ$conf.int
  newdf = data.frame('v'=rownames(confint1),confint1)
  newdf$'Lower90' = hr_90$X5..
  newdf$'Upper90' = hr_90$X95..
  newdf = rbind(newdf,rep('',dim(newdf)[2]))
  write.table(newdf,file=filename,append=TRUE,sep=',',row.names=F,col.names = TRUE)
  concordance1 = summ$concordance
  newline = data.frame('Concordance',concordance1[[1]],'se',concordance1[[2]])
  newline = rbind(newline,rep('',dim(newline)[2]))
  write.table(newline,file=filename,append=TRUE,sep=',',row.names=FALSE,col.names = FALSE)
  llt = summ$logtest
  newline = data.frame('Likelihood ratio test',llt[[1]],'df',llt[[2]],'pvalue',llt[[3]])
  newline = rbind(newline,rep('',dim(newline)[2]))
  write.table(newline,file=filename,append=TRUE,sep=',',row.names=FALSE,col.names = FALSE)
  wt = summ$waldtest
  newline = data.frame('Wald test',wt[[1]],'df',wt[[2]],'pvalue',wt[[3]])
  newline = rbind(newline,rep('',dim(newline)[2]))
  write.table(newline,file=filename,append=TRUE,sep=',',row.names=FALSE,col.names = FALSE)
  st = summ$sctest
  newline = data.frame('Score (logrank) test',st[[1]],'df',st[[2]],'pvalue',st[[3]])
  newline = rbind(newline,rep('',dim(newline)[2]))
  write.table(newline,file=filename,append=TRUE,sep=',',row.names=FALSE,col.names = FALSE)
  if(summ$used.robust){
    rob = summ$robscore
    newline = data.frame('Robust',rob[[1]],'df',rob[[2]],'pvalue',rob[[3]])
    newline = rbind(newline,rep('',dim(newline)[2]))
    write.table(newline,file=filename,append=TRUE,sep=',',row.names=FALSE,col.names = FALSE)
  }
}

# cox model
cox_model = function(df,formula_n,grp,summarytitle,factorized,weight_method){ #factorized means exposure is factor variable
  # df = tmerged_df_basic
  # formula_n = formula_basic
  # grp = var
  # summarytitle = 'Main'
  # factorized=factorized
  # weight_method='unweighted'
  # 
  if(grepl('SBP90',grp)){legends=c('SBP >= 90','SBP < 90')
  }else if(grepl('SBP_40',grp)){legends=c('SBP decrease <= 40','SBP decrease > 40')
  }else if(grepl('NPP',grp)){legends=c('No Narrow Pulse Pressure','Narrow Pulse Pressure')
  }else if(grepl('modified',grp)){legends=c('MSI<0.7','0.7<=MSI<1.3','MSI>=1.3')[na.omit(unique(df$exposure))]
  }else if(grepl('diastolic',grp)){legends=c('DSI<1.16','DSI>=1.16')
  }else{legends=c('SI<0.6','0.6<=SI<1.0','1.0<=SI<1.4')[na.omit(unique(df$exposure))]}

  if(weight_method=='ipw'){
    cox = coxph(as.formula(formula_n),data=df,weights=ipw)#weights=ow_fluid)
  }else if(weight_method=='ow'){
    cox = coxph(as.formula(formula_n),data=df,weights=ow)#weights=ow_fluid)
  }else{
    cox = coxph(as.formula(formula_n),data=df)#weights=ow_fluid)
  }
  if(factorized){
    if(weight_method=='ipw'){
      km_fit = survfit(Surv(tstart,tstop,event)~exposure,data=df,conf.int=0.95,weights=ipw)
    }else if(weight_method=='ow'){
      km_fit = survfit(Surv(tstart,tstop,event)~exposure,data=df,conf.int=0.95,weights=ow)#conf.int = 0.9,weights=ow_fluid)
    }else{
      km_fit = survfit(Surv(tstart,tstop,event)~exposure,data=df,conf.int=0.95)#conf.int = 0.9,weights=ow_fluid)
    }
    
  }else{
    if(weight_method=='ipw'){
      km_fit = survfit(Surv(tstart,tstop,event)~1,data=df,conf.int=0.95,weights=ipw)
    }else if(weight_method=='ow'){
      km_fit = survfit(Surv(tstart,tstop,event)~1,data=df,conf.int=0.95,weights=ow)#conf.int = 0.9,weights=ow_fluid)
    }else{
      km_fit = survfit(Surv(tstart,tstop,event)~1,data=df,conf.int=0.95)#conf.int = 0.9,weights=ow_fluid)
    }
  }
  # if(grp=='Low'){palettes=c(brewer.pal(6,'Set2')[1],'pink')
  # }else if(grp=='Medium'){palettes=c(brewer.pal(6,'Set2')[1],'#CB181C')
  # }else if(grp=='High'){palettes=c(brewer.pal(6,'Set2')[1],'#67000D')
  # }else {palettes=c(brewer.pal(6,'Set2')[1],'red')}
  p_km = ggsurvplot(km_fit,
                    conf.int=TRUE,
                    # conf.int.alpha=0.2,
                    legend.labs=legends,
                    ggtheme=theme_classic(),
                    data=df
                    # palette=palettes)
  )
  p_risk = ggsurvplot(km_fit,
                      conf.int=TRUE,
                      conf.int.alpha=0.2,
                      legend.labs=legends,
                      ggtheme=theme_classic(),
                      data=df,
                      # palette=palettes,
                      fun = 'event')
  if(weight_method=='ipw'){
    writesummary(cox,paste('Research Materials/BP/output/coxsummary/ipw/',grp,'_',summarytitle,'.csv',sep=''))
  }else if(weight_method=='ow'){
    writesummary(cox,paste('Research Materials/BP/output/coxsummary/ow/',grp,'_',summarytitle,'.csv',sep=''))
  }else{
    writesummary(cox,paste('Research Materials/BP/output/coxsummary/unweighted/',grp,'_',summarytitle,'.csv',sep=''))
  }
  
  returned_ls = list('model'=cox,'kmfit'=km_fit,'km_plot'=p_km,'risk_plot'=p_risk)
  return(returned_ls)
} 
# subgroup analyisis
subgroup_cox = function(data,group,volumelevel,factorized){
  # data=tmerged_df_basic
  # group='age_ow'
  # volumelevel=var
  # factorized=factorized
  if(grepl('sen',group)){
    formulax = 'Surv(tstart,tstop,event)~factor(exposure)'
  }else if(grepl('ow',group)){
    formulax =  'Surv(tstart,tstop,event)~factor(exposure)'
  }else{
    formulax = "Surv(tstart,tstop,event)~exposure + age + factor(gender) + factor(CCMI) + Abdominal_Pain_Tenderness + 
  Clinical_fluid_accumulation + Mucosal_Bleed + Lethargy + Hepatomegaly + HCT + Hypotension_for_age + Persistent_vomiting + Temperature + Pulse_rate + 
  Respiratory_rate + white_cell_count +neutrophils_polymorphs + monocytes + basophils + eosinophils +  haematocrit + haemoglobin + platelet" #neutrophils_polymorphs + monocytes + basophils + eosinophils + 
  }
  
  if(grepl('age',group)){
    df1 = data[data$age<55,]
    df2 = data[data$age>=55,]
    if(grepl('sen',group)){
      summtitles = list('55minus_sen','55plus_sen')
    }else{
      summtitles = list('55minus','55plus')
    }
    formula1 = gsub("\\+age",'',formulax)
  }else if(grepl('gender',group)){
    df1 = data[data$gender==0,]
    df2 = data[data$gender==1,]
    if(grepl('sen',group)){
      summtitles = list('Female_sen','Male_sen')
    }else{
      summtitles = list('Female','Male')
    }
    formula1 = gsub("\\+gender",'',formulax)
  }
  
  kmplots=list()
  result = list()
  count = 0
  for(sub_df in list(df1,df2)){
    count = count+1
    subresult = cox_model(sub_df,formula1,var,summtitles[[count]],factorized,'ow')
    kmplots[[count]]=subresult$km_plot
    result[[count]] = subresult

  }
  return(result)
}

dengue = read.csv('Research Materials/BP/data/dengue_bp.csv',header=TRUE)
colnames(dengue)[1] = "Study_Number"
merged_bp = read.csv('Research Materials/BP/data/merged_bp.csv',header=TRUE)
merged_bp = merged_bp[,c(2:45)]
merged_bp$CCMI = relevel(factor(merged_bp$CCMI),ref='0')
dengue$CCMI = relevel(factor(dengue$CCMI),ref='0')
dengue$SD_day = NA
for (ind in dengue$Study_Number) {
  if (is.na(dengue[dengue$Study_Number==ind,'days_to_SD'])==FALSE && dengue[dengue$Study_Number==ind,'days_to_SD']>0){
    dengue[dengue$Study_Number==ind,'SD_day']=dengue[dengue$Study_Number==ind,'days_to_SD']+1
  }else{
    dengue[dengue$Study_Number==ind,'SD_day']=dengue[dengue$Study_Number==ind,'stay']+1
  }
}


#4.	Postural hypotension = SBP drops by 20mmHg upon standing from lying, 
# or DBP drops by 10mmHg upon standing from lying 
#5.	And then the various shock indices have their own cut-offs in the literature 
# (shock index cut offs are 0.6- <1.0, 1.0-<1.4, >1.4, 
# modified shock index <0.7 and >1.3 etc – some places use slightly different. 
# HR/DBP maybe 1.16 as cut-off or there-abouts) 

merged_bp$shock_ind = merged_bp$Pulse_rate/merged_bp$SBP
merged_bp$shock_lvl = cut(merged_bp$shock_ind,
                              breaks = c(-Inf,0.6,1.0,1.4,Inf),
                              labels = c(1,2,3,4),
                              right=FALSE)
merged_bp$diastolic_shock_ind = merged_bp$Pulse_rate/merged_bp$DBP
merged_bp$diastolic_shock_lvl = cut(merged_bp$diastolic_shock_ind,
                              breaks = c(-Inf,1.16,Inf),
                              labels = c(1,2),
                              right=FALSE)
# mean arterial pressure (map) = dbp + 1/3(sbp-dbp)
merged_bp$MAP = merged_bp$DBP + 1/3*(merged_bp$SBP - merged_bp$DBP)
merged_bp$modified_shock_ind = merged_bp$Pulse_rate/merged_bp$MAP
merged_bp$modified_shock_lvl = cut(merged_bp$modified_shock_ind,
                              breaks = c(-Inf,0.7,1.3,Inf),
                              labels = c(1,2,3),
                              right=FALSE)
write_csv(merged_bp,'Research Materials/BP/data/merged_bp.csv')

# outliers of hospitalization period
outlier = boxplot.stats(dengue$stay)$out
to_drop = dengue[dengue$stay >= min(outlier),]$Study_Number
new_dengue = dengue[! dengue$Study_Number %in% to_drop,]
new_bp = merged_bp[! merged_bp$Study_Number %in% to_drop,]
new_bp = new_bp[,c(1,2,4,6,7,8,9,11:24,26,27,30:40,45:51)]

# time independent data


df_time_ind = 
  tmerge(data1 = new_dengue[,c(1,2,3,4,5,41,42,43,44,70)],
         data2 = new_dengue[,c(1,2,3,4,5,41,42,43,44,70)],
         id = Study_Number,
         event = event(SD_day,outcome))
df_time_ind$tstart = df_time_ind$tstart + 1

ctn_tmerged_df = 
  tmerge(data1=df_time_ind,
         data2=new_bp,
         id=Study_Number,
         Abdominal_Pain_Tenderness=tdc(Day,Abdominal_Pain_Tenderness),
         Clinical_fluid_accumulation=tdc(Day,Clinical_fluid_accumulation),
         Mucosal_Bleed=tdc(Day,Mucosal_Bleed),
         Lethargy=tdc(Day,Lethargy),
         Hepatomegaly=tdc(Day,Hepatomegaly),
         HCT=tdc(Day,HCT),
         Hypotension_for_age=tdc(Day,Hypotension_for_age),
         Persistent_vomiting = tdc(Day,Persistent_vomiting),
         Temperature = tdc(Day,Temperature),
         SBP=tdc(Day,SBP),
         DBP = tdc(Day,DBP),
         SBP90=tdc(Day,SBP90),
         SBP_decrease40=tdc(Day,SBP_decrease40),
         Narrow_Pulse_Pressure_bp=tdc(Day,Narrow_Pulse_Pressure_bp),
         shock_lvl = tdc(Day,shock_lvl),
         modified_shock_lvl = tdc(Day,modified_shock_lvl),
         diastolic_shock_lvl = tdc(Day,diastolic_shock_lvl),
         Pulse_rate = tdc(Day,Pulse_rate),
         Respiratory_rate = tdc(Day,Respiratory_rate),
         # Oxygen_saturation = tdc(Day,Oxygen_saturation), missing value 64.87%
         white_cell_count=tdc(Day,white.cell.count),
         neutrophils_polymorphs=tdc(Day,Neutrophils.polymorphs),
         lymphocytes=tdc(Day,lymphocytes),
         monocytes=tdc(Day,monocytes),
         basophils=tdc(Day,basophils),
         eosinophils=tdc(Day,eosinophils),
         # atypical_reactive_lymphocytes=tdc(Day,atypical.reactive.lymphocytes), missing value 63.68%
         haematocrit=tdc(Day,Haematocrit),
         haemoglobin=tdc(Day,Haemoglobin),
         platelet=tdc(Day,Platelet))

sbp90_df_basic = 
  tmerge(data1=df_time_ind,
         data2=new_bp,
         id=Study_Number,
         exposure=tdc(Day,SBP90))


narrow_df_basic = 
  tmerge(data1=df_time_ind,
         data2=new_bp,
         id=Study_Number,
         exposure=tdc(Day,Narrow_Pulse_Pressure_bp))

sbpDecrease40_df_basic = 
  tmerge(data1=df_time_ind,
         data2=new_bp,
         id=Study_Number,
         exposure=tdc(Day,SBP_decrease40))

shock_basic = tmerge(data1=df_time_ind,
                     data2=new_bp,
                     id=Study_Number,
                     exposure=tdc(Day,shock_lvl))

diastolic_shock_basic = tmerge(data1=df_time_ind,
                         data2=new_bp,
                         id=Study_Number,
                         exposure=tdc(Day,diastolic_shock_lvl))

modified_shock_basic=tmerge(data1=df_time_ind,
                            data2=new_bp,
                            id=Study_Number,
                            exposure=tdc(Day,modified_shock_lvl))


convert_exposure = function(df,exposure_name){
  if(exposure_name=='SBP_40'){exposure_name='SBP_decrease40'}
  if(exposure_name=='NPP'){exposure_name='Narrow_Pulse_Pressure_bp'}
  exposure_col = subset(df,select=exposure_name)[,1]
  print(exposure_name)
  exposure_col
  df <- df[, -which(names(df) == exposure_name)]
  df$exposure=NA
  df$exposure = exposure_col
  return(df)
}

dengue_sbp90 = read.csv('Research Materials/BP/data/dengue_sbp90.csv',header=TRUE)
dengue_sbp90$exposure = dengue_sbp90$SBP90
colnames(dengue_sbp90)[1] = "Study_Number"
dengue_sbp40 = read.csv('Research Materials/BP/data/dengue_sbp40.csv',header=TRUE)
dengue_sbp40$exposure = dengue_sbp40$SBP_decrease40
colnames(dengue_sbp40)[1] = "Study_Number"
dengue_narrow = read.csv('Research Materials/BP/data/dengue_narrow.csv',header=TRUE)
dengue_narrow$exposure = dengue_narrow$Narrow_Pulse_Pressure_bp
colnames(dengue_narrow)[1] = "Study_Number"
dengue_shock = read.csv('Research Materials/BP/data/dengue_shock.csv',header=TRUE)
dengue_shock$exposure = cut(dengue_shock$shock_ind,
                            breaks = c(-Inf,0.6,1.0,1.4,Inf),
                            labels = c(1,2,3,4),
                            right=FALSE)
colnames(dengue_shock)[1] = "Study_Number"

dengue_modified = read.csv('Research Materials/BP/data/dengue_modified.csv',header=TRUE)
dengue_modified$exposure = cut(dengue_modified$modified_shock_ind,
                               breaks = c(-Inf,0.7,1.3,Inf),
                               labels = c(1,2,3),
                               right=FALSE)
colnames(dengue_modified)[1] = "Study_Number"
dengue_diastolic = read.csv('Research Materials/BP/data/dengue_diastolic.csv',header=TRUE)
dengue_diastolic$exposure = cut(dengue_diastolic$diastolic_shock_ind,
                                breaks = c(-Inf,1.16,Inf),
                                labels = c(1,2),
                                right=FALSE)
colnames(dengue_diastolic)[1] = "Study_Number"


ipweights = function(df,level){
  ps_model = glm(factor(exposure)~age+factor(gender)+factor(CCMI)+
                   factor(Abdominal_Pain_Tenderness)+factor(Clinical_fluid_accumulation)+
                   factor(Mucosal_Bleed)+factor(Lethargy)+factor(Hepatomegaly)+
                         factor(HCT)+factor(Hypotension_for_age)+factor(Persistent_vomiting)+
                   Temperature+Pulse_rate+Respiratory_rate+white.cell.count+
                   Neutrophils.polymorphs+lymphocytes+monocytes+basophils+eosinophils+
                   Haematocrit+Haemoglobin+Platelet,family=binomial(link=logit),data=df)
  ps_scores <- predict(ps_model, type = 'response')
  ipw <- ifelse(df$exposure == 1, 
                      1 / ps_scores, 
                      1 / (1 - ps_scores))
  # cap = quantile(ipw,0.99,na.rm=TRUE)
  
  ow = ifelse(df$exposure == 1,
                    1 - ps_scores,
                    ps_scores)
  df$ipw = ipw
  # df$ipw[df$ipw > cap] = cap
  df$ow = ow
  return(df)
}

multiipw = function(data){
  fit = multinom(exposure~age+factor(gender)+factor(CCMI)+
                   factor(Abdominal_Pain_Tenderness)+factor(Clinical_fluid_accumulation)+
                   factor(Mucosal_Bleed)+factor(Lethargy)+factor(Hepatomegaly)+
                   factor(HCT)+factor(Hypotension_for_age)+factor(Persistent_vomiting)+
                   Temperature+Pulse_rate+Respiratory_rate+white.cell.count+
                   Neutrophils.polymorphs+lymphocytes+monocytes+basophils+eosinophils+
                   Haematocrit+Haemoglobin+Platelet,data=data)
  ps = predict(fit,newdata=df,'probs')
  ipw=rep(NA,nrow(data))
  for(i in 1:nrow(data)){
    actual_exposure = data$exposure[i]
    ipw[i]=1/ps[i,actual_exposure]
  }
  data$ipw = ipw
  # cap = quantile(ipw,0.99,na.rm=TRUE)
  eInv = 1/ps
  h=as.numeric(1/rowSums(eInv))
  ow = eInv * h
  selected_weights = rep(NA,nrow(data))
  for(i in 1:nrow(data)){
    actual_exposure = data$exposure[i]
    selected_weights[i] = ow[i,actual_exposure]
  }
  data$ipw = ipw
  # data$ipw[data$ipw > cap] = cap
  data$ow = selected_weights
  return(data)
}

dengue_sbp90 = ipweights(dengue_sbp90)
dengue_sbp40 = ipweights(dengue_sbp40)
dengue_narrow = ipweights(dengue_narrow)
dengue_modified=multiipw(dengue_modified)
dengue_shock = multiipw(dengue_shock)
dengue_diastolic=ipweights(dengue_diastolic)
merged_bp$ipw_90 = NA
merged_bp$ow_90 = NA
merged_bp$ipw_40 = NA
merged_bp$ow_40 = NA
merged_bp$ipw_narrow = NA
merged_bp$ow_narrow = NA
merged_bp$ipw_shock=NA
merged_bp$ow_shock=NA
merged_bp$ipw_diastolic=NA
merged_bp$ow_diastolic=NA
merged_bp$ipw_modified=NA
merged_bp$ow_modified=NA
for(sn in unique(merged_bp$Study_Number)){
  merged_bp$ipw_90[merged_bp$Study_Number==sn]=dengue_sbp90$ipw[dengue_sbp90$Study_Number==sn]
  merged_bp$ow_90[merged_bp$Study_Number==sn]=dengue_sbp90$ow[dengue_sbp90$Study_Number==sn]
  merged_bp$ipw_40[merged_bp$Study_Number==sn]=dengue_sbp40$ipw[dengue_sbp40$Study_Number==sn]
  merged_bp$ow_40[merged_bp$Study_Number==sn]=dengue_sbp40$ow[dengue_sbp40$Study_Number==sn]
  merged_bp$ipw_narrow[merged_bp$Study_Number==sn]=dengue_narrow$ipw[dengue_narrow$Study_Number==sn]
  merged_bp$ow_narrow[merged_bp$Study_Number==sn]=dengue_narrow$ow[dengue_narrow$Study_Number==sn]
  merged_bp$ipw_shock[merged_bp$Study_Number==sn]=dengue_shock$ipw[dengue_shock$Study_Number==sn]
  merged_bp$ow_shock[merged_bp$Study_Number==sn]=dengue_shock$ow[dengue_shock$Study_Number==sn]
  merged_bp$ipw_diastolic[merged_bp$Study_Number==sn]=dengue_diastolic$ipw[dengue_diastolic$Study_Number==sn]
  merged_bp$ow_diastolic[merged_bp$Study_Number==sn]=dengue_diastolic$ow[dengue_diastolic$Study_Number==sn]
  merged_bp$ipw_modified[merged_bp$Study_Number==sn]=dengue_modified$ipw[dengue_modified$Study_Number==sn]
  merged_bp$ow_modified[merged_bp$Study_Number==sn]=dengue_modified$ow[dengue_modified$Study_Number==sn]
}

for(var in c('SBP90','SBP_40','NPP','shock','diastolic','modified')){
  var='modified'
  factorized=TRUE
  var1=var
  if(var=='shock'){var1='shock_lvl'
  }else if(var=='diastolic'){var1='diastolic_shock_lvl'
  }else if(var=='modified'){var1='modified_shock_lvl'}
  tmerged_df_adjusted=convert_exposure(ctn_tmerged_df,var1)
  if(var=='SBP'){
    other_var = 'DBP'
    tmerged_df_basic = convert_exposure(ctn_tmerged_df,var1)
    factorized=FALSE
    }else if(var=='DBP'){
      other_var='SBP'
      tmerged_df_basic=convert_exposure(ctn_tmerged_df,var)
      factorized=FALSE
      }else if(var=='SBP90'){
    weighted_df = dengue_sbp90
    tmerged_df_basic=sbp90_df_basic
  }else if(var=='SBP_40'){
    tmerged_df_basic=sbpDecrease40_df_basic
    weighted_df = dengue_sbp40
  }else if(var=='NPP'){
    tmerged_df_basic=narrow_df_basic
    weighted_df = dengue_narrow
  }else if(var=='shock'){
    tmerged_df_basic=shock_basic
    weighted_df = dengue_shock
  }else if(var=='modified'){
    tmerged_df_basic=modified_shock_basic
    weighted_df = dengue_modified
  }else if(var=='diastolic'){
    tmerged_df_basic=diastolic_shock_basic
    weighted_df = dengue_diastolic
  }
  tmerged_df_basic$ipw = NA
  tmerged_df_basic$ow = NA
  tmerged_df_adjusted$ipw = NA
  tmerged_df_adjusted$ow = NA
  for(sn in tmerged_df_basic$Study_Number){
    tmerged_df_basic$ipw[tmerged_df_basic$Study_Number==sn]=weighted_df$ipw[weighted_df$Study_Number==sn]
    tmerged_df_basic$ow[tmerged_df_basic$Study_Number==sn]=weighted_df$ow[weighted_df$Study_Number==sn]
  }
  for(sn in tmerged_df_adjusted$Study_Number){
    tmerged_df_adjusted$ipw[tmerged_df_adjusted$Study_Number==sn]=weighted_df$ipw[weighted_df$Study_Number==sn]
    tmerged_df_adjusted$ow[tmerged_df_adjusted$Study_Number==sn]=weighted_df$ow[weighted_df$Study_Number==sn]
  }
  
  formula_basic = "Surv(tstart,tstop,event)~exposure"
  basic_model = cox_model(tmerged_df_basic,formula_basic,var,'Main',factorized,'unweighted')
  
  formula_adjusted = "Surv(tstart,tstop,event)~exposure + age + factor(gender) + CCMI + Abdominal_Pain_Tenderness + 
  Clinical_fluid_accumulation + Mucosal_Bleed + Lethargy + Hepatomegaly + HCT + Hypotension_for_age + Persistent_vomiting + Temperature + Pulse_rate + 
  Respiratory_rate + white_cell_count +neutrophils_polymorphs + monocytes + basophils + eosinophils +  haematocrit + haemoglobin + platelet" #neutrophils_polymorphs + monocytes + basophils + eosinophils + 
  adjusted_model = cox_model(tmerged_df_adjusted,formula_adjusted,paste(var,'adjusted',sep='_'),'Main',factorized,'unweighted')

  if(var != 'SBP' & var != 'DBP'){
    basic_model_ipw = cox_model(tmerged_df_basic,formula_basic,var,'Main',factorized,'ipw')
    basic_model_ow = cox_model(tmerged_df_basic,formula_basic,var,'Main',factorized,'ow')
    age_model_ow = subgroup_cox(tmerged_df_basic,'age_ow',var,factorized)
    gender_model_ow = subgroup_cox(tmerged_df_basic,'gender_ow',var,factorized)
    # adjusted_model_weight = cox_model(tmerged_df_adjusted,formula_adjusted,paste(var,'doublyRobust',sep='_'),'Main',factorized,'ipw')
    # zph = cox.zph(basic_model_ow$model)
    # zph
    # plot(zph[1],lwd=2,main = paste(var))
    # abline(0,0,col=1,lty=3,lwd=2)
    # abline(h=basic_model_ow$coefficients[1],col=3,lwd=2,lty=2)
    # legend('topright',legend=c('Reference line for null effect','Average hazard over time','Time-varying hazard'),lty=c(3,2,1),col=c(1,3,1),lwd=2)
    # 
  }
}

df_time_ind$sbp90_duration = NA
df_time_ind$sbp40_duration = NA
df_time_ind$npp_duration = NA
df_time_ind$sbp90_ow = NA
df_time_ind$sbp40_ow = NA
df_time_ind$npp_ow = NA

for(sn in df_time_ind$Study_Number){
  temp_df = ctn_tmerged_df[ctn_tmerged_df$Study_Number==sn,]
  df_time_ind$sbp90_duration[df_time_ind$Study_Number==sn] = sum(temp_df$SBP90,na.rm = TRUE)
  df_time_ind$sbp40_duration[df_time_ind$Study_Number==sn] = sum(temp_df$SBP_decrease40,na.rm = TRUE)
  df_time_ind$npp_duration[df_time_ind$Study_Number==sn] = sum(temp_df$Narrow_Pulse_Pressure_bp,na.rm = TRUE)
  df_time_ind$sbp90_ow[df_time_ind$Study_Number==sn] = dengue_sbp90$ow[dengue_sbp90$Study_Number==sn]
  df_time_ind$sbp40_ow[df_time_ind$Study_Number==sn] = dengue_sbp40$ow[dengue_sbp40$Study_Number==sn]
  df_time_ind$npp_ow[df_time_ind$Study_Number==sn] = dengue_narrow$ow[dengue_narrow$Study_Number==sn]
  }

coxph(Surv(SD_day,event)~sbp90_duration,data=df_time_ind,weights=sbp90_ow)
coxph(Surv(SD_day,event)~npp_duration,data=df_time_ind,weights=npp_ow)
coxph(Surv(SD_day,event)~sbp40_duration,data=df_time_ind,weights=sbp40_ow)


library(openxlsx)

# Define the base path for the files
base_path1 <- 'Research Materials/BP/output/coxsummary/unweighted'
base_path2 <- 'Research Materials/BP/output/coxsummary/ipw'
base_path3 <- 'Research Materials/BP/output/coxsummary/ow'
# List all CSV files in the specified directory
csv_files1 <- list.files(pattern = "*.csv", path = base_path1)
csv_files2 <- list.files(pattern = "*.csv", path = base_path2)
csv_files3 <- list.files(pattern = "*.csv", path = base_path3)
# Read each CSV file and store in a list of data frames
data_list1 <- map(csv_files1, ~ read.csv(file.path(base_path1, .),header=FALSE))
data_list2 <- map(csv_files2, ~ read.csv(file.path(base_path2, .),header=FALSE))
data_list3 <- map(csv_files3, ~ read.csv(file.path(base_path3, .),header=FALSE))
# Get the sheet names by removing the '.csv' extension from the file names
sheet_names1 <- sub("\\.csv$", "", csv_files1)
sheet_names2 <- sub("\\.csv$", "", csv_files2)
sheet_names3 <- sub("\\.csv$", "", csv_files3)
# Write each data frame to a separate sheet in an Excel workbook
write.xlsx(data_list1, 
           file = file.path(base_path1, 'unweighted_Cox_Model_Summary.xlsx'), 
           sheetName = sheet_names1)
write.xlsx(data_list2, 
           file = file.path(base_path2, 'ipw_Cox_Model_Summary.xlsx'), 
           sheetName = sheet_names2)
write.xlsx(data_list3, 
           file = file.path(base_path3, 'ow_Cox_Model_Summary.xlsx'), 
           sheetName = sheet_names3)

weighted_dengue = dengue
weighted_dengue$ipw_90 = dengue_sbp90$ipw
weighted_dengue$ow_90 = dengue_sbp90$ow
weighted_dengue$ipw_40 = dengue_sbp40$ipw
weighted_dengue$ow_40 = dengue_sbp40$ow
weighted_dengue$ipw_narrow = dengue_narrow$ipw
weighted_dengue$ow_narrow = dengue_narrow$ow
write_csv(dengue_sbp90,'Research Materials/BP/data/dengue_sbp90_w.csv')
write_csv(dengue_sbp40,'Research Materials/BP/data/dengue_sbp40_w.csv')
write_csv(dengue_narrow,'Research Materials/BP/data/dengue_narrow_w.csv')
write_csv(weighted_dengue,'Research Materials/BP/data/dengue_w.csv')
write_csv(merged_bp,'Research Materials/BP/data/merged_bp_w.csv')








