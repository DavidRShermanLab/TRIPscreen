#################
#
# TRIP Untreated vs. INH Processing Scripts
# Shuyi Ma, March 2018
#
# This code processes the RPKM sequence abundance values from the second tab of TableS1 in supplementary materials, and
# generates as outputs the relative log2 abundance fold change for each TFI strain upon induction relative to uninduced,
# under the presence or absence of the drug isoniazid.
#
# The input file "TableS1.csv" is the .csv version of the second tab of TableS1 from the supplementary materials of the manuscript
# The optional flag "IntermediateReportFlag" will cause the code to generate intermediate tables if the user switches this variable to "TRUE"
#
#################

library(Biobase)
library(limma)
library(Rmisc)

IntermediateReportFlag = FALSE

dataraw = read.csv('TableS1.csv',as.is = T)
TFI = colnames(dataraw)[6:ncol(dataraw)]

# These TFs contain NAN Count values or Day 0 RPM values less than 5, will be ignored
tmp = colMeans(dataraw[,TFI])
badcols = which(is.na(tmp))
tmp = colMeans(dataraw[dataraw$Day == 0,TFI])
badcols = c(badcols,which(tmp < 5))
dataraw2 = dataraw[,-which(colnames(dataraw) %in% TFI[badcols])]
TFIb = TFI[-badcols]

if (IntermediateReportFlag){write.csv(dataraw2,'TRIPrawTFIfiltered.csv')}


# Note: TFs that change significantly in log phase growth without induction. These may give harder to interpret results.
mean0 = colMeans(dataraw2[dataraw2$Day == 0,TFIb])
mean7 = colMeans(dataraw2[dataraw2$Day == 7 & dataraw2$Induced == FALSE & dataraw2$Drug == "Untreated",TFIb])
meandiff = mean7 - mean0
TFIb[meandiff < -1.5]
#[1] "Rv0047c"       "Rv0880"      

uDrug = unique(dataraw2$Drug)
uConc = unique(dataraw2$ConcxMIC)
uInd = unique(dataraw2$Induced)
uDay = unique(dataraw2$Day)
uExperiment = unique(dataraw2$Experiment)

## Ratios

# Calculate the abundances for the samples collected from Day 0
tmp1 = data.frame()
tmp2 = data.frame()
for (e in 1:length(uExperiment)){
  mean0 = colMeans(dataraw2[dataraw2$Day == 0 & dataraw2$Experiment %in% uExperiment[[e]],TFIb])

  for (d in 1:length(uDrug)){
    for (a in 1:length(uInd)){
      tmp1 = rbind(tmp1,dataraw2[dataraw2$Drug == uDrug[d] & dataraw2$Induced == uInd[a] & dataraw2$Day == 7 & dataraw2$Experiment %in% uExperiment[[e]],1:5])
      tmp2 = rbind(tmp2,sweep(dataraw2[dataraw2$Drug == uDrug[d] & dataraw2$Induced == uInd[a] & dataraw2$Day == 7 & dataraw2$Experiment %in% uExperiment[[e]],TFIb],2,mean0))
    }
  }
}
dataDeltaTimegood = cbind(tmp1,tmp2)
dataDeltaTimegood$Day = NULL

if (IntermediateReportFlag){write.csv(dataDeltaTimegood,'TRIP_DeltaTime.csv')}

## Statistical Tests comparing induced, uninduced condition for each drug/concentration.
# T-test
tmp1 = data.frame()
tmp2 = data.frame()
for (d in 1:length(uDrug)){
  tmp3 = dataDeltaTimegood[dataDeltaTimegood$Drug == uDrug[d],]
  tmp1 = rbind(tmp1,tmp3[1,1:2])
  Pv = numeric()
  for (g in 1:length(TFIb)){
    tp = t.test(eval(parse(text=paste(TFIb[g],'~ Induced'))), tmp3)
    Pv = c(Pv,tp$p.value)
  }
  Pv = p.adjust(Pv, method = "BH")
  tmp2 = rbind(tmp2,Pv)
}
dataDeltaTimepvalT = cbind(tmp1,tmp2)
colnames(dataDeltaTimepvalT)[3:ncol(dataDeltaTimepvalT)] = TFIb


#Wilcoxon ranksum test
tmp1 = data.frame()
tmp2 = data.frame()
for (d in 1:length(uDrug)){
  tmp3 = dataDeltaTimegood[dataDeltaTimegood$Drug == uDrug[d],]
  tmp1 = rbind(tmp1,tmp3[1,1:2])
  Pv = numeric()
    for (g in 1:length(TFIb)){
      wp = wilcox.test(tmp3[tmp3$Induced == TRUE,TFIb[g]],tmp3[tmp3$Induced == FALSE,TFIb[g]])
      Pv = c(Pv,wp$p.value)
    }
    Pv = p.adjust(Pv, method = "BH")
    tmp2 = rbind(tmp2,Pv)
}
dataDeltaTimepvalWilcox = cbind(tmp1,tmp2)
colnames(dataDeltaTimepvalWilcox)[3:ncol(dataDeltaTimepvalWilcox)] = TFIb

sigTFutTTest = which(dataDeltaTimepvalT[dataDeltaTimepvalT$Drug == "Untreated",TFIb] < 0.05)
sigTFutWilcox = which(dataDeltaTimepvalWilcox[dataDeltaTimepvalWilcox$Drug == "Untreated",TFIb] < 0.05)
sigTFutAll = intersect(sigTFutTTest,sigTFutWilcox)

# Z-scores
tmp1c = data.frame()
tmp2c = data.frame()
tmp1d = data.frame()
tmp2d = data.frame()

for (c in 1:length(uDrug)){
  tmp3b = dataDeltaTimegood[dataDeltaTimegood$Drug == uDrug[c],]
  sd0a = apply(tmp3b[,TFIb],2,sd)
  mean0a = apply(tmp3b[tmp3b$Induced == FALSE,TFIb],2,mean)

  tmp6=tmp3b[tmp3b$Induced == TRUE,]
  mdiffa = sweep(tmp6[,TFIb],2,mean0a)
  zsa = sweep(mdiffa,2,sd0a,FUN='/')
  
  tmp1d = rbind(tmp1d,tmp6[,1:2])
  tmp2d = rbind(tmp2d,zsa)
  
  uE = unique(tmp3b$Experiment)
  for (e in 1:length(uE)){
    tmp4 = tmp3b[tmp3b$Experiment %in% uE[[e]],]
    sd0e = apply(tmp4[,TFIb],2,sd)
    mean0 = apply(tmp4[tmp4$Induced == FALSE,TFIb],2,mean)

    tmp5 = tmp4[tmp4$Induced == TRUE,]
    mdiff = sweep(tmp5[,TFIb],2,mean0)
    zs = sweep(mdiff,2,sd0e,FUN = '/')

    tmp1c = rbind(tmp1c,tmp5[,1:2])
    tmp2c = rbind(tmp2c,zs)

  }

}
dataDeltaTimezscore = cbind(tmp1c,tmp2c)
dataDeltaTimezscoreExpAvg = cbind(tmp1d,tmp2d)
colnames(dataDeltaTimezscore)[3:ncol(dataDeltaTimezscore)] = TFIb
colnames(dataDeltaTimezscoreExpAvg)[3:ncol(dataDeltaTimezscoreExpAvg)] = TFIb

tmp1 = data.frame()
tmp2 = data.frame()
tmp1d = data.frame()
tmp2d = data.frame()
for (c in 1:length(uDrug)){
  tmp4 = dataDeltaTimezscore[dataDeltaTimezscore$Drug == uDrug[c],]
  tmp4b = dataDeltaTimezscoreExpAvg[dataDeltaTimezscoreExpAvg$Drug == uDrug[c],]
  tmp1 = rbind(tmp1,tmp4[1,1:2])
  tmp2 = rbind(tmp2,apply(tmp4[,TFIb],2,mean))
  tmp1d = rbind(tmp1d,tmp4b[1,1:2])
  tmp2d = rbind(tmp2d,apply(tmp4b[,TFIb],2,mean))
}
dataDeltaTimezscoreAvg = cbind(tmp1,tmp2)
dataDeltaTimezscoreAvgExpAvg = cbind(tmp1d,tmp2d)
colnames(dataDeltaTimezscoreAvg)[3:ncol(dataDeltaTimezscoreAvg)] = TFIb
colnames(dataDeltaTimezscoreAvgExpAvg)[3:ncol(dataDeltaTimezscoreAvgExpAvg)] = TFIb

## Ratio-ing Induced vs. Uninduced

tmp1 = data.frame()
tmp2 = data.frame()
for (d in 1:length(uDrug)){
  tmp3 = dataDeltaTimegood[dataDeltaTimegood$Drug == uDrug[d],]
  
    for (e in 1:length(uExperiment)){
      
      tmp1 = rbind(tmp1,tmp3[tmp3$Induced == TRUE & tmp3$Experiment %in% uExperiment[[e]],1:4])
      
      tmp2 = rbind(tmp2,sweep(tmp3[tmp3$Induced == TRUE & tmp3$Experiment %in% uExperiment[[e]],TFIb],2,apply(tmp3[tmp3$Induced == FALSE & tmp3$Experiment %in% uExperiment[[e]],TFIb],2,mean)))
    }
 
}
dataDtDaReps = cbind(tmp1,tmp2)
dataDtDaReps$Induced = NULL
if (IntermediateReportFlag){write.csv(dataDtDaReps,'TRIPinhReps_DeltaTime_DeltaATC.csv')}

tmp1 = data.frame()
tmp2 = data.frame()
for (c in 1:length(uDrug)){
  tmp3 = dataDtDaReps[dataDtDaReps$Drug == uDrug[c],]
  tmp1 = rbind(tmp1,tmp3[1,1:3])
  tmp2 = rbind(tmp2,apply(tmp3[,TFIb],2,mean))
}
colnames(tmp2) = TFIb
avgDtDaAll = cbind(tmp1,tmp2)

# Normalize by OD Information
dataDtDaRepsOD = dataDtDaReps
for (i in TFIb){
  dataDtDaRepsOD[,i] = dataDtDaReps[,i]/dataDtDaRepsOD[,2]
}
if (IntermediateReportFlag){write.csv(dataDtDaRepsOD,'TRIPinh2016ODnormalized_DeltaTime_DeltaATC.csv')}

tmp1 = data.frame()
tmp2 = data.frame()
for (c in uDrug){
  tmp3 = dataDtDaRepsOD[dataDtDaRepsOD$Drug == c,]
  tmp1 = rbind(tmp1,tmp3[1,1:2])
  tmp2 = rbind(tmp2,apply(tmp3[,TFIb],2,mean))
}
colnames(tmp2) = TFIb
avgDtDaOD = cbind(tmp1,tmp2)

INHvUTsummary = data.frame(Untreated = t(avgDtDaOD[avgDtDaOD$Drug == 'Untreated',TFIb]),INH = t(avgDtDaOD[avgDtDaOD$Drug == 'INH',TFIb]))
colnames(INHvUTsummary) = c('Untreated','INH')
INHvUTsummary$INHdiffUT = INHvUTsummary$INH - INHvUTsummary$Untreated

INHvUTsummary$zUT = t(dataDeltaTimezscoreAvgExpAvg[dataDeltaTimezscoreAvgExpAvg$Drug == 'Untreated',TFIb])
INHvUTsummary$zINH1 = t(dataDeltaTimezscoreAvgExpAvg[dataDeltaTimezscoreAvgExpAvg$Drug == 'INH',TFIb])
INHvUTsummary$TpUT = t(dataDeltaTimepvalT[dataDeltaTimepvalT$Drug == 'Untreated',TFIb])
INHvUTsummary$TpINH1 = t(dataDeltaTimepvalT[dataDeltaTimepvalT$Drug == 'INH',TFIb])
write.csv(INHvUTsummary,'TRIPINH1vUntreatedAvg.csv')



