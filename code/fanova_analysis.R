################### Load Libraries ###################

#libraries used for data wrangling & exploratory analysis 
library(data.table)
library(tidyverse)
library(lubridate)
library(dplyr)

#libraries for visualizations
library(ggplot2)
library(plotly)

#libraries for functional data analysis
library(fda)
library(fda.usc)

#libraries for spatial estimation 
library(RandomFields)


################### Read Data Sets ###################

#set working directory to current directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

#read in SIF data. We have already subsetted the data and labeled it according to burned or unburned in exploratory.R
sifdat<-fread(file="../data/mendocino_sif.csv") #fread is a fast function to read in data from data.table package 
dim(sifdat)
head(sifdat)

#remove NA
sifdat<-sifdat[-which(is.na(sifdat$sif)),]

################### Plot Raw Data ###################

#generic data frame 
df<-sifdat

#fire day
fireday<-paste("2018-07-27")

mytheme=theme_light()+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
mytheme2=theme_classic()+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))

#Raw SIF (original sif time series ) with dashed line representing fire
Raw_sif_plot<-ggplot(data=df, aes(x=date, y=sif, group=lonlat, color = group), alpha=0.7)+
  geom_line()+labs(title="Raw SIF Time Series", y="SIF", x="Time", color="Group")+scale_color_manual(labels = c("Burned", "Unburned"), values=c('black','lightgreen'))+
  ylim(-1,1.5)+theme_classic()+mytheme+scale_x_date(breaks = "4 months")+
  geom_vline(xintercept = as.Date(fireday), linetype = "dashed")

#Raw SIF (original sif time series ) 
Raw_sif_plot<-ggplot(data=df, aes(x=date, y=sif, group=lonlat, color = group), alpha=0.7)+
  geom_line()+labs(title="Raw SIF Time Series", y="SIF", x="Time", color="Group")+scale_color_manual(labels = c("Burned", "Unburned"), values=c('black','lightgreen'))+
  ylim(-1,1.5)+theme_classic()+mytheme+scale_x_date(breaks = "4 months")

#Raw SIF with day axis
Raw_sif_plot<-ggplot(data=df, aes(x=dayssinceapril1, y=sif, group=lonlat, color = group), alpha=0.7)+
  geom_line()+labs(title="Raw SIF Time Series", y="SIF", x="Days Since April 1, 2018", color="Group")+scale_color_manual(labels = c("Burned", "Unburned"), values=c('black','lightgreen'))+
  ylim(-1,1.5)+theme_classic()+mytheme


################### Obtain Smooth Functions ###################

#subset unburned group
df_u<-subset(df, group=="unburned")
df_b<-subset(df, group=="burned")


#Create list with individual time series 
dfU_curves<-NULL
dfB_curves<-NULL
for (i in 1:length(unique(df_u$lonlat))){
  dfU_curves[[i]]<-subset(df_u, lonlat==unique(df_u$lonlat)[i])
}

for (i in 1:length(unique(df_b$lonlat))){
  dfB_curves[[i]]<-subset(df_b, lonlat==unique(df_b$lonlat)[i])
}

length(dfU_curves)+length(dfB_curves)


#Find a common domain to smooth over. Largest range at which all the time series overal
range<-matrix(NA,nrow=length(dfU_curves)+length(dfB_curves),ncol=2)

for (i in 1:length(dfU_curves)){
  x<-dfU_curves[[i]]$dayssinceapril1
  range[i,]<-range(x, na.rm = TRUE)
}

for (i in 1:length(dfB_curves)){
  x<-dfB_curves[[i]]$dayssinceapril1
  range[i+length(dfU_curves),]<-range(x, na.rm = TRUE)
}

start<-max(range[,1]) #first time point
end<-min(range[,2]) #last time point
c(start,end)

#create new xvalues to discretize your smooth functions to
argvals<-seq(from=start, to=end, length.out=50) 
#create matrix to store rediscretized functions
fullseriesU<-matrix(NA, nrow=length(dfU_curves), ncol=length(argvals)) 
fullseriesB<-matrix(NA, nrow=length(dfB_curves), ncol=length(argvals))

#Find the optimal smoothing parameter 
#combine all time series into one list 
dfL<-c(dfU_curves,dfB_curves)

#search over all grid of lambda to minimize gcv
lambda<-seq(from=-10, to=20, by=1)
gcv<-matrix(nrow=length(lambda), ncol=length(dfL))
for (i in 1:length(lambda)){
  for (j in 1:length(dfL)){
    new_time<-dfL[[j]]$dayssinceapril1 #time
    new_series<-dfL[[j]]$sif #series 
    wt<-1/dfL[[j]]$sif_sigma^2
    basis<-create.bspline.basis(rangeval = c(start,end),nbasis=20)
    fdpar = fdPar(basis, Lfdobj=int2Lfd(2), lambda=1*10^lambda[i])
    smooth<-smooth.basis(argvals=new_time, y = new_series, wtvec=wt, fdParobj=fdpar)
    gcv[i,j]<-smooth$gcv[[1]]
  }
  print(i)
}


head(gcv)
meangcv<-rowMeans(gcv)
gcvdf<-data.frame(lambda=lambda, meangcv=meangcv)
plot(gcvdf$lambda, gcvdf$meangcv, type="l")
optlambda<-1*10^gcvdf$lambda[which.min(gcvdf$meangcv)] #1000


# smooth each group using the optimal smoothing parameter lambda 
optlambda<-1*10^4
for (i in 1:nrow(fullseriesU)){
  new_time<-dfU_curves[[i]]$dayssinceapril1 #time
  new_series<-dfU_curves[[i]]$sif #series 
  wt<-1/dfU_curves[[i]]$sif_sigma^2
  basis<-create.bspline.basis(rangeval = c(start,end),nbasis=20)
  fdpar = fdPar(basis, Lfdobj=int2Lfd(2), lambda=optlambda)
  smooth<-smooth.basis(argvals=new_time, y = new_series, wtvec=wt, fdParobj=fdpar)
  fullseriesU[i,]<-fdata(smooth$fd, argvals=argvals)$data #save the re-discretized series post-smoothing 
}

for (i in 1:nrow(fullseriesB)){
  new_time<-dfB_curves[[i]]$dayssinceapril1 #time
  new_series<-dfB_curves[[i]]$sif #series 
  wt<-1/dfB_curves[[i]]$sif_sigma^2
  basis<-create.bspline.basis(rangeval = c(start,end),nbasis=20)
  fdpar = fdPar(basis, Lfdobj=int2Lfd(2), lambda=optlambda)
  smooth<-smooth.basis(argvals=new_time, y = new_series, wtvec=wt, fdParobj=fdpar)
  fullseriesB[i,]<-fdata(smooth$fd, argvals=argvals)$data #save the re-discretized series post-smoothing 
}

#create combined long for data set for plotting
smoothdfU<-data.frame(argvals=argvals, smoothU = t(fullseriesU), meanUnburned=rowMeans(t(fullseriesU)))
smoothdfU_l<-smoothdfU %>% pivot_longer(-argvals, names_to = "variable")
smoothdfU_l$group<-"unburned"

smoothdfB<-data.frame(argvals=argvals, smoothB = t(fullseriesB), meanBurned=rowMeans(t(fullseriesB)))
smoothdfB_l<-smoothdfB %>% pivot_longer(-argvals, names_to = "variable")
smoothdfB_l$group<-"burned"
smoothdf<-rbind(smoothdfU_l, smoothdfB_l)
head(smoothdf)

#plot smooth curves along with smooth mean curves of each group
smoothplot1<-ggplot(data=smoothdf, aes(x=argvals, y=value, group=variable, col=group))+geom_line()+scale_color_manual(labels = c("burned", "unburned"), values=c('darkgrey','lightgreen'))+
  geom_line(data=subset(smoothdf, variable == "meanUnburned" | variable == "meanBurned"), aes(x=argvals, y=value, lty = variable), col="black")+
  labs(title="Smooth SIF Curves", x="Days Since April 1, 2018", y="SIF", col="Group", lty="Mean")+mytheme
ggplotly(smoothplot1)

#smooth curves with no means
smooth_plot2<-ggplot(data=smoothdf, aes(x=argvals, y=value, group=variable, col=group))+geom_line()+scale_color_manual(labels = c("burned", "unburned"), values=c('black','lightgreen'))+
  labs(title="Smooth Curves", x="Days Since April 1, 2018", y="SIF", col="Group", lty="Mean")+mytheme2+ylim(-0.5,1.5)

#smooth means plot
meanplot<-subset(smoothdf, variable == "meanUnburned" | variable == "meanBurned")
meanplot$variable[meanplot$variable=="meanUnburned"]<-"unburned"
meanplot$variable[meanplot$variable=="meanBurned"]<-"burned"

smooth_mean_plot<-ggplot()+
  geom_line(data=meanplot, aes(x=argvals, y=value, lty = variable), col="black")+
  geom_line()+labs(title="Mean Smooth Curves", x="Days Since April 1, 2018", y="SIF", lty="Mean")+mytheme2+ylim(-0.5,1.5)



#################### Functional ANOVA ########################
#accounts for spatial and temporal dependence 

#represent the smooth curves as a functional object
basis<-create.bspline.basis(rangeval = range(argvals),nbasis=20)
fdpar = fdPar(basis, Lfdobj=int2Lfd(2))

smooth1<-smooth.basis(argvals = argvals, y=t(fullseriesU), fdParobj=fdpar)
smooth2<-smooth.basis(argvals = argvals, y=t(fullseriesB), fdParobj=fdpar)
plot(smooth2) #plot smooth objects

#calculate test statistic
smean1<-mean.fd(smooth1$fd)
smean2<-mean.fd(smooth2$fd)
stat<-metric.lp(fdata(smean1, argvals=argvals), fdata(smean2, argvals=argvals))[1] #1.97
stat #7.056177 

# represent each curves by fpca (functional principal component analysis) 
pcalist1 = pca.fd(smooth1$fd, 4, centerfns=TRUE)
pcalist2 = pca.fd(smooth2$fd, 4, centerfns=TRUE)
#plot.pca.fd(pcalist1)
pcalist1$values # vector length 20 #eigenvalues 
pcalist1$varprop # 2 PCs represents 0.9914916% of variability  
pcalist1$scores #5 x 2 #scores 
sum(pcalist1$varprop) #0.97
sum(pcalist2$varprop) #1

dim(pcalist1$harmonics$coefs) #20 x 4 harmonics or eigenfunctions
pcmean1<-pcalist1$meanfd 
pcmean2<-pcalist2$meanfd

scores1 = pcalist1$scores
PC1 = pcalist1$harmonics

scores2 = pcalist2$scores
PC2 = pcalist2$harmonics


#reconstruct the curves from the fpca explansions
fullseries1<-fullseriesU
fullseries2<-fullseriesB
samplepcres1<-matrix(NA, ncol=nrow(fullseries1), nrow=50)
samplepcres2<-matrix(NA, ncol=nrow(fullseries2), nrow=50)
for (i in 1:nrow(fullseries1)){
  samplepcres1[,i]<-fdata(pcmean1+scores1[i,1]*PC1[1]+ scores1[i,2]*PC1[2] + scores1[i,3]*PC1[3] +scores1[i,4]*PC1[4], argvals=argvals)$data
}
for (i in 1:nrow(fullseries2)){
  samplepcres2[,i]<-fdata(pcmean2+scores2[i,1]*PC2[1]+ scores2[i,2]*PC2[2] + scores2[i,3]*PC2[3] + scores2[i,4]*PC2[4], argvals=argvals)$data
}


#stat
stat<-metric.lp(fdata(mean(fdata(t(samplepcres1)))$data, argvals=argvals), 
                fdata(mean(fdata(t(samplepcres2)))$data, argvals=argvals))[1] 

stat #7.95 #double check test stat (should be the same after "reconstruction" from fpca)


#calculate the sample difference between mean curves 
sample_diff<-c(mean(fdata(t(samplepcres2)))$data-mean(fdata(t(samplepcres1)))$data) #Burned - Unburned 


#obtain the sample locations 
loc1<-as.matrix(data.frame(x_coord=rep(NA,nrow(fullseries1)),y_coord=rep(NA,nrow(fullseries1))))
loc2<-as.matrix(data.frame(x_coord=rep(NA,nrow(fullseries2)),y_coord=rep(NA,nrow(fullseries2))))
for (i in 1:nrow(fullseries1)){
  loc1[i,1]<-unique(dfU_curves[[i]]$lon)
  loc1[i,2]<-unique(dfU_curves[[i]]$lat)
}
for (i in 1:nrow(fullseries2)){
  loc2[i,1]<-unique(dfB_curves[[i]]$lon)
  loc2[i,2]<-unique(dfB_curves[[i]]$lat)
}

## put together estimated scores and its locations
dta1_PC1 <- cbind(loc1, scores=scores1[,1])
dta1_PC2 <- cbind(loc1, scores=scores1[,2])
dta1_PC3 <- cbind(loc1, scores=scores1[,3])
dta1_PC4 <- cbind(loc1, scores=scores1[,4])

dta2_PC1 <- cbind(loc2, scores=scores2[,1])
dta2_PC2 <- cbind(loc2, scores=scores2[,2])
dta2_PC3 <- cbind(loc2, scores=scores2[,3])
dta2_PC4 <- cbind(loc2, scores=scores2[,4])

#visualize the scores (spatial dependence is retained in the scores
score1_1<-ggplot(data.frame(dta1_PC1), aes(x=x_coord, y=y_coord, col=scores))+geom_point()+geom_point(size=10)+scale_color_gradient(low="yellow", high="red")
score1_2<-ggplot(data.frame(dta1_PC2), aes(x=x_coord, y=y_coord, col=scores))+geom_point()+geom_point(size=10)+scale_color_gradient(low="yellow", high="red")
score1_3<-ggplot(data.frame(dta1_PC3), aes(x=x_coord, y=y_coord, col=scores))+geom_point()+geom_point(size=10)+scale_color_gradient(low="yellow", high="red")
score1_4<-ggplot(data.frame(dta1_PC4), aes(x=x_coord, y=y_coord, col=scores))+geom_point()+geom_point(size=10)+scale_color_gradient(low="yellow", high="red")

score2_1<-ggplot(data.frame(dta2_PC1), aes(x=x_coord, y=y_coord, col=scores))+geom_point()+geom_point(size=10)+scale_color_gradient(low="yellow", high="red")
score2_2<-ggplot(data.frame(dta2_PC2), aes(x=x_coord, y=y_coord, col=scores))+geom_point()+geom_point(size=10)+scale_color_gradient(low="yellow", high="red")
score2_3<-ggplot(data.frame(dta2_PC3), aes(x=x_coord, y=y_coord, col=scores))+geom_point()+geom_point(size=10)+scale_color_gradient(low="yellow", high="red")
score2_4<-ggplot(data.frame(dta2_PC4), aes(x=x_coord, y=y_coord, col=scores))+geom_point()+geom_point(size=10)+scale_color_gradient(low="yellow", high="red")


## estimate the spatial parameters of the scores 
estmodel <- RMexp(var=NA, scale=NA)
fit <- RFfit(estmodel, data=dta1_PC1)
v1_1<-fit@ml@globalvariance #var
s1_1<-fit@ml@param[1] #scale 

estmodel <- RMexp(var=NA, scale=NA)
fit <- RFfit(estmodel, data=dta1_PC2)
v1_2<-fit@ml@globalvariance #var
s1_2<-fit@ml@param[1] #scale

estmodel <- RMexp(var=NA, scale=NA)
fit <- RFfit(estmodel, data=dta1_PC3)
v1_3<-fit@ml@globalvariance #var
s1_3<-fit@ml@param[1] #scale

estmodel <- RMexp(var=NA, scale=NA)
fit <- RFfit(estmodel, data=dta1_PC4)
v1_4<-fit@ml@globalvariance #var
s1_4<-fit@ml@param[1] #scale

estmodel <- RMexp(var=NA, scale=NA)
fit <- RFfit(estmodel, data=dta2_PC1)
v2_1<-fit@ml@globalvariance #var
s2_1<-fit@ml@param[1] #scale 

estmodel <- RMexp(var=NA, scale=NA)
fit <- RFfit(estmodel, data=dta2_PC2)
v2_2<-fit@ml@globalvariance #var
s2_2<-fit@ml@param[1] #scale

estmodel <- RMexp(var=NA, scale=NA)
fit <- RFfit(estmodel, data=dta2_PC3)
v2_3<-fit@ml@globalvariance #var
s2_3<-fit@ml@param[1] #scale 

estmodel <- RMexp(var=NA, scale=NA)
fit <- RFfit(estmodel, data=dta2_PC4)
v2_4<-fit@ml@globalvariance #var
s2_4<-fit@ml@param[1] #scale

# number of bootstrap iterations (number of resample sets)
N<-2000

# empty matrices and vectors for storing 
dep_pcres_mean1<-matrix(NA, ncol=N, nrow=50)
dep_pcres_mean2<-matrix(NA, ncol=N, nrow=50)
dist<-matrix(NA, N)
argvals_tilde<-rep(NA, 50)

# spatial locations same as sample data 
x_coord1<-loc1[,1]
y_coord1<-loc1[,2]
x_coord2<-loc2[,1]
y_coord2<-loc2[,2]

#seed for reproduceability 
set.seed(30)

#functions to reconstruct curves from fpca components
fun_reconstruct1 <- function(x) {
  fdata(x[1]*PC1[1]+ x[2]*PC1[2] + x[3]*PC1[3] + x[4]*PC1[4], argvals=argvals)$data
}
fun_reconstruct2 <- function(x) {
  fdata(x[1]*PC2[1]+ x[2]*PC2[2] + x[3]*PC2[3] + x[4]*PC2[4], argvals=argvals)$data
}

# Generate resamples 
ptm_start <- proc.time()
for (j in 1:N){
  #create spatial model with desired parameter values
  scoremodel1_1 <- RMtrend(mean=0)+RMexp(var=v1_1, scale = s1_1)
  scoremodel1_2 <- RMtrend(mean=0)+RMexp(var=v1_2, scale = s1_2)
  scoremodel1_3 <- RMtrend(mean=0)+RMexp(var=v1_3, scale = s1_3)
  scoremodel1_4 <- RMtrend(mean=0)+RMexp(var=v1_4, scale = s1_4)
  
  scoremodel2_1 <- RMtrend(mean=0)+RMexp(var=v2_1, scale = s2_1)
  scoremodel2_2 <- RMtrend(mean=0)+RMexp(var=v2_2, scale = s2_2)
  scoremodel2_3 <- RMtrend(mean=0)+RMexp(var=v2_3, scale = s2_3)
  scoremodel2_4 <- RMtrend(mean=0)+RMexp(var=v2_4, scale = s2_4)
  
  #simulate spatial field of scores for 4 principal components for each group
  z1_1 <- RFsimulate(scoremodel1_1, x=x_coord1, y=y_coord1)$variable1
  z1_2 <- RFsimulate(scoremodel1_2, x=x_coord1, y=y_coord1)$variable1
  z1_3 <- RFsimulate(scoremodel1_3, x=x_coord1, y=y_coord1)$variable1
  z1_4 <- RFsimulate(scoremodel1_4, x=x_coord1, y=y_coord1)$variable1
  
  z2_1 <- RFsimulate(scoremodel2_1, x=x_coord2, y=y_coord2)$variable1
  z2_2 <- RFsimulate(scoremodel2_2, x=x_coord2, y=y_coord2)$variable1
  z2_3 <- RFsimulate(scoremodel2_3, x=x_coord2, y=y_coord2)$variable1
  z2_4 <- RFsimulate(scoremodel2_4, x=x_coord2, y=y_coord2)$variable1
  
  #data of generated scores
  scores_tilde1<-cbind(z1_1,z1_2,z1_3,z1_4)
  scores_tilde2<-cbind(z2_1,z2_2,z2_3,z2_4)
  
  #construct resample curves from the fpcs
  pcres1<-apply(scores_tilde1, 1,FUN=fun_reconstruct1)
  pcres2<-apply(scores_tilde2, 1,FUN=fun_reconstruct2)
  
  #obtain mean resample curve for each group and save
  dep_pcres_mean1[,j]<-mean(fdata(t(pcres1)))$data
  dep_pcres_mean2[,j]<-mean(fdata(t(pcres2)))$data
  
  #obtain and save monte carlo statistic (distance between mean resample curve for each group)
  dist[j] <-metric.lp(fdata(mean(fdata(t(pcres1)))$data, argvals=argvals), 
                      fdata(mean(fdata(t(pcres2)))$data, argvals=argvals))[1] 
  
  #print itteration to keep track of progress
  print(j)
}

ptm <- proc.time()-ptm_start #time to generate resamples (typically a few mins)

#calculate p-values 
pval1<-length(dist[dist>stat])/length(dist) 
pval1 #0

#create visualizations 
resdiff<--dep_pcres_mean1+dep_pcres_mean2
resdf<-data.frame(argvals=argvals, res=resdiff, sample=sample_diff)
resdflong<-resdf %>% pivot_longer(-argvals, names_to = "variable")
mytheme=theme_light()+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
mytheme2=theme_classic()+theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))

#fanova visualization with p-value annotated
fanova1<-ggplot()+geom_line(data=resdflong, aes(x=as_date(argvals, origin = "2018-04-01"), y=value, group=variable), colour="darkgray", alpha = 0.2)+
  geom_line(data=subset(resdflong, variable == "res.1"|variable == "sample"),aes(x=as_date(argvals, origin = "2018-04-01"),y=value, col=variable))+
  scale_color_manual(labels = c("Resamples", "Sample"), values=c("darkgray",'black'))+
  theme(legend.title = element_blank())+geom_vline(xintercept = as.Date(fireday), linetype = "dashed")+
  labs(title="fANOVA Visualization (Burned - Unburned)",x= "time", y = "Mean Difference", col="")+mytheme2+scale_x_date(date_breaks="2 month")+
  annotate(geom="text", x=as_date(argvals[1])+Inf, y=Inf,label=paste("global p-value =",pval1),hjust = 1, vjust = 1)+ylim(-0.5,0.5)

#fanova visualization no p-value annotated
fanova2<-ggplot()+geom_line(data=resdflong, aes(x=as_date(argvals, origin = "2018-04-01"), y=value, group=variable), colour="darkgray", alpha = 0.2)+
  geom_line(data=subset(resdflong, variable == "res.1"|variable == "sample"),aes(x=as_date(argvals, origin = "2018-04-01"),y=value, col=variable))+
  scale_color_manual(labels = c("Resamples", "Sample"), values=c("darkgray",'black'))+
  theme(legend.title = element_blank())+geom_vline(xintercept = as.Date(fireday), linetype = "dashed")+
  labs(title="fANOVA Visualization (Burned - Unburned)",x= "time", y = "Mean Difference", col="")+mytheme2+scale_x_date(date_breaks="2 month")+
  ylim(-0.5,0.5)

#fanova visualization day axis
fanova3<-ggplot()+geom_line(data=resdflong, aes(x=argvals, y=value, group=variable), colour="darkgray", alpha = 0.2)+
  geom_line(data=subset(resdflong, variable == "res.1"|variable == "sample"),aes(x=argvals,y=value, col=variable))+
  scale_color_manual(labels = c("Resamples", "Sample"), values=c("darkgray",'black'))+
  theme(legend.title = element_blank())+
  labs(title="fANOVA Visualization (Burned - Unburned)",x= "Days since April 1, 2018", y = "Mean Difference", col="")+mytheme2+
  annotate(geom="text", x=Inf, y=Inf,label=paste("global p-value =",pval1),hjust = 1, vjust = 1)+ylim(-0.5,0.5)

#sample difference curve only 
diff_curve<-ggplot()+
  geom_line(data=subset(resdflong, variable == "sample"),aes(x=argvals,y=value, col=variable))+
  scale_color_manual(labels = c("Sample"), values=c('black'))+
  theme(legend.title = element_blank())+
  labs(title="Sample Difference (Burned - Unburned)",x= "Days since April 1, 2018", y = "Mean Difference", col="")+mytheme2+
  ylim(-0.5,0.5)

# create histogram of bootstrap distances to visualize p-value calculation
d<-data.frame(dist=dist)
head(d)

# Montecarlo statistic distribution with test statistic (dashed line)
hist1<-ggplot(data=d, aes(x=dist))+geom_histogram(col="black",fill="white")+mytheme2+geom_vline(xintercept = stat, linetype=2)+xlim(0,8.5)+
  annotate(geom="text", x=Inf, y=Inf,label=paste("global p-value =",pval1),hjust = 1, vjust = 1)+labs(y="Count", x="Distance", title = "Bootstrap Distances")

#test statistic only 
hist2<-ggplot(data=d, aes(x=dist))+geom_histogram(col="white",fill="white")+mytheme2+geom_vline(xintercept = stat, linetype=2)+xlim(0,8.5)+
  annotate(geom="text", x=Inf, y=Inf,label=paste("Test Statistic"),hjust = 1, vjust = 1)+labs(y="Count", x="Distance", title = "Bootstrap Distances")
