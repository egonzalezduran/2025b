
####################### R Code for statistics and modelling ################################################################################
## Running title: Control of mitochondrial inheritance                                                                                  ###  
## by Enrique Gonzalez-Duran, Zizhen Liang, Joachim Forner, Dennis Kleinschmidt, Weiqi Wang, Liwen Jiang, Kin Pan Chung* & Ralph Bock*  ###
## 2025                                                                                                                                 ###    
## Version 15.07.25 by Enrique Gonzalez-Duran                                                                                           ###                                                                                                  
## R version 4.3.3                                                                                                                      ###  
## Max Planck Institute of Molecular Plant Physiology, Potsdam-Golm, Germany                                                            ### 
##                                                                                                                                      ###
######################################################################################################################################  ###                         

##This code is set to work in R version 4.3.3.

##List of additional packages:
# ggplot2     v. 3.5.1
# ggsignif    v. 0.6.4  
# rstudioapi  v. 0.16.0
# MASS        v. 7.3-60.0.1
# dplyr       v. 1.1.4
# multcomp    v. 1.4-25

wants <- c("ggplot2","ggsignif","rstudioapi","MASS","dplyr","multcomp") ## searches for and installs the packages
has   <- wants %in% rownames(installed.packages())
if(any(!has)) install.packages(wants[!has])

library(ggplot2)
library(ggsignif)
library(rstudioapi)
library(MASS) 
library(dplyr)
library(multcomp)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # sets the location of the file as the working directory, works in R Studio
getwd()


####### Comparison of Number of mitochondria in GC#################################

mitoGC <-  as.data.frame(matrix(data=NA, nrow=10, ncol=2))
mitoGC[,1]<-c(0,4,0,1,1,0,1,0,1,2)
mitoGC[,2]<-c(8,6,2,5,3,6,5,3,2,5)
colnames(mitoGC) <- c("pollen at 25°C","pollen at 10°C")
row.names(mitoGC) <- c("TOMO1","TOMO2","TOMO3","TOMO4","TOMO5","TOMO6","TOMO7","TOMO8","TOMO9","TOMO10")
mitoGC

t.test(mitoGC["pollen at 10°C"],mitoGC["pollen at 25°C"],alternative = "two.sided",paired=FALSE,var.equal=FALSE,conf.level = 0.95)

#Welch Two Sample t-test

#data:  mitoGC["pollen at 10°C"] and mitoGC["pollen at 25°C"]
#t = 4.7678, df = 15.272, p-value = 0.0002374
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  1.937743 5.062257
#sample estimates:
#  mean of x mean of y 
#4.5       1.0 

dfGCplot <-  as.data.frame(matrix(data=NA, nrow=20, ncol=2))
dfGCplot[1:10,2]<-c(8,6,2,5,3,6,5,3,2,5)
dfGCplot[11:20,2]<-c(0,4,0,1,1,0,1,0,1,2)
dfGCplot[,1]<- c( rep("Pollen at 10°C",10),rep("Pollen at 25°C",10))
colnames(dfGCplot) <- c("Group","Value")
dfGCplot$Group= with(dfGCplot,reorder(Group,Value,median))
dfGCplot

plot.GC <- ggplot(dfGCplot,aes(x=Group,y=Value,fill=Group))+
            geom_boxplot(outlier.size = 0, alpha=0.7, linewidth=0.5, show.legend=FALSE)+
            scale_y_continuous(limits=c(0,10),breaks = c(0,1,2,3,4,5,6,7,8,9,10),labels=c(0," ",2," ",4," ",6," ",8," ",10))+
            geom_point(aes(fill=Group),alpha = 0.9, size = 4, shape=21, show.legend=FALSE, position = position_jitter(width = 0.3,height = 0,seed=10))+
            scale_fill_manual(values=c("#faaa44","#1F78B4"))+          
              xlab(" ") +
              ylab("Number of mitochondria in GC") +
              theme(panel.background = element_blank(), 
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  text=element_text(size=12,family="sans"),
                  axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
                  axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"))+
            geom_signif(comparisons=list(c("Pollen at 10°C", "Pollen at 25°C")), annotations="***",
                  y_position = 9.3, tip_length = 0.05, vjust=0.4) 

plot.GC

pdf(file="MitoGC_R_output.pdf", height = 4, width = 5)
plot.GC
dev.off()

####### Model for effect of cold and dpd1 mutation in mitochondrial inheritance################################

PT.df <-  as.data.frame(matrix(data=NA, nrow=4, ncol=6))
PT.df[1,] <- c("nad9xDK53gh","gh","dk53",2,1105,NA)
PT.df[2,] <- c("nad9xDK53c","cold","dk53",3,386,NA)
PT.df[3,] <- c("nad9xdpd1gh","gh","dpd1",1,375,NA)
PT.df[4,] <- c("nad9xdpd1c","cold","dpd1",27,368,NA)
PT.df[,4]<- as.numeric(PT.df[,4])
PT.df[,5]<- as.numeric(PT.df[,5])
colnames(PT.df) <- c("Group","Temperature","Genotype","Positives","Total.seedlings","PT")

PT.df$Genotype = factor(PT.df$Genotype,
                        levels=unique(PT.df$Genotype))
PT.df$Temperature = factor(PT.df$Temperature,
                           levels=unique(PT.df$Temperature))
PT.df$Group = factor(PT.df$Group,
                     levels=unique(PT.df$Group))
convec1<-as.vector(as.factor(unique(PT.df[,3])))
convec1
convec2<-as.vector(as.factor(unique(PT.df[,2])))
convec2
convec3<-as.vector(as.factor(unique(PT.df[,1])))
convec3

PT.df[,6] <- PT.df$Positives / PT.df$Total.seedlings

M1 <- glm(PT ~ Group, weights = Total.seedlings, family= binomial(link= "log"),
          contrast=list(Group=contr.treatment(convec3,base=1)),data=PT.df)
summary(M1)

hpsetM1<-cftest(glht(M1), c("Groupnad9xDK53c","Groupnad9xdpd1gh","Groupnad9xdpd1c"), test= Chisqtest())
hptestedcM1<-summary(hpsetM1, test=adjusted(type="holm"))
ciM1<-confint(glht(M1),c("Groupnad9xDK53c","Groupnad9xdpd1gh","Groupnad9xdpd1c"),level=0.95,calpha = adjusted_calpha())
hptestedcM1

ciM1
ciM1df<- as.data.frame(matrix(data= ciM1[["confint"]],ncol=3,nrow=4))
colnames(ciM1df) <- c("Estimate","Upper","Lower")
row.names(ciM1df) <- c("nad9xDK53gh","nad9xDK53c","nad9xdpd1gh","nad9xdpd1c")
ciM1df_freq <- ciM1df[1:4,1:3] + ciM1df[1,1]
ciM1df_freq[1,1:3] <- ciM1df_freq[1,1:3]- ciM1df[1,1]
ciM1df_freq_log2 <- ciM1df_freq / 0.69314718056
ciM1df_freq_numeric <- 2 ^ (ciM1df_freq_log2)
ciM1df_freq_log2[,4] <- c("WT 25°C","WT 10°C","dpd1 25°C","dpd1 10°C")
ciM1df_freq_log2[,5] <- c("","ns","ns","***") 
colnames(ciM1df_freq_log2) <- c("Frequency","Upper","Lower","Group","Significant")
ciM1df_freq_log2$Group = factor(ciM1df_freq_log2$Group,
                     levels=unique(ciM1df_freq_log2$Group))
ciM1df_freq_log2
ciM1df_freq_numeric 


plot.M1 <- ggplot(ciM1df_freq_log2, aes(x = Group,y= Frequency, fill=Group))+
  scale_fill_manual(values=c("#faaa44","#1F78B4","#3b1678","#fe6800"))+
  scale_colour_manual(values=c("#faaa44","#1F78B4","#3b1678","#fe6800"))+
  scale_y_continuous(limits=c(-13,0),breaks = c(0,-2,-4,-6,-8,-10,-12))+
  geom_errorbar(aes(ymin=Lower, ymax=Upper, colour = Group),width=0,linetype=1, show.legend= FALSE, size=1)+
  geom_point(aes(fill=Group),alpha = 1, size = 4, shape=21, show.legend=FALSE) +
  geom_text(aes(label=Significant, y=-0.5))+
  xlab(" ") +
  ylab("Freq of paternal mt transmission (Log2)") +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(size=12,family="sans"),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"))
plot.M1
                                                 
pdf(file="Mitoinheritance_R_output_M1.pdf", height = 4, width = 5)
plot.M1
dev.off()                                        

####### Corrected estimate of mitochondrial transmission to account for pooling strategy##########

## the following code simulates the number of positive pools (mean and visual distribution) that could be observed at a given transminad9ion rate 

set.seed(seed = 1)

#Parameters:

Prob <- c(0.07336957,0.087) #this vector can have multiple frequencies
#the first is the empirical probability (7.33%, the second frequency in most simulations produces on average the number of pools we observed empirically. This number was found through stepwise additions of 0.001 to the empirical rate. 
Nseeds <- 370 #approximate number of seeds in which we saw 27 nad9-positive samples
Nsim <- 10000 #number of simulations
Size.cluster <- 5 #size of the pool

------------------------------------
length(Prob) #number of probabilities being tested
Ncluster <- Nseeds / Size.cluster #number of pools 

mito.data <- data.frame(matrix(data=NA,nrow=(length(Prob)*Nseeds),ncol=((Nsim)+1)))
mito.data[,1] <- rep(Prob,each=Nseeds)   

#produces a table where the first column has the probability, each row is a seedling to be simulated, and each remaining column a simulation. All columns aside from the first are empty 


for (j in 1:(length(Prob))){
  for (i in 2:(Nsim+1)) {
    mito.data[((Nseeds*(j-1))+1):(Nseeds*j),i] <- rbinom(Nseeds, 1, Prob[j])
  }
} #produces Nseeds rows per simulation (column), where each seedling (rows) can be nad9-positive or not (1 or 0), at random according to the frequency in the Prob vector

mito.after.clus <- data.frame(matrix(data=NA,nrow=(length(Prob)*Nseeds)/Size.cluster,ncol=(ncol(mito.data)-1)))
for (i in 1:Nsim) {
  for (j in 1:(Ncluster*length(Prob))) {
    if ( sum(mito.data[(((j-1)*Size.cluster)+1):(j*Size.cluster),i+1])>= 1 ){ 
      mito.after.clus[j,i] <- 1 }
    else {
      mito.after.clus[j,i] <- 0 }
  }}

#the table is condensed: 1 pool (one row of the output table) consists of 5 seedlings (rows in the input table). A 1 is assigned if there was one or more seedlings in the pool positive.  

positives.per.simulation.mito<- data.frame(matrix(data=NA,nrow=(Nsim*length(Prob)),ncol= 1))
for (j in 1:(length(Prob))){
  for (i in 1:Nsim) {
    positives.per.simulation.mito[((Nsim*(j-1))+i),1] <-sum(mito.after.clus[((j-1)*Ncluster+1):(j*Ncluster),i])
  }}
df.plot <- data.frame(prob=factor(rep(Prob,each=Nsim)),
                      pos= positives.per.simulation.mito)
colnames(df.plot)<- c("prob","pos")

#simply sums the number of positive pools per simulation, per frequency given in prom

summary.mitoprobs <- df.plot %>% 
  group_by(prob) %>%
  summarize(mean= mean(pos), sd= sd(pos))

#summarizes the table above

SimPoolPlot<- ggplot(df.plot, aes(x=pos, color= prob, fill= prob)) + 
  geom_density(alpha=.2) +
  geom_histogram(alpha=.1, binwidth = 1) +
  geom_vline(data=summary.mitoprobs, aes(xintercept=mean,color=prob),linetype="dashed") 
SimPoolPlot

#makes an histogram of the nad9-positive pools (x-axis), vs number of simulations where that number of positives are found (y-axis), color-coded depending on the proposed frequency

pdf(file="PoolSimulation.pdf", height = 4, width = 5)
SimPoolPlot
dev.off()
summary.mitoprobs 
sink(file= "positive pools expected at given transmission frequency.txt")
summary.mitoprobs 
sink()

