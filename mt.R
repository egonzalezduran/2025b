
####################### R Code for statistics and modelling ###############################################################################
## High-frequency biparental inheritance of plant mitochondria upon chilling stress and loss of a genome-degrading nuclease                                                                                                                                      ###  
## by Enrique Gonzalez-Duran, Zizhen Liang, Joachim Forner, Dennis Kleinschmidt, Weiqi Wang, Liwen Jiang, Kin Pan Chung* & Ralph Bock*  ###
## 2025                                                                                                                                 ###    
## Version 14.11.25 by Enrique Gonzalez-Duran                                                                                           ###                                                                                                  
## R version 4.3.3                                                                                                                      ###  
## Max Planck Institute of Molecular Plant Physiology, Potsdam-Golm, Germany                                                            ###                                                                                                                                       ###
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
getwd() # shows the working directory


####### A. Comparison of Number of mitochondria in GC at different temperatures#################################

#t-test
mitoGC<- as.data.frame(read.table("data_for_GC_comparison_ttest.txt", header = TRUE)) #passes data from .txt file to a data frame
mitoGC_ttest<- t.test(mitoGC["pollen_at_10C"],mitoGC["pollen_at_25C"],alternative = "two.sided",paired=FALSE,var.equal=FALSE,conf.level = 0.95) #performs t.test
mitoGC_ttest #shows result

sink(file= "GC_comparison_ttest_result.txt") #sinks result to a .txt file
mitoGC_ttest
sink()

#plot of the comparison
dfGCplot<- as.data.frame(read.table("data_for_GC_comparison_plot.txt", header = TRUE)) #passes data table from .txt file to a data frame
dfGCplot
plot.GC <- ggplot(dfGCplot,aes(x=Group,y=Value,fill=Group))+
  geom_boxplot(outlier.size = 0, alpha=0.7, linewidth=0.5, show.legend=FALSE)+
  scale_y_continuous(limits=c(0,10),breaks = c(0,1,2,3,4,5,6,7,8,9,10),labels=c(0," ",2," ",4," ",6," ",8," ",10))+
  geom_point(aes(fill=Group),alpha = 0.9, size = 4, shape=21, show.legend=FALSE, position = position_jitter(width = 0.3,height = 0,seed=10))+
  scale_fill_manual(values=c("#1F78B4","#faaa44"))+          
  xlab(" ") +
  ylab("Number of mitochondria in GC") +
  theme(panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(size=12,family="sans"),
        axis.line.x = element_line(linewidth = 0.5, linetype = "solid", colour = "black"),
        axis.line.y = element_line(linewidth = 0.5, linetype = "solid", colour = "black"))+
  geom_signif(comparisons=list(c("Pollen_at_10C", "Pollen_at_25C")), annotations="***",
              y_position = 9.3, tip_length = 0.05, vjust=0.4) 

plot.GC
pdf(file="MitoGC_plot.pdf", height = 4, width = 5)
plot.GC
dev.off()

####### B. Model for effect of cold and dpd mutation in mitochondrial inheritance################################

PT.df <- as.data.frame(read.table("RT_screen_data.txt", header = TRUE))
PT.df[,6] <- PT.df$Positives / PT.df$Total.seedlings #Fills column PT with the transmission rate (positives/total seedlings screen)
PT.df
PT.df$Group = factor(PT.df$Group,
                     levels=unique(PT.df$Group))
compvector<-as.vector(as.factor(unique(PT.df[,1])))
compvector

M1 <- glm(PT ~ Group, weights = Total.seedlings, family= binomial(link= "log"),
          contrast=list(Group=contr.treatment(compvector,base=1)),data=PT.df)
summary(M1) #Calculates the linear model for the statistics

hpsetM1<-cftest(glht(M1), c("Groupnad9xDK53c","Groupnad9xdpd1gh","Groupnad9xdpd1c"), test= Chisqtest()) #Sets the groups of hypothesis to be tested based on the M1 output. Changing the names of groups
hptestedcM1<-summary(hpsetM1, test=adjusted(type="holm")) #calculates differences between groups based on Chisq statistics, ajusted for multiple comparisons through Holm method
hptestedcM1 #show test results

sink(file= "Hypothesis_testing_results_screen.txt") ### prints the result of the hypothesis testing to file
hptestedcM1
sink()

ciM1<-confint(glht(M1),c("Groupnad9xDK53c","Groupnad9xdpd1gh","Groupnad9xdpd1c"),level=0.95,calpha = adjusted_calpha()) #calculates confidence intervals
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

ciM1df_freq_log2 #values for the plot 
ciM1df_freq_numeric #unused values, same as previous data but in numeric instead of log2

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

pdf(file="RTscreen_results.pdf", height = 4, width = 5) #Prints the results of the test of inheritance
plot.M1
dev.off()                                        

####### C. Simulation to produce a corrected estimate of mitochondrial transmission accounting for pooling strategy##########

## the following code simulates the number of positive pools (mean and visual distribution) that could be observed at a given transmission rate 

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

#simply sums the number of positive pools per simulation, per frequency given in Prob

summary.mitoprobs <- df.plot %>% 
  group_by(prob) %>%
  summarize(mean= mean(pos), sd= sd(pos))

#summarizes the table above

SimPoolPlot<- ggplot(df.plot, aes(x=pos, color= prob, fill= prob)) + 
  geom_density(alpha=.2) +
  geom_vline(data=summary.mitoprobs, aes(xintercept=mean,color=prob),linetype="dashed") 
SimPoolPlot
summary.mitoprobs 



pdf(file="PoolSimulation.pdf", height = 4, width = 5)#makes a density plot of the nad9-positive pools (x-axis), vs density of simulations where that number of positives are found (y-axis), color-coded depending on the proposed frequency
SimPoolPlot
dev.off()

sink(file= "positive pools expected at given transmission frequency.txt")
summary.mitoprobs 
sink()




