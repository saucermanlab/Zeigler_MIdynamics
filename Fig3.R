## Generates PCA plots for Figure 3
# must run fig3.m script first to generate data for use in these plots
# last updated ACZ 2.18.2020


# install.packages('ggplot2')
# install.packages('tidyverse')
# install.packages('ggthemes')
# install.packages('ggrepel')

#can skip the installation step if these are already installed
library(ggplot2)
library(tidyverse)
library(readr)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(readxl)

##
# set working directory as folder where files are saved
# setwd()

#load matrix as DF
scores=data.frame(read_csv('PCAscores.csv', col_names = T))  #load in the csv you exported in matlab, make sure the path and filename are the same for the exported file and this file
colnames(scores) <- c('scores1','scores2','scores3','labels','highlight','NodeMatch')

loads <- data.frame(read_csv('PCAload.csv', col_names = T))
colnames(loads) <- c('coeff1','coeff2','coeff3','speciesNames')

#ggplot2 help here: https://ggplot2.tidyverse.org/



#plot of scores per stimulus
#must update labels with actual %explained
ggplot(scores)+geom_point(aes(x=scores1,y=scores2,color=as.factor(highlight),alpha=0.8))+geom_text_repel(data=scores[scores$NodeMatch==1,],aes(size=25,label=labels,x=scores1,y=scores2,color=as.factor(highlight)))+
  labs(title='PCA of Stimuli and Timepoints')+ theme_classic()+scale_colour_manual(values = c("#e31a1c","#7bccc4","#3f007d","#1d91c0","#253494","black"))+
  xlim(-1.5,2.5)+ylim(-1,1)+theme(legend.position="none")+
  xlab("PC1 (81% Variance Explained)")+
  ylab("PC2 (7.8% Variance Explained)")+
  geom_hline(yintercept=0, linetype="dashed",color='#bdbdbd')+geom_vline(xintercept=0, linetype="dashed",color='#bdbdbd')

#biplot of loadings per node
ggplot(loads)+geom_point(aes(x=coeff1,y=coeff2))+geom_text_repel(aes(label=speciesNames,x=coeff1,y=coeff2))+
  xlab("PC1 (81% Variance Explained)")+
  ylab("PC2 (7.8% Variance Explained)")+
  labs(title='Loadings of Nodes')+
  theme_classic()+xlim(-1.5,2.5)+ylim(-1,1)+geom_hline(yintercept=0, linetype="dashed",color='#bdbdbd')+geom_vline(xintercept=0, linetype="dashed",color='#bdbdbd')

## showing pc1-3
#plot of scores per stimulus
#must update geom_text_repel with the indices of actual matches
ggplot(scores)+geom_point(aes(x=scores1,y=scores3,color=as.factor(highlight),alpha=0.8))+geom_text_repel(data=scores[scores$NodeMatch == 1,],aes(size=25,label=labels,x=scores1,y=scores2,color=as.factor(highlight)))+
  labs(title='PCA of Stimuli and Timepoints')+ theme_classic()+scale_colour_manual(values = c("#e31a1c","#7bccc4","#3f007d","#1d91c0","#253494","black"))+
  xlim(-1.5,2.5)+ylim(-1,1)+theme(legend.position="none")+
  xlab("PC1 (81% Variance Explained)")+
  ylab("PC3 (6.8% Variance Explained)")+
  geom_hline(yintercept=0, linetype="dashed",color='#bdbdbd')+geom_vline(xintercept=0, linetype="dashed",color='#bdbdbd')

#biplot of loadings per node
ggplot(loads)+geom_point(aes(x=coeff1,y=coeff3))+geom_text_repel(aes(label=speciesNames,x=coeff1,y=coeff3))+
  xlab("PC1 (81% Variance Explained)")+
  ylab("PC3 (6.8% Variance Explained)")+
  labs(title='Loadings of Nodes')+
  theme_classic()+xlim(-1.5,2.5)+ylim(-1,1)+geom_hline(yintercept=0, linetype="dashed",color='#bdbdbd')+geom_vline(xintercept=0, linetype="dashed",color='#bdbdbd')

## showing PC 2-3
ggplot(scores)+geom_point(aes(x=scores2,y=scores3,color=as.factor(highlight),alpha=0.8))+geom_text_repel(data=scores[scores$NodeMatch == 1,],aes(size=25,label=labels,x=scores1,y=scores2,color=as.factor(highlight)))+
  labs(title='PCA of Stimuli and Timepoints')+ theme_classic()+scale_colour_manual(values = c("#e31a1c","#7bccc4","#3f007d","#1d91c0","#253494","black"))+
  xlim(-1.7,2.5)+ylim(-1,1)+theme(legend.position="none")+
  xlab("PC2 (7.8% Variance Explained)")+
  ylab("PC3 (6.8% Variance Explained)")+
  geom_hline(yintercept=0, linetype="dashed",color='#bdbdbd')+geom_vline(xintercept=0, linetype="dashed",color='#bdbdbd')

ggplot(loads)+geom_point(aes(x=coeff2,y=coeff3))+geom_text_repel(aes(label=speciesNames,x=coeff2,y=coeff3))+
  xlab("PC2 (7.8% Variance Explained)")+
  ylab("PC3 (6.8% Variance Explained)")+
  labs(title='Loadings of Nodes')+
  theme_classic()+xlim(-1.7,2.5)+ylim(-1,1)+geom_hline(yintercept=0, linetype="dashed",color='#bdbdbd')+geom_vline(xintercept=0, linetype="dashed",color='#bdbdbd')

