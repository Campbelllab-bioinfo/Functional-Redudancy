############################################################################################
# Computation and prediction of Microbial functional redundancy
############################################################################################

## required R libraries
# general libraries
library(reshape2)
library(ggplot2)
library(patchwork)
library(stringr)
library(ggpubr)
library(tidyverse)
library(openxlsx)

# stats and FD libraries
library(ade4)
library(vegan)
library(adiv)
library(parallelDist)

# location of the script
rstudioapi::getSourceEditorContext()$path

# create a theme for the figures
theme_depth<-function() {
  theme(aspect.ratio = 1,
        panel.background  = element_rect(fill = "#edf0fa",colour = "gray15", size = 2),
        panel.grid.major = element_line(),
        panel.grid.minor = element_line(),
        legend.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        axis.title = element_text(size = 18,hjust = 0.5),
        panel.spacing = unit(1, "lines"),
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 16, face = "bold",hjust = 0)) }

#############################################################################################
# Data Import and Formating
#############################################################################################

# Info on samples
info<-read.csv("/export/lv6/user/pramond/jojy/metadata.csv", header = TRUE, sep = ";")
info$SF<-factor(info$SF, levels = c("0.2","0.8"),labels = c("Free-Living","Particle-Attached"), ordered = TRUE)

# Abundance Table
rpkg<-read.csv("/export/lv6/user/pramond/jojy/coverm_formatted.csv", sep = ";", header = TRUE, row.names = 1)
rownames(rpkg)<-gsub("-",".",rownames(rpkg)) # there was a mismatch between MAGs names in the RPKG and the trait tables due to this character

## Trait Tables from METABOLIC
# import
traits<-read.xlsx("/export/lv6/user/pramond/jojy/METABOLIC_result.xlsx", sheet = 1)
traits<-separate_longer_delim(cols= Hmm.file, data = traits, ", ") # get as rows/traits as there are hmms
# trait information table
trait.info<-traits[, c(1:10)]
trait.info$ID<-make.unique(gsub(".hmm","",trait.info$Hmm.file), sep = '_')
# format trait table
traits<-traits[,grep(".Hit.numbers",colnames(traits) )] # select numerical columns
rownames(traits)<-trait.info$ID;colnames(traits)<-gsub(".Hit.numbers","",colnames(traits)) # format trait and MAG names
traits<-as.data.frame(t(traits)) # transpose
# filter out MAGs absent from abundance table
rownames(traits)[!rownames(traits) %in% rownames(rpkg)] # these MAGs are present in the trait table but not in the abundance table
traits<-traits[rownames(traits) %in% rownames(rpkg),] # we remove them
traits<-decostand(traits, "pa") # convert to 0 and 1s
traits<-traits[,colSums(traits)>0] # remove traits absent in all MAGs
# create trait-based distance matrix
trait.dist<-parDist(as.matrix(traits), method = "euclidean")

# sanity check
dim(rpkg)
table(colnames(rpkg) %in% info$sample)
table(rownames(rpkg) %in% rownames(traits))

#########################################################################
# Compute and explore Microbial functional Diversity (potential)
#########################################################################

## Compute observed Red
uni<-adiv::uniqueness(comm = t(rpkg), dis = trait.dist, abundance = TRUE)
uni$red # coefficients N (species richness), Q (quadratic diversity), D (Simpson diversity), U=Q/D (uniqueness), R=1-U (redundancy), and Pavoine and Ricotta (2019) Ustar=(1-D)/(1-Q) (uniqueness) and Rstar=1-Ustar (redundancy); in this third data frame, coefficients are columns and communities are rows, the coefficients are thus calculated per community only.
# get info on the samples
fred<-merge(uni$red, info, by.x = "row.names", by.y = "sample")

## explore the patterns of FD
ggplot(fred, aes(x = Q, y = Rstar, size = N, color = SF, shape = Bay) )+
  labs(x = "Functional Diversity (Q)", y = "Functional Redundancy (Rstar)")+
  geom_point(alpha = 0.75)+
  scale_size(name = "# of MAGs", range = c(3,9))+
  scale_color_manual(name = "Fraction", values = c("#f59898", "#991c1c") )+
  theme_depth()
    

## explore the patterns of FD
ggplot(fred, aes(x = Q, y = Rstar, color = Salinity, shape = Bay) )+
  labs(x = "Functional Diversity (Q)", y = "Functional Redundancy (Rstar)")+
  geom_point(alpha = 0.75, size = 4)+
  #scale_color_manual(name = "Fraction" )+
  theme_depth()

## explore the patterns of FD
ggplot(fred, aes(y = Rstar, x = Bay) )+
  labs(x = "", y = "Functional Redundancy (Rstar)")+
  geom_boxplot()+
  #scale_color_manual(name = "Fraction" )+
  theme_depth()

## explore the patterns of FD
ggplot(fred, aes(y = Rstar, x = SF) )+
  labs(x = "", y = "Functional Redundancy (Rstar)")+
  geom_boxplot()+
  #scale_color_manual(name = "Fraction" )+
  theme_depth()

#########################################################################
# Drivers of functional redundancy
#########################################################################

## FRED correlation with env. variables
# divide by size-fraction
cor.fred.PA<-na.omit(data.frame(t(cor(fred[fred$SF ==  "Particle-Attached","Rstar"],fred[fred$SF ==  "Particle-Attached", 14:ncol(fred)] ))));colnames(cor.fred.PA)<-c("R2");cor.fred.PA$SF<-"Particle-Attached"
cor.fred.FL<-na.omit(data.frame(t(cor(fred[fred$SF ==  "Free-Living","Rstar"],fred[fred$SF ==  "Free-Living", 14:ncol(fred)] ))));colnames(cor.fred.FL)<-c("R2");cor.fred.FL$SF<-"Free-Living"
# merge results
cor.fred<-as.data.frame(rbind(cor.fred.PA,cor.fred.FL))
# format
cor.fred$env<-gsub("1","",rownames(cor.fred))
cor.fred<-cor.fred[cor.fred$env != c("Latitude", "Longitude"),]
cor.fred$env<-factor(cor.fred$env,rev(c("Temperature", "Salinity", "Depth", "nCells","SE","ChlA","Nitrate","Ammonium","Phosphate", "Silicate" )), ordered = TRUE)

# plot it
genv<-ggplot(cor.fred, aes(x=R2, y = env, color = R2))+
  facet_grid(.~SF)+
  geom_point(size = 4)+
  labs(title ="Correlation with Functional Redundancy (Rstar)", y = "" ,x = bquote('R'^2))+
  geom_segment(aes(yend=env), xend=0, size = 2) +
  theme_bw()+
  scale_color_gradientn(colours = c("coral3","gray85", "#011945"), guide = "none")+
  theme(axis.title = element_text(size = 18,hjust = 0.5))+
  theme(plot.title = element_text(size = 18, face = "bold",hjust = 0))+
  theme(plot.subtitle  = element_text(size = 18,hjust = 0))+
  theme(panel.spacing = unit(1, "lines"))+
  theme(strip.background = element_blank())+
  theme(strip.text =  element_text(size = 16))+
  theme(axis.text = element_text(size = 16))+
  theme(legend.text = element_text(size = 16))+
  theme(legend.title = element_text(size = 16, face = "bold"))+
  theme(legend.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA))+
  theme(plot.title = element_text(size = 16, face = "bold",hjust = 0));genv
