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

#########################################################################
# Compute and explore Expressed Microbial functional Diversity
#########################################################################

## import and format expr data
expr<-read.xlsx("/export/lv6/user/pramond/jojy/trans_new_all.xlsx")
expr<-expr[,-grep("RNA2",colnames(expr))] # I suppose RNA1 and RNA2 are just paired end reads of the same sample (?)
# => do no include duplicates and please rename samples to a common ID with table info
expr$MAG<-gsub("\\_NODE.*","",expr$Gene.ID) # removing useless info from MAGs names
expr$MAG<-gsub("_","_",expr$MAG) # formating MAGs names
expr$MAG<-gsub("-",".",expr$MAG)

## Problem with MAGs not present in RPKG (or not matching names)
# list of MAGs present in expression but not in RPKG
unique(expr[unique(expr$MAG) %in% rownames(rpkg) == FALSE, "MAG"]) 
# list of MAGs present in RPKG but not in expression
rownames(rpkg[rownames(rpkg) %in% unique(expr$MAG) == FALSE, ]) #
# => homogenize MAG names across expression and RPKG

## continue formating nevertheless
mexpr<-melt(expr[, c(grep("bam|MAG|KOS", colnames(expr)))], id.vars = c("MAG","KOS")) # long format 
mexpr$sample<-gsub("_RNA.fq.bam.sorted|_RNA1.fq.bam.sorted|R1.tr50_","",mexpr$variable) # removing useless info on samples names
mexpr<-mexpr[mexpr$value >0,-c(3)] # removing unexpressed traits

## create a list of MAG x Trait table for each sample
lsp<-list() # empty list to collect results
for (i in unique(mexpr$sample)){ # for all samples
  # make table
  mat<-dcast(formula = MAG ~ KOS ,data = mexpr[mexpr$sample == i,], value.var = c("value"), fill = 0, fun.aggregate = sum)
  rownames(mat)<-mat$MAG;mat<-mat[,-1] # format table
  lsp[[i]]<-decostand(mat, method = "pa") # convert to presence absence and save
}

# => problem of the match between samples names in expression and RPKG
# => please homogenize sample and MAG names across analyses
names(lsp) %in% colnames(rpkg)
colnames(rpkg)<-gsub("classified_|Bay|.bac.ar|_DNA|_", "", colnames(rpkg))
names(lsp)<-gsub("_", "", names(lsp))

# We can only work with those samples that are in both RPKG and expression
rpkg2<-rpkg[,colnames(rpkg) %in% names(lsp)]

# compute expressed FRed
# loop to subselect expression and compute expressed Fred per sample
expr.fred<-NULL
for (i in colnames(rpkg)){
  # rpkg
  rpsp<-rpkg[,i, drop=F]
  mags<-rownames(rpsp)[rowSums(rpsp)>0] # get MAGs present (RPKG) in the sample

  # select id for tpm
  # we can only work with samples that have an equivalent in expression data
  if(i %in% names(lsp) == FALSE ){next} 

  # select tpm and rpkg
  tpm.sp<-lsp[[i]]
  tpm.sp<-tpm.sp[rownames(tpm.sp) %in% mags,] # we work only with MAGs present (RPKG) in the sample
  tpm.sp<-tpm.sp[, colSums(tpm.sp)!= 0] # we remove MAGs absent from expression
  tpm.sp<-decostand(tpm.sp, "pa") # presence/absence of expression
  rpsp<-rpsp[rownames(rpsp) %in% rownames(tpm.sp),,drop = F] # work only with MAGs present in expression data

  # compute functional distance based on trait expressed (0/1)
  tpm.dist<-parDist(as.matrix(tpm.sp), method = "euclidean")

  # compute uniqueness (weighted by abundance in RPKG)
  uni<-cbind(adiv::uniqueness(comm = t(rpsp), dis = tpm.dist, abundance = TRUE)$red, length(mags),nrow(tpm.sp) )
  # get number of MAGs in RPKG and in expression data = those included in computation
  colnames(uni)[c(ncol(uni)-1, ncol(uni))]<-c("n.rpkg.MAGs", "n.expr.MAGs")
  expr.fred<-rbind(expr.fred,uni)# save the results

  # print progress
  print(paste("Done with ", match(i, colnames(rpkg))," / ", ncol(rpkg), sep = ""))
}

# get info on the samples where expressed FRED was computed
fred2<-fred
fred2$ID<-gsub("classified_|Bay|_DNA.bac.ar|_|.bac.ar","",fred2$Row.names)
expr.fred2<-merge(expr.fred, fred2, by.x = "row.names", by.y = "ID", all.x = TRUE, suffixes = c(".expr", ".pot"))

## explore the patterns of FD
ggplot(expr.fred2, aes(x = Q.pot, y = Q.expr, size = N.expr, color = SF, shape = Bay) )+
  labs(x = "Functional Diversity (Q)", y = "Expressed Functional Diversity (Q)")+
  geom_point(alpha = 0.75)+
  scale_size(name = "# of MAGs\n(expressing traits)", range = c(3,9))+
  scale_color_manual(name = "Fraction", values = c("#f59898", "#991c1c") )+
  theme_depth()

## explore the patterns of FD
ggplot(expr.fred2, aes(x = Rstar.pot, y = Rstar.expr, size = N.expr, color = SF, shape = Bay) )+
  labs(x = "Functional Redundancy (Rstar)", y = "Expressed Functional Redundancy (Rstar)")+
  geom_point(alpha = 0.75)+
  scale_size(name = "# of MAGs\n(expressing traits)", range = c(3,9))+
  scale_color_manual(name = "Fraction", values = c("#f59898", "#991c1c") )+
  theme_depth()



