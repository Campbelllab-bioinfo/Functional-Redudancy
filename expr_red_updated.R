#expressed functional redudancy
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

#load data
# Abundance Table
rpkg <- read.xlsx("/home/jojyj/Desktop/MTROL_MS1_final/Functional_redudancy/expressed_fr/abd.xlsx")## Import and format expression data
rownames(rpkg) <- rpkg$Genome

#load TPM values
expr2 <- read.xlsx("/home/jojyj/Desktop/MTROL_MS1_final/Functional_redudancy/expressed_fr/TPM_new_sep25.xlsx")


#Pre_processing
# list of MAGs present in expression but not in RPKG
unique(expr2[!expr2$Genome %in% rownames(rpkg), "Genome"]) 
# list of MAGs present in RPKG but not in expression
rownames(rpkg[!rownames(rpkg) %in% unique(expr2$Genome), ])

## continue formatting 
mexpr2 <- melt(expr2[, c(grep("RNA|Genome", colnames(expr2)))], id.vars = c("Genome")) # long format 
mexpr2$sample<-gsub("_RNA1|_RNA","",mexpr2$variable)
#mexpr2$sample <- mexpr2$variable  # Keep the original sample names
mexpr2 <- mexpr2[mexpr2$value > 0, ]  # removing unexpressed traits (filtered the rows)

## create a list of MAG x Trait table 
lsp2 <- list() # empty list to collect results

for (i in unique(mexpr2$sample)) { # for all samples
  # make table using dcast
  mat <- dcast(formula = Genome ~ variable, data = mexpr2[mexpr2$sample == i, ], value.var = "value", fill = 0, fun.aggregate = sum)
  
  rownames(mat) <- mat$Genome  # Set row names
  mat <- mat[, -1]  # Remove the Genome column
  
  lsp2[[i]] <- decostand(mat, method = "pa")  # Convert to presence/absence and save
}

## (working fine till here)
sample_match <- names(lsp2) %in% colnames(rpkg)  # Check for matching names
print(sample_match)
colnames(rpkg)

# Clean up RPKG column names for consistency
colnames(rpkg) <- gsub("classified_|_DNA|.bac.ar", "", colnames(rpkg))
#names(lsp2) <- gsub("_", "", names(lsp2))  # Ensure names in lsp2 are consistent
colnames(rpkg)

# Subset RPKG based on matching samples
rpkg2 <- rpkg[, colnames(rpkg) %in% names(lsp2)]  
colnames(rpkg)

# compute expressed FRed
expr.fred2<-NULL


for (i in colnames(rpkg)){
  # rpkg
  rpsp<-rpkg[,i, drop=F]
  mags<-rownames(rpsp)[rowSums(rpsp)>0] # get MAGs present (RPKG) in the sample
  
  # select id for tpm
  # we can only work with samples that have an equivalent in expression data
  if(i %in% names(lsp2) == FALSE ){next} 
  
  # select tpm and rpkg
  tpm.sp<-lsp2[[i]]
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
  expr.fred2<-rbind(expr.fred2,uni)# save the results
  
  # print progress
  print(paste("Done with ", match(i, colnames(rpkg))," / ", ncol(rpkg), sep = ""))
}

##Error in rowSums(rpsp) : 'x' must be numeric

# Get info on the samples where expressed FRED was computed
fred2 <- fred  # Assuming 'fred' is defined earlier in your script
fred2$ID <- gsub("classified_|_DNA\\.bac\\.ar|\\.bac\\.ar", "", fred2$Row.names)  # Clean IDs
expr.fred3 <- merge(expr.fred2, fred2, by.x = "row.names", by.y = "ID", all.x = TRUE, suffixes = c(".expr", ".pot"))
