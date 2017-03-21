
# This script puts together data for the trait/demographic information
# It will get 6 factor scores for the 10 traits 


# Setup -------------------------------------------------------------------

# Setup
if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))
library(plyr)
suppressMessages(library(doMC))
registerDoMC(30)
library(RColorBrewer)
library(dynamicTreeCut)

library(ggplot2)
source("/data1/famface01/command/encoding/FirstPass/plot_functions.R")
source("/data1/famface01/command/encoding/SubRois_Unfam/lm_functions.R")


# Load Trait Data -----------------------------------------------------------

# Load demographics information
base <- "/data1/famface01/analysis/encoding/12_Features"
demos <- read.csv(sprintf('%s/demographics_unfam_df.csv', base))
demos <- demos[,-1]
df.demos0 <- demos[,-1]

# Get video names
demo.vnames <- df.demos0$video # typo
demo.vnames <- sub("_fr[0-9]{3}", "", demo.vnames)

# Merge long beard and beard thing 
df.demos0$facial_hair <- revalue(df.demos0$facial_hair, c("Long beard"="Beard"))
table(df.demos0$facial_hair)

# Create big matrix with all and a smaller one without video/hair
df.demos0[,1] <- demo.vnames
df.demos    <- df.demos0[,-c(1,7)] # remove video, hair

# Load trait information
base         <- "/data1/famface01/analysis/encoding/12_Features"
traits       <- read.csv(sprintf('%s/personality_traits_grpave.csv', base))
trait.vnames <- as.character(traits[,1])
traits       <- traits[,-1]
# Reorder to match with demos
oinds       <- sapply(demo.vnames, function(x) which(trait.vnames == x))
all.equal(trait.vnames[oinds], demo.vnames)
df.traits   <- traits[oinds,]

# Save vnames
vnames <- demo.vnames

# Relevant functions
formula_to_mat <- function(formula, data) {
  mf <- model.frame(formula=formula, data=data)
  X <- model.matrix(attr(mf, "terms"), data=mf)
  return(X)
}

compile.demo.trait <- function(df.demos, df.traits, ret.intercept=TRUE) {
  ## DEMOGRAPHICS
  # Remove the effect of gender from makeup
  df.demos2 <- df.demos
  df.demos2$makeup <- lm(makeup ~ gender, data=df.demos2)$residuals
  # Center age
  df.demos2$age <- scale(df.demos2$age, scale=F, center=T)
  # Set the reference for gender to Male
  df.demos2 <- within(df.demos2, gender <- relevel(gender, ref = "Male"))
  # Set the reference for facial hair to when none
  df.demos2 <- within(df.demos2, facial_hair <- relevel(facial_hair, ref = "None"))
  # Set the reference for race to when white
  df.demos2 <- within(df.demos2, race <- relevel(race, ref = "White"))
  # Set the reference for eyes to Cannot tell
  df.demos2 <- within(df.demos2, eye <- relevel(eye, ref = "Cannot tell"))
  # Get the matrix
  mat.demos1 <- formula_to_mat(~.-1, df.demos2)
  mat.demos2 <- formula_to_mat(~., df.demos2) # note this will have the intercept
  
  ## TRAITS
  mat.traits <- scale(df.traits, scale=F, center=T)
  
  if (ret.intercept) {
    return(cbind(mat.demos2, mat.traits))
  } else {
    return(cbind(mat.demos1, mat.traits))
  }
}

# Compile the demographics and traits
mat.all <- compile.demo.trait(df.demos, df.traits)


# Load Shape/Texture Data ---------------------------------------------------

# Finally load in the shape/texture data (this is )
basedir <- "/data1/famface01/analysis/encoding/12_Features/identity_pca"
sym.shapes <- read.csv(file.path(basedir, "shape_sym_pca_scores.csv"), row.names = 1)
shape.vnames <- rownames(sym.shapes) # vnames to use
sym.shapes.eigs <- read.table(file.path(basedir, "shape_sym_pca_eigs.txt"))
sym.textures    <- read.csv(file.path(basedir, "texture_sym_pca_scores.csv"), row.names = 1)
## combine
X1 <- sym.shapes[,1:65]
X2 <- sym.textures[,1:200]
colnames(X1) <- sprintf("shape_%s", colnames(X1))
colnames(X2) <- sprintf("texture_%s", colnames(X2))
X <- cbind(X1,X2)


# Rearrange rows for shape/texture ------------------------------------------

# rearrange of the shapes here
## get the new inds
inds <- sapply(vnames, function(vname) which(shape.vnames==vname))
all.equal(shape.vnames[inds], vnames)
## set new inds
X0 <- X # save
X <- as.matrix(X[inds,])


# OLD: Reduce Dimensionality -----------------------------------------------------

### I did this before but skip to next section for where I'm at now

# So from this analysis, it seems that maybe should split up traits from other things

# Use factor analysis on demographic/traits
library(psych)

# remove intercept
mat.all2 <- mat.all[,-1]
#mat.all2 <- mat.all2[,c(1,9,15,16:25)]
head(mat.all2)

# Determine the number of components needed
fa.parallel(mat.all2) # suggests 9 comps
vss(mat.all2, 25) # suggests 18 or 21 (going with 21)
# also went with 21 over 18 because it includes middle eastern participants
# also there are local minima at 16 and 11

# Get the factors
fac.res <- psych::fa(mat.all2, nfactors=18, residuals=T, rotate='varimax', 
                     fm='minres')
#loadings(fac.res)
#residuals(fac.res)

# Save the loadings and scores, rename the columns
fac.loads  <- unclass(fac.res$loadings)
fac.scores <- fac.res$scores
colnames(fac.loads) <- sprintf("factor%02i", 1:ncol(fac.loads))
colnames(fac.scores) <- sprintf("factor%02i", 1:ncol(fac.scores))

# Plot
library(corrplot)
hc1     <- hclust(dist(fac.loads))
hc1     <- as.dendrogram(hc1)
ord.hc1 <- order.dendrogram(hc1)
loading <- fac.loads[ord.hc1,]
col <- rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                              "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                              "#4393C3", "#2166AC", "#053061"))(200))
corrplot(loading, tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=T)


# Select Trait Factors ----------------------------------------------------

mat.all2 <- mat.all[,-1]
mat.all2 <- mat.all2[,c(20:29)]
head(mat.all2)

# Determine the number of components needed
psych::fa.parallel(mat.all2) # suggests 5 factors
psych::vss(mat.all2, ncol(mat.all2)) # also suggests 5, 6, 7 or 9 factors

# Get the factors (VSS complexity 2 maximizes at 6)
fac.res <- psych::fa(mat.all2, nfactors=6, residuals=T, rotate='varimax', 
                     fm='minres')

# Save the loadings and scores, rename the columns
fac.loads  <- unclass(fac.res$loadings)
fac.scores <- fac.res$scores
colnames(fac.loads) <- sprintf("factor%02i", 1:ncol(fac.loads))
colnames(fac.scores) <- sprintf("factor%02i", 1:ncol(fac.scores))

# Plot
library(corrplot)
hc1     <- hclust(dist(fac.loads))
hc1     <- as.dendrogram(hc1)
ord.hc1 <- order.dendrogram(hc1)
loading <- fac.loads[ord.hc1,]
col <- rev(colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                              "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                              "#4393C3", "#2166AC", "#053061"))(200))
corrplot(loading, tl.col='black', tl.cex=.75, diag=T, col=col, is.corr=T)

# Assign Names
fac.names <- c("unemotional", "competent", "trustworthy", "typical", 
               "memorable", "attractive")
colnames(fac.scores) <- fac.names


# Other Features ----------------------------------------------------------

mat.all3 <- mat.all[,-1]
mat.all3 <- mat.all3[,-c(20:29)]

# split up
facial.hair <- mat.all3[,c(2:6)]
race <- mat.all3[,9:13]
eye <- mat.all3[,14:18]
misc <- mat.all3[,c(1,7,8,19)] # age, makeup, gender, glasses

# Combine it all back
all.mat <- cbind(fac.scores, misc, facial.hair, race, eye)
all.groups <- rep(c("traits", "age", "makeup", "gender", "glasses", "facial_hair", "race", "eye"), c(ncol(fac.scores), 1, 1, 1, 1, ncol(facial.hair), ncol(race), ncol(eye)))


# Plots -------------------------------------------------------------------

# Basic breakdowns of the factors
table(df.demos$gender)
table(df.demos$facial_hair)
table(df.demos$race)
table(df.demos$glasses)
table(df.demos0$hair)
table(df.demos0$eye)

# Plot all
cols <- brewer.pal(11, "RdBu")
heatmap(all.mat, col=rev(cols), scale="none", labRow=F, margins=c(7,5))
## other version
hc.cols <- hclust(dist(t(all.mat)), method="ward.D2")
hc.rows <- hclust(dist(all.mat), method="ward.D2")
heatmap(all.mat, col=rev(cols), Colv=as.dendrogram(hc.cols), Rowv=as.dendrogram(hc.rows), 
        scale="none", labRow=F, margins=c(7,5))



# Save ------------------------------------------------------------------- 

# we save the 0th demos with more stuff in it
write.csv(df.demos0, "/data1/famface01/analysis/misc/320_roi_task_activity/12_facefeats_df-demos.csv")
write.table(demo.vnames, "/data1/famface01/analysis/misc/320_roi_task_activity/12_facefeats_demo-vnames.txt", 
            row.names=F, col.names=F)

write.csv(all.mat, "/data1/famface01/analysis/misc/320_roi_task_activity/12_facefeats_all-mat.csv")
write.table(all.groups, "/data1/famface01/analysis/misc/320_roi_task_activity/12_facefeats_all-groups.txt", 
            row.names=F, col.names=F)

write.csv(X, "/data1/famface01/analysis/misc/320_roi_task_activity/12_facefeats_sym-shape-texture.csv")

