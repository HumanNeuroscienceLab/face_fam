
if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))
library(plyr)
suppressMessages(library(doMC))
registerDoMC(24)

# Base Paths
preproc.base <- "/data1/famface01/analysis/preprocessed"
roi.base     <- "/data1/famface01/analysis/roi/Functional_v3"
task.base    <- "/data1/famface01/analysis/task_activity"

load.cur <- function(subj) {
  indir <- "/data1/famface01/analysis/encoding/ShapeAnalysis/data"
  infile <- sprintf("%s/roi_n_more_%s.rda", indir, subj)
  dat <- NULL
  load(infile)
  dat
}



##

dat      <- load.cur("sub03")
trait.df <- dat$features$traits
trait.pca <- prcomp(trait.df, scale=F, retx=T)
trait.pca$rotation[,1]

base <- "/data1/famface01/analysis/encoding/12_Features"
traits.df <- read.csv(sprintf('%s/personality_traits_grpave.csv', base))
vnames <- traits.df[,1]
traits.df <- traits.df[,-1]
dim(traits.df)
res <- rpca::rpca(scale(traits.df, scale=F, center=T))
head(res$L)
cor(res$L)
cor(traits.df)

trait.pca <- prcomp(traits.df, scale=F, retx=T)

# gives much fewer!
# TODO: try ortho=T option
library(PMA)
cv.out <- SPC.cv(as.matrix(traits.df), sumabsvs=seq(1.2, sqrt(ncol(traits.df)), len=6))
print(cv.out)
plot(cv.out)
out <- SPC(as.matrix(traits.df), sumabsv=cv.out$bestsumabs, K=10)
rownames(out$v) <- colnames(traits.df)
round(out$v, 4)

library(psych)
vss(traits.df)

library(corrplot)
col <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                          "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                          "#4393C3", "#2166AC", "#053061"))(200)
corrplot(cor(traits.df), order = "hclust", tl.col='black', tl.cex=.75, diag=F, 
         col=rev(col)) 

# non-sparse version (http://web.stanford.edu/class/psych253/tutorials/FactorAnalysis.html)
res1 = factanal(traits.df, factors = 6, rotation = "varimax", na.action = na.omit, 
               scores="regression")
# very highly related
diag(cor(res1$scores, res2$scores))
res$loadings
res$rotmat

barplot(res1$uniquenesses, ylab="Uniqueness", main="Lower means more unique (less shared)")
res1$loadings[,1]

# sparse version of life
library(fanc)
ret <- fanc(as.matrix(traits.df), 5, cor.factor=T)
print(ret) #print candidates of gamma and rho
#output for fixed tuning parameters
out(ret, rho=0.1, gamma=Inf)
#select a model via model selection criterion
ret2 <- select(ret, criterion="BIC", gamma=Inf, scores=T)
#plot solution path
#plot(fit)


# NMF: https://cran.r-project.org/web/packages/NMF/vignettes/NMF-vignette.pdf
# https://www.r-bloggers.com/extracting-latent-variables-from-rating-scales-factor-analysis-vs-nonnegative-matrix-factorization/
library(NMF)
X <- as.matrix(traits.df); nres <- nmf(X, 5)
