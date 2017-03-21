
# Setup -------------------------------------------------------------------

# Path for my own packages
if (!any(.libPaths() == "/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2")) .libPaths(c("/home/zshehzad/R/x86_64-redhat-linux-gnu-library/3.2", .libPaths()))
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))

# General workhorse function
library(plyr)

# Parallelization
suppressMessages(library(doMC))
registerDoMC(24)

# If we were to use ggplot2 and other plotting needs
library(RColorBrewer)
library(ggplot2)
source("/data1/famface01/command/encoding/FirstPass/plot_functions.R")

# This gets us the simple_lm function for faster GLMs
source("/data1/famface01/command/encoding/SubRois_Unfam/lm_functions.R")

subjects <- sprintf("sub%02i", 1:6) # again skipping sub01
#indir    <- "/data1/famface01/analysis/encoding/ShapeAnalysis/data"
rnames   <- c("r.vATL","r.aFFA", "r.pFFA", "r.OFA", "r.EBA", 'r.Amyg', 
              "l.vATL", "l.FFA", "l.OFA", "l.EBA", 'l.Amyg')



# Load --------------------------------------------------------------------

## Time-Series Data
#load("/data1/famface01/analysis/misc/320_roi_task_activity/std_roi_dat.rda", verbose=T)

# Beta Data
load("/data1/famface01/analysis/misc/320_roi_task_activity/dgamma_betas.rda", verbose=T)

# Rnames
rnames2 <- sub("[.]", " ", rnames)
rnames2 <- sub("^r", "R", rnames2)
rnames2 <- sub("^l", "L", rnames2)


# Plot --------------------------------------------------------------------

library(dynamicTreeCut)
d    <- dist(grp.bs1[,,2])
dmat <- as.matrix(d)
hc   <- hclust(d, method="ward.D2")

ind  <- order.dendrogram(as.dendrogram(hc))
cl   <- cutreeDynamic(hc, distM=dmat) # this will be 5
cl   <- cutree(hc, k=6)
ucl  <- unique(cl)

rcols <- as.character(factor(cl, levels=sort(ucl), brewer.pal(length(ucl), "Set3")))
table(cl)
cols <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
seplines <- sum(grepl("^R", rnames2)) + 0.5
heatmap(grp.bs1[,,2], Rowv=as.dendrogram(hc), Colv=NA, scale="none", 
        labRow=NA, col=cols, labCol=rnames2, RowSideColors=rcols, 
        add.expr = abline(v=seplines, lwd=4, lty=3, col="grey70"))


mat <- grp.bs1[ind,,2]
bvnames2 <- bvnames[ind]
cl2 <- cl[ind]
dmat2 <- dmat[ind,]
sort(rowMeans(dmat2[cl2==1,cl2==1]))[1:10]
sort(rowMeans(cor(t(mat[cl2==1,]))), dec=T)[1:10]

cl2.vids <- lapply(unique(cl2), function(ci) names(sort(rowMeans(dmat2[cl2==ci,cl2==ci]))[1:5]))

head(bvnames[ind])
single.files <- Sys.glob(sprintf("/data1/famface01/data/stimuli/vids/single_frames/%s_*.jpg", bvnames2))

images <- llply(single.files, load.image, .progress="text")
images <- llply(images, resize, 640, 360, .progress="text") # have size be all the same
images2 <- llply(images, resize, 360, 180, .progress="text") # have size be all the same
images2 <- llply(images2, function(img) squeeze(as.array(img)), .progress="text")

plot(images[[30]])

ind.add <- cumsum(c(0, rle(cl2)$lengths[-length(rle(cl2)$lengths)]))

library(abind)

table(cl2)
unique(cl2)
ncol   <- 10
nrows  <- ceiling(table(cl2)[unique(cl2)]/ncol)
nextra <- ncol*nrows - table(cl2)[unique(cl2)]
nextra[nextra==0] <- ncol

odir <- "/data1/famface01/figures/bsclust_show-vid-pics"
blank <- squeeze(as.array(as.cimg(array(0, c(360,180,3)))))

for (ci in 1:length(unique(cl2))) {
  cat(ci, "\n")
  
  starts <- seq(1,ncol*nrows[ci],by=ncol)
  tmp  <- lapply(1:length(starts), function(i) {
    si <- starts[i] + ind.add[ci]
    if (i == length(starts)) {
      ei <- si+nextra[ci]-1
    } else {
      ei <- (si+ncol-1)
    }
    imgs <- images2[si:ei]
    if (length(imgs) < ncol) {
      for (j in 1:(ncol-nextra[ci])) {
        imgs <- append(imgs, list(blank))
      }
    }
    abind(imgs, along=1)
  })
  tmp2 <- abind(tmp, along=2)
  dim(tmp2)
  
  if (!file.exists(odir)) dir.create(odir)
  
  png(filename=file.path(odir, sprintf("clust%i-%i.png", unique(cl2)[ci], ci)), units="in", width=10, height=(10/dim(tmp2)[1] * dim(tmp2)[2]), pointsize=12, res=150)
  
  plot(as.cimg(tmp2), bty='n', axes=FALSE, ann=FALSE, xaxs="i", yaxs="i")
  title(main=sprintf("Cluster %i - %i", unique(cl2)[ci], ci))
  abline(v=seq(1,dim(tmp2)[1]+1,length.out=ncol+1), h=seq(1,dim(tmp2)[2]+1, length.out=nrows[ci]+1))
  
  dev.off()
}


# actually we have duplicates
# and they are in fact in the same cluster
which(bvnames2 %in% c("Unknown_21", "Unknown_29"))
cl2[which(bvnames2 %in% c("Unknown_21", "Unknown_29"))]


tmp <- sapply(images2, dim)
tmp[,which(tmp[1,]!=854)]
854/640

resize()


cl.single.files <- Sys.glob(sprintf("/data1/famface01/data/stimuli/vids/single_frames/%s_*.jpg", unlist(cl2.vids)))

library(imager)
img <- load.image(cl.single.files[1])
images <- lapply(cl.single.files, load.image)
plot(images[[30]])
dev.new()
dev.off()
tmp <- images[[1]] %>% as.array %>% squeeze
