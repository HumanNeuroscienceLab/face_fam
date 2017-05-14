library(plyr)
library(bigmemory)

# Load info and regressors ------------------------------------------------

# For each frame in the video, we extracted the features from each layer of the neural network
indir <- "/data1/famface01/analysis/misc/openface/layer_features"
infiles <- list.files(indir)
vidname <- sub("_fr[0-9]{3}.mat", "", infiles)
frame   <- as.integer(sub(".*fr", "", sub(".mat", "", infiles)))
vdf     <- data.frame(ind=1:length(infiles), vid=vidname, fr=frame, fn=infiles)
## select only the unfamiliar faces
base <- "/data1/famface01/analysis/encoding/12_Features"
tmp     <- read.csv(sprintf('%s/personality_traits_grpave.csv', base))
vdf     <- vdf[vdf$vid %in% as.character(tmp$X),]
vdf$ind <- 1:nrow(vdf)

# Load demographics information
base <- "/data1/famface01/analysis/encoding/12_Features"
demos <- read.csv(sprintf('%s/demographics_unfam_df.csv', base))
demos <- demos[,-1]
demo.vnames <- sub("_fr[0-9]{3}", "", demos$video)
# Reorder rows to match features
oinds       <- sapply(vdf$vid, function(x) which(demo.vnames == x))
all.equal(demo.vnames[oinds], as.character(vdf$vid))
df.demos    <- demos[oinds,-c(1:2)]

# Load trait information
base         <- "/data1/famface01/analysis/encoding/12_Features"
traits       <- read.csv(sprintf('%s/personality_traits_grpave.csv', base))
trait.vnames <- as.character(traits[,1])
traits       <- traits[,-1]
# Reorder
oinds       <- sapply(vdf$vid, function(x) which(trait.vnames == x))
all.equal(trait.vnames[oinds], as.character(vdf$vid))
df.traits   <- traits[oinds,]

# Get pose information
pbase <- "/data1/famface01/data/stimuli/vids/eight_frames/three_dee"
df.pose <- ldply(unique(vdf$vid), function(vname) {
  fpath <- sprintf("%s/%s/%s_pose.txt", pbase, vname, vname)
  tab <- read.table(fpath)
  colnames(tab) <- c("pitch", "yaw", "roll")
  tab$fd <- Rfast::dista(c(0,0,0), t(tab))
  tab$vid <- vname
  tab
}, .progress="text")
all.equal(as.character(df.pose$vid), as.character(vdf$vid))
df.pose <- df.pose[,colnames(df.pose) != "vid"]

# Get luminance information
base <- "/data1/famface01/analysis/misc/400_nnet_specificity"
df.lum <- read.csv(file.path(base, "feats_ave_luminance.csv"))
lum.vnames <- as.character(df.lum$vid)
df.lum <- df.lum[,-1]
# Reorder
oinds       <- sapply(1:nrow(vdf), function(i) {
  which((lum.vnames == vdf$vid[i]) & (df.lum$frame == vdf$fr[i]))
})
all.equal(lum.vnames[oinds], as.character(vdf$vid))
df.lum   <- df.lum[oinds,"luminance",drop=F]

# Combine
df.info <- cbind(vdf, df.pose, df.lum, df.demos, df.traits)
df.info$vid <- factor(df.info$vid)

# Save --------------------------------------------------------------------

setwd("/data1/famface01/command/misc/face_representations/322_asap_rois")
write.csv(vdf2, file="120_gender/z_df_info.csv")
