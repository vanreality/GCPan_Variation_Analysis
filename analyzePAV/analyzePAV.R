library("dplyr")
setwd("/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/script/analyzePAV")

root <- "/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"

ref.normal <- read.table(paste0(root, "/data/refPAV/ref/Normal.ref.pav"), header = TRUE)
ref.tumor  <- read.table(paste0(root, "/data/refPAV/ref/Tumor.ref.pav"), header = TRUE)
pan.normal <- read.table(paste0(root, "/data/refPAV/pan/Normal.pan.pav"), header = TRUE)
pan.tumor  <- read.table(paste0(root, "/data/refPAV/pan/Tumor.pan.pav"), header = TRUE)

hupan.normal <- read.table(paste0(root, "/data/PAV/final.normal.cds.cov"), header = TRUE)
rownames(hupan.normal) <- hupan.normal[,1]
hupan.normal <- hupan.normal[,-1]
colnames(hupan.normal) <- read.table(paste0(root, "/data/PAV/normal.individual.info"), header = FALSE)[,1]
hupan.normal <- data.frame(ifelse(hupan.normal > 0.8, 1, 0))

hupan.tumor <- read.table(paste0(root, "/data/PAV/final.tumor.cds.cov"), header = TRUE)
rownames(hupan.tumor) <- hupan.tumor[,1]
hupan.tumor <- hupan.tumor[,-1]
colnames(hupan.tumor) <- read.table(paste0(root, "/data/PAV/tumor.individual.info"), header = FALSE)[,1]
hupan.tumor <- data.frame(ifelse(hupan.tumor > 0.8, 1, 0))

# ref.normal <- ref.normal[rownames(hupan.normal),]
# ref.tumor  <- ref.tumor[rownames(hupan.tumor),]
# pan.normal <- pan.normal[rownames(hupan.normal),]
# pan.tumor  <- pan.tumor[rownames(hupan.tumor),]

# ref.normal <- ref.normal[-grep("PAR_Y",rownames(ref.normal)), ]

compare <- ref.normal - pan.normal
plot(sort(apply(compare, 2, function(x){sum(x==1)})), ylab = "MAP2REF=Prensence, MAP2PAN=Absence", xlab = "185 samples")
plot(sort(apply(compare, 2, function(x){sum(x==-1)})), ylab = "MAP2REF=Absence, MAP2PAN=Presence", xlab = "185 samples")

compare <- ref.tumor - pan.tumor
plot(apply(compare, 2, function(x){sum(x==1)}), ylab = "MAP2REF=Prensence, MAP2PAN=Absence", xlab = "185 samples")
plot(apply(compare, 2, function(x){sum(x==-1)}), ylab = "MAP2REF=Absence, MAP2PAN=Presence", xlab = "185 samples")

compare <- pan.normal - hupan.normal


