library("dplyr")
setwd("/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/script/analyzePAV")

cov1 <- read.table("/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/data/tmp/output/data/WGC013277D.sta", header = FALSE)
cov2 <- read.table("/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/data/tmp/WGC013277D/WGC013277D.bam.cov", header = FALSE)

cov1 <- cov1[,c(1, 5, 8)]
cov2 <- cov2[,c(1, 5, 6)]

cov1 <- cov1 %>% filter(V1 %in% cov2$V1)
cov1 <- cov1 %>% filter(V1 %in% unique(cov1$V1))

dup <- unique(cov1[duplicated(cov1$V1),1])
cov1 <- cov1 %>% filter(!(V1 %in% dup))
cov2 <- cov2 %>% filter(!(V1 %in% dup))

diff_gene <- cov2 %>% filter(!(V1 %in% cov1$V1)) #difference between gtf and gff
cov2 <- cov2 %>% filter(!(V1 %in% diff_gene$V1))

colnames(cov1) <- c("ID", "GENE", "CDS")
colnames(cov2) <- c("ID", "GENE", "CDS")

cov2$GENE <- round(cov2$GENE, 4)
cov2$CDS <- round(cov2$CDS, 4)

compare <- cbind(cov1, cov2)
compare$GENEDIFF <- compare[,2] - compare[,5]
compare$CDSDIFF <- compare[,3] - compare[,6]
