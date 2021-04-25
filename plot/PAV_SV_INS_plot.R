library("dplyr")
library("ggplot2")
library("ggtree")
library("reshape2")
library("ComplexHeatmap")
library("scales")

root <- "/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
setwd(paste0(root, "/script/plot"))

# <---------------------------------------------------------------- unfiltered ---------------------------------------------------------------->
summary <- read.csv(paste0(root, "/result/blastn/summary.csv"), header = TRUE, comment.char = "#")
summary <- summary[-4,]
colnames(summary)[1] <- "Type"
summary <- melt(summary, id.vars = "Type", variable.name = "SampleName", value.name = "Count")

p <- ggplot(summary,aes(x=SampleName)) + 
  theme_set(theme_bw()) + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.text.x = element_blank()) + 
  theme(axis.ticks.x = element_blank()) + 
  theme(legend.title = element_blank()) + 
  xlab("185 samples") + 
  scale_fill_discrete(name="Type", breaks=c("mapped_novel", "mapped_ref", "no_mapped"), labels=c("Mapped to novel sequence", "Mapped to reference", "No mapped")) +
  ggtitle("INS detected by SV caller(Manta) mapped to pan-genome")

p + geom_bar(aes(y=Count,fill=factor(Type)),position="stack",stat = "identity", width =1) + 
  scale_y_continuous(expand = expansion(mult=c(0,0.01)))

p + geom_bar(aes(y=Count,fill=factor(Type)),position="fill",stat = "identity", width =1) + 
  scale_y_continuous(labels = percent, expand = expansion(mult=0)) + 
  ylab("Percentage")

# <---------------------------------------------------------------- 90% coverage filtered ---------------------------------------------------------------->
summary <- read.csv(paste0(root, "/result/blastn/summary.0.9.csv"), header = TRUE, comment.char = "#")
summary <- summary[-4,]
colnames(summary)[1] <- "Type"
summary <- melt(summary, id.vars = "Type", variable.name = "SampleName", value.name = "Count")

p <- ggplot(summary,aes(x=SampleName)) + 
  theme_set(theme_bw()) + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.text.x = element_blank()) + 
  theme(axis.ticks.x = element_blank()) + 
  theme(legend.title = element_blank()) + 
  xlab("185 samples") + 
  scale_fill_discrete(name="Type", breaks=c("mapped_novel", "mapped_ref", "no_mapped"), labels=c("Mapped to novel sequence", "Mapped to reference", "No mapped")) +
  ggtitle("INS detected by SV caller(Manta) filtered with 90% coverage threshold mapped to pan-genome")

p + geom_bar(aes(y=Count,fill=factor(Type)),position="stack",stat = "identity", width =1) + 
  scale_y_continuous(expand = expansion(mult=c(0,0.01)))

p + geom_bar(aes(y=Count,fill=factor(Type)),position="fill",stat = "identity", width =1) + 
  scale_y_continuous(labels = percent, expand = expansion(mult=0)) + 
  ylab("Percentage")







# Read in blastn result file
blastn.out <- paste0(root, "/result/blastn/WGC012904D.WGC013444D/germlineSV.blastn.out.filtered")
blastn.out <- read.table(blastn.out, header = FALSE, stringsAsFactors = FALSE)
colnames(blastn.out) <- c("QSEQID", "SSEQID", "PIDENT", "LENGTH", "MISMATCH", 
                          "GAPOPEN", "QSTART", "QEND", "SSTART", "SEND", "EVALUE", "BITSCORE", "QLEN")

# filter max bitscore
max.score <- blastn.out %>% group_by(QSEQID) %>% summarise(MAX = max(BITSCORE))
blastn.out <- do.call(rbind, lapply(1:nrow(max.score), function(x) {
  subset(blastn.out, (QSEQID == unlist(max.score[x, "QSEQID"])) & (BITSCORE == unlist(max.score[x, "MAX"])))
}))
blastn.out <- blastn.out %>% filter(grepl("GC", SSEQID))

