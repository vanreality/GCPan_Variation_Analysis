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

# <---------------------------------------------------------------- 90% coverage filtered ---------------------------------------------------------------->
summary <- read.csv(paste0(root, "/result/blastn/summary.0.9.csv"), header = TRUE, comment.char = "#")

# <---------------------------------------------------------------- 90% coverage filtered and TRF masked ---------------------------------------------------------------->
summary <- read.csv(paste0(root, "/result/blastn/summary.masked.csv"), header = TRUE, comment.char = "#")

# <---------------------------------------------------------------- 90% coverage filtered and TRF masked ver 2 ---------------------------------------------------------------->
summary <- read.csv(paste0(root, "/result/blastn/summary.masked.csv"), header = TRUE, comment.char = "#")
summary <- summary[-4,]
summary <- summary[-3,]
colnames(summary)[1] <- "Type"
summary <- melt(summary, id.vars = "Type", variable.name = "SampleName", value.name = "Count")

p <- ggplot(summary,aes(x=SampleName)) + 
  theme_set(theme_bw()) + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) + 
  theme(axis.text.x = element_blank()) + 
  theme(axis.ticks.x = element_blank()) + 
  theme(legend.title = element_blank()) + 
  xlab("Samples (n=185)") + 
  scale_fill_manual(values = c("#EEB422", "#8968CD"), name="Type", breaks=c("mapped_novel", "mapped_ref"), labels=c("Novel sequence", "Reference"))

#c("#FFED6F","#BEBADA" )
#c("#EEB422", "#8968CD")
pdf("INS_A.pdf", width = 9.5, height = 2.4)
p + geom_bar(aes(y=Count,fill=factor(Type)),position="stack",stat = "identity", width =1) + 
  scale_y_continuous(expand = expansion(mult=c(0,0.01)))
dev.off()

pdf("INS_P.pdf", width = 9.5, height = 2.4)
p + geom_bar(aes(y=Count,fill=factor(Type)),position="fill",stat = "identity", width =1) + 
  scale_y_continuous(labels = percent, expand = expansion(mult=0)) + 
  ylab("Percentage")
dev.off()

p = 0
for (x in 1:ncol(summary[,-1])) {
  p <- p + summary[1, x + 1]/summary[2, x + 1]
}
p / 185

sum(summary[1,-1])/185
sum(summary[2,-1])/185

# <---------------------------------------------------------------- new gene covered by SV INS ---------------------------------------------------------------->
samples <- unlist(read.table(paste0(root, "/script/analyzeSV/185samples"), header = FALSE, stringsAsFactors = FALSE))

# GC_NRS_013394 1446 2635  GC001807
# GC_NRS_012564 8967 16379 GC000643

res <- do.call(rbind, lapply(samples, function(x) {
  blastn <- read.table(paste0(root, "/result/blastn/", x, "/germlineSV.blastn.out.filtered"))
  colnames(blastn) <- c("QSEQID", "SSEQID", "PIDENT", "LENGTH", "MISMATCH", 
                        "GAPOPEN", "QSTART", "QEND", "SSTART", "SEND", "EVALUE", "BITSCORE", "QLEN")
  blastn <- blastn %>% filter(grepl("GC_NRS_013394",SSEQID) | grepl("GC_NRS_012564",SSEQID))
  return(blastn)
}))

# <---------------------------------------------------------------- filtered insertion length density ---------------------------------------------------------------->
samples <- unlist(read.table(paste0(root, "/script/analyzeSV/185samples"), header = FALSE))

res <- do.call(cbind, lapply(samples, function(x) {
  blastn <- read.table(paste0(root, "/result/blastn/", x, "/germlineSV.blastn.out.filtered"))
  colnames(blastn) <- c("QSEQID", "SSEQID", "PIDENT", "LENGTH", "MISMATCH", 
                        "GAPOPEN", "QSTART", "QEND", "SSTART", "SEND", "EVALUE", "BITSCORE", "QLEN")
  blastn <- blastn %>% filter(grepl("GC",SSEQID))
  len <- data.frame(blastn$QLEN)
  colnames(len) <- x
  rep <- data.frame(rep(NA,292-nrow(len)))
  colnames(rep) <- x
  return(rbind(len, rep))
}))

pd <- melt(res, variable.name = "samples", value.name = "ins_length", na.rm = TRUE)

ggplot(pd, aes(x=samples, y=ins_length)) + 
  geom_boxplot(outlier.size = .1) + 
  xlab("185 samples") + 
  theme(axis.text.x = element_blank()) + 
  theme(axis.ticks.x = element_blank())

ggplot(pd, aes(x=ins_length, color=samples, group=samples)) +
  geom_density(alpha = .1) + 
  theme(legend.position = "none")










