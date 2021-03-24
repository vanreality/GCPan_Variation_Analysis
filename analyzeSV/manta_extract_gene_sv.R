library("dplyr")

args <- commandArgs(T)

# Read in manta vcf file
manta.sv.vcf <- args[1]
# manta.sv.vcf <- "/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/result/manta/WGC012904D.WGC013444D/result/results/variants/germlineSV.vcf"
cols <- colnames(read.table(pipe(paste0('grep -v "##" ', manta.sv.vcf,' | grep "#"| sed s/#//')), header = TRUE))
manta <- read.table(manta.sv.vcf, col.names = cols, stringsAsFactors = FALSE)
manta <- manta %>% filter(FILTER == "PASS")

get_sv_type <- function(x){
  if (grepl("BND", x)){
    # SVTYPE = BND
    root <- gsub(".*MATEID=", "", x)
    root <- gsub(":[01];.*", "", root)
    mate1 <- paste0(root, ":0")
    mate2 <- paste0(root, ":1")
    alt1 <- manta %>% filter(ID == mate1) %>% .$ALT
    alt2 <- manta %>% filter(ID == mate2) %>% .$ALT
    # Determine sv type based on breakpoint orientation
    if (gsub("(.*chr)|(\\:.*)", "", alt1) != gsub("(.*chr)|(\\:.*)", "", alt2)){
      sv.type <- "BND"
      
    } else if ((grepl("\\[", alt1) & grepl("\\[", alt2)) | (grepl("\\]", alt1) & grepl("\\]", alt2))){
      sv.type <- "INV"
      
    } else if (grepl("[A-Za-z.]\\[", alt1) & grepl("^\\]", alt2)){
      sv.type <- "DEL"
      
    } else if (grepl("^\\]", alt1) & grepl("[A-Za-z.]\\[", alt2)){
      sv.type <- "DUP/INS"
      
    } else{
      sv.type <- "BND"
    }
  } else{
    # SVTYPE = DEL/DUP/INS/INV
    sv.type <- gsub(".*SVTYPE=", "", x)
    sv.type <- gsub(";.*", "", sv.type)
  }
  return(sv.type)
}

manta$SVTYPE <- sapply(manta$INFO, get_sv_type)

# Extract deletions and inversions to a bed format file
gap <- manta %>% filter(SVTYPE == "DEL" | SVTYPE == "INV") %>% .[,c(1:5, 8)]
gap1 <- gap %>% filter(grepl("DEL", ID) | grepl("INV", ID))
gap1$END <- gsub("(^END=)|(;SVTYPE=.*)", "", gap1$INFO)
gap2 <- gap %>% filter(grepl("BND.*:0$", ID))
gap2$END <- gsub("(.*:)|(\\[[A-Za-z.]*)|(\\][A-Za-z.]*)", "", gap2$ALT)
gap <- rbind(gap1, gap2)
gap <- gap[,c(1, 2, 7, 3)]

for (i in 1:length(gap[,1])) {
  if(as.numeric(gap[i,2]) > as.numeric(gap[i,3])){
    tmp <- gap[i, 2]
    gap[i,2] <- gap[i, 3]
    gap[i,3] <- tmp
  }
}

# filter
# TODO dedup
# gap <- gap %>% filter(as.numeric(END) - as.numeric(POS) < 100000)
#
# gap <- gap[order(gap$CHROM, as.numeric(gap$POS)),]
#
# write.table(gap, paste0(manta.sv.vcf, ".bed"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# unfiltered
gap <- gap[order(gap$CHROM, as.numeric(gap$POS)),]

write.table(gap, paste0(manta.sv.vcf, ".unfiltered.bed"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

