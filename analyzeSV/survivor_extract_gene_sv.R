library("dplyr")

args <- commandArgs(T)

# Read in survivor merged vcf file
survivor.sv.vcf <- args[1]
control <- args[2]

# survivor.sv.vcf <- "/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/result/survivor/WGC012904D.WGC013444D/WGC012904D.WGC013444D.germline.merged.vcf"
# control <- "germline"

cols <- colnames(read.table(pipe(paste0('grep -v "##" ', survivor.sv.vcf,' | grep "#"| sed s/#//')), header = TRUE))
survivor <- read.table(survivor.sv.vcf, col.names = cols, stringsAsFactors = FALSE)
survivor <- survivor %>% filter(FILTER == "PASS")

get_sv_type <- function(x){
  # SVTYPE = DEL/DUP/INS/INV
  sv.type <- gsub(".*SVTYPE=", "", x)
  sv.type <- gsub(";.*", "", sv.type)
  return(sv.type)
}

survivor$SVTYPE <- sapply(survivor$INFO, get_sv_type)

# Extract deletions and inversions (GT : 1/1) to a bed format file
if (control == "germline") {
  survivor <- survivor[which(gsub(":.*", "", survivor[,ncol(survivor) - 2]) == "1/1"),]
} else {
  survivor <- survivor[which(gsub(":.*", "", survivor[,ncol(survivor) - 1]) == "1/1"),]
}
gap <- survivor %>% filter(SVTYPE == "DEL" | SVTYPE == "INV") %>% .[,c(1:5, 8)]
gap$END <- gsub("(.*;END=)|(;CIPOS=.*)", "", gap$INFO)
gap <- gap[,c("CHROM", "POS", "END", "ID")]


# for (i in 1:length(gap[,1])) {
#   if(as.numeric(gap[i,2]) > as.numeric(gap[i,3])){
#     tmp <- gap[i, 2]
#     gap[i,2] <- gap[i, 3]
#     gap[i,3] <- tmp
#   }
# }

# gap <- gap[order(gap$CHROM, as.numeric(gap$POS)),]

write.table(gap, paste0(survivor.sv.vcf, ".bed"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

