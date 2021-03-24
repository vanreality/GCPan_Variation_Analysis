library("dplyr")

args <- commandArgs(T)

# Read in svaba vcf file
svaba.sv.vcf <- args[1]
# svaba.sv.vcf <- "/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/result/svaba/WGC012904D.WGC013444D/WGC012904D.WGC013444D.svaba.germline.sv.vcf"
cols <- colnames(read.table(pipe(paste0('grep -v "##" ', svaba.sv.vcf,' | grep "#"| sed s/#//')), header = TRUE))
cols <- sapply(cols, function(x) gsub("(.*GATK_bam\\.)|(\\.recal.*)", "", x))
svaba <- read.table(svaba.sv.vcf, col.names = cols, stringsAsFactors = FALSE)
svaba <- svaba %>% filter(FILTER == "PASS")

get_sv_type <- function(x){
  # Find mate pair
  root <- gsub(":[12]", "", x)
  mate1 <- paste0(root, ":1")
  mate2 <- paste0(root, ":2")
  alt1 <- svaba %>% filter(ID == mate1) %>% .$ALT
  alt2 <- svaba %>% filter(ID == mate2) %>% .$ALT
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
  return(sv.type)
}

svaba$SVTYPE <- sapply(svaba$ID, get_sv_type)

# Extract deletions and inversions to a bed format file
gap <- svaba %>% filter((SVTYPE == "DEL" | SVTYPE == "INV") & grepl(":1", ID)) %>% .[,1:5]
gap$END <- gsub("(.*:)|(\\[[A-Za-z.]*)|(\\][A-Za-z.]*)", "", gap$ALT)
gap <- gap[,c(1, 2, 6, 3)]

# filter SV length < 100,000
# gap <- gap %>% filter(as.numeric(END) - as.numeric(POS) < 100,000)
#
# gap <- gap[order(gap$CHROM, gap$POS),]
#
# write.table(gap, paste0(svaba.sv.vcf, ".bed"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# unfiltered SV
gap <- gap[order(gap$CHROM, gap$POS),]

write.table(gap, paste0(svaba.sv.vcf, ".unfiltered.bed"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


# TODO rewrite bedtools intersect and coverage calculate python script

# cds.bed <- "/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/data/gff/cds.ann.bed"
# cds <- read.table(cds.bed, col.names = cbind("CHROM", "POS", "END", "ENSG", "GENE"))
# cds <- arrange(cds, CHROM, POS)

# intersect_sv_cds <- function(sv, cds) {
#   if (sv$CHROM == cds$CHROM){
    
#   }
# }

# cov <- intersect_sv_cds(gap, cds)



