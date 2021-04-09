library("dplyr")

args <- commandArgs(T)

# Read in svaba vcf file
svaba.sv.vcf <- args[1]
# svaba.sv.vcf <- "/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/result/svaba/WGC012904D.WGC013444D/WGC012904D.WGC013444D.svaba.germline.sv.vcf"
cols <- colnames(read.table(pipe(paste0('grep -v "##" ', svaba.sv.vcf,' | grep "#"| sed s/#//')), header = TRUE))
cols <- sapply(cols, function(x) gsub("(.*GATK_bam\\.)|(\\.recal.*)", "", x))
svaba <- read.table(svaba.sv.vcf, col.names = cols, stringsAsFactors = FALSE)

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
svaba <- svaba %>% filter((SVTYPE == "BND") | (gsub(".*:", "", ID) == "1"))

for (i in 1:nrow(svaba)) {
  svaba[i, "INFO"] <- gsub("BND", svaba[i, "SVTYPE"], svaba[i, "INFO"])
  if (svaba[i, "SVTYPE"] != "BND") {
    svaba[i, "INFO"] <- paste0(svaba[i, "INFO"], ";END=", gsub("(.*:)|(\\].)|(\\[.)", "", svaba[i, "ALT"]))
  }
}

cmd <- paste0("cat ", svaba.sv.vcf, " | grep '^#' > ", svaba.sv.vcf, ".converted")
system(cmd)
write.table(svaba, 
            paste0(svaba.sv.vcf, ".converted"), 
            append = TRUE, 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE, 
            sep = "\t")



