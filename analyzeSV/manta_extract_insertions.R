library("dplyr")

args <- commandArgs(T)

# Read in manta vcf file
manta.sv.vcf <- args[1]
# manta.sv.vcf <- "/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome/result/manta/WGC012906D.WGC013446D/result/results/variants/somaticSV.vcf"
cols <- colnames(read.table(pipe(paste0('grep -v "##" ', manta.sv.vcf,' | grep "#"| sed s/#//')), header = TRUE))
manta <- read.table(manta.sv.vcf, col.names = cols, stringsAsFactors = FALSE)
manta <- manta %>% filter(FILTER == "PASS")

ins <- manta %>% filter(grepl("SVTYPE=INS", INFO))
if (nrow(ins) != 0) {
  ins <- do.call(rbind,lapply(1:nrow(ins), function(x){
    if (ins[x, "ALT"] == "<INS>") {
      if (grepl("LEFT", ins[x, "INFO"])) {
        left <- gsub(".*LEFT_SVINSSEQ=", "", ins[x, "INFO"])
        left <- gsub(";.*", "", left)
        left <- cbind(ins[x, 1:2], ID = paste0(ins[x, 3], "_LEFT"), SEQ = left)
        
        right <- gsub(".*RIGHT_SVINSSEQ=", "", ins[x, "INFO"])
        right <- gsub(";.*", "", right)
        right <- cbind(ins[x, 1:2], ID = paste0(ins[x, 3], "_RIGHT"), SEQ = right)
        
        data <- data.frame(rbind(left, right), stringsAsFactors = FALSE)
        colnames(data) <- c("CHROM", "POS", "ID", "SEQ")
        return(data)
      } else {
        sequence <- gsub(".*SVINSSEQ=", "", ins[x, "INFO"])
        sequence <- gsub(";.*", "", sequence)
        
        data <- data.frame(cbind(ins[x, 1:3], sequence), stringsAsFactors = FALSE)
        colnames(data) <- c("CHROM", "POS", "ID", "SEQ")
        return(data)
      }
    } else {
      data <- data.frame(cbind(ins[x, c(1:3, 5)]), stringsAsFactors = FALSE)
      colnames(data) <- c("CHROM", "POS", "ID", "SEQ")
      return(data)
    }
  }))
  
  fasta <- do.call(rbind, lapply(1:nrow(ins), function(x) {
    INFO <- paste0(">", ins[x, "ID"], "  ", ins[x, "CHROM"], "  ", ins[x, "POS"])
    SEQ <- as.vector(ins[x, "SEQ"])
    return(rbind(INFO, SEQ))
  }))
  
  write.table(fasta, paste0(manta.sv.vcf, ".fasta"), row.names = FALSE, col.names = FALSE, quote = FALSE)
}

