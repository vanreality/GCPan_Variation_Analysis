library("dplyr")
library("ggplot2")
library("ggtree")
library("reshape2")
library("ComplexHeatmap")
library("circlize")
library("RCircos")

root <- "/lustre/home/acct-clswcc/clswcc-fsy/gastricCancer/panGenome"
setwd(paste0(root, "/script/plot"))

phenotype <- read.table(paste0(root, "/data/PAV/phenotype.txt"), header = TRUE, sep = "\t")
# z <- do.call(rbind,lapply(3:13, function(x){cbind(phen = colnames(phenotype)[x],as.data.frame(table(phenotype[[x]])))}))
namesN <- read.table(paste0(root, "/data/PAV/normal.individual.info"))
namesT <- read.table(paste0(root, "/data/PAV/tumor.individual.info"))
phenotype$NameN <- namesN %>% filter(V1 == phenotype$Normal) %>% .$V2
phenotype$NameT <- namesT %>% filter(V1 == phenotype$Tumor) %>% .$V2
cds <- read.table(paste0(root, "/data/PAV/cds.ann.bed"))
cds <- cds %>% filter(V1 == "chrY") %>% .[, c(1, 4, 5)]
cds <- unique(cds)

### SV PAV compare result
# 0 : SV    / Absence
# 1 : noSV  / Presence
# 2 : SV    / Presence    ** which indicates some variations detected by SV method but not PAV method
# 3 : noSV  / Absence     ** which indicates some variations detected by PAV method but not SV method

#<---------------------------------------------- SV PAV compare result plot ----------------------------------------------------------------->

plot_SV_PAV_compare <- function(x) {
        compare.result <- x %>% filter(!grepl("GC", Gene))
        compare.result <- compare.result[-which(compare.result$Gene %in% cds$V4),]
        data <- as.matrix(compare.result[,-1])[,phenotype$NameN]
        rownames(data) <- as.matrix(compare.result[,1])
        Heatmap(data, 
                col = structure(c("red", "white", "black", "blue"), names = c("0", "1", "2", "3")),
                name = "Compare result", 
                show_row_names = FALSE,
                show_column_names = FALSE,
                column_title = paste0(length(data[1,]), " samples"),
                row_title = paste0(length(data[,1]), " genes"),
                heatmap_legend_param = list(at = c(0, 1, 2, 3), labels = c("SV/Absence", "noSV/Presence", "SV/Presence", "noSV/Absence")),
                top_annotation = HeatmapAnnotation(gender = as.character(phenotype$Gender), 
                                                   age = phenotype$Age, 
                                                   diameter = phenotype$Diameter,
                                                   Lauren = as.character(phenotype$Lauren),
                                                   Borrman = as.character(phenotype$Borrman)))
}

#<----------- manta germline SV ----------->
compare.result <- read.csv(paste0(root, "/result/compareSVPAV/manta.germline.sv/coverage.0.8.gene.csv"))
plot_SV_PAV_compare(compare.result)

#<----------- manta somatic SV ----------->
compare.result <- read.csv(paste0(root, "/result/compareSVPAV/manta.somatic.sv/coverage.0.8.gene.csv"))
plot_SV_PAV_compare(compare.result)

#<----------- svaba germline SV ----------->
compare.result <- read.csv(paste0(root, "/result/compareSVPAV/svaba.germline.sv/coverage.0.8.gene.csv"))
plot_SV_PAV_compare(compare.result)

#<----------- svaba somatic SV ----------->
compare.result <- read.csv(paste0(root, "/result/compareSVPAV/svaba.somatic.sv/coverage.0.8.gene.csv"))
plot_SV_PAV_compare(compare.result)

#<-------------------------------------------------------- SV density and distribution plot ------------------------------------------------------------>

chrs <- paste0("chr", c(1:22,"X", "Y"))

get_manta_sv_type <- function(x){
        if (grepl("BND", x)){
                # SVTYPE = BND
                root <- gsub(".*MATEID=", "", x)
                root <- gsub(":[01];.*", "", root)
                mate1 <- paste0(root, ":0")
                mate2 <- paste0(root, ":1")
                alt1 <- vcf %>% filter(ID == mate1) %>% .$ALT
                alt2 <- vcf %>% filter(ID == mate2) %>% .$ALT
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

get_svaba_sv_type <- function(x){
        # Find mate pair
        root <- gsub(":[12]", "", x)
        mate1 <- paste0(root, ":1")
        mate2 <- paste0(root, ":2")
        alt1 <- vcf %>% filter(ID == mate1) %>% .$ALT
        alt2 <- vcf %>% filter(ID == mate2) %>% .$ALT
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

plot_SV_density <- function(x){
        tmp1 <- x %>% filter(!grepl("BND", INFO))
        tmp1$SVCHROM <- tmp1$CHROM
        tmp1$END <- gsub("(^END=)|(;SVTYPE=.*)", "", tmp1$INFO)
        tmp2 <- x %>% filter(grepl("BND", INFO))
        tmp2$SVCHROM <- gsub("(:.*)|([A-Za-z.]*\\[)|([A-Za-z.]*\\])", "", tmp2$ALT)
        tmp2$END <- gsub("(.*:)|(\\[[A-Za-z.]*)|(\\][A-Za-z.]*)", "", tmp2$ALT)
        vcf <- rbind(tmp1, tmp2)[, c(1, 2, 8, 9, 7)]
        vcf$END <- as.numeric(vcf$END)
        
        # density
        density.data <- vcf %>% filter((CHROM %in% chrs) & (CHROM == SVCHROM)) %>% .[, c(1, 2, 4, 5)]
        ggplot(density.data,aes(x = factor(CHROM, chrs), y = END - POS))+
                geom_violin()+
                geom_jitter(aes(color=SVTYPE), alpha = .8, size = .5)+
                scale_y_continuous(limits=c(-100,10000), trans = "log10")+
                theme_classic()+
                labs(x = "chromosome", y = "SV length")
}

plot_SV_circle <- function(x){
        tmp1 <- x %>% filter(!grepl("BND", INFO))
        tmp1$SVCHROM <- tmp1$CHROM
        tmp1$END <- gsub("(^END=)|(;SVTYPE=.*)", "", tmp1$INFO)
        tmp2 <- x %>% filter(grepl("BND", INFO))
        tmp2$SVCHROM <- gsub("(:.*)|([A-Za-z.]*\\[)|([A-Za-z.]*\\])", "", tmp2$ALT)
        tmp2$END <- gsub("(.*:)|(\\[[A-Za-z.]*)|(\\][A-Za-z.]*)", "", tmp2$ALT)
        vcf <- rbind(tmp1, tmp2)[, c(1, 2, 8, 9, 7)]
        vcf$END <- as.numeric(vcf$END)
        
        # circle
        data(UCSC.HG38.Human.CytoBandIdeogram)
        cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
        chr.exclude <- NULL
        tracks.inside <- 2
        tracks.outside <- 0
        RCircos.Set.Core.Components(cyto.info, chr.exclude,tracks.inside, tracks.outside)
        RCircos.List.Plot.Parameters()
        RCircos.Set.Plot.Area()     
        RCircos.Chromosome.Ideogram.Plot() 
        
        # track 1 sv number count
        histogram.data <- data.frame()
        step <- 10000000
        for (i in chrs) {
                max_length <- max(vcf %>% filter(CHROM == i) %>% .$POS)
                if(is.finite(max_length)) {
                        for (j in seq(0, max_length, step)) {
                                m <- if(j + step - 1 > max_length) max_length else j + step - 1
                                count <- vcf %>% filter((CHROM == i)&(POS > j)&(POS < m)) %>% nrow(.)
                                histogram.data <- rbind(histogram.data, cbind(i, j, m, count))
                        } 
                }
        }
        colnames(histogram.data) <- c("CHROM", "POS", "END", "COUNT")
        histogram.data$CHROM <- as.character(histogram.data$CHROM)
        histogram.data$POS <- as.numeric(as.character(histogram.data$POS))
        histogram.data$END <- as.numeric(as.character(histogram.data$END))
        histogram.data$COUNT <- as.numeric(as.character(histogram.data$COUNT))
        
        track.num <- 1
        data.col <- 4
        side <- "in"
        RCircos.Histogram.Plot(histogram.data, data.col, track.num, side)
        
        # track 2 link
        link <- vcf %>% filter((CHROM %in% chrs) & (SVCHROM %in% chrs) & (CHROM != SVCHROM))
        names(link)[3] <- "CHROM1"
        names(link)[4] <- "POS1"
        link$END <- link$POS
        link$END1 <- link$POS1
        link <- link[, c(1, 2, 6, 3, 4, 7)]
        
        track.num <- 2
        RCircos.Link.Plot(link, track.num, TRUE)
}

filter_vcf <- function(x) {
        cols <- colnames(read.table(pipe(paste0('grep -v "##" ', x,' | grep "#"| sed s/#//')), header = TRUE))
        cols <- sapply(cols, function(x) gsub("(.*GATK_bam\\.)|(\\.recal.*)", "", x))
        vcf <- read.table(x, col.names = cols, stringsAsFactors = FALSE)
        
        
        #vcf <- vcf %>% filter(FILTER == "PASS") %>% .[, c(1:5, 8)]
        
        
        return(vcf)
}

#<----------- manta germline SV ----------->
vcf.path <- paste0(root, "/result/manta/WGC012904D.WGC013444D/result/results/variants/germlineSV.vcf")
vcf <- filter_vcf(vcf.path)
vcf$SVTYPE <- sapply(vcf$INFO, get_manta_sv_type)
plot_SV_density(vcf)
plot_SV_circle(vcf)

#<----------- manta somatic SV ----------->
vcf.path <- paste0(root, "/result/manta/WGC012904D.WGC013444D/result/results/variants/somaticSV.vcf")
vcf <- filter_vcf(vcf.path)
vcf$SVTYPE <- sapply(vcf$INFO, get_manta_sv_type)
plot_SV_density(vcf)
plot_SV_circle(vcf)

#<----------- svaba germline SV ----------->
vcf.path <- paste0(root, "/result/svaba/WGC012904D.WGC013444D/WGC012904D.WGC013444D.svaba.germline.sv.vcf")
vcf <- filter_vcf(vcf.path)
vcf$SVTYPE <- sapply(vcf$ID, get_svaba_sv_type)
plot_SV_density(vcf)
plot_SV_circle(vcf)

#<----------- svaba somatic SV ----------->
vcf.path <- paste0(root, "/result/svaba/WGC012904D.WGC013444D/WGC012904D.WGC013444D.svaba.somatic.sv.vcf")
vcf <- filter_vcf(vcf.path)
vcf$SVTYPE <- sapply(vcf$ID, get_svaba_sv_type)
plot_SV_density(vcf)
plot_SV_circle(vcf)











