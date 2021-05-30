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
gene.bed <- read.table(paste0(root, "/data/PAV/gene.bed"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)

### SV PAV compare result
# 0 : SV    / Absence
# 1 : noSV  / Presence
# 2 : SV    / Presence    ** which indicates some variations detected by SV method but not PAV method
# 3 : noSV  / Absence     ** which indicates some variations detected by PAV method but not SV method

#<---------------------------------------------- SV PAV compare result plot ----------------------------------------------------------------->

phenotype_asso_genes <- c("ANKRD9",  "BHLHA9", "BPY2", "BPY2B", "BPY2C", "CDY1", "CDY1B", 
                          "EIF1AY","GPR42", "GRIN2D",  "HCN2",  "HIST1H4B",  "HLA-DRB1", 
                          "HLA-DRB5", "MAFA",  "MAP9",  "OR4C11",  "OR4P4",  "OR4S2", "PIM3",
                          "PSG4", "PSG9",  "SRY",  "TAF4",  "TAS2R43", "TCF15",
                          "TRIM48",  "TRIM64",  "UGT2B28")

SV_PAV_compare_stat <- function(x,y,z){
        compare.result <- x %>% filter(!grepl("GC", Gene))
        compare.result <- compare.result[-which(compare.result$Gene %in% cds$V4),]
        info <- gene.bed[which(gene.bed$V4 %in% compare.result$Gene), c(1, 6)]
        data <- as.matrix(compare.result[,-1])[, z]
        rownames(data) <- as.matrix(info[,2])
        return(data)
}

plot_SV_PAV_compare <- function(x, y, z) {
        data <- SV_PAV_compare_stat(x, z)
        mark_logi <- rownames(data) %in% phenotype_asso_genes
        ht <- Heatmap(data, 
                col = structure(c("#D6CD56", "#BB4F50", "#74BBE8", "#3C71A0"), names = c("0", "1", "2", "3")),
                left_annotation = rowAnnotation(a = anno_simple(as.character(rownames(data) %in% inter.genes), 
                                                                width = unit(2, 'mm'),
                                                                col = structure(c("white", "#E68C4C"), names = c("FALSE", "TRUE") )),
                                                show_annotation_name = F),
                right_annotation = rowAnnotation(c = anno_mark(at = (1:nrow(data))[mark_logi],
                                                               labels = rownames(data)[mark_logi],
                                                               link_width = unit(10, "mm"))),
                name = "Comparison", 
                show_row_names = FALSE,
                show_column_names = FALSE,
                column_dend_side = "bottom",
                column_title_side = "bottom",
                column_title = paste0(y, " (n=",length(data[1,]), ")"),
                row_title = paste0(length(data[,1]), " genes"),
                heatmap_legend_param = list(at = c(0, 1, 2, 3), labels = c("SV(-)/PAV(-)", "SV(+)/PAV(+)", "SV(-)/PAV(+)", "SV(+)/PAV(-)")))
        draw(ht, heatmap_legend_side = "right", 
             merge_legends = T,
             annotation_legend_list = Legend(labels= c("Shared genes", ""),
                                                  legend_gp = grid::gpar(fill = structure(c( "#E68C4C", "white"), names = c("TRUE", "FALSE") ))))
}

#<----------- data process --------------------->
compare.res.normal <- read.csv(paste0(root, "/result/compareSVPAV/survivor.germline.sv/coverage.0.8.gene.csv"), stringsAsFactors = FALSE)
compare.res.normal.stat <- SV_PAV_compare_stat(compare.res.normal, phenotype$NameN)

compare.res.tumor <- read.csv(paste0(root, "/result/compareSVPAV/survivor.somatic.sv/coverage.0.8.gene.csv"), stringsAsFactors = FALSE)
compare.res.tumor.stat <- SV_PAV_compare_stat(compare.res.tumor, phenotype$NameT)

inter.genes <- intersect(rownames(compare.res.normal.stat), rownames(compare.res.tumor.stat))

#<----------- survivor germline SV ----------->
title <- "Normal samples"

pdf("DEL_N.pdf", width = 11.5, height = 9)
plot_SV_PAV_compare(compare.res.normal, title, phenotype$NameN)
dev.off()

# pdf("DEL_N.pdf", width = 11, height =8.5)
# plot_SV_PAV_compare(compare.result[apply(compare.result, 1, function(x) {sum(x==3)})>54,], title, phenotype$NameN)
# dev.off()

#<----------- survivor somatic SV ----------->
title <- "Tumor samples"

pdf("DEL_T.pdf", width = 11.5, height = 9)
plot_SV_PAV_compare(compare.res.tumor, title, phenotype$NameT)
dev.off()

# pdf("DEL_T.pdf", width = 11, height =8.5)
# plot_SV_PAV_compare(compare.result[apply(compare.result, 1, function(x) {sum(x==3)})>54,], title, phenotype$NameT)
# dev.off()

# 
# #<----------- manta germline SV ----------->
# compare.result <- read.csv(paste0(root, "/result/compareSVPAV/manta.germline.sv/coverage.0.8.gene.csv"))
# plot_SV_PAV_compare(compare.result)
# 
# compare.result <- read.csv(paste0(root, "/result/compareSVPAV/manta.germline.sv/coverage.unfiltered.0.8.gene.csv"))
# plot_SV_PAV_compare(compare.result)
# 
# #<----------- manta somatic SV ----------->
# compare.result <- read.csv(paste0(root, "/result/compareSVPAV/manta.somatic.sv/coverage.0.8.gene.csv"))
# plot_SV_PAV_compare(compare.result)
# 
# compare.result <- read.csv(paste0(root, "/result/compareSVPAV/manta.somatic.sv/coverage.unfiltered.0.8.gene.csv"))
# plot_SV_PAV_compare(compare.result)
# 
# #<----------- svaba germline SV ----------->
# compare.result <- read.csv(paste0(root, "/result/compareSVPAV/svaba.germline.sv/coverage.0.8.gene.csv"))
# plot_SV_PAV_compare(compare.result)
# 
# compare.result <- read.csv(paste0(root, "/result/compareSVPAV/svaba.germline.sv/coverage.unfiltered.0.8.gene.csv"))
# plot_SV_PAV_compare(compare.result)
# 
# #<----------- svaba somatic SV ----------->
# compare.result <- read.csv(paste0(root, "/result/compareSVPAV/svaba.somatic.sv/coverage.0.8.gene.csv"))
# plot_SV_PAV_compare(compare.result)
# 
# compare.result <- read.csv(paste0(root, "/result/compareSVPAV/svaba.somatic.sv/coverage.unfiltered.0.8.gene.csv"))
# plot_SV_PAV_compare(compare.result)

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






