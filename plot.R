library(getopt)
opt_spec <- matrix(c('help', 'h', 0, 'logical', 'help manual',
                     'input', 'i', 1, 'character', 'filepath for *.BAF.txt',
                     'ratioGraph', 'r', 0, 'logical', 'if RatioGraph should be drawn',
                     'BAFGraph', 'f', 0, 'logical', 'if BAFGraph should be drawn',
                     'output', 'o', 1, 'character', 'filepath for output figures',
                     'bed', 'b', 1, 'character', 'panel bed file'),
                   byrow=TRUE, ncol=5)
opt = getopt(opt_spec, commandArgs(TRUE))
if (!is.null(opt$help)){
  cat(getopt(opt_spec, usage=TRUE))
  q(save='no', status=1)
}
stopifnot(file.exists(opt$input, opt$bed))

dir.create(opt$output)
if(is.null(opt$BAFGraph) == F){
  baflist <- list.files(path = opt$input, pattern = "BAF.txt$")
  bed <- read.delim(opt$bed,  header = F, stringsAsFactors = F)
  colnames(bed) <- c("chrom", "start", "end", "gene")
  lim <- sum(bed$end - bed$start)
  for(v in 1:length(baflist)){
    cat(paste(baflist[v], "\n"))
	patient <- unlist(strsplit(baflist[v], split = "_"))[1]
    baf <- read.delim(paste(opt$input, baflist[v], sep = "/"), header = T, stringsAsFactors = F)
    
    
    chr <- c(1:22, "X", "Y")
    bed_new <-  NULL
    chr_info <- data.frame(chr = paste("chr", chr, sep = ""), mid = NA, stringsAsFactors = F)
    baf_new <-  NULL
    for(i in 1:24){
      baf_chr <- subset(baf, Chromosome == chr[i])
      bed_chr <- subset(bed, chrom == paste("chr", chr[i], sep = ""))
      if(i == 1){
        for(j in 1:nrow(bed_chr)){
          bed_chr$end_norm[j] <- sum(bed_chr$end[1:j] - bed_chr$start[1:j])
          bed_chr$start_norm[j] <- bed_chr$end_norm[j] - (bed_chr$end[j] - bed_chr$start[j])
        }
      }else{
        for(j in 1:nrow(bed_chr)){
          bed_chr$end_norm[j] <- sum(bed_chr$end[1:j] - bed_chr$start[1:j]) + tmp
          bed_chr$start_norm[j] <- bed_chr$end_norm[j] - (bed_chr$end[j] - bed_chr$start[j])
        }
      }
      chr_info$mid[i] <- (bed_chr$start_norm[1] + bed_chr$end_norm[nrow(bed_chr)])/2
      for(j in 1:nrow(baf_chr)){
        n <- which(bed_chr$start < baf_chr$Position[j] & bed_chr$end >= baf_chr$Position[j])
        baf_chr$pos_norm[j] <- baf_chr$Position[j] - bed_chr$start[n] + bed_chr$start_norm[n]
      }
      if(is.null(bed_new)){
        bed_new <- bed_chr
      }else{
        bed_new <- rbind(bed_new, bed_chr)
      }
      if(is.null(baf_new)){
        baf_new <- baf_chr
      }else{
        baf_new <- rbind(baf_new, baf_chr)
      }
      tmp <- bed_chr$end_norm[nrow(bed_chr)]
    }
    
    tiff(filename = paste(opt$output, "/", gsub("txt", "tiff", baflist[v]), sep = ""), width = 1500, height = 400,
         units = "px", pointsize = 14, bg = "white", res = NA)
    par(mar = c(5,5.5,4,2))
    plot(baf_new$pos_norm,baf_new$BAF, xlim = c(0,lim), ylim = c(-0.1,1.1), xlab = "Chromosome", cex.main = 2, cex.axis = 1.5, cex.lab = 1.8, las = 2,
         ylab = "BAF", main = patient, pch = ".", col = colors()[1], xaxt = "n")
    #axis(1, labels = chr_info$chr, at = chr_info$mid, tck = 0, las = 2, cex.axis = 1.5, hadj = 0.75)
	axis(1, labels = gsub("chr", "", chr_info$chr), at = chr_info$mid, tck = 0, cex.axis = 1.5, hadj = 0.75, las = 2)
    #axis(2, tck = 0, las = 2, cex = 1.2, hadj = 0.75)
    tt <- which(baf_new$A == 0.5)		
    points(baf_new$pos_norm[tt], baf_new$BAF[tt], pch = 20, cex = 0.6, col = colors()[92])
    tt <- which(baf_new$A != 0.5 & baf_new$A >= 0)
    points(baf_new$pos_norm[tt], baf_new$BAF[tt], pch = 20, cex = 0.6, col = "royalblue2")
    tt <- 1
    pres <- 1
    
    if (length(baf_new$A) > 4) {
      for (j in c(2:(length(baf_new$A) - pres - 1))) {
        if (baf_new$A[j] == baf_new$A[j+pres]) {	
          tt[length(tt)+1] <- j 
        }
      }
      tt <- intersect(tt, which((baf_new$A == 0.5 | baf_new$A == 1 | baf_new$A == 0) & baf_new$A >= 0))
      points(baf_new$pos_norm[tt], baf_new$A[tt], pch = ".", col = colors()[24], cex=4)
      points(baf_new$pos_norm[tt], baf_new$B[tt], pch = ".", col = colors()[24], cex=4)	
    }
    if (length(baf_new$A) > 4) {
      for (j in c(2:(length(baf_new$A) - pres - 1))) {
        if (baf_new$A[j] == baf_new$A[j+pres]) {	
          tt[length(tt)+1] <- j 
        }
      }
      tt <- intersect(tt, which(baf_new$A != 0.5 & baf_new$A != 0 & baf_new$A != 1 & baf_new$A >= 0))
      points(baf_new$pos_norm[tt], baf_new$A[tt], pch = 20, col = colors()[453], cex=0.4)
      points(baf_new$pos_norm[tt], baf_new$B[tt], pch = 20, col = colors()[453], cex=0.4)	
    }
    
    dev.off()
  }
}


if(is.null(opt$ratioGraph) == F){
  ratiolist <- list.files(path = opt$input, pattern = "ratio.txt$")
  cnvlist <- list.files(path = opt$input, pattern = "CNVs$")
  bed <- read.delim(opt$bed,  header = F, stringsAsFactors = F)
  colnames(bed) <- c("chrom", "start", "end", "gene")
  for(v in 1:length(ratiolist)){
    cat(paste(ratiolist[v], "\n"))
    patient <- unlist(strsplit(ratiolist[v], split = "_"))[1]
    ratio <- read.delim(paste(opt$input, ratiolist[v], sep = "/"), header = T, stringsAsFactors = F)
    ratio$Ratio[which(ratio$Ratio <= 0.01)] <- 0.01
    ratio$End <- as.integer(unlist(strsplit(ratio$Gene, split = "-"))[seq(2, nrow(ratio)*2, 2)])
    ratio$Start <- ratio$Start - 1
    cnv <- read.delim(paste(opt$input, cnvlist[v], sep = "/"), header = F, stringsAsFactors = F)
    colnames(cnv)[1:4] <- c("Chromosome", "Start", "End", "CopyNumber")
    
    chr <- c(1:22, "X", "Y")
    bed_new <-  NULL
    chr_info <- data.frame(chr = paste("chr", chr, sep = ""), mid = NA, stringsAsFactors = F)
    ratio_new <-  NULL
    cnv_new <-  NULL
    for(i in 1:24){
      ratio_chr <- subset(ratio, Chromosome == chr[i])
      bed_chr <- subset(bed, chrom == paste("chr", chr[i], sep = ""))
      cnv_chr <- subset(cnv, Chromosome == chr[i])
      if(i == 1){
        for(j in 1:nrow(bed_chr)){
          bed_chr$end_norm[j] <- sum(bed_chr$end[1:j] - bed_chr$start[1:j])
          bed_chr$start_norm[j] <- bed_chr$end_norm[j] - (bed_chr$end[j] - bed_chr$start[j])
        }
      }else{
        for(j in 1:nrow(bed_chr)){
          bed_chr$end_norm[j] <- sum(bed_chr$end[1:j] - bed_chr$start[1:j]) + tmp
          bed_chr$start_norm[j] <- bed_chr$end_norm[j] - (bed_chr$end[j] - bed_chr$start[j])
        }
      }
      chr_info$mid[i] <- (bed_chr$start_norm[1] + bed_chr$end_norm[nrow(bed_chr)])/2
      for(j in 1:nrow(ratio_chr)){
        n <- which(bed_chr$start <= ratio_chr$Start[j] & bed_chr$end >= ratio_chr$Start[j])
        ratio_chr$start_norm[j] <- ratio_chr$Start[j] - bed_chr$start[n] + bed_chr$start_norm[n]
        ratio_chr$end_norm[j] <- ratio_chr$End[j] - bed_chr$start[n] + bed_chr$start_norm[n]
      }
      for(j in 1:nrow(cnv_chr)){
        n <- which(bed_chr$start <= cnv_chr$Start[j] & bed_chr$end >= cnv_chr$Start[j])
        cnv_chr$start_norm[j] <- cnv_chr$Start[j] - bed_chr$start[n] + bed_chr$start_norm[n]
        n <- which(bed_chr$start <= cnv_chr$End[j] & bed_chr$end >= cnv_chr$End[j])
        cnv_chr$end_norm[j] <- cnv_chr$End[j] - bed_chr$start[n] + bed_chr$start_norm[n]
      }
      if(is.null(bed_new)){
        bed_new <- bed_chr
      }else{
        bed_new <- rbind(bed_new, bed_chr)
      }
      if(is.null(ratio_new)){
        ratio_new <- ratio_chr
      }else{
        ratio_new <- rbind(ratio_new, ratio_chr)
      }
      if(is.null(cnv_new)){
        cnv_new <- cnv_chr
      }else{
        cnv_new <- rbind(cnv_new, cnv_chr)
      }
      tmp <- bed_chr$end_norm[nrow(bed_chr)]
    }
    
    tiff(filename = paste(opt$output, "/", gsub("txt", "tiff", ratiolist[v]), sep = ""), width = 1500, height = 400,
         units = "px", pointsize = 14, bg = "white", res = NA)
    par(mar = c(5,5.5,4,2))
	plot(ratio_new$start_norm,ratio_new$Ratio, ylim = c(log2(min(ratio_new$Ratio)), log2(max(ratio_new$Ratio))), 
         ylab = "Normalized Copy Number (log2)", xlab = "Chromosome", cex.main = 2, cex.axis = 1.5, cex.lab = 1.8, las = 2,
         main = patient, pch = ".", col = colors()[1], xaxt = "n")
    #axis(1, labels = chr_info$chr, at = chr_info$mid, tck = 0, las = 2, cex = 1.5, hadj = 0.75)
    axis(1, labels = gsub("chr", "", chr_info$chr), at = chr_info$mid, tck = 0, cex.axis = 1.5, hadj = 0.75, las = 2)
	tt <- which(ratio_new$CopyNumber == 2 )
    if(length(tt) != 0){
      for(j in 1:length(tt)){
        aa <- ratio_new$start_norm[tt[j]]:ratio_new$end_norm[tt[j]]
        #points(aa,rep(log2(ratio_new$Ratio[tt[j]]), length(aa)),pch = 20, cex = 0.1,col = colors()[88])
		points(aa,rep(log2(ratio_new$Ratio[tt[j]]), length(aa)),pch = 20, cex = 0.4,col = colors()[88])
      }
    }
    tt <- which(ratio_new$CopyNumber > 2 )
    if(length(tt) != 0){
      for(j in 1:length(tt)){
        aa <- ratio_new$start_norm[tt[j]]:ratio_new$end_norm[tt[j]]
        #points(aa,rep(log2(ratio_new$Ratio[tt[j]]), length(aa)),pch = 20, cex = 0.1,col = colors()[136])
		points(aa,rep(log2(ratio_new$Ratio[tt[j]]), length(aa)),pch = 20, cex = 0.4,col = colors()[136])
      }
    }
    tt <- which(ratio_new$CopyNumber < 2 & ratio_new$CopyNumber!= -1)
    if(length(tt) != 0){
      for(j in 1:length(tt)){
        aa <- ratio_new$start_norm[tt[j]]:ratio_new$end_norm[tt[j]]
        #points(aa,rep(log2(ratio_new$Ratio[tt[j]]), length(aa)),pch = 0, cex = 0.1,col = colors()[461])
		points(aa,rep(log2(ratio_new$Ratio[tt[j]]), length(aa)),pch = 20, cex = 0.4,col = colors()[461])
      }
    }
    tt <- which(cnv_new$CopyNumber == 2 )
    if(length(tt) != 0){
      for(j in 1:length(tt)){
        aa <- cnv_new$start_norm[tt[j]]:cnv_new$end_norm[tt[j]]
        points(aa,rep(log2(cnv_new$CopyNumber[tt[j]]/2), length(aa)),col = colors()[24], pch = ".", cex = 4)
      }
    }
    tt <- which(cnv_new$CopyNumber > 2 )
    if(length(tt) != 0){
      for(j in 1:length(tt)){
        aa <- cnv_new$start_norm[tt[j]]:cnv_new$end_norm[tt[j]]
        points(aa,rep(log2(cnv_new$CopyNumber[tt[j]]/2), length(aa)),col = colors()[24], pch = ".", cex = 4)
      }
    }
    tt <- which(cnv_new$CopyNumber < 2 )
    if(length(tt) != 0){
      for(j in 1:length(tt)){
        aa <- cnv_new$start_norm[tt[j]]:cnv_new$end_norm[tt[j]]
        points(aa,rep(log2(cnv_new$CopyNumber[tt[j]]/2), length(aa)),col = colors()[24], pch = ".", cex = 4)
      }
    }
    dev.off()
  }
}
