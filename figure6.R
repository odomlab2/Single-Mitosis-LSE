############################
#####     Figure 6     #####  
############################
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
require(changepoint)
require(pheatmap)
require(viridis)
require(scales)

# Chromosome size vector
txdb <- keepStandardChromosomes(TxDb.Mmusculus.UCSC.mm10.knownGene)
ch.vec <- seqlengths(txdb); ch.vec <- ch.vec[!(names(ch.vec) %in% c('chrY', 'chrM'))]

# Read in mutations
snv.new <- 'Path to PF1_SisterMutations.rds'
all.mutations <- readRDS(snv.new); all.mutations <- all.mutations[all.mutations$FILTER]

# Split into single samples
rds.sis <- split(all.mutations, all.mutations$SampleName)
rds.sis <- rds.sis[grep('Neg', names(rds.sis), invert=TRUE)]

# Annotate T0 SNVs
snp.match <- seq(1, length(rds.sis), by=2)
for(i in snp.match){
  sister1 <- rds.sis[[i]]; sister2 <- rds.sis[[(i+1)]]
  elementMetadata(rds.sis[[i]])$TimePoint <- rep('T1', length(rds.sis[[i]]))
  elementMetadata(rds.sis[[(i+1)]])$TimePoint <- rep('T1', length(rds.sis[[(i+1)]]))
  rds.sis[[i]]$TimePoint[queryHits(findOverlaps(sister1, sister2))] <- 'T0'
  rds.sis[[(i+1)]]$TimePoint[subjectHits(findOverlaps(sister1, sister2))] <- 'T0'
}

## Segmentation (using information from SNP pairs)
##------------------------------------------------
# Sister matching index
snp.match <- seq(1, length(rds.sis), by=2)

# Reverse complement vector 
rev.comp.vec <- c('A', 'C', 'T', 'G'); names(rev.comp.vec) <- c('T', 'G', 'A', 'C')

# Changepoint segementation, UV
seg.maps <- lapply(snp.match, function(i){
  
  message(names(rds.sis)[i])
  # Get granges from both sisters, select only UV SNPs
  sister.1 <- rds.sis[[i]]; sister.2 <- rds.sis[[(i + 1)]]
  sister.1 <- sister.1[(sister.1$REF == 'C' & sister.1$ALT == 'T') | (sister.1$REF == 'G' & sister.1$ALT == 'A') & sister.1$TimePoint=='T1']
  sister.2 <- sister.2[(sister.2$REF == 'C' & sister.2$ALT == 'T') | (sister.2$REF == 'G' & sister.2$ALT == 'A') & sister.2$TimePoint=='T1']
  
  # For using rev comp in segmentation
  sister.2.swap <- sister.2 
  
  # Reverse complement sister 2 for segmentation
  elementMetadata(sister.2.swap)$REF <- rev.comp.vec[elementMetadata(sister.2)$REF]
  elementMetadata(sister.2.swap)$ALT <- rev.comp.vec[elementMetadata(sister.2)$ALT]
  sisters.swap <- unlist(as(list(sister.1, sister.2.swap), 'GRangesList'))
  
  # Order combined sister GRanges
  sisters.swap.order <- order(as.numeric(gsub('X', 20, gsub('chr','', as.character(seqnames(sisters.swap))))), as.numeric(start(sisters.swap)), decreasing = FALSE)
  sisters.swap <- sisters.swap[sisters.swap.order]
  
  # Segment using both sister data  
  segment.types <- do.call(rbind, lapply(1:length(ch.vec), function(iii){
    
    rel.snps <- sisters.swap[as.character(seqnames(sisters.swap))==names(ch.vec)[iii]]
    is.c <- rel.snps$REF=='C'; is.g <- rel.snps$REF=='G'
    
    # Not enough SNPs to segment
    if(length(rel.snps) < 10){
      results <- data.frame('start_gen'=sum(ch.vec[1:(iii-1)]),
                            'start_chr'=1,
                            'end_gen'=sum(ch.vec[1:iii]),
                            'end_chr'=ch.vec[[iii]],
                            'skew'=NA,
                            'SCE_site_est'=NA,
                            'SCE_site_uncertainty'=NA,
                            'chrom'=names(ch.vec)[iii])
      return(results)
      next
    }
    
    # Create vector of -1 to 1
    bin.vec <- rep(0, length(rel.snps)); bin.vec[rel.snps$REF=='C'] <- 1
    tmp.dist <- values(distanceToNearest(rel.snps))$distance
    
    # Calculate fit points
    cut.points <- cpt.mean(bin.vec, penalty='AIC', method='PELT', class=FALSE,  Q = 12)
    
    # Get SCE site and ambiguity
    sce_site <- c(); sce_ambiguity <- c()
    if(length(cut.points)==1){
      sce_site <- NA; sce_ambiguity <- NA
    } else {
      for(iiii in 1:(length(cut.points))){
        if(iiii ==  length(cut.points)){ 
          sce_site <- c(sce_site, NA); sce_ambiguity <- c(sce_ambiguity, NA)
        } else {
          sce_site <- c(sce_site, as.numeric(start(rel.snps)[cut.points[iiii]]) - (as.numeric(start(rel.snps)[cut.points[iiii]]) - as.numeric(start(rel.snps)[(cut.points[iiii] - 1)]))/2)
          sce_ambiguity <- c(sce_ambiguity, (as.numeric(start(rel.snps)[cut.points[iiii]]) - as.numeric(start(rel.snps)[(cut.points[iiii] - 1)])))
        }
        
      }
    }
    
    # Skew calculations
    # For Sister 1
    sis1.snps <- sister.1[as.character(seqnames(sister.1))==names(ch.vec)[iii]]
    sis1.vec <- rep(0, length(sis1.snps)); sis1.vec[sis1.snps$REF=='C'] <- 1
    
    # For Sister 2
    sis2.snps <- sister.2[as.character(seqnames(sister.2))==names(ch.vec)[iii]]
    sis2.vec <- rep(0, length(sis2.snps)); sis2.vec[sis2.snps$REF=='C'] <- 1
    
    # Skew of bins for coloring
    if(length(cut.points)==1){
      skew.windows <- data.frame('sister1_skew'= round(mean(sis1.vec), digits=2), 'sister2_skew'=round(mean(sis2.vec), digits=2))
      
    } else {
      
      skew.windows <- data.frame(do.call(rbind, lapply(1:length(cut.points), function(h){
        
        # If chromosome represents only a single segment, take all info
        
        # Otherwise, if there is more than one segment....
        # For first segment
        if(h==1){
          tmp.seg <- GRanges(seqnames = names(ch.vec)[iii], ranges = IRanges(start=1, end= (start(rel.snps)[ cut.points[h] ] - 1) ))
        } 
        
        # For last segment
        else if(h==length(cut.points)){
          tmp.seg <- GRanges(seqnames = names(ch.vec)[iii], ranges = IRanges(start= (start(rel.snps)[ cut.points[h-1] ] ), end= as.numeric(ch.vec[iii]) ))
        }
        
        # For segments in between
        else {
          tmp.seg <- GRanges(seqnames = names(ch.vec)[iii], ranges = IRanges(start= (start(rel.snps)[ cut.points[h-1] ] ), end= (start(rel.snps)[ cut.points[h] ] - 1)  ))
          
        }
        sis1.ov <- queryHits(findOverlaps(sis1.snps, tmp.seg)); sis2.ov <- queryHits(findOverlaps(sis2.snps, tmp.seg))
        c('sister1_skew'=round(mean(sis1.vec[sis1.ov]), digits=2), 'sister2_skew'=round(mean(sis2.vec[sis2.ov]), digits=2))
        
      })))
      
    }
    
    # Position of segment start and end points, stitched genome
    if(iii==1){ current.begin = 1 } else{ current.begin <- sum(ch.vec[1:(iii-1)]) } 
    
    # Local chromosome segment
    if(length(cut.points)==1){
      chr.start <- 1
      chr.end <- as.numeric(ch.vec[iii])
    } else {
      chr.start <- c(1, (start(rel.snps[cut.points[-length(cut.points)]])-1) ) 
      chr.end <- start(rel.snps)[cut.points]
    }
    
    # Set last point to end of chromosomes
    if(chr.end[length(chr.end)] <- ch.vec[iii]){ chr.end[length(chr.end)] <- ch.vec[iii]}
    
    #Rresult data frame of start, end, skew and chromosome
    results <- data.frame('chrom'=rep(names(ch.vec)[iii], length(cut.points)),
                          'start_gen'=c(current.begin,  (current.begin + as.numeric(start(rel.snps)[cut.points[-length(cut.points)]]) )),
                          'start_chr'=chr.start,
                          'end_gen'=c( (current.begin + as.numeric(start(rel.snps)[cut.points]))),
                          'end_chr'=chr.end,
                          'SCE_site_est'=sce_site,
                          'SCE_site_uncertainty'=sce_ambiguity,
                          skew.windows)
    if(results$end_gen[nrow(results)] < sum(ch.vec[1:iii])){
      results$end_gen[nrow(results)] <- sum(ch.vec[1:iii])
    }
    return(results)
    
  }))
  return(segment.types)
})
names(seg.maps) <- paste0('pair', 1:length(seg.maps))



#######  Panel 6b  ########
#######------------########
# Example clone and chromosome
clone <- 3
chrom <- 'chr2'

# Variables for plotting
segment.colors <- colorRampPalette(c('dodgerblue4', 'dodgerblue3', 'gray60', 'gray70', 'gray80', 'goldenrod2', 'goldenrod3'))(101)
segment.index <- seq(0,1,by=0.01)
ch.lines <- c(1, unlist(lapply(1:length(ch.vec), function(x){ sum(ch.vec[1:x]) })))
ch.slice <-c(0, unlist(lapply(1:length(ch.vec), function(x){sum(ch.vec[1:x])})))
ch.slice <- ch.slice[-length(ch.slice)]; names(ch.slice) <- names(ch.vec)
pt.sz <- 1.2; clr.trns <- 0.7

# Plot Panel 6b
par(mfrow=c(3,2), mar=c(4,4,4,1))
for(chrom in c('chr17', 'chr8', 'chr2')){
  # Get points for making scatter and plot scatter, use C > T and G > A
  sis1.gr <- rds.sis[[snp.match[clone]]]; sis2.gr <- rds.sis[[(snp.match[clone]+1)]] 
  sis1.gr <- sis1.gr[((sis1.gr$REF=='C' & sis1.gr$ALT=='T') | (sis1.gr$REF == 'G' & sis1.gr$ALT=='A')) & sis1.gr$TimePoint=='T1']
  sis2.gr <- sis2.gr[((sis2.gr$REF=='C' & sis2.gr$ALT=='T') | (sis2.gr$REF == 'G' & sis2.gr$ALT=='A')) & sis2.gr$TimePoint=='T1']
  sis1.gr <- sis1.gr[as.character(seqnames(sis1.gr))==chrom]; sis2.gr <- sis2.gr[as.character(seqnames(sis2.gr))==chrom]
  
  # Distance to nearest mutation
  sis1.dist <- log10(values(distanceToNearest(sis1.gr))$distance + 1)
  sis1.dist[sis1.gr$REF=='G'] <- sis1.dist[sis1.gr$REF=='G'] * -1
  sis2.dist <- log10(values(distanceToNearest(sis2.gr))$distance + 1)
  sis2.dist[sis2.gr$REF=='G'] <- sis2.dist[sis2.gr$REF=='G'] * -1
  
  s1.is.c <- sis1.gr$REF=='C'; s2.is.c <- sis2.gr$REF=='C' 
  plot.colors.s1 <- rep('dodgerblue4', length(sis1.gr)); plot.colors.s1[s1.is.c] <- 'goldenrod'; plot.colors.s1 <- adjustcolor(plot.colors.s1, alpha.f=clr.trns)
  plot.colors.s2 <- rep('dodgerblue4', length(sis2.gr)); plot.colors.s2[s2.is.c] <- 'goldenrod'; plot.colors.s2 <- adjustcolor(plot.colors.s2, alpha.f=clr.trns)
  
  plot(start(sis1.gr)/10^6, sis1.dist, ylim=c(-10,10), col=plot.colors.s1, ylab='dist. to nearest', pch=19, cex=pt.sz, xlab='pos in mb', main=chrom, cex.axis=1, cex.lab=1.2)
  abline(h=0, col=1, lty=3)
  plot(start(sis2.gr)/10^6, sis2.dist, ylim=c(-10,10), col=plot.colors.s2, ylab='dist. to nearest', pch=19, cex=pt.sz, xlab='pos in mb', main=chrom,  yaxt='n', cex.axis=1, cex.lab=1.2)
  abline(h=0, col=1, lty=3)
}

# Clear Plot
dev.off()



#######  Panel 6c  ########
#######------------########

# Variables for plotting
pt.sz <- 0.6
clr.trns <- 0.2
par(mfrow=c(2,1), mar=c(2,4,4,4))
chr.name.pt <- unlist(lapply(2:length(ch.lines), function(x){ (ch.lines[x]-ch.lines[(x-1)])/2 + ch.lines[(x-1)] }))

# Get points for making scatter and plot scatter, use C > T and G > A
sis1.gr <- rds.sis[[snp.match[clone]]]; sis2.gr <- rds.sis[[(snp.match[clone]+1)]] 
sis1.gr <- sis1.gr[((sis1.gr$REF=='C' & sis1.gr$ALT=='T') | (sis1.gr$REF == 'G' & sis1.gr$ALT=='A')) & sis1.gr$TimePoint=='T1']
sis2.gr <- sis2.gr[((sis2.gr$REF=='C' & sis2.gr$ALT=='T') | (sis2.gr$REF == 'G' & sis2.gr$ALT=='A')) & sis2.gr$TimePoint=='T1']

# Distance to nearest mutation
sis1.dist <- log10(values(distanceToNearest(sis1.gr))$distance + 1)
sis1.dist[sis1.gr$REF=='G'] <- sis1.dist[sis1.gr$REF=='G'] * -1
sis2.dist <- log10(values(distanceToNearest(sis2.gr))$distance + 1)
sis2.dist[sis2.gr$REF=='G'] <- sis2.dist[sis2.gr$REF=='G'] * -1

# Adjust graph position to stiched together chromosomes
sis1.df <- data.frame(sis1.gr, stringsAsFactors = FALSE)
graph_pos_s1 <- unlist(lapply(1:nrow(sis1.df), function(x) { 
  if(as.character(sis1.df$seqnames[x])=='chr1') { as.numeric(sis1.df$start[x])
  } else {  sum(as.numeric(ch.vec[1:(grep(as.character(sis1.df$seqnames[x]),names(ch.vec))-1)])) + as.numeric(sis1.df$start[x] )} }))

sis2.df <- data.frame(sis2.gr, stringsAsFactors = FALSE)
graph_pos_s2 <- unlist(lapply(1:nrow(sis2.df), function(x) { 
  if(as.character(sis2.df$seqnames[x])=='chr1') { as.numeric(sis2.df$start[x])
  } else {  sum(as.numeric(ch.vec[1:(grep(as.character(sis2.df$seqnames[x]),names(ch.vec))-1)])) + as.numeric(sis2.df$start[x] )} }))

# Point color
s1.is.c <- sis1.gr$REF=='C'; s2.is.c <- sis2.gr$REF=='C' 
plot.colors.s1 <- adjustcolor(c(rep('dodgerblue4', sum(!(s1.is.c))), rep('goldenrod', sum((s1.is.c)))), alpha.f = 0.4)
plot.colors.s2 <- adjustcolor(c(rep('dodgerblue4', sum(!(s2.is.c))), rep('goldenrod', sum((s2.is.c)))), alpha.f = 0.4)

# Segment boundary box plotting variables
box.boundaries <- 11; seg.box <- c(-box.boundaries, box.boundaries, box.boundaries, -box.boundaries)

# Add skew boxes
segment.df <- seg.maps[[clone]]

# Plot Panel 6c
# Plot scatter for Sister 1
plot(c(graph_pos_s1[!(s1.is.c)], graph_pos_s1[s1.is.c]), c(sis1.dist[!(s1.is.c)], sis1.dist[s1.is.c]), ylim=c(-10,10),
     col=plot.colors.s1, ylab='dist to nearest (log10)', pch=19, cex=pt.sz, xlab='position', main='sister 1', 
     xlim=c(0,sum(ch.vec)), xaxt='n')

for(iii in 1:nrow(segment.df)){
  polygon(rep(c(segment.df$start_gen[iii], segment.df$end_gen[iii]), each=2), seg.box, border = NA,
          col=adjustcolor(segment.colors[which( as.character(segment.index) == as.character(segment.df[iii,'sister1_skew']) )], alpha=clr.trns) )
}

# Add chromosome lines and possibly SCE sites
abline(v=ch.lines[-length(ch.lines)], col='black', lwd=1.2); abline(h=0, col=1, lty=3)
mtext(gsub('chr','', names(ch.vec)), at=chr.name.pt, side=3)

# Sister 2 mutations
par(mar=c(2,4,4,4))
plot(c(graph_pos_s2[!(s2.is.c)], graph_pos_s2[s2.is.c]), c(sis2.dist[!(s2.is.c)], sis2.dist[s2.is.c]), ylim=c(-10,10),
     col=plot.colors.s2, ylab='dist. to nearest (log10)', pch=19, cex=pt.sz, xlab='position', main='sister 2', xlim=c(0,sum(ch.vec)),
     xaxt='n')

# Add skew boxes for sister 2
for(iii in 1:nrow(segment.df)){
  polygon(rep(c(segment.df$start_gen[iii], segment.df$end_gen[iii]), each=2), seg.box, border = NA,
          col=adjustcolor(segment.colors[which( as.character(segment.index) == as.character(segment.df[iii,'sister2_skew']) )], alpha=clr.trns) )
}

# Add lines for chromosomes, and possible SCE sites
abline(v=ch.lines[-length(ch.lines)], col='black', lwd=1.2); abline(h=0, col=1, lty=3)
mtext(gsub('chr','', names(ch.vec)), at=chr.name.pt, side=3)

# Clear Plot
dev.off()



#######  Panel 6d  ########
#######------------########
# Sliding genomic windows
gen.gr <- makeGRangesFromDataFrame(data.frame('chromosome'=names(ch.vec), 'start'=rep(1, length(ch.vec)), 'end'=as.numeric(ch.vec), stringsAsFactors = FALSE),
                                   start.field='start', end.field='end', seqnames.field = 'chromosome')
gt.win <- unlist(slidingWindows(gen.gr, width = 1*10^7, step = 10^6))

# Chromosome mapping boundaries
chrom.breaks <- round(scales::rescale(ch.lines,to=c(1,1000))[-c(1, (length(ch.lines)+1) )], digits=0)
map.order <- c(3,1,2,4,5,6,7) # Order clones so example clone of 6c is at top

# Scale segments
seg.mat <- do.call(rbind, lapply(map.order, function(x){
  rescaled.points <- round(scales::rescale(as.numeric(c(as.numeric(seg.maps[[x]]$start_gen), as.numeric(seg.maps[[x]]$end_gen))), to=c(1,1000)), digits=0)
  rescale.start <- rescaled.points[1:nrow(seg.maps[[x]])]; rescale.start[1] <- 0
  rescale.end <- rescaled.points[c((nrow(seg.maps[[x]])+1):(nrow(seg.maps[[x]]) * 2))] - rescale.start
  
  return( rbind(as.numeric(unlist(lapply(1:nrow(seg.maps[[x]]), function(zz) { rep(seg.maps[[x]][zz,'sister1_skew'], rescale.end[zz]) }))),
                as.numeric(unlist(lapply(1:nrow(seg.maps[[x]]), function(zz) { rep(seg.maps[[x]][zz,'sister2_skew'], rescale.end[zz]) })))) )
}))

# Calculate skew
segment.skews <- do.call(cbind, lapply(names(rds.sis), function(x){
  tmp.gr <- rds.sis[[x]]
  tmp.gr <- tmp.gr[(tmp.gr$REF=='C' & tmp.gr$ALT=='T') | (tmp.gr$REF=='G' & tmp.gr$ALT=='A') & tmp.gr$TimePoint=='T1']
  results <- rep(-1, length(tmp.gr)); results[tmp.gr$REF=='C'] <- 1
  
  snp.win.ov <- findOverlaps(gt.win, tmp.gr)
  res <- rep(NA, length(gt.win)); names(res) <- 1:length(gt.win)
  win.skews <- unlist(lapply(split(results[subjectHits(snp.win.ov)], queryHits(snp.win.ov)), mean))
  res[names(win.skews)] <- as.numeric(win.skews)
  res
}))
colnames(segment.skews) <- names(rds.sis)
rownames(segment.skews) <- 1:length(gt.win)

# Number of snps in each segment for each sample
s.in.win <- do.call(cbind, lapply(names(rds.sis), function(x){
  tmp.gr <- rds.sis[[x]]
  tmp.gr <- tmp.gr[(tmp.gr$REF=='C' & tmp.gr$ALT=='T') | (tmp.gr$REF=='G' & tmp.gr$ALT=='A') & tmp.gr$TimePoint=='T1']
  snp.win.ov <- findOverlaps(tmp.gr, gt.win)
  res.vector <- rep(0, length(gt.win)); names(res.vector) <- 1:length(res.vector)
  snps.per.window <- (table(subjectHits(snp.win.ov)))
  res.vector[names(snps.per.window)] <- as.numeric(snps.per.window)
  res.vector
}))

# Take only segments that have at least 3 SNVs
segment.skews <- segment.skews[(rowSums(s.in.win >= 3)==ncol(segment.skews)),]

# Calculate correlation between segments
seg.cor <- cor(segment.skews, )
for(i in 1:ncol(seg.cor)){seg.cor[i,i] <- NA}

# Pair Annotation
annotation <- data.frame(Pair = rep(1:7, each=2))
rownames(annotation) <- colnames(seg.cor) # check out the row names of annotation
pair.colors <- grey.colors(7); names(pair.colors) <- 1:length(pair.colors)
ann_colors = list( Pair = pair.colors )
rownames(seg.mat) <- rownames(annotation)

# Plot Panel 6d
pheatmap(seg.mat, color=segment.colors, show_colnames = FALSE, show_rownames=FALSE, border_color = NA, main='Skew Segmentation', 
         cluster_rows = FALSE, cluster_cols = FALSE, gaps_col  = chrom.breaks, gaps_row =seq(2,nrow(seg.mat), by=2), cellheight = 15, 
         annotation_row = annotation, annotation_colors = ann_colors, cellwidth=0.65)

# Clear Plot
dev.off()



#######  Panel 6e  ########
#######------------########

# Sliding genomic windows
gen.gr <- makeGRangesFromDataFrame(data.frame('chromosome'=names(ch.vec), 'start'=rep(1, length(ch.vec)), 'end'=as.numeric(ch.vec), stringsAsFactors = FALSE),
                                   start.field='start', end.field='end', seqnames.field = 'chromosome')
gt.win <- unlist(slidingWindows(gen.gr, width = 1*10^7, step = 10^6))

# Calculate skew in sliding windows
segment.skews <- do.call(cbind, lapply(names(rds.sis), function(x){
  tmp.gr <- rds.sis[[x]]
  tmp.gr <- tmp.gr[(tmp.gr$REF=='C' & tmp.gr$ALT=='T') | (tmp.gr$REF=='G' & tmp.gr$ALT=='A') & tmp.gr$TimePoint=='T1']
  results <- rep(-1, length(tmp.gr)); results[tmp.gr$REF=='C'] <- 1
  
  snp.win.ov <- findOverlaps(gt.win, tmp.gr)
  res <- rep(NA, length(gt.win)); names(res) <- 1:length(gt.win)
  win.skews <- unlist(lapply(split(results[subjectHits(snp.win.ov)], queryHits(snp.win.ov)), mean))
  res[names(win.skews)] <- as.numeric(win.skews)
  res
}))
colnames(segment.skews) <- names(rds.sis)
rownames(segment.skews) <- 1:length(gt.win)

# Number of mutations in each segment for each sample
s.in.win <- do.call(cbind, lapply(names(rds.sis), function(x){
  tmp.gr <- rds.sis[[x]]
  tmp.gr <- tmp.gr[(tmp.gr$REF=='C' & tmp.gr$ALT=='T') | (tmp.gr$REF=='G' & tmp.gr$ALT=='A') & tmp.gr$TimePoint=='T1']
  snp.win.ov <- findOverlaps(tmp.gr, gt.win)
  res.vector <- rep(0, length(gt.win)); names(res.vector) <- 1:length(res.vector)
  snps.per.window <- (table(subjectHits(snp.win.ov)))
  res.vector[names(snps.per.window)] <- as.numeric(snps.per.window)
  res.vector
}))

# Take only segments that have at least 3 mutations
segment.skews <- segment.skews[(rowSums(s.in.win >= 3)==ncol(segment.skews)),]

# Calculate correlation between segments
seg.cor <- cor(segment.skews, )
for(i in 1:ncol(seg.cor)){seg.cor[i,i] <- NA}

# Heatmap of correlation between sisters and clones
pheatmap(seg.cor, col=viridis(100), na_col = 'white', cluster_rows = FALSE, cluster_cols = FALSE, annotation_col = annotation,
         annotation_colors = ann_colors, annotation_row = annotation, show_rownames = FALSE, show_colnames = FALSE)

# Get average R squared values for variance explanation in text
pair.cors <- c()
clone.cors <- c()
first.row <- 2:nrow(seg.cor)
for(i in 1:(ncol(seg.cor)-1)){
  message(first.row[i])
  # Even numbers have no pairs correlation
  if((i %% 2) == 0) {
    clone.cors <- c(clone.cors, seg.cor[first.row[i]:nrow(seg.cor), i])
  }
  
  # Odd numbers have pairs correlation, rest clones
  else {
    if(first.row[i]==nrow(seg.cor)){
      pair.cors <- c(pair.cors, seg.cor[first.row[i],i])
    } else {
      pair.cors <- c(pair.cors, seg.cor[first.row[i],i])
      clone.cors <- c(clone.cors, seg.cor[(first.row[i] + 1):nrow(seg.cor), i])
    }
  }
}

# Scatter plot colors
jet.colors <- grDevices::colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

# Plot individual scatters
par(mfrow=c(2,1), mar=c(4,4,4,4))
smoothScatter(segment.skews[,"PF1_CCS1"], segment.skews[,"PF1_CCS2"], colramp = jet.colors, nbin=256, bandwidth=c(0.03,0.03), nrpoints=0, 
                xlim=c(-1,1), ylim=c(-1,1), xlab='sister 1 skew', ylab='sister 2 skew')
abline(0, -1, col='black', lty=3, lwd=2)

smoothScatter(segment.skews[,"PF1_CES2"], segment.skews[,"PF1_CGS1"], colramp = jet.colors, nbin=256, bandwidth=c(0.03,0.03), nrpoints=0, 
                xlim=c(-1,1), ylim=c(-1,1), xlab='clone 1 skew', ylab='clone 2 skew')
abline(0, -1, col='black', lty=3, lwd=2)

# Clear Plot
dev.off()



#######  Panel 6g  ########
#######------------########

# Vector to assign snp pairs to segmentation
snv.index <- lapply(snp.match, function(x){ return(names(rds.sis)[c(x, (x+1))])}); names(snv.index) <- names(seg.maps)

##  Assign mutations from each pair to strand (recombined or not)
##---------------------------------------------------------------
sce.map <- lapply(names(seg.maps), function(i){
  message(i)
  current.map <- seg.maps[[i]]
  
  # Select sister  pairs, use T1 UV SNVs
  sis1.snvs <- rds.sis[[snv.index[[i]][1]]]; sis1.snvs <- sis1.snvs[sis1.snvs$TimePoint=='T1']
  sis1.snvs <- sis1.snvs[(sis1.snvs$REF=='C' & sis1.snvs$ALT=='T')|(sis1.snvs$REF=='G' & sis1.snvs$ALT=='A')]
  elementMetadata(sis1.snvs)$SCE_status <- rep('mixed', length(sis1.snvs))
  
  sis2.snvs <- rds.sis[[snv.index[[i]][2]]]; sis2.snvs <- sis2.snvs[sis2.snvs$TimePoint=='T1']
  sis2.snvs <- sis2.snvs[(sis2.snvs$REF=='C' & sis2.snvs$ALT=='T')|(sis2.snvs$REF=='G' & sis2.snvs$ALT=='A')]
  elementMetadata(sis2.snvs)$SCE_status <- rep('mixed', length(sis2.snvs))
  
  # Rescale skew to -1, 1
  current.map$sister1_skew <- rescale(x = current.map$sister1_skew, newrange = c(-1,1))
  current.map$sister2_skew <- rescale(x = current.map$sister2_skew, newrange = c(-1,1))
  
  # Split on chromosomes. Splits with  > 1 row have at least one crossover event
  chrom.split <- split(current.map, current.map$chr)
  chrom.split <- chrom.split[unlist(lapply(chrom.split, nrow)) > 1]
  
  # Loop through chromsomes 
  for(ii in names(chrom.split)) {
    current.chromosome <- chrom.split[[ii]]
    
    # Loop through segments, negating last
    for(iii in 1:(nrow(current.chromosome) - 1)){
      
      # Determine mixed side
      left.side <- abs((current.chromosome$sister1_skew[iii]) - (current.chromosome$sister2_skew[(iii)]) )
      right.side <- abs((current.chromosome$sister1_skew[(iii+1)]) - (current.chromosome$sister2_skew[(iii+1)]) )
      mixed.side <- which.min(c(left.side, right.side))
      
      # Mixed side on left
      if(mixed.side == 1){
        # GRanges of 1mb window around break
        sce.window <- data.frame('chrom'=ii, 'start'=current.chromosome$start_chr[iii], 'end'=current.chromosome$SCE_site_est[iii])
        sce.window <- makeGRangesFromDataFrame(sce.window, seqnames.field = 'chrom')
        
        # Overlap vectors for both sisters
        overlap.vec1 <- rep(FALSE, length(sis1.snvs)); overlap.vec1[queryHits(findOverlaps(sis1.snvs, sce.window))] <- TRUE
        overlap.vec2 <- rep(FALSE, length(sis2.snvs)); overlap.vec2[queryHits(findOverlaps(sis2.snvs, sce.window))] <- TRUE
        
        # Assign recombined strands to sisters and extract relevant mutations
        if(current.chromosome$sister1_skew[iii+1] > current.chromosome$sister2_skew[iii+1]){
          sis1.snvs$SCE_status[(overlap.vec1 & sis1.snvs$REF=="G" & sis1.snvs$ALT=="A")] <- 'SCE'
          sis1.snvs$SCE_status[(overlap.vec1 & sis1.snvs$REF=="C" & sis1.snvs$ALT=="T")] <- 'Non_SCE'
          sis2.snvs$SCE_status[(overlap.vec2 & sis2.snvs$REF=="C" & sis2.snvs$ALT=="T")] <- 'SCE'
          sis2.snvs$SCE_status[(overlap.vec2 & sis2.snvs$REF=="G" & sis2.snvs$ALT=="A")] <- 'Non_SCE'
          
        } else{
          sis1.snvs$SCE_status[(overlap.vec1 & sis1.snvs$REF=="C" & sis1.snvs$ALT=="T")] <- 'SCE'
          sis1.snvs$SCE_status[(overlap.vec1 & sis1.snvs$REF=="G" & sis1.snvs$ALT=="A")] <- 'Non_SCE'
          sis2.snvs$SCE_status[(overlap.vec2 & sis2.snvs$REF=="G" & sis2.snvs$ALT=="A")] <- 'SCE'
          sis2.snvs$SCE_status[(overlap.vec2 & sis2.snvs$REF=="C" & sis2.snvs$ALT=="T")] <- 'Non_SCE' 
        }
        
      }
      # On right
      else{
        # GRanges of 1mb window around break
        sce.window <- data.frame('chrom'=ii, 'start'=current.chromosome$SCE_site_est[(iii)], 'end'=(current.chromosome$end_chr[(iii+1)] + 10^6) )
        sce.window <- makeGRangesFromDataFrame(sce.window, seqnames.field = 'chrom')
        
        # Overlap vectors for both sisters
        overlap.vec1 <- rep(FALSE, length(sis1.snvs)); overlap.vec1[queryHits(findOverlaps(sis1.snvs, sce.window))] <- TRUE
        overlap.vec2 <- rep(FALSE, length(sis2.snvs)); overlap.vec2[queryHits(findOverlaps(sis2.snvs, sce.window))] <- TRUE
        
        # Assign recombined strands to sisters and extract relevant mutations
        if(current.chromosome$sister1_skew[iii] > current.chromosome$sister2_skew[iii]){
          sis1.snvs$SCE_status[(overlap.vec1 & sis1.snvs$REF=="G" & sis1.snvs$ALT=="A")] <- 'SCE'
          sis1.snvs$SCE_status[(overlap.vec1 & sis1.snvs$REF=="C" & sis1.snvs$ALT=="T")] <- 'Non_SCE'
          sis2.snvs$SCE_status[(overlap.vec2 & sis2.snvs$REF=="C" & sis2.snvs$ALT=="T")] <- 'SCE'
          sis2.snvs$SCE_status[(overlap.vec2 & sis2.snvs$REF=="G" & sis2.snvs$ALT=="A")] <- 'Non_SCE' 
        } else{
          sis1.snvs$SCE_status[(overlap.vec1 & sis1.snvs$REF=="C" & sis1.snvs$ALT=="T")] <- 'SCE'
          sis1.snvs$SCE_status[(overlap.vec1 & sis1.snvs$REF=="G" & sis1.snvs$ALT=="A")] <- 'Non_SCE'
          sis2.snvs$SCE_status[(overlap.vec2 & sis2.snvs$REF=="G" & sis2.snvs$ALT=="A")] <- 'SCE'
          sis2.snvs$SCE_status[(overlap.vec2 & sis2.snvs$REF=="C" & sis2.snvs$ALT=="T")] <- 'Non_SCE' 
          
        }
      }
    }  
  }
  res <- list(sis1.snvs, sis2.snvs); names(res) <- snv.index[[i]]
  return(res)
}); names(sce.map) <- names(seg.maps)


##  Number of SCE events passing filters
##--------------------------------------
slide.window <- 10^6 + 1
step.size <- 10^5
step.number <- 200
total.size <- ( (step.size * step.number) + (slide.window) )

sce.counts <- matrix(nrow=(length(sce.map)), ncol=4); rownames(sce.counts) <- names(sce.map)
colnames(sce.counts) <- c('tot.segs','near.end','under.10mb', 'filtered')

for(x in names(sce.map)){
  current.map <- seg.maps[[x]]
  
  # Remove the segments that represent the end of the chromosomes
  current.map <- current.map[!(is.na(current.map$SCE_site_est)),]
  sce.counts[x,'tot.segs'] <- nrow(current.map)
  
  # Remove elements that are less than 1/2 window from chromosome ends
  end.chr <- (current.map$SCE_site_est >= ((total.size + 1)/2)) & (current.map$SCE_site_est < (ch.vec[current.map$chrom] - ((total.size + 1)/2)))
  current.map <- current.map[end.chr, ]
  sce.counts[x,'near.end'] <- sum(!(end.chr))
  
  # Remove elements that are less than 10mb long
  small.fragments <- (current.map$end_chr - current.map$start_chr) >= ((total.size))
  current.map <- current.map[small.fragments,]
  sce.counts[x, 'under.10mb'] <- sum(!(small.fragments)); sce.counts[x,'filtered'] <- nrow(current.map)
}

##  Retrieve mutation density around an SCE
##-----------------------------------------
sce.slide <- Reduce('+', lapply(names(sce.map), function(x){
  message(x)
  
  # Retrieve map, filter to include SCE sites that are at least 10mb from chromosome ends
  current.map <- seg.maps[[x]]
  
  # Remove the segments that represent the end of the chromosomes
  current.map <- current.map[!(is.na(current.map$SCE_site_est)),]
  
  # Remove elements below minimum size
  current.map <- current.map[(current.map$SCE_site_est - current.map$start_chr) >= (total.size + 1)/2, ]
  
  # Remove elements that are less than set size from chromosome ends
  current.map <- current.map[current.map$SCE_site_est < (ch.vec[current.map$chrom] - ((total.size + 1)/2)), ]
  
  # Create GRanges of SCE site
  gen.gr <- makeGRangesFromDataFrame(data.frame('chromosome'=current.map$chrom, 
                                                'start'= as.numeric(current.map$SCE_site_est), 
                                                'end'= as.numeric(current.map$SCE_site_est), stringsAsFactors = FALSE), 
                                     start.field='start', end.field='end', seqnames.field = 'chromosome')
  gen.gr <- resize(gen.gr, fix='center', width=total.size)

  # Order each SCE site so separated strands are on the left
  window.order <- unlist(lapply(1:length(gen.gr), function(xx){
    left.flank <- resize(gen.gr[xx], fix = 'start', width = total.size/2 )
    sce.stat <- table(subsetByOverlaps(sce.map[[x]][[1]], left.flank)$SCE_status)
    
    if(names(which.max(sce.stat))=='mixed'){
      'reverse'
    } else{
      'forward'
    }
  }))
  
  # Loop through windows, reverse window order if it is on the wrong side and tabulate mutation identities
  t(Reduce('+', lapply(1:length(window.order), function(xxx){
    gt.win <- unlist(slidingWindows(gen.gr[xxx], width = slide.window, step = step.size))  
    if(window.order[xxx]=='forward'){
      gt.win <- gt.win
    } else {
      gt.win <- rev(gt.win)
    }
    
    do.call(rbind, lapply(1:length(gt.win), function(v){
      sis1.res <- rep(0, 3); names(sis1.res) <- unique(sce.map[[1]][[1]]$SCE_status)
      sis2.res <- rep(0, 3); names(sis2.res) <- unique(sce.map[[1]][[1]]$SCE_status)
      sis1.cts <- table(subsetByOverlaps(sce.map[[x]][[1]], gt.win[v])$SCE_status)
      sis1.res[names(sis1.cts)] <- as.numeric(sis1.cts)
      
      sis2.cts <- table(subsetByOverlaps(sce.map[[x]][[2]], gt.win[v])$SCE_status)
      sis2.res[names(sis2.cts)] <- as.numeric(sis2.cts)
      
      count.res <- c(as.numeric(sis1.res), as.numeric(sis2.res))
      names(count.res) <- c(paste0('sis1_', names(sis1.res)), paste0('sis2_', names(sis2.res)))
      count.res
    }))
    
  })))
  
}))

# Plot Panel 6g
plot(1:ncol(sce.slide)-(0.5*ncol(sce.slide)), runmean( Rle((colSums(sce.slide[c(3,6),])/2)/(sum(sce.counts[,'filtered']))), k=5, endrule = 'constant'), pch=19, 
     col='darkgrey', type='l', lwd=2, ylab='mu/mb',  xlab='distance to SCE (10^5)', ylim=c(0,2.7))
  lines(1:ncol(sce.slide)-(0.5*ncol(sce.slide)), runmean( Rle(colSums(sce.slide[c(2,5),])/sum(sce.counts[,'filtered'])), k=5, endrule = 'constant'), col='coral4', lwd=2)
  lines(1:ncol(sce.slide)-(0.5*ncol(sce.slide)), runmean( Rle(colSums(sce.slide[c(1,4),])/sum(sce.counts[,'filtered'])), k=5, endrule = 'constant'), col='coral', lwd=2)
  abline(v=0, col='darkgrey', lty=3, lwd=2)
  legend('topleft', legend=c('SCE strand','non-SCE strand'), fill=c('coral4','coral'), bty='n', border=NA)
  legend('topright', legend=c('mix'), fill=c('darkgrey'), bty='n', border=NA)
