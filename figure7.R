############################
#####     Figure 7     #####  
############################
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
require(BSgenome.Mmusculus.UCSC.mm10)
require(changepoint)
require(pheatmap)
require(mixtools)
require(scales)
require(grid)

# Chromosome size vector
txdb <- keepStandardChromosomes(TxDb.Mmusculus.UCSC.mm10.knownGene)
ch.vec <- seqlengths(txdb); ch.vec <- ch.vec[!(names(ch.vec) %in% c('chrY', 'chrM'))]

# Mutations called using 99588 as background, which has higher depth and thus few X chromosome mutations
rds.file <- 'Path to F1CastB6_mutations.rds'
rds <- readRDS(rds.file) ; rds <- rds[rds$FILTER]
rds.cd <- split(rds, rds$SAMPLE)

# Tumor example
tumor.sel <- 'MCB6A_T1'
snp.ex <- rds.cd[[tumor.sel]]
snp.dist <- values(distanceToNearest(snp.ex))$distance



#######  Panel 7b  ########
#######------------########

# Mutations signature set, always as reference C or T
mutations.of.interest <- c(paste0(rep('C',3),'_',c('A','G','T')), paste0(rep('T',3),'_',c('A','C','G')))

# All possible mutation types in each trinucleotide context
mut.vec <- paste0(c(rep(paste0(rep(c('A','C','G','T'), each=4),'C', rep(c('A','C','G','T'),4)),3),
                    rep(paste0(rep(c('A','C','G','T'), each=4),'T', rep(c('A','C','G','T'),4)),3)),'_',
                  rep(mutations.of.interest, each=16))

# Result and revcomp slicing vector
results <- rep(0, length(mut.vec)); names(results) <- mut.vec 
rev.comp.vec <- c('A', 'C', 'G', 'T'); names(rev.comp.vec) <- c('T', 'G', 'C', 'A')

# Get bp change, reverse complement G and A mutations to represent everything as C or T mutations
current.gr <- snp.ex

# Retrieve trinucleotide context
current.bases <- resize(current.gr, fix = 'center', width=3)
trinucs <- Biostrings::getSeq(x=BSgenome.Mmusculus.UCSC.mm10, current.bases)
trinucs[current.gr$REF=='G' | current.gr$REF=='A'] <- reverseComplement(trinucs[current.gr$REF=='G' | current.gr$REF=='A'])

# Reverse complement G and A mutations to represent everything as a C or T reference
ref.base <- unlist(lapply(as.character(trinucs), function(x){unlist(strsplit(x,''))[[2]]}))
alt.base <- current.gr$ALT; alt.base[current.gr$REF %in% c('G','A')] <- rev.comp.vec[alt.base[current.gr$REF %in% c('G','A')]]

# Create names for mutations
trinucs <- paste0(as.character(trinucs),'_',ref.base,'_',alt.base)
trinuc.ct <- table(trinucs)
trinuc.ct <- trinuc.ct[(names(trinuc.ct) %in% names(results))]
results[names(trinuc.ct)] <- as.numeric(trinuc.ct)

# Colors
col.trinuc <- rep(c('dodgerblue1', 'black', 'firebrick2', 'darkgrey', 'limegreen', 'salmon'), each=16)

# Plot Panel 7b
par(mar=c(4,6,3,3))
barplot(results/(sum(results)), col=col.trinuc, las=2, cex.names = 0.5, ylab='frequency', 
        main='Mutation Signature', border=NA, xaxt='n', ylim=c(-0.02,0.08), cex.axis=1)
rect(0.2,-0.015,19.1,-0.005, col=col.trinuc[1], border = NA)
rect(19.2,-0.015,38.4,-0.005, col=col.trinuc[17], border = NA)
rect(38.5,-0.015,57.6,-0.005, col=col.trinuc[33], border = NA)
rect(57.7,-0.015,76.9,-0.005, col=col.trinuc[49], border = NA)
rect(77,-0.015,96.1,-0.005, col=col.trinuc[65], border = NA)
rect(96.2,-0.015,115.3,-0.005, col=col.trinuc[81], border = NA)
text(gsub('_','>',mutations.of.interest), x=c(10,seq(29,105,by=19)), y=-0.01, col='white')

# Clear Plot
dev.off()



#######  Panel 7c & g  ########
#######----------------########

# Plotting variables
segment.colors <- colorRampPalette(c('dodgerblue4','grey','goldenrod'))(101)
segment.index <- seq(0,1,by=0.01)
ch.lines <- c(1, unlist(lapply(1:length(ch.vec), function(x){ sum(ch.vec[1:x]) })))

## Segmentation 
##-------------
skew.mat <- lapply(names(rds.cd), function(z){ 

  skew.tmp <- lapply(c('C3H', 'CAST', ''), function(sp){

    do.call(rbind, lapply(1:(length(ch.vec)), function(i){

      # Get granges and data frame. Select C3H allele (check order by position)
      tmp.gr <- rds.cd[[z]][(as.character(seqnames(rds.cd[[z]]))==names(ch.vec)[i]) & (grepl(sp, rds.cd[[z]]$species_tag)) ]
      tmp.gr <- tmp.gr[(tmp.gr$REF=='A' | tmp.gr$REF=='T')]
      if(length(tmp.gr) > 10){
        
        # Create vector of -1 to 1
        bin.vec <- rep(0, length(tmp.gr)); bin.vec[tmp.gr$REF=='T'] <- 1
        tmp.dist <- values(distanceToNearest(tmp.gr))$distance
        # Calculate fit points
        cut.points <- cpt.mean(bin.vec, penalty='AIC', method='PELT', class=FALSE,  Q = 12)
        # Skew of bins for coloring
        res <- c()
        for(h in 1:length(cut.points)){
          if(h==1){ tmp <- bin.vec[1:(cut.points[h]-1)]
          res <- c(res,round(mean(tmp), digits=2) )
          next
          } 
          else if(h==length(cut.points)){
            res <- c(res, round(mean(bin.vec[cut.points[h-1]:cut.points[h]]), digits=2) )
          }
          else {
            tmp <- bin.vec[cut.points[(h-1)]:(cut.points[h]-1)]
            res <- c(res, round(mean(tmp), digits=2))
          }
        }
        # Chr specific pattern
        results.chr <- data.frame('start_chr'=c(1,  (1 + start(tmp.gr)[cut.points[-length(cut.points)]])),
                                  'end_chr'=c( start(tmp.gr)[cut.points]))
        if(results.chr$end_chr[nrow(results.chr)] < ch.vec[i]){
          results.chr$end_chr[nrow(results.chr)] <- ch.vec[i]
        }
        # position of segment start and end points for genome skew
        if(i==1){ current.begin = 0 } else{ current.begin <- sum(ch.vec[1:(i-1)]) } 
        # result data frame of start, end, skew and chromosome
        results <- data.frame('start_gen'=c(current.begin+1,  (as.numeric(current.begin) + (start(tmp.gr)[cut.points[-length(cut.points)]])) + 1),
                              'end_gen'=c( (as.numeric(current.begin) + start(tmp.gr)[cut.points])),
                              'skew'=res, 'chrom'=rep(names(ch.vec)[i], length(cut.points)))
        if(results$end_gen[nrow(results)] < sum(ch.vec[1:i])){
          results$end_gen[nrow(results)] <- sum(ch.vec[1:i])
        }
        res.final <- cbind(results.chr, results)
        return(res.final)
      }
    }))
  }); names(skew.tmp) <- c('C3H', 'CAST', 'Unsegregated')
  return(skew.tmp)
}); names(skew.mat) <- names(rds.cd)

# Adjusted genomic positions for juxtaposed chromsomes
gen.adj <- c(1, unlist(lapply(1:(length(ch.vec)-1), function(x){ sum(ch.vec[1:x])}))); names(gen.adj) <- names(ch.vec)
pt.adj <- start(snp.ex) + gen.adj[as.character(seqnames(snp.ex))]

# Point colors for C and T mutations
is.a <- snp.ex$REF == 'A'; is.t <- snp.ex$REF == 'T'
point.colors <- rep('gray', length(snp.ex)); point.colors[is.t] <- 'goldenrod'; point.colors[is.a] <- 'dodgerblue4'
point.colors <- alpha(point.colors,alpha = 0.3)

# Chromosome name locations
chr.name.pt <- unlist(lapply(2:length(ch.lines), function(x){ (ch.lines[x]-ch.lines[(x-1)])/2 + ch.lines[(x-1)] }))

# Plot dimensions
par( mar=c(1,5,1,0))
m <- matrix(c(rep(1,9),rep(2,12),rep(3,9), rep(4,9),rep(5,12),rep(6,9), rep(7,9),rep(8,12),rep(9,9)), ncol=3,byrow = TRUE)
layout(m)

# Plot panels 7c & g
for(spec.name in names(skew.mat[[tumor.sel]])){
  cs <- c()
  ifelse(spec.name == 'Unsegregated', cs <- c('C3H', 'CAST', 'Undetermined'), cs <- spec.name)
  spec.points <- (snp.ex$species_tag %in% cs)
  sk.un <- skew.mat[[tumor.sel]][[spec.name]]
  
  # A bases
  plot(pt.adj[is.a & spec.points], log10(snp.dist[(is.a & spec.points)]+1), pch=19, cex=0.4, col=point.colors[(is.a & spec.points)], frame.plot=FALSE, xaxt='n',
       ylab='dist. (log10)', xlab='', ylim=c(2,6), cex.axis=1.5, xlim=c(1,sum(ch.vec)), main=spec.name)
  abline(v=ch.lines, col=1, lty=6)
  
  
  # Blank plot
  plot(pt.adj[is.a & spec.points], log10(snp.dist[(is.a & spec.points)]+1), pch=19, cex=0.1, col='white', 
       frame.plot=FALSE, xaxt='n',yaxt='n', ylab='', xlab='', ylim=c(0,1), xlim=c(1,sum(ch.vec)))
  
  # Chromosome txt
  mtext(gsub('chr','', names(ch.vec)), at=chr.name.pt, side=3)
  
  # Add skew boxes
  for(iii in 1:nrow(sk.un)){
    if(sk.un$skew[iii] < 0.2){
      polygon(rep(c(sk.un$start_gen[iii],sk.un$end_gen[iii]), each=2), c(0.66,0.99,0.99,0.66), border = NA,
              col=segment.colors[which( as.character(segment.index) == as.character(sk.un$skew[iii]) )])
      
      
    } else if(sk.un$skew[iii] > 0.8) {
      polygon(rep(c(sk.un$start_gen[iii],sk.un$end_gen[iii]), each=2), c(0,0.33,0.33,0), border = NA,
              col=segment.colors[which( as.character(segment.index) == as.character(sk.un$skew[iii]) )])
    } else{
      polygon(rep(c(sk.un$start_gen[iii],sk.un$end_gen[iii]), each=2), c(0.33,0.66,0.66,0.33), border = NA,
              col=segment.colors[which( as.character(segment.index) == as.character(sk.un$skew[iii]) )])
      polygon(rep(c(sk.un$start_gen[iii],sk.un$end_gen[iii]), each=2), c(0.33,0.66,0.66,0.33), border = 'red',
              col=NA, lwd=2)
    }}
  abline(v=ch.lines, h=c(0,1), lty=2)
  
  # T bases
  plot(pt.adj[(is.t & spec.points)], -1 * log10(snp.dist[(is.t & spec.points)]+1), pch=19, cex=0.4, 
       col=point.colors[(is.t & spec.points)], frame.plot=FALSE,xlim=c(1,sum(ch.vec)), ylab='dist. (log10)', 
       xlab='chromosome coordinate', ylim=c(-6,-2), cex.axis=1.5, xaxt='n')
  abline(v=ch.lines, col=1, lty=6)
}

# Clear Plot
dev.off()



#######  Panel 7d  ########
#######------------########

# # Calculate read delta, add random noise
spec.dist <- unlist(lapply(1:length(rds.cd), function(cg) {
  rnd.noise <- sample(seq(0.5,1, by=0.01), length(rds.cd[[cg]]), replace = TRUE)
  (log2(rds.cd[[cg]]$ALT_C3H + rnd.noise) - log2(rds.cd[[cg]]$ALT_CAST + rnd.noise))
}))
 
# Total reads
tot.rds <- unlist(lapply(1:length(rds.cd), function(x){ rds.cd[[x]]$ALT_C3H + rds.cd[[x]]$ALT_CAST}))
 
##  Fit mixture models (haplotype determination)
##----------------------------------------------
spec.mix <- normalmixEM(spec.dist[tot.rds >= 2], lambda = .5, mu = c(-3.5, 3.5), sigma = c(1.2, 1.2))
cast.snps <- max(spec.dist[(tot.rds >= 2)][spec.mix$posterior[,'comp.1'] >= 0.99])
c3h.snps <- min(spec.dist[(tot.rds >= 2)][spec.mix$posterior[,'comp.2'] >= 0.99])

# Plot Panel 7d
plot(density(spec.mix$x), cex.axis = 1.4, cex.lab = 1.5, cex.main = 1.5,
     main = "Haplotype specific reads", xlab = "log2(C3H/CAST)", col='white')
polygon(density(spec.mix$x), col = 'darkgrey', border=NA)
abline(v=c(c3h.snps, cast.snps), lty=2, lwd=1.5)

# Clear Plot
dev.off()



#######  Panel 7e  ########
#######------------########

# Add species designation to Granges
for(i in names(rds.cd)){
  res <- rep('Undetermined', length(rds.cd[[i]])); tmp <- data.frame(rds.cd[[i]], stringsAsFactors = FALSE)
  rnd.noise <- sample(seq(0.1,0.5, by=0.01), nrow(tmp), replace = TRUE)
  snp.delta <- log2(tmp$ALT_C3H + rnd.noise)-log2(tmp$ALT_CAST + rnd.noise)
  total.reads <- sum(tmp$ALT_C3H, tmp$ALT_CAST)
  res[(snp.delta >= c3h.snps) & total.reads >= 2] <- 'C3H'
  res[(snp.delta <= cast.snps) & total.reads >= 2] <- 'CAST'
  elementMetadata(rds.cd[[i]])$species_tag <- res
}

# Data frame of number of all mutations and whether they are haplotype specific
spec.assign <- do.call(rbind, lapply(1:length(rds.cd), function(x){
  res <- c('C3H'=0, 'CAST'=0, 'Undetermined'=0)
  spec.cts <- table(rds.cd[[x]]$species_tag)
  res[names(spec.cts)] <- as.numeric(spec.cts)
  res
}))

# Plot Panel 7e
barplot(100, col='darkgrey', ylim=c(0,100), border=NA, cex.axis=1.5, ylab='% allelic assignment')
barplot(sum(colSums(spec.assign[,c(1:2)]))/sum(spec.assign) * 100, col='brown', add=TRUE, border=NA, yaxt='n')
barplot(sum(spec.assign[,1])/sum(spec.assign) * 100, col='gold4', add=TRUE, border=NA, yaxt='n')
text(0.7, 20, round( (sum(spec.assign[,1])/sum(spec.assign) * 100), digits=2 ), col='white')
text(0.7, 70, round( (sum(spec.assign[,2])/sum(spec.assign) * 100), digits=2 ), col='white')

# Clear Plot
dev.off()



#######  Panel 7f  ########
#######------------########

# Point adjustments for single chromosomes
example.chrs <- c('chr2','chr3')
chr.zoom <- as.character(seqnames(snp.ex)) %in% example.chrs
example.chr.xlim <- c(gen.adj[(min(which(names(ch.vec)%in% example.chrs)))], gen.adj[(max(which(names(ch.vec)%in% example.chrs)) + 1)])
ch.zoom.lines <- gen.adj[example.chrs]; ch.zoom.lines <- as.numeric(ch.zoom.lines[-1])

# Plot Panel 7f
par( mar=c(3,5,3,3))
m <- matrix(c(rep(rep(c(1,4,7), each=3),3),
              rep(rep(c(2,5,8), each=3),2),
              rep(rep(c(3,6,9), each=3),3)) , ncol=9, byrow = TRUE)
layout(m)
for(spec.name in names(skew.mat[[tumor.sel]])[c(2,3,1)] ){
  cs <- c()
  ifelse(spec.name == 'Unsegregated', cs <- c('CAST', 'Undetermined', 'C3H'), cs <- spec.name)
  spec.points <- (snp.ex$species_tag %in% cs)
  sk.un <- skew.mat[[tumor.sel]][[spec.name]]
  sk.un <- sk.un[(sk.un$chrom %in% example.chrs), ]
  
  # A bases
  plot(pt.adj[is.a & spec.points & chr.zoom], log10(snp.dist[(is.a & spec.points & chr.zoom)]+1), pch=19, cex=1.1, 
       col=point.colors[(is.a & spec.points & chr.zoom)], xaxt='n', main=spec.name,
       ylab='distance to nearest (log10)', xlab='', ylim=c(2,6), cex.axis=1.5, xlim=example.chr.xlim)
  abline(v=ch.zoom.lines, col=1, lty=6)
  
  if(spec.name=='CAST'){
    legend('topleft', col='dodgerblue4', pch=19, legend=c('A>N'), bty='n', cex=1.3)
  }
  
  # Blank plot
  plot(pt.adj[is.a & spec.points & chr.zoom], log10(snp.dist[(is.a & spec.points &chr.zoom)]+1), pch=19, cex=0.1, col='white', 
       frame.plot=FALSE, xaxt='n',yaxt='n', ylab='', xlab='', ylim=c(0,1), xlim=example.chr.xlim)
  
  mtext(names(ch.vec)[c(2,3)], at=chr.name.pt[c(2,3)], side=3)
  
  # Add skew boxes
  for(iii in 1:nrow(sk.un)){
    if(sk.un$skew[iii] < 0.2){
      polygon(rep(c(sk.un$start_gen[iii],sk.un$end_gen[iii]), each=2), c(0.66,0.99,0.99,0.66), border = NA,
              col=segment.colors[which( as.character(segment.index) == as.character(sk.un$skew[iii]) )])
      
      
    } else if(sk.un$skew[iii] > 0.8) {
      polygon(rep(c(sk.un$start_gen[iii],sk.un$end_gen[iii]), each=2), c(0,0.33,0.33,0), border = NA,
              col=segment.colors[which( as.character(segment.index) == as.character(sk.un$skew[iii]) )])
    } else{
      polygon(rep(c(sk.un$start_gen[iii],sk.un$end_gen[iii]), each=2), c(0.33,0.66,0.66,0.33), border = NA,
              col=segment.colors[which( as.character(segment.index) == as.character(sk.un$skew[iii]) )])
      polygon(rep(c(sk.un$start_gen[iii],sk.un$end_gen[iii]), each=2), c(0.33,0.66,0.66,0.33), border = 'red',
              col=NA, lwd=2)
    }}
  abline(v=ch.zoom.lines, col=1, lty=6)
  
  # T bases
  plot(pt.adj[(is.t & spec.points & chr.zoom)], -1 * log10(snp.dist[(is.t & spec.points & chr.zoom)]+1), pch=19, cex=1.1, 
       col=point.colors[(is.t & spec.points & chr.zoom)], xlim=example.chr.xlim, ylab='distance to nearest (log10)', 
       xlab='chromosome coordinate', ylim=c(-6,-2), cex.axis=1.5, xaxt='n')
  abline(v=ch.zoom.lines, col=1, lty=6)
  if(spec.name=='CAST'){
    legend('bottomleft', col='goldenrod', pch=19, legend=c('T>N'), bty='n', cex=1.3)
  }
}

# Clear Plot
dev.off()



#######  Panel 7h  ########
#######------------########
chrom.breaks <- round(scales::rescale(ch.lines, to=c(1,1000))[-c(1, (length(ch.lines)+1) )], digits = 0)

# Add chromosome X row to C3H which doesn't have segmentation signal
skt <- skew.mat
chx.row <- c(1,as.numeric(ch.vec['chrX']), sum(ch.vec[-length(ch.vec)]), sum(ch.vec), NA, 'chrX')
names(chx.row) <- c('start_chr', 'end_chr', 'start_gen', 'end_gen', 'skew', 'chrom')
for(x in names(skt)){
  tmp.mat <- data.frame(skt[[x]][['C3H']], stringsAsFactors = FALSE)
  tmp.mat$chrom <- as.character(tmp.mat$chrom)
  tmp.mat <- data.frame(rbind(tmp.mat, chx.row), stringsAsFactors = FALSE)
  skt[[x]][['C3H']] <- tmp.mat
}

# Create matrices of segmentation per haplotype
seg.mat <- lapply(names(skt[['MCB6A_T1']]), function(x){
  
  do.call(rbind, lapply(names(skt), function(z){
    skew.tmp <- skt[[z]][[x]]
    rescaled.points <- round(scales::rescale(as.numeric(c(as.numeric(skew.tmp$start_gen), as.numeric(skew.tmp$end_gen))), to=c(1,1000)), digits=0)
    rescale.start <- rescaled.points[1:nrow(skew.tmp)]; rescale.start[1] <- 0
    rescale.end <- rescaled.points[c((nrow(skew.tmp)+1):(nrow(skew.tmp) * 2))] - rescale.start
    
    return(as.numeric(unlist(lapply(1:nrow(skew.tmp), function(zz) {
      rep(skew.tmp[zz,'skew'], rescale.end[zz])
    }))))
  }))
}); names(seg.mat) <- names(skt[['MCB6A_T1']])

# Chromosome gap positions in heatmap
chrom.breaks <- round( scales::rescale(c(1, (ch.lines)[-c(1, (length(ch.lines)+1) )], sum(ch.vec)), to = c(1,1000)), digits=0)
chrom.breaks <- chrom.breaks[-c(1,length(chrom.breaks))]

# Plot Panel 7h heatmaps
for(i in names(seg.mat)){
  
  pheatmap(seg.mat[[i]], color=segment.colors, show_colnames = FALSE, show_rownames=FALSE, border_color = NA, main=i, 
           cluster_rows = FALSE, cluster_cols = FALSE, gaps_col  = chrom.breaks, cellheight = 15, cellwidth=0.65,
           legend_breaks = c(0.1,0.9), legend_labels = c("T>N", "A>N"))
}


