############################
#####     Figure 5     #####  
############################
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
require(RColorBrewer)
require(randtests)
require(scales)

# Chromosome size vector
txdb <- keepStandardChromosomes(TxDb.Mmusculus.UCSC.mm10.knownGene)
ch.vec <- seqlengths(txdb); ch.vec <- ch.vec[!(names(ch.vec) %in% c('chrY', 'chrM'))]

# Read in mutations
snv.new <- '/omics/groups/OE0538/internal/users/p281o/publications/single_cell_split/chip53_NewCall/filtered/PF1_SisterMutations.rds'
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



#######  Panel 5a  ########
#######------------########

##  RL20 metric
##-------------
rl20.df <- do.call(rbind, lapply(names(rds.sis), function(x){
  
  # Acute UV
  snvs.uv <- rds.sis[[x]]
  snvs.uv <- snvs.uv[(snvs.uv$REF=='C'&snvs.uv$ALT=='T') | (snvs.uv$REF=='G'&snvs.uv$ALT=='A') & snvs.uv$TimePoint=='T1']
  uv.bin <- rep(0, length(snvs.uv)); uv.bin[snvs.uv$REF=='C'] <- 1
  uv.rle <- unlist(lapply(split(uv.bin, as.character(seqnames(snvs.uv))), function(z){unlist(rle(z)$lengths) }))
  uv.rle <- uv.rle[order(uv.rle, decreasing = TRUE)]
  
  # UV rl20 and pval
  uv.rl20 <- unlist(lapply(1:length(uv.rle), function(xx){
    if(sum(uv.rle[1:xx])  >=  (sum(uv.rle) * 0.2)){
      return(uv.rle[[xx]]) }}))[1]
  uv.pval <- -log10(runs.test(uv.bin, threshold = 0.5)$p.value)
  
  # Total ROS
  snvs.ros.tot <- rds.sis[[x]]
  snvs.ros.tot <- snvs.ros.tot[(snvs.ros.tot$REF=='C'&snvs.ros.tot$ALT=='A') | (snvs.ros.tot$REF=='G'&snvs.ros.tot$ALT=='T')]
  ros.tot.bin <- rep(0, length(snvs.ros.tot)); ros.tot.bin[snvs.ros.tot$REF=='C'] <- 1
  ros.tot.rle <- unlist(lapply(split(ros.tot.bin, as.character(seqnames(snvs.ros.tot))), function(z){unlist(rle(z)$lengths) }))
  ros.tot.rle <- ros.tot.rle[order(ros.tot.rle, decreasing = TRUE)]
  
  # ROS rl20 and pval
  ros.tot.rl20 <- unlist(lapply(1:length(ros.tot.rle), function(xx){
    if(sum(ros.tot.rle[1:xx])  >=  (sum(ros.tot.rle) * 0.2)){
      return(ros.tot.rle[[xx]]) }}))[1]
  ros.tot.pval <- -log10(runs.test(ros.tot.bin, threshold = 0.5)$p.value)
  
  results <- c('UV_rl20'=uv.rl20, 'UV_pval'=uv.pval, 'ROS_total_rl20'=ros.tot.rl20, 'ROS_total_pval'=ros.tot.pval)
  return(results)
}))

# Plotting colors
bar.cols <- brewer.pal(8, 'Dark2')[c(2,1,7)]

# Plot Panel 5a
layout.matrix <- matrix(c(rep(1, 9), rep(2,3)), nrow = 3, ncol = 4, byrow = FALSE)
layout(mat=layout.matrix)
par(mar=c(5,5,4,0))
plot(rl20.df[,'ROS_total_pval'], rl20.df[,'ROS_total_rl20'], col=alpha('cyan4', alpha = 0.6), pch=19, log='y',
     ylim=c(1,50), ylab='rl20', xlab='-log10(pval)', xlim=c(0,15), cex=2.7, cex.lab=2.1, cex.axis=1.8)
abline(v=-log10(0.05), col='blue', lty=2)
abline(v=4.5, h=5.5, col='blue')
legend('topright', col = c('firebrick4', 'cyan4'), bty='n', legend=c('UV', 'ROS'),
       border=NA, cex=1.8, pch=19)

par(mar=c(5,1,4,1))
plot(order(rl20.df[,'UV_pval']), rl20.df[,'UV_rl20'], log='y', ylim=c(1,50), pch=19, col=alpha('firebrick4', alpha=0.6), cex=2.7, xaxt='n',
     yaxt='n', xlab='', ylab='', xlim=c(0,15))
abline(h=5.5, col='blue')

# Clear Plot
dev.off()



#######  Panel 5b  ########
#######------------########

# Chromosome boundary lines
ch.lines <- c(1, unlist(lapply(1:length(ch.vec), function(x){ sum(ch.vec[1:x]) })))

# Plotting variables
clone <- 5
pt.sz <- 0.35; clr.trns <- 0.4

##  UV scatter
##------------
# Subset on example clone mutations, use C > T and G > A unique to sisters
sis1.uv <- rds.sis[[snp.match[clone]]]
sis1.uv <- sis1.uv[((sis1.uv$REF=='C' & sis1.uv$ALT=='T') | (sis1.uv$REF == 'G' & sis1.uv$ALT=='A')) & sis1.uv$TimePoint=='T1']

# Calculate distance to nearest mutation
sis1.uv.dist <- log10(values(distanceToNearest(sis1.uv))$distance + 1)
sis1.uv.dist[sis1.uv$REF=='G'] <- sis1.uv.dist[sis1.uv$REF=='G'] * -1

# Adjust graph position to juxtaposed chromosomes
sis1.uv.df <- data.frame(sis1.uv, stringsAsFactors = FALSE)
graph_pos_s1 <- unlist(lapply(1:nrow(sis1.uv.df), function(x) { 
  if(as.character(sis1.uv.df$seqnames[x])=='chr1') { as.numeric(sis1.uv.df$start[x])
  } else {  sum(as.numeric(ch.vec[1:(grep(as.character(sis1.uv.df$seqnames[x]),names(ch.vec))-1)])) + as.numeric(sis1.uv.df$start[x] )} }))

# Color points based on reference base at mutation
s1.uv.c <- sis1.uv$REF=='C'
plot.colors.s1 <- adjustcolor(c(rep('dodgerblue4', sum(!(s1.uv.c))), rep('goldenrod', sum((s1.uv.c)))), alpha.f = 0.25)

# Plot scatter for Sister 1
par(mar=c(0.4,4,4,4))
plot(c(graph_pos_s1[!(s1.uv.c)], graph_pos_s1[s1.uv.c]), c(sis1.uv.dist[!(s1.uv.c)], sis1.uv.dist[s1.uv.c]), ylim=c(-10,10),
     col=plot.colors.s1, ylab='distance to nearest (log10)', pch=19, cex=pt.sz, xlab='position', main='UV', 
     xlim=c(0,sum(ch.vec)), xaxt='n')

# Add chromosome lines
abline(v=ch.lines[-length(ch.lines)], col='black', lwd=0.3); abline(h=0, col=1, lty=3)
chr.name.pt <- unlist(lapply(2:length(ch.lines), function(x){ (ch.lines[x]-ch.lines[(x-1)])/2 + ch.lines[(x-1)] }))
mtext(gsub('chr','', names(ch.vec)), at=chr.name.pt, side=3)

##  ROS Scatter
##-------------

# Subset on example clone mutations, use C > T and G > A for ROS mutations
sis1.ros.gr <- rds.sis[[snp.match[clone]]]
sis1.ros.gr <- sis1.ros.gr[((sis1.ros.gr$REF=='C' & sis1.ros.gr$ALT=='A') | (sis1.ros.gr$REF == 'G' & sis1.ros.gr$ALT=='T')) & sis1.ros.gr$TimePoint=='T1']

# Calculate distance to nearest mutation
sis1.ros.dist <- log10(values(distanceToNearest(sis1.ros.gr))$distance + 1)
sis1.ros.dist[sis1.ros.gr$REF=='G'] <- sis1.ros.dist[sis1.ros.gr$REF=='G'] * -1

# Adjust graph position to juxtaposed chromosomes
sis1.ros.df <- data.frame(sis1.ros.gr, stringsAsFactors = FALSE)
graph_pos_s1 <- unlist(lapply(1:nrow(sis1.ros.df), function(x) { 
  if(as.character(sis1.ros.df$seqnames[x])=='chr1') { as.numeric(sis1.ros.df$start[x])
  } else {  sum(as.numeric(ch.vec[1:(grep(as.character(sis1.ros.df$seqnames[x]),names(ch.vec))-1)])) + as.numeric(sis1.ros.df$start[x] )} }))

s1.ros.c <- sis1.ros.gr$REF=='C'
plot.colors.s1 <- adjustcolor(c(rep('dodgerblue4', sum(!(s1.ros.c))), rep('goldenrod', sum((s1.ros.c)))), alpha.f = 0.25)

# Plot scatter for Sister 1
par(mar=c(0.4,4,4,4))
plot(c(graph_pos_s1[!(s1.ros.c)], graph_pos_s1[s1.ros.c]), c(sis1.ros.dist[!(s1.ros.c)], sis1.ros.dist[s1.ros.c]), ylim=c(-10,10),
     col=plot.colors.s1, ylab='distance to nearest (log10)', pch=19, cex=pt.sz, xlab='position', main='ROS', 
     xlim=c(0,sum(ch.vec)), xaxt='n')

# Add chromosome lines
abline(v=ch.lines[-length(ch.lines)], col='black', lwd=0.3); abline(h=0, col=1, lty=3)
mtext(gsub('chr','', names(ch.vec)), at=chr.name.pt, side=3)

# Clear Plot
dev.off()



#######  Panel 5c,d  ########
#######--------------########

# Create sliding genomic windows
gen.gr <- makeGRangesFromDataFrame(data.frame('chromosome'=names(ch.vec), 'start'=rep(1, length(ch.vec)), 'end'=as.numeric(ch.vec), stringsAsFactors = FALSE),
                                   start.field='start', end.field='end', seqnames.field = 'chromosome')
gt.win <- unlist(slidingWindows(gen.gr, width = 1*10^7, step = 10^5))
gt.center <- resize(gt.win, fix = 'center', width = 1)

# Graph X position
graph_pos_x <- unlist(lapply(1:length(gt.center), function(x) { 
  current.window <- gt.center[x]
  if(as.character(seqnames(current.window))=='chr1') { start(current.window)
  } else {  sum(as.numeric(ch.vec[1:(grep(as.character(seqnames(current.window)),names(ch.vec))-1)])) + start(current.window)} }))

# Calculate mutation skew
segment.skews <- data.frame(do.call(rbind, lapply(names(rds.sis), function(x){
  tmp.gr <- rds.sis[[x]]
  uv.snvs <- tmp.gr[((tmp.gr$REF=='C' & tmp.gr$ALT=='T') | (tmp.gr$REF=='G' & tmp.gr$ALT=='A')) & tmp.gr$TimePoint=='T1']
  ros.snvs <- tmp.gr[((tmp.gr$REF=='C' & tmp.gr$ALT=='A') | (tmp.gr$REF=='G' & tmp.gr$ALT=='T'))]
  
  # Find overlaps for both snp types, make empty result vectors
  win.uv <- findOverlaps(gt.win, uv.snvs); win.ros <- findOverlaps(gt.win, ros.snvs)
  skew.uv <- rep(NA, length(gt.win)); uv.count <- rep(0, length(gt.win))
  skew.ros <- rep(NA, length(gt.win)); ros.count <- rep(0, length(gt.win))
  
  uv.results <- rep(-1, length(uv.snvs)); uv.results[uv.snvs$REF=='C'] <- 1
  ros.results <- rep(-1, length(ros.snvs)); ros.results[ros.snvs$REF=='C'] <- 1
  
  uv.split <- split(uv.results[subjectHits(win.uv)], queryHits(win.uv))
  skew.uv[as.numeric(names(uv.split))] <- unlist(lapply(uv.split, mean))
  uv.count[as.numeric(names(lengths(uv.split)))] <- as.numeric(lengths(uv.split))
  
  ros.split <- split(ros.results[subjectHits(win.ros)], queryHits(win.ros))
  skew.ros[as.numeric(names(ros.split))] <- unlist(lapply(ros.split, mean))
  ros.count[as.numeric(names(lengths(ros.split)))] <- as.numeric(lengths(ros.split))
  
  return(cbind(skew.uv, uv.count, skew.ros, ros.count, rep(x, length(ros.count))))
  
})), stringsAsFactors = FALSE); colnames(segment.skews) <- c('Skew_UV', 'UV_count', 'Skew_ROS', 'ROS_count', 'Sample')


# Adapt bernoulli trial to include information about C > T mutations not likely from UV
# Assumption is that the C > T from background is proportional to the C > A signal
ct.ratio <- unlist(lapply(names(rds.sis), function(x){
  tmp.gr <- rds.sis[[x]]
  t0.muts <- tmp.gr[tmp.gr$TimePoint=='T0']
  ct.muts <- t0.muts[(t0.muts$REF=='C'&t0.muts$ALT=='T') | (t0.muts$REF=='G'&t0.muts$ALT=='A')]
  ca.muts <- t0.muts[(t0.muts$REF=='C'&t0.muts$ALT=='A') | (t0.muts$REF=='G'&t0.muts$ALT=='T')]
  length(ct.muts)/(length(ct.muts)+length(ca.muts))
}))

# Proportion of CT SNVs likely coming from Non-UV 
ct.bg <- do.call(rbind, lapply(names(rds.sis), function(x){
  tmp.gr <- rds.sis[[x]]
  t0.muts <- tmp.gr[tmp.gr$TimePoint=='T1']
  ct.muts <- t0.muts[(t0.muts$REF=='C'&t0.muts$ALT=='T') | (t0.muts$REF=='G'&t0.muts$ALT=='A')]
  ca.muts <- t0.muts[(t0.muts$REF=='C'&t0.muts$ALT=='A') | (t0.muts$REF=='G'&t0.muts$ALT=='T')]
  bg.ct <- round(length(ca.muts) * mean(ct.ratio), digits=0)
  prop.tot.ct <- round(((length(ca.muts) * mean(ct.ratio))/(length(ct.muts)+length(ca.muts))), digits=2)
  c('BG_CT'=bg.ct, 'Prop_Tot_CT'=prop.tot.ct)
})); rownames(ct.bg) <- names(rds.sis)

# Create ratios trained on proportion of background mutation
mixes <- matrix(rep(c(0,0.5,0.5,1), length(rds.sis)), nrow=length(rds.sis), byrow = TRUE)
mixes[,1] <- mixes[,1]+(0.5 * ct.bg[, 'Prop_Tot_CT']) # Half of the error
mixes[,4] <- mixes[,4]-(0.5 * ct.bg[, 'Prop_Tot_CT']) # Half of the error
rownames(mixes) <- names(rds.sis)

# create skews based on background mutation
set.seed(101)
bins.per.sample <- table(segment.skews$Sample[as.numeric(segment.skews[, 'UV_count'])>= 10])
uv.exp <- scales::rescale(x = unlist(lapply(rownames(mixes), function(x){
  snvs.per.bin <- segment.skews$UV_count[(segment.skews$Sample %in% x) & as.numeric(segment.skews$UV_count) >=10]
  mix.props <- sample(mixes[x,], as.numeric(bins.per.sample[x]), replace=TRUE)
  unlist(lapply(1:length(snvs.per.bin), function(xx){
    mean(rbinom(snvs.per.bin[xx],1, mix.props[xx] ))
  }))
})), to = c(-1,1))

# Sample UV model 100 times to make a smooth plot
uv.smth <- scales::rescale(x = unlist(lapply(rownames(mixes), function(x){
  snvs.per.bin <- as.numeric(segment.skews$UV_count[(segment.skews$Sample %in% x) & as.numeric(segment.skews$UV_count) >=10])
  mix.props <- sample(mixes[x,], as.numeric(bins.per.sample[x]), replace=TRUE)
  unlist(lapply(1:length(snvs.per.bin), function(xx){
    mean(rbinom(rep(snvs.per.bin[xx],100),1, rep(mix.props[xx],100) ))
  }))
})), to = c(-1,1))

# Sample ROS model
ros.samplings <- as.numeric(segment.skews[,"ROS_count"][as.numeric(segment.skews[,'ROS_count']) >= 10])
ros.expected <- scales::rescale(x = unlist(lapply(1:length(ros.samplings), function(x){ mean(rbinom(ros.samplings[x],1, 0.5))})), to = c(-1,1))
ros.smth <- scales::rescale(x = unlist(lapply(rep(1:length(ros.samplings),100), function(x){ mean(rbinom(ros.samplings[x],1, 0.5))})), to = c(-1,1))

# Density of UV and ROS models
dens.uv <- density(uv.smth, bw=0.1)
dens.ros <- density(ros.smth, bw=0.1)

# Plot Panels 5c,d
label.size <- 1.75
par(mfrow=c(2,2))

plot(dens.uv, col='white', xlab='skew', xlim=c(-1.3,1.3), main='mutation phasing', cex.axis=label.size, cex.lab=label.size, ylab='density')
polygon(dens.uv, col = alpha('darkgrey', alpha = 0.8), border=NA)

hist(as.numeric(segment.skews[as.numeric(segment.skews[,'UV_count']) >= 10,'Skew_UV']), breaks=20, xlab='skew', ylab='density', main='UV', col=alpha('firebrick4', alpha=0.6), 
     cex.axis=label.size, cex.lab=label.size, freq=FALSE, border=NA)
hist(uv.exp, breaks=20, main="", xlab='', ylab='', col=alpha('darkgrey', alpha=0.7) , cex.axis=label.size, cex.lab=label.size, freq=FALSE, border=NA, add=TRUE)

plot(dens.ros, col='white', xlab='skew', ylab='density', xlim=c(-1.3,1.3), main='no mutation phasing', cex.axis=label.size, cex.lab=label.size)
polygon(dens.ros, col = alpha('darkgrey', alpha=0.8), border=NA)

hist(as.numeric(segment.skews[as.numeric(segment.skews[,'ROS_count']) >= 10,'Skew_ROS']), breaks=10, xlab='skew', ylab='density', main='ROS', col=alpha('cyan4', alpha=0.5), cex.axis=label.size, 
     cex.lab=label.size, freq=FALSE, border=NA)
hist(ros.expected, breaks=10, main="", xlab='', ylab='', col=alpha('darkgrey', alpha=0.7), cex.axis=label.size, 
     cex.lab=label.size, freq = FALSE, border=NA, add=TRUE)

# Clear Plot
dev.off()



#######  Panel 5e  ########
#######------------########
# Line points for comparisons to the opposite distribution
uv.ros.comp <- qqplot(as.numeric(segment.skews[as.numeric(segment.skews[,'UV_count']) >= 10,'Skew_UV']),  
                      ros.expected, plot=FALSE)
ros.uv.comp <- qqplot(as.numeric(segment.skews[as.numeric(segment.skews[,'ROS_count']) >= 10,'Skew_ROS']),  
                      as.numeric(segment.skews[as.numeric(segment.skews[,'UV_count']) >= 10,'Skew_UV']), plot=FALSE)

# Plot Panel 5e
par(mfrow=c(2,1))
qqplot(uv.exp, as.numeric(segment.skews[as.numeric(segment.skews[,'UV_count']) >= 10,'Skew_UV']),  pch=19, cex=0.8, 
       col=alpha('firebrick4', alpha=0.2), xlab='theoretical quartiles', ylab='sample quartiles', main='UV')
  points(uv.ros.comp, col=alpha('darkgrey', alpha=0.1), pch=19, cex=0.8)
  abline(0,1, col= 'black', lwd=2, lty=6)

qqplot(as.numeric(segment.skews[as.numeric(segment.skews[,'ROS_count']) >= 10,'Skew_ROS']),  ros.expected, pch=19, cex=0.8, 
       col=alpha('cyan4', alpha=0.2), xlab='theoretical quartiles', ylab='sample quartiles', main='ROS')
  points(ros.uv.comp, col=alpha('darkgrey', alpha=0.1), pch=19, cex=0.8)
  abline(0,1, col='black', lwd=2, lty=6)
