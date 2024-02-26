
############################
#####     Figure 2     #####  
############################
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicRanges)
library(ggplot2)
library(psych)



#######  Panel 2b  ########
#######------------########

# RDS file location
snv.new <- '/omics/groups/OE0538/internal/users/p281o/publications/single_cell_split/chip53_NewCall/filtered/PF1_SisterMutations.rds'

# Read in SRA table, split by sister name, remove negative controls
all.mutations <- readRDS(snv.new)
all.mutations <- all.mutations[all.mutations$FILTER]
rds.sis <- split(all.mutations, all.mutations$SampleName)
rds.sis <- rds.sis[-grep('Neg', names(rds.sis))]

# Index for referencing each second daughter of SNPs
snp.match <- seq(1, length(rds.sis), by=2)

#  Annotate T0 SNVs (by sharing between sisters)
for(i in snp.match){
  sister1 <- rds.sis[[i]]; sister2 <- rds.sis[[(i+1)]]
  elementMetadata(rds.sis[[i]])$TimePoint <- rep('T1', length(rds.sis[[i]]))
  elementMetadata(rds.sis[[(i+1)]])$TimePoint <- rep('T1', length(rds.sis[[(i+1)]]))
  rds.sis[[i]]$TimePoint[queryHits(findOverlaps(sister1, sister2))] <- 'T0'
  rds.sis[[(i+1)]]$TimePoint[subjectHits(findOverlaps(sister1, sister2))] <- 'T0'
}

# Count mutation types
mut.counts <- do.call(rbind, lapply(1:length(snp.match), function(x){
  uv.sis1 <- ((rds.sis[[snp.match[x]]]$REF=='C' & rds.sis[[snp.match[x]]]$ALT=='T') | (rds.sis[[snp.match[x]]]$REF=='G' & rds.sis[[snp.match[x]]]$ALT=='A'))
  uv.sis2 <- ((rds.sis[[(snp.match[x]+1)]]$REF=='C' & rds.sis[[(snp.match[x]+1)]]$ALT=='T') | (rds.sis[[(snp.match[x]+1)]]$REF=='G' & rds.sis[[(snp.match[x]+1)]]$ALT=='A'))
  ros.sis1 <-  ((rds.sis[[snp.match[x]]]$REF=='C' & rds.sis[[snp.match[x]]]$ALT=='A') | (rds.sis[[snp.match[x]]]$REF=='G' & rds.sis[[snp.match[x]]]$ALT=='T'))
  ros.sis2 <- ((rds.sis[[(snp.match[x]+1)]]$REF=='C' & rds.sis[[(snp.match[x]+1)]]$ALT=='A') | (rds.sis[[(snp.match[x]+1)]]$REF=='G' & rds.sis[[(snp.match[x]+1)]]$ALT=='T'))
  other.sis1 <- !(uv.sis1) & !(ros.sis1)
  other.sis2 <- !(uv.sis2) & !(ros.sis2)
  matrix(c(sum(uv.sis1), sum(uv.sis2), sum(ros.sis1), sum(ros.sis2), sum(other.sis1), sum(other.sis2)), ncol=2, byrow = TRUE)
}))

# Plot Panel 2b
par(mfrow=c(1,3)); y.limits <- c(500, 5000)
  plot(mut.counts[seq(2, nrow(mut.counts), by=3),], pch=19, col=c('dodgerblue4'), xlab='', ylab='', cex=1.7, xlim=y.limits, ylim=y.limits)
    abline(0,1,lty=3, col='black', lwd=1.5)
  plot(mut.counts[seq(1, nrow(mut.counts), by=3),], pch=19, col=c('firebrick4'), xlab='', ylab='', cex=1.7, xlim=y.limits, ylim=y.limits)
    abline(0,1,lty=3, col='black', lwd=1.5)
  plot(mut.counts[seq(3, nrow(mut.counts), by=3),], pch=19, col=c('darkgrey'), xlab='', ylab='', cex=1.7, xlim=y.limits, ylim=y.limits)
    abline(0,1,lty=3, col='black', lwd=1.5)

    

#######  Panel 2c  ########
#######------------########

# Mutational overlap between sisters and clones
ov.stats <- lapply(names(rds.sis), function(x){
  ov.cts <- do.call(cbind, lapply(names(rds.sis), function(y){
    res <- rep(0, length(rds.sis[[x]]))
    res[queryHits(findOverlaps(rds.sis[[x]], rds.sis[[y]], minoverlap = 0))] <- 1
    res
  })); colnames(ov.cts) <- names(rds.sis)
  ov.cts <- ov.cts[,grep(x, colnames(ov.cts), invert=TRUE)]
  ov.cts
}); names(ov.stats) <- names(rds.sis)

# Variable to assign sisters
overlap.pair <- names(ov.stats)[c(2,1,4,3,6,5,8,7,10,9,12,11,14,13)]

# Number of mutations unique (none), sister shared (sisters) or clone shared (clones)
overlap.stats <- do.call(rbind, lapply(1:length(ov.stats), function(x){
  message(x)
  tmp.df <- data.frame(ov.stats[[x]], stringsAsFactors = FALSE)
  
  c("none" = sum(rowSums(tmp.df)==0)/nrow(tmp.df), 
    "sisters"=sum(tmp.df[,overlap.pair[x]])/nrow(tmp.df), 
    "clones"=sum(rowSums(tmp.df[,grep(overlap.pair[x], colnames(tmp.df), invert=TRUE)])>0)/nrow(tmp.df))
  
})); overlap.stats <- round(overlap.stats, digits=3) * 100

# data frame for ggplot2
data <- data.frame( name=rep(colnames(overlap.stats), each=nrow(overlap.stats)), value=as.numeric(overlap.stats) )

# Plot Panel 2c
bar.cols <- rep('darkgrey',3)
p <- ggplot(data, aes(x=reorder(name, value), y=value)) + ggtitle('') + ylab('% SNV overlap') + xlab('') 
p + geom_bar(position = 'dodge', stat = 'summary', fun='mean', width=0.6, fill=bar.cols) +
  coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.background = element_blank(), axis.line = element_line(colour = "black"),
                       axis.text.y = element_text( size = 16)) +
  scale_y_continuous(limits = c(0, 100)) +
  geom_jitter(aes(x = name), shape = 21, fill='black', size=1.5, width=0.15)



#######  Panel 2d  ########
#######------------########

# Color vector for reference base
col.trinuc <- rep(c('dodgerblue1', 'black', 'firebrick2', 'darkgrey', 'limegreen', 'salmon'), each=16)

# Mutations signature set, always as reference C or T
mutations.of.interest <- c(paste0(rep('C',3),'_',c('A','G','T')), paste0(rep('T',3),'_',c('A','C','G')))

# All possible mutation types in each trinucleotide context
mut.vec <- paste0(c(rep(paste0(rep(c('A','C','G','T'), each=4),'C', rep(c('A','C','G','T'),4)),3),
                    rep(paste0(rep(c('A','C','G','T'), each=4),'T', rep(c('A','C','G','T'),4)),3)),'_',
                  rep(mutations.of.interest, each=16))

# List of mutations that are shared and unique to each sister
not.unique <- unlist(as(lapply(names(rds.sis)[snp.match], function(x){rds.sis[[x]][rds.sis[[x]]$TimePoint=='T0']}), "GRangesList"))
all.unique <- unlist(as(lapply(names(rds.sis), function(x){rds.sis[[x]][rds.sis[[x]]$TimePoint=='T1']}), "GRangesList"))
snp.list <- list('Unique'=all.unique, 'Shared'=not.unique)

# Plot Panel 2d
par(mfrow=c(2,1))
for(i in 1:length(snp.list)){
  
  # Result and Reverse Complement slicing vector
  results <- rep(0, length(mut.vec)); names(results) <- mut.vec 
  rev.comp.vec <- c('A', 'C', 'G', 'T'); names(rev.comp.vec) <- c('T', 'G', 'C', 'A')
  
  # Base pair change at site 
  current.gr <- snp.list[[i]]
  
  # Retrieve trinucleotide context
  current.bases <- resize(current.gr, fix = 'center', width=3)
  trinucs <- Biostrings::getSeq(x=BSgenome.Mmusculus.UCSC.mm10, current.bases)
  
  # Reverse complement G and A mutations to represent everything as a C or T reference
  trinucs[current.gr$REF=='G' | current.gr$REF=='A'] <- reverseComplement(trinucs[current.gr$REF=='G' | current.gr$REF=='A'])
  ref.base <- unlist(lapply(as.character(trinucs), function(x){unlist(strsplit(x,''))[[2]]}))
  alt.base <- current.gr$ALT; alt.base[current.gr$REF %in% c('G','A')] <- rev.comp.vec[alt.base[current.gr$REF %in% c('G','A')]]
  
  # Create names for mutations to match signature position
  trinucs <- paste0(as.character(trinucs),'_',ref.base,'_',alt.base)
  trinuc.ct <- table(trinucs)
  trinuc.ct <- trinuc.ct[(names(trinuc.ct) %in% names(results))]
  results[names(trinuc.ct)] <- as.numeric(trinuc.ct)
  
  # Plot Panel 2d
  barplot(results/sum(results), col=col.trinuc, las=2, cex.names = 0.5, ylab='frequency',
          main=paste0(names(snp.list)[i],' Mutations'), border=NA, xaxt='n')
}



#######  Panel 2e  ########
#######------------########

# Cosmic SBS signatures version 3.2
cosmic.sigs <- read.delim('/omics/groups/OE0538/internal/users/p281o/data_published/annotation/mm10/GRCm38_based/COSMIC_v3.2_SBS_mm10.txt',
                          sep='\t', stringsAsFactors = FALSE)

# Reorder mutations to match the mutation signature as above in Panel 2d
tri.changes <- unlist(lapply(cosmic.sigs$Type, function(x){
  bases.tmp <- unlist(strsplit(x,''))
  paste0( paste(bases.tmp[c(1,3,7)], collapse=''),'_',bases.tmp[3],'_', bases.tmp[5])
}))
rownames(cosmic.sigs) <- tri.changes; cosmic.sigs <- cosmic.sigs[,-1]
cosmic.sigs <- cosmic.sigs[mut.vec,]

# Plot Panel 2e
par(mfrow=c(2,1))
# SBS7a, UV
barplot(cosmic.sigs[,'SBS7a'], col=col.trinuc, las=2, cex.names=0.5, ylab='frequency', main='SBS7a', border=NA, xaxt='n')
# SBS18, ROS
barplot(cosmic.sigs[,'SBS18'], col=col.trinuc, las=2, cex.names=0.5, ylab='frequency', main='SBS18', border=NA, xaxt='n')



#######  Panels 2f,g,h  ########
#######-----------------########

# Selection vector for UV or oxidation damage (based on C>T or C>A mutations respectively)
uv.snvs <- lapply(names(rds.sis), function(x){ (rds.sis[[x]]$REF=='C' & rds.sis[[x]]$ALT=='T') | (rds.sis[[x]]$REF=='G' & rds.sis[[x]]$ALT=='A')})
ox.snvs <- lapply(names(rds.sis), function(x){ (rds.sis[[x]]$REF=='C' & rds.sis[[x]]$ALT=='A') | (rds.sis[[x]]$REF=='G' & rds.sis[[x]]$ALT=='T')})
names(uv.snvs) <- names(rds.sis); names(ox.snvs) <- names(rds.sis)

# Mutation categories of unique UV, Unique ROS and Shared 
new.uv <- unlist(lapply(names(rds.sis), function(x){  rds.sis[[x]]$VAF[!(rds.sis[[x]]$TimePoint=='T0') & uv.snvs[[x]]] }))
new.ox <- unlist(lapply(names(rds.sis), function(x){ rds.sis[[x]]$VAF[!(rds.sis[[x]]$TimePoint=='T0') & ox.snvs[[x]]] }))
shared.vafs <- unlist(lapply(names(rds.sis), function(x){ rds.sis[[x]]$VAF[(rds.sis[[x]]$TimePoint=='T0') ] }))

# VAF plot x-axis labels
tick.marks <- seq(0,1,by=0.25)

# Plot Panels 2f, 2g and 2h
par(mfrow=c(1,3))

# Shared SNV vafs
dens.fix <- density(shared.vafs, bw=0.025)  
plot(dens.fix, col='white', xlab='Variant Allele Frequency', xlim=c(-0.05,1.05), main='Shared sister SNVs', ylim=c(0,5), xaxt='n')
axis(side=1, at=tick.marks, labels=FALSE, tck=-0.06)
polygon(dens.fix, col = 'darkgrey', border=NA)

# Ox
dens.ox <- density(new.ox, bw=0.025)
plot(dens.ox, col='white', xlab='Variant Allele Frequency', xlim=c(-0.05,1.05), main='C>A sister unique', ylim=c(0,4), xaxt='n')
axis(side=1, at=tick.marks, labels=FALSE, tck=-0.06)
polygon(dens.ox, col = 'dodgerblue4', border=NA)

# UV
dens.uv <- density(new.uv, bw=0.025)
plot(dens.uv, col='white', xlab='Variant Allele Frequency', xlim=c(-0.05,1.05), main='C>T sister unique', ylim=c(0,3.5), xaxt='n')
axis(side=1, at=tick.marks, labels=FALSE, tck=-0.06)
polygon(dens.uv, col = 'firebrick3', border=NA)



#######  Panels 2i  ########
#######-------------########

# Calculate pearsons median skewness per sister and mutation type
pear.skew <-do.call(rbind, lapply(1:length(rds.sis), function(x){
  
  uv.stats <- data.frame(describe(rds.sis[[x]]$VAF[ rds.sis[[x]]$TimePoint=='T1' & uv.snvs[[x]]]))[,c('mean','median','sd')]
  ros.stats <- data.frame(describe(rds.sis[[x]]$VAF[ rds.sis[[x]]$TimePoint=='T1' & ox.snvs[[x]]]))[,c('mean','median','sd')]
  share.stats <- data.frame(describe(rds.sis[[x]]$VAF[ rds.sis[[x]]$TimePoint=='T0']))[,c('mean','median','sd')]
  
  uv.p <- (3 * (uv.stats[,'mean'] - uv.stats[,'median']))/uv.stats[,'sd']
  ros.p <- (3 * (ros.stats[,'mean'] - ros.stats[,'median']))/ros.stats[,'sd']
  share.p <- (3 * (share.stats[,'mean'] - share.stats[,'median']))/share.stats[,'sd']
  
  return(c('Share'=share.p, 'UV'=uv.p, 'ROS'=ros.p))
}))

##  Plot Panel 2i
boxplot(pear.skew, ylim=c(-1,1), col=c('darkgrey', 'firebrick3','dodgerblue4'), horizontal=TRUE)
abline(v=0, col=1, lwd=2, lty=3)


