
############################
#####     Figure 4     #####  
############################
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
library(RColorBrewer)
library(Repitools)
library(regioneR)
library(ggplot2)
library(QuasR)
library(Gviz)

# RDS file location
snv.new <- '/omics/groups/OE0538/internal/users/p281o/publications/single_cell_split/chip53_NewCall/filtered/PF1_SisterMutations.rds'
all.mutations <- readRDS(snv.new); all.mutations <- all.mutations[all.mutations$FILTER]

# Split into single samples
sister.split <- split(all.mutations, all.mutations$SampleName)
sister.split <- sister.split[grep('Neg', names(sister.split), invert=TRUE)]

# Index for referencing mitotic sister of each pair
snp.match <- seq(1, length(sister.split), by=2)

# Annotate T0 mutations by shared representation in mitotic sisters
for(i in snp.match){
  sister1 <- sister.split[[i]]; sister2 <- sister.split[[(i+1)]]
  elementMetadata(sister.split[[i]])$TimePoint <- rep('T1', length(sister.split[[i]]))
  elementMetadata(sister.split[[(i+1)]])$TimePoint <- rep('T1', length(sister.split[[(i+1)]]))
  sister.split[[i]]$TimePoint[queryHits(findOverlaps(sister1, sister2))] <- 'T0'
  sister.split[[(i+1)]]$TimePoint[subjectHits(findOverlaps(sister1, sister2))] <- 'T0'
}

## Add trinucleotide context to speed up rate calculations
for(tri in 1: length(sister.split)){
  elementMetadata(sister.split[[tri]])$TriN <- as.character(as.character(BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm10, resize(sister.split[[tri]], fix='center', width=3))))
}

# All unique mutations
all.unique <- unlist(as(lapply(snp.match, function(x){
  sis1 <- sister.split[[x]]; sis1 <- sis1[sis1$TimePoint=='T1']
  sis2 <- sister.split[[(x+1)]]
  unlist(as(list(sis1, sis2), 'GRangesList'))
}),'GRangesList'))
names(all.unique) <- 1:length(all.unique)

#  RNA data
rna.cts <- readRDS('/omics/groups/OE0538/internal/users/p281o/publications/single_cell_split/Genic_and_RNA/PF1_RNA_counts.rds')
gene.dist <- 1000; gene.width <- 1000
rna.cts <- rna.cts[(values(distanceToNearest(rna.cts))$distance >= gene.dist) & width(rna.cts) >= gene.width]

# Expression binning, take average tpm from 3 replicates
ps <- 0.1 # pseudocount
avg.tpm <- rowMeans(data.frame(rna.cts, stringsAsFactors=FALSE)[,c('RNA_rep1', 'RNA_rep2', 'RNA_rep3')])

# Set minimum expression at 0
some.exp <- log2(avg.tpm + ps) >= 0
mean.exp <- log2(avg.tpm + ps)

# Bin genes by expression in TPM
exp.bins <- rep(1, length(rna.cts))
expressed.bins <- cut(mean.exp[some.exp], breaks=quantile(mean.exp[some.exp], probs=seq(0,1,1/3)), labels = FALSE)
expressed.bins[is.na(expressed.bins)] <- 1
exp.bins[some.exp] <- expressed.bins + 1
elementMetadata(rna.cts)$Exp_Bin <- exp.bins

# Colors for adding bin annotation below
gene.colors <- brewer.pal(9, 'Blues')[c(2,4,7,9)]; names(gene.colors) <- 1:4



#######  Panel 4a  ########
#######------------########

# tile Genome
txdb <- keepStandardChromosomes(TxDb.Mmusculus.UCSC.mm10.knownGene)
ch.vec <- seqlengths(txdb); ch.vec <- ch.vec[!(names(ch.vec) %in% c('chrY', 'chrM'))]
chr.win <- tileGenome(seqlengths = ch.vec, tilewidth = 10^5, cut.last.tile.in.chrom = TRUE)

# Bam Files for RNA and DNA
rna.bams <- read.delim('/omics/groups/OE0538/internal/users/p281o/projects/lce_mechanism/sequencing/seqID_27232/samples.txt', sep='\t', stringsAsFactors = FALSE)
rna.bam <- rna.bams[3,'FileName']
atac.bams <- read.delim('/omics/groups/OE0538/internal/users/p281o/projects/lce_mechanism/sequencing/seqID_30724/samples_quasr.txt', sep='\t', stringsAsFactors = FALSE)
atac.bam <- atac.bams[1,'FileName']

# Using Chek2 window from first verion of panel
gene.window <- chr.win
current.window <- subjectHits(findOverlaps(rna.cts[grep('Chek2', rna.cts$GENENAME)], gene.window))

# RNA track
chr <- as.character(seqnames(gene.window[current.window])); gen <- 'mm10'; names(gen) <- chr
rna.track <- DataTrack(range = rna.bam, genome = "mm10", type = "l", name = "Coverage", window = -1, chromosome = chr)

# Gene track for window
grtrack <- GeneRegionTrack(txdb, genome = gen, chromosome = chr, name = "UCSC genes",stacking = 'squish', col='dodgerblue4', 
                           fill='dodgerblue4', background.title='dodgerblue4')

# Ideogram track
itrack <- IdeogramTrack(genome = gen, chromosome = chr)

# Genome axis track
gtrack <- GenomeAxisTrack()

# RNA track
rna.track <- AlignmentsTrack(range=rna.bam, start = start(gene.window[current.window]), end=end(gene.window[current.window]), isPaired=TRUE,
                             name='RNA', type='coverage', col='chocolate4', fill='chocolate4', coverageHeight=0.01, background.title='chocolate4', 
                             ylim=c(0,30))
# ATAC track
atac.track <- AlignmentsTrack(range=atac.bam, start = start(gene.window[current.window]), end=end(gene.window[current.window]), isPaired=TRUE,
                              name='ATAC', type='coverage', col='forestgreen', fill='forestgreen', coverageHeight=0.01, 
                              background.title='forestgreen', ylim=c(0,1000))

# Gene Model track
gene.colors <- brewer.pal(9, 'Blues')[c(2,4,7,9)]; names(gene.colors) <- 1:4
win.cols <- gene.colors[as.numeric(rna.cts$Exp_Bin[queryHits(findOverlaps(rna.cts, gene.window[current.window] ))])]
gmTrack <- AnnotationTrack(rna.cts, start(gene.window[current.window]), end=end(gene.window[current.window]), 
                           chromosome = chr, id=rna.cts$GENENAME, background.title='darkslategray4', width=2,
                           genome = "mm10", name = "Gene model", fill= win.cols , cex=0.9)

# Mutation track
gene.mutations <- subsetByOverlaps(all.unique, gene.window[current.window])
mut.track <- AnnotationTrack(gene.mutations, start(gene.window[current.window]), end=end(gene.window[current.window]), 
                             chromosome = chr, background.title='black', width=2,
                             genome = "mm10", name = "Mut", fill='black', cex=0.9)

# Plot Panel 4a
plotTracks(list(itrack, gtrack, rna.track, atac.track, grtrack, gmTrack, mut.track), from = start(gene.window[current.window]), 
           to = end(gene.window[current.window]), sizes = c(0.05,0.1,0.25,0.25,0.125,0.125, 0.1), featureAnnotation='id', cex.title=1 )



#######  Panel 4b  ########
#######------------########
layout(matrix(c(rep(1,3), rep(2,6)), nrow=3, byrow = TRUE))
par(mar=c(1,4,2,1))
  hist(mean.exp[(values(distanceToNearest(rna.cts))$distance >= 1000)], main='', xlab='', breaks=40, col='darkslategray4', border=NA, xaxt='n', ylim=c(5000,10000))
    abline(v=c(0,quantile(mean.exp[some.exp], probs=seq(0,1,1/3))), col=2, lty=2, lwd=2)
  hist(mean.exp[(values(distanceToNearest(rna.cts))$distance >= 1000)], main='', xlab='', breaks=40, col='darkslategray4', border=NA, xaxt='n', ylim=c(-100,1500))
    axis(3)
  abline(v=c(0,quantile(mean.exp[some.exp], probs=seq(0,1,1/3))), col=2, lty=2, lwd=2)
  rect(c(min(mean.exp), c(0,quantile(mean.exp[some.exp], probs=seq(0,1,1/3)))[-c(1,5)]), 
     rep(-100,4), quantile(mean.exp[some.exp], probs=seq(0,1,1/3)), rep(-5,4), col=gene.colors, border=NA, xpd=TRUE)

  

#######  Panel 4c  ########
#######------------########

# Tile Genome, retain regions with 95% mappability
gt <- tileGenome(seqlengths = ch.vec, tilewidth = 10^4, cut.last.tile.in.chrom = TRUE)
gt.mapp <- mappabilityCalc(gt, BSgenome.Mmusculus.UCSC.mm10)
analyzeable.gen <- gt[gt.mapp >= 0.95]

#  Genomic Trinucleotide Frequencies
bg.tri <- colSums(trinucleotideFrequency(BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm10, analyzeable.gen)))
# Relevant trinucleotides are at C or G bases (both UV and ROS happen at these bases)
bg.tri <- bg.tri[unlist(lapply(names(bg.tri), function(x){ (unlist(strsplit(x,''))[2] %in% c('C','G')) }))]

# Minimum coverage variable
min.cov <- 10

# Genomic rates for UV, all ROS, T1 and T0 ROS
snp.rt.all <- do.call(rbind, lapply(1:length(sister.split), function(x){
  all.snps <- sister.split[[x]]; all.snps <- all.snps[all.snps$coverage >= min.cov] 
  
  # Mutation identities
  all.uv <- all.snps[(((all.snps$REF=='C' & all.snps$ALT=='T') | (all.snps$REF=='G' & all.snps$ALT=='A')) & all.snps$TimePoint=='T1')]
  all.ros <- all.snps[(((all.snps$REF=='C' & all.snps$ALT=='A') | (all.snps$REF=='G' & all.snps$ALT=='T')))]

  # UV rate
  uv <- colSums(trinucleotideFrequency(BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm10, resize(all.uv, fix='center', width=3))))
  uv <- uv[ names(bg.tri) ]
  uv.res <-(weighted.mean(uv/bg.tri, bg.tri/sum(bg.tri)) * 10^6)
  
  # ROS all rate
  ros <- colSums(trinucleotideFrequency(BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm10, resize(all.ros, fix='center', width=3))))
  ros <- ros[ names(bg.tri) ]
  ros.res <-(weighted.mean(ros/bg.tri, bg.tri/sum(bg.tri)) * 10^6)
  
  return(c('UV'=uv.res , 'ROS_All'=ros.res))
})); rownames(snp.rt.all) <- names(sister.split)


# Rates in genes figure
gene.dist <- 1000; gene.width <- 1000
min.cov <- 10

# Genic Trinucleotide frequencies
relevant.genes <- rna.cts[(values(distanceToNearest(rna.cts))$distance >= gene.dist) & width(rna.cts) >= gene.width]

# Trinucleotide counts at C and G bases in relevant genic bin
gene.tri <- lapply(unique(exp.bins), function(eb){
  genes.in.bin <- relevant.genes[relevant.genes$Exp_Bin==eb]
  genes.res <- colSums(trinucleotideFrequency(BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm10, genes.in.bin)))
  genes.res <- genes.res[unlist(lapply(names(genes.res), function(ti){ (unlist(strsplit(ti,''))[2] %in% c('C','G')) }))]
  genes.res
}); names(gene.tri) <- unique(exp.bins)
total.genic.bg <- colSums(do.call(rbind, gene.tri))/sum(colSums(do.call(rbind, gene.tri)))

# Calculate mutation rates for UV and ROS in genic bins
genic.rates <- do.call(rbind, lapply(names(gene.tri), function(x){
  message(x)
  
  current.genes <- relevant.genes[relevant.genes$Exp_Bin==as.numeric(x)]
  
  do.call(rbind, lapply(names(sister.split), function(xx){
    
    # use both pairs, remove redundant snps, filter on at least 20 reads
    sister <- sister.split[[xx]]; sister <- sister[sister$coverage >= min.cov]
    gene.ov <- findOverlaps(sister, current.genes, ignore.strand=TRUE)
    is.in.gene <- rep(FALSE, length(sister)); is.in.gene[queryHits(gene.ov)] <- TRUE
    
    # Mutation identities
    all.uv <- (((sister$REF=='C' & sister$ALT=='T') | (sister$REF=='G' & sister$ALT=='A')) & sister$TimePoint=='T1')
    all.ros <- (((sister$REF=='C' & sister$ALT=='A') | (sister$REF=='G' & sister$ALT=='T')))

    # UV rate and delta
    uv.res <- rep(0, length(gene.tri[[x]])); names(uv.res) <- names(gene.tri[[x]])
    uv.tri <- table(sister$TriN[is.in.gene & all.uv]); uv.tri <- uv.tri[names(uv.tri) %in% names(uv.res)]
    uv.res[names(uv.tri)] <- as.numeric(uv.tri)
    uv.rate <-(weighted.mean(uv.res/gene.tri[[x]], total.genic.bg) * 10^6)
    uv.delta <- uv.rate/snp.rt.all[xx,'UV']
    
    # ROS rate and delta
    ros.res <- rep(0, length(gene.tri[[x]])); names(ros.res) <- names(gene.tri[[x]])
    ros.tri <- table(sister$TriN[is.in.gene & all.ros]); ros.tri <- ros.tri[names(ros.tri) %in% names(ros.res)]
    ros.res[names(ros.tri)] <- as.numeric(ros.tri)
    ros.rate <-(weighted.mean(ros.res/gene.tri[[x]], total.genic.bg) * 10^6)
    ros.delta <- ros.rate/snp.rt.all[xx,'ROS_All']
    
    c('UV_rate'=uv.rate, "UV_delta"=uv.delta, 'ROS_rate'=ros.rate, 'ROS_Delta'=ros.delta)
  })) 
}))

# Combine rates, add TPM strata tag
all.rates <- c(genic.rates[,'UV_delta'], genic.rates[,'ROS_Delta'])
strata <- rep(c(seq(1,8,by=2), seq(2,8,by=2)), each=14)

# Barplot with error bars
uv.mean <- unlist(lapply(split((genic.rates[,'UV_delta']), rep(1:4, each=14)), mean))
uv.sd <- unlist(lapply(split((genic.rates[,'UV_delta']), rep(1:4, each=14)), sd))
ros.mean <- unlist(lapply(split((genic.rates[,'ROS_Delta']), rep(1:4, each=14)), mean))
ros.sd <- unlist(lapply(split((genic.rates[,'ROS_Delta']), rep(1:4, each=14)), sd))

# Individual points
d.pt <- data.frame('source'=rep(c('UV','ROS'), each=nrow(genic.rates)), 'bin'=rep(1:4, each=14), 
                   'rate'=c(genic.rates[,'UV_delta'], genic.rates[,'ROS_Delta']) )

# Mean by sample, bin and source
dmg.df <- data.frame('damage'=rep(c('UV','ROS'), each=4), 'bin'=as.factor(rep(1:4,2)),
                     'mean'=c(uv.mean, ros.mean),'sd'=c(uv.sd,ros.sd))

# Plot Panel 4c
p <- ggplot(dmg.df, aes(x=bin, y=mean, fill=damage)) 

p + geom_bar(stat="identity", color=NA, position=position_dodge()) +
  geom_point(data=d.pt, aes(x=bin,y=rate,fill=source), position = position_jitterdodge(jitter.width=0.25),size=2,
             color=grey.colors(10)[3] ) +
  geom_errorbar(aes(ymin=mean - 2*sd, ymax=mean + 2*sd), width=.2, position=position_dodge(.9)) +
  labs(title="", x="TPM bin", y = "rate/background") + theme_classic() + scale_fill_manual(values=c('turquoise3','sienna3')) +
  scale_y_continuous(name ="genic rate/background rate\n(mu/mb)", limits=c(0,1.25)) +
  theme(axis.text.x = element_text(face="bold", color="black", size=14), axis.text.y = element_text(face="bold", color="black", size=14),
        axis.title=element_text(size=14,face="bold"))



#######  Panel 4e  ########
#######------------########

# Iterate by expression bin
genic.exp.all <- do.call(rbind, lapply(names(gene.tri), function(z){
  
  # Genes in current bin
  bin.genes <- relevant.genes[relevant.genes$Exp_Bin == z]
  
  # Loop by map
  genic.snps <- do.call(rbind, lapply(names(sister.split), function(i){
    
    # All mutations for current pair, meeting minimum coverage threshold
    all.snps <- sister.split[[i]]
    all.snps <- all.snps[all.snps$coverage >= min.cov]
    
    # All overlapping mutations in genes from bin
    snp.gene.ov <- findOverlaps(all.snps, bin.genes)
    snp.ov <- all.snps[queryHits(snp.gene.ov)]; gene.ov <- bin.genes[subjectHits(snp.gene.ov)]
    
    # Selective vectors for UV, ROS, and shared mutations
    is.uv <- ((snp.ov$REF=='C' & snp.ov$ALT=='T') | (snp.ov$REF=='G' & snp.ov$ALT=='A')) & snp.ov$TimePoint=='T1'
    is.ros <- ((snp.ov$REF=='C' & snp.ov$ALT=='A') | (snp.ov$REF=='G' & snp.ov$ALT=='T'))
    is.t0 <- snp.ov$TimePoint=='T0'
    
    # Template UV rate
    tuv.res <- rep(0, length(gene.tri[[z]])); names(tuv.res) <- names(gene.tri[[z]])
    temp.uv <- (as.character(strand(gene.ov)) == '+' & is.uv & snp.ov$REF=='G') |
      (as.character(strand(gene.ov)) == '-' & is.uv & snp.ov$REF=='C')
    uv.temp.tri <- table(snp.ov$TriN[temp.uv])
    tuv.res[names(uv.temp.tri)] <- as.numeric(uv.temp.tri)
    uv.temp.rate <-(weighted.mean(tuv.res/gene.tri[[z]], bg.tri/sum(bg.tri)) * 10^6)

    # Non-template UV rate
    ntuv.res <- rep(0, length(gene.tri[[z]])); names(ntuv.res) <- names(gene.tri[[z]])
    nontemp.uv <- (as.character(strand(gene.ov)) == '+' & is.uv & snp.ov$REF=='C') |
      (as.character(strand(gene.ov)) == '-' & is.uv & snp.ov$REF=='G')
    uv.nontemp.tri <- table(snp.ov$TriN[nontemp.uv])
    ntuv.res[names(uv.nontemp.tri)] <- as.numeric(uv.nontemp.tri)
    uv.nontemp.rate <-(weighted.mean(ntuv.res/gene.tri[[z]], bg.tri/sum(bg.tri)) * 10^6)

    # Template ROS rate (all)
    tros.res <- rep(0, length(gene.tri[[z]])); names(tros.res) <- names(gene.tri[[z]])
    temp.ros <- (as.character(strand(gene.ov)) == '+' & is.ros & snp.ov$REF=='C') |
      (as.character(strand(gene.ov)) == '-' & is.ros & snp.ov$REF=='G')
    ros.temp.tri <- table(snp.ov$TriN[temp.ros])
    tros.res[names(ros.temp.tri)] <- as.numeric(ros.temp.tri)
    ros.temp.rate <-(weighted.mean(tros.res/gene.tri[[z]], bg.tri/sum(bg.tri)) * 10^6)
    
    # Non-template ROS rate (all)
    ntros.res <- rep(0, length(gene.tri[[z]])); names(ntros.res) <- names(gene.tri[[z]])
    nontemp.ros <- (as.character(strand(gene.ov)) == '+' & is.ros & snp.ov$REF=='G' ) |
      (as.character(strand(gene.ov)) == '-' & is.ros & snp.ov$REF=='C' )
    ros.nontemp.tri <- table(snp.ov$TriN[nontemp.ros])
    ntros.res[names(ros.nontemp.tri)] <- as.numeric(ros.nontemp.tri)
    ros.nontemp.rate <-(weighted.mean(ntros.res/gene.tri[[z]], bg.tri/sum(bg.tri)) * 10^6)
 

    return(list('Template_UV' = (uv.temp.rate/snp.rt.all[i,'UV']) * 2, 
                'NonTemplate_UV' = (uv.nontemp.rate/snp.rt.all[i,'UV']) * 2,
                'Template_ROS' = (ros.temp.rate/snp.rt.all[i,'ROS_All']) * 2, 
                'NonTemplate_ROS'= (ros.nontemp.rate/snp.rt.all[i,'ROS_All']) * 2, 
                'Num_Genes' = length(bin.genes)))
  }))
  
  return(genic.snps)
}))

# Plot Panel 4e
par(mfrow=c(1,2))
template.color <- alpha('black', alpha=0.6)
nontemplate.color <- alpha('darkgrey', alpha=0.55)

# UV mutations
stripchart(unlist(genic.exp.all[,'Template_UV']) ~ rep(1:4, each=14), method='jitter', vertical=TRUE, pch=19, col=template.color,
           main='UV mutations', ylab='mutations/mb (temp/non-temp)', xaxt='n', ylim=c(0,1.5), cex=1.5)
stripchart( unlist(genic.exp.all[,'NonTemplate_UV']) ~ rep(1:4, each=14), method='jitter', vertical=TRUE, pch=19, col=nontemplate.color,
            main='', ylab='mutations/mb (temp/non-temp)', xaxt='n', add=TRUE, cex=1.5)
abline(h=c(1), col='black', lty=3)

# ROS mutations
stripchart(unlist(genic.exp.all[,'Template_ROS']) ~ rep(1:4, each=14), method='jitter', vertical=TRUE, pch=19, col=template.color,
           main='ROS mutations', ylab='mutations/mb (temp/non-temp)', xaxt='n', ylim=c(0,1.5), cex=1.5)
stripchart( unlist(genic.exp.all[,'NonTemplate_ROS']) ~ rep(1:4, each=14), method='jitter', vertical=TRUE, pch=19, col=nontemplate.color,
            main='', ylab='mutations/mb (temp/non-temp)', xaxt='n', add=TRUE, cex=1.5)
abline(h=c(1), col='black', lty=3)



#######  Panel 4f  ########
#######------------########

# ATAC peaks file
atac.peaks <- readRDS('/omics/groups/OE0538/internal/users/p281o/publications/single_cell_split/ATAC/ATAC_Peaks_PF1.rds')

# Trinucleotide frequencies centered on C and G bases in ATAC peaks
peak.tri <- colSums(trinucleotideFrequency(BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm10, atac.peaks)))
peak.tri <- peak.tri[unlist(lapply(names(peak.tri), function(x){ (unlist(strsplit(x,''))[2] %in% c('C','G')) }))]

# Selective vector for ROS or UV mutations
is.uv <- (((all.unique$REF=='C' & all.unique$ALT=='T') | (all.unique$REF=='G' & all.unique$ALT=='A')) & all.unique$TimePoint=='T1')
is.ros <- (((all.unique$REF=='C' & all.unique$ALT=='A') | (all.unique$REF=='G' & all.unique$ALT=='T')))

# Tag to identify sister pairs
sister.tags <- unique(substr( unique(all.unique$SampleName), 1, nchar(unique(all.unique$SampleName)) - 1 ))

# Empty vector to keep result structure
res <- rep(0, length(bg.tri)); names(res) <- names(bg.tri)

# Iterate through sister pairs
peak.rates <- do.call(rbind, lapply(sister.tags, function(x){
  
  # Select mutations for sister pair
  uv.muts <- all.unique[(is.uv & grepl(x, all.unique$SampleName) )]
  ros.muts <- all.unique[(is.ros & grepl(x, all.unique$SampleName) )]
  
  # Overlap with ATAC peaks
  atac.uv <- findOverlaps(uv.muts, atac.peaks)
  atac.ros <- findOverlaps(ros.muts, atac.peaks)
  
  # Genomic rates
  uv.bg <- res; uv.bg[names(table(uv.muts$TriN))] <- as.numeric(table(uv.muts$TriN))
  uv.bg <- (weighted.mean(uv.bg/bg.tri, bg.tri/sum(bg.tri))*10^6)
  ros.bg <- res; ros.bg[names(table(ros.muts$TriN))] <- as.numeric(table(ros.muts$TriN))
  ros.bg <- (weighted.mean(ros.bg/bg.tri, bg.tri/sum(bg.tri))*10^6)
  
  # ATAC peak rates
  uv.pk <- res; uv.pk[names(table(uv.muts$TriN[queryHits(atac.uv)]))] <- as.numeric(table(uv.muts$TriN[queryHits(atac.uv)]))
  uv.pk <- (weighted.mean(uv.pk/peak.tri, peak.tri/sum(peak.tri))*10^6)
  ros.pk <- res; ros.pk[names(table(ros.muts$TriN[queryHits(atac.ros)]))] <- as.numeric(table(ros.muts$TriN[queryHits(atac.ros)]))
  ros.pk <- (weighted.mean(ros.pk/peak.tri, peak.tri/sum(peak.tri))*10^6)
  
  return(c('UV_genome'=uv.bg, 'UV_ATAC'=uv.pk, 'ROS_genome'=ros.bg, 'ROS_ATAC'=ros.pk))
  
}))

# Data frame of means and standard deviation
peak.deltas <- data.frame('damage'=c('UV','ROS'), 'mean'=c(mean(peak.rates[,'UV_ATAC']/peak.rates[,'UV_genome']), mean(peak.rates[,'ROS_ATAC']/peak.rates[,'ROS_genome']) ),
                          'sd'=c(sd(peak.rates[,'UV_ATAC']/peak.rates[,'UV_genome']), sd(peak.rates[,'ROS_ATAC']/peak.rates[,'ROS_genome']) ))

# Data frame of individual points
pk.pt <- data.frame('source'=rep(c('ROS', 'UV'), each=7), 
                    'rate'=c(peak.rates[,'ROS_ATAC']/peak.rates[,'ROS_genome'], peak.rates[,'UV_ATAC']/peak.rates[,'UV_genome']))

# Plot Panel 4f
p <-ggplot(peak.deltas, aes(x=damage, y=mean, fill=damage))

p + geom_bar(stat="identity", color=NA, position=position_dodge()) +
  geom_point(data=pk.pt, aes(x=source,y=rate,fill=source), position = position_jitterdodge(jitter.width=0.4),size=1.7,
             color=grey.colors(10)[3] ) +
  geom_errorbar(aes(ymin=mean - 2*sd, ymax=mean + 2*sd), width=.2, position=position_dodge(.9)) +
  labs(title="", x="damage", y = "rate/background") + theme_classic() + scale_fill_manual(values=c('turquoise3','sienna3')) +
  scale_y_continuous(name ="peak rate/background rate\n(mu/mb)", limits=c(0,1)) +
  theme(axis.text.x = element_text(face="bold", color="black", size=14), axis.text.y = element_text(face="bold", color="black", size=14),
        axis.title=element_text(size=14,face="bold"))



#######  Panel 4g  ########
#######------------########

# Create random SNP positions, with same number of SNPs on each chromosome as control
set.seed(101)
rnd.positions <- randomizeRegions(all.unique, genome=BSgenome.Mmusculus.UCSC.mm10, allow.overlaps = FALSE, per.chromosome = TRUE)
seqlevels(rnd.positions) <- seqlevelsInUse(rnd.positions)
names(rnd.positions) <- 1:length(rnd.positions)
              
# QuasR object
proj <- qAlign(sampleFile = '/omics/groups/OE0538/internal/users/p281o/projects/lce_mechanism/sequencing/seqID_30724/samples_quasr.txt',
  genome = 'BSgenome.Mmusculus.UCSC.mm10', 
  paired = 'fr',
  cacheDir = '/omics/groups/OE0538/internal/users/p281o/tmp')
              
# Profile read depth around SNPs
window.size <- 10001
              
# Table with 10kb profile around all mutations
snp <- qProfile(proj, query=all.unique, upstream=window.size, clObj = cluObj)
atac.snps <- snp$replicate1
              
# Table with 10kb profile around randomized positions of overlapping genes for each expression category
rnd <- qProfile(proj, query=rnd.positions, upstream=window.size, clObj = cluObj)
rnd.snps <- rnd$replicate1

# Distinguish UV and ROS mutations
uv.muts <- ((all.unique$REF=='C' & all.unique$ALT=='T') | (all.unique$REF=='G' & all.unique$ALT=='A')) & all.unique$TimePoint=='T1'
ros.muts <- ((all.unique$REF=='G' & all.unique$ALT=='T') | (all.unique$REF=='C' & all.unique$ALT=='A'))

# Smooth signal
smth.win <- 201
uv.signal <- as.vector(runmean(Rle(log2(colMeans(atac.snps[uv.muts,])+1)), k=smth.win, endrule = 'constant'))
ros.signal <- as.vector(runmean(Rle(log2(colMeans(atac.snps[ros.muts,])+1)), k=smth.win, endrule = 'constant'))
rnd.signal <- as.vector(runmean(Rle(log2(colMeans(rnd.snps[uv.muts,])+1)), k=smth.win, endrule = 'constant'))

# Plot Panel 4g
line.width <- 3
plot(-window.size:window.size, uv.signal, type='l', ylim=c(0.1,0.2), xlab='pos rel. to SNP', ylab='mean signal', lwd=line.width, 
     col='firebrick3', cex.axis=1.5, cex.lab=2)
lines(-window.size:window.size, ros.signal, type='l', col='dodgerblue4', lwd=line.width)
lines(-window.size:window.size, rnd.signal, type='l', col='darkgrey', lwd=line.width)
legend('bottomright', legend=c('UV', 'ROS', 'random'), fill=c('firebrick3', 'dodgerblue4', 'darkgrey'), bty='n', border=NA, cex=2)
abline(v=0, col='darkgrey', lty=3, lwd=2)

