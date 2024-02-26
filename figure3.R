############################
#####     Figure 3     #####  
############################
require(GenomicRanges)
require(Rsamtools)
require(scales)

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

# Create GRanges of all unique mutations (eliminate shared T0 mutations from one sister)
all.unique <- unlist(as(lapply(snp.match, function(x){
  sis1 <- sister.split[[x]]; sis1 <- sis1[sis1$TimePoint=='T1']
  sis2 <- sister.split[[(x+1)]]
  unlist(as(list(sis1, sis2), 'GRangesList'))
}),'GRangesList'))

# Split mutations back into single samples
rds.sis <- split(all.unique, all.unique$SampleName)

## Multiallelism at CC dinucleotide mutations
##-------------------------------------------
sis.id <- data.frame(do.call(rbind, lapply(rds.sis, function(x){
  
  # Select sister unique mutations adjacent to eachother
  rel.snps <- x[x$TimePoint=='T1']
  rel.dist <- distanceToNearest(rel.snps)
  dist.mat <- rel.dist[values(rel.dist)$distance == 0]
  
  # Remove duplicated entries
  f.mat <- data.frame("V1"=queryHits(dist.mat), "V2"=subjectHits(dist.mat), 'dist_indx'=1:length(dist.mat))
  f.mat[1:2] <- t( apply(f.mat[1:2], 1, sort) )
  dist.mat <- dist.mat[(1:length(dist.mat) %in% (f.mat[!duplicated(f.mat[1:2]),'dist_indx']))]
  
  # Remove trinucleotide substitutions
  all.positions <- c(queryHits(dist.mat), subjectHits(dist.mat))
  dist.mat <- dist.mat[ !((queryHits(dist.mat) %in% all.positions[duplicated(all.positions)]) | (subjectHits(dist.mat) %in% all.positions[duplicated(all.positions)]))]
  
  # Return data frame
  snp.one <- paste0(rel.snps$REF[queryHits(dist.mat)], '>', rel.snps$ALT[queryHits(dist.mat)])
  snp.two <- paste0(rel.snps$REF[subjectHits(dist.mat)], '>', rel.snps$ALT[subjectHits(dist.mat)])
  dual.mut <- paste0(rel.snps$REF[queryHits(dist.mat)], rel.snps$REF[subjectHits(dist.mat)], '>', rel.snps$ALT[queryHits(dist.mat)], rel.snps$ALT[subjectHits(dist.mat)])
  mut.pairs <- data.frame(cbind('Mutation_One'=snp.one, 'Mutation_Two'=snp.two), "Dual_Mutation"=dual.mut, 
                          'VAF1'= rel.snps$VAF[queryHits(dist.mat)], 'VAF2' = rel.snps$VAF[subjectHits(dist.mat)],
                          'chrom'=as.character(seqnames(rel.snps[queryHits(dist.mat)])), 'start'=start(rel.snps[queryHits(dist.mat)]), 'end'=end(rel.snps[subjectHits(dist.mat)]))
  mut.pairs$distance <- values(dist.mat)$distance

  return(mut.pairs)

})), stringsAsFactors=FALSE)

# Subset on only UV dual mutations
uv.dual <- sis.id[sis.id$Dual_Mutation=='CC>TT'| sis.id$Dual_Mutation=='GG>AA',]
uv.dual$SisterName <- gsub('\\.\\d+', '', rownames(uv.dual))
uv.gr <- makeGRangesFromDataFrame(uv.dual, seqnames.field = 'chrom', start.field = 'start', end.field = 'end', keep.extra.columns = TRUE)


##  Extracting Reads over dual mutations  
##----------------------------------------
# Chip 53
ref.key <- read.delim('/omics/groups/OE0538/internal/users/p281o/publications/single_cell_split/sra_transfer/Sample_key.txt', sep='\t', stringsAsFactors = FALSE)
ref.slice <- ref.key$SRA_names; names(ref.slice) <- ref.key$Name
df.53 <- read.delim('/omics/groups/OE0538/internal/users/p281o/publications/single_cell_split/chip53_NewCall/snv_input.txt', sep='\t', stringsAsFactors = FALSE, header=FALSE)
bams.53 <- df.53[,2]; names(bams.53) <- as.character(c(ref.slice[gsub('_mf_unq.bam','', basename(df.53[1:8,2]))], ref.slice[18:19]))

# Chip 43
df.43 <- read.delim('/omics/groups/OE0538/internal/users/p281o/projects/lce_mechanism/sequencing/seqID_24871/snv_input.txt', sep='\t', stringsAsFactors = FALSE, header=FALSE)
bams.43 <- df.43[-7,1]
slice.43 <- ref.key$SRA_names; names(slice.43) <- ref.key$ASID
names(bams.43) <- slice.43[gsub('AS-(\\d+)_sort.*','\\1', basename(bams.43))]

# All bams
bam.files <- c(bams.53, bams.43)

# Relevant info to retrieve for overlapping reads
what <- c("strand", "pos", "qwidth", "seq", "cigar")
rev.comp.vec <- c("A"="T","C"="G","G"="C","T"="A")

# Read evidence for each dual mutation
all.snvs <- lapply(1:length(uv.gr), function(z){
  
  # Current dual mutation
  which <- uv.gr[z]
  param <- ScanBamParam(which=which, what=what, reverseComplement = TRUE, flag=scanBamFlag(isPaired=TRUE, isProperPair = TRUE), mapqFilter = 20)
  bamFile <- bam.files[[uv.gr$SisterName[z]]]
  bam <- scanBam(bamFile, param=param)
  
  # Get insertion and deletion info
  cigar <- bam[[1]]$cigar
  has.insert <- grepl('I', cigar)
  has.deletion <- grepl('D', cigar)
  
  # Test that read covers both bases
  both.bases <- which((bam[[1]]$pos <= start(which)) & (bam[[1]]$pos + bam[[1]]$qwidth >= (end(which)+1)))
  
  # Iterate through each set of reads and elminate entries with an insertion or deletion
  do.call(rbind, lapply(both.bases, function(x){
    matches <- unlist(strsplit(cigar[x],'\\D+'))
    ids <- unlist(strsplit(cigar[x],'\\d+'))[-1]
    base.positions <- (start(uv.gr[z])-bam[[1]]$pos[x] + 1):(end(uv.gr[z])-bam[[1]]$pos[x] + 1)
    
    # Compensate for insertions in reads
    if(has.insert[x]){
      
      insertion.pos <- unlist(lapply(which(ids %in% 'I'), function(i){ sum(as.numeric(matches[1:(i-1)]))+1}))
      insertion.size <- as.numeric(matches[(which(ids%in%'I'))])
      insertion.adj <- 0
      
      for(ii in 1:length(insertion.pos)){
        # Cases where insertion affects read position
        if((base.positions[1] + insertion.size[ii] > bam[[1]]$qwidth[x]) ){
          base.positions <- NA
        } 
        if( (base.positions[1] > insertion.pos[ii]) & (base.positions[2]+insertion.size[ii] < bam[[1]]$qwidth[x]) ){
          insertion.adj <- insertion.adj + insertion.size[ii]
        } 
        if( (base.positions[2] < insertion.pos[ii]) ){
          insertion.adj <- insertion.adj
        } 
        base.positions <- base.positions + insertion.adj
      }
    }
    
    # Compensate for deletions in reads
    if(has.deletion[x]){
      
      deletion.pos <- unlist(lapply(which(ids %in% 'D'), function(d){ sum(as.numeric(matches[1:(d-1)]))+1}))
      deletion.size <- as.numeric(matches[(which(ids%in%'D'))])
      deletion.adj <- 0
      
      for(dd in 1:length(deletion.pos)){
        
        if((base.positions[2] < deletion.pos[dd])){
          deletion.adj <- 0 
        } 
        
        if((base.positions[1] > deletion.pos[dd]) & (base.positions[1] - deletion.size[dd] > 0) ){
          deletion.adj <- deletion.adj - deletion.size[dd]
        } 
        # Cases where deletion affects mutation position
        if( (base.positions[1] >= deletion.pos[dd]) & (deletion.pos[dd] + deletion.size[dd] > base.positions[1]) ){
          base.positions <- NA
        } 
        base.positions <- base.positions + deletion.adj
      }
    }
    
    if(bam[[1]]$strand[[x]]=='+'){
      current.read <- as.character(bam[[1]]$seq[[x]])
    } else{
      current.read <- as.character(reverseComplement(bam[[1]]$seq[[x]]))
    }
    
    # Return sequences from modified base position
    if( sum(is.na(base.positions))>0 ){ res <- 'XX'} else{ res <- paste(unlist(strsplit(current.read,''))[base.positions], collapse='') }
    res
  }))
})

# Select multiallelic sites under the criteria: At least 2 alternate alleles with >= 3 reads as evidence
is.mult <- (unlist(lapply(all.snvs, function(x){
  res <- table(x)
  res <- res[as.numeric(res)>= 3]
  length(res) > 2
})))
mu.id <- all.snvs[(is.mult)]



#######  Panel 3a  ########
#######------------########

# Bases to flank
window.flank <- 10

# Example sites of biallelic and multiallelic mutations
dual.alleles <- c('Biallelic'=64, 'Multiallelic'=66)

# Plot biallelic and multiallelic examples
for(dp.ex in 1:length(dual.alleles)){
  message(dp.ex)
  
  # Rsamtools parameters
  which <- uv.gr[dual.alleles[[dp.ex]]]
  param <- ScanBamParam(which=which, what=what, reverseComplement = TRUE, flag=scanBamFlag(isPaired=TRUE, isProperPair = TRUE), mapqFilter = 20)
  bamFile <- bam.files[[uv.gr$SisterName[dual.alleles[[dp.ex]] ]]]
  bam <- scanBam(bamFile, param=param)
  
  # Get insertion and deletion info
  cigar <- bam[[1]]$cigar
  has.insert <- grepl('I', cigar)
  has.deletion <- grepl('D', cigar)
  
  # Base vector
  b.vec <- c('N'=0, 'A'=2, 'C'=3,'T'=4,'G'=5)
  
  # Get read events
  read.mat <- do.call(rbind, lapply(1:length(bam[[1]]$seq) , function(x){
    #message(x)
    
    # Result vector
    res <- rep(0, ((window.flank * 2) + 2) )
    
    # Adjust for insertions or deletions
    matches <- unlist(strsplit(cigar[x],'\\D+'))
    ids <- unlist(strsplit(cigar[x],'\\d+'))[-1]
    base.positions <- (start(uv.gr[dual.alleles[[dp.ex]]])-bam[[1]]$pos[x] + 1):(end(uv.gr[dual.alleles[[dp.ex]]])-bam[[1]]$pos[x] + 1)
    width.adjust <- 0
    read.length <- bam[[1]]$qwidth[x]
    
    # Check if insertion
    if(has.insert[x]){
      
      insertion.pos <- unlist(lapply(which(ids %in% 'I'), function(i){ sum(as.numeric(matches[1:(i-1)]))+1}))
      insertion.size <- as.numeric(matches[(which(ids%in%'I'))])
      insertion.adj <- 0
      
      for(ii in 1:length(insertion.pos)){
        # Cases where insertion shifts outside of read info
        if((base.positions[1] + insertion.size[ii] > bam[[1]]$qwidth[x]) ){
          base.positions <- c(NA,NA)
        } 
        if( (base.positions[1] > insertion.pos[ii]) & (base.positions[2]+insertion.size[ii] < bam[[1]]$qwidth[x]) ){
          insertion.adj <- insertion.adj + insertion.size[ii]
        } 
        if( (base.positions[2] < insertion.pos[ii]) ){
          insertion.adj <- insertion.adj
        } 
        base.positions <- base.positions + insertion.adj
      }
    }
    
    # Check if deletion
    if(has.deletion[x]){
      
      deletion.pos <- unlist(lapply(which(ids %in% 'D'), function(d){ sum(as.numeric(matches[1:(d-1)]))+1}))
      deletion.size <- as.numeric(matches[(which(ids%in%'D'))])
      deletion.adj <- 0
      
      for(dd in 1:length(deletion.pos)){
        
        if((base.positions[2] < deletion.pos[dd])){
          deletion.adj <- 0 
        } 
        
        if((base.positions[1] > deletion.pos[dd]) & (base.positions[1] - deletion.size[dd] > 0) ){
          deletion.adj <- deletion.adj - deletion.size[dd]
        } 
        # Cases where deletion shifts below read coverage
        if( (base.positions[1] >= deletion.pos[dd]) & (deletion.pos[dd] + deletion.size[dd] > base.positions[1]) ){
          base.positions <- c(NA,NA)
        } 
        base.positions <- base.positions + deletion.adj
      }
    }
    
    if(bam[[1]]$strand[[x]]=='+'){
      current.read <- as.character(bam[[1]]$seq[[x]])
    } else{
      current.read <- as.character(reverseComplement(bam[[1]]$seq[[x]]))
    }
    
    
    if( sum(is.na(base.positions))>0 ){ 
      res <- res
    } else{ 
      base.changes <- unlist(strsplit(current.read,''))[(base.positions[1]):(base.positions[2])] 
      
      # If read extends beyon 3 prime
      if( (read.length + width.adjust) - base.positions[2] > 0){
        
        last.base <- min(c(2*window.flank + 2, ((window.flank +2) + (read.length + width.adjust) - base.positions[2])))
        res[((window.flank) + 3):last.base] <- 1
      } 
      # If second base in mutation is last in read
      if((read.length + width.adjust) == base.positions[2]){
        res[(window.flank +3):(2*window.flank + 2)] <- 0
      }
      
      # If read extends before 5 prime
      if( base.positions[1]==1){
        
        res[1:window.flank] <- 0
      } else{
        bases.covered <- min(c(10, base.positions[1]))
        first.base <- 11 - bases.covered
        res[first.base:window.flank] <- 1
      }
      # make base change color at mutation site
      res[c(window.flank+c(1,2))] <- as.numeric(b.vec[base.changes])
    }
    res
  }))
  
  # Order by read start
  row.order <- order(unlist(lapply(1:nrow(read.mat), function(x){ min(which(read.mat[x,]>=1))})), decreasing=FALSE)
  read.mat <- read.mat[row.order,]
  
  # Plot
  pheatmap(read.mat, cluster_rows = FALSE, cluster_cols = FALSE, color = c('white', 'grey', 'green', 'blue','red','orange'), 
           border_color = NA, legend=FALSE, main = names(dual.alleles)[dp.ex], filename = paste0(names(dual.alleles)[dp.ex],'_Figure3a.pdf'))
}



#######  Panel 3b  ########
#######------------########

# Assign VAF of each base relative to replication fork direction
first.fork <- rep(NA, nrow(uv.dual)); first.fork[uv.dual$Dual=='CC>TT'] <- uv.dual$VAF2[uv.dual$Dual=='CC>TT']; first.fork[uv.dual$Dual=='GG>AA'] <- uv.dual$VAF1[uv.dual$Dual=='GG>AA']
second.fork <- rep(NA, nrow(uv.dual)); second.fork[uv.dual$Dual=='CC>TT'] <- uv.dual$VAF1[uv.dual$Dual=='CC>TT']; second.fork[uv.dual$Dual=='GG>AA'] <- uv.dual$VAF2[uv.dual$Dual=='GG>AA']

# Plot Panel 3b
plot(first.fork, second.fork, col=alpha(ifelse(is.mult,'red','black'),alpha = 0.7), pch=19, xlab='cC VAF', ylab='Cc VAF', cex=1.3)
legend('topleft', legend=c('biallelic', 'multiallelic'), pch=19, col=c(1,2), bty='n')



#######  Panel 3c  ########
#######------------########

# Identities of dual mutants
mut.cat <- table(sis.id$Dual)
dual.pyrimidines <- c('TC','CC','CT','TT','GA','GG','AG','AA')
bar.col <- rep('darkgrey', length(mut.cat))
bar.col <- rep('darkgrey', length(mut.cat)); bar.col[unlist(lapply(names(mut.cat), function(x){unlist(strsplit(x, '>'))[[1]]})) %in% dual.pyrimidines] <- 'firebrick4'
barplot(mut.cat, col=bar.col, las=2)

# What are read representations across duals with different VAFs
most.common.vars <- do.call(rbind, lapply(mu.id, function(x){
  read.id <- table(x)/sum(as.numeric(table(x)))
  read.id <- read.id[order(as.numeric(read.id), decreasing=TRUE)]
  read.id <- read.id[!(names(read.id) %in% c('GG','CC'))]
  c(as.numeric(read.id)[1:2], names(read.id)[1:2])
}))

# Plot Panel 3c
boxplot(list(as.numeric(most.common.vars[,c(1)])* 100,as.numeric(most.common.vars[,c(2)])* 100), 
        ylim=c(0,50), pch=19, notch=TRUE, cex=0.4, names=c('1st', '2nd'), ylab='VAF')



#######  Panel 3d  ########
#######------------########

# Vectors of most common first and second alternate allleles
first.alternate <- round(table(most.common.vars[,3])/sum(table(most.common.vars[,3])), digits=2)
first.alternate <- first.alternate[order(first.alternate, decreasing=TRUE)]
second.alternate <- round(table(most.common.vars[,4])/sum(table(most.common.vars[,4])), digits=2)
second.alternate <- second.alternate[order(second.alternate, decreasing=TRUE)]

# Modify to reflect C to T Change
first.alternate <- c('CC>TT'=sum(first.alternate[c('AA','TT')]), 'CC>CT'=sum(first.alternate[c('CT','AG')]), 'CC>TC'=sum(first.alternate[c('TC','GA')]),'other'=1-sum(first.alternate[c('AA','TT','CT','AG','TC','GA' )]))
second.alternate <- c('CC>TT'=sum(second.alternate[c('AA','TT')]), 'CC>CT'=sum(second.alternate[c('CT','AG')]), 'CC>TC'=sum(second.alternate[c('TC','GA')]),'other'=1-sum(second.alternate[c('AA','TT','CT','AG','TC','GA' )]))

# Color vector for mutation identitiy
all.subs <- unique( c(names(first.alternate), names(second.alternate)) )
bar.colors <- rainbow(length(all.subs)); names(bar.colors) <- all.subs

# Vector for dual mutation identity, variant allele 1 followed by variant allele 2
dual_mut <- table(paste0(most.common.vars[,3],'>', most.common.vars[,4]))
dual_mut2 <- c('CT>TT'=sum(dual_mut[c('CT>TT', 'AG>AA')]), 'TT>CT'=sum(dual_mut[c('TT>CT','AA>AG')]), 'TC>TT'=sum(dual_mut[c('TC>TT', 'GA>AA')]), 
               'TT>TC'=sum(dual_mut[c('TT>TC', 'AA>GA')]), 'other'=sum(dual_mut)-sum(dual_mut[c('CT>TT', 'AG>AA','TT>CT','AA>AG','TC>TT', 'GA>AA','TT>TC', 'AA>GA')]))


# Plot Panel 3d
barplot(rev(dual_mut2), horiz = TRUE, xlab='# of multiallelic sites')
