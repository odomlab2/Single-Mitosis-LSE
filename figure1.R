############################
#####     Figure 1     #####  
############################
require(apcluster)
require(pheatmap)


#######  Panel 1c  ########
#######------------########

# Set working directory where you would like the files to be downloaded, or alternatively update file paths below
# Download GitHub respository, extract
download.file(url = "https://github.com/odomlab2/Single-Mitosis-LSE/archive/refs/heads/main.zip",
              destfile = 'Single_Mitosis_LSE_GitHub.zip')
unzip(zipfile = 'Single_Mitosis_LSE_GitHub.zip', exdir = './', )
file.remove('Single_Mitosis_LSE_GitHub.zip')

# Path to tables Phenomex cell measurements
bl.path <- 'Single-Mitosis-LSE-main/Phenomex_Tables'

# Load time points 
tps <- lapply(list.files(bl.path,pattern='TP', full.names = TRUE), read.delim, sep='\t', stringsAsFactors= FALSE)
names(tps) <- unlist(lapply(list.files(bl.path,pattern='TP'), function(x){unlist(strsplit(x,'_'))[[1]]}))

# Select PenIDs seen in all timepoints
pen.ids <- unlist(lapply(tps, function(x) { unique(x$PenId) }))
all.timepoints <- as.numeric(names(table(pen.ids))[as.numeric(table(pen.ids))==length(tps)])

# Remove PenIDs where 
# 1) there is more than one cell in timepoint 0 
# 2) Cells are very small or very large
sc.size <- unique(tps[[1]]$PenId[tps[[1]]$CellsInPen == 1 & tps[[1]]$DiameterMicrons >= 6 & tps[[1]]$DiameterMicrons <= 10])
all.timepoints <- all.timepoints[(all.timepoints %in% sc.size)]

# Remove columns with no added information
col.prune <- colnames(tps[[1]])[unlist(lapply(1:ncol(tps[[1]]), function(x) length(table(tps[[1]][,x])) > 1 ))]
tpf <- lapply(names(tps), function(x){ 
  tmp <- tps[[x]][(tps[[x]]$PenId %in% all.timepoints), col.prune]
  tmp <- tmp[order(tmp$View, tmp$PenId),]
  tmp })

# Cell counts per Timepoint
cell.counts <- do.call(cbind, lapply(tpf, function(x){ lengths(split(x$CellsInPen, x$PenId))/4}))

# Filter pens where there was a doubling and later decrease in cell number (n = 5)
prolif.filt <- unlist(lapply(1:nrow(cell.counts), function(x){ identical(order(cell.counts[x,]), 1:ncol(cell.counts))}))
cell.counts <- cell.counts[prolif.filt,]

# Filter Pens where there was a jump from 1 to 3 cells (n = 1, counting error)
cell.jump <- unlist(lapply(1:nrow(cell.counts), function(x){ cts.tmp <- cell.counts[x,] 
max(unlist(lapply(2:length(cts.tmp), function(z){ cts.tmp[z] - cts.tmp[(z-1)] }))) > 1 }))
cell.counts <- cell.counts[!(cell.jump),]

# Normalize all fluorophore signals to background, scale from 0 - 1000
tp.sc <- lapply(tpf, function(x){
  filt.x <- x[x$PenId %in% rownames(cell.counts), ]
  fluoro <- split(filt.x, filt.x$Cube)
  fl.cts <- do.call(cbind, lapply(fluoro, function(z){
    sc.signal <- z$MedianBrightness - z$MedianBackgroundBrightness
    sc.signal[sc.signal < 0] <- 0
    scales::rescale(sc.signal, to=c(0,1000)) }))
  res <- data.frame('View'=fluoro[[1]]$View, 'PenId'=fluoro[[1]]$PenId, 'TargetIndex'=fluoro[[1]]$TargetIndex, 
                    'CellsInPen'=fluoro[[1]]$CellsInPen, fl.cts[,-grep('OEP',colnames(fl.cts))], stringsAsFactors=FALSE)
  res
}); names(tp.sc) <- names(tps)

# Take mean delta of divided cells
delta.mat <- matrix(rep(NA, length(cell.counts)), ncol=ncol(cell.counts), nrow=nrow(cell.counts))
rownames(delta.mat) <- rownames(cell.counts); colnames(delta.mat) <- names(tp.sc)
for(i in names(tp.sc)){
  delta.sig <- log2(tp.sc[[i]]$TRED + 1)-log2(tp.sc[[i]]$FITC + 1)
  names(delta.sig) <- tp.sc[[i]]$PenId
  by.pen.id <- lapply(split(delta.sig, names(delta.sig)), mean)
  delta.mat[names(by.pen.id), i] <- as.numeric(by.pen.id)
}

# Assign clusters to cells based on FUCCI signal across the time course
g1.mat <- delta.mat[rowSums(is.na(delta.mat))==0,]
is.cycling <- (as.numeric(apply(g1.mat, 1, var)) >= 4)
g1.mat <- g1.mat[is.cycling, ]
f.cl <- apcluster(negDistMat(r=2), g1.mat, q=0.01)
cell.clusters <- f.cl@clusters
g1.mat <- g1.mat[as.numeric(unlist(cell.clusters)), ]

# Data frame of cluster signal
cluster.id <- data.frame('cluster'=rep(1:length(cell.clusters), lengths(cell.clusters)), stringsAsFactors = FALSE)
rownames(cluster.id) <- rownames(g1.mat)

# Subselect on G1, S and G2/M clusters
gaps.sel <- lengths(cell.clusters)[c(5,10)]; gaps.sel[2] <- sum(gaps.sel)
cc.sel <- c(names(cell.clusters[[5]]), names(cell.clusters[[10]]), names(cell.clusters[[4]]))
mat.sel <- g1.mat[cc.sel,]

# Plot Panel 1c
pheatmap(mat.sel, cluster_rows=FALSE, cluster_cols = FALSE, show_rownames = FALSE, na_col = 'grey', cellwidth = 20, 
         color = colorRampPalette(c('forestgreen','khaki', 'firebrick2'))(100), border=NA, gaps_row = gaps.sel, gaps_col=1, 
         cellheight=2, border_color=NA, legend_breaks = c(-5, 0, 5), legend_labels = c("G2/M", "S", "G1"))


