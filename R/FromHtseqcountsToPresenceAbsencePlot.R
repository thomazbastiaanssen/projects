###IMPORTANT###
#Always run line 1:92 before using this script, as it contains custom functions. 

###set up and load
options(stringsAsFactors = F)

### load heatmap libraries
library("RColorBrewer");
library("gplots");
library("Biobase");
library("biomaRt");
library("DESeq2");

### load custom functions and wrapper
#DESeqprep takes a data.frame structure in the format of the htseqcount output containing only the two conditions that are to be compared.

DESeqPrep <- function(x){
  condition  <- factor(gsub("[0-9]*$", "", colnames(x)))
  colData    <- data.frame(condition , row.names = colnames(x))
  levels     <- as.vector(unique(condition))
  return(DESeqDataSetFromMatrix(countData = x, 
                                colData   = colData,
                                design    = formula(~ condition))
  )
}

### 
#extractResults gathers the results from the statistical analysis performed by DESeq2 into a dataframe.

extractResults <- function(relevantComparisons){
  outputdf  <-data.frame(matrix(NA, nrow = length(row.names(relevantComparisons[[1]])), ncol = 2*length(relevantComparisons)))
  colnumber <- 1
  for(comparison in relevantComparisons){
    prepdds      = DESeqPrep(comparison)                
    tempdds      = DESeq(prepdds)
    tempres      = results(tempdds)
    outputdf[,colnumber]                                          = tempres$padj
    colnames(outputdf)[colnumber]                                 = paste("padj", paste0(unique(gsub("[0-9]*$", "", colnames(comparison))), collapse=""))
    outputdf[,(colnumber + length(relevantComparisons))]          = tempres$log2FoldChange
    colnames(outputdf)[(colnumber + length(relevantComparisons))] = paste0(unique(gsub("[0-9]*$", "", colnames(comparison))), collapse="")
    colnumber = colnumber +1
  }
  rownames(outputdf) <- rownames(relevantComparisons[[1]])
  return(outputdf)
}

### 
#formatForPlot formats the output extractResults to a presence absence matrix ready for plotting

formatForPlot <- function(resultdf = resultdf, cutoff = 0.1){
  uplist <- list()
  dolist <- list()
  comparison = 1
  
  for(comparison in 1:(ncol(resultdf)/2)){
    uplist[[comparison]] <- row.names(resultdf)[intersect(which(resultdf[(ncol(resultdf)/2)+comparison] >0), which(resultdf[comparison] < cutoff))]
    names(uplist)[comparison] <- paste("UP",names(resultdf[comparison +9]), sep=".")
    dolist[[comparison]] <- row.names(resultdf)[intersect(which(resultdf[(ncol(resultdf)/2)+comparison] <0), which(resultdf[comparison] < cutoff))]
    names(dolist)[comparison] <- paste("DO",names(resultdf[comparison +9]), sep=".")
  }
  totlist <- append(uplist, dolist)
  
  totmat <- sapply(totlist, function(x){as.numeric(sort(unique(unlist(totlist))) %in% x)}) 
  rownames(totmat) <- sort(unique(unlist(totlist))) 
  
  upmat <- mymat[, 1:length(uplist) ]
  domat <- mymat[,(1+length(uplist)):(2*length(uplist))]*-1
  totmat <- as.matrix(data.frame(upmat, domat))
  outputmat <- (totmat[,1:9]+ totmat[,10:18])
  colnames(outputmat) <- names(resultdf[1:(ncol(resultdf)/2) +9])
  return(outputmat)
}

###
#plotPA is a wrapper to plot presence/absence plots. The filter can be used to exclude columns from the graph. 

plotPA <-function(PAmatrix = PAmatrix, filter = c(1, 2, 3, 4, 5, 7, 6, 8, 9), title = "", col = redgreen(75) ){
  return(
    heatmap.2(
      PAmatrix[,filter][which(rowSums(PAmatrix != 0)>1),], 
      margin     = c(8, 4),
      col        = col, 
      trace      = "none", 
      dendrogram = "none", 
      Colv       = F, 
      Rowv       = T, 
      main       = title 
    )
  )
}

##################################################
##################################################
### You may start modifying the script from here.
### load the htseqcount output file, make sure the path is correct.
setwd("C:/Users/bioscgen3/Desktop/Thomaz/R/Raw")
allcounts <- read.delim("htseqcount.txt",  header=T, row.names=1)


### define the separate conditions from the htseqcount output file 
conp <- allcounts[,1 :4 ]
cons <- allcounts[,5 :12]
gfp  <- allcounts[,13:16]
gfrp <- allcounts[,17:20]
gfrs <- allcounts[,21:28]
gfs  <- allcounts[,29:40]

###define the comparisons
#p vs s:
countsconpcons <- data.frame(conp , cons)
countsgfpgfs   <- data.frame(gfp  , gfs )
countsgfrpgfrs <- data.frame(gfrp , gfrs)

#con vs gf:
countsconpgfp  <- data.frame(conp , gfp )
countsconsgfs  <- data.frame(cons , gfs )

#con or gf vs gfr:
countsconsgfrs  <- data.frame(cons, gfrs)
countsconpgfrp  <- data.frame(conp, gfrp)
countsgfpgfrp   <- data.frame(gfp , gfrp)
countsgfsgfrs   <- data.frame(gfs , gfrs)


###Create a vector with the names of all relevant comparisons as defined above.
relevantComparisons <- list(countsconpcons, countsgfpgfs, countsgfrpgfrs, countsconpgfp, countsconsgfs, countsconsgfrs, countsconpgfrp, countsgfpgfrp, countsgfsgfrs)

###This will take a couple of minutes
resultdf <- extractResults(relevantComparisons = relevantComparisons)

###Generate a matrix to be plotted, adjust the desired adjusted p-value using the 'cutoff' argument.
PAmatrix <- formatForPlot(resultdf = resultdf, cutoff = 0.1)

###Plot the presence/absence plot!
plotPA(PAmatrix, title ="untitled" )


