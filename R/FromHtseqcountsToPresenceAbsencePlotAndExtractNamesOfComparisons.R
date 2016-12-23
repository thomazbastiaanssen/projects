###IMPORTANT###
#Always run line 1:150 before using this script, as it contains custom functions. 

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
  
  upmat <- totmat[, 1:length(uplist) ]
  domat <- totmat[,(1+length(uplist)):(2*length(uplist))]*-1
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
      PAmatrix[,filter][which(abs(rowSums(PAmatrix[,filter]))!=0),], 
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

###getNamesOfComparisons is a function that produces dataframes of all genes that are relevant between a comparison on the PA gene. 
#The function takes either one or a vector with two integers representing the column number you want to select. 
#Its function depends on the mode. 
#"similar" gives you the genes that are upregulated and downregulated in both comparisons in the vector.
#"different" gives you the genes that are upregulated or downregulated in in one but not the other vector.
#"simple" gives you the genes that are upregulated and downregulated in the column specified.
getNamesOfComparisons <- function(mode = "similar", x){
  if     (mode == "similar"){
    
    outlist <- list("upregulated in both"  = names(which( PAmatrix[,x[1]] == PAmatrix[,x[2]] & PAmatrix[,x[1]]== 1  )),
                    "downregulated in both"= names(which( PAmatrix[,x[1]] == PAmatrix[,x[2]] & PAmatrix[,x[1]]==-1 )))
    
    maxlen = max(length(outlist[[1]]),
                 length(outlist[[2]]))
    
    outdf = data.frame("upregulated in both"   = c(outlist[[1]], rep(NA, maxlen - length(outlist[[1]]))), 
                       "downregulated in both" = c(outlist[[2]], rep(NA, maxlen - length(outlist[[2]]))))
    
    
    return(outdf)
  }
  else if(mode == "different"){
    outlist  <- list("upregulated in nr1"   = names(which( PAmatrix[,x[1]] != PAmatrix[,x[2]] & PAmatrix[,x[1]]== 1  )), 
                     "downregulated in nr1" = names(which( PAmatrix[,x[1]] != PAmatrix[,x[2]] & PAmatrix[,x[1]]==-1  )), 
                     "upregulated in nr2"   = names(which( PAmatrix[,x[1]] != PAmatrix[,x[2]] & PAmatrix[,x[2]]== 1  )), 
                     "downregulated in nr2" = names(which( PAmatrix[,x[1]] != PAmatrix[,x[2]] & PAmatrix[,x[2]]==-1 ))) 
    
    maxlen = max(length(outlist[[1]]),
                 length(outlist[[2]]),
                 length(outlist[[3]]),
                 length(outlist[[4]]))
    
    outdf = data.frame("upregulated in nr1"   = c(outlist[[1]], rep(NA, maxlen - length(outlist[[1]]))), 
                       "downregulated in nr1" = c(outlist[[2]], rep(NA, maxlen - length(outlist[[2]]))),
                       "upregulated in nr2"   = c(outlist[[3]], rep(NA, maxlen - length(outlist[[3]]))),
                       "downregulated in nr2" = c(outlist[[4]], rep(NA, maxlen - length(outlist[[4]])))
                              )
    
    return(outdf)
  }
  else if(mode == "simple"){
    outlist      <- list("upregulated"      = names(which( PAmatrix[,x[1]] == 1  )),
                         "downregulated"    = names(which( PAmatrix[,x[1]] ==-1 )))
    
    maxlen = max(length(outlist[[1]]),
                 length(outlist[[2]]))
    
    outdf = data.frame("upregulated"   = c(outlist[[1]], rep(NA, maxlen - length(outlist[[1]]))), 
                       "downregulated" = c(outlist[[2]], rep(NA, maxlen - length(outlist[[2]]))))
    return(outdf)
  }
  else{
    print("Error, please select a mode. Options are: similar, different or simple (between quotes)")
  }
}

##################################################
##################################################
### You may start modifying the script from here.
### load the htseqcount output file, make sure the path is correct.
#setwd("C:/Users/bioscgen3/Desktop/Thomaz/R/Raw")
setwd("C:/Users/Thomaz/Desktop/School/Master/R/Raw")
getwd()
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
plotPA(PAmatrix, title = "Presence absence plot for all comparisons")

###Calibrate your plot

#Modify selection only if you want to plot a specific comparison. 
#select what comparison you want here, leave the ones you dont care about at "<2".
#To show genes that are not upregulated in a certain comparison, Use "<1".
#To show genes that are not downregulated in a certain comparison, Use "<1".
#You can also use >comparison1 == comparison2 to find genes that are similarly regulated across comparisons.
#If you are working with a different dataset, make sure the number of columns is equal to the number of argments in 'selection'!
selection <- which    (PAmatrix[,1]  <2  &
                       PAmatrix[,2]  <2  & 
                       PAmatrix[,3]  <2  & 
                       PAmatrix[,4]  <2  & 
                       PAmatrix[,5]  <2  & 
                       PAmatrix[,6]  <2  & 
                       PAmatrix[,7]  <2  & 
                       PAmatrix[,8]  <2  & 
                       PAmatrix[,9]  <2  
)

###exclude a number from this factor to exclude the corresponding column from the PA plot. 
showTheseColumns <- c(1, 2, 4, 5)

###plot the adjusted presence/absence plot using the arguments 'selection' and 'showTheseColumns' defined above.
plotPA(PAmatrix[selection,], filter = showTheseColumns, title = "untitled")


###define for what column or columns you want to retrieve the genenames. 1 signifies the leftmost column. 
x <- c(1, 2)

###generate a dataframe with all genenames in the comparison defined in x. 
#Set mode to "similar", "different" or "simple"
#See line 92 for extra documentation
selectedNames <- getNamesOfComparisons(mode = "simple", x = 2)

#Write the genenames to a file:
write.csv(selectedNames, file="namesOfSelectedGenes.csv", row.names = FALSE, na="")

#To save the heatmap as an eps file:
outputFile = paste("heatmap_",title = "untitledPA",".eps", sep = "") 
dev.copy2eps(file = outputFile)
