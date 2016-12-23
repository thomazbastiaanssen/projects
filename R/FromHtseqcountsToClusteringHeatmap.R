###IMPORTANT###
#Always run line 1:75 before using this script, as it contains custom functions. 

###set up and load
options(stringsAsFactors = F)

### load heatmap libraries
library("RColorBrewer");
library("gplots");
library("Biobase");
library("biomaRt");
library("DESeq2");

### load custom wrapper functions
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


#DEplot1 takes one mandatory object: vsdMat; This should be produced using the DESeq2 function 'assay' and contains the Z-scores of the reads. 
#The optional objects are 'filter'; an index for the relevant genes for the plot (strongly recommended)
#groupcols; this defines the colour of the legend on top of the plot
#Don't forget to add a title with > title = "mytitle"!

DEplot1 <- function(vsdMat    = vsdMat,
                    filter    = sigED, 
                    title     = "",  
                    groupcols = brewer.pal(8, "Dark2")[1:length(unique(colData(dds)$condition))]){
  
  return(heatmap.2(vsdMat[sigED,], 
                   col          = col, trace="none", 
                   main         = title, margin=c(10, 10), 
                   ColSideColors= groupcols[colData(dds)$condition], 
                   scale        = "row"))
  
}

#DEplot2 takes two mandatory objects: carpet; this is the default name for the output of the DEplot1 wrapper.
#vsdMat; This should be produced using the DESeq2 function 'assay' and contains the Z-scores of the reads.
#The optional objects are 'filter'; an index for the relevant genes for the plot (strongly recommended)
#groupcols; this defines the colour of the legend on top of the plot
#Don't forget to add a title with > title = "mytitle"!

DEplot2 <- function(x         = carpet, 
                    vsdMat    = vsdMat, 
                    controls  = controls, 
                    filter    = sigED, 
                    title     = "", 
                    groupcols = brewer.pal(8, "Dark2")[1:length(unique(colData(dds)$condition))]){
  
  tempsorted <- names(sort(rowMeans(t(carpet$carpet)[,controls])))
  
  return(heatmap.2(data.matrix(vsdMat[sigED,])[tempsorted,], 
                   col = rev(col), 
                   trace="none", 
                   main = title, 
                   margin=c(10, 10), 
                   ColSideColors=groupcols[colData(dds)$condition], 
                   scale="row", 
                   Rowv=F, 
                   dendrogram = "col"  
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



### Run DESeq. make sure DESeqPrep has the right input file. 
### You most likely want to run the following code up until the last line for each of the relevant comparisons. 

dds      <- DESeqPrep(countsconpcons)               #DESeqPrep is a wrapper function, it is defined earlier in this document. 
dds      <- DESeq(dds)                              #DESeq will take a couple of seconds.
vsd      <- varianceStabilizingTransformation(dds)  #this one too, but it's faster.
vsdMat   <- assay(vsd)
res      <- results(dds)
sigED    <- which(res$padj <0.1 )


###Define your colours
col <- colorRampPalette(brewer.pal(10, 'RdBu'))(100)                          #this controls the colourscale used to signify expression levels.
groupcols <- brewer.pal(8, "Dark2")[1:length(unique(colData(dds)$condition))] #this controls the colour of the legend above the plot.

#Produce the unsorted heatmap

carpet   <- DEplot1(vsdMat, title="untitled")

#indicate the location of the control groups on the heatmap in a vector format (leftmost = 1)

controls <- c(1, 2, 3, 4)

#Produce the sorted heatmap

DEplot2(carpet, vsdMat, title="untitled", controls= controls)

#To save the heatmap as an eps file:
outputFile = paste("heatmap_",title = "untitled",".eps", sep = "") 
dev.copy2eps(file = outputFile)