# beaveR

# Installation

`devtools::install_github("NPSDC/beaveR")`

## Example of how to process parse the output of TreeTerminus and obtain fixed groups from it

### Creating a treeSummarizedExperiment object
```{r}
clustFile <- file.path(treeTermDir, 'cluster_nwk.txt') #path to file output by TreeTerminus
quantDir <- file.path(dir, 'out_sal') # path to directory containing Salmon quantified files for an RNASeq experiment
samples <- as.vector(outer(c(1:6), c(1,2), function(x,y) paste(x,y,sep='_')))
quantFiles <- file.path(quantDir, samples, 'quant.sf')
coldata <- data.frame(files=quantFiles, names=samples, condition=factor(rep(1:2, each=6)))
tse <- buildTSE(treeTermFile = clustFile, coldata = coldata)
print(tse)
```

### Find the summary information associated with a node
```{r}
node <- 300
nodeInf <- findNodeInformation(tse, node=node, type='children')
print(nodeInf)
```

### Find the groups associated with the objective function that minimizes mean infRV and height
```{r}
tree <- TreeSummarizedExperiment::rowTree(tse)
gamma <- 0.1
descSize <- sapply(phangorn::Descendants(tree, seq(nrow(tse))), length)
metric <- (mcols(tse)[['meanInfRV']] + ape::node.depth(tree, 2)*gamma) *descSize
objS <- solveForOptimalCut(tse, metVec = metric, type = 'min')
print(objS[['optVal']])
print(length(objS[['cut']]))
```

### Find the groups associated with the objective function that maximizes weighted log fold change
```{r}
tree <- TreeSummarizedExperiment::rowTree(tse)
lfc <- getScaledLFC(tse, "condition")
descSize <- sapply(phangorn::Descendants(tree, seq(nrow(tse))), length)
metric <- lfc/mcols(tse)[['meanInfRV']]*descSize
objS <- solveForOptimalCut(tse, metVec = metric, type = 'max')
print(objS[['optVal']])
print(length(objS[['cut']]))
```
