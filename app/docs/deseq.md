# Example Prep for data object
## Load processing libraries
```
library(DESeq2)
library(tidyr)
library(tidyverse)
library(dplyr)
```

## Load the merged gene counts file and metadata
These should both be provided to prep the object. The merged gene counts file will come from featureCounts currently. Future updates will include tximport for Salmon data as well. 
Metadata should be provided by the investigator in relation to the specific project. 
```
#Load the file and organize it
filecontents <- read.csv('./merged_gene_counts.txt', header=T, sep='\t', row.names=1)
counts <- filecontents[, -c(1)]
counts <- counts[,order(names(counts))]

#Create a quick table to parse the ensembl gene to gene symbol
genetable <- subset(filecontents, select=c(1))
names(genetable)[1] <- 'Geneid'
genetable$Gene <- rownames(genetable)

#Load the metadata and organize as well
filecontents = read.csv('./meta.csv', header=T, sep=',', row.names=1)
sampletable <- filecontents[order(row.names(filecontents)),]
countdata <- as.matrix(counts)[, colnames(counts) %in% rownames(sampletable)]

#Set your factor for comparisons
sampletable$condition <- factor(sampletable$Group)

#Run DESeq to create the dds object. You can customize the design here if it is multifactorial. 
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=sampletable, design=~ condition)
dds <- DESeq(dds)
```

## Create the normalized assay data
Pull the VST normalized data from the dds object. You can also do this with rld if you prefer.
```
vst <- vst(dds,blind=TRUE)

mat <- assay(vst)
#If you need to regress out a batch effect
#mat <- limma::removeBatchEffect(mat, vst$batch)
#assay(vst) <- mat
```

## Set up the cross comparisons and run those real quick
```
#Define the list of objects
comparisons <- list()

#Set up the cross comparisons based on your design. This example is for the demo data. 
compares <- list(
    list('Cell1Fungi','Control'),
    list('Cell2Fungi','Control'),
    list('Cell1Bug','Control'),
    list('Cell2Bug','Control')
)

#Loop through each and create a nice output dataframe. 
for(i in compares){
    out_res <- results(dds, contrast=c('condition',i[[1]], i[[2]]), alpha=0.05)
    shrink <- lfcShrink(dds, contrast=c('condition',i[[1]], i[[2]]),res=out_res, type="normal")
    comp_out <- shrink[order(shrink$padj),]
    comp_out <- merge(as.data.frame(comp_out), as.data.frame(mat),
    by='row.names', sort=FALSE)
    names(comp_out)[1] <- 'Gene'

    DataIn <- merge(as.data.frame(comp_out), as.data.frame(genetable), by="Gene", sort=FALSE)
    DataIn <- DataIn %>% relocate(Geneid, .after = Gene)
    comparisons[[paste(i[[1]],'_vs_',i[[2]], sep='')]] <- DataIn
}
```

## Save it as an output image to load.
```
save.image('data.RData')
```