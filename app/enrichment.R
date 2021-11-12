observe({
    updateSelectInput(session, "EnrichComp", choices=c(names(comparisons),'None'), selected='None')
})

# Gene Set Enrichment 
############################################################################

### Table part
GSEAReactive <- reactive({
    validate(need(input$EnrichComp!='None', message = "Please choose the comparison."))
    withProgress(message='Running Gene Set Enrichment',
        detail='Please stand by...',
        {
            shiny::setProgress(value = 0.15, detail = "Loading Reference Database")

            library(organism, character.only=TRUE)

            shiny::setProgress(value = 0.45, detail = "Pulling comparison data")

            DataSet=comparisons[[input$EnrichComp]]

            original_gene_list <- DataSet$log2FoldChange
            names(original_gene_list) <- DataSet$Gene
            gene_list<-na.omit(original_gene_list)
            gene_list = sort(gene_list, decreasing = TRUE)

            ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
            dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
            df2 = DataSet[DataSet$Gene %in% dedup_ids$ENSEMBL,]
            df2$Entrez = dedup_ids$ENTREZID
            kegg_gene_list <- df2$log2FoldChange
            names(kegg_gene_list) <- df2$Entrez
            kegg_gene_list<-na.omit(kegg_gene_list)
            kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

            shiny::setProgress(value = 0.65, detail = "Executing Gene Set Analysis")
            if(input$GSESet=='GO'){
                gse <- gseGO(geneList=gene_list, 
                    ont =input$GSEOnt, 
                    keyType = "ENSEMBL", 
                    nPerm = as.numeric(input$NPerm),
                    pvalueCutoff = 0.05,
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    verbose = TRUE, 
                    OrgDb = organism, 
                    pAdjustMethod = input$GSEPVal)

                return(gse)
            }
            else if(input$GSESet=='kegg'){
                gse <- gseKEGG(geneList=kegg_gene_list,
                    organism=korg,
                    nPerm = as.numeric(input$NPerm),
                    pvalueCutoff = 0.05,
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    verbose = TRUE,
                    pAdjustMethod = input$GSEPVal
                )

                return(gse)
            }
            else if(input$GSESet=='msigdb'){
                gene_sets=msigdbr(species=msig_org, category=input$GSEMSig)
                gene_sets$gs_name = gsub("_", " ", gene_sets$gs_name, fixed = TRUE)
                msigdbr_t2g = gene_sets %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()

                gse <- GSEA(kegg_gene_list, 
                    TERM2GENE=msigdbr_t2g, 
                    nPerm = as.numeric(input$NPerm),
                    pvalueCutoff = 0.05,
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    verbose = TRUE,
                    pAdjustMethod = input$GSEPVal
                )

                return(gse)
            }
        }
    )

})


output$gsea_table <- renderDataTable({
    DataIn <- GSEAReactive()
    DT::datatable(as.data.frame(DataIn), style = "bootstrap", options=list(pageLength = 15,scrollX=TRUE))
})

### Plot part
observe({
    output$gseaplot <- renderPlot({
        DataIn <- GSEAReactive()

        updateSelectizeInput(session,'GSEAPathwayChoice', choices=DataIn@result$Description)
        if(input$GSESet=='kegg'){
            updateSelectizeInput(session, 'GSEAKeggChoice', choices=DataIn@result$Description)
        }
        
        if(input$GSEAPlotType=='dotplot'){
            plot <- dotplot(DataIn, showCategory=as.numeric(input$NCategories), split=".sign", font.size=as.numeric(input$GSEFontSize)) + facet_grid(.~.sign) + scale_y_discrete(labels=function(x) str_wrap(x, width=as.numeric(input$GSEYWidth)))
            plot
        }
        else if(input$GSEAPlotType=='ridgechart'){
            plot <- ridgeplot(DataIn, showCategory = as.numeric(input$NCategories)) + labs(x = "enrichment distribution") + theme(text = element_text(size=as.numeric(input$GSEFontSize)))
            plot
        }
        else if(input$GSEAPlotType=='emap'){
            plot <- emapplot(DataIn, showCategory = as.numeric(input$NCategories), font.size=as.numeric(input$GSEFontSize))
            plot
        }
    }, height=input$EnrichHeight, width=input$EnrichWidth)
})

output$DownloadGSEATable <- downloadHandler(
    filename=function(){
        paste(input$EnrichComp,'_',input$GSESet,'_gsea','.csv',sep='')
    },
    content = function(file){
        write.csv(GSEAReactive(), file)
    }
)

output$DownloadGSEAPlot <- downloadHandler(
    filename=function(){
        paste(input$EnrichComp,'_',input$GSESet,'_gsea.',input$DownGSEAFormat,sep='')
    },
    content=function(file){   
        DataIn <- GSEAReactive()
        if(input$GSEAPlotType=='dotplot'){
            plot <- dotplot(DataIn, showCategory=as.numeric(input$NCategories), split=".sign", font.size=as.numeric(input$GSEFontSize)) + facet_grid(.~.sign) + scale_y_discrete(labels=function(x) str_wrap(x, width=as.numeric(input$GSEYWidth)))
        }
        else if(input$GSEAPlotType=='ridgechart'){
            plot <- ridgeplot(DataIn, showCategory = as.numeric(input$NCategories)) + labs(x = "enrichment distribution") + theme(text = element_text(size=as.numeric(input$GSEFontSize)))
        }
        else if(input$GSEAPlotType=='emap'){
            plot <- emapplot(DataIn, showCategory = as.numeric(input$NCategories), font.size=as.numeric(input$GSEFontSize))
        }
        if(input$DownGSEAFormat=='jpeg'){
            jpeg(file, height=input$EnrichHeight, width=input$EnrichWidth)
            print(plot)
            dev.off()
        }
        if(input$DownGSEAFormat=='png'){
            png(file, height=input$EnrichHeight, width=input$EnrichWidth)
            print(plot)
            dev.off()
        }
        if(input$DownGSEAFormat=='tiff'){
            tiff(file, height=input$EnrichHeight, width=input$EnrichWidth)
            print(plot)
            dev.off()
        }
    }
)

# ORA
############################################################################
ORAReactive <- reactive({
    validate(need(input$EnrichComp!='None', message = "Please choose the comparison."))
    withProgress(message='Running Over Representation',
        detail='Please stand by...',
        {
            shiny::setProgress(value = 0.15, detail = "Loading Reference Database")

            library(organism, character.only=TRUE)

            shiny::setProgress(value = 0.45, detail = "Pulling comparison data")

            DataSet=comparisons[[input$EnrichComp]]

            if(input$Enrichpval =='padj'){
                DataSet = subset(DataSet, padj < 0.05)
            }
            else if (input$Enrichpval=='pvalue'){
                DataSet = subset(DataSet, pvalue < 0.05)
            }
            
            original_gene_list <- DataSet$log2FoldChange
            names(original_gene_list) <- DataSet$Gene
            gene_list<-na.omit(original_gene_list)
            gene_list = sort(gene_list, decreasing = TRUE)
            genes <- names(gene_list)

            ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
            dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
            df2 = DataSet[DataSet$Gene %in% dedup_ids$ENSEMBL,]
            df2$Entrez = dedup_ids$ENTREZID
            kegg_gene_list <- df2$log2FoldChange
            names(kegg_gene_list) <- df2$Entrez
            kegg_gene_list<-na.omit(kegg_gene_list)
            kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
            kegg_genes <- names(kegg_gene_list)

            shiny::setProgress(value = 0.65, detail = "Executing Over Representation Analysis")
            if(input$ORASet=='GO'){
                ora <- enrichGO(gene=genes, 
                    ont =input$ORAOnt, 
                    keyType = "ENSEMBL", 
                    pvalueCutoff = 0.05,
                    OrgDb = organism, 
                    pAdjustMethod = input$ORAPVal)

                return(ora)
            }
            else if(input$ORASet=='kegg'){
                ora <- enrichKEGG(gene=kegg_genes,
                    organism=korg,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = input$ORAPVal
                )

                return(ora)
            }
            else if(input$ORASet=='msigdb'){
                gene_sets=msigdbr(species=msig_org, category=input$ORAMSig)
                gene_sets$gs_name = gsub("_", " ", gene_sets$gs_name, fixed = TRUE)
                msigdbr_t2g = gene_sets %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()

                ora <- enricher(kegg_genes, 
                    TERM2GENE=msigdbr_t2g, 
                    pvalueCutoff = 0.05,
                    pAdjustMethod = input$ORAPVal
                )

                return(ora)
            }
        }
    )

})


output$ora_table <- renderDataTable({
    DataIn <- ORAReactive()
    DT::datatable(as.data.frame(DataIn), style = "bootstrap", options=list(pageLength = 15,scrollX=TRUE))
})

### Plot part
observe({
    output$oraplot <- renderPlot({
        DataIn <- ORAReactive()

        if(input$ORASet=='kegg'){
            updateSelectizeInput(session, 'ORAKeggChoice', choices=DataIn@result$Description)
        }
        if(input$ORAPlotType=='dotplot'){
            plot <- dotplot(DataIn, showCategory=as.numeric(input$ORANCategories), font.size=as.numeric(input$ORAFontSize)) + scale_y_discrete(labels=function(x) str_wrap(x, width=as.numeric(input$ORAYWidth)))
            plot
        }
        else if(input$ORAPlotType=='barchart'){
            plot <- barplot(DataIn, showCategory=as.numeric(input$ORANCategories),font.size=as.numeric(input$ORAFontSize))
            plot
        }
        else if(input$ORAPlotType=='emap'){
            plot <- emapplot(DataIn, showCategory = as.numeric(input$ORANCategories), font.size=as.numeric(input$ORAFontSize))
            plot
        }
    }, height=input$EnrichHeight, width=input$EnrichWidth)
})

output$DownloadORATable <- downloadHandler(
    filename=function(){
        paste(input$EnrichComp,'_',input$ORASet,'_ORA','.csv',sep='')
    },
    content = function(file){
        write.csv(ORAReactive(), file)
    }
)

output$DownloadORAPlot <- downloadHandler(
    filename=function(){
        paste(input$EnrichComp,'_',input$ORASet,'_ORA.',input$DownORAFormat,sep='')
    },
    content=function(file){   
        DataIn <- ORAReactive()
        if(input$ORAPlotType=='dotplot'){
            plot <- dotplot(DataIn, showCategory=as.numeric(input$ORANCategories), font.size=as.numeric(input$ORAFontSize)) + scale_y_discrete(labels=function(x) str_wrap(x, width=as.numeric(input$ORAYWidth)))
        }
        else if(input$ORAPlotType=='barchart'){
            plot <- barplot(DataIn, showCategory=as.numeric(input$ORANCategories),font.size=as.numeric(input$ORAFontSize))
        }
        else if(input$ORAPlotType=='emap'){
            plot <- emapplot(DataIn, showCategory = as.numeric(input$ORANCategories), font.size=as.numeric(input$ORAFontSize))
        }
        if(input$DownORAFormat=='jpeg'){
            jpeg(file, height=input$EnrichHeight, width=input$EnrichWidth)
            print(plot)
            dev.off()
        }
        if(input$DownORAFormat=='png'){
            png(file, height=input$EnrichHeight, width=input$EnrichWidth)
            print(plot)
            dev.off()
        }
        if(input$DownORAFormat=='tiff'){
            tiff(file, height=input$EnrichHeight, width=input$EnrichWidth)
            print(plot)
            dev.off()
        }
    }
)

# Single Pathways
############################################################################

observe({
    output$gsea_soloplot <- renderPlot({
        DataIn <- GSEAReactive()
        DF <- as.data.frame(DataIn)
        index <- which(DF$Description==input$GSEAPathwayChoice)

        plot <- gseaplot(DataIn, geneSetID=as.numeric(index), by="runningScore", title=input$GSEAPathwayChoice)
        plot
        
    }, height=input$EnrichHeight, width=input$EnrichWidth)
})

output$DownloadGSEASoloPlot <- downloadHandler(
    filename=function(){
        paste(input$EnrichComp,input$GSEAPathwayChoice,'.',input$DownGSEASoloFormat,sep='')
    },
    content=function(file){   
        DataIn <- GSEAReactive()
        DF <- as.data.frame(DataIn)
        index <- which(DF$Description==input$GSEAPathwayChoice)

        plot <- gseaplot(DataIn, geneSetID=as.numeric(index), by="runningScore", title=input$GSEAPathwayChoice)

        if(input$DownGSEASoloFormat=='jpeg'){
            jpeg(file, height=input$EnrichHeight, width=input$EnrichWidth)
            print(plot)
            dev.off()
        }
        if(input$DownGSEASoloFormat=='png'){
            png(file, height=input$EnrichHeight, width=input$EnrichWidth)
            print(plot)
            dev.off()
        }
        if(input$DownGSEASoloFormat=='tiff'){
            tiff(file, height=input$EnrichHeight, width=input$EnrichWidth)
            print(plot)
            dev.off()
        }
    }
)

KEGGReactiveGSEA <- eventReactive(input$MakePathviewGSEA, {
    validate(need(input$GSESet=='kegg', message = "You must select KEGG pathways."))

    withProgress(message='Running Pathview',
        detail='Please stand by...',
        {
            shiny::setProgress(value = 0.15, detail = "Loading databases")

            shiny::setProgress(value = 0.35, detail = "Loading comparisons")

            DataSet=comparisons[[input$EnrichComp]]

            shiny::setProgress(value = 0.65, detail = "Pulling genelists")
            original_gene_list <- DataSet$log2FoldChange
            names(original_gene_list) <- DataSet$Gene
            gene_list<-na.omit(original_gene_list)
            gene_list = sort(gene_list, decreasing = TRUE)

            ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
            dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
            df2 = DataSet[DataSet$Gene %in% dedup_ids$ENSEMBL,]
            df2$Entrez = dedup_ids$ENTREZID
            kegg_gene_list <- df2$log2FoldChange
            names(kegg_gene_list) <- df2$Entrez
            kegg_gene_list<-na.omit(kegg_gene_list)
            kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

            DataIn <- GSEAReactive()
            DF <- as.data.frame(DataIn)
            index <- which(DF$Description==input$GSEAKeggChoice)
            print(input$GSEAKeggChoice)
            keggid <- DF[as.numeric(index), "ID"]
            print(keggid)

            shiny::setProgress(value = 0.85, detail = "Making plots")
            plot <- pathview(gene.data=kegg_gene_list, 
                pathway.id=keggid,
                species=korg,
                high=list(gene='red', cpd='red'),
                low=list(gene='blue',cpd='blue')
                )
            file.copy(paste0(keggid,".pathview.png"),paste0("tempimage"))
            TempValue$ImagePath=paste0(keggid,".pathview.png")
            return(list(
                src = paste0("tempimage"),
                filetype = "image/png",
                alt = "Pathview Image"
            ))
    })
})

TempValue<- reactiveValues(ImagePath=NULL)

observe({
    output$gsea_keggplot <- renderImage({
        KEGGReactiveGSEA()
    }, deleteFile=TRUE)
})

output$DownloadGSEAKEGGPlot <- downloadHandler(
    filename=function(){
        TempValue$ImagePath
    },
    content=function(file){   
        print(paste0(getwd(),'/',TempValue$ImagePath))
        file.copy(paste0(getwd(),'/',TempValue$ImagePath), file)
    }
)


KEGGReactiveORA <- eventReactive(input$MakePathviewORA, {
    validate(need(input$ORASet=='kegg', message = "You must select KEGG pathways."))

    withProgress(message='Running Pathview',
        detail='Please stand by...',
        {
            shiny::setProgress(value = 0.15, detail = "Loading databases")

            shiny::setProgress(value = 0.35, detail = "Loading comparisons")

            DataSet=comparisons[[input$EnrichComp]]

            shiny::setProgress(value = 0.65, detail = "Pulling genelists")
            original_gene_list <- DataSet$log2FoldChange
            names(original_gene_list) <- DataSet$Gene
            gene_list<-na.omit(original_gene_list)
            gene_list = sort(gene_list, decreasing = TRUE)

            ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
            dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
            df2 = DataSet[DataSet$Gene %in% dedup_ids$ENSEMBL,]
            df2$Entrez = dedup_ids$ENTREZID
            kegg_gene_list <- df2$log2FoldChange
            names(kegg_gene_list) <- df2$Entrez
            kegg_gene_list<-na.omit(kegg_gene_list)
            kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

            DataIn <- ORAReactive()
            DF <- as.data.frame(DataIn)
            index <- which(DF$Description==input$ORAKeggChoice)
            keggid <- DF[as.numeric(index), "ID"]

            shiny::setProgress(value = 0.85, detail = "Making plots")
            plot <- pathview(gene.data=kegg_gene_list, 
                pathway.id=keggid,
                species=korg,
                high=list(gene='red', cpd='red'),
                low=list(gene='blue',cpd='blue')
                )
            file.copy(paste0(keggid,".pathview.png"),paste0("tempimage"))
            TempValue$ImagePathORA=paste0(keggid,".pathview.png")
            return(list(
                src = paste0("tempimage"),
                filetype = "image/png",
                alt = "Pathview Image"
            ))
    })
})

observe({
    output$ora_keggplot <- renderImage({
        KEGGReactiveORA()
    }, deleteFile=TRUE)
})

output$DownloadGSEAORAPlot <- downloadHandler(
    filename=function(){
        TempValue$ImagePathORA
    },
    content=function(file){   
        print(paste0(getwd(),'/',TempValue$ImagePathORA))
        file.copy(paste0(getwd(),'/',TempValue$ImagePathORA), file)
    }
)