observe({
    GenesList=genetable$Geneid
    updateSelectizeInput(session, "HSelectedGenes", choices=GenesList, server=TRUE, options = list(maxOptions = 50))
})
observe({
    samplenames=colnames(mat)
    updateSelectizeInput(session, 'HSamples', choices=samplenames, server=TRUE)
})
observe({
    updateSelectInput(session, 'HGroups', choices=c("None",colnames(sampletable)), selected='None')
})

static_heatmap_plotter <- reactive({
    withProgress(message='Pulling heatmap data',
        detail='Please stand by...',{
    memory.limit(size=25000)

    cluster_r=TRUE
    if (input$HRowClust=='FALSE'){
        cluster_r=FALSE
    }
    cluster_c=TRUE
    if (input$HColClust=='FALSE'){
        cluster_c=FALSE
    }

    shownames=FALSE
    if (input$HShowNames=='TRUE'){
        shownames=TRUE
    }

    my_colors = brewer.pal(n = 11, name = "RdBu")
    my_colors = colorRampPalette(my_colors)(50)
    my_colors = rev(my_colors)
    
    shiny::setProgress(value = 0.55, detail = "Pulling genes")

    if(input$HSubType=='sigdif'){
        DataSetIn=comparisons

        sig_genes <- c()

        if (input$HSigP=='padj'){
            for (i in names(DataSetIn)){
                dtmp <- as.data.frame(DataSetIn[[i]])
                gtmp <- dtmp %>% filter(padj < 0.05)
                glist <- as.list(gtmp['Gene'])
                for (item in glist$Gene){
                    sig_genes <- c(sig_genes, item)
                }
            }
        }
        if (input$HSigP=='pvalue'){
            for (i in names(DataSetIn)){
                dtmp <- as.data.frame(DataSetIn[[i]])
                gtmp <- dtmp %>% filter(pvalue < 0.05)
                glist <- as.list(gtmp['Gene'])
                for (item in glist$Gene){
                    sig_genes <- c(sig_genes, item)
                }
            }
        }

        sig_genes <- unique(sig_genes)

        if(input$HSamp=='TRUE'){
            validate(need(length(input$HSamples)>2, message = "Please choose at least 2 samples."))
            DataIn <- subset(mat, TRUE, c(input$HSamples))
        }
        else{
            DataIn=mat
        }

        if(input$HGrouping=='TRUE'){
            heads <- list()
            opts <- c()
            colchoice=input$HGroups
            for(i in unique(sampletable[[colchoice]])){
                new <- unique(list(subset(as.data.frame(sampletable), eval(parse(text=colchoice))==i)$SampleName))
                opts <- c(opts, i)
                
                heads[[paste(i)]] <- list(paste(i), c(levels(factor(new[[1]]))))
            }

            temp <- as.data.frame(mat)
            for(i in heads){
                temp[[paste(i[[1]])]] <- rowMeans(temp[,c(i[[2]])])
            }

            temp <- temp[,colnames(temp) %in% opts]
            DataIn=temp
        }

        DTmpIn <- DataIn[rownames(DataIn) %in% sig_genes,]

        plot <- pheatmap(DTmpIn, cluster_rows=cluster_r, cluster_cols=cluster_c, color=my_colors,
            scale=input$HScale, show_rownames=shownames, fontsize_col=input$HXsize, angle_col=input$Hang)
        
        return(plot)
    }
    if(input$HSubType=='Select'){
        validate(need(length(input$HSelectedGenes)>1, message = "Please choose at least 2 genes."))
        chosen_genes <- c(input$HSelectedGenes)

        if(input$HSamp=='TRUE'){
            validate(need(length(input$HSamples)>2, message = "Please choose at least 2 samples."))
            DataIn <- subset(mat, TRUE, c(input$HSamples))
        }
        else{
            DataIn=mat
        }

        if(input$HGrouping=='TRUE'){
            heads <- list()
            opts <- c()
            colchoice=input$HGroups
            for(i in unique(sampletable[[colchoice]])){
                new <- unique(list(subset(as.data.frame(sampletable), eval(parse(text=colchoice))==i)$SampleName))
                opts <- c(opts, i)
                
                heads[[paste(i)]] <- list(paste(i), c(levels(factor(new[[1]]))))
            }

            temp <- as.data.frame(mat)
            for(i in heads){
                temp[[paste(i[[1]])]] <- rowMeans(temp[,c(i[[2]])])
            }

            temp <- temp[,colnames(temp) %in% opts]
            DataIn=temp
        }

        DataSet <- merge(as.data.frame(DataIn), as.data.frame(genetable), by="row.names", sort=FALSE)
        DataSet <- subset(DataSet, select=-c(Gene,Row.names))
        DataSet <- distinct(DataSet, Geneid, .keep_all=TRUE)
        DataSet <- DataSet %>% remove_rownames %>% column_to_rownames(var="Geneid")

        chosen_genes <- unique(chosen_genes)
        DTmpIn <- DataSet[rownames(DataSet) %in% chosen_genes,]
        plot <- pheatmap(DTmpIn, cluster_rows=cluster_r, cluster_cols=cluster_c, color=my_colors,
            scale=input$HScale, show_rownames=shownames, fontsize_col=input$HXsize, angle_col=input$Hang)
        
        return(plot)
    }

    })
})

observe({
    output$heatmap_static <- renderPlot({
        static_heatmap_plotter()
    },height=input$HHeight, width=input$HWidth)
})

output$DownloadHS <- downloadHandler(
    filename=function(){
        paste('heatmap',input$DownHSFormat,sep='.')
    },
    content=function(file){   
        if(input$DownHSFormat=='jpeg'){
            jpeg(file, height=input$CHeight, width=input$CWidth)
            print(static_heatmap_plotter())
            dev.off()
        }
        if(input$DownHSFormat=='png'){
            png(file, height=input$CHeight, width=input$CWidth)
            print(static_heatmap_plotter())
            dev.off()
        }
        if(input$DownHSFormat=='tiff'){
            tiff(file, height=input$CHeight, width=input$CWidth, res=1000)
            print(static_heatmap_plotter())
            dev.off()
        }
    }
)