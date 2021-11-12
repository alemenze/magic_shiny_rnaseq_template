observe({
    GenesList=genetable$Geneid
    updateSelectizeInput(session, "gene_select", choices=GenesList, server=TRUE, options = list(maxOptions = 50))   
})

observe({
    updateSelectInput(session, "GEGroup", choices=colnames(sampletable))
})

# Box
############################################################################
box_plotter <- reactive({
    validate(need(length(input$gene_select)>0, message = "Please choose at least 1 gene."))
    validate(need(length(input$GEGroup)>0, message = "Please choose the grouping."))

    DataIn=dds

    TmpDF=genetable
    index <- which(TmpDF$Geneid==input$gene_select)
    ens_gene_select <- TmpDF[as.numeric(index), "Gene"]

    filtered <- t(log2((counts(DataIn[ens_gene_select, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
           merge(colData(DataIn), ., by="row.names") %>%
           gather(gene, expression, (ncol(.)-length(ens_gene_select)+1):ncol(.))
    
    p <- ggplot(filtered, aes_string(input$GEGroup, "expression", fill=input$GEGroup))
    p <- p + geom_boxplot() + facet_wrap(~gene, scales="free_y")
    
    p <- p + xlab(' ') + ylab('Normalized log2 Expression') + theme(
      plot.margin = unit(c(1,1,1,1), "cm"),
      axis.text.x = element_text(angle = as.numeric(input$GEang),size=as.numeric(input$GELabelSize)),
      axis.text.y = element_text(size=as.numeric(input$GELabelSize)),
      plot.title=element_text(size=as.numeric(input$GETitleSize)),
      legend.position=input$GELegendPos) + ggtitle(input$gene_select)
    
    p
})

observe({
    output$box_static <- renderPlot({
        box_plotter()
    }, height=input$GEHeight, width=input$GEWidth)
})

output$DownloadBox <- downloadHandler(
    filename=function(){
        paste('boxplot',input$DownBoxFormat,sep='.')
    },
    content=function(file){   
        if(input$DownBoxFormat=='jpeg'){
            jpeg(file, height=input$GEHeight, width=input$GEWidth)
            print(box_plotter())
            dev.off()
        }
        if(input$DownBoxFormat=='png'){
            png(file, height=input$GEHeight, width=input$GEWidth)
            print(box_plotter())
            dev.off()
        }
        if(input$DownBoxFormat=='tiff'){
            tiff(file, height=input$GEHeight, width=input$GEWidth)
            print(box_plotter())
            dev.off()
        }
    }
)

# Violin
############################################################################
violin_plotter <- reactive({
    validate(need(length(input$gene_select)>0, message = "Please choose at least 1 gene."))
    validate(need(length(input$GEGroup)>0, message = "Please choose the grouping."))

    DataIn=dds

    TmpDF=genetable
    index <- which(TmpDF$Geneid==input$gene_select)
    ens_gene_select <- TmpDF[as.numeric(index), "Gene"]

    filtered <- t(log2((counts(DataIn[ens_gene_select, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
           merge(colData(DataIn), ., by="row.names") %>%
           gather(gene, expression, (ncol(.)-length(ens_gene_select)+1):ncol(.))
    
    p <- ggplot(filtered, aes_string(input$GEGroup, "expression", fill=input$GEGroup))
    p <- p + geom_violin(trim = FALSE) + facet_wrap(~gene, scales="free_y")
    
    p <- p + xlab(' ') + ylab('Normalized log2 Expression') + theme(
      plot.margin = unit(c(1,1,1,1), "cm"),
      axis.text.x = element_text(angle = as.numeric(input$GEang),size=as.numeric(input$GELabelSize)),
      axis.text.y = element_text(size=as.numeric(input$GELabelSize)),
      plot.title=element_text(size=as.numeric(input$GETitleSize)),
      legend.position=input$GELegendPos) + ggtitle(input$gene_select)
    
    p
})

observe({
    output$violin_static <- renderPlot({
        violin_plotter()
    }, height=input$GEHeight, width=input$GEWidth)
})

output$DownloadViolin <- downloadHandler(
    filename=function(){
        paste('violinplot',input$DownViolinFormat,sep='.')
    },
    content=function(file){   
        if(input$DownViolinFormat=='jpeg'){
            jpeg(file, height=input$GEHeight, width=input$GEWidth)
            print(violin_plotter())
            dev.off()
        }
        if(input$DownViolinFormat=='png'){
            png(file, height=input$GEHeight, width=input$GEWidth)
            print(violin_plotter())
            dev.off()
        }
        if(input$DownViolinFormat=='tiff'){
            tiff(file, height=input$GEHeight, width=input$GEWidth)
            print(violin_plotter())
            dev.off()
        }
    }
)