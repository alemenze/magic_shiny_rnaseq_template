# Venn
############################################################################

observe({
    updateSelectInput(session, "VennComp1", choices=c(names(comparisons),'None'), selected='None')
    updateSelectInput(session, "VennComp2", choices=c(names(comparisons),'None'), selected='None')
    updateSelectInput(session, "VennComp3", choices=c(names(comparisons),'None'), selected='None')
    updateSelectInput(session, "VennComp4", choices=c(names(comparisons),'None'), selected='None')
    updateSelectInput(session, "VennComp5", choices=c(names(comparisons),'None'), selected='None')
})

static_venn_plotter <- reactive({
    DataSetIn=comparisons

    validate(need(input$VennComp1 !='None', message='You must select at least one group'))
    

    sigger <- function(list_choice){
        sig_genes <- c()

        DataSet=DataSetIn[[paste(list_choice)]]

        if (input$Vennpval=='padj'){
            dtmp <- as.data.frame(DataSet)
            gtmp <- dtmp %>% filter(padj < 0.05)
            glist <- as.list(gtmp['Gene'])
            for (item in glist$Gene){
                sig_genes <- c(sig_genes, item)
            }   
        }
        if (input$Vennpval=='pvalue'){
            dtmp <- as.data.frame(DataSet)
            gtmp <- dtmp %>% filter(pvalue < 0.05)
            glist <- as.list(gtmp['Gene'])
            for (item in glist$Gene){
                sig_genes <- c(sig_genes, item)
            }
        }

        sig_genes <- unique(sig_genes)
        return(sig_genes)
    }

    vennlist <- list()
    venncolors <- list()

    if(input$VennComp1!='None'){
        Venn1 <- sigger(input$VennComp1)
        vennlist[[input$VennComp1]] <- Venn1
        venncolors[[input$VennComp1]] <- input$color1
    }
    if(input$VennComp2!='None'){
        Venn2 <- sigger(input$VennComp2)
        vennlist[[input$VennComp2]] <- Venn2
        venncolors[[input$VennComp2]] <- input$color2
    }
    if(input$VennComp3!='None'){
        Venn3 <- sigger(input$VennComp3)
        vennlist[[input$VennComp3]] <- Venn3
        venncolors[[input$VennComp3]] <- input$color3
    }
    if(input$VennComp4!='None'){
        Venn4 <- sigger(input$VennComp4)
        vennlist[[input$VennComp4]] <- Venn4
        venncolors[[input$VennComp4]] <- input$color4
    }
    if(input$VennComp5!='None'){
        Venn5 <- sigger(input$VennComp5)
        vennlist[[input$VennComp5]] <- Venn5
        venncolors[[input$VennComp5]] <- input$color5
    }

    fillcolors <- unlist(venncolors)

    vennplot <- venn.diagram(x=vennlist,
        fill=fillcolors,
        lwd=input$VennLineThick,
        lty=input$VennLineType,
        cat.cex=input$VennLabelFontSize,
        cex=input$VennFontSize,
        alpha=input$VennAlpha,
        filename=NULL
    )
    return(vennplot)
})

observe({
    output$venn_static <- renderPlot({
        grid.draw(static_venn_plotter())
    }, height=input$VennHeight, width=input$VennWidth)
})

output$DownloadVennS <- downloadHandler(
    filename=function(){
        paste('venn',input$DownVennSFormat,sep='.')
    },
    content=function(file){   
        if(input$DownVennSFormat=='jpeg'){
            jpeg(file, height=input$VennHeight, width=input$VennWidth)
            print(grid.draw(static_venn_plotter()))
            dev.off()
        }
        if(input$DownVennSFormat=='png'){
            png(file, height=input$VennHeight, width=input$VennWidth)
            print(grid.draw(static_venn_plotter()))
            dev.off()
        }
        if(input$DownVennSFormat=='tiff'){
            tiff(file, height=input$VennHeight, width=input$VennWidth)
            print(grid.draw(static_venn_plotter()))
            dev.off()
        }
    }
)

# UpSet
############################################################################
observe({
    updateSelectizeInput(session, 'upset_compares', choices=names(comparisons), server=TRUE)  
})

upset_plotter <- reactive({
    validate(need(length(input$upset_compares)>1, message = "Please choose at least 2 comparisons."))
    DataSetIn=comparisons

    sigger <- function(list_choice){
        sig_genes <- c()

        DataSet=DataSetIn[[list_choice]]

        if (input$Upsetpval=='padj'){
            dtmp <- as.data.frame(DataSet)
            gtmp <- dtmp %>% filter(padj < 0.05)
            glist <- as.list(gtmp['Gene'])
            for (item in glist$Gene){
                sig_genes <- c(sig_genes, item)
            }   
        }
        if (input$Upsetpval=='pvalue'){
            dtmp <- as.data.frame(DataSet)
            gtmp <- dtmp %>% filter(pvalue < 0.05)
            glist <- as.list(gtmp['Gene'])
            for (item in glist$Gene){
                sig_genes <- c(sig_genes, item)
            }
        }

        sig_genes <- unique(sig_genes)
        return(sig_genes)
    }

    vennlist <- list()

    for(item in input$upset_compares){
        vennlist[[paste(item)]] <- sigger(item)
    }

    upset(fromList(vennlist), order.by='freq',
        text.scale=c(input$UpSetintsizetitle, input$UpSetintsizeticklabels, input$UpSetsetsizetitle,
            input$UpSetsetsizeticklabels, input$UpSetsetnames, input$UpSetnumbersonbar),
        point.size=input$UpSetPointSize, 
        line.size=input$UpSetLineSize,
        mainbar.y.label='Gene intersections', sets.x.label='Significant genes',
        )

})


observe({
    output$upset_static <- renderPlot({
        upset_plotter()
    }, height=input$VennHeight, width=input$VennWidth)
})

output$DownloadUpset <- downloadHandler(
    filename=function(){
        paste('upset',input$DownUpsetFormat,sep='.')
    },
    content=function(file){   
        if(input$DownUpsetFormat=='jpeg'){
            jpeg(file, height=input$VennHeight, width=input$VennWidth)
            print(upset_plotter())
            dev.off()
        }
        if(input$DownUpsetFormat=='png'){
            png(file, height=input$VennHeight, width=input$VennWidth)
            print(upset_plotter())
            dev.off()
        }
        if(input$DownUpsetFormat=='tiff'){
            tiff(file, height=input$VennHeight, width=input$VennWidth)
            print(upset_plotter())
            dev.off()
        }
    }
)