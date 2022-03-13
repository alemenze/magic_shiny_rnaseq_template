observe({
    GenesList=genetable$Geneid
    updateSelectInput(session, "VolComp", choices=names(comparisons))
    updateSelectizeInput(session, "VSelectedGenes", choices=GenesList, server=TRUE, options = list(maxOptions = 50))
})

static_volcano_plotter <- reactive({
    DataSet=comparisons[[paste(input$VolComp)]]

    if (input$VLabelType=='None'){
        plot <- EnhancedVolcano(DataSet,
            lab='',
            title=input$VolComp,
            x='log2FoldChange',
            y=input$Vpval,
            legendPosition = input$VLegendPos,
            legendLabels=c('Not sig.','Log2FC',paste(input$Vpval),paste(input$Vpval, ' & log2FC', sep='')),
            pCutoff=0.05, 
            FCcutoff=1, 
            pointSize=input$VpointSize,
            labSize=input$VLabSize,
            legendLabSize=input$VLegLabSize,
            legendIconSize=input$VLegIconSize
            )
    
        return(plot)
    }
    if (input$VLabelType != "Select"){
        plot <- EnhancedVolcano(DataSet,
            lab=as.character(DataSet[[input$VLabelType]]),
            title=input$VolComp,
            x='log2FoldChange',
            y=input$Vpval,
            legendPosition = input$VLegendPos,
            legendLabels=c('Not sig.','Log2FC',paste(input$Vpval),paste(input$Vpval, ' & log2FC', sep='')),
            pCutoff=0.05, 
            FCcutoff=1, 
            pointSize=input$VpointSize,
            labSize=input$VLabSize,
            legendLabSize=input$VLegLabSize,
            legendIconSize=input$VLegIconSize
            )
    
        return(plot)
    } else {
        plot <- EnhancedVolcano(DataSet,
            lab=as.character(DataSet$Geneid),
            title=input$VolComp,
            selectLab=c(input$VSelectedGenes),
            x='log2FoldChange',
            y=input$Vpval,
            legendPosition = input$VLegendPos,
            legendLabels=c('Not sig.','Log2FC',paste(input$Vpval),paste(input$Vpval, ' & log2FC', sep='')),
            pCutoff=0.05, 
            FCcutoff=1, 
            pointSize=input$VpointSize,
            labSize=input$VLabSize,
            legendLabSize=input$VLegLabSize,
            legendIconSize=input$VLegIconSize
            )
    
        return(plot)
    }
    
})

observe({
    output$vol_static <- renderPlot({
        static_volcano_plotter()
    }, height=input$VHeight, width=input$VWidth)
})

output$DownloadVS <- downloadHandler(
    filename=function(){
        paste(input$VolComp,input$DownVSFormat,sep='.')
    },
    content=function(file){   
        if(input$DownVSFormat=='jpeg'){
            jpeg(file, height=input$VHeight, width=input$VWidth)
            print(static_volcano_plotter())
            dev.off()
        }
        if(input$DownVSFormat=='png'){
            png(file, height=input$VHeight, width=input$VWidth)
            print(static_volcano_plotter())
            dev.off()
        }
        if(input$DownVSFormat=='tiff'){
            tiff(file, height=input$VHeight, width=input$VWidth, res=1000)
            print(static_volcano_plotter())
            dev.off()
        }
    }
)
