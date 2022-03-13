# PCA Plots 
############################################################################

observe({
    updateSelectInput(session, "PCAColor", choices=colnames(sampletable), selected='Group')
    updateSelectInput(session, "PCAShape", choices=c("None",colnames(sampletable)), selected='None')
})

pcaplotter <- reactive({
    p <- pca(mat, metadata=sampletable, removeVar=0.1)
    if (input$PCAShape=="None"){
        pcaplot <- biplot(p, colby=input$PCAColor, 
            legendPosition=input$PCALegendPos, legendLabSize=input$PCALegendLabel, encircle=input$Ellipsoid, legendIconSize=input$PCALegendIcon, pointSize=input$PCAPointSize, labSize=input$PCALabelSize)
    } else {
        pcaplot <- biplot(p, colby=input$PCAColor, shape=input$PCAShape,
            legendPosition=input$PCALegendPos, legendLabSize=input$PCALegendLabel, encircle=input$Ellipsoid, legendIconSize=input$PCALegendIcon, pointSize=input$PCAPointSize, labSize=input$PCALabelSize)
    }
    return(pcaplot)
})

observe({
    output$pcaplot <- renderPlot({
        pcaplotter()
    }, height=input$CHeight, width=input$CWidth)
})

output$DownloadPCA <- downloadHandler(
    filename=function(){
        paste('PCA',input$DownPCAFormat,sep='.')
    },
    content=function(file){   
        if(input$DownPCAFormat=='jpeg'){
            jpeg(file, height=input$CHeight, width=input$CWidth)
            print(pcaplotter())
            dev.off()
        }
        if(input$DownPCAFormat=='png'){
            png(file, height=input$CHeight, width=input$CWidth)
            print(pcaplotter())
            dev.off()
        }
        if(input$DownPCAFormat=='tiff'){
            tiff(file, height=input$CHeight, width=input$CWidth, res=1000)
            print(pcaplotter())
            dev.off()
        }
    }
)

# Distance matrix
############################################################################
distanceplotter <- reactive({
    sampleDists <- as.matrix(dist(t(mat)))

    anno_cols <- which(names(sampletable)==input$PCAColor)
    annotations=sampletable[,anno_cols, drop=F]

    colors <- colorRampPalette(rev(brewer.pal(9, "Greys")))(30)

    plot <- pheatmap(sampleDists, annotation_row=annotations, annotation_col=annotations, col=colors,
        fontsize_col=input$DMXsize, fontsize_row=input$DMYsize,angle_col=input$DMang
    )
    return(plot)
})

observe({
    output$distplot <- renderPlot({
        distanceplotter()
    }, height=input$CHeight, width=input$CWidth)
})

output$DownloadDM <- downloadHandler(
    filename=function(){
        paste('distance_matrix',input$DownDMFormat,sep='.')
    },
    content=function(file){   
        if(input$DownDMFormat=='jpeg'){
            jpeg(file, height=input$CHeight, width=input$CWidth)
            print(distanceplotter())
            dev.off()
        }
        if(input$DownDMFormat=='png'){
            png(file, height=input$CHeight, width=input$CWidth)
            print(distanceplotter())
            dev.off()
        }
        if(input$DownDMFormat=='tiff'){
            tiff(file, height=input$CHeight, width=input$CWidth, res=1000)
            print(distanceplotter())
            dev.off()
        }
    }
)

# Eigencorplots
############################################################################
observe({
    updateSelectizeInput(session, "EigenMetavars", choices=colnames(sampletable))
})
eigencorplotter <- reactive({
    
    p <- pca(mat, metadata=sampletable, removeVar=0.1)
    
    validate(need(length(input$EigenMetavars)>=2, message='You must select at least 2 metadata columns'))
    
    eplot <- eigencorplot(p, 
        metavars=input$EigenMetavars,
        cexCorval=input$ELabSize,
        cexLabX=input$EXLabSize,
        cexLabY=input$EYLabSize,
        colCorval = 'white',
        fontCorval = 2,
        posLab = 'bottomleft',
        rotLabX = 45,
        posColKey = input$ELegendPos,
        cexLabColKey = input$ELabColSize,
        scale = TRUE,
        main = 'PC1-10 metadata correlations',
        colFrame = 'white',
        plotRsquared = input$ERsquare,
        corFUN = input$EFun,
        corUSE = 'pairwise.complete.obs',
        corMultipleTestCorrection = input$ECorrect,
        signifSymbols = c('***', '**', '*', ''),
        signifCutpoints = c(0, 0.0001, 0.001, 0.01, 1))
    return(eplot)
})

observe({
    output$eigencorplotout <- renderPlot({
        eigencorplotter()
    }, height=input$CHeight, width=input$CWidth)
})

output$DownloadEigen <- downloadHandler(
    filename=function(){
        paste('eigencorplot',input$DownEigenFormat,sep='.')
    },
    content=function(file){   
        if(input$DownEigenFormat=='jpeg'){
            jpeg(file, height=input$CHeight, width=input$CWidth)
            print(eigencorplotter())
            dev.off()
        }
        if(input$DownEigenFormat=='png'){
            png(file, height=input$CHeight, width=input$CWidth)
            print(eigencorplotter())
            dev.off()
        }
        if(input$DownEigenFormat=='tiff'){
            tiff(file, height=input$CHeight, width=input$CWidth, res=1000)
            print(eigencorplotter())
            dev.off()
        }
    }
)