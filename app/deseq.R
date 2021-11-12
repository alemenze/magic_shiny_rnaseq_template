# Table Output
############################################################################

observe({
    updateSelectInput(session, "DESeqTable", choices=names(comparisons))
})

DESeqTableReactive <- reactive({
    datain <- comparisons[[paste(input$DESeqTable)]]

    return(datain)
})

output$deseq_counts <- renderDataTable({
    DT::datatable(DESeqTableReactive(), style='bootstrap',options=list(pageLength = 15,scrollX=TRUE))
})

output$DESeqDownload <- downloadHandler(
    filename=function(){
        paste(paste(input$DESeqTable),'.csv',sep='')
    },
    content = function(file){
        write.csv(DESeqTableReactive(), file)
    }
)

output$TotalSig <- renderText({
    length(rownames(subset(DESeqTableReactive(), padj < 0.05)))
})

output$UpSig <- renderText({
    length(rownames(subset(DESeqTableReactive(), (padj < 0.05 & log2FoldChange > 0))))
})

output$DownSig <- renderText({
    length(rownames(subset(DESeqTableReactive(), (padj < 0.05 & log2FoldChange < 0))))
})
