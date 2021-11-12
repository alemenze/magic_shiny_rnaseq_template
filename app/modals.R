# DESeq2 processing
########################################################################################################
observeEvent(input$deseq_processing, {
    showModal(modalDialog(
        column(12,includeMarkdown("docs/deseq.md"), align='center', hr()),
        easyClose = TRUE,
        footer = modalButton("Close"),
        size='l'
    ))    
})