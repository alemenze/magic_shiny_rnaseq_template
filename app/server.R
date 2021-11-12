function(input, output, session) {
    load('data.RData')

    options(shiny.maxRequestSize=30*1024^2)
    source('ui.R',local=TRUE)
    source('modals.R',local=TRUE)
    source('deseq.R',local=TRUE)
    source('clustering.R',local=TRUE)
    source('volcanoes.R',local=TRUE)
    source('heatmaps.R', local=TRUE)
    source('geplots.R', local=TRUE)
    source('venn.R', local=TRUE)
    source('enrichment.R', local=TRUE)
}