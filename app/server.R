function(input, output, session) {
    #load('data.RData')

    options(shiny.maxRequestSize=30*1024^2)
    source('ui.R',local=TRUE)
}