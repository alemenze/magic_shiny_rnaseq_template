library(shiny)
require(shinyjs)
library(shinythemes)
require(shinycssloaders)

library(tidyr)
library(tidyverse)
library(dplyr)
library(DT)
library(colourpicker)
library(RColorBrewer)
library(ggplot2)

library(DESeq2)
library(PCAtools) #need dev version for ellipses
library('pheatmap')
library('EnhancedVolcano')
library(VennDiagram)
library(UpSetR)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(msigdbr)
library(stringr)
library(pathview)

tagList(
    tags$head(
        #includeHTML(("www/GA.html")), If including google analytics
        tags$style(type = 'text/css','.navbar-brand{display:none;}')
    ),
    fluidPage(theme = shinytheme('yeti'),
            windowTitle = "Bulk RNAseq Analysis",
            titlePanel(
                fluidRow(
                column(2, tags$a(href='http://www.bioinformagic.io/', tags$img(height =75 , src = "MaGIC_Icon_0f344c.svg")), align = 'center'), 
                column(10, fluidRow(
                  column(10, h1(strong('RNA Seq Interactive Visualization Tool'), align = 'center')),
                  column(10, h2(strong('Project'), align = 'center'))))
                ),
                windowTitle = "Bulk RNAseq Analysis" ),
                tags$style(type='text/css', '.navbar{font-size:20px;}'),
                tags$style(type='text/css', '.nav-tabs{padding-bottom:20px;}'),
                tags$head(tags$style(".modal-dialog{ width:1300px}")),

        navbarPage(title ="", id='NAVTABS',

        ## Intro Page
##########################################################################################################################################################
            tabPanel('Introduction',
                fluidRow(
                  column(2,
                  ),
                  column(8,
                      column(12, align = "center", 
                          style="margin-bottom:25px;",
                          h2("Introduction")),
                      hr(),
                      column(12, align="center",
                          style="margin-bottom:50px;", 
                          column(2,),
                          column(8,markdown("Welcome to the RNAseq data explorer by [the Molecular and Genomics Informatics Core](http://www.bioinformagic.io).
                            This tool is designed to utilize bulk RNA seq data processed by the Core's pipelines and enable you to visualize the data however you please. 
                            Each instance is custimizable and can have additional modules added. Please contact the Core with any specific requests. 
                          ")),
                          column(2,)       
                      )
                  ),
                  column(2,
                  )
                ),

                fluidRow(
                  column(2,
                  ),
                  column(8,
                      column(12, align = "center", 
                          style="margin-bottom:25px;",
                          h2("Primary processing")),
                      hr(),
                      column(12, align="center",
                          style="margin-bottom:50px;", 
                          column(2,),
                          column(8,markdown("
                          Primary alignment is performed via [STAR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/) along with associated quality trimming and QC assessment. 

                          #Add more details here

                          ```
                          code block
                          ```
                          ")),
                          column(2,)       
                      )
                  ),
                  column(2,
                  )
                ),

                fluidRow(
                  column(2,
                  ),
                  column(8,
                      column(12, align = "center", 
                          style="margin-bottom:25px;",
                          h2("Secondary processing")),
                      hr(),
                      column(12, align="center",
                          style="margin-bottom:50px;", 
                          column(2,),
                          column(8,markdown("Secondary processing is performed via [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) in R. 
                            #More details on this too
                            "),
                            actionButton('deseq_processing','DESeq2 Example',class='btn btn-info',style="margin-top:15px;")
                            ),
                          column(2,)       
                      )
                  ),

                  column(2,
                  )
                )
               
            ),

        ## Comparison Tables
##########################################################################################################################################################
            tabPanel('Comparison Tables',
                fluidRow(
                    column(3,
                        wellPanel(
                            conditionalPanel("input.Tables=='DESeq Counts'", align='center',
                                h2('Comparison Data tables', align='center'),
                                selectInput("DESeqTable", label='Select table for viewing', choices=NULL),
                                fluidRow(
                                    h3('Significant genes in comparison', align = "center")                                    
                                ),
                                fluidRow(
                                    column(4, 
                                        h4('Total: '), 
                                        textOutput('TotalSig')
                                    ),
                                    column(4, 
                                        h4('Up: '),
                                        textOutput('UpSig')
                                    ),
                                    column(4,
                                        h4('Down: '),
                                        textOutput('DownSig')
                                    )
                                )
                            )
                        )
                    ),
                    column(9,
                        tabsetPanel(id='Tables',
                            tabPanel(title='DESeq Counts', hr(),
                                withSpinner(type=6, color='#5bc0de',
                                    dataTableOutput('deseq_counts')
                                ),
                                fluidRow(
                                    column(12, align='center',downloadButton('DESeqDownload', 'Download the Table'))
                                )
                            )
                        ),
                    )
                )
            ),

        ## Clustering Page
##########################################################################################################################################################
            tabPanel('Clustering Plots',
                fluidRow(
                    column(3,
                        wellPanel(
                                conditionalPanel("input.ClusteringPlots=='PCA' || input.ClusteringPlots=='Distance Matrix'",
                                    selectInput("PCAColor", label='Color By', choices=NULL)
                                ),
                                conditionalPanel("input.ClusteringPlots=='PCA'",
                                    selectInput("PCAShape", label='Shape By', choices=NULL),
                                    radioButtons("Ellipsoid",label='Draw Ellipsoids (requires >3 samples/group)',inline=TRUE, choices=c('No'='FALSE','Yes'='TRUE'), selected='FALSE'),
                                    radioButtons("PCALegendPos", label="Legend Position", inline=TRUE, 
                                        choices=c("Top"='top',"Bottom"='bottom','Right'='right','left'='left'), selected="top"),
                                    sliderInput("PCAPointSize","Point Size: ", min=1, max=30, step=1, value=6),
                                    sliderInput("PCALabelSize","Label Size: ", min=1, max=30, step=1, value=5),
                                    sliderInput("PCALegendIcon","Legend Icon Size", min=1, max=30, step=1, value=8),
                                    sliderInput("PCALegendLabel","Legend Label Size", , min=1, max=30, step=1, value=12)
                                ),
                                conditionalPanel("input.ClusteringPlots=='Distance Matrix'",
                                    sliderInput("DMXsize", 'X-Axis Label Size', min=1, max=30, step=1, value=10),
                                    sliderInput("DMYsize", 'Y-Axis Label Size', min=1, max=30, step=1, value=10),
                                    radioButtons("DMang", label='X-Axis Angle', inline=TRUE,
                                        choices=c('0'=0,'45'=45,'90'=90,'270'=270,'315'=315), selected=45)
                                ),
                                conditionalPanel("input.ClusteringPlots=='Eigencorplots'",
                                    selectizeInput("EigenMetavars", label='Select the metadata variables', choices=NULL, multiple=TRUE),
                                    radioButtons("EFun", label='Correlation Function', inline=TRUE,
                                        choices=c("Pearson"='pearson','Spearman'='spearman','Kendall'='kendall'), selected='pearson'),
                                    radioButtons("ERsquare", label="R or R squared", inline=TRUE, 
                                        choices=c("R"='FALSE',"R Squared"='TRUE'), selected="FALSE"),
                                    radioButtons("ECorrect", label="Multiple Testing Correction", inline=TRUE, 
                                        choices=c('False Discovery Rate'='fdr',"Bonferroni"='bonferroni','None'='none'), selected='fdr'),
                                    radioButtons("ELegendPos", label="Legend Position", inline=TRUE, 
                                        choices=c("Top"='top',"Bottom"='bottom','Right'='right','left'='left'), selected="top"),
                                    sliderInput("ELabSize","Value Label size: ", min=1, max=30, step=1, value=2),
                                    sliderInput("ELabColSize","Legend Label Size: ", min=1, max=30, step=1, value=2),
                                    sliderInput("EXLabSize","X-axis Label Size: ", min=1, max=30, step=1, value=2),
                                    sliderInput("EYLabSize","Y-axis Label Size: ", min=1, max=30, step=1, value=2)                                    
                                ),
                                sliderInput('CHeight', label='Plot Heights: ', min=50, max=2000, step=10, value=800),
                                sliderInput('CWidth', label='Plot Widths: ',  min=50, max=2000, step=10, value=800)
                            )    
                        ),
                    column(9,
                        tabsetPanel(id='ClusteringPlots',
                            tabPanel(title='PCA', hr(), 
                                withSpinner(type=6, color='#5bc0de',
                                        plotOutput("pcaplot", height='100%')
                                ),
                                fluidRow(align='center',style="margin-top:25px;",
                                    uiOutput('PCADownload')
                                )
                            ),
                            tabPanel(title='Distance Matrix',hr(), 
                                withSpinner(type=6, color='#5bc0de',
                                        plotOutput("distplot", height='100%')
                                ),
                                fluidRow(align='center',style="margin-top:25px;",
                                    uiOutput('DMDownload')
                                )
                            ),
                            tabPanel(title='Eigencorplots',hr(), 
                                withSpinner(type=6, color='#5bc0de',  
                                    plotOutput('eigencorplotout', height='100%')
                                ),
                                fluidRow(align='center',style="margin-top:25px;",
                                    uiOutput('EigenDownload')
                                )
                            )
        
                        )
                    )
                )
            ),

        ## Volcano Page
##########################################################################################################################################################
            tabPanel('Volcano Plots',
                fluidRow(
                    column(3,
                        wellPanel(
                                conditionalPanel("input.VolcanoPlots=='Volcanoes'",
                                    selectInput("VolComp", label="Choose Comparison", choices=NULL),
                                    radioButtons("Vpval", label="FDR adjusted P or p-value", inline=TRUE, 
                                        choices=c("FDR adjusted P"='padj',"P-value"='pvalue'), selected="padj"),
                                    radioButtons("VLabelType", label="Select Label Type", inline=TRUE, 
                                        choices=c("Ensembl Gene"='Gene',"Gene Symbol"='Geneid', 'Selected Genes'='Select','None'='None'), selected="None"),
                                    radioButtons("VLegendPos", label="Legend Position", inline=TRUE, 
                                        choices=c("Top"='top',"Bottom"='bottom','Right'='right','left'='left'), selected="top"),
                                    sliderInput("VpointSize","Point size: ", min=1, max=30, step=1, value=5),
                                    sliderInput("VLabSize","Label size: ", min=1, max=30, step=1, value=4),
                                    sliderInput("VLegLabSize","Legend Label size: ", min=1, max=30, step=1, value=10),
                                    sliderInput("VLegIconSize","Legend Icon size: ", min=1, max=30, step=1, value=5),
                                    conditionalPanel("input.VLabelType=='Select'",
                                        selectizeInput("VSelectedGenes", "Please list genes of choice", choices = NULL, multiple=TRUE, options=list(placeholder='Search'))
                                        )                                    
                                ),
                                sliderInput('VHeight', label='Plot Heights: ', min=50, max=2000, step=10, value=800),
                                sliderInput('VWidth', label='Plot Widths: ',  min=50, max=2000, step=10, value=800)
                            )    
                        ),
                    column(9,
                        tabsetPanel(id='VolcanoPlots',
                            tabPanel(title='Volcanoes', hr(), 
                                withSpinner(type=6, color='#5bc0de',  
                                    plotOutput("vol_static", height='100%')
                                ),
                                fluidRow(align='center',style="margin-top:25px;",
                                    uiOutput('VolcanoDownload')
                                )
                            )
                        )
                    )
                )
            ),

        ## Heatmaps
##########################################################################################################################################################
            tabPanel('Heatmap Plots',
                fluidRow(
                    column(3,
                        wellPanel(
                                conditionalPanel("input.HeatPlots=='Heatmap'",
                                    radioButtons("HSubType", label="Select heatmap subsetting", inline=TRUE, 
                                        choices=c("Significantly different"='sigdif','Select your own'='Select'), selected="sigdif"),
                                    conditionalPanel("input.HSubType=='sigdif'",
                                        radioButtons('HSigP', label='Significance based on: ', inline=TRUE,
                                            choices=c('P Value'='pvalue','FDR P'='padj'), selected='padj')
                                        ),
                                    conditionalPanel("input.HSubType=='Select'",
                                        selectizeInput("HSelectedGenes", "Please list genes of choice", choices = NULL, multiple=TRUE, options=list(placeholder='Search'))
                                        ),
                                    radioButtons('HRowClust', label='Cluster rows: ', inline=TRUE,
                                        choices=c('True'='TRUE','False'='FALSE'), selected='TRUE'),
                                    radioButtons('HColClust', label='Cluster columns: ', inline=TRUE,
                                        choices=c('True'='TRUE','False'='FALSE'), selected='TRUE'),
                                    radioButtons('HScale', label='Scale By: ', inline=TRUE,
                                        choices=c('Row'='row','Column'='column','None'='none'), selected='row'),
                                    radioButtons('HShowNames', label='Show the rownames: ', inline=TRUE,
                                        choices=c('No'='FALSE','Yes'='TRUE'), selected='FALSE'),
                                    sliderInput("HXsize", 'X-Axis Label Size', min=1, max=30, step=1, value=10),
                                    radioButtons("Hang", label='X-Axis Angle', inline=TRUE,
                                        choices=c('0'=0,'45'=45,'90'=90,'270'=270,'315'=315), selected=45)
                                ),
                                sliderInput('HHeight', label='Plot Heights: ', min=50, max=2000, step=10, value=800),
                                sliderInput('HWidth', label='Plot Widths: ',  min=50, max=2000, step=10, value=800)
                            )    
                        ),
                    column(9,
                        tabsetPanel(id='HeatPlots',
                            tabPanel(title='Heatmap', hr(), 
                                withSpinner(type=6, color='#5bc0de',  
                                    plotOutput('heatmap_static', height='100%')
                                ),
                                fluidRow(align='center',style="margin-top:25px;",
                                    uiOutput('HeatmapDownload')
                                )
                            )
                        )
                    )
                )
            ),
        
        ## Expression Plots
##########################################################################################################################################################
            tabPanel('Gene Expression Plots',
                fluidRow(
                    column(3,
                        wellPanel(
                                selectizeInput("gene_select", label="Name of gene",
                                    choices = NULL, multiple=FALSE, options=list(placeholder='Search')),
                                selectInput("GEGroup", label='Group by', choices=NULL),
                                sliderInput("GELabelSize", 'Label Size', min=1, max=30, step=1, value=5),
                                sliderInput("GETitleSize", 'Title Size', min=1, max=30, step=1, value=10),
                                radioButtons("GEang", label='X-Axis Angle', inline=TRUE,
                                        choices=c('0'=0,'45'=45,'90'=90,'270'=270,'315'=315), selected=45),
                                radioButtons("GELegendPos", label="Legend Position", inline=TRUE, 
                                        choices=c("Top"='top',"Bottom"='bottom','Right'='right','left'='left'), selected="bottom"),
                                sliderInput('GEHeight', label='Plot Heights: ', min=50, max=2000, step=10, value=800),
                                sliderInput('GEWidth', label='Plot Widths: ',  min=50, max=2000, step=10, value=800)
                            )    
                        ),
                    column(9,
                        tabsetPanel(id='GEPlots',
                            tabPanel('Box Plots', hr(), 
                                withSpinner(type=6, color='#5bc0de',  
                                    plotOutput('box_static', height='100%')
                                ),
                                fluidRow(align='center',style="margin-top:25px;",
                                    uiOutput('BoxPlotDownload')
                                )
                            ),
                            tabPanel('Violin Plots', hr(), 
                                withSpinner(type=6, color='#5bc0de',  
                                    plotOutput('violin_static', height='100%')
                                ),
                                fluidRow(align='center',style="margin-top:25px;",
                                    uiOutput('ViolinPlotDownload')
                                )
                            )
                        )
                    )
                )
            ),

        ## Venn Diagrams
##########################################################################################################################################################
            tabPanel('Venn Diagrams',
                fluidRow(
                    column(3,
                        wellPanel(
                                conditionalPanel("input.VennPlots=='Venn Diagrams'",
                                    column(6, selectInput("VennComp1", label="Choose Comparison1", choices=NULL)),
                                    column(6, colourInput("color1", "Pick color1", "#0000FF",palette = "limited")),
                                    column(6, selectInput("VennComp2", label="Choose Comparison2", choices=NULL)),
                                    column(6, colourInput("color2", "Pick color2", "#00FF00",palette = "limited")),
                                    column(6, selectInput("VennComp3", label="Choose Comparison3", choices=NULL)),
                                    column(6, colourInput("color3", "Pick color3", "#FFFF00",palette = "limited")),
                                    column(6, selectInput("VennComp4", label="Choose Comparison4", choices=NULL)),
                                    column(6, colourInput("color4", "Pick color4", "#7FFFD4",palette = "limited")),
                                    column(6, selectInput("VennComp5", label="Choose Comparison5", choices=NULL)),
                                    column(6, colourInput("color5", "Pick color5", "#FF0000",palette = "limited")),
                                    radioButtons("Vennpval", label="FDR adjusted P or p-value", inline=TRUE, 
                                        choices=c("FDR adjusted P"='padj',"P-value"='pvalue'), selected="padj"),
                                    sliderInput("VennFontSize","Font size: ", min=1, max=30, step=1, value=5),
                                    sliderInput("VennLabelFontSize","Label Font size: ", min=1, max=10, step=0.1, value=1),
                                    sliderInput("VennLineThick","Line thickness: ", min=1, max=30, step=1, value=1),
                                    sliderInput("VennLineType","Line Type: ", min=1, max=10, step=1, value=1),
                                    sliderInput("VennAlpha","Color Opacity: ", min=0, max=1, value=0.5)                                    
                                ),
                                conditionalPanel("input.VennPlots=='UpSet Plots'",
                                    selectizeInput("upset_compares", label="Choose comparisons",
                                        choices = NULL, multiple=TRUE, options=list(placeholder='Search')),
                                    radioButtons("Upsetpval", label="FDR adjusted P or p-value", inline=TRUE, 
                                        choices=c("FDR adjusted P"='padj',"P-value"='pvalue'), selected="padj"),
                                    sliderInput("UpSetintsizetitle","Integer title size: ", min=0.1, max=5, step=0.1, value=2),
                                    sliderInput("UpSetintsizeticklabels","Integer tick size: ", min=0.1, max=5, step=0.1, value=2),
                                    sliderInput("UpSetsetsizetitle","Set title size: ", min=0.1, max=5, step=0.1, value=2),
                                    sliderInput("UpSetsetsizeticklabels","Set tick size: ", min=0.1, max=5, step=0.1, value=2),
                                    sliderInput("UpSetsetnames","Names size: ", min=0.1, max=5, step=0.1, value=3),
                                    sliderInput("UpSetnumbersonbar","Numbers on bar size: ", min=0.1, max=5, step=0.1, value=2),
                                    sliderInput("UpSetPointSize","Point size: ", min=0.1, max=10, step=0.1, value=6), 
                                    sliderInput("UpSetLineSize","Line size: ", min=0.1, max=10, step=0.1, value=4)                           
                                ),
                                sliderInput('VennHeight', label='Plot Heights: ', min=50, max=2000, step=10, value=800),
                                sliderInput('VennWidth', label='Plot Widths: ',  min=50, max=2000, step=10, value=800)  
                            )     
                        ),
                    column(9,
                        tabsetPanel(id='VennPlots',
                            tabPanel(title='Venn Diagrams', hr(), 
                                withSpinner(type=6, color='#5bc0de',  
                                    plotOutput("venn_static", height='100%')
                                ),
                                fluidRow(align='center',style="margin-top:25px;",
                                    uiOutput('VennDownload')
                                )
                            ),
                            tabPanel(title='UpSet Plots', hr(), 
                                withSpinner(type=6, color='#5bc0de',  
                                    plotOutput("upset_static", height='100%')
                                ),
                                fluidRow(align='center',style="margin-top:25px;",
                                    uiOutput('UpSetDownload')
                                )
                            )
                        )
                    )
                )
            ),

        ## Pathways Page
##########################################################################################################################################################
            tabPanel('Enrichment Pathways',
                fluidRow(
                    column(3,
                        wellPanel(
                                conditionalPanel("input.Enrichments=='Gene Set Enrichment' || input.Enrichments=='Over Representation'",
                                    selectInput("EnrichComp", label='Choose Enrichment Comparison', choices=NULL)
                                ),

                                conditionalPanel("input.Enrichments=='Gene Set Enrichment'",
                                    h3('GSEA Settings: '),
                                    radioButtons("GSESet", label="Selected gene set: ", inline=TRUE, 
                                        choices=c("Gene ontology"='GO',"KEGG"='kegg','MSigDB'='msigdb'), selected="GO"),
                                    conditionalPanel("input.GSESet=='GO'",
                                        radioButtons("GSEOnt", label='Select Ontology Source: ', inline=TRUE,
                                        choices=c('All'='ALL','Biological Processes'='BP','Molecular Functions'='MF','Cellular Components'='CC'), selected='ALL')
                                    ),
                                    conditionalPanel("input.GSESet=='msigdb'",
                                        radioButtons("GSEMSig", label='Select Ontology Source: ', inline=TRUE,
                                        choices=c('H'='H','C1'='C1','C2'='C2','C3'='C3','C4'='C4','C5'='C5','C6'='C6','C7'='C7'), selected='H')
                                    ), 
                                    sliderInput("NPerm","Number of Permutations: ", min=500, max=10000, step=500, value=1000),      
                                    radioButtons("GSEPVal", label="P value correction: ", inline=TRUE, 
                                        choices=c("None"='none','Bonferroni'='bonferroni','FDR'='fdr'), selected="none"),
                                    h3('Plot Settings: '),
                                    radioButtons("GSEAPlotType", label="Plot Type: ", inline=TRUE, 
                                        choices=c("Dot Plot"='dotplot','Ridge Chart'='ridgechart','Enrichment Map'='emap'), selected="dotplot"), 
                                    sliderInput("NCategories","Number of pathways to show: ", min=5, max=100, step=5, value=10),
                                    sliderInput("GSEFontSize","Font Size: ", min=1, max=30, step=1, value=10),
                                    conditionalPanel("input.GSEAPlotType=='dotplot'",
                                        sliderInput("GSEYWidth","Pathway Label Width: ", min=5, max=100, step=5, value=40)
                                    ),
                                    radioButtons("ShowGSEATable", label="Show GSEA Table: ", inline=TRUE, 
                                        choices=c("No"='no','Yes'='yes'), selected="no"),                  
                                ),

                                conditionalPanel("input.Enrichments=='Over Representation'",
                                    h3('Over Representation Settings: '),
                                    radioButtons("Enrichpval", label="FDR adjusted P or p-value for input selection", inline=TRUE, 
                                        choices=c("FDR adjusted P"='padj',"P-value"='pvalue'), selected="padj"),
                                    radioButtons("ORASet", label="Selected gene set: ", inline=TRUE, 
                                        choices=c("Gene ontology"='GO',"KEGG"='kegg','MSigDB'='msigdb'), selected="GO"),
                                    conditionalPanel("input.ORASet=='GO'",
                                        radioButtons("ORAOnt", label='Select Ontology Source: ', inline=TRUE,
                                        choices=c('All'='ALL','Biological Processes'='BP','Molecular Functions'='MF','Cellular Components'='CC'), selected='ALL')
                                    ),
                                    conditionalPanel("input.ORASet=='msigdb'",
                                        radioButtons("ORAMSig", label='Select Ontology Source: ', inline=TRUE,
                                        choices=c('H'='H','C1'='C1','C2'='C2','C3'='C3','C4'='C4','C5'='C5','C6'='C6','C7'='C7'), selected='H')
                                    ), 
                                    radioButtons("ORAPVal", label="P value correction: ", inline=TRUE, 
                                        choices=c("None"='none','Bonferroni'='bonferroni','FDR'='fdr'), selected="none"),
                                    h3('Plot Settings: '),
                                    radioButtons("ORAPlotType", label="Plot Type: ", inline=TRUE, 
                                        choices=c("Dot Plot"='dotplot','Bar Chart'='barchart','Enrichment Map'='emap'), selected="dotplot"), 
                                    sliderInput("ORANCategories","Number of pathways to show: ", min=5, max=100, step=5, value=10),
                                    sliderInput("ORAFontSize","Font Size: ", min=1, max=30, step=1, value=10),
                                    conditionalPanel("input.ORAPlotType=='dotplot'",
                                        sliderInput("ORAYWidth","Pathway Label Width: ", min=5, max=100, step=5, value=40)
                                    ),
                                    radioButtons("ShowORATable", label="Show Over Representation Table: ", inline=TRUE, 
                                        choices=c("No"='no','Yes'='yes'), selected="no"),                       
                                ),
                                
                                conditionalPanel("input.Enrichments=='Selected Pathways'",
                                    h4('Data will be based on the last run GSE/ORA analysis'),
                                    radioButtons("SelectPlotType", label='Select type of plot', inline=TRUE,
                                        choices=c('GSEA Plot'='gsea','GSEA Kegg Plot'='gsea_kegg','ORA Kegg Plot'='ora_kegg'), selected='gsea'),
                                    conditionalPanel("input.SelectPlotType=='gsea'",
                                        selectizeInput("GSEAPathwayChoice", label='Choose GSEA Pathway of interest', choices=NULL, 
                                            multiple=F, options=list(placeholder='Enter pathway name'))
                                    ),
                                    conditionalPanel("input.SelectPlotType=='gsea_kegg'",
                                        selectizeInput("GSEAKeggChoice", label='Choose GSEA KEGG Pathway of interest', choices=NULL, 
                                            multiple=F, options=list(placeholder='Enter pathway name')),
                                        actionButton("MakePathviewGSEA","Make your Pathview plot", class='btn btn-info')
                                    ),
                                    conditionalPanel("input.SelectPlotType=='ora_kegg'",
                                        selectizeInput("ORAKeggChoice", label='Choose ORA KEGG Pathway of interest', choices=NULL, 
                                            multiple=F, options=list(placeholder='Enter pathway name')),
                                        actionButton("MakePathviewORA","Make your Pathview plot", class='btn btn-info')
                                    )
                                    
                                ),
                                conditionalPanel("input.Enrichments=='Gene Set Enrichment' || input.Enrichments=='Over Representation'",
                                    sliderInput('EnrichHeight', label='Plot Heights: ', min=50, max=2000, step=10, value=800),
                                    sliderInput('EnrichWidth', label='Plot Widths: ',  min=50, max=2000, step=10, value=800)
                                )
                            )    
                        ),
                    column(9,
                        tabsetPanel(id='Enrichments',
                            tabPanel(title='Gene Set Enrichment',hr(),
                                fluidRow(
                                    withSpinner(type=6, color='#5bc0de',  
                                        plotOutput("gseaplot", height='100%')
                                    )
                                ),
                                fluidRow(align='center',style="margin-top:25px;",
                                    uiOutput('GSEAPlotDownload')
                                ),
                                fluidRow(style="margin-top:25px;",
                                    conditionalPanel("input.ShowGSEATable=='yes'",
                                        h3('Gene Set Enrichment Table'),
                                        withSpinner(type=6, color='#5bc0de',  
                                            dataTableOutput("gsea_table")
                                        )
                                    )
                                ),
                                fluidRow(align='center',style="margin-top:25px;",
                                    uiOutput('GSEATableDownload')
                                )
                            ),

                            tabPanel(title='Over Representation',hr(),
                                fluidRow(
                                    withSpinner(type=6, color='#5bc0de',  
                                        plotOutput("oraplot", height='100%')
                                    )
                                ),
                                fluidRow(align='center',style="margin-top:25px;",
                                    uiOutput('ORAPlotDownload')
                                ),
                                fluidRow(style="margin-top:25px;",
                                    conditionalPanel("input.ShowORATable=='yes'",
                                        h3('Over Representation Analysis Table'),
                                        withSpinner(type=6, color='#5bc0de',  
                                            dataTableOutput("ora_table")
                                        )
                                    )
                                ),
                                fluidRow(align='center',style="margin-top:25px;",
                                    uiOutput('ORATableDownload')
                                )
                            ),

                            tabPanel(title='Selected Pathways',hr(),
                                conditionalPanel("input.SelectPlotType=='gsea'",
                                    fluidRow(
                                        withSpinner(type=6, color='#5bc0de',  
                                            plotOutput('gsea_soloplot', height='100%')
                                        )
                                    ),
                                    fluidRow(align='center',style="margin-top:25px;",
                                        uiOutput('GSEASoloDownload')
                                    ),
                                ),
                                conditionalPanel("input.SelectPlotType=='gsea_kegg'",
                                    fluidRow(
                                        withSpinner(type=6, color='#5bc0de',  
                                            plotOutput('gsea_keggplot', height='100%')
                                        )
                                    ),
                                    fluidRow(align='center',style="margin-top:25px;",
                                        uiOutput('GSEAKEGGDownload')
                                    )
                                ),
                                conditionalPanel("input.SelectPlotType=='ora_kegg'",
                                    fluidRow(
                                        withSpinner(type=6, color='#5bc0de',  
                                            plotOutput('ora_keggplot', height='100%')
                                        )
                                    ),
                                    fluidRow(align='center',style="margin-top:25px;",
                                        uiOutput('ORAKEGGDownload')
                                    )
                                )
                            )                             
                        )
                    )
                )
            ),

        ## Page Template
##########################################################################################################################################################
            tabPanel('Page',
                fluidRow(
                    column(3,
                        wellPanel(
                            h2('test', align='center')
                            # 

                        )
                    ),
                    column(9,
                        tabsetPanel(id='Plot',
                            tabPanel(title='output plot', hr()
                            
                            )
                        )
                    )
                )
            ),
            
        ## Footer
##########################################################################################################################################################
            tags$footer(
                wellPanel(
                    fluidRow(
                        column(4, align='center',
                        tags$a(href="https://github.com/alemenze/magic_shiny_rnaseq_template", icon("github", "fa-3x")),
                        tags$h4('GitHub to submit issues/requests')
                        ),
                        column(4, align='center',
                        tags$a(href="http://www.bioinformagic.io/", icon("magic", "fa-3x")),
                        tags$h4('Bioinfor-MaGIC Home Page')
                        ),
                        column(4, align='center',
                        tags$a(href="https://alemenze.github.io/", icon("address-card", "fa-3x")),
                        tags$h4("Developer's Page")
                        )
                    ),
                    fluidRow(
                        column(12, align='center',
                            HTML('<a href="https://www.youtube.com/watch?v=dQw4w9WgXcQ">
                            <p>&copy; 
                                <script language="javascript" type="text/javascript">
                                var today = new Date()
                                var year = today.getFullYear()
                                document.write(year)
                                </script>
                            </p>
                            </a>
                            ')
                        )
                    ) 
                )
            )
        )#Ends navbarPage,
    )#Ends fluidpage
)#Ends tagList
