library(shiny)
require(shinyjs)
library(shinythemes)
require(shinycssloaders)

library(tidyr)
library(tidyverse)
library(dplyr)
library(DT)
library(rhandsontable)
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
            windowTitle = "scRNA",
            titlePanel(
                fluidRow(
                column(2, tags$a(href='http://www.bioinformagic.io/', tags$img(height =75 , src = "MaGIC_Icon_0f344c.svg")), align = 'center'), 
                column(10, fluidRow(
                  column(10, h1(strong('RNA Seq Interactive Visualization Tool'), align = 'center')),
                  column(10, h2(strong('Project'), align = 'center'))))
                ),
                windowTitle = "scRNA" ),
                tags$style(type='text/css', '.navbar{font-size:20px;}'),
                tags$style(type='text/css', '.nav-tabs{padding-bottom:20px;}'),
                tags$head(tags$style(".modal-dialog{ width:1300px}")),

        navbarPage(title ="", id='NAVTABS',

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
