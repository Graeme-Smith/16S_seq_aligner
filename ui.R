source("chooser.R")
library(markdown)
library(DT)
library(Biostrings)
library(ggtree)
library(phangorn)
library(msa)
library(ggplot2)

seq_df <-
  read.csv(
    "shiny_app_data.csv",
    header = TRUE,
    sep = ","
  )
genus_names <- unique(seq_df$genus)

navbarPage(
  "Query 16S rRNA Database",
  tabPanel(
    "Select Sequences",
    titlePanel("Select Sequences"),
    chooserInput(
      "mychooser",
      "Available 16S Sequences",
      "Selected 16S Sequences",
      genus_names,
      c(),
      size = 10,
      multiple = TRUE
    ),
    actionButton("runAll", label="Run Analysis"),
    p("Click the button to calculate Multiple Sequence Alignments and phylogeny for selected sequences.")
  ),
  tabPanel("Data",
           titlePanel("Selected data"),
           DT::dataTableOutput("table")
           ),
  
  tabPanel("Multiple Sequence Alignment",
           verbatimTextOutput("selection"),
           downloadButton('downloadPDF')
           # tags$iframe(style="height:600px; width:100%", 
           #             src="shiny_placeholder.pdf")
  ),
  tabPanel("Distance Matrix",
           titlePanel("Visualise Distance Matrix"),
           mainPanel(plotOutput("dm_heatmap"),
                     DT::dataTableOutput("dm_table"))
  ),
  tabPanel("Phylogenetic Tree",
           titlePanel("Visualise Phylogentic tree"),
           sidebarLayout(
             sidebarPanel(
               actionButton("pdfButton", "Download PDF")
             ),
             mainPanel(plotOutput("treePlot"))
           )),

  tabPanel("Summary Metrics",
           "Place holder for summary metrics"
  ),
  tabPanel("Help",
           "Place html help file here"
  )
)
