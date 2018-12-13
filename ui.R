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
RV <- reactiveValues(data = seq_df)

navbarPage(
  "Query 16S rRNA Database",
  tabPanel(
    "Select Sequences",
    chooserInput(
      "mychooser",
      "Available 16S Sequences",
      "Selected 16S Sequences",
      genus_names,
      c(),
      size = 10,
      multiple = TRUE
    ),
    downloadButton("goButton", "Run Analysis"),
    p(
      "Click the button to calculate Multiple Sequence Alignments and phylogeny for selected sequences."
    )
  ),
  tabPanel("Data",
           DT::dataTableOutput("table")),
  tabPanel("Multiple Sequence Alignment",
           verbatimTextOutput("selection"),
           tags$iframe(style="height:600px; width:100%", 
                       src="shiny_placeholder.pdf")
                       #src="shiny_placeholder.pdf")
  ),
  tabPanel("Phylogenetic Tree",
           sidebarLayout(
             sidebarPanel(
               radioButtons("plotType", "Plot Type",
                            c("ggTree" = "p", "MSA" = "l")),
               tags$head(tags$script(src = "message-handler.js")),
               actionButton("pdfButton", "Download PDF")
             ),
             mainPanel(plotOutput("plot"))
           )),
  tabPanel("Summary Metrics",
           "Place holder for summary metrics"
  ),
  tabPanel("Help",
           "Place html help file here"
  )
)
