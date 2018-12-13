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

aligneR <- function(selected_df) {
  # session$sendCustomMessage(type = 'testmessage', message = type(selected_df))
  dna2 <- DNAStringSet(selected_df$full_seq)
  names(dna2) <-paste0(selected_df$genus_species, "_", selected_df$gi_number)
  selected_alignment <- msa(dna2)
  # my_align <- msaConvert(selected_alignment, "phangorn::phyDat")
  # tmpFile <- tempfile(pattern="msa", tmpdir=".", fileext=".pdf")
  # msa_pdf <- msaPrettyPrint(my_align, file="temp_msa.pdf", output="pdf",
  #                showNames="left", showNumbering="none", showLogo="top",
  #                showConsensus="bottom", logoColors="rasmol",
  #                verbose=FALSE, askForOverwrite=FALSE)
  # dm <- dist.ml(my_align)
  # treeUPGMA <- upgma(dm)
  # my_tree <- ggtree(treeUPGMA) + geom_tiplab() + ggplot2::xlim(0, 0.005)
  return(selected_alignment)
}

tree_makeR <- function(alignment) {
  my_align <- msaConvert(alignment, "phangorn::phyDat")
  dm <- dist.ml(my_align)
  treeUPGMA <- upgma(dm)
  my_tree <- ggtree(treeUPGMA) + geom_tiplab() + ggplot2::xlim(0, 0.005)
  return(my_tree)
}


RV <- reactiveValues(data = seq_df)
library(Biostrings)

function(input, output, session) {
  output$selection <- renderPrint({
    input$mychooser[2]
  })
  output$table <- DT::renderDataTable({
    DT::datatable(RV$data[seq_df$genus %in% unlist(input$mychooser[2]),])
  })
  output$plot <- renderPlot({
    plot(cars, type = input$plotType)
  })
  output$summary <- renderPrint({
    summary(cars)
  })
  
  
  
  # MA <- observeEvent(input$goButton, {
  #  aligneR(RV$data[seq_df$genus %in% unlist(input$mychooser[2]),])
  # })
  
  output$goButton = downloadHandler(
    filename = 'myreport.pdf',
    content = function(file) {
      msaPrettyPrint(
        aligneR(RV$data[seq_df$genus %in% unlist(input$mychooser[2]),])
        , file = 'myreport.pdf'
        , output="pdf"
        , showNames="left"
        , showLogo="top"
        , consensusColor="BlueRed"
        , logoColors="accessible area"
        , askForOverwrite=FALSE)
      file.rename("myreport.pdf", file) # move pdf to file for downloading
    },
    contentType = 'application/pdf'
  )
  
}