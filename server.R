library(markdown)
library(DT)
library(Biostrings)
library(ggtree)
library(phangorn)
library(msa)
library(ggplot2)
library(reshape2)

# 
seq_df <-
  read.csv(
    "shiny_app_data.csv",
    header = TRUE,
    sep = ","
  )

# Function to perform multiple alignment
aligneR <- function(selected_df) {
  selected_dna <- DNAStringSet(selected_df$full_seq)
  names(selected_dna) <-paste0(selected_df$genus_species, "_", selected_df$gi_number)
  selected_alignment <- msa(selected_dna, type = "DNA")
  return(selected_alignment)
}

# Function to return distance matrix
dist_matrixR <- function(alignment){
  my_align <- msaConvert(alignment, "phangorn::phyDat")
  dm <- dist.ml(my_align) 
  return(dm)
}

# Function to return tree
tree_makeR <- function(dm) {
  treeUPGMA <- upgma(dm)
  my_tree <- ggtree(treeUPGMA) + geom_tiplab()  + xlim(0, 0.6)
  return(my_tree)
}

RV <- reactiveValues(data = seq_df)

# Interactive elements ###########################################################
function(input, output, session) {
  
  output$selection <- renderPrint({
    input$mychooser[2]
  })
  
  output$table <- DT::renderDataTable({
    DT::datatable(RV$data[seq_df$genus %in% unlist(input$mychooser[2]),])
  })

# Analysis ###########################################################
  

  
  pp <- eventReactive(input$runAll, {
    message("calculating MSA...")
    seq_alignment <- aligneR(RV$data[seq_df$genus %in% unlist(input$mychooser[2]),])
    message("Finished calculating MSA...")
    message("calculating distance matrix...")
    distance_matrix <- dist_matrixR(seq_alignment)
    message("Finished calculating distance matrix...")
    message("calculating tree...")
    phylo_tree <- tree_makeR(distance_matrix)
    message("Finished calculating tree...")
    message(is.ggplot(phylo_tree))
    return(phylo_tree)
  }) 
 
  output$treePlot <- renderPlot({
    pp()
  })
  
  bb <- eventReactive(input$runAll, {
    message("calculating MSA...")
    seq_alignment <- aligneR(RV$data[seq_df$genus %in% unlist(input$mychooser[2]),])
    message("Finished calculating MSA...")
    message("calculating distance matrix...")
    distance_matrix <- dist_matrixR(seq_alignment)
    message("Finished calculating distance matrix...")
    return(distance_matrix)
  }) 
  
  output$dm_heatmap <- renderPlot({
    library(lattice) 
    levelplot(as.matrix(bb()), col.regions = heat.colors(100)[length(heat.colors(100)):1], scales=list(x=list(rot=90)))
  })
  
  
  output$dm_table <- DT::renderDataTable({
    DT::datatable(as.matrix(bb()))
  })
  # output$goButton = downloadHandler(
  #   filename = 'myreport.pdf',
  #   content = function(file) {
  #     msaPrettyPrint(
  #       aligneR(RV$data[seq_df$genus %in% unlist(input$mychooser[2]),])
  #       , file = 'myreport.pdf'
  #       , output="pdf"
  #       , showNames="left"
  #       , showLogo="top"
  #       , consensusColor="BlueRed"
  #       , logoColors="accessible area"
  #       , askForOverwrite=FALSE)
  #     file.rename("myreport.pdf", file) # move pdf to file for downloading
  #   },
  #   contentType = 'application/pdf'
  # )
  
}