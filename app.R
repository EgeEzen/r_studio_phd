# renv::install("dplyr")
# renv::install("ggplot2")
# renv::install("ggthemr")
# renv::install("readxl")
# renv::install("shiny")
# renv::install("tidyr")
# renv::install("rsconnect")
# renv::install("ggtext")
# renv::restore()

#> shiny::runApp("~/Desktop/PhD/2024 Second Term/shiny app rna seq all/")
#> rsconnect::deployApp("~/Desktop/PhD/2024 Second Term/shiny app rna seq all/",appFileManifest = "manifest.txt")

library(shiny)
library(ggplot2)
library(readxl)
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggtext)
library(envalysis)
library(ggsignif)

# read the data 
# the list goes like nrm_counts_copy,mat_copy,results_vs_df, metadata
healthy_vs_ra_data <- readRDS("healthy_vs_ra_shiny.rds")
tnf_vs_nontnf_data <- readRDS("tnf_vs_nontnf_shinyapp.rds")

gene_list <- rownames(healthy_vs_ra_data[[3]])

# UI and server code follow below

# UI code
ui <- fluidPage(
  titlePanel("Gene Expression - RNA-Sequencing Cultured SFs"),
  tabsetPanel(
    tabPanel(
      title="RA vs Healthy",
      sidebarLayout(
        sidebarPanel(
          selectizeInput("healthy_gene", "Select Gene Name:", choices = NULL, multiple = FALSE, 
                         options = list(placeholder = "Type to search gene...")),
          actionButton("healthy_update", "Draw Plot"),
          tags$hr(),
          div(
            actionButton("toggleSettings", "Settings for Download", icon = icon("cogs")),
            conditionalPanel(
              condition = "input.toggleSettings % 2 == 1",  # Toggle visibility based on button click
              numericInput("healthy_width", "Plot Width (inches)", value = 4, min = 3, max = 20, step = 0.5),
              numericInput("healthy_height", "Plot Height (inches)", value = 4, min = 3, max = 20, step = 0.5)
            )
          ),
          downloadButton("download_healthyPlot", "Download Plot as PNG"),
          tags$hr(),
          downloadButton("download_healthyData", "Download Data as CSV") 
        ),
        mainPanel(
          plotOutput("healthyPlot")
        )
      )
    ),
    
    tabPanel(
      title = "TNF vs Non-TNF Stimulation",
      sidebarLayout(
        sidebarPanel(
          selectizeInput("tnf_gene", "Select Gene Name:", choices = NULL, multiple = FALSE, 
                         options = list(placeholder = "Type to search gene...")),
          actionButton("tnf_update", "Draw Plot"),
          tags$hr(),
          div(
            actionButton("toggleSettings", "Settings for Download", icon = icon("cogs")),
            conditionalPanel(
              condition = "input.toggleSettings % 2 == 1",  # Toggle visibility based on button click
              numericInput("tnf_width", "Plot Width (inches)", value = 4, min = 3, max = 20, step = 0.5),
              numericInput("tnf_height", "Plot Height (inches)", value = 4, min = 3, max = 20, step = 0.5)
            )
          ),
          downloadButton("download_tnfPlot", "Download Violin Plot as PNG"),
          tags$hr(),
          downloadButton("download_tnfData", "Download Data as CSV")
        ),
        mainPanel(
          plotOutput("tnfPlot")
        )
      )
    )
  ),
  tags$hr(), 
  fluidRow(
    column(12,
           tags$footer(
             HTML("<p style='text-align: center; color: #666;'>
          This Shiny app visualizes gene expression data  <br>
        </p>
        <p style='text-align: center; font-size: 0.8em;'>
          2025 Ege Ezen / University of Zurich - Department of Rheumatology / Caroline Ospelt Lab <br>
          Last updated: January 2025.
        </p>"),
             style = "margin-top: 20px; padding: 10px; background-color: #f8f8f8;"
             
           )
    )
  )
)

# Server code
server <- function(input, output, session) {
  
  # --- Healthy Plot ---
  
  # Initialize server-side selectize input
  updateSelectizeInput(session, "healthy_gene", choices = gene_list, server = TRUE)
  
  observeEvent(input$healthy_update, {
    gene_of_interest <- input$healthy_gene
    nrm_counts_copy <- healthy_vs_ra_data[[1]]
    results_vs_df <- healthy_vs_ra_data[[3]]
    
    #prepare df for dotplot
    nrm_counts_filtered <- nrm_counts_copy[rownames(nrm_counts_copy) == gene_of_interest,,drop=FALSE]
    nrm_counts_filtered <- as.data.frame(t(nrm_counts_filtered))
    nrm_counts_filtered <- rownames_to_column(nrm_counts_filtered, var = "Sample") 
    nrm_counts_filtered$Condition <- ifelse(grepl("Healthy",nrm_counts_filtered$Sample), "Healthy", "RA") 
    colnames(nrm_counts_filtered) <- c("Sample","Counts","Condition") 
    checking_condition <- (results_vs_df[rownames(results_vs_df) == gene_of_interest,,drop=FALSE ])$padj
    
    if (nrow(nrm_counts_filtered) == 0) {
      output$healthyPlot <- renderPlot({
        ggplot() + 
          geom_text(aes(x = 1, y = 1, label = "Gene not found!"), size = 6, color = "red") +
          theme_void()
      })
    } else {
      # Render the dot plot 
      
      output$healthyPlot <- renderPlot({
        if (checking_condition <= 0.05){ # with significance
          ggplot(nrm_counts_filtered, aes(x = Condition, y = Counts, fill = Condition)) +
            geom_violin(trim = TRUE, alpha = 0.8) + 
            geom_jitter(width = 0.2, size = 1, alpha = 0.6, color = "black") +
            labs(title = paste0("Gene Expression of\n", gene_of_interest," in Healthy vs RA"), 
                 x = "Condition", y = "Normalized Counts") +
            theme_publish() +
            ylim(NA, max(nrm_counts_filtered$Counts) * 1.2)+
            theme(axis.text.x = element_text(angle = 45, hjust = 1),  
                  legend.position = "none") +
            geom_signif(comparisons = list(c("Healthy", "RA")), # TNF non TNF or Healthy Ra
                               annotations = "*",  # Star for significance
                               y_position = max(nrm_counts_filtered$Counts) * 1.05, 
                               tip_length = 0.02, textsize = 7)
        }
        else {#without significance
          
          ggplot(nrm_counts_filtered, aes(x = Condition, y = Counts, fill = Condition)) +
            geom_violin(trim = TRUE, alpha = 0.8) + 
            geom_jitter(width = 0.2, size = 1, alpha = 0.6, color = "black") +
            labs(title = paste0("Gene Expression of\n", gene_of_interest," in Healthy vs RA"), 
                 x = "Condition", y = "Normalized Counts") +
            theme_publish() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),  
                  legend.position = "none")  
          
        }
        
      })
    }
    
    # Create a download handler for the data
    output$download_healthyData <- downloadHandler(
      filename = function() {
    paste0("dotplot_of_", input$healthy_gene, "_healthy_ra.csv")
      },
      content = function(file) {
        write.csv(nrm_counts_filtered, file, row.names = FALSE)
      }
    )
  })
  
  output$download_healthyPlot <- downloadHandler(
    filename = function() {
      paste0("dotplot_of_",input$healthy_gene,"_healthy_ra.png")
    },
    content = function(file) {
      ggsave(file, plot = last_plot(), device = "png", dpi = 300, bg = "transparent", 
             width = input$healthy_width, height = input$healthy_height)
    }
  )
  
  
  # --- TNF vs non TNF part ---
  
  # Initialize server-side selectize input
  updateSelectizeInput(session, "tnf_gene", choices = gene_list, server = TRUE)
  
  observeEvent(input$tnf_update, {
    gene_of_interest <- input$tnf_gene
    nrm_counts_copy <- tnf_vs_nontnf_data[[1]]
    results_vs_df <- tnf_vs_nontnf_data[[3]]
    
    #prepare df for dotplot
    nrm_counts_filtered <- nrm_counts_copy[rownames(nrm_counts_copy) == gene_of_interest,,drop=FALSE]
    nrm_counts_filtered <- as.data.frame(t(nrm_counts_filtered))
    nrm_counts_filtered <- rownames_to_column(nrm_counts_filtered, var = "Sample") 
    nrm_counts_filtered$Condition <- ifelse(grepl("TNF",nrm_counts_filtered$Sample), "TNF", "non-TNF") 
    colnames(nrm_counts_filtered) <- c("Sample","Counts","Condition") 
    checking_condition <- (results_vs_df[rownames(results_vs_df) == gene_of_interest,,drop=FALSE ])$padj
    
    if (nrow(nrm_counts_filtered) == 0) {
      output$tnfPlot <- renderPlot({
        ggplot() + 
          geom_text(aes(x = 1, y = 1, label = "Gene not found!"), size = 6, color = "red") +
          theme_void()
      })
    } else {
      # Render the dot plot 
      
      output$tnfPlot <- renderPlot({
        if (checking_condition <= 0.05){ # with significance
          ggplot(nrm_counts_filtered, aes(x = Condition, y = Counts, fill = Condition)) +
            geom_violin(trim = TRUE, alpha = 0.8) + 
            geom_jitter(width = 0.2, size = 1, alpha = 0.6, color = "black") +
            labs(title = paste0("Gene Expression of\n", gene_of_interest," in \nTNF vs non-TNF in RA"), 
                 x = "Condition", y = "Normalized Counts") +
            theme_publish() +
            ylim(NA, max(nrm_counts_filtered$Counts) * 1.2)+
            theme(axis.text.x = element_text(angle = 45, hjust = 1),  
                  legend.position = "none") +
            geom_signif(comparisons = list(c("non-TNF", "TNF")), # TNF non TNF or Healthy Ra
                        annotations = "*",  # Star for significance
                        y_position = max(nrm_counts_filtered$Counts) * 1.05, 
                        tip_length = 0.02, textsize = 7)
        }
        else {#without significance
          
          ggplot(nrm_counts_filtered, aes(x = Condition, y = Counts, fill = Condition)) +
            geom_violin(trim = TRUE, alpha = 0.8) + 
            geom_jitter(width = 0.2, size = 1, alpha = 0.6, color = "black") +
            labs(title = paste0("Gene Expression of\n", gene_of_interest," in \nTNF vs non-TNF in RA"), 
                 x = "Condition", y = "Normalized Counts") +
            theme_publish() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),  
                  legend.position = "none")  
          
        }
        
      })
    }
    
    # Create a download handler for the data
    output$download_tnfData <- downloadHandler(
      filename = function() {
        paste0("dotplot_of_", input$tnf_gene, "_tnf_vs_non_tnf.csv")
      },
      content = function(file) {
        write.csv(nrm_counts_filtered, file, row.names = FALSE)
      }
    )
  })

  output$download_tnfPlot <- downloadHandler(
    filename = function() {
      paste0("dotplot_of_",input$tnf_gene,"_tnf_vs_non_tnf.png")
    },
    content = function(file) {
      ggsave(file, plot = last_plot(), device = "png", dpi = 300, bg = "transparent", 
             width = input$tnf_width, height = input$tnf_height)
    }
  )
  
}

# Run the app
shinyApp(ui = ui, server = server)

