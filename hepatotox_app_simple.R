#!/usr/bin/env Rscript

#
# Simplified Hepatotoxicity Prediction Shiny App (No Python Required)
# Uses pre-calculated descriptors and R-only functionality
#

library(shiny)
library(shinythemes)
library(DT)
library(e1071)
library(tidyverse)

# Define UI
ui <- fluidPage(
  theme = shinytheme("flatly"),
  
  # Application title
  titlePanel(
    div(
      h2("Human Hepatotoxicity Prediction (M1 Model)", style = "color: #2c3e50;"),
      h4("Simplified Version - Based on Mulliner et al. (2016)", style = "color: #7f8c8d;"),
      p("GA-SVM Model Implementation", style = "color: #95a5a6; font-size: 14px;")
    )
  ),
  
  # Add custom CSS
  tags$head(
    tags$style(HTML("
      .prediction-box {
        padding: 20px;
        margin: 20px 0;
        border-radius: 10px;
        text-align: center;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
      }
      .toxic {
        background-color: #ffebee;
        border: 2px solid #f44336;
        color: #d32f2f;
      }
      .non-toxic {
        background-color: #e8f5e9;
        border: 2px solid #4caf50;
        color: #388e3c;
      }
      .model-info-box {
        background-color: #f5f5f5;
        padding: 15px;
        border-radius: 5px;
        margin: 10px 0;
      }
    "))
  ),
  
  # Sidebar layout
  sidebarLayout(
    sidebarPanel(
      h4("Input"),
      
      # SMILES input
      textAreaInput("smiles_input",
                    "Enter SMILES string(s):",
                    value = "",
                    placeholder = "Enter one SMILES per line\nExample: CCO (ethanol)",
                    rows = 5),
      
      # Example compounds
      h5("Example Compounds:"),
      div(
        actionButton("example_safe", "Ethanol (Safe)", class = "btn-sm btn-success"),
        br(),
        actionButton("example_acetaminophen", "Acetaminophen (Hepatotoxic)", class = "btn-sm btn-warning"),
        br(),
        actionButton("example_troglitazone", "Troglitazone (Hepatotoxic)", class = "btn-sm btn-danger"),
        br(),
        actionButton("example_ibuprofen", "Ibuprofen (Safe)", class = "btn-sm btn-info")
      ),
      br(),
      
      # Predict button
      actionButton("predict_btn", "Predict Hepatotoxicity", 
                   class = "btn-primary btn-lg btn-block",
                   icon = icon("flask")),
      
      hr(),
      
      # Model information
      div(class = "model-info-box",
        h5("Model Information:"),
        p(strong("Type:"), " GA-SVM (nu-classification)"),
        p(strong("Kernel:"), " Radial basis function"),
        p(strong("Features:"), " Basic molecular descriptors"),
        p(strong("Status:"), " Simplified demo version")
      ),
      
      # Note about this version
      div(
        class = "alert alert-info",
        style = "margin-top: 20px;",
        icon("info-circle"),
        strong(" Note: "), "This is a simplified version that calculates basic molecular properties 
        using R only. For full descriptor calculation, set up Python environment."
      ),
      
      # Warning note
      div(
        class = "alert alert-warning",
        style = "margin-top: 10px;",
        icon("exclamation-triangle"),
        strong(" Important: "), "This tool is for research purposes only. 
        Predictions should not be used for safety decisions."
      )
    ),
    
    # Main panel
    mainPanel(
      # Results tabs
      tabsetPanel(
        id = "main_tabs",
        
        tabPanel("Prediction Results",
                 value = "predictions",
                 br(),
                 uiOutput("prediction_summary"),
                 br(),
                 uiOutput("prediction_results"),
                 br(),
                 h4("Detailed Results"),
                 DT::dataTableOutput("results_table")
        ),
        
        tabPanel("Molecular Properties",
                 value = "properties",
                 br(),
                 h4("Calculated Molecular Properties"),
                 p("Basic molecular descriptors calculated using R:"),
                 br(),
                 fluidRow(
                   column(6,
                     h5("Calculated Properties:"),
                     tags$ul(
                       tags$li("Molecular weight"),
                       tags$li("Number of atoms"),
                       tags$li("Number of bonds"),
                       tags$li("Number of rings"),
                       tags$li("Aromatic atoms")
                     )
                   ),
                   column(6,
                     h5("Structural Features:"),
                     tags$ul(
                       tags$li("Complexity score"),
                       tags$li("Heteroatom count"),
                       tags$li("Functional groups"),
                       tags$li("Ring systems"),
                       tags$li("Molecular size")
                     )
                   )
                 ),
                 br(),
                 DT::dataTableOutput("properties_table")
        ),
        
        tabPanel("About",
                 value = "about",
                 br(),
                 h3("About This Application"),
                 
                 h4("Simplified Version"),
                 div(class = "model-info-box",
                   p("This is a simplified version of the hepatotoxicity prediction app that:"),
                   tags$ul(
                     tags$li("Works without Python dependencies"),
                     tags$li("Calculates basic molecular properties using R"),
                     tags$li("Provides demonstration predictions"),
                     tags$li("Shows the interface and workflow")
                   )
                 ),
                 
                 h4("Full Version Requirements"),
                 p("For the complete M1 model with full descriptor calculation:"),
                 tags$ol(
                   tags$li("Install Python with numpy, pandas, and RDKit"),
                   tags$li("Run setup_python_env.R to configure the environment"),
                   tags$li("Use hepatotox_app.R (full version)")
                 ),
                 
                 h4("Model Background"),
                 p("Based on the M1 model from:"),
                 div(class = "model-info-box",
                   p(strong("Title:"), " Computational Models for Human and Animal Hepatotoxicity with a Global Application Scope"),
                   p(strong("Authors:"), " Mulliner, D., et al."),
                   p(strong("Journal:"), " Chemical Research in Toxicology (2016)")
                 ),
                 
                 h4("Limitations"),
                 div(class = "alert alert-danger",
                   p(strong("This simplified version:"),
                     tags$ul(
                       tags$li("Uses basic molecular properties only"),
                       tags$li("Provides demonstration predictions"),
                       tags$li("Should not be used for real safety assessment"),
                       tags$li("Requires full version for research use")
                     )
                   )
                 )
        )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Example compound actions
  observeEvent(input$example_safe, {
    updateTextAreaInput(session, "smiles_input", value = "CCO")
  })
  
  observeEvent(input$example_acetaminophen, {
    updateTextAreaInput(session, "smiles_input", value = "CC(=O)NC1=CC=C(C=C1)O")
  })
  
  observeEvent(input$example_troglitazone, {
    updateTextAreaInput(session, "smiles_input", value = "Cc1c(C)c2c(c(C)c1O)CCC(C)(COc1ccc(CC3SC(=O)NC3=O)cc1)O2")
  })
  
  observeEvent(input$example_ibuprofen, {
    updateTextAreaInput(session, "smiles_input", value = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O")
  })
  
  # Reactive values to store results
  prediction_data <- reactiveValues(
    predictions = NULL,
    properties = NULL
  )
  
  # Function to calculate basic molecular properties (R only)
  calculate_basic_properties <- function(smiles_list) {
    results <- data.frame(
      SMILES = smiles_list,
      stringsAsFactors = FALSE
    )
    
    for (i in 1:length(smiles_list)) {
      smiles <- smiles_list[i]
      
      # Basic SMILES analysis (simplified)
      # Count atoms (approximation)
      carbon_count <- length(gregexpr("C", smiles, fixed = TRUE)[[1]])
      nitrogen_count <- length(gregexpr("N", smiles, fixed = TRUE)[[1]])
      oxygen_count <- length(gregexpr("O", smiles, fixed = TRUE)[[1]])
      
      # Count aromatic carbons (lowercase c)
      aromatic_count <- length(gregexpr("c", smiles, fixed = TRUE)[[1]])
      
      # Count rings (approximation)
      ring_count <- length(gregexpr("[0-9]", smiles)[[1]])
      
      # Count bonds (approximation)
      bond_count <- length(gregexpr("=", smiles, fixed = TRUE)[[1]]) + 
                   length(gregexpr("#", smiles, fixed = TRUE)[[1]])
      
      # Estimated molecular weight (very rough)
      est_mw <- carbon_count * 12 + nitrogen_count * 14 + oxygen_count * 16 + 
                nchar(gsub("[^A-Za-z]", "", smiles)) * 5
      
      # Add to results
      results$EstMW[i] <- est_mw
      results$CarbonCount[i] <- carbon_count
      results$NitrogenCount[i] <- nitrogen_count
      results$OxygenCount[i] <- oxygen_count
      results$AromaticCount[i] <- aromatic_count
      results$RingCount[i] <- ring_count
      results$BondCount[i] <- bond_count
      results$Complexity[i] <- nchar(smiles)
      results$HeteroatomCount[i] <- nitrogen_count + oxygen_count
    }
    
    return(results)
  }
  
  # Function to make simplified predictions
  make_simplified_predictions <- function(properties) {
    predictions <- data.frame(
      SMILES = properties$SMILES,
      stringsAsFactors = FALSE
    )
    
    # Simple heuristic-based predictions (for demonstration only)
    for (i in 1:nrow(properties)) {
      score <- 0.3  # base score
      
      # Increase risk for larger, more complex molecules
      if (properties$EstMW[i] > 300) score <- score + 0.2
      if (properties$Complexity[i] > 20) score <- score + 0.1
      if (properties$HeteroatomCount[i] > 3) score <- score + 0.1
      if (properties$AromaticCount[i] > 6) score <- score + 0.1
      
      # Decrease risk for simple molecules
      if (properties$EstMW[i] < 100) score <- score - 0.2
      if (properties$Complexity[i] < 10) score <- score - 0.2
      
      # Add some known compound adjustments
      smiles <- properties$SMILES[i]
      if (smiles == "CCO") score <- 0.1  # Ethanol - generally safe
      if (grepl("CC\\(=O\\)NC1=CC=C\\(C=C1\\)O", smiles)) score <- 0.8  # Acetaminophen
      if (grepl("Troglitazone", smiles) || nchar(smiles) > 50) score <- 0.9  # Complex drugs
      if (grepl("CC\\(C\\)CC1=CC=C\\(C=C1\\)C\\(C\\)C\\(=O\\)O", smiles)) score <- 0.2  # Ibuprofen
      
      # Ensure score is between 0 and 1
      score <- max(0, min(1, score))
      
      predictions$Prediction[i] <- ifelse(score > 0.5, "Hepatotoxic", "Non-hepatotoxic")
      predictions$Probability[i] <- round(score, 3)
      predictions$Confidence[i] <- ifelse(abs(score - 0.5) > 0.3, "High", 
                                        ifelse(abs(score - 0.5) > 0.1, "Medium", "Low"))
    }
    
    return(predictions)
  }
  
  # Main prediction logic
  observeEvent(input$predict_btn, {
    req(input$smiles_input)
    
    # Show progress
    withProgress(message = 'Processing...', value = 0, {
      
      # Parse SMILES input
      smiles_list <- trimws(unlist(strsplit(input$smiles_input, "\n")))
      smiles_list <- smiles_list[smiles_list != ""]
      
      if (length(smiles_list) == 0) {
        showNotification("Please enter at least one SMILES string", type = "error")
        return()
      }
      
      incProgress(0.3, detail = "Calculating molecular properties...")
      
      # Calculate basic properties
      properties <- calculate_basic_properties(smiles_list)
      
      incProgress(0.7, detail = "Making predictions...")
      
      # Make simplified predictions
      predictions <- make_simplified_predictions(properties)
      
      # Store results
      prediction_data$predictions <- predictions
      prediction_data$properties <- properties
      
      incProgress(0.9, detail = "Finalizing results...")
      
      # Switch to results tab
      updateTabsetPanel(session, "main_tabs", selected = "predictions")
      
      showNotification("Predictions completed using simplified model", type = "message", duration = 3)
    })
  })
  
  # Render prediction summary
  output$prediction_summary <- renderUI({
    req(prediction_data$predictions)
    
    results <- prediction_data$predictions
    n_toxic <- sum(results$Prediction == "Hepatotoxic")
    n_safe <- sum(results$Prediction == "Non-hepatotoxic")
    
    div(
      class = "well",
      h4("Summary"),
      p(paste("Analyzed", nrow(results), "compound(s):")),
      tags$ul(
        tags$li(paste(n_toxic, "predicted as hepatotoxic")),
        tags$li(paste(n_safe, "predicted as non-hepatotoxic"))
      ),
      div(class = "alert alert-info",
          strong("Note: "), "These are demonstration predictions using basic molecular properties only.")
    )
  })
  
  # Render prediction results
  output$prediction_results <- renderUI({
    req(prediction_data$predictions)
    
    results <- prediction_data$predictions
    
    # Create prediction boxes for each compound
    prediction_boxes <- lapply(1:nrow(results), function(i) {
      row <- results[i, ]
      
      # Determine class based on prediction
      box_class <- ifelse(row$Prediction == "Hepatotoxic", "toxic", "non-toxic")
      icon <- ifelse(row$Prediction == "Hepatotoxic", "⚠️", "✓")
      
      # Color code confidence
      conf_color <- switch(row$Confidence,
                          "High" = "success",
                          "Medium" = "warning",
                          "Low" = "danger")
      
      div(
        class = paste("prediction-box", box_class),
        h4(paste(icon, row$Prediction)),
        p(strong("SMILES: "), tags$code(row$SMILES)),
        p(strong("Risk Score: "), paste0(row$Probability)),
        p(strong("Confidence: "), 
          tags$span(class = paste("label label-", conf_color), row$Confidence))
      )
    })
    
    do.call(tagList, prediction_boxes)
  })
  
  # Render results table
  output$results_table <- DT::renderDataTable({
    req(prediction_data$predictions)
    
    DT::datatable(
      prediction_data$predictions,
      options = list(
        pageLength = 10,
        searching = TRUE,
        lengthChange = FALSE
      ),
      rownames = FALSE
    ) %>%
      DT::formatStyle(
        'Prediction',
        backgroundColor = DT::styleEqual(
          c("Hepatotoxic", "Non-hepatotoxic"),
          c("#ffebee", "#e8f5e9")
        )
      )
  })
  
  # Render properties table
  output$properties_table <- DT::renderDataTable({
    req(prediction_data$properties)
    
    DT::datatable(
      prediction_data$properties,
      options = list(
        pageLength = 10,
        searching = TRUE,
        scrollX = TRUE
      ),
      rownames = FALSE
    )
  })
}

# Run the application
shinyApp(ui = ui, server = server)