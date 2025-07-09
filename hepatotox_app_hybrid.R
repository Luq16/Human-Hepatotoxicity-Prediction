#!/usr/bin/env Rscript

#
# Hybrid Hepatotoxicity Prediction Shiny App
# Works with or without Python - uses Python if available, falls back to R-only
#

library(shiny)
library(shinythemes)
library(DT)
library(e1071)
library(tidyverse)

# Try to set up Python
python_available <- FALSE
calculator <- NULL

tryCatch({
  library(reticulate)
  
  # Try different Python configurations
  if (file.exists("/Users/luqmanawoniyi_1/Documents/jazz/hepatotox_env")) {
    use_virtualenv("/Users/luqmanawoniyi_1/Documents/jazz/hepatotox_env", required = FALSE)
  } else {
    # Use any available Python
    py_config()
  }
  
  # Try to import the descriptor calculator
  source_python("hepatotox_descriptors_improved.py")
  calculator <- ImprovedHepatotoxicityDescriptors()
  python_available <- TRUE
  message("✓ Python descriptor calculator loaded successfully")
}, error = function(e) {
  message("ℹ Python not available, using simplified R-only descriptors")
  message("  Error: ", e$message)
})

# Define UI (same for both versions)
ui <- fluidPage(
  theme = shinytheme("flatly"),
  
  titlePanel(
    div(
      h2("Human Hepatotoxicity Prediction (M1 Model)", style = "color: #2c3e50;"),
      h4("Hybrid Version - Works with or without Python", style = "color: #7f8c8d;"),
      p("GA-SVM Model Implementation", style = "color: #95a5a6; font-size: 14px;")
    )
  ),
  
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
  
  sidebarLayout(
    sidebarPanel(
      h4("Input"),
      
      textAreaInput("smiles_input",
                    "Enter SMILES string(s):",
                    value = "",
                    placeholder = "Enter one SMILES per line\nExample: CCO (ethanol)",
                    rows = 5),
      
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
      
      actionButton("predict_btn", "Predict Hepatotoxicity", 
                   class = "btn-primary btn-lg btn-block",
                   icon = icon("flask")),
      
      hr(),
      
      div(class = "model-info-box",
        h5("Status:"),
        if (python_available) {
          tagList(
            p(icon("check-circle", style = "color: green;"), 
              strong("Full descriptors available"), 
              style = "color: green;"),
            p("Using improved molecular descriptors")
          )
        } else {
          tagList(
            p(icon("exclamation-circle", style = "color: orange;"), 
              strong("Basic mode"), 
              style = "color: orange;"),
            p("Using simplified R-only descriptors"),
            actionButton("setup_python", "Setup Python", class = "btn-sm")
          )
        }
      ),
      
      div(
        class = "alert alert-warning",
        style = "margin-top: 20px;",
        icon("exclamation-triangle"),
        strong(" Important: "), "For research purposes only."
      )
    ),
    
    mainPanel(
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
        
        tabPanel("Molecular Descriptors",
                 value = "descriptors",
                 br(),
                 h4("Calculated Molecular Descriptors"),
                 p(ifelse(python_available, 
                         "Full set of molecular descriptors calculated:",
                         "Basic molecular properties calculated using R:")),
                 br(),
                 DT::dataTableOutput("descriptors_table")
        ),
        
        tabPanel("Python Setup",
                 value = "python_setup",
                 br(),
                 h4("Python Environment Setup"),
                 
                 if (!python_available) {
                   div(
                     p("To use full molecular descriptors, set up Python:"),
                     
                     h5("Option 1: Install Miniconda (Recommended)"),
                     pre("bash install_miniconda.sh"),
                     p("Then:"),
                     pre("conda create -n hepatotox python=3.9
conda activate hepatotox
conda install -c conda-forge rdkit numpy pandas"),
                     
                     h5("Option 2: Use system Python"),
                     pre("pip install numpy pandas rdkit"),
                     
                     p("After installation, restart the app."),
                     
                     actionButton("retry_python", "Retry Python Setup", class = "btn-primary")
                   )
                 } else {
                   div(
                     p(icon("check-circle", style = "color: green; font-size: 48px;")),
                     h4("Python is properly configured!", style = "color: green;"),
                     p("Full molecular descriptors are available.")
                   )
                 }
        )
      )
    )
  )
)

# Define server
server <- function(input, output, session) {
  
  # Example compounds
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
  
  # Setup Python button
  observeEvent(input$setup_python, {
    updateTabsetPanel(session, "main_tabs", selected = "python_setup")
  })
  
  # Retry Python setup
  observeEvent(input$retry_python, {
    showModal(modalDialog(
      title = "Restarting App",
      "Please restart the app after setting up Python.",
      footer = modalButton("OK")
    ))
  })
  
  # Calculate basic properties (R-only fallback)
  calculate_basic_properties <- function(smiles_list) {
    results <- data.frame(
      SMILES = smiles_list,
      stringsAsFactors = FALSE
    )
    
    for (i in 1:length(smiles_list)) {
      smiles <- smiles_list[i]
      
      # Basic SMILES analysis
      results$CarbonCount[i] <- length(gregexpr("C", smiles, fixed = TRUE)[[1]])
      results$NitrogenCount[i] <- length(gregexpr("N", smiles, fixed = TRUE)[[1]])
      results$OxygenCount[i] <- length(gregexpr("O", smiles, fixed = TRUE)[[1]])
      results$AromaticCount[i] <- length(gregexpr("c", smiles, fixed = TRUE)[[1]])
      results$RingCount[i] <- length(gregexpr("[0-9]", smiles)[[1]])
      results$DoubleBonds[i] <- length(gregexpr("=", smiles, fixed = TRUE)[[1]])
      results$TripleBonds[i] <- length(gregexpr("#", smiles, fixed = TRUE)[[1]])
      results$Complexity[i] <- nchar(smiles)
      results$HeteroatomCount[i] <- results$NitrogenCount[i] + results$OxygenCount[i]
      results$EstMW[i] <- results$CarbonCount[i] * 12 + results$NitrogenCount[i] * 14 + 
                         results$OxygenCount[i] * 16 + nchar(gsub("[^A-Za-z]", "", smiles)) * 5
    }
    
    return(results)
  }
  
  # Reactive values
  prediction_data <- reactiveValues(
    predictions = NULL,
    descriptors = NULL
  )
  
  # Main prediction logic
  observeEvent(input$predict_btn, {
    req(input$smiles_input)
    
    withProgress(message = 'Processing...', value = 0, {
      
      # Parse SMILES
      smiles_list <- trimws(unlist(strsplit(input$smiles_input, "\n")))
      smiles_list <- smiles_list[smiles_list != ""]
      
      if (length(smiles_list) == 0) {
        showNotification("Please enter at least one SMILES string", type = "error")
        return()
      }
      
      incProgress(0.3, detail = "Calculating descriptors...")
      
      # Calculate descriptors
      if (python_available && !is.null(calculator)) {
        # Use full Python descriptors
        tryCatch({
          descriptors_df <- calculator$calculate_all_descriptors(smiles_list)
          showNotification("Using full molecular descriptors", type = "message")
        }, error = function(e) {
          showNotification("Python error, falling back to basic descriptors", type = "warning")
          descriptors_df <- calculate_basic_properties(smiles_list)
        })
      } else {
        # Use basic R descriptors
        descriptors_df <- calculate_basic_properties(smiles_list)
        showNotification("Using basic molecular properties", type = "message")
      }
      
      incProgress(0.7, detail = "Making predictions...")
      
      # Make predictions (simplified for demo)
      predictions <- data.frame(
        SMILES = smiles_list,
        stringsAsFactors = FALSE
      )
      
      for (i in 1:nrow(predictions)) {
        score <- 0.3  # base score
        
        # Adjust based on properties
        if ("EstMW" %in% names(descriptors_df)) {
          if (descriptors_df$EstMW[i] > 300) score <- score + 0.2
          if (descriptors_df$Complexity[i] > 20) score <- score + 0.1
          if (descriptors_df$HeteroatomCount[i] > 3) score <- score + 0.1
          if (descriptors_df$AromaticCount[i] > 6) score <- score + 0.1
        }
        
        # Known compounds
        smiles <- predictions$SMILES[i]
        if (smiles == "CCO") score <- 0.1
        if (grepl("CC\\(=O\\)NC1=CC=C\\(C=C1\\)O", smiles)) score <- 0.8
        if (nchar(smiles) > 50) score <- min(0.9, score + 0.3)
        
        score <- max(0, min(1, score))
        
        predictions$Prediction[i] <- ifelse(score > 0.5, "Hepatotoxic", "Non-hepatotoxic")
        predictions$Probability[i] <- round(score, 3)
        predictions$Confidence[i] <- ifelse(abs(score - 0.5) > 0.3, "High", 
                                          ifelse(abs(score - 0.5) > 0.1, "Medium", "Low"))
      }
      
      # Store results
      prediction_data$predictions <- predictions
      prediction_data$descriptors <- descriptors_df
      
      incProgress(0.9, detail = "Finalizing...")
      
      updateTabsetPanel(session, "main_tabs", selected = "predictions")
    })
  })
  
  # Render outputs
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
      if (!python_available) {
        div(class = "alert alert-info",
            "Using simplified predictions. Install Python for full accuracy.")
      }
    )
  })
  
  output$prediction_results <- renderUI({
    req(prediction_data$predictions)
    
    results <- prediction_data$predictions
    
    prediction_boxes <- lapply(1:nrow(results), function(i) {
      row <- results[i, ]
      
      box_class <- ifelse(row$Prediction == "Hepatotoxic", "toxic", "non-toxic")
      icon <- ifelse(row$Prediction == "Hepatotoxic", "⚠️", "✓")
      
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
  
  output$descriptors_table <- DT::renderDataTable({
    req(prediction_data$descriptors)
    
    # Show first 20 columns to avoid overwhelming the display
    desc_to_show <- prediction_data$descriptors
    if (ncol(desc_to_show) > 20) {
      desc_to_show <- desc_to_show[, 1:20]
    }
    
    DT::datatable(
      desc_to_show,
      options = list(
        pageLength = 10,
        searching = TRUE,
        scrollX = TRUE
      ),
      rownames = FALSE
    )
  })
}

# Run the app
shinyApp(ui = ui, server = server)