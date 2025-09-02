#!/usr/bin/env Rscript

#
# Human Hepatotoxicity Prediction Shiny App
# Based on M1 model from Mulliner et al. (2016)
# "Computational Models for Human and Animal Hepatotoxicity with a Global Application Scope"
# Chemical Research in Toxicology
#

library(shiny)
library(shinythemes)
library(DT)
library(reticulate)
library(e1071)
library(tidyverse)
library(plotly)

# Configure Python environment
use_python("/usr/bin/python3", required = FALSE)

# Source the Python descriptor calculator
source_python("hepatotox_descriptors_improved.py")

# Define UI
ui <- fluidPage(
  theme = shinytheme("flatly"),
  
  # Application title
  titlePanel(
    div(
      h2("Human Hepatotoxicity Prediction (M1 Model)", style = "color: #2c3e50;"),
      h4("Based on Mulliner et al. (2016) - Chemical Research in Toxicology", style = "color: #7f8c8d;"),
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
      .descriptor-panel {
        max-height: 400px;
        overflow-y: auto;
      }
      .model-info-box {
        background-color: #f5f5f5;
        padding: 15px;
        border-radius: 5px;
        margin: 10px 0;
      }
      .descriptor-highlight {
        background-color: #fff3cd;
        padding: 2px 4px;
        border-radius: 3px;
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
      
      # Example compounds with known hepatotoxicity
      h5("Example Compounds:"),
      div(
        actionButton("example_safe", "Ethanol (Safe at low doses)", class = "btn-sm btn-success"),
        br(),
        actionButton("example_acetaminophen", "Acetaminophen (Dose-dependent)", class = "btn-sm btn-warning"),
        br(),
        actionButton("example_troglitazone", "Troglitazone (Hepatotoxic)", class = "btn-sm btn-danger"),
        br(),
        actionButton("example_diclofenac", "Diclofenac (NSAID)", class = "btn-sm btn-info")
      ),
      br(),
      
      # Model selection (for future use when multiple models available)
      selectInput("model_version",
                  "Model Version:",
                  choices = c("M1 - Human Hepatotoxicity" = "M1"),
                  selected = "M1"),
      
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
        p(strong("Features:"), " Optimized subset from 100+ descriptors"),
        p(strong("Performance:"), " High accuracy on external validation")
      ),
      
      # Warning note
      div(
        class = "alert alert-warning",
        style = "margin-top: 20px;",
        icon("exclamation-triangle"),
        strong(" Important: "), "This tool is for research purposes only. 
        Predictions should not be used as the sole basis for safety decisions. 
        Always consult with toxicology experts and conduct appropriate testing."
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
        
        tabPanel("Molecular Descriptors",
                 value = "descriptors",
                 br(),
                 h4("Calculated Molecular Descriptors"),
                 p("The M1 model uses an optimized subset of these descriptors selected by genetic algorithm:"),
                 br(),
                 # Descriptor categories
                 fluidRow(
                   column(6,
                     h5("Key Descriptor Categories:"),
                     tags$ul(
                       tags$li("Pharmacophore descriptors (PP, PD, PA, PL)"),
                       tags$li("MACCS fingerprint keys"),
                       tags$li("MOE-style descriptors"),
                       tags$li("ADMET-related descriptors")
                     )
                   ),
                   column(6,
                     h5("Important Features:"),
                     tags$ul(
                       tags$li("Molecular size and complexity"),
                       tags$li("Lipophilicity (logP-related)"),
                       tags$li("H-bond donors/acceptors"),
                       tags$li("Aromatic features")
                     )
                   )
                 ),
                 br(),
                 div(class = "descriptor-panel",
                     DT::dataTableOutput("descriptors_table")
                 )
        ),
        
        tabPanel("Model Details",
                 value = "model_details",
                 br(),
                 h4("GA-SVM Model Implementation"),
                 
                 fluidRow(
                   column(6,
                     h5("Support Vector Machine Parameters:"),
                     div(class = "model-info-box",
                       p(strong("SVM Type:"), " nu-classification"),
                       p(strong("Kernel:"), " Radial basis function (RBF)"),
                       p(strong("Nu:"), " 0.7"),
                       p(strong("Gamma:"), " 0.01"),
                       p(strong("Tolerance:"), " 0.01"),
                       p(strong("Epsilon:"), " 0.1")
                     )
                   ),
                   column(6,
                     h5("Genetic Algorithm Parameters:"),
                     div(class = "model-info-box",
                       p(strong("Population Size:"), " 1000"),
                       p(strong("Mutation Rate:"), " 0.01"),
                       p(strong("Fitness Metric:"), " AUC (Area Under ROC Curve)"),
                       p(strong("Feature Selection:"), " Optimized descriptor subset"),
                       p(strong("Cross-validation:"), " n-fold")
                     )
                   )
                 ),
                 
                 h5("Model Training Process:"),
                 p("The M1 model was developed using:"),
                 tags$ol(
                   tags$li("Large dataset of compounds with known human hepatotoxicity data"),
                   tags$li("Genetic algorithm to select optimal descriptor subset from 100+ molecular descriptors"),
                   tags$li("Support vector machine with RBF kernel for classification"),
                   tags$li("Extensive cross-validation and external validation")
                 ),
                 
                 h5("Applicability Domain:"),
                 p("The model is applicable to drug-like organic compounds. Performance may be reduced for:"),
                 tags$ul(
                   tags$li("Very large molecules (MW > 1000 Da)"),
                   tags$li("Inorganic compounds"),
                   tags$li("Compounds very dissimilar to the training set")
                 )
        ),
        
        tabPanel("About",
                 value = "about",
                 br(),
                 h3("About the M1 Hepatotoxicity Model"),
                 
                 h4("Publication"),
                 div(class = "model-info-box",
                   p(strong("Title:"), " Computational Models for Human and Animal Hepatotoxicity with a Global Application Scope"),
                   p(strong("Authors:"), " Mulliner, D., Schmidt, F., Stolte, M., Spirkl, H. P., Czich, A., & Amberg, A."),
                   p(strong("Journal:"), " Chemical Research in Toxicology"),
                   p(strong("Year:"), " 2016"),
                   p(strong("DOI:"), " 10.1021/acs.chemrestox.5b00465")
                 ),
                 
                 h4("Model Performance"),
                 p("The M1 model showed excellent performance in the original publication:"),
                 tags$ul(
                   tags$li("High accuracy in cross-validation"),
                   tags$li("Good performance on external validation set"),
                   tags$li("Broad applicability across chemical space")
                 ),
                 
                 h4("Implementation Notes"),
                 p("This Shiny app implementation:"),
                 tags$ul(
                   tags$li("Uses the GA-SVM framework from the original publication"),
                   tags$li("Calculates all molecular descriptors using RDKit"),
                   tags$li("Would apply the optimized M1 model for predictions (model file required)"),
                   tags$li("Provides interpretable results with confidence scores")
                 ),
                 
                 h4("Limitations and Disclaimers"),
                 div(class = "alert alert-danger",
                   p(strong("Important:"), " This implementation is for research and educational purposes only."),
                   p("Key limitations:"),
                   tags$ul(
                     tags$li("Predictions are probabilistic and should not be the sole basis for safety decisions"),
                     tags$li("Model performance depends on chemical similarity to training set"),
                     tags$li("Does not account for dose, duration, or individual susceptibility factors"),
                     tags$li("Should be used as part of a comprehensive safety assessment")
                   )
                 ),
                 
                 h4("Contact"),
                 p("For questions about this implementation, please refer to the original publication.")
        )
      )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Initialize the descriptor calculator
  calculator <- NULL
  
  # Try to initialize the Python calculator
  tryCatch({
    calculator <<- ImprovedHepatotoxicityDescriptors()
    showNotification("Descriptor calculator initialized successfully", type = "success", duration = 3)
  }, error = function(e) {
    showNotification(paste("Error initializing descriptor calculator:", e$message), type = "error")
  })
  
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
  
  observeEvent(input$example_diclofenac, {
    updateTextAreaInput(session, "smiles_input", value = "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl")
  })
  
  # Reactive values to store results
  prediction_data <- reactiveValues(
    predictions = NULL,
    descriptors = NULL,
    model_loaded = FALSE
  )
  
  # Function to load M1 model (supports both GA-SVM and quick models)
  load_m1_model <- function() {
    model_file <- "M1_hepatotoxicity_model.RData"
    if (file.exists(model_file)) {
      load(model_file)
      if (exists("model")) {
        return(model)
      }
    }
    return(NULL)
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
      
      incProgress(0.2, detail = "Validating SMILES...")
      
      # Calculate descriptors
      tryCatch({
        if (!is.null(calculator)) {
          incProgress(0.3, detail = "Calculating molecular descriptors...")
          
          # Calculate descriptors using Python
          descriptors_df <- calculator$calculate_all_descriptors(smiles_list)
          
          incProgress(0.6, detail = "Applying M1 model...")
          
          # Try to load the actual M1 model
          m1_model <- load_m1_model()
          
          if (!is.null(m1_model)) {
            # Use actual model for predictions
            showNotification("Using trained M1 model for predictions", type = "success", duration = 3)
            
            # Select only the descriptors used by the model
            available_descriptors <- intersect(m1_model$DESCRIPTORS, names(descriptors_df))
            missing_descriptors <- setdiff(m1_model$DESCRIPTORS, names(descriptors_df))
            
            if (length(missing_descriptors) > 0) {
              showNotification(paste("Warning: Missing descriptors:", paste(missing_descriptors, collapse = ", ")), 
                             type = "warning", duration = 5)
            }
            
            # Prepare model input data
            model_input <- descriptors_df[, available_descriptors, drop = FALSE]
            
            # Handle missing values (replace with 0)
            model_input[is.na(model_input)] <- 0
            
            # Make predictions using the SVM model
            if ("svm_model" %in% names(m1_model)) {
              # Quick model format
              predictions_raw <- predict(m1_model$svm_model, model_input, probability = TRUE)
              prob_scores <- attr(predictions_raw, "probabilities")[, "1"]
            } else {
              # GA-SVM format (assumes model is the SVM object directly)
              predictions_raw <- predict(m1_model, model_input, probability = TRUE)
              prob_scores <- attr(predictions_raw, "probabilities")[, "1"]
            }
            
            predictions <- data.frame(
              SMILES = descriptors_df$SMILES,
              Prediction = ifelse(as.numeric(as.character(predictions_raw)) == 1, "Hepatotoxic", "Non-hepatotoxic"),
              Probability = round(prob_scores, 3),
              Confidence = ifelse(abs(prob_scores - 0.5) > 0.3, "High", 
                                ifelse(abs(prob_scores - 0.5) > 0.1, "Medium", "Low")),
              stringsAsFactors = FALSE
            )
          } else {
            # Mock predictions for demonstration
            # In practice, this would use the actual trained M1 model
            showNotification("Note: Using demonstration mode - actual M1 model file not found", 
                           type = "warning", duration = 5)
            
            # Create mock predictions based on some descriptor values
            # This is just for demonstration and should be replaced with actual model
            mock_toxicity_score <- apply(descriptors_df[, c("SMR_C", "CACO2", "A_HEAVY_C")], 1, function(x) {
              # Simple mock scoring based on molecular properties
              score <- 0.3  # base score
              if (!is.na(x[1]) && x[1] > 50) score <- score + 0.2  # High molecular refractivity
              if (!is.na(x[2]) && x[2] < -1) score <- score + 0.2  # Low permeability
              if (!is.na(x[3]) && x[3] > 30) score <- score + 0.1  # Large molecule
              return(min(score + runif(1, -0.1, 0.1), 0.95))
            })
            
            predictions <- data.frame(
              SMILES = descriptors_df$SMILES,
              Prediction = ifelse(mock_toxicity_score > 0.5, "Hepatotoxic", "Non-hepatotoxic"),
              Probability = round(mock_toxicity_score, 3),
              Confidence = ifelse(abs(mock_toxicity_score - 0.5) > 0.3, "High", 
                                ifelse(abs(mock_toxicity_score - 0.5) > 0.1, "Medium", "Low")),
              stringsAsFactors = FALSE
            )
          }
          
          # Store results
          prediction_data$predictions <- predictions
          prediction_data$descriptors <- descriptors_df
          
          incProgress(0.9, detail = "Finalizing results...")
          
          # Switch to results tab
          updateTabsetPanel(session, "main_tabs", selected = "predictions")
          
        } else {
          showNotification("Descriptor calculator not initialized", type = "error")
        }
        
      }, error = function(e) {
        showNotification(paste("Error during prediction:", e$message), type = "error")
      })
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
      )
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
      icon <- ifelse(row$Prediction == "Hepatotoxic", "WARNING", "OK")
      
      # Color code confidence
      conf_color <- switch(row$Confidence,
                          "High" = "success",
                          "Medium" = "warning",
                          "Low" = "danger")
      
      div(
        class = paste("prediction-box", box_class),
        h4(paste(icon, row$Prediction)),
        p(strong("SMILES: "), tags$code(row$SMILES)),
        p(strong("Probability Score: "), paste0(row$Probability)),
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
        lengthChange = FALSE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      rownames = FALSE,
      extensions = 'Buttons'
    ) %>%
      DT::formatStyle(
        'Prediction',
        backgroundColor = DT::styleEqual(
          c("Hepatotoxic", "Non-hepatotoxic"),
          c("#ffebee", "#e8f5e9")
        )
      ) %>%
      DT::formatStyle(
        'Confidence',
        backgroundColor = DT::styleEqual(
          c("High", "Medium", "Low"),
          c("#d4edda", "#fff3cd", "#f8d7da")
        )
      )
  })
  
  # Render descriptors table
  output$descriptors_table <- DT::renderDataTable({
    req(prediction_data$descriptors)
    
    # Select subset of descriptors to display
    # In practice, you would show only the descriptors used by the M1 model
    important_descriptors <- c("SMILES", "A_HEAVY_C", "SMR_C", "CACO2", 
                              "PP10", "PD1", "PA7", "PL4", "PL5",
                              "MKEY14", "MKEY25", "MKEY56", "MKEY83",
                              "PEOE_VSA5_C", "SLOGP_VSA4_C", "CHIRAL_U_C",
                              "HL2", "A", "LGD7_5", "P_FU4", "DRDRDR")
    
    descriptors_to_show <- prediction_data$descriptors[, intersect(names(prediction_data$descriptors), 
                                                                   important_descriptors)]
    
    # Round numeric columns
    numeric_cols <- sapply(descriptors_to_show, is.numeric)
    descriptors_to_show[numeric_cols] <- round(descriptors_to_show[numeric_cols], 4)
    
    DT::datatable(
      descriptors_to_show,
      options = list(
        pageLength = 10,
        searching = TRUE,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),
      rownames = FALSE,
      extensions = 'Buttons'
    )
  })
}

# Run the application
shinyApp(ui = ui, server = server)