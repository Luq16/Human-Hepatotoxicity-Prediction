# Human Hepatotoxicity Prediction Shiny App

This Shiny app implements the M1 human hepatotoxicity prediction model based on Mulliner et al. (2016).

## Requirements

### R Packages
```R
install.packages(c("shiny", "shinythemes", "DT", "reticulate", "e1071", "tidyverse", "plotly"))

# For GA-SVM model:
install.packages(c("genalg", "ggplot2", "reshape", "pROC", "grid", "plyr"))
```

### Python Requirements
```bash
pip install numpy pandas rdkit
```

## Running the App

1. Make sure you have both R and Python installed with the required packages

2. Run the Shiny app:
```R
shiny::runApp("hepatotox_app.R")
```

Or from command line:
```bash
R -e "shiny::runApp('hepatotox_app.R')"
```

## How It Works

1. **Input**: Enter SMILES strings (one per line) in the text area
2. **Descriptor Calculation**: The app uses the Python script `hepatotox_descriptors.py` to calculate molecular descriptors
3. **Prediction**: The app would use the trained M1 GA-SVM model to predict hepatotoxicity
4. **Output**: Results show prediction (Hepatotoxic/Non-hepatotoxic), probability score, and confidence level

## Model Information

The M1 model is a GA-SVM (Genetic Algorithm - Support Vector Machine) model that:
- Uses nu-classification SVM with RBF kernel
- Features selected by genetic algorithm from 100+ molecular descriptors
- Trained on human hepatotoxicity data

### Key Parameters:
- SVM type: nu-classification
- Kernel: Radial basis function
- Nu: 0.7
- Gamma: 0.01

## Note on Model File

The actual trained M1 model file (`M1_hepatotoxicity_model.RData`) is required for real predictions. Without it, the app runs in demonstration mode with mock predictions.

To use the actual model:
1. Train the model using the GA-SVM scripts in `tx5b00465_si_002/`
2. Save the model as `M1_hepatotoxicity_model.RData`
3. Place it in the same directory as the app

## Example Compounds

The app includes example compounds with known hepatotoxicity status:
- Ethanol (Safe at therapeutic doses)
- Acetaminophen (Dose-dependent hepatotoxicity)
- Troglitazone (Withdrawn due to hepatotoxicity)
- Diclofenac (NSAID with potential hepatotoxicity)

## Disclaimer

This tool is for research purposes only. Predictions should not be used as the sole basis for safety decisions. Always consult with toxicology experts and conduct appropriate testing.

## Reference

Mulliner, D., Schmidt, F., Stolte, M., Spirkl, H. P., Czich, A., & Amberg, A. (2016). Computational models for human and animal hepatotoxicity with a global application scope. Chemical research in toxicology, 29(5), 757-767.