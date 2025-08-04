# Feature Selection Options for Hepatotoxicity Prediction

## Overview

Feature selection has been implemented as **optional** steps that can be applied at different stages of the modeling pipeline. This gives you flexibility to choose the level of filtering based on your needs.

## Two Types of Feature Selection

### 1. QSAR Molecular Filters (Pre-descriptor calculation)
**Purpose**: Filter out molecules that are problematic for QSAR modeling  
**Based on**: Pipeline Pilot workflow from the original study  
**When to use**: Before descriptor calculation to remove unsuitable molecules

**Filters applied**:
- Molecular weight limits (100-800 Da)
- LogP limits (< 5.5)
- Polar surface area (< 150 Ų)
- H-bond donors/acceptors limits
- Rotatable bonds (< 12)
- Ring count limits
- Inorganic compounds
- Large hydrocarbons
- Polyphenols

### 2. Statistical Descriptor Filters (Post-descriptor calculation)
**Purpose**: Remove problematic descriptors for machine learning  
**Based on**: Statistical properties of calculated descriptors  
**When to use**: After descriptor calculation, before GA-SVM training

**Filters applied**:
- Constant features (same value for all molecules)
- Low variance features (variance < 0.01)
- Highly correlated features (correlation > 0.95)

## Usage Options

### Basic Usage (No filtering)
```python
from hepatotox_descriptors_improved import ImprovedHepatotoxicityDescriptors

calculator = ImprovedHepatotoxicityDescriptors()
descriptors = calculator.calculate_all_descriptors(smiles_list)
```

### With Molecular Preprocessing Only
```python
descriptors = calculator.calculate_all_descriptors(smiles_list, preprocess=True)
```

### With QSAR Molecular Filters
```python
descriptors = calculator.calculate_all_descriptors(smiles_list, 
                                                  preprocess=True, 
                                                  apply_qsar_filters=True)
```

### With Statistical Descriptor Filters
```python
# First calculate descriptors
descriptors = calculator.calculate_all_descriptors(smiles_list)

# Then apply statistical filters
filtered_descriptors, stats = calculator.filter_descriptors(descriptors, 
                                                           apply_statistical_filters=True)
```

### Full Pipeline (All filters)
```python
# 1. Calculate with molecular filters
descriptors = calculator.calculate_all_descriptors(smiles_list, 
                                                  preprocess=True, 
                                                  apply_qsar_filters=True)

# 2. Apply statistical filters
filtered_descriptors, stats = calculator.filter_descriptors(descriptors, 
                                                           apply_statistical_filters=True)
```

## GA-SVM Integration

The **GA-SVM algorithm already performs feature selection** as its main purpose. Therefore:

- **QSAR filters** are most useful to remove problematic molecules before training
- **Statistical filters** can help by removing obviously problematic descriptors
- **GA selection** will still select the optimal subset from remaining descriptors

## Recommendations

### For Research/Production Use:
1. **Always use**: `preprocess=True` (molecular standardization)
2. **Recommended**: `apply_qsar_filters=True` (removes problematic molecules)
3. **Optional**: Statistical descriptor filters (can improve GA-SVM efficiency)

### For Quick Testing:
1. Use basic calculation without filters
2. Apply only preprocessing if needed

### For Training the M1 Model:
```python
# Recommended approach for M1 model training
descriptors = calculator.calculate_all_descriptors(compounds_smiles, 
                                                  preprocess=True, 
                                                  apply_qsar_filters=True)

# Optional: Apply statistical filters to reduce descriptor space
filtered_descriptors, stats = calculator.filter_descriptors(descriptors, 
                                                           apply_statistical_filters=True)

# Then use with GA-SVM training scripts
```

## Benefits of Optional Design

1. **Flexibility**: Choose filters based on your dataset and goals
2. **Compatibility**: Works with original descriptors if no filtering needed
3. **Performance**: Skip expensive filtering for small datasets
4. **Research**: Compare results with/without different filter combinations
5. **Debugging**: Easier to isolate issues when filters can be disabled

## Error Handling

All feature selection is gracefully handled:
- Missing modules: Falls back to basic calculation with warnings
- Failed filters: Reports issues and continues with available data
- Empty results: Provides meaningful error messages

The system will work even if the `feature_selection.py` module is not available, ensuring backward compatibility.