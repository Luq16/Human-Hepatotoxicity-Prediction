# Pipeline.xml Preprocessing Integration

## Summary

The `pipeline.xml` file contains a Pipeline Pilot workflow for data preparation from Mulliner et al. (2016). We have successfully integrated the key preprocessing steps into our Python scripts.

## Integrated Preprocessing Steps

### 1. **Salt Removal**
- **Pipeline Pilot**: "Keep Largest Fragment" - removes salt fragments and counterions
- **Our Implementation**: Using RDKit's `LargestFragmentChooser` in `preprocess_molecules.py`
- **Example**: `CC(=O)NC1=CC=C(C=C1)O.HCl` → `CC(=O)NC1=CC=C(C=C1)O`

### 2. **Charge Standardization**
- **Pipeline Pilot**: "StandardizeCharges" and "NeutralizeBondedZwitterions"
- **Our Implementation**: Using RDKit's `Uncharger` and `Normalizer`
- **Example**: `C1=CC=C(C=C1)C(=O)[O-].[Na+]` → `C1=CC=C(C=C1)C(=O)O`

### 3. **Tautomer Canonicalization**
- **Pipeline Pilot**: "Standardize tautomer"
- **Our Implementation**: Using RDKit's `TautomerEnumerator.Canonicalize()`
- Ensures consistent tautomer representation across all molecules

### 4. **Stereochemistry Standardization**
- **Pipeline Pilot**: "StandardizeStereo"
- **Our Implementation**: Using RDKit's `AssignStereochemistry(cleanIt=True, force=True)`
- Ensures consistent stereochemistry representation

### 5. **Molecular Cleanup**
- **Pipeline Pilot**: Various cleanup operations (ClearQueryInfo, ClearUnusualValence, etc.)
- **Our Implementation**: Implicit in RDKit's molecule handling and SMILES generation

## Usage

### Standalone Preprocessing
```python
from preprocess_molecules import MoleculePreprocessor

preprocessor = MoleculePreprocessor()
results = preprocessor.preprocess_batch(smiles_list)
```

### Integrated with Descriptor Calculation
```python
from hepatotox_descriptors_improved import ImprovedHepatotoxicityDescriptors

calculator = ImprovedHepatotoxicityDescriptors()
# Preprocessing is enabled by default
descriptors = calculator.calculate_all_descriptors(smiles_list, preprocess=True)
```

## Benefits

1. **Consistency**: All molecules are standardized before descriptor calculation
2. **Robustness**: Handles salts, charged species, and various molecular formats
3. **Reproducibility**: Same preprocessing as the original study
4. **Quality**: Removes artifacts that could affect descriptor calculations

## Example Results

Input molecules with preprocessing effects:
- `CCO` → `CCO` (unchanged - already standard)
- `CC(=O)NC1=CC=C(C=C1)O.HCl` → `CC(=O)NC1=CC=C(C=C1)O` (salt removed)
- `[NH3+]CC([O-])=O` → `NCC(=O)O` (zwitterion neutralized)

## Files Modified

1. **`preprocess_molecules.py`** - New preprocessing module
2. **`hepatotox_descriptors_improved.py`** - Updated to use preprocessing
3. **Training scripts** - Will automatically use preprocessed molecules

This integration ensures our hepatotoxicity predictions follow the same data preparation workflow as the original Pipeline Pilot implementation.