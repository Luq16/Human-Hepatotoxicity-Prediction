#!/usr/bin/env python3
"""
Feature Selection for Hepatotoxicity Prediction
Based on Pipeline Pilot QSAR filters and statistical methods

This module implements:
1. QSAR-based molecular filters (from pipeline.xml)
2. Statistical feature selection (variance, correlation)
3. Pre-GA descriptor filtering
"""

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import warnings
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler
from scipy.stats import pearsonr
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class QSARMolecularFilter:
    """
    Implements QSAR-based molecular filtering from the Pipeline Pilot workflow.
    Filters out molecules that may be problematic for QSAR modeling.
    """
    
    def __init__(self):
        """Initialize with default QSAR filter thresholds."""
        # Default thresholds from drug-like/lead-like filters
        self.thresholds = {
            'MW_high': 800,        # Molecular weight upper limit
            'MW_low': 100,         # Molecular weight lower limit
            'RotBonds_high': 12,   # Rotatable bonds upper limit
            'ALogP_high': 5.5,     # LogP upper limit
            'PolarSurfaceArea_high': 150,  # PSA upper limit
            'NHacceptors_high': 12,  # H-bond acceptors upper limit
            'NHdonors_high': 6,     # H-bond donors upper limit
            'NCarbons_low': 3,      # Minimum carbon count
            'NRings_high': 6,       # Maximum ring count
            'NRings_low': 0,        # Minimum ring count
            'NStereoatoms_high': 8  # Maximum stereo centers
        }
    
    def evaluate_molecule(self, smiles):
        """
        Evaluate a single molecule against QSAR filters.
        
        Args:
            smiles (str): SMILES string
            
        Returns:
            dict: Filter results with reasons and values
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'pass': False, 'reasons': ['Invalid_SMILES'], 'values': [smiles]}
        
        reasons = []
        values = []
        
        try:
            # Calculate molecular properties
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            tpsa = Descriptors.TPSA(mol)
            hba = Descriptors.NumHAcceptors(mol)
            hbd = Descriptors.NumHDonors(mol)
            rotbonds = Descriptors.NumRotatableBonds(mol)
            rings = Descriptors.RingCount(mol)
            carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
            stereo_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
            
            # Check inorganic content (non-organic elements)
            organic_elements = {1, 6, 7, 8, 9, 15, 16, 17, 35, 53}  # H, C, N, O, F, P, S, Cl, Br, I
            inorganic_atoms = [atom for atom in mol.GetAtoms() 
                             if atom.GetAtomicNum() not in organic_elements]
            
            # Check if it's a hydrocarbon (only C and H)
            elements = set(atom.GetAtomicNum() for atom in mol.GetAtoms())
            is_hydrocarbon = elements.issubset({1, 6})
            
            # Apply filters
            if mw > self.thresholds['MW_high'] or mw < self.thresholds['MW_low']:
                reasons.append('MW')
                values.append(mw)
            
            if logp > self.thresholds['ALogP_high']:
                reasons.append('ALogP')
                values.append(logp)
            
            if tpsa > self.thresholds['PolarSurfaceArea_high']:
                reasons.append('PSA')
                values.append(tpsa)
            
            if hba > self.thresholds['NHacceptors_high']:
                reasons.append('Acceptors')
                values.append(hba)
            
            if hbd > self.thresholds['NHdonors_high']:
                reasons.append('Donors')
                values.append(hbd)
            
            if rotbonds > self.thresholds['RotBonds_high']:
                reasons.append('RotatableBonds')
                values.append(rotbonds)
            
            if rings > self.thresholds['NRings_high'] or rings < self.thresholds['NRings_low']:
                reasons.append('Rings')
                values.append(rings)
            
            if carbons < self.thresholds['NCarbons_low']:
                reasons.append('NCarbons')
                values.append(carbons)
            
            if stereo_centers > self.thresholds['NStereoatoms_high']:
                reasons.append('Stereoatoms')
                values.append(stereo_centers)
            
            if inorganic_atoms:
                reasons.append('Inorganic')
                values.append('TRUE')
            
            if is_hydrocarbon and carbons > 8:  # Large hydrocarbons
                reasons.append('Hydrocarbon')
                values.append('TRUE')
            
            # Check for polyphenol (multiple phenol groups)
            phenol_pattern = Chem.MolFromSmarts('c1ccccc1O')
            if phenol_pattern:
                phenol_matches = len(mol.GetSubstructMatches(phenol_pattern))
                if phenol_matches > 2:
                    reasons.append('Polyphenol')
                    values.append(phenol_matches)
            
        except Exception as e:
            logger.error(f"Error evaluating molecule {smiles}: {e}")
            reasons.append('Calculation_Error')
            values.append(str(e))
        
        return {
            'pass': len(reasons) == 0,
            'reasons': reasons,
            'values': values,
            'smiles': smiles
        }
    
    def filter_molecules(self, smiles_list):
        """
        Filter a list of molecules.
        
        Args:
            smiles_list (list): List of SMILES strings
            
        Returns:
            dict: Filtering results
        """
        results = []
        for smiles in smiles_list:
            result = self.evaluate_molecule(smiles)
            results.append(result)
        
        passed = [r for r in results if r['pass']]
        failed = [r for r in results if not r['pass']]
        
        return {
            'total': len(smiles_list),
            'passed': len(passed),
            'failed': len(failed),
            'pass_rate': len(passed) / len(smiles_list) if smiles_list else 0,
            'results': results,
            'passed_smiles': [r['smiles'] for r in passed],
            'failed_smiles': [r['smiles'] for r in failed]
        }


class StatisticalFeatureSelector:
    """
    Statistical feature selection methods for descriptor filtering.
    """
    
    def __init__(self, variance_threshold=0.01, correlation_threshold=0.95):
        """
        Initialize with thresholds for feature selection.
        
        Args:
            variance_threshold (float): Minimum variance for features
            correlation_threshold (float): Maximum correlation between features
        """
        self.variance_threshold = variance_threshold
        self.correlation_threshold = correlation_threshold
        self.selected_features = None
        self.feature_stats = {}
    
    def remove_low_variance_features(self, X, feature_names):
        """
        Remove features with low variance.
        
        Args:
            X (array): Feature matrix
            feature_names (list): Feature names
            
        Returns:
            tuple: (filtered_X, filtered_feature_names, removed_features)
        """
        selector = VarianceThreshold(threshold=self.variance_threshold)
        X_filtered = selector.fit_transform(X)
        
        # Get selected feature indices
        selected_idx = selector.get_support()
        filtered_features = [feature_names[i] for i in range(len(feature_names)) if selected_idx[i]]
        removed_features = [feature_names[i] for i in range(len(feature_names)) if not selected_idx[i]]
        
        self.feature_stats['variance_removed'] = removed_features
        self.feature_stats['variance_kept'] = filtered_features
        
        logger.info(f"Removed {len(removed_features)} low-variance features")
        
        return X_filtered, filtered_features, removed_features
    
    def remove_correlated_features(self, X, feature_names):
        """
        Remove highly correlated features.
        
        Args:
            X (array): Feature matrix
            feature_names (list): Feature names
            
        Returns:
            tuple: (filtered_X, filtered_feature_names, removed_features)
        """
        # Calculate correlation matrix
        corr_matrix = np.corrcoef(X.T)
        
        # Find highly correlated pairs
        to_remove = set()
        corr_pairs = []
        
        for i in range(len(feature_names)):
            for j in range(i + 1, len(feature_names)):
                if abs(corr_matrix[i, j]) > self.correlation_threshold:
                    # Remove the feature with lower variance
                    var_i = np.var(X[:, i])
                    var_j = np.var(X[:, j])
                    
                    if var_i < var_j:
                        to_remove.add(i)
                    else:
                        to_remove.add(j)
                    
                    corr_pairs.append({
                        'feature1': feature_names[i],
                        'feature2': feature_names[j],
                        'correlation': corr_matrix[i, j],
                        'removed': feature_names[i] if var_i < var_j else feature_names[j]
                    })
        
        # Remove correlated features
        keep_idx = [i for i in range(len(feature_names)) if i not in to_remove]
        X_filtered = X[:, keep_idx]
        filtered_features = [feature_names[i] for i in keep_idx]
        removed_features = [feature_names[i] for i in to_remove]
        
        self.feature_stats['correlation_pairs'] = corr_pairs
        self.feature_stats['correlation_removed'] = removed_features
        self.feature_stats['correlation_kept'] = filtered_features
        
        logger.info(f"Removed {len(removed_features)} highly correlated features")
        
        return X_filtered, filtered_features, removed_features
    
    def remove_constant_features(self, X, feature_names):
        """
        Remove features that are constant across all samples.
        
        Args:
            X (array): Feature matrix
            feature_names (list): Feature names
            
        Returns:
            tuple: (filtered_X, filtered_feature_names, removed_features)
        """
        # Find constant features
        constant_features = []
        keep_idx = []
        
        for i in range(X.shape[1]):
            if np.all(X[:, i] == X[0, i]):  # All values are the same
                constant_features.append(feature_names[i])
            else:
                keep_idx.append(i)
        
        X_filtered = X[:, keep_idx]
        filtered_features = [feature_names[i] for i in keep_idx]
        
        self.feature_stats['constant_removed'] = constant_features
        self.feature_stats['constant_kept'] = filtered_features
        
        logger.info(f"Removed {len(constant_features)} constant features")
        
        return X_filtered, filtered_features, constant_features
    
    def select_features(self, X, feature_names):
        """
        Apply all statistical feature selection methods.
        
        Args:
            X (array): Feature matrix
            feature_names (list): Feature names
            
        Returns:
            tuple: (filtered_X, filtered_feature_names, selection_stats)
        """
        logger.info(f"Starting feature selection with {X.shape[1]} features")
        
        # 1. Remove constant features
        X, feature_names, _ = self.remove_constant_features(X, feature_names)
        
        # 2. Remove low variance features
        X, feature_names, _ = self.remove_low_variance_features(X, feature_names)
        
        # 3. Remove highly correlated features
        X, feature_names, _ = self.remove_correlated_features(X, feature_names)
        
        self.selected_features = feature_names
        
        logger.info(f"Feature selection complete: {len(feature_names)} features selected")
        
        return X, feature_names, self.feature_stats


def prefilter_descriptors_for_ga(descriptor_df, target_column=None):
    """
    Pre-filter descriptors before GA-SVM optimization.
    
    Args:
        descriptor_df (DataFrame): Descriptor dataframe
        target_column (str): Name of target column (if present)
        
    Returns:
        tuple: (filtered_df, selection_stats)
    """
    logger.info("Pre-filtering descriptors for GA-SVM...")
    
    # Separate features from target and ID columns
    exclude_cols = ['ID', 'SMILES']
    if target_column:
        exclude_cols.append(target_column)
    
    feature_cols = [col for col in descriptor_df.columns if col not in exclude_cols]
    
    # Extract feature matrix
    X = descriptor_df[feature_cols].values
    
    # Handle missing values
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    
    # Apply statistical feature selection
    selector = StatisticalFeatureSelector()
    X_filtered, selected_features, stats = selector.select_features(X, feature_cols)
    
    # Reconstruct dataframe
    filtered_df = descriptor_df[exclude_cols].copy()
    for i, feature in enumerate(selected_features):
        filtered_df[feature] = X_filtered[:, i]
    
    logger.info(f"Descriptor pre-filtering complete: {len(selected_features)}/{len(feature_cols)} features kept")
    
    return filtered_df, stats


def main():
    """Example usage of feature selection."""
    
    # Example molecules for QSAR filtering
    test_molecules = [
        "CCO",  # Simple, should pass
        "CC(=O)NC1=CC=C(C=C1)O",  # Acetaminophen, should pass
        "C" * 50,  # Very large hydrocarbon, should fail
        "[Na+].[Cl-]",  # Inorganic, should fail
        "C1=CC=C(C=C1)O" * 5,  # Polyphenol, should fail
        "CC(C)(C)C(C)(C)C(C)(C)C(C)(C)C(C)(C)C(C)(C)C(C)(C)C(C)(C)C",  # High MW, should fail
    ]
    
    print("QSAR Molecular Filtering Example")
    print("=" * 50)
    
    # Test molecular filtering
    qsar_filter = QSARMolecularFilter()
    filter_results = qsar_filter.filter_molecules(test_molecules)
    
    print(f"Total molecules: {filter_results['total']}")
    print(f"Passed filters: {filter_results['passed']} ({filter_results['pass_rate']:.1%})")
    print(f"Failed filters: {filter_results['failed']}")
    
    print("\nFailed molecules:")
    for result in filter_results['results']:
        if not result['pass']:
            print(f"  {result['smiles'][:50]}... - {', '.join(result['reasons'])}")
    
    # Example descriptor filtering
    print("\nStatistical Feature Selection Example")
    print("=" * 50)
    
    # Create example descriptor matrix
    np.random.seed(42)
    n_samples, n_features = 100, 50
    X = np.random.randn(n_samples, n_features)
    
    # Add some problematic features
    X[:, 10] = 0  # Constant feature
    X[:, 11] = X[:, 0] + 0.01 * np.random.randn(n_samples)  # Highly correlated
    X[:, 12] = 0.001 * np.random.randn(n_samples)  # Low variance
    
    feature_names = [f"Feature_{i}" for i in range(n_features)]
    
    selector = StatisticalFeatureSelector()
    X_filtered, selected_features, stats = selector.select_features(X, feature_names)
    
    print(f"Original features: {len(feature_names)}")
    print(f"Selected features: {len(selected_features)}")
    print(f"Removed constant: {len(stats['constant_removed'])}")
    print(f"Removed low variance: {len(stats['variance_removed'])}")
    print(f"Removed correlated: {len(stats['correlation_removed'])}")


if __name__ == "__main__":
    main()