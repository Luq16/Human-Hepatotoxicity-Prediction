#!/usr/bin/env python3
"""
Molecular Preprocessing for Hepatotoxicity Prediction
Based on Pipeline Pilot workflow preprocessing steps from Mulliner et al. (2016)

This module provides standardization and preprocessing functions for molecules
before descriptor calculation.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class MoleculePreprocessor:
    """
    Preprocesses molecules according to the Pipeline Pilot workflow steps.
    """
    
    def __init__(self):
        """Initialize the preprocessor with standard settings."""
        # Initialize salt remover with common salts
        self.salt_remover = SaltRemover.SaltRemover()
        
        # Initialize uncharger for neutralization
        self.uncharger = rdMolStandardize.Uncharger()
        
        # Initialize tautomer canonicalizer
        self.tautomer_canon = rdMolStandardize.TautomerEnumerator()
        
        # Initialize normalizer
        self.normalizer = rdMolStandardize.Normalizer()
        
        # Initialize fragment remover
        self.fragment_remover = rdMolStandardize.LargestFragmentChooser()
        
    def preprocess_molecule(self, smiles):
        """
        Preprocess a single molecule through all standardization steps.
        
        Args:
            smiles (str): Input SMILES string
            
        Returns:
            str: Standardized SMILES string or None if preprocessing fails
        """
        try:
            # Step 1: Parse SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.warning(f"Could not parse SMILES: {smiles}")
                return None
            
            # Step 2: Remove salts/counterions (keep largest fragment)
            logger.debug("Removing salts...")
            mol = self.fragment_remover.choose(mol)
            
            # Step 3: Normalize the molecule
            logger.debug("Normalizing molecule...")
            mol = self.normalizer.normalize(mol)
            
            # Step 4: Neutralize charges where possible
            logger.debug("Neutralizing charges...")
            mol = self.uncharger.uncharge(mol)
            
            # Step 5: Canonicalize tautomer
            logger.debug("Canonicalizing tautomer...")
            mol = self.tautomer_canon.Canonicalize(mol)
            
            # Step 6: Remove explicit hydrogens (for consistent representation)
            mol = Chem.RemoveHs(mol)
            
            # Step 7: Standardize stereochemistry representation
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
            
            # Return standardized SMILES
            return Chem.MolToSmiles(mol, isomericSmiles=True)
            
        except Exception as e:
            logger.error(f"Error preprocessing {smiles}: {str(e)}")
            return None
    
    def preprocess_batch(self, smiles_list, verbose=False):
        """
        Preprocess a batch of molecules.
        
        Args:
            smiles_list (list): List of SMILES strings
            verbose (bool): If True, show progress
            
        Returns:
            list: List of tuples (original_smiles, standardized_smiles, success)
        """
        results = []
        
        for i, smiles in enumerate(smiles_list):
            if verbose and i % 100 == 0:
                logger.info(f"Processing molecule {i+1}/{len(smiles_list)}")
            
            standardized = self.preprocess_molecule(smiles)
            success = standardized is not None
            
            results.append((smiles, standardized, success))
        
        # Summary statistics
        n_success = sum(1 for _, _, success in results if success)
        logger.info(f"Successfully preprocessed {n_success}/{len(smiles_list)} molecules")
        
        return results
    
    def get_preprocessing_stats(self, results):
        """
        Get statistics about preprocessing results.
        
        Args:
            results (list): Output from preprocess_batch
            
        Returns:
            dict: Statistics about preprocessing
        """
        stats = {
            'total': len(results),
            'success': sum(1 for _, _, success in results if success),
            'failed': sum(1 for _, _, success in results if not success),
            'changed': 0,
            'unchanged': 0
        }
        
        # Count how many molecules were changed
        for original, standardized, success in results:
            if success:
                if original != standardized:
                    stats['changed'] += 1
                else:
                    stats['unchanged'] += 1
        
        stats['success_rate'] = stats['success'] / stats['total'] if stats['total'] > 0 else 0
        stats['change_rate'] = stats['changed'] / stats['success'] if stats['success'] > 0 else 0
        
        return stats


def preprocess_for_descriptors(smiles_list):
    """
    Convenience function to preprocess molecules before descriptor calculation.
    
    Args:
        smiles_list (list): List of SMILES strings
        
    Returns:
        list: List of standardized SMILES strings (None for failed molecules)
    """
    preprocessor = MoleculePreprocessor()
    results = preprocessor.preprocess_batch(smiles_list)
    
    # Extract standardized SMILES
    standardized_smiles = []
    for original, standardized, success in results:
        if success:
            standardized_smiles.append(standardized)
        else:
            # Keep original if preprocessing failed
            standardized_smiles.append(original)
            logger.warning(f"Using original SMILES for failed preprocessing: {original}")
    
    return standardized_smiles


def main():
    """Example usage of the preprocessor."""
    
    # Example molecules including salts and charged species
    test_smiles = [
        "CCO",  # Ethanol (simple, no change expected)
        "CC(=O)NC1=CC=C(C=C1)O",  # Acetaminophen
        "CC(=O)NC1=CC=C(C=C1)O.HCl",  # Acetaminophen HCl salt
        "C1=CC=C(C=C1)C(=O)[O-].[Na+]",  # Sodium benzoate
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",  # Caffeine
        "[NH3+]CC([O-])=O",  # Zwitterionic amino acid
        "INVALID_SMILES",  # Invalid (should fail)
    ]
    
    print("Molecular Preprocessing Example")
    print("=" * 50)
    
    # Initialize preprocessor
    preprocessor = MoleculePreprocessor()
    
    # Process molecules
    results = preprocessor.preprocess_batch(test_smiles, verbose=True)
    
    # Display results
    print("\nPreprocessing Results:")
    print("-" * 50)
    for original, standardized, success in results:
        status = "✓" if success else "✗"
        if success:
            changed = " (changed)" if original != standardized else " (unchanged)"
            print(f"{status} {original} → {standardized}{changed}")
        else:
            print(f"{status} {original} → FAILED")
    
    # Show statistics
    stats = preprocessor.get_preprocessing_stats(results)
    print("\nStatistics:")
    print(f"Total molecules: {stats['total']}")
    print(f"Successfully processed: {stats['success']} ({stats['success_rate']:.1%})")
    print(f"Failed: {stats['failed']}")
    print(f"Changed during preprocessing: {stats['changed']} ({stats['change_rate']:.1%})")


if __name__ == "__main__":
    main()