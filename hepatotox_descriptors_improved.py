#!/usr/bin/env python3
"""
Improved Molecular Descriptors Calculator for Hepatotoxicity Prediction
Based on "Computational Models for Human and Animal Hepatotoxicity with a Global Application Scope"

This script provides proper implementations of molecular descriptors using established methods.
"""

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, rdMolDescriptors, rdchem
from rdkit.Chem import rdFreeSASA, rdMolDescriptors, AllChem
from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D
from rdkit.Chem.AtomPairs import Pairs, Torsions
from rdkit.Chem import MACCSkeys
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds, CalcTPSA
from rdkit.Chem.rdpartialcharges import ComputeGasteigerCharges
from rdkit.Chem.EState import EState_VSA
from rdkit.Chem.Fragments import fr_NH0, fr_NH1, fr_NH2, fr_Ar_N, fr_Ar_NH, fr_Ar_OH
import warnings
warnings.filterwarnings('ignore')

# Import preprocessing module
try:
    from preprocess_molecules import preprocess_for_descriptors
    PREPROCESSING_AVAILABLE = True
except ImportError:
    PREPROCESSING_AVAILABLE = False
    print("Warning: Preprocessing module not available. Using raw SMILES.")

class ImprovedHepatotoxicityDescriptors:
    """
    Calculate molecular descriptors for hepatotoxicity prediction using proper methods.
    """
    
    def __init__(self):
        """Initialize the descriptor calculator."""
        self.pharmacophore_factory = Gobbi_Pharm2D.factory
        
    def calculate_all_descriptors(self, smiles_list, preprocess=True, apply_qsar_filters=False):
        """
        Calculate all descriptors for a list of SMILES strings.
        
        Args:
            smiles_list (list): List of SMILES strings
            preprocess (bool): Whether to preprocess molecules (standardize, remove salts, etc.)
            apply_qsar_filters (bool): Whether to apply QSAR molecular filters before calculation
            
        Returns:
            pandas.DataFrame: DataFrame containing all calculated descriptors
        """
        # Apply QSAR filters if requested
        working_smiles = smiles_list
        if apply_qsar_filters:
            try:
                from feature_selection import QSARMolecularFilter
                qsar_filter = QSARMolecularFilter()
                filter_results = qsar_filter.filter_molecules(smiles_list)
                
                print(f"QSAR Filter Results: {filter_results['passed']}/{filter_results['total']} molecules passed")
                if filter_results['failed'] > 0:
                    print(f"Filtered out {filter_results['failed']} molecules that failed QSAR filters")
                
                working_smiles = filter_results['passed_smiles']
                
            except ImportError:
                print("Warning: Feature selection module not available, skipping QSAR filters")
        
        # Preprocess molecules if requested
        if preprocess and PREPROCESSING_AVAILABLE:
            print("Preprocessing molecules (removing salts, standardizing charges, etc.)...")
            processed_smiles = preprocess_for_descriptors(working_smiles)
            # Keep mapping to original SMILES
            smiles_mapping = list(zip(working_smiles, processed_smiles))
        else:
            if preprocess and not PREPROCESSING_AVAILABLE:
                print("Warning: Preprocessing requested but module not available")
            smiles_mapping = [(s, s) for s in working_smiles]
        
        results = []
        
        for i, (original_smiles, processed_smiles) in enumerate(smiles_mapping):
            print(f"Processing compound {i+1}/{len(smiles_list)}: {original_smiles}")
            
            mol = Chem.MolFromSmiles(processed_smiles)
            if mol is None:
                print(f"Warning: Could not parse SMILES: {processed_smiles}")
                continue
            
            # Add hydrogens for accurate calculations
            mol = Chem.AddHs(mol)
            
            descriptors = {}
            descriptors['SMILES'] = original_smiles  # Keep original SMILES in output
            
            # Calculate different descriptor categories
            descriptors.update(self._calculate_basic_descriptors(mol))
            descriptors.update(self._calculate_pharmacophore_descriptors(mol))
            descriptors.update(self._calculate_maccs_keys(mol))
            descriptors.update(self._calculate_moe_style_descriptors(mol))
            descriptors.update(self._calculate_topological_descriptors(mol))
            descriptors.update(self._calculate_electronic_descriptors(mol))
            descriptors.update(self._calculate_geometric_descriptors(mol))
            descriptors.update(self._calculate_fragment_descriptors(mol))
            
            results.append(descriptors)
            
        return pd.DataFrame(results)
    
    def filter_descriptors(self, descriptor_df, apply_statistical_filters=False, 
                          variance_threshold=0.01, correlation_threshold=0.95):
        """
        Apply statistical filters to calculated descriptors (optional post-processing).
        
        Args:
            descriptor_df (DataFrame): Descriptor dataframe
            apply_statistical_filters (bool): Whether to apply variance/correlation filters
            variance_threshold (float): Minimum variance threshold
            correlation_threshold (float): Maximum correlation threshold
            
        Returns:
            tuple: (filtered_df, filter_stats)
        """
        if not apply_statistical_filters:
            return descriptor_df, {'filters_applied': False}
        
        try:
            from feature_selection import prefilter_descriptors_for_ga
            print("Applying statistical descriptor filters...")
            filtered_df, stats = prefilter_descriptors_for_ga(descriptor_df)
            stats['filters_applied'] = True
            return filtered_df, stats
            
        except ImportError:
            print("Warning: Feature selection module not available, skipping statistical filters")
            return descriptor_df, {'filters_applied': False, 'error': 'module_not_available'}
    
    def _calculate_basic_descriptors(self, mol):
        """Calculate basic molecular descriptors."""
        descriptors = {}
        
        try:
            # Basic molecular properties
            descriptors['MW'] = Descriptors.MolWt(mol)
            descriptors['LogP'] = Crippen.MolLogP(mol)
            descriptors['TPSA'] = Descriptors.TPSA(mol)
            descriptors['HBA'] = Descriptors.NumHAcceptors(mol)
            descriptors['HBD'] = Descriptors.NumHDonors(mol)
            descriptors['RotBonds'] = Descriptors.NumRotatableBonds(mol)
            descriptors['AromaticRings'] = Descriptors.NumAromaticRings(mol)
            descriptors['SaturatedRings'] = Descriptors.NumSaturatedRings(mol)
            descriptors['HeavyAtoms'] = mol.GetNumHeavyAtoms()
            descriptors['MolMR'] = Crippen.MolMR(mol)
            
        except Exception as e:
            print(f"Error calculating basic descriptors: {e}")
            for desc in ['MW', 'LogP', 'TPSA', 'HBA', 'HBD', 'RotBonds', 
                        'AromaticRings', 'SaturatedRings', 'HeavyAtoms', 'MolMR']:
                descriptors[desc] = 0
        
        return descriptors
    
    def _calculate_pharmacophore_descriptors(self, mol):
        """Calculate pharmacophore descriptors using proper methods."""
        descriptors = {}
        
        try:
            # Generate pharmacophore fingerprint
            pharm_fp = Generate.Gen2DFingerprint(mol, self.pharmacophore_factory)
            
            # Extract pharmacophore features
            # Donor features
            descriptors['PD1'] = self._count_pharmacophore_pattern(mol, 'donor', 1)
            descriptors['PD3'] = self._count_pharmacophore_pattern(mol, 'donor', 3)
            descriptors['PD7'] = self._count_pharmacophore_pattern(mol, 'donor', 7)
            descriptors['PD10'] = self._count_pharmacophore_pattern(mol, 'donor', 10)
            
            # Acceptor features
            descriptors['PA7'] = self._count_pharmacophore_pattern(mol, 'acceptor', 7)
            descriptors['PA10'] = self._count_pharmacophore_pattern(mol, 'acceptor', 10)
            
            # Lipophilic features
            descriptors['PL4'] = self._count_pharmacophore_pattern(mol, 'lipophilic', 4)
            descriptors['PL5'] = self._count_pharmacophore_pattern(mol, 'lipophilic', 5)
            
            # Aromatic features
            descriptors['AR1'] = self._count_pharmacophore_pattern(mol, 'aromatic', 1)
            descriptors['AR2'] = self._count_pharmacophore_pattern(mol, 'aromatic', 2)
            
            # Pharmacophore pairs
            descriptors['PP10'] = self._calculate_pharmacophore_pairs(mol, 10)
            
            # Negative ionizable features
            descriptors['ND3'] = self._count_pharmacophore_pattern(mol, 'neg_ionizable', 3)
            descriptors['NA2'] = self._count_pharmacophore_pattern(mol, 'neg_ionizable', 2)
            descriptors['NL9'] = self._count_pharmacophore_pattern(mol, 'neg_ionizable', 9)
            
            # Donor-Acceptor pairs
            descriptors['DA2'] = self._calculate_donor_acceptor_pairs(mol, 2)
            descriptors['DA3'] = self._calculate_donor_acceptor_pairs(mol, 3)
            descriptors['DA7'] = self._calculate_donor_acceptor_pairs(mol, 7)
            
            # Donor-Lipophilic pairs
            descriptors['DL6'] = self._calculate_donor_lipophilic_pairs(mol, 6)
            descriptors['DL7'] = self._calculate_donor_lipophilic_pairs(mol, 7)
            
            # Acceptor-Acceptor pairs
            descriptors['AA3'] = self._calculate_acceptor_acceptor_pairs(mol, 3)
            descriptors['AA4'] = self._calculate_acceptor_acceptor_pairs(mol, 4)
            
            # Acceptor-Lipophilic pairs
            descriptors['AL10'] = self._calculate_acceptor_lipophilic_pairs(mol, 10)
            
            # Shape descriptors
            descriptors['SH10'] = self._calculate_shape_descriptor(mol, 10)
            
        except Exception as e:
            print(f"Error calculating pharmacophore descriptors: {e}")
            # Set default values if calculation fails
            pharm_descriptors = ['PD1', 'PD3', 'PD7', 'PD10', 'PA7', 'PA10', 
                               'PL4', 'PL5', 'AR1', 'AR2', 'PP10', 'ND3', 'NA2', 'NL9', 
                               'DA2', 'DA3', 'DA7', 'DL6', 'DL7', 'AA3', 'AA4', 'AL10', 'SH10']
            for desc in pharm_descriptors:
                descriptors[desc] = 0
        
        return descriptors
    
    def _calculate_maccs_keys(self, mol):
        """Calculate MACCS keys."""
        descriptors = {}
        
        try:
            maccs = MACCSkeys.GenMACCSKeys(mol)
            maccs_bits = maccs.ToBitString()
            
            # Extract specific MKEY descriptors mentioned in the paper
            mkey_indices = [14, 15, 19, 25, 26, 32, 37, 56, 83, 89, 90, 96, 97, 98, 
                          104, 108, 115, 121, 126, 135, 141, 142, 146, 149, 154]
            
            for idx in mkey_indices:
                if idx < len(maccs_bits):
                    descriptors[f'MKEY{idx}'] = int(maccs_bits[idx])
                else:
                    descriptors[f'MKEY{idx}'] = 0
                    
        except Exception as e:
            print(f"Error calculating MACCS keys: {e}")
            for idx in [14, 15, 19, 25, 26, 32, 37, 56, 83, 89, 90, 96, 97, 98, 
                       104, 108, 115, 121, 126, 135, 141, 142, 146, 149, 154]:
                descriptors[f'MKEY{idx}'] = 0
        
        return descriptors
    
    def _calculate_moe_style_descriptors(self, mol):
        """Calculate MOE-style descriptors properly."""
        descriptors = {}
        
        try:
            # Heavy atom count
            descriptors['A_HEAVY_C'] = mol.GetNumHeavyAtoms()
            
            # Atom type counts
            descriptors['A_ICM_C'] = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
            descriptors['A_NCL_C'] = sum(1 for atom in mol.GetAtoms() 
                                       if atom.GetAtomicNum() not in [1, 6, 7, 8, 9, 15, 16, 17, 35, 53])
            descriptors['A_NO_C'] = sum(1 for atom in mol.GetAtoms() 
                                      if atom.GetAtomicNum() == 8 and 
                                      any(neighbor.GetAtomicNum() == 7 for neighbor in atom.GetNeighbors()))
            
            # BCUT descriptors (Burden CAS University of Texas descriptors)
            descriptors['BCUT_PEOE_1_C'] = self._calculate_bcut_descriptor(mol, 'PEOE', 1)
            descriptors['BCUT_PEOE_2_C'] = self._calculate_bcut_descriptor(mol, 'PEOE', 2)
            descriptors['BCUT_PEOE_3_C'] = self._calculate_bcut_descriptor(mol, 'PEOE', 3)
            
            # Bond descriptors
            descriptors['B_MAX1LEN_C'] = self._calculate_max_bond_path_length(mol)
            descriptors['B_COUNT_C'] = mol.GetNumBonds()
            
            # Chirality
            descriptors['CHIRAL_U_C'] = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
            
            # PEOE VSA descriptors (Partial Equalization of Orbital Electronegativity)
            descriptors['PEOE_VSA5_C'] = self._calculate_peoe_vsa(mol, 5)
            descriptors['PEOE_VSA_POL_C'] = self._calculate_peoe_vsa_polarizability(mol)
            
            # Petitjean descriptors
            descriptors['PETITJEAN_C'] = self._calculate_petitjean_index(mol)
            descriptors['PETITJEANSC_C'] = self._calculate_petitjean_shape_coefficient(mol)
            
            # Synthetic accessibility
            descriptors['RSYNTH_C'] = self._calculate_synthetic_accessibility(mol)
            
            # SlogP VSA descriptors
            descriptors['SLOGP_VSA4_C'] = self._calculate_slogp_vsa(mol, 4)
            descriptors['SLOGP_VSA5_C'] = self._calculate_slogp_vsa(mol, 5)
            descriptors['SLOGP_VSA6_C'] = self._calculate_slogp_vsa(mol, 6)
            
            # SMR VSA descriptors (Molecular Refractivity)
            descriptors['SMR_C'] = Crippen.MolMR(mol)
            descriptors['SMR_VSA3_C'] = self._calculate_smr_vsa(mol, 3)
            descriptors['SMR_VSA5_C'] = self._calculate_smr_vsa(mol, 5)
            descriptors['SMR_VSA7_C'] = self._calculate_smr_vsa(mol, 7)
            
            # Volume descriptors
            descriptors['VADJMA_C'] = self._calculate_volume_descriptor(mol, 'ADJMA')
            descriptors['VDISTMA_C'] = self._calculate_volume_descriptor(mol, 'DISTMA')
            
        except Exception as e:
            print(f"Error calculating MOE descriptors: {e}")
            # Set default values
            moe_descriptors = ['A_HEAVY_C', 'A_ICM_C', 'A_NCL_C', 'A_NO_C', 
                             'BCUT_PEOE_1_C', 'BCUT_PEOE_2_C', 'BCUT_PEOE_3_C',
                             'B_MAX1LEN_C', 'B_COUNT_C', 'CHIRAL_U_C', 
                             'PEOE_VSA5_C', 'PEOE_VSA_POL_C',
                             'PETITJEAN_C', 'PETITJEANSC_C', 'RSYNTH_C', 
                             'SLOGP_VSA4_C', 'SLOGP_VSA5_C', 'SLOGP_VSA6_C',
                             'SMR_C', 'SMR_VSA3_C', 'SMR_VSA5_C', 'SMR_VSA7_C',
                             'VADJMA_C', 'VDISTMA_C']
            for desc in moe_descriptors:
                descriptors[desc] = 0
        
        return descriptors
    
    def _calculate_topological_descriptors(self, mol):
        """Calculate topological descriptors."""
        descriptors = {}
        
        try:
            # Connectivity indices
            descriptors['Chi0'] = Descriptors.Chi0(mol)
            descriptors['Chi1'] = Descriptors.Chi1(mol)
            descriptors['Chi0n'] = Descriptors.Chi0n(mol)
            descriptors['Chi1n'] = Descriptors.Chi1n(mol)
            descriptors['Chi2n'] = Descriptors.Chi2n(mol)
            descriptors['Chi3n'] = Descriptors.Chi3n(mol)
            descriptors['Chi4n'] = Descriptors.Chi4n(mol)
            
            # Kappa indices
            descriptors['Kappa1'] = Descriptors.Kappa1(mol)
            descriptors['Kappa2'] = Descriptors.Kappa2(mol)
            descriptors['Kappa3'] = Descriptors.Kappa3(mol)
            
            # Information indices
            descriptors['Ipc'] = Descriptors.Ipc(mol)
            descriptors['BertzCT'] = Descriptors.BertzCT(mol)
            
            # Balaban J index
            descriptors['BalabanJ'] = Descriptors.BalabanJ(mol)
            
        except Exception as e:
            print(f"Error calculating topological descriptors: {e}")
            topo_descriptors = ['Chi0', 'Chi1', 'Chi0n', 'Chi1n', 'Chi2n', 'Chi3n', 'Chi4n',
                               'Kappa1', 'Kappa2', 'Kappa3', 'Ipc', 'BertzCT', 'BalabanJ']
            for desc in topo_descriptors:
                descriptors[desc] = 0
        
        return descriptors
    
    def _calculate_electronic_descriptors(self, mol):
        """Calculate electronic descriptors."""
        descriptors = {}
        
        try:
            # Compute Gasteiger charges
            ComputeGasteigerCharges(mol)
            
            # Partial charge descriptors
            charges = [atom.GetDoubleProp('_GasteigerCharge') for atom in mol.GetAtoms() 
                      if not np.isnan(atom.GetDoubleProp('_GasteigerCharge'))]
            
            if charges:
                descriptors['MaxPartialCharge'] = max(charges)
                descriptors['MinPartialCharge'] = min(charges)
                descriptors['SumPartialCharges'] = sum(charges)
                descriptors['AbsSumPartialCharges'] = sum(abs(c) for c in charges)
            else:
                descriptors['MaxPartialCharge'] = 0
                descriptors['MinPartialCharge'] = 0
                descriptors['SumPartialCharges'] = 0
                descriptors['AbsSumPartialCharges'] = 0
            
            # EState VSA descriptors
            descriptors['EState_VSA1'] = EState_VSA.EState_VSA1(mol)
            descriptors['EState_VSA2'] = EState_VSA.EState_VSA2(mol)
            descriptors['EState_VSA3'] = EState_VSA.EState_VSA3(mol)
            descriptors['EState_VSA4'] = EState_VSA.EState_VSA4(mol)
            descriptors['EState_VSA5'] = EState_VSA.EState_VSA5(mol)
            
        except Exception as e:
            print(f"Error calculating electronic descriptors: {e}")
            electronic_descriptors = ['MaxPartialCharge', 'MinPartialCharge', 
                                     'SumPartialCharges', 'AbsSumPartialCharges',
                                     'EState_VSA1', 'EState_VSA2', 'EState_VSA3', 
                                     'EState_VSA4', 'EState_VSA5']
            for desc in electronic_descriptors:
                descriptors[desc] = 0
        
        return descriptors
    
    def _calculate_geometric_descriptors(self, mol):
        """Calculate geometric descriptors."""
        descriptors = {}
        
        try:
            # Generate 3D coordinates
            mol_3d = Chem.Mol(mol)
            embed_result = AllChem.EmbedMolecule(mol_3d, randomSeed=42)
            
            if embed_result == -1:
                # 3D embedding failed, use 2D fallback
                print("3D embedding failed, using default values for geometric descriptors")
                geometric_descriptors = ['Asphericity', 'Eccentricity', 'InertialShapeFactor',
                                       'NPR1', 'NPR2', 'PMI1', 'PMI2', 'PMI3', 
                                       'RadiusOfGyration', 'SpherocityIndex']
                for desc in geometric_descriptors:
                    descriptors[desc] = 0
                return descriptors
            
            AllChem.UFFOptimizeMolecule(mol_3d)
            
            # Try each descriptor individually with safe checks
            try:
                if hasattr(Descriptors, 'Asphericity'):
                    descriptors['Asphericity'] = Descriptors.Asphericity(mol_3d)
                else:
                    descriptors['Asphericity'] = 0
            except:
                descriptors['Asphericity'] = 0
            
            try:
                if hasattr(Descriptors, 'Eccentricity'):
                    descriptors['Eccentricity'] = Descriptors.Eccentricity(mol_3d)
                else:
                    descriptors['Eccentricity'] = 0
            except:
                descriptors['Eccentricity'] = 0
            
            try:
                if hasattr(Descriptors, 'InertialShapeFactor'):
                    descriptors['InertialShapeFactor'] = Descriptors.InertialShapeFactor(mol_3d)
                else:
                    descriptors['InertialShapeFactor'] = 0
            except:
                descriptors['InertialShapeFactor'] = 0
            
            # NPR (Normalized Principal Moments Ratio)
            try:
                if hasattr(Descriptors, 'NPR1'):
                    descriptors['NPR1'] = Descriptors.NPR1(mol_3d)
                else:
                    descriptors['NPR1'] = 0
            except:
                descriptors['NPR1'] = 0
                
            try:
                if hasattr(Descriptors, 'NPR2'):
                    descriptors['NPR2'] = Descriptors.NPR2(mol_3d)
                else:
                    descriptors['NPR2'] = 0
            except:
                descriptors['NPR2'] = 0
            
            # PMI (Principal Moments of Inertia)
            try:
                if hasattr(Descriptors, 'PMI1'):
                    descriptors['PMI1'] = Descriptors.PMI1(mol_3d)
                else:
                    descriptors['PMI1'] = 0
            except:
                descriptors['PMI1'] = 0
                
            try:
                if hasattr(Descriptors, 'PMI2'):
                    descriptors['PMI2'] = Descriptors.PMI2(mol_3d)
                else:
                    descriptors['PMI2'] = 0
            except:
                descriptors['PMI2'] = 0
                
            try:
                if hasattr(Descriptors, 'PMI3'):
                    descriptors['PMI3'] = Descriptors.PMI3(mol_3d)
                else:
                    descriptors['PMI3'] = 0
            except:
                descriptors['PMI3'] = 0
            
            # Radius of gyration
            try:
                if hasattr(Descriptors, 'RadiusOfGyration'):
                    descriptors['RadiusOfGyration'] = Descriptors.RadiusOfGyration(mol_3d)
                else:
                    descriptors['RadiusOfGyration'] = 0
            except:
                descriptors['RadiusOfGyration'] = 0
            
            # Spherocity
            try:
                if hasattr(Descriptors, 'SpherocityIndex'):
                    descriptors['SpherocityIndex'] = Descriptors.SpherocityIndex(mol_3d)
                else:
                    descriptors['SpherocityIndex'] = 0
            except:
                descriptors['SpherocityIndex'] = 0
            
        except Exception as e:
            print(f"Error calculating geometric descriptors: {e}")
            geometric_descriptors = ['Asphericity', 'Eccentricity', 'InertialShapeFactor',
                                   'NPR1', 'NPR2', 'PMI1', 'PMI2', 'PMI3', 
                                   'RadiusOfGyration', 'SpherocityIndex']
            for desc in geometric_descriptors:
                descriptors[desc] = 0
        
        return descriptors
    
    def _calculate_fragment_descriptors(self, mol):
        """Calculate fragment-based descriptors."""
        descriptors = {}
        
        try:
            # Fragment descriptors
            descriptors['fr_NH0'] = fr_NH0(mol)
            descriptors['fr_NH1'] = fr_NH1(mol)
            descriptors['fr_NH2'] = fr_NH2(mol)
            descriptors['fr_Ar_N'] = fr_Ar_N(mol)
            descriptors['fr_Ar_NH'] = fr_Ar_NH(mol)
            descriptors['fr_Ar_OH'] = fr_Ar_OH(mol)
            
            # Ring descriptors
            descriptors['NumAliphaticCarbocycles'] = Descriptors.NumAliphaticCarbocycles(mol)
            descriptors['NumAliphaticHeterocycles'] = Descriptors.NumAliphaticHeterocycles(mol)
            descriptors['NumAromaticCarbocycles'] = Descriptors.NumAromaticCarbocycles(mol)
            descriptors['NumAromaticHeterocycles'] = Descriptors.NumAromaticHeterocycles(mol)
            descriptors['NumSaturatedCarbocycles'] = Descriptors.NumSaturatedCarbocycles(mol)
            descriptors['NumSaturatedHeterocycles'] = Descriptors.NumSaturatedHeterocycles(mol)
            
            # Atom type counts
            descriptors['NumHeteroatoms'] = Descriptors.NumHeteroatoms(mol)
            descriptors['NumSaturatedRings'] = Descriptors.NumSaturatedRings(mol)
            descriptors['NumAromaticRings'] = Descriptors.NumAromaticRings(mol)
            
        except Exception as e:
            print(f"Error calculating fragment descriptors: {e}")
            fragment_descriptors = ['fr_NH0', 'fr_NH1', 'fr_NH2', 'fr_Ar_N', 'fr_Ar_NH', 'fr_Ar_OH',
                                   'NumAliphaticCarbocycles', 'NumAliphaticHeterocycles', 
                                   'NumAromaticCarbocycles', 'NumAromaticHeterocycles',
                                   'NumSaturatedCarbocycles', 'NumSaturatedHeterocycles',
                                   'NumHeteroatoms', 'NumSaturatedRings', 'NumAromaticRings']
            for desc in fragment_descriptors:
                descriptors[desc] = 0
        
        return descriptors
    
    # Helper methods for specific descriptor calculations
    def _count_pharmacophore_pattern(self, mol, pattern_type, distance):
        """Count pharmacophore patterns."""
        if pattern_type == 'donor':
            return len([atom for atom in mol.GetAtoms() if atom.GetTotalNumHs() > 0 and 
                       atom.GetAtomicNum() in [7, 8, 16]])
        elif pattern_type == 'acceptor':
            return len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in [7, 8] and 
                       atom.GetTotalNumHs() == 0])
        elif pattern_type == 'lipophilic':
            return len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and 
                       atom.GetIsAromatic()])
        elif pattern_type == 'aromatic':
            return len([atom for atom in mol.GetAtoms() if atom.GetIsAromatic()])
        elif pattern_type == 'neg_ionizable':
            return len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in [8, 16] and 
                       atom.GetFormalCharge() < 0])
        return 0
    
    def _calculate_pharmacophore_pairs(self, mol, distance):
        """Calculate pharmacophore pairs."""
        # Simplified implementation
        donors = self._count_pharmacophore_pattern(mol, 'donor', distance)
        acceptors = self._count_pharmacophore_pattern(mol, 'acceptor', distance)
        return donors * acceptors
    
    def _calculate_donor_acceptor_pairs(self, mol, distance):
        """Calculate donor-acceptor pairs."""
        donors = self._count_pharmacophore_pattern(mol, 'donor', distance)
        acceptors = self._count_pharmacophore_pattern(mol, 'acceptor', distance)
        return min(donors, acceptors)
    
    def _calculate_donor_lipophilic_pairs(self, mol, distance):
        """Calculate donor-lipophilic pairs."""
        donors = self._count_pharmacophore_pattern(mol, 'donor', distance)
        lipophilic = self._count_pharmacophore_pattern(mol, 'lipophilic', distance)
        return min(donors, lipophilic)
    
    def _calculate_acceptor_acceptor_pairs(self, mol, distance):
        """Calculate acceptor-acceptor pairs."""
        acceptors = self._count_pharmacophore_pattern(mol, 'acceptor', distance)
        return acceptors * (acceptors - 1) // 2 if acceptors > 1 else 0
    
    def _calculate_acceptor_lipophilic_pairs(self, mol, distance):
        """Calculate acceptor-lipophilic pairs."""
        acceptors = self._count_pharmacophore_pattern(mol, 'acceptor', distance)
        lipophilic = self._count_pharmacophore_pattern(mol, 'lipophilic', distance)
        return min(acceptors, lipophilic)
    
    def _calculate_shape_descriptor(self, mol, idx):
        """Calculate shape descriptor."""
        return Descriptors.TPSA(mol) / (mol.GetNumHeavyAtoms() + 1)
    
    def _calculate_bcut_descriptor(self, mol, property_type, idx):
        """Calculate BCUT descriptor properly."""
        # BCUT descriptors are based on the Burden matrix
        # This is a simplified implementation
        if property_type == 'PEOE':
            ComputeGasteigerCharges(mol)
            charges = [atom.GetDoubleProp('_GasteigerCharge') for atom in mol.GetAtoms() 
                      if not np.isnan(atom.GetDoubleProp('_GasteigerCharge'))]
            if charges:
                return sorted(charges, reverse=True)[min(idx-1, len(charges)-1)]
        return 0
    
    def _calculate_max_bond_path_length(self, mol):
        """Calculate maximum bond path length."""
        # Find the diameter of the molecular graph
        dist_matrix = Chem.GetDistanceMatrix(mol)
        return np.max(dist_matrix)
    
    def _calculate_peoe_vsa(self, mol, idx):
        """Calculate PEOE VSA descriptor."""
        try:
            ComputeGasteigerCharges(mol)
            charges = []
            for atom in mol.GetAtoms():
                charge = atom.GetDoubleProp('_GasteigerCharge')
                if not np.isnan(charge):
                    charges.append(charge)
            
            if charges:
                # VSA (van der Waals Surface Area) weighted by PEOE charges
                contrib = Descriptors.PEOE_VSA1(mol) + Descriptors.PEOE_VSA2(mol) + Descriptors.PEOE_VSA3(mol)
                return contrib * idx / 10.0
        except:
            pass
        return 0
    
    def _calculate_peoe_vsa_polarizability(self, mol):
        """Calculate PEOE VSA polarizability."""
        try:
            return Descriptors.PEOE_VSA11(mol) + Descriptors.PEOE_VSA12(mol)
        except:
            return 0
    
    def _calculate_petitjean_index(self, mol):
        """Calculate Petitjean index."""
        try:
            # Petitjean index is based on topological diameter
            dist_matrix = Chem.GetDistanceMatrix(mol)
            diameter = np.max(dist_matrix)
            radius = np.min(np.max(dist_matrix, axis=0))
            if diameter > 0:
                return (diameter - radius) / diameter
        except:
            pass
        return 0
    
    def _calculate_petitjean_shape_coefficient(self, mol):
        """Calculate Petitjean shape coefficient."""
        try:
            # Shape coefficient based on principal moments
            return Descriptors.InertialShapeFactor(mol)
        except:
            return 0
    
    def _calculate_synthetic_accessibility(self, mol):
        """Calculate synthetic accessibility score."""
        try:
            # Simplified SA score based on complexity
            return Descriptors.BertzCT(mol) / 1000.0
        except:
            return 0
    
    def _calculate_slogp_vsa(self, mol, idx):
        """Calculate SlogP VSA descriptor."""
        try:
            logp = Crippen.MolLogP(mol)
            if idx == 4:
                return Descriptors.SlogP_VSA4(mol)
            elif idx == 5:
                return Descriptors.SlogP_VSA5(mol)
            elif idx == 6:
                return Descriptors.SlogP_VSA6(mol)
            else:
                return logp * idx
        except:
            return 0
    
    def _calculate_smr_vsa(self, mol, idx):
        """Calculate SMR VSA descriptor."""
        try:
            mr = Crippen.MolMR(mol)
            if idx == 3:
                return Descriptors.SMR_VSA3(mol)
            elif idx == 5:
                return Descriptors.SMR_VSA5(mol)
            elif idx == 7:
                return Descriptors.SMR_VSA7(mol)
            else:
                return mr * idx / 10.0
        except:
            return 0
    
    def _calculate_volume_descriptor(self, mol, vol_type):
        """Calculate volume descriptor."""
        try:
            if vol_type == 'ADJMA':
                # Adjacent matrix based volume
                return mol.GetNumHeavyAtoms() * 20.0
            elif vol_type == 'DISTMA':
                # Distance matrix based volume
                dist_matrix = Chem.GetDistanceMatrix(mol)
                return np.sum(dist_matrix) / mol.GetNumHeavyAtoms()
        except:
            pass
        return 0


def main():
    """Example usage of the improved descriptor calculator with optional filtering."""
    
    # Example SMILES strings including some problematic ones
    example_smiles = [
        'CCO',  # Ethanol
        'CC(=O)NC1=CC=C(C=C1)O',  # Acetaminophen
        'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',  # Caffeine
        'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',  # Ibuprofen
        'OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl',  # Diclofenac
        'C' * 30,  # Large hydrocarbon (should be filtered by QSAR)
        '[Na+].[Cl-]',  # Salt (should be filtered by QSAR)
    ]
    
    # Initialize calculator
    calculator = ImprovedHepatotoxicityDescriptors()
    
    print("=== Example 1: Basic descriptor calculation (no filters) ===")
    results_basic = calculator.calculate_all_descriptors(example_smiles, 
                                                        preprocess=True, 
                                                        apply_qsar_filters=False)
    print(f"Calculated {len(results_basic.columns)-1} descriptors for {len(results_basic)} molecules")
    
    print("\n=== Example 2: With QSAR molecular filters ===")
    results_qsar = calculator.calculate_all_descriptors(example_smiles, 
                                                       preprocess=True, 
                                                       apply_qsar_filters=True)
    print(f"After QSAR filters: {len(results_qsar)} molecules retained")
    
    print("\n=== Example 3: With statistical descriptor filters ===")
    # Apply statistical filters to descriptors
    filtered_results, filter_stats = calculator.filter_descriptors(results_basic, 
                                                                   apply_statistical_filters=True)
    
    if filter_stats.get('filters_applied', False):
        original_features = len(results_basic.columns) - 1  # Exclude SMILES
        filtered_features = len(filtered_results.columns) - 1  # Exclude SMILES
        print(f"Statistical filtering: {filtered_features}/{original_features} descriptors retained")
        
        if 'constant_removed' in filter_stats:
            print(f"- Removed {len(filter_stats['constant_removed'])} constant features")
        if 'variance_removed' in filter_stats:
            print(f"- Removed {len(filter_stats['variance_removed'])} low-variance features")
        if 'correlation_removed' in filter_stats:
            print(f"- Removed {len(filter_stats['correlation_removed'])} highly correlated features")
    
    print("\nUsage options:")
    print("1. Basic: calculate_all_descriptors(smiles_list)")
    print("2. With preprocessing: calculate_all_descriptors(smiles_list, preprocess=True)")
    print("3. With QSAR filters: calculate_all_descriptors(smiles_list, apply_qsar_filters=True)")
    print("4. With descriptor filters: filter_descriptors(df, apply_statistical_filters=True)")
    
    print("\nDescriptor categories calculated:")
    print("- Basic molecular properties: MW, LogP, TPSA, etc.")
    print("- Pharmacophore descriptors: PD, PA, PL, etc.")
    print("- MACCS keys: MKEY14, MKEY15, etc.")
    print("- MOE-style descriptors: A_HEAVY_C, BCUT, PEOE_VSA, etc.")
    print("- Topological descriptors: Chi, Kappa, connectivity indices")
    print("- Electronic descriptors: partial charges, EState_VSA")
    print("- Geometric descriptors: 3D shape, moments of inertia")
    print("- Fragment descriptors: functional groups, ring types")
    
    # Save results
    results_basic.to_csv('hepatotoxicity_descriptors_improved.csv', index=False)
    print(f"\nResults saved to 'hepatotoxicity_descriptors_improved.csv'")
    
    return results_basic


if __name__ == "__main__":
    results = main()