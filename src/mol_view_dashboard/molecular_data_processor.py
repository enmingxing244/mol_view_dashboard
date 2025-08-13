#!/usr/bin/env python3
"""
Molecular Data Processor
Handles loading, processing, and analysis of molecular data including:
- CSV data loading and validation
- Molecular descriptor calculation  
- Chemical space analysis (PCA + t-SNE)
- Data preparation for visualization
"""

import pandas as pd
import numpy as np
import logging
from typing import Dict, List, Tuple, Optional, Any
from pathlib import Path
import base64
from io import BytesIO

from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# RDKit imports
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, Lipinski, QED, Draw
    from rdkit.Chem import rdMolDescriptors, rdFingerprintGenerator
    from rdkit.Contrib.SA_Score import sascorer
    from rdkit import DataStructs
except ImportError as e:
    raise ImportError(f"RDKit is required. Install with: conda install -c conda-forge rdkit. Error: {e}")


class MolecularDataProcessor:
    """Processes molecular data for analysis and visualization"""
    
    def __init__(self, config_manager):
        """
        Initialize molecular data processor
        
        Args:
            config_manager: Configuration manager instance
        """
        self.config = config_manager
        self.logger = logging.getLogger(__name__)
        self.df = None
        self.fingerprint_matrix = None
        self.pca_coords = None
        self.tsne_coords = None
        self.pca_model = None
        
    def load_and_process_data(self, csv_file: str) -> pd.DataFrame:
        """
        Load CSV data and process molecular structures
        
        Args:
            csv_file: Path to CSV file containing SMILES and properties
            
        Returns:
            Processed DataFrame with valid molecules
        """
        self.logger.info(f"Loading data from: {csv_file}")
        
        # Load CSV file
        try:
            self.df = pd.read_csv(csv_file)
            self.logger.info(f"Original data shape: {self.df.shape}")
            self.logger.info(f"Columns: {list(self.df.columns)}")
        except Exception as e:
            raise ValueError(f"Error loading CSV file: {e}")
        
        # Validate required columns
        smiles_col = self.config.get('input.smiles_column', 'SMILES')
        if smiles_col not in self.df.columns:
            raise ValueError(f"Required SMILES column '{smiles_col}' not found in CSV")
        
        # Handle optional Name column
        name_col = self.config.get('input.name_column', 'Name')
        if name_col not in self.df.columns:
            self.df[name_col] = [f"Compound_{i+1}" for i in range(len(self.df))]
            self.logger.info(f"Generated {name_col} column")
        
        # Map Name to title for internal consistency
        self.df['title'] = self.df[name_col]
        
        # Create molecule objects from SMILES
        self.logger.info("Creating molecule objects from SMILES...")
        self.df['Mol'] = self.df[smiles_col].apply(
            lambda x: Chem.MolFromSmiles(x) if pd.notna(x) else None
        )
        
        # Report failed SMILES
        failed_smiles = self.df[self.df['Mol'].isna()]
        if len(failed_smiles) > 0:
            self.logger.warning(f"Failed to parse {len(failed_smiles)} SMILES")
            for idx, row in failed_smiles.head().iterrows():
                self.logger.warning(f"  {idx}: {row[smiles_col]}")
        
        # Filter valid molecules
        self.df = self.df[self.df['Mol'].notna()].reset_index(drop=True)
        self.logger.info(f"Valid molecules: {len(self.df)}")
        
        return self.df
    
    def calculate_molecular_descriptors(self) -> pd.DataFrame:
        """
        Calculate comprehensive molecular descriptors
        
        Returns:
            DataFrame with calculated descriptors
        """
        if self.df is None:
            raise ValueError("No molecular data loaded. Call load_and_process_data() first.")
        
        if not self.config.get('analysis.calculate_descriptors', True):
            self.logger.info("Descriptor calculation disabled in configuration")
            return self.df
        
        self.logger.info(f"Calculating molecular descriptors for {len(self.df)} compounds...")
        
        descriptors = []
        
        for i, mol in enumerate(self.df['Mol']):
            if i % 100 == 0:
                self.logger.info(f"  Processing molecule {i+1}/{len(self.df)}")
            
            try:
                desc = {
                    # Basic properties
                    'MW': Descriptors.MolWt(mol),
                    'LogP': Crippen.MolLogP(mol),
                    'TPSA': rdMolDescriptors.CalcTPSA(mol),
                    'HBA': rdMolDescriptors.CalcNumHBA(mol),
                    'HBD': rdMolDescriptors.CalcNumHBD(mol),
                    'RotBonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
                    'NumRings': rdMolDescriptors.CalcNumRings(mol),
                    
                    # Drug-likeness
                    'QED': QED.qed(mol),
                    'SAscore': sascorer.calculateScore(mol),
                }
                
            except Exception as e:
                self.logger.warning(f"Error calculating descriptors for molecule {i}: {e}")
                desc = {key: np.nan for key in [
                    'MW', 'LogP', 'TPSA', 'HBA', 'HBD', 'RotBonds', 'NumRings', 'QED', 'SAscore'
                ]}
            
            descriptors.append(desc)
        
        # Add descriptors to dataframe
        desc_df = pd.DataFrame(descriptors)
        for col in desc_df.columns:
            self.df[col] = desc_df[col]
        
        # For testing purposes, keep all compounds even if some descriptors failed
        self.logger.info(f"Final dataset: {len(self.df)} compounds (with some potential missing descriptors)")
        
        # Log descriptor summary
        self.logger.info("Descriptor summary:")
        summary_cols = ['MW', 'LogP', 'TPSA', 'QED', 'SAscore', 'HBA', 'HBD', 'RotBonds', 'NumRings']
        available_cols = [col for col in summary_cols if col in self.df.columns]
        if available_cols:
            self.logger.info(f"\n{self.df[available_cols].describe()}")
        
        return self.df
    
    def calculate_molecular_fingerprints(self) -> np.ndarray:
        """
        Calculate molecular fingerprints for chemical space analysis
        
        Returns:
            Fingerprint matrix (n_compounds x n_bits)
        """
        if self.df is None:
            raise ValueError("No molecular data loaded. Call load_and_process_data() first.")
        
        self.logger.info(f"Calculating molecular fingerprints for {len(self.df)} compounds...")
        
        # Get fingerprint configuration
        fp_config = self.config.get('analysis.fingerprints', {})
        fp_type = fp_config.get('type', 'morgan')
        radius = fp_config.get('radius', 2)
        n_bits = fp_config.get('n_bits', 2048)
        
        # Generate fingerprints based on type
        if fp_type.lower() == 'morgan':
            generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
        elif fp_type.lower() == 'rdkit':
            generator = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=n_bits)
        else:
            self.logger.warning(f"Unknown fingerprint type '{fp_type}', using Morgan")
            generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
        
        fingerprints = []
        
        for i, mol in enumerate(self.df['Mol']):
            if i % 100 == 0:
                self.logger.info(f"  Processing molecule {i+1}/{len(self.df)}")
            
            fp = generator.GetFingerprint(mol)
            fp_array = np.zeros((n_bits,))
            DataStructs.ConvertToNumpyArray(fp, fp_array)
            fingerprints.append(fp_array)
        
        self.fingerprint_matrix = np.array(fingerprints)
        self.logger.info(f"Fingerprint matrix shape: {self.fingerprint_matrix.shape}")
        
        return self.fingerprint_matrix
    
    def perform_chemical_space_analysis(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Perform PCA and t-SNE analysis on molecular fingerprints
        
        Returns:
            Tuple of (pca_coords, tsne_coords, pca_explained_variance)
        """
        if self.fingerprint_matrix is None:
            self.calculate_molecular_fingerprints()
        
        if len(self.df) == 0:
            self.logger.warning("No compounds available for chemical space analysis")
            return np.array([]), np.array([]), np.array([])
        
        if len(self.df) < 2:
            self.logger.warning("Insufficient compounds for chemical space analysis (need at least 2)")
            return np.array([]), np.array([]), np.array([])
        
        self.logger.info("Performing chemical space analysis...")
        
        # Standardize features
        scaler = StandardScaler()
        fp_scaled = scaler.fit_transform(self.fingerprint_matrix)
        
        # PCA analysis
        pca_config = self.config.get('analysis.chemical_space.pca', {})
        if pca_config.get('enabled', True):
            self.logger.info("  Running PCA...")
            n_components = pca_config.get('n_components', 2)
            
            self.pca_model = PCA(n_components=n_components, random_state=42)
            self.pca_coords = self.pca_model.fit_transform(fp_scaled)
            
            explained_variance = self.pca_model.explained_variance_ratio_
            self.logger.info(f"  PCA explained variance: {explained_variance}")
            self.logger.info(f"  Total explained variance: {sum(explained_variance):.3f}")
            
            # Add PCA coordinates to dataframe
            for i in range(n_components):
                self.df[f'PCA_{i+1}'] = self.pca_coords[:, i]
        
        # t-SNE analysis
        tsne_config = self.config.get('analysis.chemical_space.tsne', {})
        if tsne_config.get('enabled', True) and len(self.df) >= 10:
            self.logger.info("  Running t-SNE...")
            n_components = tsne_config.get('n_components', 2)
            perplexity = min(tsne_config.get('perplexity', 30), len(self.df)//4)
            random_state = tsne_config.get('random_state', 42)
            
            tsne = TSNE(
                n_components=n_components, 
                random_state=random_state, 
                perplexity=perplexity
            )
            self.tsne_coords = tsne.fit_transform(fp_scaled)
            
            # Add t-SNE coordinates to dataframe
            for i in range(n_components):
                self.df[f'tSNE_{i+1}'] = self.tsne_coords[:, i]
        else:
            if len(self.df) < 10:
                self.logger.info(f"  t-SNE skipped: requires at least 10 compounds (found {len(self.df)})")
            else:
                self.logger.info("  t-SNE disabled in configuration")
        
        self.logger.info("Chemical space analysis completed")
        
        return (
            self.pca_coords, 
            self.tsne_coords if self.tsne_coords is not None else np.array([]), 
            self.pca_model.explained_variance_ratio_ if self.pca_model else np.array([])
        )
    
    def mol_to_base64_png(self, mol, size=(200, 140)) -> str:
        """
        Convert molecule to base64 encoded PNG image
        
        Args:
            mol: RDKit molecule object
            size: Image size tuple (width, height)
            
        Returns:
            Base64 encoded image string
        """
        if mol is None:
            return ""
        
        try:
            img = Draw.MolToImage(mol, size=size)
            buffered = BytesIO()
            img.save(buffered, format="PNG")
            img_str = base64.b64encode(buffered.getvalue()).decode()
            return f"data:image/png;base64,{img_str}"
        except Exception as e:
            self.logger.warning(f"Error generating structure image: {e}")
            return ""
    
    def prepare_visualization_data(self) -> List[Dict[str, Any]]:
        """
        Prepare data for visualization dashboard
        
        Returns:
            List of dictionaries containing compound data for visualization
        """
        if self.df is None:
            raise ValueError("No molecular data processed")
        
        self.logger.info("Preparing visualization data...")
        
        data_points = []
        smiles_col = self.config.get('input.smiles_column', 'SMILES')
        
        for idx, row in self.df.iterrows():
            if idx % 50 == 0:
                self.logger.info(f"  Processing compound {idx+1}/{len(self.df)}")
            
            # Generate structure image
            img_base64 = self.mol_to_base64_png(row['Mol'])
            
            # Base data point
            data_point = {
                'id': int(idx),
                'title': str(row.get('title', f"Compound_{idx+1}")),
                'smiles': str(row[smiles_col]),
            }
            
            # Add molecular descriptors
            descriptor_cols = [
                'MW', 'LogP', 'TPSA', 'HBA', 'HBD', 'QED', 'SAscore', 
                'RotBonds', 'NumRings'
            ]
            
            for col in descriptor_cols:
                if col in row and pd.notna(row[col]):
                    data_point[col.lower()] = float(row[col])
            
            # Add chemical space coordinates
            if 'PCA_1' in row and pd.notna(row['PCA_1']):
                data_point['pca_x'] = float(row['PCA_1'])
                data_point['pca_y'] = float(row['PCA_2'])
            
            if 'tSNE_1' in row and pd.notna(row['tSNE_1']):
                data_point['tsne_x'] = float(row['tSNE_1'])
                data_point['tsne_y'] = float(row['tSNE_2'])
            
            # Add custom property columns
            property_cols = self.config.get('input.property_columns', []) or []
            for prop_col in property_cols:
                if prop_col in row and pd.notna(row[prop_col]):
                    # Normalize property name for JavaScript
                    prop_key = prop_col.lower().replace(' ', '_').replace('-', '_')
                    try:
                        data_point[prop_key] = float(row[prop_col])
                    except (ValueError, TypeError):
                        data_point[prop_key] = str(row[prop_col])
            
            # Add structure image
            data_point['image'] = img_base64
            
            data_points.append(data_point)
        
        self.logger.info(f"Visualization data prepared for {len(data_points)} compounds")
        return data_points
    
    def get_available_properties(self) -> List[str]:
        """
        Get list of available properties for plotting
        
        Returns:
            List of property column names
        """
        if self.df is None:
            return []
        
        # Standard molecular descriptors
        standard_props = [
            'MW', 'LogP', 'TPSA', 'HBA', 'HBD', 'QED', 'SAscore',
            'RotBonds', 'NumRings'
        ]
        
        # Custom properties from input
        custom_props = self.config.get('input.property_columns', []) or []
        
        # Filter to only include properties that exist in the dataframe
        available_props = []
        for prop in standard_props + custom_props:
            if prop in self.df.columns:
                available_props.append(prop)
        
        return available_props
    
    def get_data_summary(self) -> Dict[str, Any]:
        """
        Get summary statistics of the dataset
        
        Returns:
            Dictionary with dataset summary information
        """
        if self.df is None:
            return {}
        
        summary = {
            'total_compounds': len(self.df),
            'available_properties': self.get_available_properties(),
            'has_pca': 'PCA_1' in self.df.columns,
            'has_tsne': 'tSNE_1' in self.df.columns,
            'pca_variance': self.pca_model.explained_variance_ratio_.tolist() if self.pca_model else [],
        }
        
        # Add property ranges for the most common descriptors
        common_props = ['MW', 'LogP', 'TPSA', 'QED']
        property_ranges = {}
        
        for prop in common_props:
            if prop in self.df.columns:
                property_ranges[prop] = {
                    'min': float(self.df[prop].min()),
                    'max': float(self.df[prop].max()),
                    'mean': float(self.df[prop].mean()),
                    'std': float(self.df[prop].std())
                }
        
        summary['property_ranges'] = property_ranges
        
        return summary