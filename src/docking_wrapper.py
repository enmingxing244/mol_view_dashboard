#!/usr/bin/env python3
"""
Molecular Docking Wrapper
Handles AutoDock Vina integration for molecular docking analysis
- Prepares ligands from SMILES
- Runs Vina docking calculations
- Processes docking results
- Integrates results with visualization data
"""

import os
import subprocess
import pandas as pd
import numpy as np
import logging
from typing import Dict, List, Tuple, Optional, Any
from pathlib import Path
import tempfile
import shutil

# RDKit imports for molecule preparation
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolAlign
    from rdkit.Chem.rdMolAlign import AlignMol
except ImportError as e:
    logging.warning(f"RDKit not available for docking: {e}")


class DockingError(Exception):
    """Custom exception for docking-related errors"""
    pass


class VinaDockingWrapper:
    """Wrapper for AutoDock Vina molecular docking"""
    
    def __init__(self, config_manager):
        """
        Initialize docking wrapper
        
        Args:
            config_manager: Configuration manager instance
        """
        self.config = config_manager
        self.logger = logging.getLogger(__name__)
        self.docking_results = []
        self.temp_dir = None
        
        # Check if docking is enabled
        if not self.config.is_docking_enabled():
            self.logger.info("Docking is disabled in configuration")
            return
        
        # Validate docking requirements
        self._validate_docking_setup()
    
    def _validate_docking_setup(self) -> None:
        """
        Validate that all required files and tools are available for docking
        
        Raises:
            DockingError: If docking setup is invalid
        """
        # Check for required files
        protein_pdb = self.config.get('docking.protein_pdb')
        vina_config = self.config.get('docking.vina_config')
        
        if not protein_pdb or not Path(protein_pdb).exists():
            raise DockingError(f"Protein PDB file not found: {protein_pdb}")
        
        if not vina_config or not Path(vina_config).exists():
            raise DockingError(f"Vina configuration file not found: {vina_config}")
        
        # Check for Vina executable
        if not self._check_vina_executable():
            raise DockingError("AutoDock Vina executable not found in PATH")
        
        self.logger.info("Docking setup validation completed")
    
    def _check_vina_executable(self) -> bool:
        """
        Check if AutoDock Vina is available
        
        Returns:
            True if Vina executable is found, False otherwise
        """
        try:
            result = subprocess.run(['vina', '--help'], 
                                 capture_output=True, 
                                 text=True, 
                                 timeout=10)
            return result.returncode == 0
        except (subprocess.SubprocessError, FileNotFoundError):
            return False
    
    def prepare_ligands(self, df: pd.DataFrame) -> List[str]:
        """
        Prepare ligand files from SMILES strings
        
        Args:
            df: DataFrame containing SMILES and molecule data
            
        Returns:
            List of prepared ligand file paths
        """
        if not self.config.is_docking_enabled():
            return []
        
        self.logger.info(f"Preparing ligands for {len(df)} compounds...")
        
        # Create temporary directory for ligand files
        self.temp_dir = tempfile.mkdtemp(prefix="mol_docking_")
        self.logger.info(f"Using temporary directory: {self.temp_dir}")
        
        ligand_files = []
        smiles_col = self.config.get('input.smiles_column', 'SMILES')
        
        for idx, row in df.iterrows():
            if idx % 10 == 0:
                self.logger.info(f"  Preparing ligand {idx+1}/{len(df)}")
            
            try:
                mol = row['Mol']
                if mol is None:
                    continue
                
                # Add hydrogens and generate 3D coordinates
                mol_h = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol_h, randomSeed=42)
                AllChem.MMFFOptimizeMolecule(mol_h)
                
                # Write to SDF file
                ligand_file = Path(self.temp_dir) / f"ligand_{idx:04d}.sdf"
                writer = Chem.SDWriter(str(ligand_file))
                writer.write(mol_h)
                writer.close()
                
                ligand_files.append(str(ligand_file))
                
            except Exception as e:
                self.logger.warning(f"Failed to prepare ligand {idx}: {e}")
                continue
        
        self.logger.info(f"Successfully prepared {len(ligand_files)} ligand files")
        return ligand_files
    
    def run_vina_docking(self, ligand_files: List[str]) -> List[Dict[str, Any]]:
        """
        Run AutoDock Vina docking for prepared ligands
        
        Args:
            ligand_files: List of ligand file paths
            
        Returns:
            List of docking results dictionaries
        """
        if not self.config.is_docking_enabled() or not ligand_files:
            return []
        
        self.logger.info(f"Running Vina docking for {len(ligand_files)} ligands...")
        
        # Get docking configuration
        protein_pdb = self.config.get('docking.protein_pdb')
        vina_config = self.config.get('docking.vina_config')
        output_dir = Path(self.config.get('docking.output_dir', 'docking_results'))
        
        # Create output directory
        output_dir.mkdir(exist_ok=True)
        
        # Docking parameters
        params = self.config.get('docking.parameters', {})
        exhaustiveness = params.get('exhaustiveness', 8)
        num_modes = params.get('num_modes', 9)
        energy_range = params.get('energy_range', 3)
        
        docking_results = []
        
        for i, ligand_file in enumerate(ligand_files):
            if i % 5 == 0:
                self.logger.info(f"  Docking ligand {i+1}/{len(ligand_files)}")
            
            try:
                # Extract ligand index from filename
                ligand_idx = int(Path(ligand_file).stem.split('_')[1])
                
                # Output file for docking result
                output_file = output_dir / f"result_{ligand_idx:04d}.pdbqt"
                log_file = output_dir / f"log_{ligand_idx:04d}.txt"
                
                # Build Vina command
                vina_cmd = [
                    'vina',
                    '--config', str(vina_config),
                    '--receptor', str(protein_pdb),
                    '--ligand', str(ligand_file),
                    '--out', str(output_file),
                    '--log', str(log_file),
                    '--exhaustiveness', str(exhaustiveness),
                    '--num_modes', str(num_modes),
                    '--energy_range', str(energy_range)
                ]
                
                # Run Vina
                result = subprocess.run(vina_cmd, 
                                      capture_output=True, 
                                      text=True, 
                                      timeout=300)  # 5 minute timeout per ligand
                
                if result.returncode == 0:
                    # Parse docking results
                    binding_energy = self._parse_vina_output(str(log_file))
                    
                    docking_result = {
                        'ligand_index': ligand_idx,
                        'binding_energy': binding_energy,
                        'output_file': str(output_file),
                        'log_file': str(log_file),
                        'success': True
                    }
                    
                else:
                    self.logger.warning(f"Vina failed for ligand {ligand_idx}: {result.stderr}")
                    docking_result = {
                        'ligand_index': ligand_idx,
                        'binding_energy': np.nan,
                        'success': False,
                        'error': result.stderr
                    }
                
                docking_results.append(docking_result)
                
            except subprocess.TimeoutExpired:
                self.logger.warning(f"Vina timeout for ligand {i}")
                docking_results.append({
                    'ligand_index': i,
                    'binding_energy': np.nan,
                    'success': False,
                    'error': 'Timeout'
                })
                
            except Exception as e:
                self.logger.warning(f"Error docking ligand {i}: {e}")
                docking_results.append({
                    'ligand_index': i,
                    'binding_energy': np.nan,
                    'success': False,
                    'error': str(e)
                })
        
        successful_dockings = sum(1 for r in docking_results if r.get('success', False))
        self.logger.info(f"Completed docking: {successful_dockings}/{len(ligand_files)} successful")
        
        self.docking_results = docking_results
        return docking_results
    
    def _parse_vina_output(self, log_file: str) -> float:
        """
        Parse Vina log file to extract binding energy
        
        Args:
            log_file: Path to Vina log file
            
        Returns:
            Best binding energy (kcal/mol)
        """
        try:
            with open(log_file, 'r') as f:
                lines = f.readlines()
            
            # Look for the results table
            for i, line in enumerate(lines):
                if 'mode' in line and 'affinity' in line:
                    # Next line should contain the best result
                    if i + 1 < len(lines):
                        result_line = lines[i + 1].strip()
                        parts = result_line.split()
                        if len(parts) >= 2:
                            try:
                                return float(parts[1])  # Binding affinity
                            except ValueError:
                                pass
            
            # If table format not found, look for other patterns
            for line in lines:
                if 'Estimated Free Energy of Binding' in line:
                    parts = line.split()
                    for j, part in enumerate(parts):
                        if 'kcal/mol' in part and j > 0:
                            try:
                                return float(parts[j-1])
                            except ValueError:
                                pass
            
            self.logger.warning(f"Could not parse binding energy from {log_file}")
            return np.nan
            
        except Exception as e:
            self.logger.warning(f"Error parsing Vina output {log_file}: {e}")
            return np.nan
    
    def integrate_docking_results(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Integrate docking results into the main dataframe
        
        Args:
            df: Main dataframe with molecular data
            
        Returns:
            DataFrame with docking results added
        """
        if not self.config.is_docking_enabled() or not self.docking_results:
            # Add empty docking columns
            df['binding_energy'] = np.nan
            df['docking_success'] = False
            return df
        
        self.logger.info("Integrating docking results...")
        
        # Create a mapping from ligand index to docking results
        docking_map = {r['ligand_index']: r for r in self.docking_results}
        
        # Add docking columns
        binding_energies = []
        docking_success = []
        
        for idx in range(len(df)):
            if idx in docking_map:
                result = docking_map[idx]
                binding_energies.append(result.get('binding_energy', np.nan))
                docking_success.append(result.get('success', False))
            else:
                binding_energies.append(np.nan)
                docking_success.append(False)
        
        df['binding_energy'] = binding_energies
        df['docking_success'] = docking_success
        
        # Log summary
        successful_count = sum(docking_success)
        if successful_count > 0:
            valid_energies = [e for e in binding_energies if not np.isnan(e)]
            if valid_energies:
                self.logger.info(f"Docking integration complete:")
                self.logger.info(f"  Successful dockings: {successful_count}/{len(df)}")
                self.logger.info(f"  Binding energy range: {min(valid_energies):.2f} to {max(valid_energies):.2f} kcal/mol")
                self.logger.info(f"  Mean binding energy: {np.mean(valid_energies):.2f} kcal/mol")
        
        return df
    
    def prepare_docking_visualization_data(self, df: pd.DataFrame) -> List[Dict[str, Any]]:
        """
        Prepare docking data for visualization
        
        Args:
            df: DataFrame with docking results
            
        Returns:
            List of docking visualization data
        """
        if not self.config.is_docking_enabled():
            return []
        
        viz_data = []
        output_dir = Path(self.config.get('docking.output_dir', 'docking_results'))
        
        for idx, row in df.iterrows():
            if row.get('docking_success', False):
                # Find corresponding output files
                result_file = output_dir / f"result_{idx:04d}.pdbqt"
                
                if result_file.exists():
                    viz_data.append({
                        'compound_id': idx,
                        'compound_name': row.get('title', f"Compound_{idx}"),
                        'binding_energy': row.get('binding_energy', np.nan),
                        'pose_file': str(result_file),
                        'smiles': row.get(self.config.get('input.smiles_column', 'SMILES'), '')
                    })
        
        self.logger.info(f"Prepared docking visualization data for {len(viz_data)} compounds")
        return viz_data
    
    def cleanup(self) -> None:
        """Clean up temporary files"""
        if self.temp_dir and Path(self.temp_dir).exists():
            try:
                shutil.rmtree(self.temp_dir)
                self.logger.info(f"Cleaned up temporary directory: {self.temp_dir}")
            except Exception as e:
                self.logger.warning(f"Failed to clean up temporary directory: {e}")
    
    def __del__(self):
        """Destructor to ensure cleanup"""
        self.cleanup()
    
    def get_docking_summary(self) -> Dict[str, Any]:
        """
        Get summary of docking results
        
        Returns:
            Dictionary with docking summary statistics
        """
        if not self.config.is_docking_enabled() or not self.docking_results:
            return {'enabled': False}
        
        successful_results = [r for r in self.docking_results if r.get('success', False)]
        binding_energies = [r['binding_energy'] for r in successful_results 
                          if not np.isnan(r.get('binding_energy', np.nan))]
        
        summary = {
            'enabled': True,
            'total_compounds': len(self.docking_results),
            'successful_dockings': len(successful_results),
            'success_rate': len(successful_results) / len(self.docking_results) if self.docking_results else 0,
        }
        
        if binding_energies:
            summary.update({
                'best_binding_energy': min(binding_energies),
                'worst_binding_energy': max(binding_energies),
                'mean_binding_energy': np.mean(binding_energies),
                'std_binding_energy': np.std(binding_energies)
            })
        
        return summary