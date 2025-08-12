#!/usr/bin/env python3
"""
Molecular Docking Wrapper
Handles AutoDock Vina integration with MGLTools preparation
- Prepares receptors using MGLTools (clean water/ions, add charges)
- Prepares ligands from SMILES using MGLTools prepare_ligand4.py
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
    """Wrapper for AutoDock Vina molecular docking with MGLTools preparation"""
    
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
        
        # MGLTools paths
        self.mgltools_path = "/Users/enmingxing/Projects/mol_view_dashboard/packages/mgltools_1.5.7_MacOS-X/installed"
        self.mgl_python = f"{self.mgltools_path}/bin/python"
        self.utilities_path = f"{self.mgltools_path}/MGLToolsPckgs/AutoDockTools/Utilities24"
        
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
        
        if not protein_pdb or not Path(protein_pdb).exists():
            raise DockingError(f"Protein PDB file not found: {protein_pdb}")
        
        # Check for MGLTools installation
        if not Path(self.mgl_python).exists():
            raise DockingError(f"MGLTools Python not found: {self.mgl_python}")
        
        if not Path(self.utilities_path).exists():
            raise DockingError(f"MGLTools utilities not found: {self.utilities_path}")
        
        # Check for Vina executable
        if not self._check_vina_executable():
            vina_path = self.config.get('docking.vina_executable')
            if vina_path:
                raise DockingError(f"AutoDock Vina executable not found at: {vina_path}")
            else:
                raise DockingError("AutoDock Vina executable not found. Please specify vina_executable in configuration.")
        
        # Check for required docking parameters
        center = self.config.get('docking.binding_site.center')
        size = self.config.get('docking.binding_site.size')
        
        if not center or len(center) != 3:
            raise DockingError("Binding site center (x, y, z) must be specified in configuration")
        
        if not size or len(size) != 3:
            raise DockingError("Binding site size (x, y, z) must be specified in configuration")
        
        self.logger.info("Docking setup validation completed")
    
    def _check_vina_executable(self) -> bool:
        """
        Check if AutoDock Vina is available
        
        Returns:
            True if Vina executable is found, False otherwise
        """
        # Get custom vina executable path from config
        vina_executable = self.config.get('docking.vina_executable', 'vina')
        
        try:
            result = subprocess.run([vina_executable, '--help'], 
                                 capture_output=True, 
                                 text=True, 
                                 timeout=10)
            return result.returncode == 0
        except (subprocess.SubprocessError, FileNotFoundError):
            return False
    
    def _get_vina_executable(self) -> str:
        """
        Get the Vina executable path
        
        Returns:
            Path to Vina executable
        """
        return self.config.get('docking.vina_executable', 'vina')
    
    def prepare_receptor(self) -> str:
        """
        Prepare receptor using MGLTools - clean water/ions and convert to PDBQT
        
        Returns:
            Path to prepared receptor PDBQT file
        """
        if not self.config.is_docking_enabled():
            return ""
        
        protein_pdb = self.config.get('docking.protein_pdb')
        self.logger.info(f"Preparing receptor from: {protein_pdb}")
        
        # Create temp directory if not exists
        if not self.temp_dir:
            self.temp_dir = tempfile.mkdtemp(prefix="mol_docking_")
        
        # Output PDBQT file path
        receptor_pdbqt = Path(self.temp_dir) / "receptor.pdbqt"
        
        # MGLTools prepare_receptor4.py command
        prepare_receptor_script = f"{self.utilities_path}/prepare_receptor4.py"
        
        # Get receptor preparation options from config
        receptor_prep = self.config.get('docking.mgltools.receptor_prep', {})
        remove_ligands = receptor_prep.get('remove_ligands', True)
        
        cmd = [
            self.mgl_python,
            prepare_receptor_script,
            '-r', str(protein_pdb),
            '-o', str(receptor_pdbqt),
            '-A', 'hydrogens',  # Add hydrogens
            '-U', 'waters',     # Remove waters
            '-U', 'nonstdres'   # Remove non-standard residues
        ]
        
        # Add ligand removal option if enabled
        if remove_ligands:
            cmd.extend(['-U', 'lps'])  # Remove ligands and cofactors
            self.logger.info("Receptor preparation will remove: waters, non-standard residues, and native ligands/cofactors")
        else:
            self.logger.info("Receptor preparation will remove: waters and non-standard residues (keeping native ligands)")
        
        self.logger.info(f"Running: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            
            if result.returncode == 0 and receptor_pdbqt.exists():
                self.logger.info(f"Receptor prepared successfully: {receptor_pdbqt}")
                return str(receptor_pdbqt)
            else:
                raise DockingError(f"Receptor preparation failed: {result.stderr}")
                
        except subprocess.TimeoutExpired:
            raise DockingError("Receptor preparation timed out")
        except Exception as e:
            raise DockingError(f"Error preparing receptor: {e}")
    
    def prepare_ligands(self, df: pd.DataFrame) -> List[str]:
        """
        Prepare ligand files from SMILES using MGLTools
        
        Args:
            df: DataFrame containing SMILES and molecule data
            
        Returns:
            List of prepared ligand PDBQT file paths
        """
        if not self.config.is_docking_enabled():
            return []
        
        self.logger.info(f"Preparing ligands for {len(df)} compounds using MGLTools...")
        
        # Create temporary directory for ligand files
        if not self.temp_dir:
            self.temp_dir = tempfile.mkdtemp(prefix="mol_docking_")
        self.logger.info(f"Using temporary directory: {self.temp_dir}")
        
        ligand_pdbqt_files = []
        smiles_col = self.config.get('input.smiles_column', 'SMILES')
        prepare_ligand_script = f"{self.utilities_path}/prepare_ligand4.py"
        
        for idx, row in df.iterrows():
            if idx % 10 == 0:
                self.logger.info(f"  Preparing ligand {idx+1}/{len(df)}")
            
            try:
                mol = row['Mol']
                if mol is None:
                    continue
                
                # Step 1: Generate 3D structure using RDKit
                mol_h = Chem.AddHs(mol)
                
                # Try to embed molecule with multiple attempts
                embed_result = AllChem.EmbedMolecule(mol_h, randomSeed=42)
                if embed_result == -1:
                    # Try with different parameters
                    embed_result = AllChem.EmbedMolecule(mol_h, randomSeed=42, useRandomCoords=True)
                    if embed_result == -1:
                        self.logger.warning(f"Failed to embed 3D coordinates for ligand {idx}")
                        continue
                
                # Optimize geometry
                try:
                    AllChem.MMFFOptimizeMolecule(mol_h)
                except Exception as e:
                    self.logger.warning(f"MMFF optimization failed for ligand {idx}: {e}")
                    # Continue anyway, the embedded coordinates might still work
                
                # Step 2: Write to PDB file using working directory approach
                compound_name = row.get('title', f"ligand_{idx:04d}")
                # Sanitize compound name for filename
                safe_name = "".join(c for c in compound_name if c.isalnum() or c in ('-', '_'))[:20]
                if not safe_name:
                    safe_name = f"ligand_{idx:04d}"
                
                ligand_pdb = Path(self.temp_dir) / f"{safe_name}.pdb"
                ligand_pdbqt = Path(self.temp_dir) / f"{safe_name}.pdbqt"
                
                # Write PDB file
                pdb_block = Chem.MolToPDBBlock(mol_h)
                if not pdb_block or pdb_block.strip() == "":
                    self.logger.warning(f"Empty PDB block generated for ligand {idx}")
                    continue
                    
                with open(ligand_pdb, 'w') as f:
                    f.write(pdb_block)
                
                if not ligand_pdb.exists() or ligand_pdb.stat().st_size == 0:
                    self.logger.warning(f"Failed to create PDB file for ligand {idx}")
                    continue
                
                # Step 3: Convert PDB to PDBQT using MGLTools (run from temp directory)
                original_dir = os.getcwd()
                try:
                    os.chdir(self.temp_dir)
                    
                    cmd = [
                        self.mgl_python,
                        prepare_ligand_script,
                        '-l', f"{safe_name}.pdb",  # Use relative path
                        '-o', f"{safe_name}.pdbqt",
                        '-A', 'hydrogens'
                    ]
                    
                    result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
                    
                    if result.returncode == 0 and ligand_pdbqt.exists():
                        ligand_pdbqt_files.append(str(ligand_pdbqt))
                        self.logger.info(f"Successfully prepared ligand {idx}: {ligand_pdbqt}")
                    else:
                        self.logger.warning(f"Failed to prepare ligand {idx} with MGLTools: {result.stderr}")
                        continue
                        
                finally:
                    os.chdir(original_dir)
                
            except Exception as e:
                self.logger.warning(f"Failed to prepare ligand {idx}: {e}")
                continue
        
        self.logger.info(f"Successfully prepared {len(ligand_pdbqt_files)} ligand PDBQT files")
        return ligand_pdbqt_files
    
    def run_vina_docking(self, ligand_files: List[str], receptor_file: str) -> List[Dict[str, Any]]:
        """
        Run AutoDock Vina docking for prepared ligands
        
        Args:
            ligand_files: List of ligand PDBQT file paths
            receptor_file: Path to prepared receptor PDBQT file
            
        Returns:
            List of docking results dictionaries
        """
        if not self.config.is_docking_enabled() or not ligand_files:
            return []
        
        self.logger.info(f"Running Vina docking for {len(ligand_files)} ligands...")
        
        # Get docking configuration
        output_dir = Path(self.config.get('docking.output_dir', 'docking_results'))
        
        # Create output directory
        output_dir.mkdir(exist_ok=True)
        
        # Docking parameters
        params = self.config.get('docking.parameters', {})
        exhaustiveness = params.get('exhaustiveness', 8)
        num_modes = params.get('num_modes', 9)
        energy_range = params.get('energy_range', 3)
        
        # Binding site configuration
        center = self.config.get('docking.binding_site.center')
        size = self.config.get('docking.binding_site.size')
        
        docking_results = []
        
        for i, ligand_file in enumerate(ligand_files):
            if i % 5 == 0:
                self.logger.info(f"  Docking ligand {i+1}/{len(ligand_files)}")
            
            try:
                # Use enumeration index as ligand identifier
                ligand_idx = i
                
                # Get compound name from filename for better organization
                compound_name = Path(ligand_file).stem
                
                # Output file for docking result
                output_file = output_dir / f"result_{ligand_idx:04d}_{compound_name}.pdbqt"
                log_file = output_dir / f"log_{ligand_idx:04d}_{compound_name}.txt"
                
                # Build Vina command with binding site parameters
                vina_executable = self._get_vina_executable()
                vina_cmd = [
                    vina_executable,
                    '--receptor', str(receptor_file),
                    '--ligand', str(ligand_file),
                    '--out', str(output_file),
                    '--log', str(log_file),
                    '--center_x', str(center[0]),
                    '--center_y', str(center[1]),
                    '--center_z', str(center[2]),
                    '--size_x', str(size[0]),
                    '--size_y', str(size[1]),
                    '--size_z', str(size[2]),
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
                    # Skip the units line and get the first data line
                    if i + 2 < len(lines):
                        # Skip the "|  (kcal/mol) | rmsd l.b.| rmsd u.b." line
                        # Skip the "-----+------------+----------+----------" line  
                        # Get the actual first result line
                        result_line = lines[i + 3].strip()
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