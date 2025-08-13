#!/usr/bin/env python3
"""
PDB Cleaner - Removes alternative conformations and problematic residues

This script processes PDB files to resolve issues that can cause problems
with MGLTools receptor preparation:

1. Alternative conformations (A, B, C conformers) - keeps highest occupancy
2. Missing atoms in residues
3. Non-standard residues
4. Waters and ions (optional)

Usage:
    python pdb_cleaner.py input.pdb output.pdb [--remove-waters] [--remove-ions]
"""

import os
import sys
import argparse
import logging
from pathlib import Path
from typing import Dict, List, Set, Tuple
from collections import defaultdict


class PDBCleaner:
    """Clean PDB files for molecular docking preparation"""
    
    def __init__(self):
        self.logger = logging.getLogger(__name__)
        
        # Standard amino acid residues
        self.standard_residues = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS',
            'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
            'TYR', 'VAL'
        }
        
        # Common ions and cofactors to potentially remove
        self.common_ions = {
            'CA', 'MG', 'ZN', 'FE', 'MN', 'CU', 'NA', 'K', 'CL', 'SO4',
            'PO4', 'CO3', 'NO3', 'BR', 'I', 'F'
        }
        
    def clean_pdb(self, input_file: str, output_file: str, 
                  remove_waters: bool = True, remove_ions: bool = False,
                  remove_ligands: bool = False) -> None:
        """
        Clean PDB file by resolving alternative conformations and other issues
        
        Args:
            input_file: Path to input PDB file
            output_file: Path to output cleaned PDB file
            remove_waters: Remove water molecules (HOH)
            remove_ions: Remove common ions and salts
            remove_ligands: Remove non-protein ligands (HETATM except ions/waters)
        """
        self.logger.info(f"Cleaning PDB file: {input_file}")
        
        # Read and parse PDB file
        with open(input_file, 'r') as f:
            lines = f.readlines()
        
        # Process alternative conformations
        cleaned_lines = self._resolve_alternative_conformations(lines)
        
        # Filter out unwanted records
        final_lines = self._filter_records(cleaned_lines, remove_waters, 
                                         remove_ions, remove_ligands)
        
        # Write cleaned PDB
        with open(output_file, 'w') as f:
            f.writelines(final_lines)
        
        self.logger.info(f"Cleaned PDB saved to: {output_file}")
        
        # Report cleaning statistics
        self._report_cleaning_stats(lines, final_lines)
    
    def _resolve_alternative_conformations(self, lines: List[str]) -> List[str]:
        """
        Resolve alternative conformations by selecting the highest occupancy
        
        Args:
            lines: List of PDB file lines
            
        Returns:
            List of cleaned lines with alternative conformations resolved
        """
        # Group all related lines (ATOM, HETATM, ANISOU) by residue and atom
        atom_groups = defaultdict(list)
        other_lines = []
        
        i = 0
        while i < len(lines):
            line = lines[i]
            
            if line.startswith(('ATOM', 'HETATM')):
                # Parse atom record
                chain_id = line[21:22].strip()
                res_num = line[22:26].strip()
                res_name = line[17:20].strip()
                atom_name = line[12:16].strip()
                alt_loc = line[16:17].strip()
                
                # Create unique key for this atom position
                key = (chain_id, res_num, res_name, atom_name)
                
                # Collect this atom line and any following ANISOU line
                atom_block = [line]
                
                # Check if next line is ANISOU for this atom
                if (i + 1 < len(lines) and 
                    lines[i + 1].startswith('ANISOU') and
                    lines[i + 1][6:11] == line[6:11]):  # Same atom serial number
                    atom_block.append(lines[i + 1])
                    i += 1  # Skip the ANISOU line in main loop
                
                atom_groups[key].append((atom_block, alt_loc))
                
            else:
                other_lines.append(line)
            
            i += 1
        
        # Select best conformation for each atom
        cleaned_atom_blocks = []
        alt_conf_resolved = 0
        
        for key, atom_list in atom_groups.items():
            if len(atom_list) == 1:
                # No alternative conformations - just add the block
                cleaned_atom_blocks.extend(atom_list[0][0])
            else:
                # Multiple conformations - select highest occupancy
                best_block = self._select_best_conformation_block(atom_list)
                cleaned_atom_blocks.extend(best_block)
                alt_conf_resolved += len(atom_list) - 1
        
        if alt_conf_resolved > 0:
            self.logger.info(f"Resolved {alt_conf_resolved} alternative conformations")
        
        # Combine other lines with cleaned atom blocks
        # Keep original order by inserting atoms at first ATOM/HETATM position
        result_lines = []
        atom_lines_inserted = False
        
        for line in other_lines:
            if not atom_lines_inserted and line.startswith(('ATOM', 'HETATM')):
                result_lines.extend(cleaned_atom_blocks)
                atom_lines_inserted = True
            elif not line.startswith(('ATOM', 'HETATM', 'ANISOU')):
                result_lines.append(line)
        
        # If no ATOM/HETATM lines in other_lines, append at end
        if not atom_lines_inserted:
            result_lines.extend(cleaned_atom_blocks)
        
        return result_lines
    
    def _select_best_conformation_block(self, atom_list: List[Tuple[List[str], str]]) -> List[str]:
        """
        Select the best conformation block from alternative conformations
        
        Args:
            atom_list: List of (block_lines, alt_loc) tuples where block_lines contains ATOM and optionally ANISOU
            
        Returns:
            Best block of PDB lines with alt_loc cleared
        """
        # Parse occupancies and select highest
        best_block = atom_list[0][0]
        best_occupancy = 0.0
        best_alt_loc = atom_list[0][1]
        
        for block_lines, alt_loc in atom_list:
            # Get occupancy from the ATOM/HETATM line (first line in block)
            atom_line = block_lines[0]
            try:
                occupancy = float(atom_line[54:60])
                if occupancy > best_occupancy:
                    best_occupancy = occupancy
                    best_block = block_lines
                    best_alt_loc = alt_loc
            except (ValueError, IndexError):
                # If occupancy can't be parsed, keep first one
                continue
        
        # Clear alternative location indicator from all lines in the block
        cleaned_block = []
        for line in best_block:
            if best_alt_loc and line.startswith(('ATOM', 'HETATM', 'ANISOU')):
                # Clear alt_loc (position 16)
                cleaned_line = line[:16] + ' ' + line[17:]
                cleaned_block.append(cleaned_line)
            else:
                cleaned_block.append(line)
        
        return cleaned_block
    
    def _filter_records(self, lines: List[str], remove_waters: bool,
                       remove_ions: bool, remove_ligands: bool) -> List[str]:
        """
        Filter out unwanted records based on options
        
        Args:
            lines: List of PDB lines
            remove_waters: Remove water molecules
            remove_ions: Remove ions
            remove_ligands: Remove non-protein ligands
            
        Returns:
            Filtered list of lines
        """
        filtered_lines = []
        removed_counts = defaultdict(int)
        
        for line in lines:
            keep_line = True
            
            if line.startswith('HETATM'):
                res_name = line[17:20].strip()
                
                # Check for waters
                if remove_waters and res_name in ('HOH', 'WAT', 'H2O', 'TIP', 'SOL'):
                    keep_line = False
                    removed_counts['waters'] += 1
                
                # Check for ions
                elif remove_ions and res_name in self.common_ions:
                    keep_line = False
                    removed_counts['ions'] += 1
                
                # Check for ligands (everything else that's not protein)
                elif remove_ligands and res_name not in self.standard_residues:
                    # Don't remove if it's water or ion (handled above)
                    if not (res_name in ('HOH', 'WAT', 'H2O', 'TIP', 'SOL') or 
                           res_name in self.common_ions):
                        keep_line = False
                        removed_counts['ligands'] += 1
            
            if keep_line:
                filtered_lines.append(line)
        
        # Report what was removed
        for category, count in removed_counts.items():
            if count > 0:
                self.logger.info(f"Removed {count} {category}")
        
        return filtered_lines
    
    def _report_cleaning_stats(self, original_lines: List[str], 
                              cleaned_lines: List[str]) -> None:
        """Report statistics about the cleaning process"""
        
        original_atoms = sum(1 for line in original_lines 
                           if line.startswith(('ATOM', 'HETATM')))
        cleaned_atoms = sum(1 for line in cleaned_lines 
                          if line.startswith(('ATOM', 'HETATM')))
        
        self.logger.info(f"Original atoms: {original_atoms}")
        self.logger.info(f"Cleaned atoms: {cleaned_atoms}")
        self.logger.info(f"Atoms removed: {original_atoms - cleaned_atoms}")


def main():
    """Main function for command line usage"""
    parser = argparse.ArgumentParser(
        description="Clean PDB files for molecular docking preparation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic cleaning - resolve alt conformations only
  python pdb_cleaner.py input.pdb output.pdb
  
  # Remove waters and ions
  python pdb_cleaner.py input.pdb output.pdb --remove-waters --remove-ions
  
  # Remove everything except protein
  python pdb_cleaner.py input.pdb output.pdb --remove-waters --remove-ions --remove-ligands
        """
    )
    
    parser.add_argument('input_pdb', help='Input PDB file')
    parser.add_argument('output_pdb', help='Output cleaned PDB file')
    parser.add_argument('--remove-waters', action='store_true',
                       help='Remove water molecules (HOH, WAT, etc.)')
    parser.add_argument('--remove-ions', action='store_true',
                       help='Remove common ions (CA, MG, ZN, etc.)')
    parser.add_argument('--remove-ligands', action='store_true',
                       help='Remove non-protein ligands (keeps only standard amino acids)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Setup logging
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format='%(levelname)s: %(message)s'
    )
    
    # Validate input file
    if not Path(args.input_pdb).exists():
        print(f"Error: Input file {args.input_pdb} not found")
        sys.exit(1)
    
    # Create output directory if needed
    output_dir = Path(args.output_pdb).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Clean the PDB file
    cleaner = PDBCleaner()
    try:
        cleaner.clean_pdb(
            args.input_pdb,
            args.output_pdb,
            remove_waters=args.remove_waters,
            remove_ions=args.remove_ions,
            remove_ligands=args.remove_ligands
        )
        print(f"Successfully cleaned PDB file: {args.output_pdb}")
    except Exception as e:
        print(f"Error cleaning PDB file: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()