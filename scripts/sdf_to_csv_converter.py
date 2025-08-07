#!/usr/bin/env python3
"""
SDF to CSV Converter
Converts SDF files to CSV format with canonicalized SMILES and molecular properties
- Reads single SDF file or multiple SDF files from a directory
- Extracts molecular structures and properties
- Generates canonicalized SMILES strings
- Outputs CSV with SMILES column and all available properties
"""

import pandas as pd
import numpy as np
import argparse
import sys
import os
from pathlib import Path

# RDKit imports
try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
except ImportError as e:
    print(f"Error: RDKit is required. Install with: conda install -c conda-forge rdkit")
    print(f"Specific error: {e}")
    sys.exit(1)


def read_sdf_file(sdf_file):
    """Read molecules from an SDF file and extract properties"""
    print(f"Reading SDF file: {sdf_file}")
    
    try:
        # Read molecules from SDF
        supplier = Chem.SDMolSupplier(str(sdf_file))
        molecules_data = []
        sdf_basename = Path(sdf_file).stem  # Get filename without extension
        
        for i, mol in enumerate(supplier):
            if mol is None:
                print(f"  Warning: Could not parse molecule {i+1}")
                continue
            
            # Generate canonicalized SMILES
            try:
                canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
            except Exception as e:
                print(f"  Warning: Could not generate SMILES for molecule {i+1}: {e}")
                continue
            
            # Extract all properties from the molecule
            mol_data = {
                'SMILES': canonical_smiles,
                'Name': f"{sdf_basename}_{i+1}",  # SDF basename + number starting from 1
                'Molecule_Index': i + 1
            }
            
            # Get all properties from the SDF
            prop_names = mol.GetPropNames()
            for prop_name in prop_names:
                try:
                    prop_value = mol.GetProp(prop_name)
                    # Try to convert to numeric if possible
                    try:
                        prop_value = float(prop_value)
                    except (ValueError, TypeError):
                        # Keep as string if not numeric
                        pass
                    mol_data[prop_name] = prop_value
                except Exception as e:
                    print(f"  Warning: Could not extract property '{prop_name}' for molecule {i+1}: {e}")
                    mol_data[prop_name] = None
            
            molecules_data.append(mol_data)
        
        print(f"  Successfully processed {len(molecules_data)} molecules")
        return molecules_data
        
    except Exception as e:
        print(f"Error reading SDF file {sdf_file}: {e}")
        return []


def process_single_sdf(sdf_file, output_file=None):
    """Process a single SDF file"""
    sdf_path = Path(sdf_file)
    
    if not sdf_path.exists():
        raise FileNotFoundError(f"SDF file not found: {sdf_file}")
    
    if not output_file:
        output_file = sdf_path.with_suffix('.csv')
    
    # Read molecules from SDF
    molecules_data = read_sdf_file(sdf_path)
    
    if not molecules_data:
        print(f"No valid molecules found in {sdf_file}")
        return
    
    # Create DataFrame
    df = pd.DataFrame(molecules_data)
    
    # Reorder columns to put SMILES and Name first
    columns = ['SMILES', 'Name'] + [col for col in df.columns if col not in ['SMILES', 'Name']]
    df = df[columns]
    
    # Save to CSV
    df.to_csv(output_file, index=False)
    
    print(f"✅ Converted {len(df)} molecules from {sdf_file} to {output_file}")
    print(f"📊 Columns: {list(df.columns)}")
    
    return df


def process_multiple_sdf(input_directory, output_file=None, merge_files=True):
    """Process multiple SDF files from a directory"""
    input_path = Path(input_directory)
    
    if not input_path.exists() or not input_path.is_dir():
        raise NotADirectoryError(f"Directory not found: {input_directory}")
    
    # Find all SDF files
    sdf_files = list(input_path.glob("*.sdf")) + list(input_path.glob("*.SDF"))
    
    if not sdf_files:
        print(f"No SDF files found in {input_directory}")
        return
    
    print(f"Found {len(sdf_files)} SDF files:")
    for sdf_file in sdf_files:
        print(f"  - {sdf_file.name}")
    
    all_molecules_data = []
    
    for sdf_file in sdf_files:
        molecules_data = read_sdf_file(sdf_file)
        
        # Add source file information
        for mol_data in molecules_data:
            mol_data['Source_File'] = sdf_file.stem
        
        all_molecules_data.extend(molecules_data)
        
        # If not merging, save individual CSV files
        if not merge_files:
            individual_output = sdf_file.with_suffix('.csv')
            df_individual = pd.DataFrame(molecules_data)
            
            columns = ['SMILES', 'Name'] + [col for col in df_individual.columns if col not in ['SMILES', 'Name']]
            df_individual = df_individual[columns]
            df_individual.to_csv(individual_output, index=False)
            print(f"  → Saved individual file: {individual_output}")
    
    if merge_files and all_molecules_data:
        # Create merged DataFrame
        df_merged = pd.DataFrame(all_molecules_data)
        
        # Reorder columns
        columns = ['SMILES', 'Name', 'Source_File'] + [col for col in df_merged.columns if col not in ['SMILES', 'Name', 'Source_File']]
        df_merged = df_merged[columns]
        
        # Save merged file
        if not output_file:
            output_file = input_path / "merged_molecules.csv"
        
        df_merged.to_csv(output_file, index=False)
        
        print(f"\n✅ Merged {len(df_merged)} molecules from {len(sdf_files)} files to {output_file}")
        print(f"📊 Columns: {list(df_merged.columns)}")
        print(f"📁 Files processed: {[f.stem for f in sdf_files]}")
        
        return df_merged


def main():
    parser = argparse.ArgumentParser(description="Convert SDF files to CSV with canonicalized SMILES")
    parser.add_argument('input', help='Input SDF file or directory containing SDF files')
    parser.add_argument('--output', '-o', help='Output CSV file (default: same name with .csv extension)')
    parser.add_argument('--individual', action='store_true', 
                       help='When processing directory, save individual CSV files instead of merging')
    parser.add_argument('--no-merge', action='store_true',
                       help='When processing directory, only save individual files (no merged file)')
    
    args = parser.parse_args()
    
    try:
        input_path = Path(args.input)
        
        if input_path.is_file():
            # Process single SDF file
            if not args.input.lower().endswith(('.sdf', '.SDF')):
                print("Warning: Input file does not have .sdf extension")
            
            process_single_sdf(args.input, args.output)
            
        elif input_path.is_dir():
            # Process multiple SDF files from directory
            merge_files = not args.no_merge
            
            if args.individual:
                # Save both individual and merged files
                process_multiple_sdf(args.input, args.output, merge_files=True)
            else:
                # Default behavior based on --no-merge flag
                process_multiple_sdf(args.input, args.output, merge_files=merge_files)
        else:
            raise FileNotFoundError(f"Input path does not exist: {args.input}")
    
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    print(f"\n🎉 Conversion completed successfully!")
    print(f"📋 Usage tips:")
    print(f"   • SMILES column contains canonicalized SMILES strings")
    print(f"   • All SDF properties are preserved in the CSV")
    print(f"   • Use the CSV with molecular_property_visualizer.py for analysis")


if __name__ == "__main__":
    main()