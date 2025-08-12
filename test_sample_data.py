#!/usr/bin/env python3
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import tempfile
from pathlib import Path

# Load the first few compounds from the sample file
df = pd.read_csv('/Users/enmingxing/Projects/mol_view_dashboard/examples/sample_compounds.csv')
print('Sample data:')
print(df.head(3)[['SMILES', 'Name']])

# Test the first compound
smiles = df.iloc[0]['SMILES']
name = df.iloc[0]['Name']

print(f'\nTesting {name}: {smiles}')

mol = Chem.MolFromSmiles(smiles)
if mol:
    mol_h = Chem.AddHs(mol)
    embed_result = AllChem.EmbedMolecule(mol_h, randomSeed=42)
    print(f'Embed result: {embed_result}')
    
    if embed_result != -1:
        AllChem.MMFFOptimizeMolecule(mol_h)
        pdb_block = Chem.MolToPDBBlock(mol_h)
        print(f'PDB block length: {len(pdb_block) if pdb_block else 0}')
        
        # Test creating file
        temp_dir = tempfile.mkdtemp()
        ligand_pdb = Path(temp_dir) / 'test_ligand.pdb'
        with open(ligand_pdb, 'w') as f:
            f.write(pdb_block)
        
        print(f'Created: {ligand_pdb}')
        print(f'Exists: {ligand_pdb.exists()}')
        print(f'Size: {ligand_pdb.stat().st_size if ligand_pdb.exists() else 0}')
    else:
        print('❌ Embed failed')
else:
    print('❌ Mol creation failed')