#!/usr/bin/env python3
"""
Debug script to test ligand preparation
"""
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import tempfile
from pathlib import Path

# Test with a simple molecule
test_smiles = ["CCO", "c1ccccc1", "CC(=O)O"]  # ethanol, benzene, acetic acid
test_names = ["Ethanol", "Benzene", "Acetic_acid"]

print("Testing ligand preparation...")

for i, (smiles, name) in enumerate(zip(test_smiles, test_names)):
    print(f"\n=== Testing {name} (SMILES: {smiles}) ===")
    
    # Step 1: Create molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"❌ Failed to create molecule from SMILES")
        continue
    
    # Step 2: Add hydrogens
    mol_h = Chem.AddHs(mol)
    print(f"✅ Added hydrogens: {mol_h.GetNumAtoms()} atoms")
    
    # Step 3: Generate 3D coordinates
    embed_result = AllChem.EmbedMolecule(mol_h, randomSeed=42)
    if embed_result == -1:
        print(f"⚠️ First embed failed, trying with random coords...")
        embed_result = AllChem.EmbedMolecule(mol_h, randomSeed=42, useRandomCoords=True)
        if embed_result == -1:
            print(f"❌ Failed to embed 3D coordinates")
            continue
    
    print(f"✅ Generated 3D coordinates: {mol_h.GetNumConformers()} conformers")
    
    # Step 4: Optimize
    try:
        AllChem.MMFFOptimizeMolecule(mol_h)
        print(f"✅ MMFF optimization successful")
    except Exception as e:
        print(f"⚠️ MMFF optimization failed: {e}")
    
    # Step 5: Generate PDB
    try:
        pdb_block = Chem.MolToPDBBlock(mol_h)
        if pdb_block and pdb_block.strip():
            print(f"✅ PDB block generated: {len(pdb_block)} characters")
            
            # Write to file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
                f.write(pdb_block)
                pdb_file = f.name
            
            pdb_path = Path(pdb_file)
            if pdb_path.exists() and pdb_path.stat().st_size > 0:
                print(f"✅ PDB file created: {pdb_file} ({pdb_path.stat().st_size} bytes)")
                
                # Show first few lines
                with open(pdb_file) as f:
                    lines = f.readlines()[:5]
                    print("📄 PDB content preview:")
                    for line in lines:
                        print(f"   {line.strip()}")
            else:
                print(f"❌ PDB file creation failed")
        else:
            print(f"❌ Empty PDB block generated")
    except Exception as e:
        print(f"❌ PDB generation failed: {e}")

print("\nDone!")