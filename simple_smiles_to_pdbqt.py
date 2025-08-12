#!/usr/bin/env python3
"""
Simple SMILES to PDBQT converter using direct PDB format
"""
import os
import subprocess
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

def smiles_to_pdbqt_simple(smiles, output_pdbqt, name="UNK"):
    """
    Convert SMILES to PDBQT using simple PDB intermediate
    """
    print(f"Converting: {smiles} -> {name}")
    
    # MGLTools paths
    mgltools_path = "/Users/enmingxing/Projects/mol_view_dashboard/packages/mgltools_1.5.7_MacOS-X/installed"
    mgl_python = f"{mgltools_path}/bin/python"
    prepare_ligand_script = f"{mgltools_path}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
    
    # Create working directory
    work_dir = Path("temp_ligand_prep")
    work_dir.mkdir(exist_ok=True)
    
    try:
        # Step 1: Create molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"❌ Invalid SMILES: {smiles}")
            return False
        
        # Step 2: Generate 3D structure
        mol_h = Chem.AddHs(mol)
        embed_result = AllChem.EmbedMolecule(mol_h, randomSeed=42)
        if embed_result == -1:
            embed_result = AllChem.EmbedMolecule(mol_h, randomSeed=42, useRandomCoords=True)
            if embed_result == -1:
                print(f"❌ Failed to generate 3D coordinates")
                return False
        
        try:
            AllChem.MMFFOptimizeMolecule(mol_h)
        except:
            pass  # Continue with unoptimized structure
        
        # Step 3: Write PDB file
        pdb_file = work_dir / f"{name}.pdb"
        pdb_block = Chem.MolToPDBBlock(mol_h)
        
        with open(pdb_file, 'w') as f:
            f.write(pdb_block)
        
        if not pdb_file.exists() or pdb_file.stat().st_size == 0:
            print(f"❌ Failed to create PDB file")
            return False
        
        print(f"✅ Created PDB: {pdb_file} ({pdb_file.stat().st_size} bytes)")
        
        # Step 4: Run MGLTools from the working directory
        original_dir = Path.cwd()
        os.chdir(work_dir)
        
        try:
            cmd = [
                mgl_python,
                prepare_ligand_script,
                '-l', f"{name}.pdb",  # Use relative path
                '-o', f"{name}.pdbqt",
                '-A', 'hydrogens'
            ]
            
            print(f"Running in {work_dir}: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
            
            pdbqt_file = Path(f"{name}.pdbqt")
            
            if result.returncode == 0 and pdbqt_file.exists():
                # Move PDBQT to final location
                final_pdbqt = original_dir / output_pdbqt
                final_pdbqt.parent.mkdir(exist_ok=True)
                pdbqt_file.rename(final_pdbqt)
                print(f"✅ Created PDBQT: {final_pdbqt}")
                return True
            else:
                print(f"❌ MGLTools failed:")
                print(f"Return code: {result.returncode}")
                print(f"STDERR: {result.stderr}")
                print(f"STDOUT: {result.stdout}")
                return False
                
        finally:
            os.chdir(original_dir)
            
    except Exception as e:
        print(f"❌ Error: {e}")
        return False

if __name__ == "__main__":
    # Test
    test_data = [
        ("CCO", "Ethanol"),
        ("c1ccccc1", "Benzene")
    ]
    
    output_dir = Path("test_pdbqt_simple")
    output_dir.mkdir(exist_ok=True)
    
    for smiles, name in test_data:
        output_file = output_dir / f"{name}.pdbqt"
        smiles_to_pdbqt_simple(smiles, output_file, name)
        print()