#!/usr/bin/env python3
"""
Convert SMILES strings directly to docking-ready PDBQT files
"""
import sys
import subprocess
import tempfile
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

def smiles_to_pdbqt(smiles, output_pdbqt, name="UNK"):
    """
    Convert SMILES string to PDBQT file ready for docking
    
    Args:
        smiles: SMILES string
        output_pdbqt: Output PDBQT file path
        name: Compound name
    
    Returns:
        True if successful, False otherwise
    """
    print(f"Converting SMILES: {smiles} -> {output_pdbqt}")
    
    # MGLTools paths
    mgltools_path = "/Users/enmingxing/Projects/mol_view_dashboard/packages/mgltools_1.5.7_MacOS-X/installed"
    mgl_python = f"{mgltools_path}/bin/python"
    prepare_ligand_script = f"{mgltools_path}/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py"
    
    try:
        # Step 1: Create molecule from SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"❌ Failed to create molecule from SMILES: {smiles}")
            return False
        
        # Step 2: Add hydrogens and generate 3D coordinates
        mol_h = Chem.AddHs(mol)
        embed_result = AllChem.EmbedMolecule(mol_h, randomSeed=42)
        if embed_result == -1:
            # Try with random coordinates
            embed_result = AllChem.EmbedMolecule(mol_h, randomSeed=42, useRandomCoords=True)
            if embed_result == -1:
                print(f"❌ Failed to generate 3D coordinates for: {smiles}")
                return False
        
        # Step 3: Optimize geometry
        try:
            AllChem.MMFFOptimizeMolecule(mol_h)
            print(f"✅ 3D structure optimized")
        except Exception as e:
            print(f"⚠️ MMFF optimization failed: {e}, continuing anyway")
        
        # Step 4: Write to SDF file (MGLTools prefers SDF over PDB for small molecules)
        temp_dir = tempfile.mkdtemp()
        sdf_file = Path(temp_dir) / f"{name}.sdf"
        
        writer = Chem.SDWriter(str(sdf_file))
        mol_h.SetProp("_Name", name)
        writer.write(mol_h)
        writer.close()
        
        # Verify SDF file was created
        if not sdf_file.exists() or sdf_file.stat().st_size == 0:
            print(f"❌ Failed to create SDF file")
            return False
        
        print(f"✅ Created SDF file: {sdf_file} ({sdf_file.stat().st_size} bytes)")
        
        # Step 5: Convert SDF to PDBQT using MGLTools
        cmd = [
            mgl_python,
            prepare_ligand_script,
            '-l', str(sdf_file),
            '-o', str(output_pdbqt),
            '-A', 'hydrogens',  # Add hydrogens
            '-U', 'waters'      # Remove waters
        ]
        
        print(f"Running: {' '.join(cmd)}")
        print(f"SDF file exists before MGLTools: {sdf_file.exists()}")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        # Clean up temp directory
        import shutil
        shutil.rmtree(temp_dir)
        
        if result.returncode == 0 and Path(output_pdbqt).exists():
            print(f"✅ Successfully created PDBQT: {output_pdbqt}")
            return True
        else:
            print(f"❌ MGLTools failed: {result.stderr}")
            return False
            
    except Exception as e:
        print(f"❌ Error: {e}")
        return False

if __name__ == "__main__":
    # Test with sample SMILES
    test_smiles = [
        ("CCO", "Ethanol"),
        ("c1ccccc1", "Benzene"),
        ("CC(=O)O", "Acetic_acid"),
        ("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "Ibuprofen")
    ]
    
    output_dir = Path("test_pdbqt_output")
    output_dir.mkdir(exist_ok=True)
    
    success_count = 0
    for smiles, name in test_smiles:
        output_file = output_dir / f"{name}.pdbqt"
        if smiles_to_pdbqt(smiles, output_file, name):
            success_count += 1
        print()
    
    print(f"Successfully converted {success_count}/{len(test_smiles)} molecules")
    
    if success_count > 0:
        print(f"\nPDBQT files created in: {output_dir}")
        for pdbqt_file in output_dir.glob("*.pdbqt"):
            print(f"  {pdbqt_file.name} ({pdbqt_file.stat().st_size} bytes)")