#!/usr/bin/env python3
"""
3D Structure Viewer Integration
Handles conversion of docking results for web visualization using 3Dmol.js
"""

import os
import json
import logging
from pathlib import Path
from typing import Dict, List, Any, Optional
import base64


class StructureViewer:
    """Handles 3D structure visualization for docking results"""
    
    def __init__(self, config_manager):
        """
        Initialize structure viewer
        
        Args:
            config_manager: Configuration manager instance
        """
        self.config = config_manager
        self.logger = logging.getLogger(__name__)
        
    def convert_pdbqt_to_pdb(self, pdbqt_file: str, output_pdb: str = None) -> str:
        """
        Convert PDBQT file to PDB format for 3Dmol.js
        
        Args:
            pdbqt_file: Path to PDBQT file
            output_pdb: Optional output PDB file path
            
        Returns:
            PDB content as string
        """
        try:
            with open(pdbqt_file, 'r') as f:
                pdbqt_content = f.read()
            
            # Convert PDBQT to PDB by removing AutoDock-specific lines
            pdb_lines = []
            for line in pdbqt_content.split('\n'):
                # Keep coordinate lines and essential PDB records
                if (line.startswith(('ATOM', 'HETATM', 'END', 'CONECT', 'MASTER', 'HEADER')) 
                    and not line.startswith(('ROOT', 'ENDROOT', 'BRANCH', 'ENDBRANCH', 'TORSDOF'))):
                    # Remove AutoDock specific columns (charges, atom types)
                    if line.startswith(('ATOM', 'HETATM')):
                        # Keep standard PDB format (first 66 characters)
                        pdb_line = line[:66].ljust(66) + '\n'
                        pdb_lines.append(pdb_line.rstrip())
                    else:
                        pdb_lines.append(line)
            
            pdb_content = '\n'.join(pdb_lines)
            
            # Write to file if requested
            if output_pdb:
                with open(output_pdb, 'w') as f:
                    f.write(pdb_content)
                self.logger.info(f"Converted {pdbqt_file} to {output_pdb}")
            
            return pdb_content
            
        except Exception as e:
            self.logger.error(f"Error converting {pdbqt_file} to PDB: {e}")
            return ""
    
    def prepare_receptor_for_viewing(self, receptor_pdbqt: str) -> str:
        """
        Prepare receptor structure for 3D viewing
        
        Args:
            receptor_pdbqt: Path to receptor PDBQT file
            
        Returns:
            PDB content as string
        """
        return self.convert_pdbqt_to_pdb(receptor_pdbqt)
    
    def prepare_docking_results_for_viewing(self, docking_data: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Prepare all docking results for 3D viewing
        
        Args:
            docking_data: List of docking result dictionaries
            
        Returns:
            Dictionary with structure data for visualization
        """
        if not docking_data:
            return {}
        
        output_dir = Path(self.config.get('docking.output_dir', 'docking_results'))
        
        # Find receptor file (should be in temp directory, but let's create a copy)
        receptor_files = list(output_dir.glob('*receptor*.pdbqt'))
        receptor_pdb_content = ""
        
        if not receptor_files:
            # Try to recreate receptor if needed
            self.logger.warning("No receptor file found in docking results")
        else:
            # Use the first receptor file found
            receptor_pdb_content = self.convert_pdbqt_to_pdb(str(receptor_files[0]))
        
        # Prepare ligand structures
        structure_data = {
            'receptor': receptor_pdb_content,
            'ligands': {},
            'binding_energies': {}
        }
        
        for result in docking_data:
            compound_id = result.get('compound_id', 0)
            compound_name = result.get('compound_name', f"Compound_{compound_id}")
            pose_file = result.get('pose_file', '')
            binding_energy = result.get('binding_energy', 0.0)
            
            if pose_file and Path(pose_file).exists():
                # Convert ligand pose to PDB
                ligand_pdb = self.convert_pdbqt_to_pdb(pose_file)
                structure_data['ligands'][compound_id] = {
                    'name': compound_name,
                    'pdb_content': ligand_pdb,
                    'binding_energy': binding_energy
                }
                structure_data['binding_energies'][compound_id] = binding_energy
        
        self.logger.info(f"Prepared {len(structure_data['ligands'])} structures for viewing")
        return structure_data
    
    def encode_structure_data(self, structure_data: Dict[str, Any]) -> str:
        """
        Encode structure data for embedding in HTML
        
        Args:
            structure_data: Structure data dictionary
            
        Returns:
            Base64 encoded JSON string
        """
        try:
            json_data = json.dumps(structure_data)
            encoded_data = base64.b64encode(json_data.encode()).decode()
            return encoded_data
        except Exception as e:
            self.logger.error(f"Error encoding structure data: {e}")
            return ""
    
    def generate_3dmol_html(self, structure_data: Dict[str, Any]) -> str:
        """
        Generate HTML for 3Dmol.js viewer
        
        Args:
            structure_data: Structure data dictionary
            
        Returns:
            HTML content for 3D viewer
        """
        # Encode data for JavaScript
        encoded_data = self.encode_structure_data(structure_data)
        
        html_content = f"""
        <div id="structure-viewer-container">
            <div id="structure-viewer-header">
                <h3>🧬 3D Structure Viewer</h3>
                <div id="structure-controls">
                    <select id="ligand-selector">
                        <option value="">Select compound...</option>
                    </select>
                    <button id="center-view">Center View</button>
                    <button id="toggle-surface">Toggle Surface</button>
                </div>
            </div>
            <div id="structure-viewer" style="width: 100%; flex: 1; min-height: 400px; background: white; border: 1px solid #ddd; border-radius: 4px; margin: 0 auto; display: flex; align-items: center; justify-content: center;"></div>
            <div id="structure-info">
                <div id="binding-energy-display"></div>
            </div>
        </div>
        
        <script>
            // Structure data (base64 encoded)
            const structureDataEncoded = '{encoded_data}';
            let structureData = {{}};
            
            // Decode structure data
            try {{
                const decodedData = atob(structureDataEncoded);
                structureData = JSON.parse(decodedData);
            }} catch (e) {{
                console.error('Error decoding structure data:', e);
            }}
            
            // Initialize 3Dmol viewer
            let viewer = null;
            let currentLigand = null;
            let showSurface = false;
            
            function initializeViewer() {{
                const element = $('#structure-viewer');
                const config = {{ backgroundColor: 'white' }};
                viewer = $3Dmol.createViewer(element, config);
                
                // Load receptor if available
                if (structureData.receptor) {{
                    viewer.addModel(structureData.receptor, 'pdb');
                    viewer.setStyle({{model: 0}}, {{
                        cartoon: {{ 
                            color: 'spectrum',
                            opacity: 0.8
                        }}
                    }});
                }}
                
                // Populate ligand selector
                const selector = document.getElementById('ligand-selector');
                Object.keys(structureData.ligands || {{}}).forEach(ligandId => {{
                    const ligand = structureData.ligands[ligandId];
                    const option = document.createElement('option');
                    option.value = ligandId;
                    option.textContent = `${{ligand.name}} (${{ligand.binding_energy.toFixed(1)}} kcal/mol)`;
                    selector.appendChild(option);
                }});
                
                viewer.zoomTo();
                viewer.render();
            }}
            
            function showLigand(ligandId) {{
                if (!viewer || !structureData.ligands[ligandId]) return;
                
                // Remove previous ligand
                if (currentLigand !== null) {{
                    viewer.removeModel(viewer.getModel(1));
                }}
                
                // Add new ligand
                const ligand = structureData.ligands[ligandId];
                viewer.addModel(ligand.pdb_content, 'pdb');
                
                // Style the ligand - green carbons with heteroatom coloring
                viewer.setStyle({{model: 1}}, {{
                    stick: {{ 
                        colorscheme: {{ 
                            'C': 'green',
                            'N': 'blue', 
                            'O': 'red', 
                            'S': 'yellow',
                            'P': 'orange',
                            'F': 'lightblue',
                            'Cl': 'lightgreen',
                            'Br': 'darkred',
                            'I': 'purple'
                        }}, 
                        radius: 0.25 
                    }}
                }});
                
                // Update binding energy display
                const energyDisplay = document.getElementById('binding-energy-display');
                energyDisplay.innerHTML = `
                    <strong>${{ligand.name}}</strong><br>
                    Binding Energy: <span style="color: #4CAF50;">${{ligand.binding_energy.toFixed(1)}} kcal/mol</span>
                `;
                
                currentLigand = ligandId;
                viewer.zoomTo();
                viewer.render();
            }}
            
            function centerView() {{
                if (viewer) {{
                    viewer.zoomTo();
                    viewer.render();
                }}
            }}
            
            function toggleSurface() {{
                if (!viewer) return;
                
                showSurface = !showSurface;
                
                if (showSurface) {{
                    viewer.addSurface($3Dmol.VDW, {{
                        opacity: 0.4,
                        color: 'lightgray'
                    }}, {{model: 0}});
                    document.getElementById('toggle-surface').textContent = 'Hide Surface';
                }} else {{
                    viewer.setSurfaceMaterialStyle(0, {{opacity: 0}});
                    document.getElementById('toggle-surface').textContent = 'Show Surface';
                }}
                
                viewer.render();
            }}
            
            // Event listeners
            document.getElementById('ligand-selector').addEventListener('change', function(e) {{
                if (e.target.value) {{
                    showLigand(e.target.value);
                }}
            }});
            
            document.getElementById('center-view').addEventListener('click', centerView);
            document.getElementById('toggle-surface').addEventListener('click', toggleSurface);
            
            // Initialize when 3Dmol is ready
            if (typeof $3Dmol !== 'undefined') {{
                initializeViewer();
            }} else {{
                document.addEventListener('DOMContentLoaded', function() {{
                    setTimeout(initializeViewer, 1000);
                }});
            }}
        </script>
        
        <style>
            #structure-viewer-container {{
                border: 1px solid #ddd;
                border-radius: 8px;
                padding: 15px;
                background: white;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                width: 100%;
                height: 100%;
                display: flex;
                flex-direction: column;
                align-items: center;
            }}
            
            #structure-viewer-header {{
                display: flex;
                justify-content: space-between;
                align-items: center;
                margin-bottom: 10px;
            }}
            
            #structure-controls {{
                display: flex;
                gap: 10px;
                align-items: center;
            }}
            
            #structure-viewer {{
                position: relative;
                margin: 0 auto !important;
                display: flex !important;
                align-items: center !important;
                justify-content: center !important;
            }}
            
            #structure-viewer canvas {{
                margin: 0 auto !important;
                display: block !important;
            }}
            
            #ligand-selector {{
                padding: 5px 10px;
                border: 1px solid #ccc;
                border-radius: 4px;
                min-width: 200px;
            }}
            
            #structure-controls button {{
                padding: 5px 12px;
                background: #2196F3;
                color: white;
                border: none;
                border-radius: 4px;
                cursor: pointer;
                font-size: 12px;
            }}
            
            #structure-controls button:hover {{
                background: #1976D2;
            }}
            
            #structure-info {{
                margin-top: 10px;
                padding: 10px;
                background: #f5f5f5;
                border-radius: 4px;
                font-family: monospace;
            }}
            
            #binding-energy-display {{
                color: #333;
            }}
        </style>
        """
        
        return html_content