#!/usr/bin/env python3
"""
Comprehensive Molecular Analysis Dashboard
Creates multiple interactive scatter plots with t-SNE visualization
- Property plots: MW vs LogP, TPSA vs QED, etc.
- t-SNE molecular space visualization
- Synchronized highlighting across all plots
- Chemical structure display on hover
"""

import pandas as pd
import numpy as np
import base64
from io import BytesIO
import json
import argparse
import sys
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler

# RDKit imports
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, Lipinski, QED, Draw
    from rdkit.Chem import rdMolDescriptors, rdFingerprintGenerator
    from rdkit.Contrib.SA_Score import sascorer
except ImportError as e:
    print(f"Error: RDKit is required. Install with: conda install -c conda-forge rdkit")
    print(f"Specific error: {e}")
    sys.exit(1)


def mol_to_base64_png(mol, size=(300, 300)):
    """Convert molecule to base64 encoded PNG"""
    if mol is None:
        return ""
    
    try:
        img = Draw.MolToImage(mol, size=size)
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode()
        return f"data:image/png;base64,{img_str}"
    except Exception as e:
        print(f"Error generating structure: {e}")
        return ""


def load_and_process_data(input_file, properties_from_input=None):
    """Load CSV and process molecular data"""
    print(f"Loading data from: {input_file}")
    
    df = pd.read_csv(input_file)
    print(f"Original data shape: {df.shape}")
    print(f"Columns: {list(df.columns)}")
    
    # Check required columns
    required_cols = ['SMILES']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    # Handle optional columns with defaults
    if 'Name' not in df.columns:
        df['Name'] = [f"Compound_{i+1}" for i in range(len(df))]
    
    # Map Name to title for internal consistency
    df['title'] = df['Name']
    
    # Handle custom property columns
    if properties_from_input:
        for i, prop in enumerate(properties_from_input):
            if prop not in df.columns:
                print(f"Warning: '{prop}' column not found in CSV")
                df[f'custom_property_{i}'] = np.nan
            else:
                df[f'custom_property_{i}'] = df[prop]
                print(f"Using '{prop}' as custom property {i}")
    # If no custom properties specified, don't look for docking score either
    
    print("Creating molecule objects from SMILES...")
    df['Mol'] = df['SMILES'].apply(lambda x: Chem.MolFromSmiles(x) if pd.notna(x) else None)
    
    # Show failed SMILES
    failed_smiles = df[df['Mol'].isna()]
    if len(failed_smiles) > 0:
        print(f"Failed to parse {len(failed_smiles)} SMILES:")
        for idx, row in failed_smiles.head().iterrows():
            print(f"  {idx}: {row['SMILES']}")
    
    # Filter valid molecules
    df = df[df['Mol'].notna()].reset_index(drop=True)
    print(f"Valid molecules: {len(df)}")
    
    return df


def calculate_comprehensive_descriptors(df):
    """Calculate comprehensive molecular descriptors"""
    print(f"Calculating molecular descriptors for {len(df)} compounds...")
    
    descriptors = []
    
    for i, mol in enumerate(df['Mol']):
        if i % 100 == 0:
            print(f"  Processing molecule {i+1}/{len(df)}")
        
        try:
            desc = {
                # Basic properties
                'MW': Descriptors.MolWt(mol),
                'LogP': Crippen.MolLogP(mol),
                'TPSA': rdMolDescriptors.CalcTPSA(mol),
                'HBA': rdMolDescriptors.CalcNumHBA(mol),
                'HBD': rdMolDescriptors.CalcNumHBD(mol),
                'RotBonds': rdMolDescriptors.CalcNumRotatableBonds(mol),
                'AromaticRings': rdMolDescriptors.CalcNumAromaticRings(mol),
                'HeavyAtoms': mol.GetNumHeavyAtoms(),
                
                # Lipinski descriptors
                'MolMR': Crippen.MolMR(mol),
                
                # Drug-likeness
                'QED': QED.qed(mol),
                'SAscore': sascorer.calculateScore(mol),
            }
            
        except Exception as e:
            print(f"Error calculating descriptors for molecule {i}: {e}")
            desc = {key: np.nan for key in [
                'MW', 'LogP', 'TPSA', 'HBA', 'HBD', 'RotBonds', 'AromaticRings', 
                'HeavyAtoms', 'MolMR', 'QED', 'SAscore'
            ]}
        
        descriptors.append(desc)
    
    # Add descriptors to dataframe
    desc_df = pd.DataFrame(descriptors)
    for col in desc_df.columns:
        df[col] = desc_df[col]
    
    # Remove rows with critical missing values
    critical_cols = ['MW', 'LogP', 'TPSA', 'QED']
    df = df.dropna(subset=critical_cols).reset_index(drop=True)
    
    print(f"Final dataset: {len(df)} compounds with complete descriptors")
    print("Descriptor summary:")
    print(df[['MW', 'LogP', 'TPSA', 'QED', 'SAscore', 'r_i_docking_score']].describe())
    
    return df


def calculate_tsne_coordinates(df):
    """Calculate t-SNE coordinates using molecular fingerprints"""
    print("Calculating t-SNE coordinates...")
    
    # Generate molecular fingerprints
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fingerprints = []
    
    for mol in df['Mol']:
        fp = generator.GetFingerprint(mol)
        fp_array = np.zeros((1,))
        Chem.DataStructs.ConvertToNumpyArray(fp, fp_array)
        fingerprints.append(fp_array)
    
    fp_matrix = np.array(fingerprints)
    print(f"Fingerprint matrix shape: {fp_matrix.shape}")
    
    # Standardize features for t-SNE
    scaler = StandardScaler()
    fp_scaled = scaler.fit_transform(fp_matrix)
    
    # Run t-SNE
    tsne = TSNE(n_components=2, random_state=42, perplexity=min(30, len(df)//4))
    tsne_coords = tsne.fit_transform(fp_scaled)
    
    df['tSNE_X'] = tsne_coords[:, 0]
    df['tSNE_Y'] = tsne_coords[:, 1]
    
    print("t-SNE calculation completed")
    return df


def create_comprehensive_dashboard(df, output_file="molecular_dashboard.html", properties_from_input=None):
    """Create comprehensive molecular analysis dashboard"""
    
    print(f"Creating comprehensive dashboard with {len(df)} compounds")
    
    # Generate structure images and prepare data
    print("Generating molecular structures...")
    data_points = []
    
    for idx, row in df.iterrows():
        if idx % 50 == 0:
            print(f"  Generating structure {idx+1}/{len(df)}")
        
        img_base64 = mol_to_base64_png(row['Mol'])
        
        # Convert all numeric values to float for JSON serialization
        data_point = {
            'id': int(idx),
            'title': str(row['title']),
            'smiles': str(row['SMILES']),
            'mw': float(row['MW']),
            'logp': float(row['LogP']),
            'tpsa': float(row['TPSA']),
            'hba': int(row['HBA']),
            'hbd': int(row['HBD']),
            'qed': float(row['QED']),
            'sascore': float(row['SAscore']),
            'rotbonds': int(row['RotBonds']),
            'tsne_x': float(row['tSNE_X']),
            'tsne_y': float(row['tSNE_Y']),
            'image': img_base64
        }
        
        # Add custom properties if provided
        if properties_from_input:
            for i, prop in enumerate(properties_from_input):
                prop_col = f'custom_property_{i}'
                if prop_col in row and pd.notna(row[prop_col]):
                    data_point[prop_col] = float(row[prop_col])
        data_points.append(data_point)
    
    # Convert to JSON
    data_json = json.dumps(data_points, indent=2)
    
    # Generate dynamic plot configurations based on properties
    plot_configs = []
    
    # Add default plots
    plot_configs.extend([
        {
            'id': 'plot1',
            'title': 'Molecular Weight vs LogP',
            'x_prop': 'mw',
            'y_prop': 'logp', 
            'color_prop': None,
            'x_label': 'Molecular Weight (Da)',
            'y_label': 'LogP',
            'color_label': None
        },
        {
            'id': 'plot2', 
            'title': 'TPSA vs QED',
            'x_prop': 'tpsa',
            'y_prop': 'qed',
            'color_prop': 'sascore',
            'x_label': 'TPSA (Ų)',
            'y_label': 'QED',
            'color_label': 'SA Score'
        },
        {
            'id': 'plot3',
            'title': 'Drug-likeness: HBA vs HBD', 
            'x_prop': 'hba',
            'y_prop': 'hbd',
            'color_prop': 'qed',
            'x_label': 'H-Bond Acceptors',
            'y_label': 'H-Bond Donors', 
            'color_label': 'QED'
        },
        {
            'id': 'plot4',
            'title': 't-SNE Molecular Space',
            'x_prop': 'tsne_x',
            'y_prop': 'tsne_y',
            'color_prop': None, 
            'x_label': 't-SNE Dimension 1',
            'y_label': 't-SNE Dimension 2',
            'color_label': None
        }
    ])
    
    # Add custom property plots
    if properties_from_input:
        num_props = len(properties_from_input)
        
        if num_props % 2 == 1:  # Odd number - use first as color, rest as scatter plots
            # Use first property as color for default plots
            color_prop = 'custom_property_0'
            color_label = properties_from_input[0]
            plot_configs[0]['color_prop'] = color_prop
            plot_configs[0]['color_label'] = color_label
            plot_configs[3]['color_prop'] = color_prop
            plot_configs[3]['color_label'] = color_label
            
            # Create scatter plots for remaining properties (pairs)
            for i in range(1, num_props, 2):
                if i + 1 < num_props:
                    plot_id = f'custom_plot_{i}'
                    x_prop = f'custom_property_{i}'
                    y_prop = f'custom_property_{i+1}'
                    x_label = properties_from_input[i]
                    y_label = properties_from_input[i+1]
                    
                    plot_configs.append({
                        'id': plot_id,
                        'title': f'{x_label} vs {y_label}',
                        'x_prop': x_prop,
                        'y_prop': y_prop,
                        'color_prop': color_prop,
                        'x_label': x_label,
                        'y_label': y_label,
                        'color_label': color_label
                    })
        else:  # Even number - create scatter plots for all pairs
            for i in range(0, num_props, 2):
                if i + 1 < num_props:
                    plot_id = f'custom_plot_{i}'
                    x_prop = f'custom_property_{i}'
                    y_prop = f'custom_property_{i+1}'
                    x_label = properties_from_input[i]
                    y_label = properties_from_input[i+1]
                    
                    plot_configs.append({
                        'id': plot_id,
                        'title': f'{x_label} vs {y_label}',
                        'x_prop': x_prop,
                        'y_prop': y_prop,
                        'color_prop': None,
                        'x_label': x_label,
                        'y_label': y_label,
                        'color_label': None
                    })
    
    # Generate custom plots HTML
    custom_plots_html = ""
    if len(plot_configs) > 4:
        custom_plots_html = "\n"
        # Group custom plots by rows of 2
        custom_configs = plot_configs[4:]
        for i in range(0, len(custom_configs), 2):
            custom_plots_html += '    <div class="custom-plots">\n'
            
            # First plot in the row
            config = custom_configs[i]
            custom_plots_html += f'        <div class="plot-container">\n'
            custom_plots_html += f'            <div class="plot-title">{config["title"]}</div>\n'
            custom_plots_html += f'            <div id="{config["id"]}"></div>\n'
            custom_plots_html += f'        </div>\n'
            
            # Second plot in the row (if exists)
            if i + 1 < len(custom_configs):
                config = custom_configs[i + 1]
                custom_plots_html += f'        <div class="plot-container">\n'
                custom_plots_html += f'            <div class="plot-title">{config["title"]}</div>\n'
                custom_plots_html += f'            <div id="{config["id"]}"></div>\n'
                custom_plots_html += f'        </div>\n'
            else:
                # Empty spacer if odd number of plots
                custom_plots_html += '        <div class="plot-container" style="visibility: hidden;"></div>\n'
            
            # Spacer for structure panel column
            custom_plots_html += '        <div class="spacer"></div>\n'
            custom_plots_html += '    </div>\n'
    
    plot_configs_json = json.dumps(plot_configs)
    
    # Create HTML dashboard
    html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Comprehensive Molecular Analysis Dashboard</title>
    <script src="https://d3js.org/d3.v7.min.js"></script>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 0;
            background: #f5f5f5;
            min-height: 100vh;
        }}
        
        .header {{
            text-align: center;
            padding: 20px;
            color: #333;
            background: white;
            border-bottom: 1px solid #ddd;
        }}
        
        .header h1 {{
            margin: 0;
            font-size: 2em;
            font-weight: normal;
        }}
        
        .stats {{
            background: #f8f9fa;
            padding: 15px;
            border: 1px solid #ddd;
            margin: 20px auto;
            max-width: 800px;
            color: #666;
            text-align: center;
        }}
        
        .dashboard {{
            display: grid;
            grid-template-columns: 1fr 1fr 400px;
            gap: 20px;
            max-width: 1800px;
            margin: 0 auto;
            padding: 20px;
        }}
        
        .custom-plots {{
            display: grid;
            grid-template-columns: 1fr 1fr 400px;
            gap: 20px;
            max-width: 1800px;
            margin: 0 auto;
            padding: 0 20px 20px 20px;
        }}
        
        .spacer {{
            /* Empty spacer for the structure panel column */
        }}
        
        .plot-container {{
            background: white;
            border: 1px solid #ddd;
            padding: 20px;
            position: relative;
        }}
        
        .plot-title {{
            text-align: center;
            margin-bottom: 15px;
            font-size: 14px;
            font-weight: normal;
            color: #333;
        }}
        
        .axis {{
            stroke: #666;
        }}
        
        .axis text {{
            fill: #333;
            font-size: 12px;
            font-weight: normal;
        }}
        
        .axis-label {{
            fill: #333;
            font-size: 12px;
            font-weight: normal;
        }}
        
        .grid-line {{
            stroke: #f0f0f0;
            stroke-dasharray: 2,2;
        }}
        
        .dot {{
            stroke: none;
            cursor: pointer;
            transition: all 0.2s;
            opacity: 0.6;
        }}
        
        .dot:hover {{
            stroke: #333;
            stroke-width: 1;
            r: 4;
            opacity: 0.9;
        }}
        
        .dot.highlighted {{
            stroke: #e74c3c;
            stroke-width: 2.5;
            r: 5;
            opacity: 0.9;
        }}
        
        .structure-panel {{
            background: white;
            border: 1px solid #ccc;
            padding: 25px;
            box-shadow: 0 4px 8px rgba(0,0,0,0.1);
            display: block;
            max-height: calc(100vh - 200px);
            overflow-y: auto;
            position: sticky;
            top: 20px;
        }}
        
        .structure-panel h3 {{
            margin: 0 0 20px 0;
            color: #333;
            border-bottom: 1px solid #ddd;
            padding-bottom: 10px;
            font-size: 18px;
            font-weight: normal;
        }}
        
        .structure-image {{
            text-align: center;
            margin: 20px 0;
            background: #f8f9fa;
            border: 1px solid #ddd;
            padding: 15px;
        }}
        
        .structure-image img {{
            max-width: 100%;
        }}
        
        .properties-grid {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 10px;
            margin: 15px 0;
        }}
        
        .property-box {{
            background: #f8f9fa;
            padding: 12px;
            border: 1px solid #ddd;
            border-left: 3px solid #666;
        }}
        
        .property-label {{
            font-weight: normal;
            color: #333;
            font-size: 11px;
            margin-bottom: 3px;
        }}
        
        .property-value {{
            color: #555;
            font-family: Arial, sans-serif;
            font-size: 11px;
            font-weight: normal;
        }}
        
        .smiles-box {{
            background: #f8f9fa;
            padding: 12px;
            border: 1px solid #ddd;
            margin: 15px 0;
            border-left: 3px solid #666;
        }}
        
        .smiles-box .property-label {{
            color: #333;
        }}
        
        .smiles-box .property-value {{
            word-break: break-all;
            font-size: 10px;
            font-weight: normal;
        }}
        
        .close-btn {{
            position: absolute;
            top: 20px;
            right: 25px;
            background: #666;
            color: white;
            border: none;
            width: 30px;
            height: 30px;
            cursor: pointer;
            font-size: 16px;
            font-weight: normal;
        }}
        
        .close-btn:hover {{
            background: #333;
        }}
        
        .tooltip {{
            position: absolute;
            padding: 10px 15px;
            background: rgba(0,0,0,0.9);
            color: white;
            border-radius: 8px;
            font-size: 12px;
            pointer-events: none;
            opacity: 0;
            transition: opacity 0.2s;
            z-index: 1001;
        }}
        
        .legend {{
            position: absolute;
            bottom: 10px;
            right: 10px;
            background: rgba(255,255,255,0.9);
            padding: 10px;
            border-radius: 5px;
            font-size: 11px;
            border: 1px solid #ddd;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Molecular Analysis Dashboard</h1>
        <div class="stats">
            <strong>Dataset:</strong> {len(df)} compounds | 
            <strong>Valid structures:</strong> {len(data_points)} | 
            <strong>Properties calculated:</strong> MW, LogP, TPSA, QED, SA Score
        </div>
    </div>
    
    <div class="dashboard">
        <div class="plot-container">
            <div class="plot-title">Molecular Weight vs LogP</div>
            <div id="plot1"></div>
        </div>
        <div class="plot-container">
            <div class="plot-title">TPSA vs QED</div>
            <div id="plot2"></div>
        </div>
        <div class="structure-panel" id="structurePanel">
            <button class="close-btn" onclick="hideStructure()">X</button>
            <h3 id="compoundName">Compound Name</h3>
            <div class="structure-image" id="structureImage">
                <p>Hover over a point to see the molecular structure</p>
            </div>
            <div class="properties-grid">
                <div class="property-box">
                    <div class="property-label">Molecular Weight</div>
                    <div class="property-value" id="prop-mw">--</div>
                </div>
                <div class="property-box">
                    <div class="property-label">LogP</div>
                    <div class="property-value" id="prop-logp">--</div>
                </div>
                <div class="property-box">
                    <div class="property-label">TPSA</div>
                    <div class="property-value" id="prop-tpsa">--</div>
                </div>
                <div class="property-box">
                    <div class="property-label">QED</div>
                    <div class="property-value" id="prop-qed">--</div>
                </div>
                <div class="property-box">
                    <div class="property-label">SA Score</div>
                    <div class="property-value" id="prop-sascore">--</div>
                </div>
                <div class="property-box">
                    <div class="property-label" id="prop-docking-label">Additional Property</div>
                    <div class="property-value" id="prop-docking">--</div>
                </div>
            </div>
            <div class="smiles-box">
                <div class="property-label">SMILES</div>
                <div class="property-value" id="compoundSmiles">--</div>
            </div>
        </div>
        <div class="plot-container">
            <div class="plot-title">Drug-likeness: HBA vs HBD</div>
            <div id="plot3"></div>
        </div>
        <div class="plot-container">
            <div class="plot-title">t-SNE Molecular Space</div>
            <div id="plot4"></div>
        </div>
    </div>
    
{custom_plots_html}
    
    <div class="tooltip" id="tooltip"></div>

    <script>
        console.log('=== Comprehensive Molecular Dashboard Starting ===');
        
        // Data and configurations
        const data = {data_json};
        const plotConfigs = {plot_configs_json};
        
        console.log('Loaded', data.length, 'compounds');
        console.log('Plot configurations:', plotConfigs);
        
        // Global variables
        let currentHighlighted = null;
        const plots = {{}};
        
        // Color scales
        const colorScales = {{}};
        
        // Initialize plots
        plotConfigs.forEach(config => {{
            createPlot(config);
        }});
        
        function createPlot(config) {{
            console.log('Creating plot:', config.title);
            
            // Dimensions
            const margin = {{top: 30, right: 100, bottom: 80, left: 100}};
            const width = 600 - margin.left - margin.right;
            const height = 480 - margin.top - margin.bottom;
            
            // Create SVG
            const svg = d3.select(`#${{config.id}}`)
                .append("svg")
                .attr("width", width + margin.left + margin.right)
                .attr("height", height + margin.top + margin.bottom);
                
            const g = svg.append("g")
                .attr("transform", `translate(${{margin.left}},${{margin.top}})`);
            
            // Scales
            const xScale = d3.scaleLinear()
                .domain(d3.extent(data, d => d[config.x_prop]))
                .nice()
                .range([0, width]);
                
            const yScale = d3.scaleLinear()
                .domain(d3.extent(data, d => d[config.y_prop]))
                .nice()
                .range([height, 0]);
            
            // Color scale (only if color property exists)
            let colorScale = null;
            let colorExtent = null;
            if (config.color_prop && data.some(d => d.hasOwnProperty(config.color_prop) && d[config.color_prop] !== null)) {{
                colorExtent = d3.extent(data, d => d[config.color_prop]);
                colorScale = d3.scaleSequential(d3.interpolateViridis)
                    .domain(colorExtent);
                colorScales[config.id] = colorScale;
            }}
            
            // Add grid
            g.selectAll(".grid-line.vertical")
                .data(xScale.ticks())
                .enter().append("line")
                .attr("class", "grid-line vertical")
                .attr("x1", d => xScale(d))
                .attr("x2", d => xScale(d))
                .attr("y1", 0)
                .attr("y2", height);
                
            g.selectAll(".grid-line.horizontal")
                .data(yScale.ticks())
                .enter().append("line")
                .attr("class", "grid-line horizontal")
                .attr("x1", 0)
                .attr("x2", width)
                .attr("y1", d => yScale(d))
                .attr("y2", d => yScale(d));
            
            // Add axes
            g.append("g")
                .attr("class", "axis")
                .attr("transform", `translate(0,${{height}})`)
                .call(d3.axisBottom(xScale));
                
            g.append("g")
                .attr("class", "axis")
                .call(d3.axisLeft(yScale));
            
            // Add axis labels
            g.append("text")
                .attr("class", "axis-label")
                .attr("transform", `translate(${{width/2}}, ${{height + 45}})`)
                .style("text-anchor", "middle")
                .text(config.x_label);
                
            g.append("text")
                .attr("class", "axis-label")
                .attr("transform", "rotate(-90)")
                .attr("y", 0 - margin.left + 25)
                .attr("x", 0 - (height / 2))
                .style("text-anchor", "middle")
                .text(config.y_label);
            
            // Add color legend (only if color scale exists)
            if (colorScale && colorExtent) {{
                const legend = g.append("g")
                    .attr("class", "legend")
                    .attr("transform", `translate(${{width + 10}}, 20)`);
                    
                const legendHeight = 100;
                const legendWidth = 20;
                
                // Create gradient for legend
                const gradient = svg.append("defs")
                    .append("linearGradient")
                    .attr("id", `gradient-${{config.id}}`)
                    .attr("x1", "0%")
                    .attr("y1", "100%")
                    .attr("x2", "0%")
                    .attr("y2", "0%");
                    
                const numStops = 10;
                for (let i = 0; i <= numStops; i++) {{
                    const offset = i / numStops;
                    const value = colorExtent[0] + offset * (colorExtent[1] - colorExtent[0]);
                    gradient.append("stop")
                        .attr("offset", `${{offset * 100}}%`)
                        .attr("stop-color", colorScale(value));
                }}
                
                legend.append("rect")
                    .attr("width", legendWidth)
                    .attr("height", legendHeight)
                    .style("fill", `url(#gradient-${{config.id}})`);
                    
                legend.append("text")
                    .attr("x", legendWidth + 5)
                    .attr("y", 0)
                    .style("font-size", "11px")
                    .text(colorExtent[1].toFixed(2));
                    
                legend.append("text")
                    .attr("x", legendWidth + 5)
                    .attr("y", legendHeight)
                    .style("font-size", "11px")
                    .text(colorExtent[0].toFixed(2));
                    
                legend.append("text")
                    .attr("x", legendWidth + 5)
                    .attr("y", legendHeight + 15)
                    .style("font-size", "10px")
                    .text(config.color_label);
            }}
            
            // Add dots
            const dots = g.selectAll(".dot")
                .data(data)
                .enter().append("circle")
                .attr("class", "dot")
                .attr("cx", d => xScale(d[config.x_prop]))
                .attr("cy", d => yScale(d[config.y_prop]))
                .attr("r", 2.5)
                .attr("fill", d => colorScale ? colorScale(d[config.color_prop]) : "#69b3a2")
                .on("mouseover", function(event, d) {{
                    highlightCompound(d.id);
                    showTooltip(event, d, config);
                    showStructure(d);
                }})
                .on("mouseout", function() {{
                    hideTooltip();
                }})
                .on("click", function(event, d) {{
                    showStructure(d);
                }});
            
            // Store plot reference
            plots[config.id] = {{
                svg: svg,
                g: g,
                dots: dots,
                config: config
            }};
            
            console.log(`Plot ${{config.id}} created with ${{dots.size()}} points`);
        }}
        
        function highlightCompound(compoundId) {{
            // Remove previous highlights
            if (currentHighlighted !== null) {{
                Object.values(plots).forEach(plot => {{
                    plot.dots.classed("highlighted", false);
                }});
            }}
            
            // Add new highlights
            Object.values(plots).forEach(plot => {{
                plot.dots.classed("highlighted", d => d.id === compoundId);
            }});
            
            currentHighlighted = compoundId;
        }}
        
        function showTooltip(event, compound, config) {{
            const tooltip = d3.select("#tooltip");
            tooltip.style("opacity", 1)
                .html(`
                    <strong>${{compound.title}}</strong><br>
                    ${{config.x_label}}: ${{compound[config.x_prop].toFixed(2)}}<br>
                    ${{config.y_label}}: ${{compound[config.y_prop].toFixed(2)}}${{config.color_prop && compound[config.color_prop] !== undefined ? '<br>' + config.color_label + ': ' + compound[config.color_prop].toFixed(2) : ''}}
                `)
                .style("left", (event.pageX + 10) + "px")
                .style("top", (event.pageY - 10) + "px");
        }}
        
        function hideTooltip() {{
            d3.select("#tooltip").style("opacity", 0);
        }}
        
        function showStructure(compound) {{
            console.log('Showing structure for:', compound.title);
            
            document.getElementById('compoundName').textContent = compound.title;
            
            if (compound.image) {{
                document.getElementById('structureImage').innerHTML = 
                    `<img src="${{compound.image}}" alt="Molecular Structure" />`;
            }} else {{
                document.getElementById('structureImage').innerHTML = 
                    '<p style="color: #999;">No structure available</p>';
            }}
            
            // Update properties
            document.getElementById('prop-mw').textContent = compound.mw.toFixed(1) + ' Da';
            document.getElementById('prop-logp').textContent = compound.logp.toFixed(2);
            document.getElementById('prop-tpsa').textContent = compound.tpsa.toFixed(1) + ' Ų';
            document.getElementById('prop-qed').textContent = compound.qed.toFixed(3);
            document.getElementById('prop-sascore').textContent = compound.sascore.toFixed(2);
            
            // Update custom property if available (use first custom property)
            let customPropFound = false;
            for (let i = 0; i < 10; i++) {{ // Check up to 10 custom properties
                if (compound[`custom_property_${{i}}`] !== undefined) {{
                    document.getElementById('prop-docking-label').textContent = '{properties_from_input[0] if properties_from_input else "Custom Property"}';
                    document.getElementById('prop-docking').textContent = compound[`custom_property_${{i}}`].toFixed(2);
                    customPropFound = true;
                    break;
                }}
            }}
            if (!customPropFound) {{
                document.getElementById('prop-docking-label').textContent = 'Additional Property';
                document.getElementById('prop-docking').textContent = '--';
            }}
            document.getElementById('compoundSmiles').textContent = compound.smiles;
            
            // Panel is always visible, no need to show/hide
        }}
        
        function hideStructure() {{
            // Clear highlights
            if (currentHighlighted !== null) {{
                Object.values(plots).forEach(plot => {{
                    plot.dots.classed("highlighted", false);
                }});
                currentHighlighted = null;
            }}
            // Reset to default content
            document.getElementById('compoundName').textContent = 'Compound Name';
            document.getElementById('structureImage').innerHTML = '<p>Hover over a point to see the molecular structure</p>';
            document.getElementById('prop-mw').textContent = '--';
            document.getElementById('prop-logp').textContent = '--';
            document.getElementById('prop-tpsa').textContent = '--';
            document.getElementById('prop-qed').textContent = '--';
            document.getElementById('prop-sascore').textContent = '--';
            document.getElementById('prop-docking').textContent = '--';
            document.getElementById('compoundSmiles').textContent = '--';
        }}
        
        // Make hideStructure globally available
        window.hideStructure = hideStructure;
        
        console.log('=== Dashboard setup complete ===');
        console.log('Hover over dots to see synchronized highlighting!');
        
    </script>
</body>
</html>
"""
    
    # Save HTML file
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"\n🎉 Comprehensive dashboard saved to: {output_file}")
    print(f"📊 Features:")
    print(f"   • 4 interactive scatter plots with synchronized highlighting")
    print(f"   • t-SNE molecular space visualization") 
    print(f"   • Color-coded third dimension on each plot")
    print(f"   • Comprehensive molecular descriptors")
    print(f"   • Structure panel with detailed properties")
    print(f"🚀 Open in your browser and hover over points!")


def main():
    parser = argparse.ArgumentParser(description="Comprehensive Molecular Analysis Dashboard")
    parser.add_argument('input', help='Input CSV file with SMILES, Name columns')
    parser.add_argument('--output', '-o', default='molecular_dashboard.html', 
                       help='Output HTML file (default: molecular_dashboard.html)')
    parser.add_argument('--property_from_input', nargs='+', help='Column names from CSV to use for plots (multiple values allowed)')
    
    args = parser.parse_args()
    
    try:
        # Load and process data
        df = load_and_process_data(args.input, args.property_from_input)
        
        # Calculate comprehensive descriptors
        df = calculate_comprehensive_descriptors(df)
        
        # Calculate t-SNE coordinates
        df = calculate_tsne_coordinates(df)
        
        # Create dashboard
        print("Creating comprehensive dashboard...")
        create_comprehensive_dashboard(df, args.output, args.property_from_input);
        
        print(f"\n✅ Success! Open {args.output} in your web browser.");
        print(f"🎯 Interactive features:");
        print(f"   • Hover over any point to highlight it across all plots");
        print(f"   • View molecular structures and detailed properties");
        print(f"   • Explore chemical space with t-SNE visualization");
        
    except Exception as e:
        print(f"❌ Error: {e}");
        import traceback;
        traceback.print_exc();
        sys.exit(1);


if __name__ == "__main__":
    main();