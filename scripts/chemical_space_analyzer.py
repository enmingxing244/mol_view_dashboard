#!/usr/bin/env python3
"""
Chemical Space Analysis Dashboard
Analyzes multiple CSV files with molecular data and creates PCA/t-SNE visualization
- Reads SMILES and Name columns from multiple CSV files
- Associates each file with a user-defined label
- Performs PCA and t-SNE on molecular fingerprints
- Creates interactive dashboard with compound source labeling
"""

import pandas as pd
import numpy as np
import base64
from io import BytesIO
import json
import argparse
import sys
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# RDKit imports
try:
    from rdkit import Chem
    from rdkit.Chem import Draw, rdFingerprintGenerator
    from rdkit import DataStructs
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


def load_multiple_csv_files(csv_files, labels):
    """Load multiple CSV files and combine them with source labels"""
    if len(csv_files) != len(labels):
        raise ValueError("Number of CSV files must match number of labels")
    
    combined_data = []
    
    for csv_file, label in zip(csv_files, labels):
        print(f"Loading data from: {csv_file} (label: {label})")
        
        df = pd.read_csv(csv_file)
        print(f"  Shape: {df.shape}")
        print(f"  Columns: {list(df.columns)}")
        
        # Check required columns
        required_cols = ['SMILES']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns in {csv_file}: {missing_cols}")
        
        # Handle optional Name column
        if 'Name' not in df.columns:
            df['Name'] = [f"Compound_{i+1}" for i in range(len(df))]
        
        # Add source label
        df['Source'] = label
        df['File'] = csv_file
        
        # Create molecule objects
        df['Mol'] = df['SMILES'].apply(lambda x: Chem.MolFromSmiles(x) if pd.notna(x) else None)
        
        # Filter valid molecules
        valid_df = df[df['Mol'].notna()].reset_index(drop=True)
        print(f"  Valid molecules: {len(valid_df)}")
        
        combined_data.append(valid_df)
    
    # Combine all dataframes
    combined_df = pd.concat(combined_data, ignore_index=True)
    print(f"\nTotal combined dataset: {len(combined_df)} compounds from {len(csv_files)} sources")
    
    return combined_df


def calculate_molecular_fingerprints(df):
    """Calculate molecular fingerprints for all compounds"""
    print(f"Calculating molecular fingerprints for {len(df)} compounds...")
    
    # Generate molecular fingerprints
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fingerprints = []
    
    for i, mol in enumerate(df['Mol']):
        if i % 100 == 0:
            print(f"  Processing molecule {i+1}/{len(df)}")
        
        fp = generator.GetFingerprint(mol)
        fp_array = np.zeros((2048,))
        DataStructs.ConvertToNumpyArray(fp, fp_array)
        fingerprints.append(fp_array)
    
    fp_matrix = np.array(fingerprints)
    print(f"Fingerprint matrix shape: {fp_matrix.shape}")
    
    return fp_matrix


def perform_dimensionality_reduction(fp_matrix):
    """Perform PCA and t-SNE on fingerprint matrix"""
    print("Performing dimensionality reduction...")
    
    # Standardize features
    scaler = StandardScaler()
    fp_scaled = scaler.fit_transform(fp_matrix)
    
    # PCA
    print("  Running PCA...")
    pca = PCA(n_components=2, random_state=42)
    pca_coords = pca.fit_transform(fp_scaled)
    
    print(f"  PCA explained variance: {pca.explained_variance_ratio_[0]:.3f}, {pca.explained_variance_ratio_[1]:.3f}")
    print(f"  Total explained variance: {sum(pca.explained_variance_ratio_):.3f}")
    
    # t-SNE
    print("  Running t-SNE...")
    perplexity = min(30, len(fp_scaled)//4)
    tsne = TSNE(n_components=2, random_state=42, perplexity=perplexity)
    tsne_coords = tsne.fit_transform(fp_scaled)
    
    return pca_coords, tsne_coords, pca.explained_variance_ratio_


def create_chemical_space_dashboard(df, pca_coords, tsne_coords, pca_variance, output_file="chemical_space_dashboard.html"):
    """Create chemical space analysis dashboard"""
    
    print(f"Creating chemical space dashboard with {len(df)} compounds")
    
    # Add coordinates to dataframe
    df['PCA_X'] = pca_coords[:, 0]
    df['PCA_Y'] = pca_coords[:, 1]
    df['tSNE_X'] = tsne_coords[:, 0]
    df['tSNE_Y'] = tsne_coords[:, 1]
    
    # Generate structure images and prepare data
    print("Generating molecular structures...")
    data_points = []
    
    for idx, row in df.iterrows():
        if idx % 50 == 0:
            print(f"  Generating structure {idx+1}/{len(df)}")
        
        img_base64 = mol_to_base64_png(row['Mol'])
        
        data_point = {
            'id': int(idx),
            'name': str(row['Name']),
            'smiles': str(row['SMILES']),
            'source': str(row['Source']),
            'file': str(row['File']),
            'pca_x': float(row['PCA_X']),
            'pca_y': float(row['PCA_Y']),
            'tsne_x': float(row['tSNE_X']),
            'tsne_y': float(row['tSNE_Y']),
            'image': img_base64
        }
        data_points.append(data_point)
    
    # Convert to JSON
    data_json = json.dumps(data_points, indent=2)
    
    # Get unique sources for color mapping
    sources = df['Source'].unique()
    source_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
                     '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    
    # Create HTML dashboard
    html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Chemical Space Analysis Dashboard</title>
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
            max-width: 1000px;
            color: #666;
            text-align: center;
        }}
        
        .dashboard {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 20px;
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
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
            position: fixed;
            top: 20px;
            right: 20px;
            width: 400px;
            background: white;
            border: 1px solid #ccc;
            padding: 25px;
            box-shadow: 0 4px 8px rgba(0,0,0,0.1);
            display: none;
            z-index: 1000;
            max-height: 85vh;
            overflow-y: auto;
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
            top: 10px;
            right: 10px;
            background: rgba(255,255,255,0.9);
            padding: 10px;
            border: 1px solid #ddd;
            font-size: 11px;
        }}
        
        .legend-item {{
            display: flex;
            align-items: center;
            margin: 3px 0;
        }}
        
        .legend-color {{
            width: 12px;
            height: 12px;
            margin-right: 5px;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Chemical Space Analysis Dashboard</h1>
        <div class="stats">
            <strong>Dataset:</strong> {len(df)} compounds from {len(sources)} sources | 
            <strong>Sources:</strong> {', '.join(sources)}
        </div>
    </div>
    
    <div class="dashboard">
        <div class="plot-container">
            <div class="plot-title">PCA Analysis (PC1: {pca_variance[0]:.1%}, PC2: {pca_variance[1]:.1%})</div>
            <div id="pca-plot"></div>
        </div>
        <div class="plot-container">
            <div class="plot-title">t-SNE Analysis</div>
            <div id="tsne-plot"></div>
        </div>
    </div>
    
    <div class="structure-panel" id="structurePanel">
        <button class="close-btn" onclick="hideStructure()">X</button>
        <h3 id="compoundName">Compound Name</h3>
        <div class="structure-image" id="structureImage">
            <p>Hover over a point to see the molecular structure</p>
        </div>
        <div class="properties-grid">
            <div class="property-box">
                <div class="property-label">Compound Name</div>
                <div class="property-value" id="prop-name">--</div>
            </div>
            <div class="property-box">
                <div class="property-label">Source</div>
                <div class="property-value" id="prop-source">--</div>
            </div>
            <div class="property-box">
                <div class="property-label">File</div>
                <div class="property-value" id="prop-file">--</div>
            </div>
            <div class="property-box">
                <div class="property-label">ID</div>
                <div class="property-value" id="prop-id">--</div>
            </div>
        </div>
        <div class="smiles-box">
            <div class="property-label">SMILES</div>
            <div class="property-value" id="compoundSmiles">--</div>
        </div>
    </div>
    
    <div class="tooltip" id="tooltip"></div>

    <script>
        console.log('=== Chemical Space Analysis Dashboard Starting ===');
        
        // Data
        const data = {data_json};
        const sources = {json.dumps(list(sources))};
        const sourceColors = {json.dumps(source_colors[:len(sources)])};
        
        console.log('Loaded', data.length, 'compounds from', sources.length, 'sources');
        
        // Create color mapping
        const colorMap = {{}};
        sources.forEach((source, i) => {{
            colorMap[source] = sourceColors[i % sourceColors.length];
        }});
        
        // Global variables
        let currentHighlighted = null;
        const plots = {{}};
        
        // Create plots
        createPlot({{
            id: 'pca-plot',
            title: 'PCA Analysis',
            x_prop: 'pca_x',
            y_prop: 'pca_y',
            x_label: 'PC1 ({pca_variance[0]:.1%})',
            y_label: 'PC2 ({pca_variance[1]:.1%})'
        }});
        
        createPlot({{
            id: 'tsne-plot',
            title: 't-SNE Analysis',
            x_prop: 'tsne_x',
            y_prop: 'tsne_y',
            x_label: 't-SNE Dimension 1',
            y_label: 't-SNE Dimension 2'
        }});
        
        function createPlot(config) {{
            console.log('Creating plot:', config.title);
            
            // Dimensions
            const margin = {{top: 30, right: 120, bottom: 80, left: 100}};
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
            
            // Add legend
            const legend = g.append("g")
                .attr("class", "legend")
                .attr("transform", `translate(${{width + 10}}, 20)`);
                
            sources.forEach((source, i) => {{
                const legendItem = legend.append("g")
                    .attr("class", "legend-item")
                    .attr("transform", `translate(0, ${{i * 20}})`);
                    
                legendItem.append("circle")
                    .attr("class", "legend-color")
                    .attr("r", 6)
                    .attr("fill", colorMap[source]);
                    
                legendItem.append("text")
                    .attr("x", 15)
                    .attr("y", 0)
                    .attr("dy", "0.35em")
                    .style("font-size", "10px")
                    .text(source);
            }});
            
            // Add dots
            const dots = g.selectAll(".dot")
                .data(data)
                .enter().append("circle")
                .attr("class", "dot")
                .attr("cx", d => xScale(d[config.x_prop]))
                .attr("cy", d => yScale(d[config.y_prop]))
                .attr("r", 2.5)
                .attr("fill", d => colorMap[d.source])
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
                    <strong>${{compound.name}}</strong><br>
                    Source: ${{compound.source}}<br>
                    ${{config.x_label}}: ${{compound[config.x_prop].toFixed(2)}}<br>
                    ${{config.y_label}}: ${{compound[config.y_prop].toFixed(2)}}
                `)
                .style("left", (event.pageX + 10) + "px")
                .style("top", (event.pageY - 10) + "px");
        }}
        
        function hideTooltip() {{
            d3.select("#tooltip").style("opacity", 0);
        }}
        
        function showStructure(compound) {{
            console.log('Showing structure for:', compound.name);
            
            document.getElementById('compoundName').textContent = compound.name;
            
            if (compound.image) {{
                document.getElementById('structureImage').innerHTML = 
                    `<img src="${{compound.image}}" alt="Molecular Structure" />`;
            }} else {{
                document.getElementById('structureImage').innerHTML = 
                    '<p style="color: #999;">No structure available</p>';
            }}
            
            // Update properties
            document.getElementById('prop-name').textContent = compound.name;
            document.getElementById('prop-source').textContent = compound.source;
            document.getElementById('prop-file').textContent = compound.file.split('/').pop();
            document.getElementById('prop-id').textContent = compound.id;
            document.getElementById('compoundSmiles').textContent = compound.smiles;
            
            document.getElementById('structurePanel').style.display = 'block';
        }}
        
        function hideStructure() {{
            document.getElementById('structurePanel').style.display = 'none';
            // Clear highlights
            if (currentHighlighted !== null) {{
                Object.values(plots).forEach(plot => {{
                    plot.dots.classed("highlighted", false);
                }});
                currentHighlighted = null;
            }}
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
    
    print(f"\n🎉 Chemical space dashboard saved to: {output_file}")
    print(f"📊 Features:")
    print(f"   • PCA analysis with explained variance: {pca_variance[0]:.1%} + {pca_variance[1]:.1%} = {sum(pca_variance):.1%}")
    print(f"   • t-SNE molecular space visualization")
    print(f"   • {len(sources)} data sources with color coding")
    print(f"   • Interactive compound highlighting across plots")
    print(f"   • Source labeling and compound metadata")
    print(f"🚀 Open in your browser and explore the chemical space!")


def main():
    parser = argparse.ArgumentParser(description="Chemical Space Analysis Dashboard")
    parser.add_argument('--csv_files', nargs='+', required=True,
                       help='Input CSV files with SMILES and Name columns')
    parser.add_argument('--labels', nargs='+', required=True,
                       help='Labels for each CSV file (same order as csv_files)')
    parser.add_argument('--output', '-o', default='chemical_space_dashboard.html',
                       help='Output HTML file (default: chemical_space_dashboard.html)')
    
    args = parser.parse_args()
    
    try:
        # Load and combine data from multiple CSV files
        df = load_multiple_csv_files(args.csv_files, args.labels)
        
        # Calculate molecular fingerprints
        fp_matrix = calculate_molecular_fingerprints(df)
        
        # Perform dimensionality reduction
        pca_coords, tsne_coords, pca_variance = perform_dimensionality_reduction(fp_matrix)
        
        # Create dashboard
        print("Creating chemical space dashboard...")
        create_chemical_space_dashboard(df, pca_coords, tsne_coords, pca_variance, args.output)
        
        print(f"\n✅ Success! Open {args.output} in your web browser.")
        print(f"🎯 Interactive features:")
        print(f"   • Hover over any point to highlight it across both plots")
        print(f"   • View molecular structures and source information")
        print(f"   • Color-coded data sources with legend")
        print(f"   • PCA and t-SNE chemical space exploration")
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()