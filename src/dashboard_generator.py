#!/usr/bin/env python3
"""
Enhanced Dashboard Generator
Creates interactive HTML dashboard with:
- User-configurable property plots
- Chemical space visualization (PCA + t-SNE)
- Synchronized highlighting across all plots
- Docking pose visualization
- Professional styling with color bars
"""

import json
import logging
from typing import Dict, List, Any, Optional
from pathlib import Path
from structure_viewer import StructureViewer


class DashboardGenerator:
    """Generates enhanced interactive HTML dashboard"""
    
    def __init__(self, config_manager):
        """
        Initialize dashboard generator
        
        Args:
            config_manager: Configuration manager instance
        """
        self.config = config_manager
        self.logger = logging.getLogger(__name__)
        self.structure_viewer = StructureViewer(config_manager)
    
    def generate_dashboard(self, 
                         data_points: List[Dict[str, Any]], 
                         data_summary: Dict[str, Any],
                         docking_data: Optional[List[Dict[str, Any]]] = None) -> str:
        """
        Generate complete HTML dashboard
        
        Args:
            data_points: List of compound data dictionaries
            data_summary: Dataset summary information
            docking_data: Optional docking visualization data
            
        Returns:
            HTML content string
        """
        self.logger.info(f"Generating dashboard for {len(data_points)} compounds")
        
        # Prepare data for JavaScript
        data_json = json.dumps(data_points, indent=2)
        summary_json = json.dumps(data_summary, indent=2)
        docking_json = json.dumps(docking_data or [], indent=2)
        
        # Check if docking is enabled
        docking_enabled = self.config.is_docking_enabled() and docking_data
        
        # Prepare 3D structure data if docking is available
        structure_viewer_html = ""
        if docking_enabled:
            structure_data = self.structure_viewer.prepare_docking_results_for_viewing(docking_data)
            structure_viewer_html = self.structure_viewer.generate_3dmol_html(structure_data)
        
        # Get available properties for plot configuration
        available_props = data_summary.get('available_properties', [])
        properties_json = json.dumps(available_props)
        
        # Get visualization configuration
        viz_config = self.config.get('visualization', {})
        title = viz_config.get('title', 'Molecular Visualization and Analysis')
        color_scheme = viz_config.get('style.color_scheme', 'viridis')
        
        # Generate HTML content
        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <script src="https://d3js.org/d3.v7.min.js"></script>
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    {"<script src='https://unpkg.com/3dmol@latest/build/3Dmol-min.js'></script>" if docking_enabled else ""}
    {self._generate_css()}
</head>
<body>
    {self._generate_header(title, data_summary)}
    {self._generate_tabs(docking_enabled)}
    {self._generate_main_tab(data_summary)}
    {self._generate_docking_tab(structure_viewer_html) if docking_enabled else ""}
    {self._generate_structure_panel()}
    <div class="tooltip" id="tooltip"></div>
    
    <script>
        // Global data and configuration
        const data = {data_json};
        const summary = {summary_json};
        const dockingData = {docking_json};
        const availableProperties = {properties_json};
        const colorScheme = '{color_scheme}';
        const dockingEnabled = {str(docking_enabled).lower()};
        
        {self._generate_javascript()}
    </script>
</body>
</html>"""
        
        return html_content
    
    def _generate_css(self) -> str:
        """Generate CSS styles for the dashboard"""
        return """
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', 'Roboto', sans-serif;
            background: #ffffff;
            color: #333;
            line-height: 1.6;
        }
        
        .header {
            background: #ffffff;
            border-bottom: 2px solid #e5e5e5;
            padding: 1.5rem 2rem;
            text-align: center;
        }
        
        .header h1 {
            font-size: 2.2rem;
            font-weight: 300;
            color: #2c3e50;
            margin-bottom: 0.5rem;
        }
        
        .stats {
            background: #f8f9fa;
            padding: 1rem 2rem;
            border: 1px solid #dee2e6;
            border-radius: 4px;
            margin: 1rem auto;
            max-width: 1000px;
            font-size: 0.95rem;
            color: #495057;
        }
        
        .tabs {
            display: flex;
            background: #f8f9fa;
            border-bottom: 1px solid #dee2e6;
            padding: 0 2rem;
        }
        
        .tab-button {
            background: none;
            border: none;
            padding: 1rem 2rem;
            cursor: pointer;
            font-size: 1rem;
            color: #6c757d;
            border-bottom: 3px solid transparent;
            transition: all 0.3s ease;
        }
        
        .tab-button.active {
            color: #2c3e50;
            border-bottom-color: #3498db;
            background: #ffffff;
        }
        
        .tab-button:hover {
            color: #2c3e50;
            background: #e9ecef;
        }
        
        .tab-content {
            display: none;
            padding: 2rem;
        }
        
        .tab-content.active {
            display: block;
        }
        
        .controls {
            background: #f8f9fa;
            border: 1px solid #dee2e6;
            border-radius: 4px;
            padding: 1.5rem;
            margin-bottom: 2rem;
        }
        
        .controls h3 {
            margin-bottom: 1rem;
            color: #2c3e50;
            font-weight: 400;
        }
        
        .control-group {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 1rem;
            margin-bottom: 1rem;
        }
        
        .control-item {
            display: flex;
            flex-direction: column;
        }
        
        .control-item label {
            font-weight: 500;
            margin-bottom: 0.5rem;
            color: #495057;
        }
        
        .control-item select {
            padding: 0.5rem;
            border: 1px solid #ced4da;
            border-radius: 4px;
            background: #ffffff;
            font-size: 0.9rem;
        }
        
        .dashboard {
            display: grid;
            grid-template-columns: 1fr 1fr;
            grid-template-rows: 1fr 1fr;
            gap: 1.5rem;
            max-width: 1400px;
            margin: 0 auto;
            height: 800px;
        }
        
        .plot-container {
            background: #ffffff;
            border: 1px solid #dee2e6;
            border-radius: 4px;
            padding: 1.5rem;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        }
        
        .plot-title {
            text-align: center;
            margin-bottom: 1rem;
            font-size: 1.1rem;
            font-weight: 500;
            color: #2c3e50;
        }
        
        .axis {
            stroke: #6c757d;
            stroke-width: 1;
        }
        
        .axis text {
            fill: #495057;
            font-size: 11px;
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
        }
        
        .axis-label {
            fill: #2c3e50;
            font-size: 12px;
            font-weight: 500;
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
        }
        
        .grid-line {
            stroke: #f1f3f4;
            stroke-width: 1;
        }
        
        .dot {
            stroke: none;
            cursor: pointer;
            transition: all 0.2s ease;
            opacity: 0.7;
        }
        
        .dot:hover {
            stroke: #2c3e50;
            stroke-width: 2;
            r: 4.5;
            opacity: 1;
        }
        
        .dot.highlighted {
            stroke: #e74c3c;
            stroke-width: 2.5;
            r: 5;
            opacity: 1;
        }
        
        .colorbar {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
        }
        
        .colorbar text {
            fill: #495057;
            font-size: 10px;
        }
        
        .colorbar-title {
            fill: #2c3e50;
            font-size: 11px;
            font-weight: 500;
        }
        
        .structure-panel {
            background: #ffffff;
            border: 1px solid #dee2e6;
            border-radius: 4px;
            padding: 1rem;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
            overflow: hidden;
            display: flex;
            flex-direction: column;
            height: 100%;
        }
        
        .structure-panel h3 {
            margin-bottom: 0.5rem;
            color: #2c3e50;
            font-weight: 500;
            font-size: 1rem;
            border-bottom: 1px solid #dee2e6;
            padding-bottom: 0.3rem;
            flex-shrink: 0;
        }
        
        .structure-image {
            text-align: center;
            margin: 0.5rem 0;
            background: #f8f9fa;
            border: 1px solid #dee2e6;
            border-radius: 4px;
            padding: 0.5rem;
            height: 140px;
            display: flex;
            align-items: center;
            justify-content: center;
            flex-shrink: 0;
        }
        
        .structure-image img {
            max-width: 100%;
            max-height: 120px;
            border-radius: 4px;
        }
        
        .structure-image .placeholder {
            color: #6c757d;
            font-style: italic;
            font-size: 0.85rem;
        }
        
        .properties-grid {
            display: grid;
            grid-template-columns: 1fr 1fr 1fr;
            gap: 0.4rem;
            margin: 0.5rem 0;
            flex: 1;
            overflow-y: auto;
        }
        
        .property-box {
            background: #f8f9fa;
            padding: 0.4rem;
            border: 1px solid #dee2e6;
            border-radius: 3px;
            border-left: 2px solid #3498db;
        }
        
        .property-label {
            font-weight: 500;
            color: #2c3e50;
            font-size: 0.7rem;
            margin-bottom: 0.1rem;
            line-height: 1.2;
        }
        
        .property-value {
            color: #495057;
            font-size: 0.75rem;
            font-family: 'SF Mono', Monaco, 'Cascadia Code', monospace;
            line-height: 1.2;
        }
        
        .smiles-box {
            background: #f8f9fa;
            padding: 0.5rem;
            border: 1px solid #dee2e6;
            border-radius: 4px;
            margin: 0.5rem 0 0 0;
            border-left: 3px solid #3498db;
            flex-shrink: 0;
        }
        
        .smiles-box .property-value {
            word-break: break-all;
            font-size: 0.65rem;
            line-height: 1.3;
        }
        
        .close-btn {
            position: absolute;
            top: 1rem;
            right: 1rem;
            background: #6c757d;
            color: white;
            border: none;
            width: 32px;
            height: 32px;
            border-radius: 4px;
            cursor: pointer;
            font-size: 16px;
            transition: background 0.2s ease;
        }
        
        .close-btn:hover {
            background: #495057;
        }
        
        .tooltip {
            position: absolute;
            padding: 0.75rem 1rem;
            background: rgba(33, 37, 41, 0.95);
            color: white;
            border-radius: 4px;
            font-size: 0.85rem;
            pointer-events: none;
            opacity: 0;
            transition: opacity 0.2s ease;
            z-index: 1000;
            max-width: 250px;
        }
        
        .docking-viewer {
            width: 100%;
            height: 500px;
            border: 1px solid #dee2e6;
            border-radius: 4px;
            background: #f8f9fa;
        }
        
        .docking-controls {
            margin-bottom: 1rem;
        }
        
        .compound-list {
            max-height: 400px;
            overflow-y: auto;
            border: 1px solid #dee2e6;
            border-radius: 4px;
        }
        
        .compound-item {
            padding: 0.75rem;
            border-bottom: 1px solid #dee2e6;
            cursor: pointer;
            transition: background 0.2s ease;
        }
        
        .compound-item:hover {
            background: #f8f9fa;
        }
        
        .compound-item.selected {
            background: #e3f2fd;
            border-left: 3px solid #2196f3;
        }
        
        .compound-name {
            font-weight: 500;
            margin-bottom: 0.25rem;
        }
        
        .binding-energy {
            font-size: 0.85rem;
            color: #6c757d;
        }
        
        .energy-good {
            color: #28a745;
            font-weight: 500;
        }
        
        .energy-moderate {
            color: #ffc107;
            font-weight: 500;
        }
        
        .energy-poor {
            color: #dc3545;
            font-weight: 500;
        }
        
        @media (max-width: 1200px) {
            .dashboard {
                max-width: 1000px;
                height: 700px;
            }
        }
        
        @media (max-width: 900px) {
            .dashboard {
                grid-template-columns: 1fr;
                grid-template-rows: auto auto auto auto;
                height: auto;
                gap: 1rem;
            }
            
            .structure-panel {
                order: 2;
            }
        }
        
        @media (max-width: 768px) {
            .control-group {
                grid-template-columns: 1fr;
            }
            
            .dashboard {
                padding: 1rem;
            }
        }
    </style>"""
    
    def _generate_header(self, title: str, data_summary: Dict[str, Any]) -> str:
        """Generate header section"""
        total_compounds = data_summary.get('total_compounds', 0)
        available_props = data_summary.get('available_properties', [])
        has_pca = data_summary.get('has_pca', False)
        has_tsne = data_summary.get('has_tsne', False)
        
        analysis_types = []
        if has_pca:
            analysis_types.append("PCA")
        if has_tsne:
            analysis_types.append("t-SNE")
        
        docking_summary = ""
        if self.config.is_docking_enabled():
            docking_summary = " | <strong>Docking:</strong> Enabled"
        
        return f"""
    <div class="header">
        <h1>{title}</h1>
        <div class="stats">
            <strong>Dataset:</strong> {total_compounds} compounds | 
            <strong>Properties:</strong> {len(available_props)} available | 
            <strong>Chemical Space:</strong> {', '.join(analysis_types) if analysis_types else 'None'}{docking_summary}
        </div>
    </div>"""
    
    def _generate_tabs(self, docking_enabled: bool) -> str:
        """Generate tab navigation"""
        docking_tab = """
        <button class="tab-button" onclick="showTab('docking')">Docking Analysis</button>""" if docking_enabled else ""
        
        return f"""
    <div class="tabs">
        <button class="tab-button active" onclick="showTab('main')">Property & Chemical Space Analysis</button>{docking_tab}
    </div>"""
    
    def _generate_main_tab(self, data_summary: Dict[str, Any]) -> str:
        """Generate main analysis tab"""
        available_props = data_summary.get('available_properties', [])
        
        # Default plot configuration
        default_config = self.config.get('visualization.property_plots.primary_plot', {})
        default_x = default_config.get('x_axis', available_props[0] if available_props else 'MW')
        default_y = default_config.get('y_axis', available_props[1] if len(available_props) > 1 else 'LogP')
        default_color = default_config.get('color_by', available_props[2] if len(available_props) > 2 else 'QED')
        
        return f"""
    <div id="main" class="tab-content active">
        <div class="controls">
            <h3>Plot Configuration</h3>
            <div class="control-group">
                <div class="control-item">
                    <label for="x-axis-select">X-Axis Property:</label>
                    <select id="x-axis-select" onchange="updatePropertyPlot()">
                        {self._generate_property_options(available_props, default_x)}
                    </select>
                </div>
                <div class="control-item">
                    <label for="y-axis-select">Y-Axis Property:</label>
                    <select id="y-axis-select" onchange="updatePropertyPlot()">
                        {self._generate_property_options(available_props, default_y)}
                    </select>
                </div>
                <div class="control-item">
                    <label for="color-select">Color Property:</label>
                    <select id="color-select" onchange="updatePropertyPlot()">
                        <option value="">None</option>
                        {self._generate_property_options(available_props, default_color)}
                    </select>
                </div>
            </div>
        </div>
        
        <div class="dashboard">
            <!-- Row 1, Col 1: Property plot -->
            <div class="plot-container">
                <div class="plot-title" id="property-plot-title">Property Visualization</div>
                <div id="property-plot"></div>
            </div>
            
            <!-- Row 1, Col 2: Structure panel -->
            <div class="structure-panel" id="structurePanel">
                <h3 id="compoundName">Molecular Details</h3>
                <div class="structure-image" id="structureImage">
                    <div class="placeholder">Hover over a data point to view molecular structure</div>
                </div>
                <div class="properties-grid" id="propertiesGrid">
                    <!-- Properties will be populated dynamically -->
                </div>
                <div class="smiles-box">
                    <div class="property-label">SMILES</div>
                    <div class="property-value" id="compoundSmiles">--</div>
                </div>
            </div>
            
            <!-- Row 2, Col 1: PCA plot -->
            <div class="plot-container">
                <div class="plot-title">PCA Chemical Space</div>
                <div id="pca-plot"></div>
            </div>
            
            <!-- Row 2, Col 2: t-SNE plot -->
            <div class="plot-container">
                <div class="plot-title">t-SNE Chemical Space</div>
                <div id="tsne-plot"></div>
            </div>
        </div>
    </div>"""
    
    def _generate_docking_tab(self, structure_viewer_html: str = "") -> str:
        """Generate docking analysis tab"""
        return f"""
    <div id="docking" class="tab-content">
        <div class="dashboard">
            <div class="plot-container" style="grid-column: span 2;">
                <div class="plot-title">3D Protein-Ligand Visualization</div>
                {structure_viewer_html if structure_viewer_html else '<div style="padding: 2rem; text-align: center; color: #6c757d;">3D visualization not available</div>'}
            </div>
            <div class="plot-container">
                <div class="plot-title">Binding Energy Distribution</div>
                <div id="binding-energy-plot"></div>
            </div>
            <div class="plot-container">
                <div class="plot-title">Top Compounds</div>
                <div class="compound-list" id="compound-list">
                    <!-- Will be populated dynamically -->
                </div>
            </div>
        </div>
    </div>"""
    
    def _generate_structure_panel(self) -> str:
        """Generate floating structure panel (hidden by default)"""
        return ""  # Structure panel is now integrated into the main layout
    
    def _generate_property_options(self, properties: List[str], selected: str = "") -> str:
        """Generate option elements for property selection"""
        options = []
        for prop in properties:
            selected_attr = 'selected' if prop == selected else ''
            options.append(f'<option value="{prop}" {selected_attr}>{prop}</option>')
        return '\n'.join(options)
    
    def _generate_javascript(self) -> str:
        """Generate JavaScript code for interactive functionality"""
        return """
        console.log('=== Enhanced Molecular Analysis Dashboard ===');
        console.log('Loaded data:', data.length, 'compounds');
        console.log('Available properties:', availableProperties);
        console.log('Docking enabled:', dockingEnabled);
        
        // Global variables
        let currentHighlighted = null;
        let plots = {};
        let propertyPlotConfig = {
            x: document.getElementById('x-axis-select').value,
            y: document.getElementById('y-axis-select').value,
            color: document.getElementById('color-select').value
        };
        
        // Color scales
        const colorScales = {};
        
        // Initialize dashboard
        initializeDashboard();
        
        function initializeDashboard() {
            console.log('Initializing dashboard...');
            
            // Create initial plots
            createPropertyPlot();
            createPCAPlot();
            createTSNEPlot();
            
            if (dockingEnabled) {
                populateDockingData();
            }
            
            console.log('Dashboard initialization complete');
        }
        
        function showTab(tabName) {
            // Hide all tab contents
            document.querySelectorAll('.tab-content').forEach(tab => {
                tab.classList.remove('active');
            });
            
            // Remove active class from all tab buttons
            document.querySelectorAll('.tab-button').forEach(btn => {
                btn.classList.remove('active');
            });
            
            // Show selected tab content
            document.getElementById(tabName).classList.add('active');
            
            // Activate corresponding tab button
            event.target.classList.add('active');
        }
        
        function createPropertyPlot() {
            const config = {
                id: 'property-plot',
                title: 'Property Visualization',
                x_prop: propertyPlotConfig.x.toLowerCase(),
                y_prop: propertyPlotConfig.y.toLowerCase(),
                color_prop: propertyPlotConfig.color ? propertyPlotConfig.color.toLowerCase() : null,
                x_label: propertyPlotConfig.x,
                y_label: propertyPlotConfig.y,
                color_label: propertyPlotConfig.color
            };
            
            createPlot(config);
            
            // Update title
            const title = propertyPlotConfig.color ? 
                `${propertyPlotConfig.x} vs ${propertyPlotConfig.y} (colored by ${propertyPlotConfig.color})` :
                `${propertyPlotConfig.x} vs ${propertyPlotConfig.y}`;
            document.getElementById('property-plot-title').textContent = title;
        }
        
        function createPCAPlot() {
            if (!summary.has_pca) {
                document.getElementById('pca-plot').innerHTML = '<p style="text-align: center; color: #6c757d; padding: 2rem;">PCA analysis not available</p>';
                return;
            }
            
            const variance = summary.pca_variance || [0, 0];
            const config = {
                id: 'pca-plot',
                title: 'PCA Chemical Space',
                x_prop: 'pca_x',
                y_prop: 'pca_y',
                color_prop: propertyPlotConfig.color ? propertyPlotConfig.color.toLowerCase() : null,
                x_label: `PC1 (${(variance[0] * 100).toFixed(1)}%)`,
                y_label: `PC2 (${(variance[1] * 100).toFixed(1)}%)`,
                color_label: propertyPlotConfig.color
            };
            
            createPlot(config);
        }
        
        function createTSNEPlot() {
            if (!summary.has_tsne) {
                document.getElementById('tsne-plot').innerHTML = '<p style="text-align: center; color: #6c757d; padding: 2rem;">t-SNE analysis not available</p>';
                return;
            }
            
            const config = {
                id: 'tsne-plot',
                title: 't-SNE Chemical Space',
                x_prop: 'tsne_x',
                y_prop: 'tsne_y',
                color_prop: propertyPlotConfig.color ? propertyPlotConfig.color.toLowerCase() : null,
                x_label: 't-SNE Dimension 1',
                y_label: 't-SNE Dimension 2',
                color_label: propertyPlotConfig.color
            };
            
            createPlot(config);
        }
        
        function createPlot(config) {
            console.log('Creating plot:', config.title);
            
            // Clear existing plot
            d3.select(`#${config.id}`).selectAll("*").remove();
            
            // Check if required data exists
            const validData = data.filter(d => 
                d.hasOwnProperty(config.x_prop) && 
                d.hasOwnProperty(config.y_prop) &&
                !isNaN(d[config.x_prop]) && 
                !isNaN(d[config.y_prop])
            );
            
            if (validData.length === 0) {
                d3.select(`#${config.id}`)
                    .append("div")
                    .style("text-align", "center")
                    .style("color", "#6c757d")
                    .style("padding", "2rem")
                    .text(`No valid data for ${config.x_label} vs ${config.y_label}`);
                return;
            }
            
            // Uniform dimensions for all plots in 2x2 grid
            const margin = {top: 20, right: 120, bottom: 80, left: 80};
            const width = 500 - margin.left - margin.right;
            const height = 350 - margin.top - margin.bottom;
            
            // Create SVG
            const svg = d3.select(`#${config.id}`)
                .append("svg")
                .attr("width", width + margin.left + margin.right)
                .attr("height", height + margin.top + margin.bottom);
                
            const g = svg.append("g")
                .attr("transform", `translate(${margin.left},${margin.top})`);
            
            // Scales
            const xExtent = d3.extent(validData, d => d[config.x_prop]);
            const yExtent = d3.extent(validData, d => d[config.y_prop]);
            
            const xScale = d3.scaleLinear()
                .domain(xExtent)
                .nice()
                .range([0, width]);
                
            const yScale = d3.scaleLinear()
                .domain(yExtent)
                .nice()
                .range([height, 0]);
            
            // Color scale
            let colorScale = null;
            let colorExtent = null;
            if (config.color_prop && validData.some(d => d.hasOwnProperty(config.color_prop) && !isNaN(d[config.color_prop]))) {
                colorExtent = d3.extent(validData, d => d[config.color_prop]);
                colorScale = d3.scaleSequential(d3.interpolateViridis)
                    .domain(colorExtent);
                colorScales[config.id] = colorScale;
            }
            
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
                .attr("transform", `translate(0,${height})`)
                .call(d3.axisBottom(xScale));
                
            g.append("g")
                .attr("class", "axis")
                .call(d3.axisLeft(yScale));
            
            // Add axis labels
            g.append("text")
                .attr("class", "axis-label")
                .attr("transform", `translate(${width/2}, ${height + 45})`)
                .style("text-anchor", "middle")
                .text(config.x_label);
                
            g.append("text")
                .attr("class", "axis-label")
                .attr("transform", "rotate(-90)")
                .attr("y", -margin.left + 20)
                .attr("x", -height / 2)
                .style("text-anchor", "middle")
                .text(config.y_label);
            
            // Add color bar
            if (colorScale && colorExtent) {
                const colorbarWidth = 20;
                const colorbarHeight = 200;
                const colorbarX = width + 20;
                const colorbarY = (height - colorbarHeight) / 2;
                
                // Create gradient
                const gradient = svg.append("defs")
                    .append("linearGradient")
                    .attr("id", `gradient-${config.id}`)
                    .attr("x1", "0%")
                    .attr("y1", "100%")
                    .attr("x2", "0%")
                    .attr("y2", "0%");
                
                const numStops = 20;
                for (let i = 0; i <= numStops; i++) {
                    const offset = i / numStops;
                    const value = colorExtent[0] + offset * (colorExtent[1] - colorExtent[0]);
                    gradient.append("stop")
                        .attr("offset", `${offset * 100}%`)
                        .attr("stop-color", colorScale(value));
                }
                
                const colorbar = g.append("g")
                    .attr("class", "colorbar")
                    .attr("transform", `translate(${colorbarX}, ${colorbarY})`);
                
                colorbar.append("rect")
                    .attr("width", colorbarWidth)
                    .attr("height", colorbarHeight)
                    .style("fill", `url(#gradient-${config.id})`)
                    .style("stroke", "#dee2e6")
                    .style("stroke-width", 1);
                
                // Color bar scale
                const colorbarScale = d3.scaleLinear()
                    .domain(colorExtent)
                    .range([colorbarHeight, 0]);
                
                colorbar.append("g")
                    .attr("transform", `translate(${colorbarWidth}, 0)`)
                    .call(d3.axisRight(colorbarScale).ticks(5).tickSize(3));
                
                // Color bar title
                colorbar.append("text")
                    .attr("class", "colorbar-title")
                    .attr("transform", "rotate(-90)")
                    .attr("x", -colorbarHeight / 2)
                    .attr("y", -10)
                    .style("text-anchor", "middle")
                    .text(config.color_label);
            }
            
            // Add dots
            const dots = g.selectAll(".dot")
                .data(validData)
                .enter().append("circle")
                .attr("class", "dot")
                .attr("cx", d => xScale(d[config.x_prop]))
                .attr("cy", d => yScale(d[config.y_prop]))
                .attr("r", 3)
                .attr("fill", d => colorScale ? colorScale(d[config.color_prop]) : "#3498db")
                .on("mouseover", function(event, d) {
                    highlightCompound(d.id);
                    showTooltip(event, d, config);
                    showStructure(d);
                })
                .on("mouseout", function() {
                    hideTooltip();
                })
                .on("click", function(event, d) {
                    showStructure(d);
                });
            
            // Store plot reference
            plots[config.id] = {
                svg: svg,
                g: g,
                dots: dots,
                config: config,
                validData: validData
            };
            
            console.log(`Plot ${config.id} created with ${validData.length} points`);
        }
        
        function updatePropertyPlot() {
            propertyPlotConfig = {
                x: document.getElementById('x-axis-select').value,
                y: document.getElementById('y-axis-select').value,
                color: document.getElementById('color-select').value || null
            };
            
            createPropertyPlot();
            createPCAPlot();  // Update PCA plot color
            createTSNEPlot(); // Update t-SNE plot color
        }
        
        function highlightCompound(compoundId) {
            // Remove previous highlights
            if (currentHighlighted !== null) {
                Object.values(plots).forEach(plot => {
                    if (plot.dots) {
                        plot.dots.classed("highlighted", false);
                    }
                });
            }
            
            // Add new highlights
            Object.values(plots).forEach(plot => {
                if (plot.dots) {
                    plot.dots.classed("highlighted", d => d.id === compoundId);
                }
            });
            
            currentHighlighted = compoundId;
        }
        
        function showTooltip(event, compound, config) {
            const tooltip = d3.select("#tooltip");
            const colorInfo = config.color_prop && compound[config.color_prop] !== undefined ? 
                `<br>${config.color_label}: ${compound[config.color_prop].toFixed(2)}` : '';
            
            tooltip.style("opacity", 1)
                .html(`
                    <strong>${compound.title}</strong><br>
                    ${config.x_label}: ${compound[config.x_prop].toFixed(2)}<br>
                    ${config.y_label}: ${compound[config.y_prop].toFixed(2)}${colorInfo}
                `)
                .style("left", (event.pageX + 10) + "px")
                .style("top", (event.pageY - 10) + "px");
        }
        
        function hideTooltip() {
            d3.select("#tooltip").style("opacity", 0);
        }
        
        function showStructure(compound) {
            console.log('Showing structure for:', compound.title);
            
            document.getElementById('compoundName').textContent = compound.title;
            
            // Update structure image
            const structureImage = document.getElementById('structureImage');
            if (compound.image) {
                structureImage.innerHTML = `<img src="${compound.image}" alt="Molecular Structure" />`;
            } else {
                structureImage.innerHTML = '<div class="placeholder">No structure available</div>';
            }
            
            // Update properties grid
            const propertiesGrid = document.getElementById('propertiesGrid');
            const propertyBoxes = [
                {label: 'MW', value: compound.mw ? `${compound.mw.toFixed(1)}` : '--'},
                {label: 'LogP', value: compound.logp ? compound.logp.toFixed(2) : '--'},
                {label: 'TPSA', value: compound.tpsa ? `${compound.tpsa.toFixed(0)}` : '--'},
                {label: 'HBA', value: compound.hba !== undefined ? compound.hba : '--'},
                {label: 'HBD', value: compound.hbd !== undefined ? compound.hbd : '--'},
                {label: 'RotBonds', value: compound.rotbonds !== undefined ? compound.rotbonds : '--'},
                {label: 'Rings', value: compound.numrings !== undefined ? compound.numrings : '--'},
                {label: 'QED', value: compound.qed ? compound.qed.toFixed(3) : '--'},
                {label: 'SA Score', value: compound.sascore ? compound.sascore.toFixed(2) : '--'}
            ];
            
            // Add docking results if available
            if (dockingEnabled && compound.binding_energy !== undefined && !isNaN(compound.binding_energy)) {
                propertyBoxes.push({
                    label: 'Binding Energy', 
                    value: `${compound.binding_energy.toFixed(2)} kcal/mol`
                });
            }
            
            propertiesGrid.innerHTML = propertyBoxes.map(prop => `
                <div class="property-box">
                    <div class="property-label">${prop.label}</div>
                    <div class="property-value">${prop.value}</div>
                </div>
            `).join('');
            
            // Update SMILES
            document.getElementById('compoundSmiles').textContent = compound.smiles;
        }
        
        function populateDockingData() {
            if (!dockingEnabled || !dockingData.length) return;
            
            console.log('Populating docking data...');
            
            // Populate compound list
            const compoundList = document.getElementById('compound-list');
            compoundList.innerHTML = dockingData
                .sort((a, b) => a.binding_energy - b.binding_energy)
                .map(compound => {
                    const energyClass = compound.binding_energy < -8 ? 'energy-good' :
                                      compound.binding_energy < -6 ? 'energy-moderate' : 'energy-poor';
                    
                    return `
                        <div class="compound-item" onclick="selectDockingCompound(${compound.compound_id})">
                            <div class="compound-name">${compound.compound_name}</div>
                            <div class="binding-energy ${energyClass}">
                                ${compound.binding_energy.toFixed(2)} kcal/mol
                            </div>
                        </div>
                    `;
                }).join('');
                
            // Create binding energy distribution plot
            createBindingEnergyPlot();
        }
        
        function selectDockingCompound(compoundId) {
            // Update selector
            document.getElementById('compound-select').value = compoundId;
            
            // Highlight in list
            document.querySelectorAll('.compound-item').forEach(item => {
                item.classList.remove('selected');
            });
            event.target.closest('.compound-item').classList.add('selected');
            
            // Load pose
            loadDockingPose();
        }
        
        function createBindingEnergyPlot() {
            if (!dockingEnabled || !dockingData.length) return;
            
            console.log('Creating binding energy distribution plot...');
            
            // Clear existing plot
            d3.select("#binding-energy-plot").selectAll("*").remove();
            
            const margin = {top: 20, right: 20, bottom: 40, left: 60};
            const width = 400 - margin.left - margin.right;
            const height = 300 - margin.top - margin.bottom;
            
            const svg = d3.select("#binding-energy-plot")
                .append("svg")
                .attr("width", width + margin.left + margin.right)
                .attr("height", height + margin.top + margin.bottom);
                
            const g = svg.append("g")
                .attr("transform", `translate(${margin.left},${margin.top})`);
            
            // Prepare data
            const energies = dockingData.map(d => d.binding_energy).sort((a, b) => a - b);
            const bins = d3.bin().thresholds(10)(energies);
            
            const xScale = d3.scaleLinear()
                .domain(d3.extent(energies))
                .range([0, width]);
                
            const yScale = d3.scaleLinear()
                .domain([0, d3.max(bins, d => d.length)])
                .range([height, 0]);
            
            // Add bars
            g.selectAll(".bar")
                .data(bins)
                .enter().append("rect")
                .attr("class", "bar")
                .attr("x", d => xScale(d.x0))
                .attr("width", d => Math.max(0, xScale(d.x1) - xScale(d.x0) - 1))
                .attr("y", d => yScale(d.length))
                .attr("height", d => height - yScale(d.length))
                .attr("fill", "#3498db")
                .attr("opacity", 0.7);
            
            // Add axes
            g.append("g")
                .attr("transform", `translate(0,${height})`)
                .call(d3.axisBottom(xScale));
                
            g.append("g")
                .call(d3.axisLeft(yScale));
            
            // Add labels
            g.append("text")
                .attr("transform", `translate(${width/2}, ${height + 35})`)
                .style("text-anchor", "middle")
                .style("font-size", "12px")
                .text("Binding Energy (kcal/mol)");
                
            g.append("text")
                .attr("transform", "rotate(-90)")
                .attr("y", -40)
                .attr("x", -height / 2)
                .style("text-anchor", "middle")
                .style("font-size", "12px")
                .text("Frequency");
        }
        
        // Make functions globally available
        window.showTab = showTab;
        window.updatePropertyPlot = updatePropertyPlot;
        window.selectDockingCompound = selectDockingCompound;
        """
    
    def save_dashboard(self, html_content: str, output_file: str) -> None:
        """
        Save HTML dashboard to file
        
        Args:
            html_content: Generated HTML content
            output_file: Output file path
        """
        try:
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            
            with open(output_path, 'w', encoding='utf-8') as f:
                f.write(html_content)
            
            self.logger.info(f"Dashboard saved to: {output_file}")
            
        except Exception as e:
            self.logger.error(f"Error saving dashboard: {e}")
            raise