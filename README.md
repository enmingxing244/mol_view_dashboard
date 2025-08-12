# Molecular Visualization and Analysis Tool

A comprehensive Python-based tool for molecular property analysis, chemical space visualization, and optional molecular docking. This enhanced toolkit creates interactive HTML dashboards with synchronized highlighting across multiple plot types including property visualizations, PCA, and t-SNE chemical space analysis.

## 🚀 Features

### Core Analysis Capabilities
- **Molecular Property Calculation**: Comprehensive descriptor calculation (MW, LogP, TPSA, QED, SA Score, etc.)
- **Chemical Space Analysis**: PCA and t-SNE visualization of molecular fingerprints
- **Interactive Property Plots**: User-configurable scatter plots with selectable X/Y axes and color coding
- **Synchronized Highlighting**: Mouse interactions highlight compounds across all plots simultaneously

### Advanced Features
- **Molecular Docking**: Optional AutoDock Vina integration for protein-ligand docking
- **Professional Styling**: Clean white backgrounds, color bars, and scientific visualization standards
- **YAML Configuration**: Flexible configuration system for customizing all analysis parameters
- **Export Capabilities**: Export processed data and analysis results

### Visualization Features
- **Interactive HTML Dashboard**: Modern web-based interface with tabs and controls
- **Structure Display**: Molecular structure visualization on hover/click
- **Docking Pose Viewer**: 3D visualization of docking poses (when docking enabled)
- **Responsive Design**: Works on desktop and mobile devices

## 📦 Installation

### Prerequisites
- Python 3.8 or higher
- RDKit (chemistry toolkit)
- AutoDock Vina (optional, for docking)

### Install Dependencies

```bash
# Using conda (recommended for RDKit)
conda create -n mol_analysis python=3.9
conda activate mol_analysis
conda install -c conda-forge rdkit pandas numpy scikit-learn pyyaml

# Or using pip
pip install -r requirements.txt
```

### Install AutoDock Vina (Optional)
For molecular docking functionality:

```bash
# On Ubuntu/Debian
sudo apt-get install autodock-vina

# On macOS with Homebrew
brew install autodock-vina

# Or download from: https://vina.scripps.edu/downloads/
```

## 📁 Repository Structure

```
mol_view_dashboard/
├── src/                              # Enhanced analysis modules
│   ├── config_manager.py            # YAML configuration handling
│   ├── molecular_data_processor.py  # Data loading and analysis
│   ├── docking_wrapper.py           # AutoDock Vina integration
│   └── dashboard_generator.py       # HTML dashboard creation
├── scripts/                          # Legacy scripts (still functional)
│   ├── sdf_to_csv_converter.py      # Convert SDF files to CSV
│   ├── molecular_property_visualizer.py  # Single dataset analysis
│   └── chemical_space_analyzer.py   # Multi-dataset comparison
├── run_analysis.py                   # Main entry point (NEW)
├── config_template.yaml             # Configuration template
├── requirements.txt                 # Python dependencies
└── README.md                        # This file
```

## 🔧 Usage

## Quick Start

### Basic Usage (NEW Enhanced Interface)

```bash
# Analyze a CSV file with SMILES
python run_analysis.py compounds.csv

# Specify output file
python run_analysis.py compounds.csv --output my_analysis.html

# Include specific property columns
python run_analysis.py compounds.csv --properties MW LogP TPSA QED
```

### Advanced Usage with Configuration

```bash
# Generate sample configuration file
python run_analysis.py --generate-config

# Use custom configuration
python run_analysis.py compounds.csv --config my_config.yaml

# Enable molecular docking
python run_analysis.py compounds.csv --protein protein.pdb --vina-config vina.conf
```

## Configuration

### YAML Configuration File
The tool uses YAML configuration files for detailed customization. Generate a template:

```bash
python run_analysis.py --generate-config
```

### Example Configuration Sections

#### Input Configuration
```yaml
input:
  csv_file: "data/compounds.csv"
  smiles_column: "SMILES"
  name_column: "Name"
  property_columns:
    - "MW"
    - "LogP"
    - "TPSA"
    - "IC50"
```

#### Docking Configuration
```yaml
docking:
  enabled: true
  protein_pdb: "data/protein.pdb"
  vina_config: "data/vina.conf"
  parameters:
    exhaustiveness: 8
    num_modes: 9
```

## Legacy Scripts (Still Functional)

### 1. Convert SDF to CSV

Convert single SDF file:
```bash
python scripts/sdf_to_csv_converter.py input.sdf --output molecules.csv
```

Convert multiple SDF files from directory:
```bash
python scripts/sdf_to_csv_converter.py /path/to/sdf_directory/
```

### 2. Single Dataset Analysis (Legacy)

```bash
python scripts/molecular_property_visualizer.py data.csv --output dashboard.html
```

### 3. Multi-Dataset Comparison (Legacy)

```bash
python scripts/chemical_space_analyzer.py \
    --csv_files dataset1.csv dataset2.csv dataset3.csv \
    --labels "ChEMBL" "PubChem" "Custom" \
    --output comparison_dashboard.html
```

## 📊 Enhanced Dashboard Features

### Main Analysis Tab
- **Property Visualization**: Configurable scatter plot with dropdown selectors for X-axis, Y-axis, and color coding
- **PCA Plot**: Principal component analysis of molecular fingerprints with explained variance
- **t-SNE Plot**: t-distributed stochastic neighbor embedding for chemical space visualization
- **Structure Panel**: Molecular structure display with comprehensive property details

### Docking Tab (if enabled)
- **3D Pose Viewer**: Interactive visualization of docking poses
- **Binding Energy Plot**: Distribution of binding energies
- **Compound Rankings**: List of compounds sorted by binding affinity

### Interactive Features
- **Synchronized Highlighting**: Hover over any point to highlight the same compound across all plots
- **Structure Display**: Click or hover to view molecular structures and properties
- **Configurable Plots**: Change X/Y axes and color coding on the fly
- **Professional Styling**: Clean white backgrounds, color bars, and scientific styling

### Molecular Properties Calculated
- **Basic**: Molecular Weight, LogP, TPSA
- **Drug-likeness**: QED, Synthetic Accessibility Score
- **Topology**: H-bond acceptors/donors, Rotatable bonds
- **Structural**: Aromatic rings, Heavy atoms, Fraction Csp3
- **Custom**: User-provided properties from CSV
- **Docking**: Binding energies and pose files (if enabled)

## 💡 Example Workflows

### Workflow 1: Drug Discovery Pipeline
```bash
# 1. Convert docking results from SDF to CSV
python scripts/sdf_to_csv_converter.py docking_results.sdf

# 2. Visualize with docking scores
python scripts/molecular_property_visualizer.py docking_results.csv \
    --property_from_input "docking_score" "binding_energy"
```

### Workflow 2: Chemical Space Comparison
```bash
# Compare different compound libraries
python scripts/chemical_space_analyzer.py \
    --csv_files approved_drugs.csv experimental.csv natural_products.csv \
    --labels "FDA_Approved" "Experimental" "Natural_Products"
```

### Workflow 3: SAR Analysis
```bash
# Structure-Activity Relationship analysis
python scripts/molecular_property_visualizer.py sar_data.csv \
    --property_from_input "IC50" "LogIC50" "Activity_Class"
```

## 📋 Input Data Requirements

### CSV Format
Required columns:
- `SMILES`: Canonical SMILES strings
- `Name` (optional): Compound names/IDs

Optional columns:
- Any numerical properties for visualization
- Custom property names (specified via command line)

Example CSV structure:
```csv
SMILES,Name,MW,LogP,Activity
CCO,Ethanol,46.07,-0.31,0.5
CC(C)O,Isopropanol,60.10,0.05,0.8
```

### SDF Format
- Standard MDL SDF format
- Properties stored as SDF tags
- Multiple molecules per file supported

## 🎨 Customization

### Color Schemes
- Default: Viridis colormap for continuous properties
- Source-based: Distinct colors for different datasets
- Customizable through CSS in generated HTML

### Plot Configuration
- Adjustable plot dimensions and margins
- Configurable axis labels and titles
- Flexible property selection for X, Y, and color dimensions

## 🔬 Technical Details

### Molecular Descriptors
- **Fingerprints**: Morgan fingerprints (radius=2, 2048 bits)
- **Dimensionality Reduction**: PCA and t-SNE with optimal parameters
- **Standardization**: Feature scaling for consistent analysis

### Performance
- Optimized for datasets up to 10,000 compounds
- Efficient molecular structure generation
- Memory-conscious fingerprint calculations
- Progress tracking for large datasets

## 🤝 Contributing

Contributions are welcome! Areas for improvement:
- Additional visualization types
- New molecular descriptors
- Performance optimizations
- Export capabilities (PNG, SVG, PDF)

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🐛 Troubleshooting

### Common Issues

**RDKit Installation**:
- Use conda for easiest RDKit installation
- Ensure correct environment activation
- Check RDKit version compatibility

**Memory Issues**:
- Reduce dataset size for large files
- Use chunked processing for >10k molecules
- Consider increasing system memory

**Browser Compatibility**:
- Modern browsers required (Chrome, Firefox, Safari, Edge)
- JavaScript must be enabled
- Local file access may require server

## 📞 Support

For questions, issues, or feature requests, please open an issue in the GitHub repository.

---

**Keywords**: cheminformatics, molecular visualization, chemical space, drug discovery, RDKit, interactive dashboards, SMILES, SDF, t-SNE, PCA