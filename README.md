# Molecular Visualization Dashboard

An open-source repository for interactive visualization of small molecules with property annotations. This toolkit provides comprehensive molecular analysis capabilities including chemical space exploration, property visualization, and interactive dashboards.

## 🚀 Features

### Core Capabilities
- **SDF to CSV Conversion**: Convert SDF files to CSV format with canonicalized SMILES
- **Property Visualization**: Create interactive scatter plots of molecular properties
- **Chemical Space Analysis**: Explore molecular diversity using PCA and t-SNE
- **Interactive Dashboards**: Generate HTML dashboards with synchronized highlighting
- **Structure Display**: View 2D molecular structures on hover/click
- **Multi-dataset Support**: Analyze and compare multiple datasets simultaneously

### Visualization Types
- **Property Correlation Plots**: MW vs LogP, TPSA vs QED, HBA vs HBD
- **Chemical Space Mapping**: t-SNE and PCA molecular space visualization
- **Custom Property Analysis**: Use your own molecular properties for visualization
- **Multi-source Comparison**: Color-code and compare compounds from different sources

## 📦 Installation

### Prerequisites
- Python 3.7+
- RDKit (chemical informatics toolkit)

### Option 1: Using Conda (Recommended)
```bash
# Create conda environment with RDKit
conda create -n molview python=3.8
conda activate molview
conda install -c conda-forge rdkit

# Install other dependencies
pip install pandas numpy scikit-learn
```

### Option 2: Using pip
```bash
# Create virtual environment
python -m venv molview
source molview/bin/activate  # On Windows: molview\Scripts\activate

# Install all dependencies
pip install -r requirements.txt
```

### Clone Repository
```bash
git clone https://github.com/your-username/mol_view_dashboard.git
cd mol_view_dashboard
```

## 📁 Repository Structure

```
mol_view_dashboard/
├── scripts/
│   ├── sdf_to_csv_converter.py      # Convert SDF files to CSV
│   ├── molecular_property_visualizer.py  # Single dataset analysis
│   └── chemical_space_analyzer.py   # Multi-dataset comparison
├── requirements.txt                 # Python dependencies
├── LICENSE                         # License file
└── README.md                       # This file
```

## 🔧 Usage

### 1. Convert SDF to CSV

Convert single SDF file:
```bash
python scripts/sdf_to_csv_converter.py input.sdf --output molecules.csv
```

Convert multiple SDF files from directory:
```bash
python scripts/sdf_to_csv_converter.py /path/to/sdf_directory/
```

**Features:**
- Extracts canonicalized SMILES strings
- Preserves all SDF properties
- Handles multiple files with source tracking
- Generates both individual and merged CSV files

### 2. Single Dataset Analysis

Analyze molecular properties from CSV file:
```bash
python scripts/molecular_property_visualizer.py data.csv --output dashboard.html
```

With custom properties:
```bash
python scripts/molecular_property_visualizer.py data.csv \
    --property_from_input "Docking_Score" "Binding_Affinity" \
    --output custom_dashboard.html
```

**Generated Plots:**
- Molecular Weight vs LogP
- TPSA vs QED (Drug-likeness)
- H-Bond Acceptors vs H-Bond Donors
- t-SNE Molecular Space
- Custom property correlations

### 3. Multi-Dataset Comparison

Compare multiple datasets with chemical space analysis:
```bash
python scripts/chemical_space_analyzer.py \
    --csv_files dataset1.csv dataset2.csv dataset3.csv \
    --labels "ChEMBL" "PubChem" "Custom" \
    --output comparison_dashboard.html
```

**Features:**
- PCA analysis with explained variance
- t-SNE chemical space visualization
- Source-based color coding
- Interactive compound highlighting
- Cross-plot synchronization

## 📊 Dashboard Features

### Interactive Elements
- **Hover Effects**: View molecular structures and properties on hover
- **Synchronized Highlighting**: Compound highlighting across all plots
- **Color Coding**: Third dimension visualization using color scales
- **Property Panel**: Detailed molecular information display
- **Responsive Design**: Works across different screen sizes

### Molecular Properties Calculated
- **Basic**: Molecular Weight, LogP, TPSA
- **Drug-likeness**: QED, Synthetic Accessibility Score
- **Topology**: H-bond acceptors/donors, Rotatable bonds
- **Structural**: Aromatic rings, Heavy atoms
- **Custom**: User-provided properties from CSV

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