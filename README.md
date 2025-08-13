# Molecular View Dashboard

A comprehensive interactive molecular visualization and docking analysis platform that transforms chemical compound data into insightful visual dashboards.

## What is this codebase for?

This tool provides **end-to-end molecular analysis** for chemists, biologists, and computational researchers:

- **Analyze molecular properties** from SMILES strings
- **Visualize chemical space** using PCA and t-SNE
- **Perform molecular docking** with AutoDock Vina + MGLTools
- **Generate interactive dashboards** with synchronized 3D visualizations

## Key Features

🧪 **Property Analysis**: Calculate 9+ molecular descriptors (MW, LogP, TPSA, QED, etc.)  
🧬 **Chemical Space**: PCA and t-SNE dimensionality reduction (requires ≥10 compounds)  
⚗️ **Molecular Docking**: AutoDock Vina integration with MGLTools preparation  
🎯 **3D Visualization**: Interactive protein-ligand viewer with cartoon representation  
📊 **Dashboard**: Professional HTML interface with synchronized highlighting  
🔧 **Configuration**: Flexible YAML-based setup  

## Installation

### Option 1: Conda (Recommended)

```bash
# Create environment from file
conda env create -f environment.yaml
conda activate mol_view_dashboard

# Install the package
pip install -e .
```

### Option 2: Manual Setup

```bash
# Create conda environment
conda create -n mol_view_dashboard python=3.11
conda activate mol_view_dashboard

# Install dependencies
conda install -c conda-forge rdkit pandas numpy scikit-learn pyyaml
pip install -r requirements.txt

# Install package
pip install -e .
```

### External Tools (for docking)

- **AutoDock Vina**: [Download here](https://vina.scripps.edu/)
- **MGLTools**: [Download here](http://mgltools.scripps.edu/)

## Quick Start

**Run the example:**

```bash
python run_analysis.py examples/sample_compounds.csv --config examples/config.yaml
```

**View results:** Open `examples/molecular_docking_mgltools.html` in your browser

## Usage

### Basic Analysis (No Docking)

```bash
python run_analysis.py your_molecules.csv --config config.yaml
```

### With Molecular Docking

```bash
python run_analysis.py compounds.csv --config docking_config.yaml --verbose
```

## Configuration Guide

### Basic Configuration (`config.yaml`)

```yaml
input:
  csv_file: "molecules.csv"
  smiles_column: "SMILES"
  name_column: "Name"
  property_columns: []  # Optional: additional CSV columns

analysis:
  calculate_descriptors: true
  chemical_space:
    pca:
      enabled: true
      n_components: 2
    tsne:
      enabled: true      # Requires ≥10 compounds
      perplexity: 30

visualization:
  output_file: "dashboard.html"
  title: "Molecular Analysis Dashboard"
  property_plots:
    primary_plot:
      x_axis: "MW"
      y_axis: "LogP" 
      color_by: "docking_score"  # or "TPSA", "QED", etc.
```

### Docking Configuration

```yaml
docking:
  enabled: true
  protein_pdb: "examples/4csv.pdb"
  vina_executable: "/path/to/vina" # the path to vina excecutable
  mgltools_path: "/path/to/mgltools/installed" # where teh mgltools was installed
  
  binding_site:
    center: [-21.23, 8.03, -0.93]  # Binding pocket center
    size: [20.0, 20.0, 20.0]       # Search box size

  parameters:
    exhaustiveness: 8
    num_modes: 9
    energy_range: 3
```

## Examples Directory

The `examples/` directory contains:

### 📁 `sample_compounds.csv`
```csv
SMILES,Name
ClC1=CC(Cl)=C(OC)C=C1NC3=C(C#N)C=NC2=CC(OCCCN4CCN(C)CC4)=C(OC)C=C23,Bosutinib
CC1=NC(NC2=NC=C(C(NC3=C(C=CC=C3Cl)C)=O)S2)=CC(N4CCN(CC4)CCO)=N1,Dasatinib
...
```
**Purpose**: Example kinase inhibitors for testing the complete workflow

### 📁 `4csv.pdb`
**Purpose**: Crystal structure of a kinase protein for docking studies

### 📁 `config.yaml`
**Purpose**: Complete configuration example with docking enabled

## What the Dashboard Shows

### 🔬 **Main Analysis Tab**
- **Property Plot**: Interactive scatter plot with configurable axes
- **PCA Plot**: Chemical space visualization showing compound clusters
- **t-SNE Plot**: Non-linear dimensionality reduction (if ≥10 compounds)

### ⚗️ **Docking Results Tab** (if enabled)
- **3D Protein-Ligand Viewer**: Interactive molecular visualization
- **Docking Score Distribution**: Histogram of binding energies
- **Top Compounds**: Ranked list by docking score
- **Structure Selector**: Switch between different ligand poses

### 📊 **Features**
- **Synchronized highlighting**: Click any point to highlight across all plots
- **Compound details**: Hover for molecular properties
- **Professional styling**: Color scales, legends, and responsive design

## Output Files

- **HTML Dashboard**: Interactive visualization (main output)
- **CSV Results**: Processed data with calculated descriptors and docking scores
- **Structure Files**: PDBQT files for protein and ligands (in `docking_results/`)

## Requirements

- **Python**: 3.11+
- **Core**: pandas, numpy, scikit-learn, pyyaml
- **Chemistry**: RDKit (conda recommended)
- **Optional**: AutoDock Vina, MGLTools (for docking)

## Troubleshooting

**"t-SNE skipped"**: Normal for <10 compounds  
**"RDKit not found"**: Install via `conda install -c conda-forge rdkit`  
**Docking fails**: Check vina_executable and mgltools_path in config  

## License

MIT License