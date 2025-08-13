"""
Molecular View Dashboard

An interactive molecular visualization and docking analysis dashboard that provides:
- Interactive property plots with configurable axes
- Chemical space visualization (PCA + t-SNE)
- Synchronized highlighting across all plots
- Molecular structure display with 3D visualization
- Professional styling with color bars
- Molecular docking results and pose visualization
"""

__version__ = "1.0.0"
__author__ = "Your Name"
__email__ = "your.email@example.com"

# Import main classes for easy access
from .config_manager import ConfigManager, ConfigurationError
from .molecular_data_processor import MolecularDataProcessor
from .docking_wrapper import VinaDockingWrapper, DockingError
from .dashboard_generator import DashboardGenerator
from .structure_viewer import StructureViewer

__all__ = [
    "ConfigManager",
    "ConfigurationError", 
    "MolecularDataProcessor",
    "VinaDockingWrapper",
    "DockingError",
    "DashboardGenerator", 
    "StructureViewer",
]