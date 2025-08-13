#!/usr/bin/env python3
"""
Configuration Manager for Molecular Analysis Pipeline
Handles YAML configuration loading, validation, and default settings
"""

import yaml
import os
from pathlib import Path
from typing import Dict, Any, Optional
import logging

class ConfigurationError(Exception):
    """Custom exception for configuration-related errors"""
    pass

class ConfigManager:
    """Manages configuration loading and validation for molecular analysis"""
    
    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize configuration manager
        
        Args:
            config_path: Path to YAML configuration file
        """
        self.config = None
        self.config_path = config_path
        self.logger = logging.getLogger(__name__)
        
        if config_path:
            self.load_config(config_path)
        else:
            self.load_default_config()
    
    def load_config(self, config_path: str) -> Dict[str, Any]:
        """
        Load configuration from YAML file
        
        Args:
            config_path: Path to YAML configuration file
            
        Returns:
            Loaded configuration dictionary
            
        Raises:
            ConfigurationError: If configuration file cannot be loaded or is invalid
        """
        try:
            config_file = Path(config_path)
            if not config_file.exists():
                raise ConfigurationError(f"Configuration file not found: {config_path}")
            
            with open(config_file, 'r', encoding='utf-8') as f:
                self.config = yaml.safe_load(f)
            
            self.logger.info(f"Configuration loaded from: {config_path}")
            self._validate_config()
            return self.config
            
        except yaml.YAMLError as e:
            raise ConfigurationError(f"Error parsing YAML configuration: {e}")
        except Exception as e:
            raise ConfigurationError(f"Error loading configuration: {e}")
    
    def load_default_config(self) -> Dict[str, Any]:
        """Load default configuration settings"""
        self.config = {
            'input': {
                'csv_file': None,
                'smiles_column': 'SMILES',
                'name_column': 'Name',
                'property_columns': []
            },
            'analysis': {
                'calculate_descriptors': True,
                'chemical_space': {
                    'pca': {'enabled': True, 'n_components': 2},
                    'tsne': {'enabled': True, 'n_components': 2, 'perplexity': 30, 'random_state': 42}
                },
                'fingerprints': {
                    'type': 'morgan',
                    'radius': 2,
                    'n_bits': 2048
                }
            },
            'docking': {
                'enabled': False,
                'protein_pdb': None,
                'vina_config': None,
                'parameters': {
                    'exhaustiveness': 8,
                    'num_modes': 9,
                    'energy_range': 3
                },
                'output_dir': 'docking_results'
            },
            'visualization': {
                'output_file': 'molecular_analysis_dashboard.html',
                'title': 'Molecular Visualization and Analysis',
                'property_plots': {
                    'primary_plot': {
                        'x_axis': 'MW',
                        'y_axis': 'LogP',
                        'color_by': 'QED'
                    }
                },
                'chemical_space_plots': {
                    'pca': {'title': 'PCA Chemical Space', 'show_variance': True},
                    'tsne': {'title': 't-SNE Chemical Space'}
                },
                'docking_plots': {
                    'pose_viewer': 'ngl',
                    'show_interactions': True
                },
                'style': {
                    'color_scheme': 'viridis',
                    'background': 'white',
                    'plot_style': 'clean',
                    'show_colorbars': True,
                    'synchronized_highlighting': True
                }
            },
            'export': {
                'formats': ['html'],
                'export_data': True,
                'export_file': 'analysis_results.csv',
                'export_docking': True
            },
            'performance': {
                'n_jobs': -1,
                'chunk_size': 1000,
                'cache_calculations': True,
                'cache_dir': '.cache'
            },
            'advanced': {
                'custom_descriptors': {'enabled': False},
                'external_tools': {
                    'rdkit_version': 'latest',
                    'sklearn_random_state': 42
                },
                'debug_mode': False,
                'log_level': 'INFO',
                'log_file': 'analysis.log'
            }
        }
        
        self.logger.info("Default configuration loaded")
        return self.config
    
    def _validate_config(self) -> None:
        """
        Validate configuration structure and required fields
        
        Raises:
            ConfigurationError: If configuration is invalid
        """
        if not self.config:
            raise ConfigurationError("No configuration loaded")
        
        # Check required sections
        required_sections = ['input', 'analysis', 'visualization']
        for section in required_sections:
            if section not in self.config:
                raise ConfigurationError(f"Missing required configuration section: {section}")
        
        # Validate input section
        input_config = self.config.get('input', {})
        if not input_config.get('csv_file'):
            self.logger.warning("No CSV file specified in configuration")
        
        # Validate file paths if they exist
        csv_file = input_config.get('csv_file')
        if csv_file and not Path(csv_file).exists():
            self.logger.warning(f"CSV file not found: {csv_file}")
        
        # Validate docking configuration if enabled
        docking_config = self.config.get('docking', {})
        if docking_config.get('enabled', False):
            protein_pdb = docking_config.get('protein_pdb')
            vina_config = docking_config.get('vina_config')
            binding_site = docking_config.get('binding_site')
            
            if not protein_pdb:
                raise ConfigurationError("Docking enabled but no protein PDB file specified")
            
            # Check that either vina_config OR binding_site is specified
            if not vina_config and not binding_site:
                raise ConfigurationError("Docking enabled but neither 'vina_config' file nor 'binding_site' parameters specified")
            
            # Validate binding site parameters if using inline configuration
            if binding_site:
                center = binding_site.get('center')
                size = binding_site.get('size')
                if not center or len(center) != 3:
                    raise ConfigurationError("Binding site 'center' must be a list of 3 coordinates [x, y, z]")
                if not size or len(size) != 3:
                    raise ConfigurationError("Binding site 'size' must be a list of 3 dimensions [x, y, z]")
            
            # Check file existence
            if not Path(protein_pdb).exists():
                self.logger.warning(f"Protein PDB file not found: {protein_pdb}")
            if vina_config and not Path(vina_config).exists():
                self.logger.warning(f"Vina configuration file not found: {vina_config}")
        
        self.logger.info("Configuration validation completed")
    
    def get(self, key_path: str, default: Any = None) -> Any:
        """
        Get configuration value using dot notation
        
        Args:
            key_path: Dot-separated key path (e.g., 'analysis.chemical_space.pca.enabled')
            default: Default value if key not found
            
        Returns:
            Configuration value or default
        """
        if not self.config:
            return default
        
        keys = key_path.split('.')
        value = self.config
        
        for key in keys:
            if isinstance(value, dict) and key in value:
                value = value[key]
            else:
                return default
        
        return value
    
    def set(self, key_path: str, value: Any) -> None:
        """
        Set configuration value using dot notation
        
        Args:
            key_path: Dot-separated key path
            value: Value to set
        """
        if not self.config:
            self.config = {}
        
        keys = key_path.split('.')
        current = self.config
        
        for key in keys[:-1]:
            if key not in current:
                current[key] = {}
            current = current[key]
        
        current[keys[-1]] = value
    
    def update_from_args(self, args) -> None:
        """
        Update configuration from command line arguments
        
        Args:
            args: Parsed command line arguments
        """
        if hasattr(args, 'input') and args.input:
            self.set('input.csv_file', args.input)
        
        if hasattr(args, 'output') and args.output:
            self.set('visualization.output_file', args.output)
        
        if hasattr(args, 'protein_pdb') and args.protein_pdb:
            self.set('docking.protein_pdb', args.protein_pdb)
            self.set('docking.enabled', True)
        
        if hasattr(args, 'vina_config') and args.vina_config:
            self.set('docking.vina_config', args.vina_config)
            self.set('docking.enabled', True)
        
        if hasattr(args, 'no_docking') and args.no_docking:
            self.set('docking.enabled', False)
        
        if hasattr(args, 'property_columns') and args.property_columns:
            self.set('input.property_columns', args.property_columns)
    
    def save_config(self, output_path: str) -> None:
        """
        Save current configuration to YAML file
        
        Args:
            output_path: Path to save configuration file
        """
        try:
            with open(output_path, 'w', encoding='utf-8') as f:
                yaml.dump(self.config, f, default_flow_style=False, indent=2)
            self.logger.info(f"Configuration saved to: {output_path}")
        except Exception as e:
            raise ConfigurationError(f"Error saving configuration: {e}")
    
    def get_full_config(self) -> Dict[str, Any]:
        """Get the complete configuration dictionary"""
        return self.config or {}
    
    def is_docking_enabled(self) -> bool:
        """Check if molecular docking is enabled"""
        return self.get('docking.enabled', False)
    
    def get_input_file(self) -> Optional[str]:
        """Get the input CSV file path"""
        return self.get('input.csv_file')
    
    def get_output_file(self) -> str:
        """Get the output HTML file path"""
        return self.get('visualization.output_file', 'molecular_analysis_dashboard.html')
    
    def get_property_columns(self) -> list:
        """Get the list of property columns to analyze"""
        return self.get('input.property_columns', [])