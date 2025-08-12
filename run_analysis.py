#!/usr/bin/env python3
"""
Molecular Visualization and Analysis Tool
Main entry point for comprehensive molecular analysis including:
- Property visualization with user-configurable axes
- Chemical space analysis (PCA + t-SNE) 
- Optional molecular docking with AutoDock Vina
- Interactive HTML dashboard with synchronized highlighting
- Professional styling and color bars
"""

import sys
import logging
import argparse
from pathlib import Path
import traceback

# Add src directory to Python path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

try:
    from config_manager import ConfigManager, ConfigurationError
    from molecular_data_processor import MolecularDataProcessor
    from docking_wrapper import VinaDockingWrapper, DockingError
    from dashboard_generator import DashboardGenerator
except ImportError as e:
    print(f"Error importing modules: {e}")
    print("Make sure all required dependencies are installed:")
    print("  conda install -c conda-forge rdkit pandas numpy scikit-learn")
    sys.exit(1)


def setup_logging(log_level: str = 'INFO', log_file: str = None) -> None:
    """
    Setup logging configuration
    
    Args:
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR)
        log_file: Optional log file path
    """
    # Configure logging format
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    date_format = '%Y-%m-%d %H:%M:%S'
    
    # Setup logging handlers
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        handlers.append(logging.FileHandler(log_file))
    
    # Configure logging
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format=log_format,
        datefmt=date_format,
        handlers=handlers
    )
    
    # Reduce noise from some libraries
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.getLogger('PIL').setLevel(logging.WARNING)


def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments
    
    Returns:
        Parsed arguments namespace
    """
    parser = argparse.ArgumentParser(
        description="Molecular Visualization and Analysis Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic analysis with CSV file
  python run_analysis.py input.csv
  
  # Use custom configuration file
  python run_analysis.py input.csv --config config.yaml
  
  # Enable docking with protein structure
  python run_analysis.py input.csv --protein protein.pdb --vina-config vina.conf
  
  # Specify custom properties and output
  python run_analysis.py input.csv --properties MW LogP TPSA --output results.html
  
  # Generate sample configuration file
  python run_analysis.py --generate-config
        """
    )
    
    # Input options
    parser.add_argument('input', nargs='?', help='Input CSV file with SMILES and properties')
    parser.add_argument('--config', '-c', help='YAML configuration file')
    parser.add_argument('--output', '-o', help='Output HTML file')
    
    # Property options
    parser.add_argument('--properties', nargs='+', 
                       help='Property columns to include in analysis')
    parser.add_argument('--smiles-column', default='SMILES',
                       help='Name of SMILES column (default: SMILES)')
    parser.add_argument('--name-column', default='Name',
                       help='Name of compound name column (default: Name)')
    
    # Docking options
    parser.add_argument('--protein-pdb', help='Protein PDB file for docking')
    parser.add_argument('--vina-config', help='AutoDock Vina configuration file')
    parser.add_argument('--no-docking', action='store_true',
                       help='Disable docking even if configured')
    
    # Analysis options
    parser.add_argument('--no-descriptors', action='store_true',
                       help='Skip molecular descriptor calculation')
    parser.add_argument('--no-pca', action='store_true',
                       help='Skip PCA analysis')
    parser.add_argument('--no-tsne', action='store_true',
                       help='Skip t-SNE analysis')
    
    # Output options
    parser.add_argument('--export-data', action='store_true',
                       help='Export processed data to CSV')
    parser.add_argument('--export-file', default='analysis_results.csv',
                       help='Export CSV filename')
    
    # Utility options
    parser.add_argument('--generate-config', action='store_true',
                       help='Generate sample configuration file and exit')
    parser.add_argument('--validate-config', action='store_true',
                       help='Validate configuration file and exit')
    
    # Debug options
    parser.add_argument('--debug', action='store_true',
                       help='Enable debug logging')
    parser.add_argument('--log-file', help='Write log to file')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose output')
    
    return parser.parse_args()


def generate_sample_config() -> None:
    """Generate sample configuration file"""
    print("Generating sample configuration file...")
    
    config_template_path = Path(__file__).parent / 'config_template.yaml'
    output_path = Path('molecular_analysis_config.yaml')
    
    if config_template_path.exists():
        import shutil
        shutil.copy(config_template_path, output_path)
        print(f"✅ Sample configuration saved to: {output_path}")
        print("📝 Edit this file to customize your analysis settings")
    else:
        print("❌ Configuration template not found")
        sys.exit(1)


def validate_configuration(config_path: str) -> None:
    """
    Validate configuration file
    
    Args:
        config_path: Path to configuration file
    """
    print(f"Validating configuration file: {config_path}")
    
    try:
        config = ConfigManager(config_path)
        print("✅ Configuration file is valid")
        print(f"📊 Input file: {config.get_input_file()}")
        print(f"📈 Output file: {config.get_output_file()}")
        print(f"🧪 Docking enabled: {config.is_docking_enabled()}")
        print(f"🔬 Properties: {len(config.get_property_columns())} specified")
        
    except ConfigurationError as e:
        print(f"❌ Configuration validation failed: {e}")
        sys.exit(1)


def main():
    """Main analysis pipeline"""
    args = parse_arguments()
    
    # Handle utility commands
    if args.generate_config:
        generate_sample_config()
        return
    
    if args.validate_config:
        if not args.config:
            print("❌ Please specify configuration file with --config")
            sys.exit(1)
        validate_configuration(args.config)
        return
    
    # Validate required arguments
    if not args.input:
        print("❌ Please specify input CSV file")
        print("Use 'python run_analysis.py --help' for usage information")
        sys.exit(1)
    
    # Setup logging
    log_level = 'DEBUG' if args.debug else 'INFO'
    setup_logging(log_level, args.log_file)
    logger = logging.getLogger(__name__)
    
    logger.info("=== Molecular Visualization and Analysis Tool ===")
    logger.info(f"Input file: {args.input}")
    
    try:
        # Initialize configuration
        logger.info("Loading configuration...")
        if args.config:
            config = ConfigManager(args.config)
            logger.info(f"Configuration loaded from: {args.config}")
        else:
            config = ConfigManager()
            logger.info("Using default configuration")
        
        # Update configuration from command line arguments
        config.update_from_args(args)
        
        # Override specific settings from arguments
        if args.no_descriptors:
            config.set('analysis.calculate_descriptors', False)
        if args.no_pca:
            config.set('analysis.chemical_space.pca.enabled', False)
        if args.no_tsne:
            config.set('analysis.chemical_space.tsne.enabled', False)
        
        # Validate input file
        input_file = Path(args.input)
        if not input_file.exists():
            raise FileNotFoundError(f"Input file not found: {args.input}")
        
        # Initialize data processor
        logger.info("Initializing molecular data processor...")
        processor = MolecularDataProcessor(config)
        
        # Load and process molecular data
        logger.info("Loading and processing molecular data...")
        df = processor.load_and_process_data(str(input_file))
        
        # Calculate molecular descriptors
        if config.get('analysis.calculate_descriptors', True):
            logger.info("Calculating molecular descriptors...")
            df = processor.calculate_molecular_descriptors()
        
        # Perform chemical space analysis
        logger.info("Performing chemical space analysis...")
        pca_coords, tsne_coords, pca_variance = processor.perform_chemical_space_analysis()
        
        # Initialize docking if enabled
        docking_data = None
        if config.is_docking_enabled():
            logger.info("Initializing molecular docking...")
            try:
                docking_wrapper = VinaDockingWrapper(config)
                
                # Prepare ligands
                ligand_files = docking_wrapper.prepare_ligands(df)
                
                if ligand_files:
                    # Run docking
                    docking_results = docking_wrapper.run_vina_docking(ligand_files)
                    
                    # Integrate results
                    df = docking_wrapper.integrate_docking_results(df)
                    
                    # Prepare visualization data
                    docking_data = docking_wrapper.prepare_docking_visualization_data(df)
                    
                    logger.info(f"Docking completed for {len(docking_data)} compounds")
                
            except DockingError as e:
                logger.warning(f"Docking failed: {e}")
                logger.warning("Continuing without docking analysis")
                config.set('docking.enabled', False)
            
            except Exception as e:
                logger.error(f"Unexpected docking error: {e}")
                logger.warning("Continuing without docking analysis")
                config.set('docking.enabled', False)
        
        # Prepare visualization data
        logger.info("Preparing visualization data...")
        data_points = processor.prepare_visualization_data()
        data_summary = processor.get_data_summary()
        
        # Generate dashboard
        logger.info("Generating interactive dashboard...")
        dashboard = DashboardGenerator(config)
        html_content = dashboard.generate_dashboard(data_points, data_summary, docking_data)
        
        # Save dashboard
        output_file = config.get_output_file()
        dashboard.save_dashboard(html_content, output_file)
        
        # Export data if requested
        if args.export_data or config.get('export.export_data', False):
            export_file = args.export_file or config.get('export.export_file', 'analysis_results.csv')
            logger.info(f"Exporting processed data to: {export_file}")
            
            # Remove Mol column before export (not serializable)
            export_df = df.drop(columns=['Mol'], errors='ignore')
            export_df.to_csv(export_file, index=False)
        
        # Print summary
        logger.info("=== Analysis Complete ===")
        logger.info(f"📊 Processed {len(data_points)} compounds")
        logger.info(f"🧪 {len(data_summary.get('available_properties', []))} properties analyzed")
        logger.info(f"📈 Dashboard saved to: {output_file}")
        
        if data_summary.get('has_pca', False):
            pca_variance = data_summary.get('pca_variance', [])
            if pca_variance:
                total_variance = sum(pca_variance) * 100
                logger.info(f"🔬 PCA explained variance: {total_variance:.1f}%")
        
        if config.is_docking_enabled() and docking_data:
            logger.info(f"⚗️  Docking: {len(docking_data)} successful poses")
        
        print(f"\n🎉 Analysis completed successfully!")
        print(f"📊 Open {output_file} in your web browser to explore the results")
        print(f"🔍 Features available:")
        print(f"   • Interactive property plots with configurable axes")
        print(f"   • Chemical space visualization (PCA + t-SNE)")
        print(f"   • Synchronized highlighting across all plots")
        print(f"   • Molecular structure display")
        print(f"   • Professional styling with color bars")
        
        if config.is_docking_enabled() and docking_data:
            print(f"   • Molecular docking results and pose visualization")
        
    except ConfigurationError as e:
        logger.error(f"Configuration error: {e}")
        sys.exit(1)
    
    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        sys.exit(1)
    
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        if args.debug or args.verbose:
            logger.error("Full traceback:")
            logger.error(traceback.format_exc())
        sys.exit(1)


if __name__ == "__main__":
    main()