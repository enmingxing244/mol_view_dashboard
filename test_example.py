#!/usr/bin/env python3
"""
Test Example for Molecular Visualization and Analysis Tool
Creates sample data and validates the enhanced analysis pipeline
"""

import pandas as pd
import numpy as np
import sys
from pathlib import Path

# Add src directory to Python path
sys.path.insert(0, str(Path(__file__).parent / 'src'))

def create_sample_data():
    """Create sample molecular data for testing"""
    print("Creating sample data...")
    
    # Sample SMILES strings for common molecules
    sample_molecules = [
        ("CCO", "Ethanol"),
        ("CC(C)O", "Isopropanol"),
        ("c1ccccc1", "Benzene"),
        ("CC(=O)O", "Acetic acid"),
        ("CCCCCCCCO", "Octanol"),
        ("Cc1ccccc1", "Toluene"),
        ("CCN(CC)CC", "Triethylamine"),
        ("CC(C)(C)O", "tert-Butanol"),
        ("c1ccc2ccccc2c1", "Naphthalene"),
        ("CC(=O)Oc1ccccc1C(=O)O", "Aspirin"),
        ("CC1=CC=C(C=C1)C(C)C", "p-Cymene"),
        ("CCCCCCCCCCCCCCCC(=O)O", "Palmitic acid"),
        ("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "Ibuprofen"),
        ("CC1=C(C(=O)N(N1C)C)C", "Antipyrine"),
        ("CC(C)NCC(C1=CC(=C(C=C1)O)CO)O", "Salbutamol"),
    ]
    
    # Create DataFrame
    data = []
    for i, (smiles, name) in enumerate(sample_molecules):
        # Add some synthetic properties for demonstration
        mw = np.random.normal(200, 100)  # Molecular weight
        logp = np.random.normal(2, 1.5)  # LogP
        tpsa = np.random.normal(50, 30)  # TPSA
        ic50 = np.random.lognormal(0, 1)  # IC50 (nM)
        binding_score = np.random.normal(-7, 2)  # Binding score
        
        data.append({
            'SMILES': smiles,
            'Name': name,
            'MW': max(50, mw),  # Ensure positive MW
            'LogP': logp,
            'TPSA': max(0, tpsa),  # Ensure positive TPSA
            'IC50_nM': ic50,
            'Binding_Score': binding_score
        })
    
    df = pd.DataFrame(data)
    
    # Save to CSV
    output_file = 'sample_compounds.csv'
    df.to_csv(output_file, index=False)
    print(f"Sample data saved to: {output_file}")
    print(f"Data shape: {df.shape}")
    print(f"Columns: {list(df.columns)}")
    
    return output_file

def test_configuration():
    """Test configuration management"""
    print("\nTesting configuration management...")
    
    try:
        from config_manager import ConfigManager
        
        # Test default configuration
        config = ConfigManager()
        print("✅ Default configuration loaded successfully")
        
        # Test configuration access
        input_file = config.get('input.csv_file')
        output_file = config.get_output_file()
        docking_enabled = config.is_docking_enabled()
        
        print(f"   Input file: {input_file}")
        print(f"   Output file: {output_file}")
        print(f"   Docking enabled: {docking_enabled}")
        
        return True
        
    except Exception as e:
        print(f"❌ Configuration test failed: {e}")
        return False

def test_data_processor():
    """Test molecular data processor"""
    print("\nTesting molecular data processor...")
    
    try:
        from config_manager import ConfigManager
        from molecular_data_processor import MolecularDataProcessor
        
        # Create sample data
        csv_file = create_sample_data()
        
        # Initialize processor
        config = ConfigManager()
        config.set('input.csv_file', csv_file)
        config.set('input.property_columns', ['MW', 'LogP', 'TPSA', 'IC50_nM', 'Binding_Score'])
        
        processor = MolecularDataProcessor(config)
        
        # Load and process data
        df = processor.load_and_process_data(csv_file)
        print(f"✅ Data loaded successfully: {len(df)} compounds")
        
        # Calculate descriptors
        df = processor.calculate_molecular_descriptors()
        print(f"✅ Descriptors calculated successfully")
        
        # Perform chemical space analysis
        pca_coords, tsne_coords, pca_variance = processor.perform_chemical_space_analysis()
        print(f"✅ Chemical space analysis completed")
        print(f"   PCA explained variance: {pca_variance}")
        
        # Prepare visualization data
        data_points = processor.prepare_visualization_data()
        print(f"✅ Visualization data prepared: {len(data_points)} compounds")
        
        # Get data summary
        summary = processor.get_data_summary()
        print(f"✅ Data summary generated")
        print(f"   Available properties: {len(summary.get('available_properties', []))}")
        
        return True
        
    except Exception as e:
        print(f"❌ Data processor test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_dashboard_generator():
    """Test dashboard generation"""
    print("\nTesting dashboard generator...")
    
    try:
        from config_manager import ConfigManager
        from molecular_data_processor import MolecularDataProcessor
        from dashboard_generator import DashboardGenerator
        
        # Create sample data and process it
        csv_file = create_sample_data()
        
        config = ConfigManager()
        config.set('input.csv_file', csv_file)
        config.set('input.property_columns', ['MW', 'LogP', 'TPSA', 'IC50_nM', 'Binding_Score'])
        config.set('visualization.output_file', 'test_dashboard.html')
        
        processor = MolecularDataProcessor(config)
        df = processor.load_and_process_data(csv_file)
        df = processor.calculate_molecular_descriptors()
        processor.perform_chemical_space_analysis()
        
        data_points = processor.prepare_visualization_data()
        data_summary = processor.get_data_summary()
        
        # Generate dashboard
        dashboard = DashboardGenerator(config)
        html_content = dashboard.generate_dashboard(data_points, data_summary)
        
        # Save dashboard
        output_file = config.get_output_file()
        dashboard.save_dashboard(html_content, output_file)
        
        print(f"✅ Dashboard generated successfully: {output_file}")
        print(f"   HTML content length: {len(html_content)} characters")
        
        # Verify file exists
        if Path(output_file).exists():
            print(f"✅ Dashboard file created successfully")
            return True
        else:
            print(f"❌ Dashboard file not found")
            return False
        
    except Exception as e:
        print(f"❌ Dashboard generator test failed: {e}")
        import traceback
        traceback.print_exc()
        return False

def test_main_interface():
    """Test main run_analysis.py interface"""
    print("\nTesting main interface...")
    
    try:
        import subprocess
        
        # Create sample data
        csv_file = create_sample_data()
        
        # Test basic usage
        cmd = [
            'python', 'run_analysis.py', csv_file,
            '--output', 'test_main_output.html',
            '--properties', 'MW', 'LogP', 'TPSA',
            '--no-docking'  # Disable docking for testing
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        
        if result.returncode == 0:
            print("✅ Main interface test passed")
            print("   Output:")
            for line in result.stdout.split('\n')[-10:]:  # Show last 10 lines
                if line.strip():
                    print(f"   {line}")
            
            # Check if output file was created
            if Path('test_main_output.html').exists():
                print("✅ Output file created successfully")
                return True
            else:
                print("❌ Output file not found")
                return False
        else:
            print("❌ Main interface test failed")
            print("   Error output:")
            print(result.stderr)
            return False
            
    except subprocess.TimeoutExpired:
        print("❌ Main interface test timed out")
        return False
    except Exception as e:
        print(f"❌ Main interface test failed: {e}")
        return False

def cleanup_test_files():
    """Clean up test files"""
    print("\nCleaning up test files...")
    
    test_files = [
        'sample_compounds.csv',
        'test_dashboard.html',
        'test_main_output.html',
        'analysis_results.csv',
        'molecular_analysis_config.yaml'
    ]
    
    for file in test_files:
        if Path(file).exists():
            Path(file).unlink()
            print(f"   Removed: {file}")

def main():
    """Run all tests"""
    print("=== Testing Molecular Visualization and Analysis Tool ===")
    
    tests = [
        ("Configuration Management", test_configuration),
        ("Data Processor", test_data_processor),
        ("Dashboard Generator", test_dashboard_generator),
        ("Main Interface", test_main_interface),
    ]
    
    results = []
    
    for test_name, test_func in tests:
        print(f"\n{'='*60}")
        print(f"Running: {test_name}")
        print('='*60)
        
        try:
            success = test_func()
            results.append((test_name, success))
        except Exception as e:
            print(f"❌ {test_name} failed with exception: {e}")
            results.append((test_name, False))
    
    # Print summary
    print(f"\n{'='*60}")
    print("TEST SUMMARY")
    print('='*60)
    
    passed = 0
    for test_name, success in results:
        status = "✅ PASSED" if success else "❌ FAILED"
        print(f"{test_name:<30} {status}")
        if success:
            passed += 1
    
    print(f"\nOverall: {passed}/{len(results)} tests passed")
    
    if passed == len(results):
        print("🎉 All tests passed! The system is working correctly.")
    else:
        print("⚠️  Some tests failed. Please check the error messages above.")
    
    # Cleanup
    cleanup_test_files()
    
    return passed == len(results)

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)