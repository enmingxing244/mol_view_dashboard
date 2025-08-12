#!/usr/bin/env python3
"""
Quick script to display analysis results
"""

import pandas as pd
import numpy as np

def show_results():
    """Display analysis results summary"""
    
    # Read the processed data
    df = pd.read_csv('processed_data.csv')
    
    print("=== MOLECULAR ANALYSIS RESULTS ===")
    print(f"📊 Total compounds analyzed: {len(df)}")
    print()
    
    print("📋 INPUT DATA OVERVIEW:")
    print(f"   Original CSV columns: SMILES, Name, Activity_uM, Binding_Score")
    print(f"   Sample compounds: {', '.join(df['Name'].head(3).tolist())}...")
    print()
    
    print("🧪 CALCULATED RDKit DESCRIPTORS:")
    descriptor_cols = ['MW', 'LogP', 'HBA', 'HBD', 'NumRings', 'QED', 'SAscore']
    
    for col in descriptor_cols:
        if col in df.columns:
            values = df[col].dropna()
            if len(values) > 0:
                print(f"   {col:<12}: {values.min():.2f} - {values.max():.2f} (mean: {values.mean():.2f})")
            else:
                print(f"   {col:<12}: No valid data")
    print()
    
    print("🔬 CHEMICAL SPACE ANALYSIS:")
    if 'PCA_1' in df.columns and 'PCA_2' in df.columns:
        print("   ✅ PCA analysis completed (2 components)")
        print(f"   PCA range X: {df['PCA_1'].min():.2f} to {df['PCA_1'].max():.2f}")
        print(f"   PCA range Y: {df['PCA_2'].min():.2f} to {df['PCA_2'].max():.2f}")
    else:
        print("   ❌ PCA analysis not available")
        
    if 'tSNE_1' in df.columns and 'tSNE_2' in df.columns:
        print("   ✅ t-SNE analysis completed (2 components)")
        print(f"   t-SNE range X: {df['tSNE_1'].min():.1f} to {df['tSNE_1'].max():.1f}")
        print(f"   t-SNE range Y: {df['tSNE_2'].min():.1f} to {df['tSNE_2'].max():.1f}")
    else:
        print("   ❌ t-SNE analysis not available")
    print()
    
    print("📈 TOP COMPOUNDS BY DIFFERENT METRICS:")
    
    # Top compounds by activity
    if 'Activity_uM' in df.columns:
        top_active = df.nsmallest(3, 'Activity_uM')[['Name', 'Activity_uM', 'MW', 'LogP', 'QED']]
        print("   Most Active (lowest Activity_uM):")
        for _, row in top_active.iterrows():
            print(f"      {row['Name']:<15} Activity: {row['Activity_uM']:.1f} μM, MW: {row['MW']:.1f}, LogP: {row['LogP']:.2f}, QED: {row['QED']:.3f}")
    
    # Top compounds by binding score
    if 'Binding_Score' in df.columns:
        top_binding = df.nsmallest(3, 'Binding_Score')[['Name', 'Binding_Score', 'MW', 'LogP', 'QED']]
        print("   Best Binding Score (most negative):")
        for _, row in top_binding.iterrows():
            print(f"      {row['Name']:<15} Score: {row['Binding_Score']:.1f}, MW: {row['MW']:.1f}, LogP: {row['LogP']:.2f}, QED: {row['QED']:.3f}")
    
    # Top drug-like compounds
    if 'QED' in df.columns:
        top_druglike = df.nlargest(3, 'QED')[['Name', 'QED', 'MW', 'LogP', 'SAscore']]
        print("   Most Drug-like (highest QED):")
        for _, row in top_druglike.iterrows():
            print(f"      {row['Name']:<15} QED: {row['QED']:.3f}, MW: {row['MW']:.1f}, LogP: {row['LogP']:.2f}, SA: {row['SAscore']:.2f}")
    print()
    
    print("📊 AVAILABLE VISUALIZATIONS:")
    print("   🎯 Interactive property plots with configurable X/Y axes")
    print("   🗺️  PCA chemical space visualization")  
    print("   🌐 t-SNE chemical space visualization")
    print("   🔬 Molecular structure display on hover")
    print("   🎨 Color-coded properties with color bars")
    print("   🔗 Synchronized highlighting across all plots")
    print()
    
    print("🌐 TO VIEW RESULTS:")
    print("   1. Open 'results.html' in your web browser")
    print("   2. Use the dropdown menus to select X-axis, Y-axis, and color properties")
    print("   3. Hover over data points to see molecular structures")
    print("   4. Click points to pin structure display")
    print("   5. Notice how highlighting works across all three plots simultaneously")
    print()
    
    print("📁 OUTPUT FILES CREATED:")
    print(f"   📄 results.html - Interactive dashboard ({220278/1024:.0f} KB)")
    print(f"   📊 processed_data.csv - All data with calculated descriptors")
    print(f"   📝 sample_compounds.csv - Original input data")

if __name__ == "__main__":
    show_results()