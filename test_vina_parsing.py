#!/usr/bin/env python3
"""
Test Vina output parsing
"""

def parse_vina_output(log_file: str) -> float:
    """Parse Vina log file to extract binding energy"""
    try:
        with open(log_file, 'r') as f:
            lines = f.readlines()
        
        print(f"Total lines: {len(lines)}")
        
        # Look for the results table
        for i, line in enumerate(lines):
            if 'mode' in line and 'affinity' in line:
                print(f"Found header at line {i}: {line.strip()}")
                # Skip the units line and get the first data line
                if i + 3 < len(lines):
                    print(f"Line {i+1}: {lines[i+1].strip()}")
                    print(f"Line {i+2}: {lines[i+2].strip()}")
                    result_line = lines[i + 3].strip()
                    print(f"Result line: '{result_line}'")
                    parts = result_line.split()
                    print(f"Parts: {parts}")
                    if len(parts) >= 2:
                        try:
                            energy = float(parts[1])
                            print(f"Extracted energy: {energy}")
                            return energy
                        except ValueError as e:
                            print(f"ValueError: {e}")
                            pass
        
        print("No results table found")
        return float('nan')
        
    except Exception as e:
        print(f"Error: {e}")
        return float('nan')

# Test with actual log file
log_file = "docking_results/log_0000_Ethanol.txt"
print(f"Testing with: {log_file}")
result = parse_vina_output(log_file)
print(f"Final result: {result}")