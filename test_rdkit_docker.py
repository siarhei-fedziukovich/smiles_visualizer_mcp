#!/usr/bin/env python3
"""
Test script to verify RDKit installation in Docker container
"""

import sys
import os

def test_rdkit_import():
    """Test RDKit import and basic functionality"""
    try:
        print("Testing RDKit import...")
        from rdkit import Chem
        print("✓ RDKit imported successfully")
        
        # Test basic functionality
        mol = Chem.MolFromSmiles("CCO")
        if mol is not None:
            print("✓ RDKit can parse SMILES strings")
            smiles = Chem.MolToSmiles(mol)
            print(f"  Example: CCO -> {smiles}")
        else:
            print("✗ RDKit cannot parse SMILES strings")
            return False
            
        # Test drawing functionality
        try:
            from rdkit.Chem import Draw
            print("✓ RDKit drawing module imported")
        except ImportError as e:
            print(f"✗ RDKit drawing module import failed: {e}")
            return False
            
        return True
        
    except ImportError as e:
        print(f"✗ RDKit import failed: {e}")
        return False

def test_matplotlib():
    """Test matplotlib in headless environment"""
    try:
        print("\nTesting matplotlib...")
        import matplotlib
        print(f"✓ Matplotlib version: {matplotlib.__version__}")
        print(f"✓ Backend: {matplotlib.get_backend()}")
        
        import matplotlib.pyplot as plt
        print("✓ Matplotlib.pyplot imported successfully")
        
        # Test basic plotting
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3], [1, 4, 2])
        plt.close()
        print("✓ Matplotlib can create plots")
        
        return True
        
    except Exception as e:
        print(f"✗ Matplotlib test failed: {e}")
        return False

def test_other_dependencies():
    """Test other dependencies"""
    dependencies = [
        ("numpy", "NumPy"),
        ("pandas", "Pandas"),
        ("networkx", "NetworkX"),
        ("plotly", "Plotly"),
        ("PIL", "Pillow"),
        ("mcp", "MCP SDK")
    ]
    
    print("\nTesting other dependencies...")
    all_passed = True
    
    for module, name in dependencies:
        try:
            __import__(module)
            print(f"✓ {name} imported successfully")
        except ImportError as e:
            print(f"✗ {name} import failed: {e}")
            all_passed = False
    
    return all_passed

def main():
    """Main test function"""
    print("=" * 50)
    print("RDKit Docker Container Test")
    print("=" * 50)
    
    # Set matplotlib backend
    os.environ['MPLBACKEND'] = 'Agg'
    
    # Run tests
    rdkit_ok = test_rdkit_import()
    matplotlib_ok = test_matplotlib()
    other_ok = test_other_dependencies()
    
    print("\n" + "=" * 50)
    print("Test Results:")
    print(f"RDKit: {'✓ PASS' if rdkit_ok else '✗ FAIL'}")
    print(f"Matplotlib: {'✓ PASS' if matplotlib_ok else '✗ FAIL'}")
    print(f"Other Dependencies: {'✓ PASS' if other_ok else '✗ FAIL'}")
    
    if rdkit_ok and matplotlib_ok and other_ok:
        print("\n🎉 All tests passed! Container is ready.")
        sys.exit(0)
    else:
        print("\n❌ Some tests failed. Check the output above.")
        sys.exit(1)

if __name__ == "__main__":
    main()
