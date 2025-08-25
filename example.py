#!/usr/bin/env python3
"""
Example usage of SMILES Visualizer MCP Server
"""

import asyncio
import json
import base64
from pathlib import Path

# Example SMILES strings
EXAMPLE_SMILES = {
    "ethanol": "CCO",
    "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "biphenyl": "C1=CC=C(C=C1)C2=CC=CC=C2"
}

async def main():
    """Main example function"""
    print("SMILES Visualizer MCP Server Example")
    print("=" * 50)
    print()
    
    # This is a demonstration of how to use the MCP server
    # In a real scenario, you would connect to the running server
    
    print("Example SMILES strings:")
    for name, smiles in EXAMPLE_SMILES.items():
        print(f"  {name}: {smiles}")
    print()
    
    print("Available tools:")
    tools = [
        "validate_smiles",
        "get_molecular_info", 
        "visualize_rdkit",
        "visualize_network",
        "visualize_plotly",
        "visualize_custom_matplotlib",
        "compare_visualizations",
        "batch_visualize"
    ]
    
    for tool in tools:
        print(f"  - {tool}")
    print()
    
    print("Example usage patterns:")
    print()
    
    # Example 1: Validate SMILES
    print("1. Validate SMILES:")
    print("   await validate_smiles('CCO')")
    print("   # Returns: {'valid': True, 'message': 'Valid SMILES', 'canonical_smiles': 'CCO'}")
    print()
    
    # Example 2: Get molecular information
    print("2. Get molecular information:")
    print("   await get_molecular_info('CCO')")
    print("   # Returns molecular weight, atom count, properties, etc.")
    print()
    
    # Example 3: Create RDKit visualization
    print("3. Create RDKit visualization:")
    print("   await visualize_rdkit('CCO', size='400,300')")
    print("   # Returns base64 encoded PNG image")
    print()
    
    # Example 4: Create network visualization
    print("4. Create network visualization:")
    print("   await visualize_network('CCO', layout='spring')")
    print("   # Returns base64 encoded PNG image with network layout")
    print()
    
    # Example 5: Create interactive Plotly visualization
    print("5. Create interactive Plotly visualization:")
    print("   await visualize_plotly('CCO')")
    print("   # Returns HTML content for interactive chart")
    print()
    
    # Example 6: Compare all visualizations
    print("6. Compare all visualizations:")
    print("   await compare_visualizations('CCO')")
    print("   # Returns all visualization types and molecular info")
    print()
    
    # Example 7: Batch processing
    print("7. Batch processing:")
    print("   smiles_list = ['CCO', 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C']")
    print("   await batch_visualize(smiles_list, visualization_type='rdkit')")
    print("   # Returns visualizations for all molecules")
    print()
    
    print("Integration with MCP clients:")
    print()
    print("The server can be integrated with any MCP-compatible client:")
    print("- Claude Desktop")
    print("- VS Code with MCP extension")
    print("- Custom MCP clients")
    print()
    
    print("Client configuration example:")
    print("""
{
  "mcpServers": {
    "smiles-visualizer": {
      "command": "python",
      "args": ["path/to/smiles_visualizer_mcp/server.py", "--host", "127.0.0.1", "--port", "8001"],
      "env": {
        "MCP_HOST": "127.0.0.1",
        "MCP_PORT": "8001"
      }
    }
  }
}
""")
    
    print("Output formats:")
    print("- Images: PNG format, base64 encoded")
    print("- Interactive visualizations: HTML (Plotly)")
    print("- Data: JSON format")
    print()
    
    print("Error handling:")
    print("- Invalid SMILES strings return error messages")
    print("- Missing dependencies are handled gracefully")
    print("- All errors are returned as structured JSON")
    print()
    
    print("To run the server:")
    print("  python server.py --host 127.0.0.1 --port 8001")
    print()
    
    print("To test the server:")
    print("  python test_server.py")
    print()

if __name__ == "__main__":
    asyncio.run(main())
