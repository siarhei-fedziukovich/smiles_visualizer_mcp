#!/usr/bin/env python3
"""
Test script for SMILES Visualizer MCP Server
"""

import asyncio
import json
import requests
import time
from typing import Dict, Any

class SMILESVisualizerTester:
    """Test client for SMILES Visualizer MCP Server"""
    
    def __init__(self, base_url: str = "http://127.0.0.1:8080"):
        self.base_url = base_url
        self.session = requests.Session()
    
    def test_tool(self, tool_name: str, arguments: Dict[str, Any]) -> Dict[str, Any]:
        """Test a specific tool"""
        url = f"{self.base_url}/tools/{tool_name}/call"
        
        try:
            response = self.session.post(
                url,
                json={"arguments": arguments},
                headers={"Content-Type": "application/json"},
                timeout=30
            )
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            return {"error": f"Request failed: {str(e)}"}
    
    def test_validate_smiles(self):
        """Test SMILES validation"""
        print("Testing SMILES validation...")
        
        test_cases = [
            ("CCO", "Valid ethanol"),
            ("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "Valid ibuprofen"),
            ("INVALID_SMILES", "Invalid SMILES"),
            ("", "Empty SMILES")
        ]
        
        for smiles, description in test_cases:
            print(f"  Testing {description}: {smiles}")
            result = self.test_tool("validate_smiles", {"smiles": smiles})
            print(f"    Result: {json.dumps(result, indent=2)}")
            print()
    
    def test_get_molecular_info(self):
        """Test molecular information retrieval"""
        print("Testing molecular information retrieval...")
        
        test_smiles = "CCO"  # Ethanol
        print(f"  Testing molecular info for: {test_smiles}")
        result = self.test_tool("get_molecular_info", {"smiles": test_smiles})
        print(f"    Result: {json.dumps(result, indent=2)}")
        print()
    
    def test_visualize_rdkit(self):
        """Test RDKit visualization"""
        print("Testing RDKit visualization...")
        
        test_smiles = "CCO"  # Ethanol
        print(f"  Testing RDKit visualization for: {test_smiles}")
        result = self.test_tool("visualize_rdkit", {
            "smiles": test_smiles,
            "size": "400,300"
        })
        
        if "content" in result and len(result["content"]) > 0:
            content = result["content"][0]
            if "text" in content:
                try:
                    data = json.loads(content["text"])
                    if "success" in data and data["success"]:
                        print(f"    Success: Generated {data['format']} image ({data['size']})")
                        print(f"    Base64 length: {len(data['image_base64'])} characters")
                    else:
                        print(f"    Error: {data.get('error', 'Unknown error')}")
                except json.JSONDecodeError:
                    print(f"    Error: Invalid JSON response")
            else:
                print(f"    Error: No text content in response")
        else:
            print(f"    Error: No content in response")
        print()
    
    def test_visualize_network(self):
        """Test network visualization"""
        print("Testing network visualization...")
        
        test_smiles = "CCO"  # Ethanol
        print(f"  Testing network visualization for: {test_smiles}")
        result = self.test_tool("visualize_network", {
            "smiles": test_smiles,
            "layout": "spring"
        })
        
        if "content" in result and len(result["content"]) > 0:
            content = result["content"][0]
            if "text" in content:
                try:
                    data = json.loads(content["text"])
                    if "success" in data and data["success"]:
                        print(f"    Success: Generated {data['format']} image")
                        print(f"    Layout: {data['layout']}")
                        print(f"    Nodes: {data['num_nodes']}, Edges: {data['num_edges']}")
                        print(f"    Base64 length: {len(data['image_base64'])} characters")
                    else:
                        print(f"    Error: {data.get('error', 'Unknown error')}")
                except json.JSONDecodeError:
                    print(f"    Error: Invalid JSON response")
            else:
                print(f"    Error: No text content in response")
        else:
            print(f"    Error: No content in response")
        print()
    
    def test_visualize_plotly(self):
        """Test Plotly visualization"""
        print("Testing Plotly visualization...")
        
        test_smiles = "CCO"  # Ethanol
        print(f"  Testing Plotly visualization for: {test_smiles}")
        result = self.test_tool("visualize_plotly", {"smiles": test_smiles})
        
        if "content" in result and len(result["content"]) > 0:
            content = result["content"][0]
            if "text" in content:
                try:
                    data = json.loads(content["text"])
                    if "success" in data and data["success"]:
                        print(f"    Success: Generated {data['format']} visualization")
                        print(f"    HTML length: {len(data['html_content'])} characters")
                    else:
                        print(f"    Error: {data.get('error', 'Unknown error')}")
                except json.JSONDecodeError:
                    print(f"    Error: Invalid JSON response")
            else:
                print(f"    Error: No text content in response")
        else:
            print(f"    Error: No content in response")
        print()
    
    def test_compare_visualizations(self):
        """Test comparison of all visualization types"""
        print("Testing comparison of all visualization types...")
        
        test_smiles = "CCO"  # Ethanol
        print(f"  Testing all visualizations for: {test_smiles}")
        result = self.test_tool("compare_visualizations", {"smiles": test_smiles})
        
        if "content" in result and len(result["content"]) > 0:
            content = result["content"][0]
            if "text" in content:
                try:
                    data = json.loads(content["text"])
                    if "error" not in data:
                        print(f"    Success: Generated all visualization types")
                        print(f"    Available types: {list(data.keys())}")
                        
                        for viz_type, viz_data in data.items():
                            if viz_type == "molecular_info":
                                print(f"      {viz_type}: Molecular information available")
                            elif "success" in viz_data and viz_data["success"]:
                                print(f"      {viz_type}: Success")
                            else:
                                print(f"      {viz_type}: Error - {viz_data.get('error', 'Unknown error')}")
                    else:
                        print(f"    Error: {data.get('error', 'Unknown error')}")
                except json.JSONDecodeError:
                    print(f"    Error: Invalid JSON response")
            else:
                print(f"    Error: No text content in response")
        else:
            print(f"    Error: No content in response")
        print()
    
    def test_batch_visualize(self):
        """Test batch visualization"""
        print("Testing batch visualization...")
        
        test_smiles_list = ["CCO", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"]
        print(f"  Testing batch visualization for {len(test_smiles_list)} molecules")
        result = self.test_tool("batch_visualize", {
            "smiles_list": test_smiles_list,
            "visualization_type": "rdkit"
        })
        
        if "content" in result and len(result["content"]) > 0:
            content = result["content"][0]
            if "text" in content:
                try:
                    data = json.loads(content["text"])
                    if "error" not in data:
                        print(f"    Success: Processed {len(data)} molecules")
                        for mol_name, mol_data in data.items():
                            if "error" in mol_data:
                                print(f"      {mol_name}: Error - {mol_data['error']}")
                            else:
                                print(f"      {mol_name}: Success")
                    else:
                        print(f"    Error: {data.get('error', 'Unknown error')}")
                except json.JSONDecodeError:
                    print(f"    Error: Invalid JSON response")
            else:
                print(f"    Error: No text content in response")
        else:
            print(f"    Error: No content in response")
        print()
    
    def run_all_tests(self):
        """Run all tests"""
        print("=" * 60)
        print("SMILES Visualizer MCP Server Test Suite")
        print("=" * 60)
        print()
        
        # Check if server is running
        try:
            response = self.session.get(f"{self.base_url}/tools", timeout=5)
            if response.status_code == 200:
                print("✓ Server is running and responding")
                print()
            else:
                print("✗ Server is not responding correctly")
                return
        except requests.exceptions.RequestException:
            print("✗ Cannot connect to server. Make sure it's running on", self.base_url)
            print("  Start the server with: python server.py")
            return
        
        # Run tests
        self.test_validate_smiles()
        self.test_get_molecular_info()
        self.test_visualize_rdkit()
        self.test_visualize_network()
        self.test_visualize_plotly()
        self.test_compare_visualizations()
        self.test_batch_visualize()
        
        print("=" * 60)
        print("Test suite completed!")
        print("=" * 60)

def main():
    """Main function"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Test SMILES Visualizer MCP Server")
    parser.add_argument(
        "--url", "-u",
        default="http://127.0.0.1:8001",
        help="Server URL (default: http://127.0.0.1:8001)"
    )
    
    args = parser.parse_args()
    
    tester = SMILESVisualizerTester(args.url)
    tester.run_all_tests()

if __name__ == "__main__":
    main()
