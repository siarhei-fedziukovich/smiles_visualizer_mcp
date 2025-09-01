#!/usr/bin/env python3
"""
SMILES Visualizer MCP Server
A Model Context Protocol server for SMILES molecular visualization
"""

import asyncio
import argparse
import json
import logging
import os
import sys
import base64
from typing import Any, Dict, List, Optional, Sequence
from pathlib import Path

from mcp.server.fastmcp import FastMCP
from mcp.types import TextContent, ImageContent

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Configure matplotlib for headless environment
import matplotlib
matplotlib.use('Agg')

# Import visualization modules
try:
    from rdkit import Chem
    from rdkit.Chem import Draw, AllChem, rdMolDescriptors
    from rdkit.Chem.Draw import rdMolDraw2D
    RDKIT_AVAILABLE = True
except ImportError as e:
    RDKIT_AVAILABLE = False
    logger.warning(f"RDKit not available. RDKit visualizations will be disabled. ex: {e}")

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from matplotlib.patches import FancyBboxPatch
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    logger.warning("Matplotlib not available. Matplotlib visualizations will be disabled.")

try:
    import networkx as nx
    NETWORKX_AVAILABLE = True
except ImportError:
    NETWORKX_AVAILABLE = False
    logger.warning("NetworkX not available. Network visualizations will be disabled.")

try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    logger.warning("Plotly not available. Plotly visualizations will be disabled.")

import numpy as np
import pandas as pd
from PIL import Image
import io

class SMILESVisualizerMCP:
    """SMILES Visualizer MCP Server implementation using FastMCP"""
    
    def __init__(self, output_dir: str = "output"):
        self.mcp = FastMCP("SMILES Visualizer")
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.setup_tools()
    
    def setup_tools(self):
        """Setup all SMILES visualization tools using FastMCP decorators"""
        
        @self.mcp.tool()
        async def validate_smiles(smiles: str) -> str:
            """Validate a SMILES string using RDKit"""
            if not RDKIT_AVAILABLE:
                return json.dumps({"valid": True, "message": "RDKit not available, assuming valid"})
            
            try:
                mol = Chem.MolFromSmiles(smiles)
                is_valid = mol is not None
                return json.dumps({
                    "valid": is_valid,
                    "message": "Valid SMILES" if is_valid else "Invalid SMILES",
                    "canonical_smiles": Chem.MolToSmiles(mol, canonical=True) if is_valid else None
                })
            except Exception as e:
                return json.dumps({"valid": False, "message": f"Error: {str(e)}"})

        @self.mcp.tool()
        async def get_molecular_info(smiles: str) -> str:
            """Get detailed molecular information from SMILES"""
            if not RDKIT_AVAILABLE:
                return json.dumps({"error": "RDKit not available"})
            
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    return json.dumps({"error": "Invalid SMILES"})
                
                info = {
                    "molecular_weight": rdMolDescriptors.CalcExactMolWt(mol),
                    "num_atoms": mol.GetNumAtoms(),
                    "num_bonds": mol.GetNumBonds(),
                    "num_rings": rdMolDescriptors.CalcNumRings(mol),
                    "formula": rdMolDescriptors.CalcMolFormula(mol),
                    "smiles": Chem.MolToSmiles(mol),
                    "canonical_smiles": Chem.MolToSmiles(mol, canonical=True),
                    "logp": rdMolDescriptors.CalcCrippenDescriptors(mol)[0],
                    "molar_refractivity": rdMolDescriptors.CalcCrippenDescriptors(mol)[1],
                    "tpsa": rdMolDescriptors.CalcTPSA(mol),
                    "num_rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
                    "num_hbd": rdMolDescriptors.CalcNumHBD(mol),
                    "num_hba": rdMolDescriptors.CalcNumHBA(mol)
                }
                return json.dumps(info, indent=2)
            except Exception as e:
                return json.dumps({"error": str(e)})

        @self.mcp.tool()
        async def visualize_rdkit(smiles: str, size: str = "200,200") -> list:
            """Create RDKit 2D molecular visualization and return as image content"""
            if not RDKIT_AVAILABLE:
                return [TextContent(type="text", text="RDKit is required for this visualization")]
            
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    return [TextContent(type="text", text=f"Invalid SMILES string: {smiles}")]
                
                # Parse size
                width, height = map(int, size.split(','))
                
                # Create 2D visualization
                img = Draw.MolToImage(mol, size=(width, height))
                
                # Convert to base64
                buffer = io.BytesIO()
                img.save(buffer, format='PNG')
                img_base64 = base64.b64encode(buffer.getvalue()).decode()
                
                return [
                    TextContent(type="text", text=f"RDKit visualization for {smiles} ({width}x{height})"),
                    ImageContent(
                        type="image",
                        data=img_base64,
                        mimeType="image/png"
                    )
                ]
            except Exception as e:
                logger.error(f"Error in visualize_rdkit for SMILES '{smiles}': {str(e)}", exc_info=True)
                return [TextContent(type="text", text=f"Error: {str(e)}")]

        @self.mcp.tool()
        async def visualize_network(smiles: str, layout: str = "spring") -> list:
            """Create network graph visualization of molecular structure"""
            if not NETWORKX_AVAILABLE or not RDKIT_AVAILABLE:
                return [TextContent(type="text", text="NetworkX and RDKit are required for this visualization")]
            
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    return [TextContent(type="text", text=f"Invalid SMILES string: {smiles}")]
                
                # Create graph
                G = nx.Graph()
                
                # Add nodes (atoms)
                for atom in mol.GetAtoms():
                    atom_idx = atom.GetIdx()
                    atom_symbol = atom.GetSymbol()
                    G.add_node(atom_idx, symbol=atom_symbol, atomic_num=atom.GetAtomicNum())
                
                # Add edges (bonds)
                for bond in mol.GetBonds():
                    begin_idx = bond.GetBeginAtomIdx()
                    end_idx = bond.GetEndAtomIdx()
                    bond_type = bond.GetBondType()
                    G.add_edge(begin_idx, end_idx, bond_type=bond_type)
                
                # Create visualization
                plt.figure(figsize=(12, 8))
                
                if layout == "spring":
                    pos = nx.spring_layout(G, k=1, iterations=50)
                elif layout == "circular":
                    pos = nx.circular_layout(G)
                elif layout == "random":
                    pos = nx.random_layout(G)
                else:
                    pos = nx.spring_layout(G, k=1, iterations=50)
                
                # Draw nodes
                nx.draw_networkx_nodes(G, pos, node_color='lightblue', 
                                      node_size=1000, alpha=0.7)
                
                # Draw edges
                nx.draw_networkx_edges(G, pos, alpha=0.5, edge_color='gray')
                
                # Add atom labels
                labels = {node: G.nodes[node]['symbol'] for node in G.nodes()}
                nx.draw_networkx_labels(G, pos, labels, font_size=12, font_weight='bold')
                
                plt.title(f"Molecular Network Graph: {smiles}")
                plt.axis('off')
                
                # Save to buffer and convert to base64
                buffer = io.BytesIO()
                plt.savefig(buffer, dpi=300, bbox_inches='tight', format='PNG')
                plt.close()
                
                img_base64 = base64.b64encode(buffer.getvalue()).decode()
                
                return [
                    TextContent(type="text", text=f"Network visualization for {smiles} ({layout} layout, {G.number_of_nodes()} nodes, {G.number_of_edges()} edges)"),
                    ImageContent(
                        type="image",
                        data=img_base64,
                        mimeType="image/png"
                    )
                ]
            except Exception as e:
                logger.error(f"Error in visualize_network for SMILES '{smiles}' with layout '{layout}': {str(e)}", exc_info=True)
                return [TextContent(type="text", text=f"Error: {str(e)}")]

        @self.mcp.tool()
        async def visualize_plotly(smiles: str) -> list:
            """Create interactive Plotly visualization and return as JSON"""
            if not PLOTLY_AVAILABLE or not RDKIT_AVAILABLE:
                return [TextContent(type="text", text="Plotly and RDKit are required for this visualization")]
            
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    return [TextContent(type="text", text=f"Invalid SMILES string: {smiles}")]
                
                mol_info = json.loads(await get_molecular_info(smiles))
                
                # Create subplots
                fig = make_subplots(
                    rows=2, cols=2,
                    subplot_titles=('Molecular Structure', 'Molecular Properties', 
                                  'Atom Distribution', 'Bond Types'),
                    specs=[[{"type": "scatter"}, {"type": "bar"}],
                           [{"type": "pie"}, {"type": "bar"}]]
                )
                
                # 1. Molecular structure (2D coordinates)
                AllChem.Compute2DCoords(mol)
                conf = mol.GetConformer()
                
                x_coords = []
                y_coords = []
                atom_labels = []
                
                for atom in mol.GetAtoms():
                    pos = conf.GetAtomPosition(atom.GetIdx())
                    x_coords.append(pos.x)
                    y_coords.append(pos.y)
                    atom_labels.append(atom.GetSymbol())
                
                fig.add_trace(
                    go.Scatter(x=x_coords, y=y_coords, mode='markers+text',
                              text=atom_labels, textposition="middle center",
                              marker=dict(size=20, color='lightblue'),
                              name='Atoms'),
                    row=1, col=1
                )
                
                # Add bonds
                for bond in mol.GetBonds():
                    begin_pos = conf.GetAtomPosition(bond.GetBeginAtomIdx())
                    end_pos = conf.GetAtomPosition(bond.GetEndAtomIdx())
                    
                    fig.add_trace(
                        go.Scatter(x=[begin_pos.x, end_pos.x], 
                                  y=[begin_pos.y, end_pos.y],
                                  mode='lines', line=dict(color='gray', width=2),
                                  showlegend=False),
                        row=1, col=1
                    )
                
                # 2. Molecular properties
                if 'error' not in mol_info:
                    properties = ['Molecular Weight', 'Number of Atoms', 'Number of Bonds', 'Number of Rings']
                    values = [mol_info['molecular_weight'], mol_info['num_atoms'], 
                             mol_info['num_bonds'], mol_info['num_rings']]
                    
                    fig.add_trace(
                        go.Bar(x=properties, y=values, name='Properties'),
                        row=1, col=2
                    )
                
                # 3. Atom distribution
                atom_counts = {}
                for atom in mol.GetAtoms():
                    symbol = atom.GetSymbol()
                    atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
                
                fig.add_trace(
                    go.Pie(labels=list(atom_counts.keys()), 
                          values=list(atom_counts.values()),
                          name='Atoms'),
                    row=2, col=1
                )
                
                # 4. Bond types
                bond_counts = {}
                for bond in mol.GetBonds():
                    bond_type = str(bond.GetBondType())
                    bond_counts[bond_type] = bond_counts.get(bond_type, 0) + 1
                
                fig.add_trace(
                    go.Bar(x=list(bond_counts.keys()), 
                          y=list(bond_counts.values()),
                          name='Bonds'),
                    row=2, col=2
                )
                
                fig.update_layout(
                    title=f"Interactive Molecular Analysis: {smiles}",
                    height=800,
                    showlegend=True
                )
                
                # Return as JSON data
                plotly_json = fig.to_json()
                return [
                    TextContent(type="text", text=f"Interactive Plotly visualization for {smiles} (JSON format)"),
                    ImageContent(
                        type="image",
                        data=plotly_json,
                        mimeType="application/vnd.plotly.v1+json"
                    )
                ]
            except Exception as e:
                logger.error(f"Error in visualize_plotly for SMILES '{smiles}': {str(e)}", exc_info=True)
                return [TextContent(type="text", text=f"Error: {str(e)}")]

        @self.mcp.tool()
        async def visualize_custom_matplotlib(smiles: str) -> list:
            """Create custom matplotlib visualization with molecular properties"""
            if not MATPLOTLIB_AVAILABLE or not RDKIT_AVAILABLE:
                return [TextContent(type="text", text="Matplotlib and RDKit are required for this visualization")]
            
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    return [TextContent(type="text", text=f"Invalid SMILES string: {smiles}")]
                
                mol_info = json.loads(await get_molecular_info(smiles))
                
                # Create figure with subplots
                fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
                
                # 1. Molecular structure (RDKit image)
                img = Draw.MolToImage(mol, size=(300, 300))
                ax1.imshow(img)
                ax1.set_title('Molecular Structure')
                ax1.axis('off')
                
                # 2. Molecular properties
                if 'error' not in mol_info:
                    properties = ['MW', 'Atoms', 'Bonds', 'Rings']
                    values = [mol_info['molecular_weight'], mol_info['num_atoms'], 
                             mol_info['num_bonds'], mol_info['num_rings']]
                    
                    bars = ax2.bar(properties, values, color=['skyblue', 'lightgreen', 'lightcoral', 'gold'])
                    ax2.set_title('Molecular Properties')
                    ax2.set_ylabel('Count/Weight')
                    
                    # Add value labels on bars
                    for bar, value in zip(bars, values):
                        height = bar.get_height()
                        ax2.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                                f'{value:.1f}', ha='center', va='bottom')
                
                # 3. Atom distribution
                atom_counts = {}
                for atom in mol.GetAtoms():
                    symbol = atom.GetSymbol()
                    atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
                
                if atom_counts:
                    ax3.pie(atom_counts.values(), labels=atom_counts.keys(), autopct='%1.1f%%')
                    ax3.set_title('Atom Distribution')
                
                # 4. Bond types
                bond_counts = {}
                for bond in mol.GetBonds():
                    bond_type = str(bond.GetBondType())
                    bond_counts[bond_type] = bond_counts.get(bond_type, 0) + 1
                
                if bond_counts:
                    ax4.bar(bond_counts.keys(), bond_counts.values(), color='lightblue')
                    ax4.set_title('Bond Types')
                    ax4.set_ylabel('Count')
                
                plt.suptitle(f'Molecular Analysis: {smiles}', fontsize=16)
                plt.tight_layout()
                
                # Save to buffer and convert to base64
                buffer = io.BytesIO()
                plt.savefig(buffer, dpi=300, bbox_inches='tight', format='PNG')
                plt.close()
                
                img_base64 = base64.b64encode(buffer.getvalue()).decode()
                
                return [
                    TextContent(type="text", text=f"Custom matplotlib visualization for {smiles}"),
                    ImageContent(
                        type="image",
                        data=img_base64,
                        mimeType="image/png"
                    )
                ]
            except Exception as e:
                logger.error(f"Error in visualize_custom_matplotlib for SMILES '{smiles}': {str(e)}", exc_info=True)
                return [TextContent(type="text", text=f"Error: {str(e)}")]

        @self.mcp.tool()
        async def compare_visualizations(smiles: str, plotly_format: str = "html") -> list:
            """Generate all visualization types for comparison"""
            # Validate SMILES first
            validation = json.loads(await validate_smiles(smiles))
            if not validation.get("valid", False):
                return [TextContent(type="text", text=f"Invalid SMILES: {smiles}")]
            
            # Get molecular info
            mol_info = json.loads(await get_molecular_info(smiles))
            
            # Generate visualizations
            viz_methods = [
                ("rdkit", visualize_rdkit),
                ("network", visualize_network),
                ("plotly", visualize_plotly),
                ("custom_matplotlib", visualize_custom_matplotlib)
            ]
            
            content_items = [
                TextContent(type="text", text=f"Comparison of visualization types for {smiles}"),
                TextContent(type="text", text=f"Molecular Info: {json.dumps(mol_info, indent=2)}")
            ]
            
            for viz_type, method in viz_methods:
                try:
                    if viz_type == "plotly":
                        result = await method(smiles, plotly_format)
                    else:
                        result = await method(smiles)
                    content_items.append(TextContent(type="text", text=f"\n--- {viz_type.upper()} VISUALIZATION ---"))
                    content_items.extend(result)
                except Exception as e:
                    logger.error(f"Error in compare_visualizations for {viz_type} visualization of SMILES '{smiles}': {str(e)}", exc_info=True)
                    content_items.append(TextContent(type="text", text=f"\n--- {viz_type.upper()} VISUALIZATION (ERROR) ---"))
                    content_items.append(TextContent(type="text", text=f"Error: {str(e)}"))
            
            return content_items

        @self.mcp.tool()
        async def batch_visualize(smiles_list: List[str], visualization_type: str = "rdkit", plotly_format: str = "html") -> list:
            """Generate visualizations for multiple SMILES strings"""
            if not smiles_list:
                return [TextContent(type="text", text="No SMILES strings provided")]
            
            content_items = [
                TextContent(type="text", text=f"Batch visualization of {len(smiles_list)} molecules using {visualization_type}")
            ]
            
            for i, smiles in enumerate(smiles_list):
                try:
                    if visualization_type == "rdkit":
                        result = await visualize_rdkit(smiles)
                    elif visualization_type == "network":
                        result = await visualize_network(smiles)
                    elif visualization_type == "plotly":
                        result = await visualize_plotly(smiles, plotly_format)
                    elif visualization_type == "custom_matplotlib":
                        result = await visualize_custom_matplotlib(smiles)
                    else:
                        return [TextContent(type="text", text=f"Unknown visualization type: {visualization_type}")]
                    
                    content_items.append(TextContent(type="text", text=f"\n--- MOLECULE {i+1}: {smiles} ---"))
                    content_items.extend(result)
                except Exception as e:
                    logger.error(f"Error in batch_visualize for molecule {i+1} with SMILES '{smiles}' using {visualization_type}: {str(e)}", exc_info=True)
                    content_items.append(TextContent(type="text", text=f"\n--- MOLECULE {i+1}: {smiles} (ERROR) ---"))
                    content_items.append(TextContent(type="text", text=f"Error: {str(e)}"))
            
            return content_items

async def run_http_streamable_server(smiles_server: SMILESVisualizerMCP, host: str, port: int):
    """Run the MCP server using HTTP Streamable transport"""
    logger.info(f"Starting SMILES Visualizer MCP Server with HTTP Streamable transport on {host}:{port}")
    
    import uvicorn
    app = smiles_server.mcp.streamable_http_app()
    config = uvicorn.Config(app, host=host, port=port, log_level="info")
    server = uvicorn.Server(config)
    await server.serve()

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="SMILES Visualizer MCP Server - A Model Context Protocol server for molecular visualization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with HTTP Streamable transport (default)
  python server.py --host 127.0.0.1 --port 8001

  # Run with different host and port
  python server.py --host 0.0.0.0 --port 8002

Environment Variables:
  MCP_HOST: Host for HTTP Streamable transport (default: 127.0.0.1)
  MCP_PORT: Port for HTTP Streamable transport (default: 8001)
  OUTPUT_DIR: Directory for saving visualizations (default: output)
        """
    )
    
    parser.add_argument(
        "--host", "-H",
        default=os.getenv("MCP_HOST", "127.0.0.1"),
        help="Host for HTTP Streamable transport (default: 127.0.0.1)"
    )
    
    parser.add_argument(
        "--port", "-p",
        type=int,
        default=int(os.getenv("MCP_PORT", "8080")),
        help="Port for HTTP Streamable transport (default: 8080)"
    )
    
    parser.add_argument(
        "--output-dir", "-o",
        default=os.getenv("OUTPUT_DIR", "output"),
        help="Directory for saving visualizations (default: output)"
    )
    
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging"
    )
    
    parser.add_argument(
        "--version",
        action="version",
        version="SMILES Visualizer MCP Server 1.0.0"
    )
    
    return parser.parse_args()

async def main():
    """Main function to run the MCP server"""
    args = parse_arguments()
    
    # Configure logging level
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    logger.info(f"SMILES Visualizer MCP Server Configuration:")
    logger.info(f"  Host: {args.host}")
    logger.info(f"  Port: {args.port}")
    logger.info(f"  Output Directory: {args.output_dir}")
    logger.info(f"  RDKit Available: {RDKIT_AVAILABLE}")
    logger.info(f"  Matplotlib Available: {MATPLOTLIB_AVAILABLE}")
    logger.info(f"  NetworkX Available: {NETWORKX_AVAILABLE}")
    logger.info(f"  Plotly Available: {PLOTLY_AVAILABLE}")
    
    try:
        # Create and run the server
        smiles_server = SMILESVisualizerMCP(output_dir=args.output_dir)
        await run_http_streamable_server(smiles_server, args.host, args.port)
    except KeyboardInterrupt:
        logger.info("Server stopped by user")
    except Exception as e:
        logger.error(f"Server error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    asyncio.run(main())
