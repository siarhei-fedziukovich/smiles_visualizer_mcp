# SMILES Visualizer MCP Server

A Model Context Protocol (MCP) server for molecular visualization using SMILES (Simplified Molecular Input Line Entry System) strings. This server provides multiple visualization approaches for chemical structures including RDKit, NetworkX, Plotly, and custom matplotlib visualizations.

## Features

- **Multiple Visualization Types**: RDKit 2D structures, network graphs, interactive Plotly charts, and custom matplotlib visualizations
- **Molecular Information**: Detailed molecular properties and descriptors
- **SMILES Validation**: Built-in validation using RDKit
- **Batch Processing**: Process multiple SMILES strings at once
- **HTTP Streamable Transport**: Modern MCP transport for easy integration
- **Base64 Image Output**: Images returned as base64 strings for easy embedding

## Available Tools

### Core Tools

1. **`validate_smiles`** - Validate SMILES strings using RDKit
2. **`get_molecular_info`** - Get detailed molecular information and properties
3. **`visualize_rdkit`** - Create RDKit 2D molecular visualizations
4. **`visualize_network`** - Create network graph visualizations
5. **`visualize_plotly`** - Create interactive Plotly visualizations
6. **`visualize_custom_matplotlib`** - Create custom matplotlib visualizations
7. **`compare_visualizations`** - Generate all visualization types for comparison
8. **`batch_visualize`** - Process multiple SMILES strings

### Molecular Properties

The server calculates various molecular properties including:
- Molecular weight
- Number of atoms, bonds, and rings
- Molecular formula
- LogP and molar refractivity
- Topological polar surface area (TPSA)
- Number of rotatable bonds
- Hydrogen bond donors/acceptors

## Installation

### Prerequisites

- Python 3.8 or higher
- RDKit (for molecular processing)
- Matplotlib (for custom visualizations)
- NetworkX (for network graphs)
- Plotly (for interactive visualizations)

### Setup

1. Clone or download the project
2. Install dependencies:

```bash
pip install -r requirements.txt
```

Or using uv (recommended):

```bash
uv pip install -r requirements.txt
```

## Quick Start

### Option 1: Direct Installation

1. **Install dependencies**
   ```bash
   pip install -r requirements.txt
   ```

2. **Run the server**
   ```bash
   python server.py --host 127.0.0.1 --port 8080 --verbose
   ```

### Option 2: Using uv (Recommended)

1. **Install dependencies with uv**
   ```bash
   uv pip install -r requirements.txt
   ```

2. **Run the server**
   ```bash
   python server.py --host 127.0.0.1 --port 8080 --verbose
   ```

### Option 3: Development Mode (Linux/macOS)

1. **Make entrypoint script executable**
   ```bash
   chmod +x dev_entrypoint.sh
   ```

2. **Run with development entrypoint**
   ```bash
   ./dev_entrypoint.sh
   ```

### Option 4: Docker

1. **Build and run with Docker Compose**
   ```bash
   docker-compose up -d
   ```

2. **Or build manually**
   ```bash
   docker build -t smiles-visualizer-mcp .
   docker run -p 8080:8080 smiles-visualizer-mcp
   ```

3. **Run with custom environment variables**
   ```bash
   docker run -p 8080:8080 \
     -e MCP_HOST=0.0.0.0 \
     -e MCP_PORT=8080 \
     -e VERBOSE=true \
     -e OUTPUT_DIR=/app/output \
     smiles-visualizer-mcp
   ```

## Usage

### Running the Server

Start the MCP server with HTTP Streamable transport:

```bash
python server.py --host 127.0.0.1 --port 8080
```

### Command Line Options

- `--host, -H`: Host address (default: 127.0.0.1)
- `--port, -p`: Port number (default: 8080)
- `--output-dir, -o`: Output directory for files (default: output)
- `--verbose, -v`: Enable verbose logging
- `--version`: Show version information

### Environment Variables

- `MCP_HOST`: Host for HTTP Streamable transport
- `MCP_PORT`: Port for HTTP Streamable transport
- `OUTPUT_DIR`: Directory for saving visualizations

## Tool Examples

### Validate SMILES

```python
# Validate a SMILES string
result = await validate_smiles("CCO")
# Returns: {"valid": true, "message": "Valid SMILES", "canonical_smiles": "CCO"}
```

### Get Molecular Information

```python
# Get detailed molecular information
info = await get_molecular_info("CCO")
# Returns molecular weight, atom count, properties, etc.
```

### Create RDKit Visualization

```python
# Create RDKit 2D visualization
result = await visualize_rdkit("CCO", size="400,300")
# Returns base64 encoded PNG image
```

### Create Network Visualization

```python
# Create network graph
result = await visualize_network("CCO", layout="spring")
# Returns base64 encoded PNG image with network layout
```

### Create Interactive Plotly Visualization

```python
# Create interactive visualization (HTML format)
result = await visualize_plotly("CCO")
# Returns HTML content for interactive chart

# Create interactive visualization (JSON format)
result = await visualize_plotly("CCO", output_format="json")
# Returns JSON data for programmatic use
```

### Compare All Visualizations

```python
# Generate all visualization types (HTML format for Plotly)
results = await compare_visualizations("CCO")
# Returns all visualization types and molecular info

# Generate all visualization types (JSON format for Plotly)
results = await compare_visualizations("CCO", plotly_format="json")
# Returns all visualization types with JSON format for Plotly
```

### Batch Processing

```python
# Process multiple molecules
smiles_list = ["CCO", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"]
results = await batch_visualize(smiles_list, visualization_type="rdkit")
# Returns visualizations for all molecules

# Process multiple molecules with Plotly JSON format
results = await batch_visualize(smiles_list, visualization_type="plotly", plotly_format="json")
# Returns Plotly visualizations in JSON format for all molecules
```

## Example Molecules

The server works with various types of molecules:

- **Simple molecules**: `CCO` (ethanol)
- **Drug molecules**: `CC(C)CC1=CC=C(C=C1)C(C)C(=O)O` (ibuprofen)
- **Complex structures**: `CN1C=NC2=C1C(=O)N(C(=O)N2C)C` (caffeine)
- **Aromatic compounds**: `C1=CC=C(C=C1)C2=CC=CC=C2` (biphenyl)

## Integration with MCP Clients

This server can be integrated with any MCP-compatible client such as:

- Claude Desktop
- VS Code with MCP extension
- Custom MCP clients

### Client Configuration

Add the server to your MCP client configuration:

```json
{
  "mcpServers": {
    "smiles-visualizer": {
      "command": "python",
      "args": ["path/to/smiles_visualizer_mcp/server.py", "--host", "127.0.0.1", "--port", "8080"],
      "env": {
        "MCP_HOST": "127.0.0.1",
        "MCP_PORT": "8080"
      }
    }
  }
}
```

## Output Formats

### Images
- **Format**: PNG
- **Encoding**: Base64
- **Usage**: Can be embedded in HTML, displayed in applications, or saved to files

### Interactive Visualizations
- **Format**: HTML (Plotly) or JSON (Plotly)
- **Features**: Zoom, pan, hover information, interactive elements
- **JSON Format**: Raw Plotly figure data for programmatic use

### Data
- **Format**: JSON
- **Content**: Molecular properties, validation results, error messages

## Error Handling

The server includes comprehensive error handling:

- Invalid SMILES strings
- Missing dependencies
- Processing errors
- Network/graph generation issues

All errors are returned as structured JSON responses with descriptive messages.

## Dependencies

### Required
- `mcp[cli]` - MCP Python SDK
- `rdkit-pypi` - Molecular processing
- `matplotlib` - Custom visualizations
- `networkx` - Network graphs
- `plotly` - Interactive visualizations
- `numpy` - Numerical operations
- `pandas` - Data manipulation
- `pillow` - Image processing
- `uvicorn` - ASGI server
- `fastapi` - Web framework

### Optional
- `seaborn` - Enhanced plotting (if available)

## Development

### Project Structure

```
smiles_visualizer_mcp/
├── server.py          # Main MCP server implementation
├── requirements.txt   # Python dependencies
└── README.md         # This file
```

### Adding New Visualizations

To add new visualization types:

1. Add the tool decorator to the `setup_tools` method
2. Implement the visualization logic
3. Return results in the expected JSON format
4. Update the `compare_visualizations` method to include the new type

### Testing

Test the server with example SMILES strings:

```bash
# Start the server
python server.py

# Test with curl (in another terminal)
curl -X POST http://127.0.0.1:8080/tools/validate_smiles/call \
  -H "Content-Type: application/json" \
  -d '{"arguments": {"smiles": "CCO"}}'

## License

This project is open source and available under the MIT License.

## Contributing

Contributions are welcome! Please feel free to submit pull requests or open issues for bugs and feature requests.

## Support

For issues and questions:
1. Check the error messages in the server logs
2. Verify all dependencies are installed correctly
3. Ensure SMILES strings are valid
4. Check that the required libraries (RDKit, matplotlib, etc.) are available

## Related Projects

- [RDKit](https://www.rdkit.org/) - Open-source cheminformatics toolkit
- [Model Context Protocol](https://modelcontextprotocol.io/) - Protocol for AI context
- [MCP Python SDK](https://github.com/modelcontextprotocol/python-sdk) - Official Python implementation
