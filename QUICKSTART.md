# Quick Start Guide - SMILES Visualizer MCP Server

This guide will help you get the SMILES Visualizer MCP Server up and running quickly.

## Prerequisites

- Python 3.8 or higher
- pip or uv package manager

## Installation

### Option 1: Direct Installation

1. **Clone or download the project**
   ```bash
   git clone <repository-url>
   cd smiles_visualizer_mcp
   ```

2. **Install dependencies**
   ```bash
   pip install -r requirements.txt
   ```

### Option 2: Using uv (Recommended)

1. **Install uv if you haven't already**
   ```bash
   pip install uv
   ```

2. **Install dependencies with uv**
   ```bash
   uv pip install -r requirements.txt
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

## Running the Server

### Basic Usage

Start the server with default settings:
```bash
python server.py
```

The server will start on `http://127.0.0.1:8080`

### Custom Configuration

Start with custom host and port:
```bash
python server.py --host 0.0.0.0 --port 8002
```

### Environment Variables

You can also use environment variables:
```bash
export MCP_HOST=0.0.0.0
export MCP_PORT=8002
python server.py
```

## Testing the Server

### Run the Test Suite

```bash
python test_server.py
```

### Manual Testing

Test with curl:
```bash
# Test SMILES validation
curl -X POST http://127.0.0.1:8080/tools/validate_smiles/call \
  -H "Content-Type: application/json" \
  -d '{"arguments": {"smiles": "CCO"}}'

# Test molecular info
curl -X POST http://127.0.0.1:8080/tools/get_molecular_info/call \
  -H "Content-Type: application/json" \
  -d '{"arguments": {"smiles": "CCO"}}'

# Test RDKit visualization
curl -X POST http://127.0.0.1:8080/tools/visualize_rdkit/call \
  -H "Content-Type: application/json" \
  -d '{"arguments": {"smiles": "CCO", "size": "400,300"}}'
```

## Example Usage

### Python Client Example

```python
import requests
import json
import base64

# Server URL
base_url = "http://127.0.0.1:8080"

# Test SMILES validation
def validate_smiles(smiles):
    response = requests.post(
        f"{base_url}/tools/validate_smiles/call",
        json={"arguments": {"smiles": smiles}},
        headers={"Content-Type": "application/json"}
    )
    return response.json()

# Test molecular info
def get_molecular_info(smiles):
    response = requests.post(
        f"{base_url}/tools/get_molecular_info/call",
        json={"arguments": {"smiles": smiles}},
        headers={"Content-Type": "application/json"}
    )
    return response.json()

# Test visualization
def visualize_molecule(smiles):
    response = requests.post(
        f"{base_url}/tools/visualize_rdkit/call",
        json={"arguments": {"smiles": smiles, "size": "400,300"}},
        headers={"Content-Type": "application/json"}
    )
    return response.json()

# Example usage
if __name__ == "__main__":
    # Test with ethanol
    smiles = "CCO"
    
    print("Validating SMILES...")
    validation = validate_smiles(smiles)
    print(json.dumps(validation, indent=2))
    
    print("\nGetting molecular info...")
    info = get_molecular_info(smiles)
    print(json.dumps(info, indent=2))
    
    print("\nCreating visualization...")
    viz = visualize_molecule(smiles)
    print(json.dumps(viz, indent=2))
```

## Integration with MCP Clients

### Claude Desktop

Add to your Claude Desktop configuration:

```json
{
  "mcpServers": {
    "smiles-visualizer": {
      "command": "python",
      "args": ["/path/to/smiles_visualizer_mcp/server.py"],
      "env": {
        "MCP_HOST": "127.0.0.1",
        "MCP_PORT": "8080"
      }
    }
  }
}
```

### VS Code

Add to your VS Code settings:

```json
{
  "mcp.servers": {
    "smiles-visualizer": {
      "command": "python",
      "args": ["/path/to/smiles_visualizer_mcp/server.py"],
      "env": {
        "MCP_HOST": "127.0.0.1",
        "MCP_PORT": "8080"
      }
    }
  }
}
```

## Available Tools

The server provides these tools:

1. **`validate_smiles`** - Validate SMILES strings
2. **`get_molecular_info`** - Get molecular properties
3. **`visualize_rdkit`** - Create RDKit 2D visualizations
4. **`visualize_network`** - Create network graphs
5. **`visualize_plotly`** - Create interactive visualizations
6. **`visualize_custom_matplotlib`** - Create custom matplotlib plots
7. **`compare_visualizations`** - Generate all visualization types
8. **`batch_visualize`** - Process multiple molecules

## Example Molecules

Try these SMILES strings:

- `CCO` - Ethanol
- `CC(C)CC1=CC=C(C=C1)C(C)C(=O)O` - Ibuprofen
- `CN1C=NC2=C1C(=O)N(C(=O)N2C)C` - Caffeine
- `CC(=O)OC1=CC=CC=C1C(=O)O` - Aspirin
- `C1=CC=C(C=C1)C2=CC=CC=C2` - Biphenyl

## Troubleshooting

### Common Issues

1. **Port already in use**
   ```bash
   # Use a different port
   python server.py --port 8002
   ```

2. **Dependencies not installed**
   ```bash
   # Reinstall dependencies
   pip install -r requirements.txt
   ```

3. **RDKit installation issues**
   ```bash
   # Try conda installation
   conda install -c conda-forge rdkit
   ```

4. **Permission denied**
   ```bash
   # Use a different host
   python server.py --host 127.0.0.1
   ```

### Logs

Enable verbose logging:
```bash
python server.py --verbose
```

### Health Check

Check if the server is running:
```bash
curl http://127.0.0.1:8080/tools
```

## Next Steps

1. **Explore the tools** - Try different visualization types
2. **Test with your molecules** - Use your own SMILES strings
3. **Integrate with your workflow** - Connect to your preferred MCP client
4. **Customize** - Modify the server for your specific needs

## Support

- Check the main README.md for detailed documentation
- Run `python example.py` for usage examples
- Use `python test_server.py` to verify functionality
