#!/bin/bash
set -e

# Development entrypoint script
echo "SMILES Visualizer MCP Server - Development Mode"
echo "================================================"

# Default values for development
MCP_HOST=${MCP_HOST:-"127.0.0.1"}
MCP_PORT=${MCP_PORT:-"8080"}
OUTPUT_DIR=${OUTPUT_DIR:-"./output"}
VERBOSE=${VERBOSE:-"true"}

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Set Python path for development
export PYTHONPATH=".:${PYTHONPATH}"

# Set matplotlib backend for headless environment
export MPLBACKEND=Agg

echo "Development Configuration:"
echo "  Host: $MCP_HOST"
echo "  Port: $MCP_PORT"
echo "  Output Directory: $OUTPUT_DIR"
echo "  Verbose: $VERBOSE"
echo "  Python Path: $PYTHONPATH"
echo ""

# Check dependencies
echo "Checking dependencies..."
python -c "import mcp; print('✓ MCP SDK')" || echo "✗ MCP SDK missing"
python -c "from rdkit import Chem; print('✓ RDKit')" || echo "✗ RDKit missing"
python -c "import matplotlib; print('✓ Matplotlib')" || echo "✗ Matplotlib missing"
python -c "import plotly; print('✓ Plotly')" || echo "✗ Plotly missing"
python -c "import networkx; print('✓ NetworkX')" || echo "✗ NetworkX missing"
echo ""

# Build server arguments
SERVER_ARGS=(
    "--host" "$MCP_HOST"
    "--port" "$MCP_PORT"
)

if [ "$VERBOSE" = "true" ]; then
    SERVER_ARGS+=("--verbose")
fi

echo "Starting server with arguments: ${SERVER_ARGS[@]}"
echo ""

# Start the server directly (no background process for development)
exec python server.py "${SERVER_ARGS[@]}"
