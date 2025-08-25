#!/bin/bash
set -e

# Function to handle graceful shutdown
cleanup() {
    echo "Received shutdown signal, cleaning up..."
    if [ ! -z "$SERVER_PID" ]; then
        kill -TERM "$SERVER_PID" 2>/dev/null || true
        wait "$SERVER_PID" 2>/dev/null || true
    fi
    exit 0
}

# Set up signal handlers
trap cleanup SIGTERM SIGINT

# Default values
MCP_HOST=${MCP_HOST:-"0.0.0.0"}
MCP_PORT=${MCP_PORT:-"8080"}
OUTPUT_DIR=${OUTPUT_DIR:-"/app/output"}
VERBOSE=${VERBOSE:-"false"}

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Set Python path
export PYTHONPATH="/app:${PYTHONPATH}"

# Check if required packages are installed
echo "Checking dependencies..."
python -c "import mcp; print('✓ MCP SDK available')" || {
    echo "✗ MCP SDK not available"
    exit 1
}

python -c "from rdkit import Chem; print('✓ RDKit available')" || {
    echo "✗ RDKit not available"
    exit 1
}

python -c "import matplotlib; print('✓ Matplotlib available')" || {
    echo "✗ Matplotlib not available"
    exit 1
}

python -c "import plotly; print('✓ Plotly available')" || {
    echo "✗ Plotly not available"
    exit 1
}

python -c "import networkx; print('✓ NetworkX available')" || {
    echo "✗ NetworkX not available"
    exit 1
}

echo "All dependencies available ✓"

# Build server arguments
SERVER_ARGS=(
    "--host" "$MCP_HOST"
    "--port" "$MCP_PORT"
)

if [ "$VERBOSE" = "true" ]; then
    SERVER_ARGS+=("--verbose")
fi

# Start the server
echo "Starting SMILES Visualizer MCP Server..."
echo "Host: $MCP_HOST"
echo "Port: $MCP_PORT"
echo "Output Directory: $OUTPUT_DIR"
echo "Verbose: $VERBOSE"
echo "Server arguments: ${SERVER_ARGS[@]}"
echo ""

# Start the server in background
python server.py "${SERVER_ARGS[@]}" &
SERVER_PID=$!

# Wait for server to start
echo "Waiting for server to start..."
sleep 5

# Check if server is running
if kill -0 "$SERVER_PID" 2>/dev/null; then
    echo "✓ Server started successfully (PID: $SERVER_PID)"
    
    # Health check
    echo "Performing health check. query health endpoint..."
    if curl -f "http://$MCP_HOST:$MCP_PORT/health" >/dev/null 2>&1; then
        echo "✓ Health check passed"
    else
        echo "✗ Health check failed"
        exit 1
    fi
    
    echo ""
    echo "SMILES Visualizer MCP Server is ready!"
    echo "Available endpoints:"
    echo "  - Health: http://$MCP_HOST:$MCP_PORT/health"
    echo ""
    echo "Press Ctrl+C to stop the server"
    
    # Wait for the server process
    wait "$SERVER_PID"
else
    echo "✗ Failed to start server"
    exit 1
fi
