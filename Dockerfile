FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    libc6-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements first for better caching
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy application files
COPY server.py .
COPY test_server.py .
COPY example.py .
COPY README.md .
COPY entrypoint.sh .

# Create output directory
RUN mkdir -p output

# Make entrypoint script executable
RUN chmod +x entrypoint.sh

# Expose port
EXPOSE 8080

# Set environment variables
ENV MCP_HOST=0.0.0.0
ENV MCP_PORT=8080
ENV OUTPUT_DIR=/app/output
ENV VERBOSE=false

# Use entrypoint script
ENTRYPOINT ["./entrypoint.sh"]
