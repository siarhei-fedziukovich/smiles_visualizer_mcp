FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Set noninteractive installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    libc6-dev \
    libxrender1 \
    libxext6 \
    libx11-6 \
    libxrandr2 \
    libxinerama1 \
    libxi6 \
    libxss1 \
    libxtst6 \
    libxfixes3 \
    libxcb1 \
    libxcb-render0 \
    libxcb-shape0 \
    libxcb-xfixes0 \
    libxcb-keysyms1 \
    libxcb-icccm4 \
    libxcb-image0 \
    libxcb-shm0 \
    libxcb-util1 \
    libxcb-randr0 \
    libxcb-xinerama0 \
    libxcb-xkb1 \
    libxkbcommon0 \
    libxkbcommon-x11-0 \
    libfontconfig1 \
    libfreetype6 \
    libpng16-16 \
    libjpeg62-turbo \
    libtiff5 \
    libwebp6 \
    libharfbuzz0b \
    libfribidi0 \
    libexpat1 \
    libzstd1 \
    liblz4-1 \
    libbz2-1.0 \
    liblzma5 \
    libgomp1 \
    libgcc-s1 \
    libstdc++6 \
    libc6 \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements first for better caching
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy application files
COPY server.py .
COPY test_server.py .
COPY example.py .
COPY test_rdkit_docker.py .
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
ENV MPLBACKEND=Agg

# Use entrypoint script
ENTRYPOINT ["./entrypoint.sh"]
