# Use Python 3.10 slim image for smaller size
FROM python:3.10-slim

# Set working directory
WORKDIR /app

# Copy requirements first for better caching
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy application source code
COPY ribostruct/ ./ribostruct/
COPY requirements.txt .

# Create jobs directory
RUN mkdir -p jobs

# Expose port 8000
EXPOSE 8000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD python -c "import os, requests; port = os.environ.get('PORT', '8000'); requests.get(f'http://localhost:{port}/health')" || exit 1

# Run the application (Render dynamically assigns $PORT)
CMD sh -c "uvicorn ribostruct.web.server:app --host 0.0.0.0 --port ${PORT:-8000}"
