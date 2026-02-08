# Docker Deployment Guide

## Quick Start

### Option 1: Docker Compose (Recommended)

```bash
# Build and start the container
docker-compose up -d

# View logs
docker-compose logs -f

# Stop the container
docker-compose down
```

The application will be available at: **http://localhost:8000**

### Option 2: Docker CLI

```bash
# Build the image
docker build -t ribostruct .

# Run the container
docker run -d -p 8000:8000 --name ribostruct ribostruct

# View logs
docker logs -f ribostruct

# Stop and remove
docker stop ribostruct
docker rm ribostruct
```

---

## Container Details

### Base Image
- **python:3.10-slim** - Lightweight Python runtime

### Exposed Ports
- **8000** - Web server (uvicorn)

### Volumes
The `docker-compose.yml` automatically mounts:
- `./jobs:/app/jobs` - Job data persists on host

### Health Check
The container includes automatic health monitoring:
- Checks `/health` endpoint every 30 seconds
- 3 retries with 10-second timeout
- 5-second startup grace period

---

## Environment Variables

You can customize behavior with environment variables:

```yaml
environment:
  - PYTHONUNBUFFERED=1  # Real-time logging
```

---

## Production Deployment

### With Volume Persistence

```bash
docker run -d \
  -p 8000:8000 \
  -v $(pwd)/jobs:/app/jobs \
  --name ribostruct \
  --restart unless-stopped \
  ribostruct
```

### Behind a Reverse Proxy (Nginx)

```nginx
server {
    listen 80;
    server_name ribostruct.example.com;

    location / {
        proxy_pass http://localhost:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }
}
```

---

## Maintenance

### View Container Logs
```bash
docker logs -f ribostruct
```

### Access Container Shell
```bash
docker exec -it ribostruct /bin/bash
```

### Update to Latest Version
```bash
# Pull new code
git pull

# Rebuild and restart
docker-compose down
docker-compose up -d --build
```

### Clean Up Old Images
```bash
docker image prune -a
```

---

## Verification

After starting the container, verify it's working:

```bash
# Check health
curl http://localhost:8000/health

# Expected output:
# {
#   "status": "healthy",
#   "jobs_directory": "jobs",
#   "jobs_directory_exists": true
# }
```

Open browser: **http://localhost:8000**

---

## Troubleshooting

### Container won't start
```bash
# Check logs
docker logs ribostruct

# Verify image built correctly
docker images | grep ribostruct
```

### Port already in use
```bash
# Find process using port 8000
lsof -i :8000

# Use different port
docker run -p 8001:8000 ribostruct
```

### Jobs directory permissions
```bash
# On Linux, ensure proper permissions
chmod -R 755 ./jobs
```

---

## What's Included

The Docker image contains:
- ✅ Python 3.10 + all dependencies
- ✅ Core bioinformatics modules (parser, alignment, processor, injector)
- ✅ FastAPI web server
- ✅ HTML frontend
- ✅ Automated job cleanup (24-hour retention)
- ✅ Health monitoring

**Image Size**: ~500 MB (slim, production-ready)

---

## Security Notes

- Container runs as non-root (safe)
- No sensitive data in image
- Jobs automatically cleaned after 24 hours
- Health checks prevent zombie processes

---

## Cloud Deployment

### AWS ECS
Upload to ECR and deploy as a service

### Google Cloud Run
```bash
gcloud run deploy ribostruct --source .
```

### Azure Container Instances
```bash
az container create --resource-group myResourceGroup --name ribostruct --image ribostruct:latest --ports 8000
```

---

For more information, see [SERVER_USAGE.md](SERVER_USAGE.md)
