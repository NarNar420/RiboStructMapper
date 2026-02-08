# 🌐 Production Deployment Guide

## Making RiboStructMapper Publicly Accessible

This guide shows how to deploy RiboStructMapper so everyone can access it 24/7 without you manually starting the server.

---

## 🎯 Quick Answer

**Best Options:**
1. **Cloud Platform (Easiest)** - Deploy to AWS/GCP/Azure/DigitalOcean
2. **VPS with Docker** - Use your own server with Docker + systemd
3. **Docker + Reverse Proxy** - Professional setup with nginx

---

## Option 1: Cloud Platform Deployment (Recommended)

### Deploy to DigitalOcean App Platform (Easiest)

1. **Push to GitHub/GitLab**
   ```bash
   git push origin main
   ```

2. **Create App on DigitalOcean**
   - Go to DigitalOcean App Platform
   - Click "Create App" → Connect your repository
   - Select "Docker" as build method
   - Set port to `8000`
   - Click "Deploy"

3. **Done!**
   - You get: `https://your-app.ondigitalocean.app`
   - Auto-restarts on crashes
   - Automatic HTTPS

**Cost:** ~$5-12/month

### Deploy to Render.com (Free Tier Available)

1. **Connect Repository**
   - Sign up at render.com
   - "New Web Service" → Connect repository

2. **Configure**
   ```yaml
   Build Command: docker build -t ribostruct .
   Start Command: docker run -p 8000:8000 ribostruct
   ```

3. **Deploy**
   - Free tier available (with limitations)
   - Automatic HTTPS
   - URL: `https://your-app.onrender.com`

---

## Option 2: VPS Deployment (Your Own Server)

### Prerequisites
- Ubuntu/Debian server (AWS EC2, DigitalOcean Droplet, etc.)
- Domain name (optional but recommended)
- SSH access

### Step 1: Server Setup

**Install Docker:**
```bash
# Update system
sudo apt update && sudo apt upgrade -y

# Install Docker
curl -fsSL https://get.docker.com -o get-docker.sh
sudo sh get-docker.sh

# Install Docker Compose
sudo apt install docker-compose -y
```

### Step 2: Deploy Application

**Clone and run:**
```bash
# Clone your repository
git clone <your-repo-url>
cd RiboStructMapper

# Start with Docker Compose
sudo docker-compose up -d

# Enable auto-start on reboot
sudo docker update --restart=always $(sudo docker ps -q)
```

### Step 3: Configure Firewall

```bash
# Allow HTTP/HTTPS
sudo ufw allow 80/tcp
sudo ufw allow 443/tcp
sudo ufw allow 8000/tcp  # For direct access
sudo ufw enable
```

### Step 4: Setup Reverse Proxy (nginx)

**Install nginx:**
```bash
sudo apt install nginx -y
```

**Create configuration:**
```bash
sudo nano /etc/nginx/sites-available/ribostruct
```

**Add this configuration:**
```nginx
server {
    listen 80;
    server_name your-domain.com;  # Replace with your domain

    client_max_body_size 100M;  # Allow large file uploads

    location / {
        proxy_pass http://localhost:8000;
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection 'upgrade';
        proxy_set_header Host $host;
        proxy_cache_bypass $http_upgrade;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
    }
}
```

**Enable and restart:**
```bash
sudo ln -s /etc/nginx/sites-available/ribostruct /etc/nginx/sites-enabled/
sudo nginx -t
sudo systemctl restart nginx
```

### Step 5: Setup HTTPS (Free SSL)

```bash
# Install Certbot
sudo apt install certbot python3-certbot-nginx -y

# Get SSL certificate
sudo certbot --nginx -d your-domain.com

# Auto-renewal is configured automatically
```

**Done!** Your app is now at: `https://your-domain.com`

---

## Option 3: Systemd Service (Auto-Start)

If you want the server to automatically start on boot without Docker:

**Create service file:**
```bash
sudo nano /etc/systemd/system/ribostruct.service
```

**Add configuration:**
```ini
[Unit]
Description=RiboStructMapper Web Application
After=network.target

[Service]
Type=simple
User=your-username
WorkingDirectory=/home/your-username/RiboStructMapper
Environment="PATH=/home/your-username/RiboStructMapper/.venv/bin"
ExecStart=/home/your-username/RiboStructMapper/.venv/bin/uvicorn ribostruct.web.server:app --host 0.0.0.0 --port 8000
Restart=always
RestartSec=10

[Install]
WantedBy=multi-user.target
```

**Enable and start:**
```bash
sudo systemctl daemon-reload
sudo systemctl enable ribostruct
sudo systemctl start ribostruct
sudo systemctl status ribostruct
```

---

## 🔒 Security Checklist

Before going public, ensure:

- [ ] **HTTPS enabled** (via Let's Encrypt or cloud platform)
- [ ] **Firewall configured** (only ports 80, 443 open)
- [ ] **File upload limits** set (prevent abuse)
- [ ] **Rate limiting** enabled (prevent DoS)
- [ ] **Regular backups** of job data
- [ ] **Monitoring** setup (uptime checks)

### Add Rate Limiting (nginx)

Add to your nginx config:
```nginx
limit_req_zone $binary_remote_addr zone=ribostruct:10m rate=10r/m;

location /submit_job {
    limit_req zone=ribostruct burst=5;
    proxy_pass http://localhost:8000;
    # ... other proxy settings
}
```

---

## 🌍 Domain Setup

### Point Domain to Server

**For VPS:**
1. Get your server's IP address: `curl ifconfig.me`
2. In your domain registrar (GoDaddy, Namecheap, etc.):
   - Add **A record**: `@` → Your server IP
   - Add **A record**: `www` → Your server IP
3. Wait 5-60 minutes for DNS propagation

**For Cloud Platforms:**
- They provide a URL automatically
- Or add custom domain in their settings

---

## 📊 Monitoring

### Check if Running

```bash
# Docker
sudo docker ps

# Systemd
sudo systemctl status ribostruct

# Check logs
sudo docker logs $(sudo docker ps -q --filter ancestor=ribostruct)
# OR
sudo journalctl -u ribostruct -f
```

### Setup Uptime Monitoring

Use free services like:
- UptimeRobot (https://uptimerobot.com)
- Pingdom
- StatusCake

Monitor: `https://your-domain.com/health`

---

## 💰 Cost Comparison

| Option | Cost | Difficulty | Best For |
|--------|------|------------|----------|
| **Render.com (Free)** | Free | ⭐ Easy | Testing/Personal |
| **DigitalOcean App** | $5-12/mo | ⭐ Easy | Small teams |
| **VPS + Docker** | $5-20/mo | ⭐⭐ Medium | Full control |
| **AWS/GCP** | $10-50/mo | ⭐⭐⭐ Hard | Enterprise |

---

## 🚀 Quick Deploy Commands

### For Docker on VPS:
```bash
# One-time setup
git clone <repo>
cd RiboStructMapper
sudo docker-compose up -d
sudo docker update --restart=always $(sudo docker ps -q)

# Updates
git pull
sudo docker-compose down
sudo docker-compose up -d --build
```

### For Cloud Platforms:
```bash
# Just push to GitHub
git push origin main
# Platform auto-deploys
```

---

## 🆘 Troubleshooting

**Can't access from outside?**
→ Check firewall: `sudo ufw status`

**Server stops after reboot?**
→ Enable auto-restart: `sudo docker update --restart=always <container-id>`

**Out of memory?**
→ Upgrade to larger instance or add swap space

**SSL not working?**
→ Check domain points to server: `nslookup your-domain.com`

---

## ✅ Final Checklist

- [ ] Application deployed and accessible
- [ ] HTTPS/SSL configured
- [ ] Auto-restart enabled (Docker/systemd)
- [ ] Domain configured (if using one)
- [ ] Firewall rules set
- [ ] Monitoring setup
- [ ] Backup strategy in place
- [ ] Test file upload (large files)
- [ ] Test processing with real data

---

**Your app is now accessible 24/7 from anywhere in the world!** 🌍
