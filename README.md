# 🧬 RiboStructMapper

**Map ribosome density data onto protein structures for visualization**

RiboStructMapper is a bioinformatics web application that integrates ribosome profiling data with protein structures. It processes genomic sequences, ribosome density profiles, and PDB structures to generate visualization-ready output files with density scores encoded in the B-factor column.

[![Python](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/)
[![FastAPI](https://img.shields.io/badge/FastAPI-0.128-green.svg)](https://fastapi.tiangolo.com/)
[![Docker](https://img.shields.io/badge/docker-ready-blue.svg)](https://www.docker.com/)

---

## ✨ Features

- 🧬 **Genomic Data Parsing** - Extracts CDS from FASTA + GTF annotations
- 🔬 **Structure Processing** - Parses PDB files and extracts amino acid sequences
- 🧮 **Sequence Alignment** - Global alignment with gap handling
- 📊 **Density Aggregation** - Nucleotide-to-amino acid mapping (mean/max/sum/median)
- ⚙️ **Offset Correction** - Ribosome position phase adjustment
- 💉 **B-Factor Injection** - Writes density scores to PDB for visualization
- 🌐 **Web Interface** - Modern dark-mode UI for easy file uploads
- 📡 **REST API** - Programmatic access for automation
- 🧹 **Auto Cleanup** - Automatic deletion of jobs older than 24 hours
- 🐳 **Docker Ready** - One-command deployment anywhere

---

## 🚀 Quick Start

### Option 1: Docker (Recommended)

```bash
# Using Docker Compose
docker-compose up -d

# Or using Docker CLI
docker build -t ribostruct .
docker run -p 8000:8000 ribostruct
```

**Access the application**: http://localhost:8000

### Option 2: Manual Installation

```bash
# Install dependencies
pip install -r requirements.txt

# Run the server
uvicorn server:app --reload
```

**Access the application**: http://127.0.0.1:8000

---

## 📋 Requirements

### Input Files

1. **PDB File** (`.pdb`) - Protein structure in PDB format
2. **FASTA File** (`.fasta`) - Genomic sequence
3. **GTF File** (`.gtf`) - Gene annotations
4. **bedGraph File** (`.bedgraph`) - Ribosome density data
5. **Offset Values** - Comma-separated integers (e.g., `0,-12,-15`)

### System Requirements

- **Python**: 3.10+
- **Memory**: 2GB+ recommended
- **Disk**: ~500MB for Docker image

---

## 🌐 Web Interface

The web interface provides an intuitive way to process your data:

1. Open http://localhost:8000 in your browser
2. Upload the 4 required files (PDB, FASTA, GTF, bedGraph)
3. Enter offset values (e.g., `0,-12`)
4. Click "Process Files"
5. Wait for processing (status updates every 2 seconds)
6. Download results as ZIP file

---

## 📡 API Reference

### POST `/submit_job`

Submit a new processing job.

**Request:**
```bash
curl -X POST "http://localhost:8000/submit_job" \
  -F "pdb_file=@protein.pdb" \
  -F "fasta_file=@genome.fasta" \
  -F "gtf_file=@annotation.gtf" \
  -F "density_file=@density.bedgraph" \
  -F "offsets=0,-12,-15"
```

**Response:**
```json
{
  "job_id": "a3b8d1b6-0b3b-4b1a-9c1a-1a2b3c4d5e6f",
  "status": "queued",
  "message": "Job created and queued for processing. Job ID: a3b8d1b6..."
}
```

### GET `/status/{job_id}`

Check job status.

**Request:**
```bash
curl "http://localhost:8000/status/a3b8d1b6-0b3b-4b1a-9c1a-1a2b3c4d5e6f"
```

**Response:**
```json
{
  "job_id": "a3b8d1b6-0b3b-4b1a-9c1a-1a2b3c4d5e6f",
  "status": "completed",
  "message": null
}
```

**Status Values:**
- `uploaded` - Files uploaded successfully
- `queued` - Job queued for processing
- `processing` - Pipeline running
- `completed` - Processing finished successfully
- `failed` - Error occurred (see message)

### GET `/download/{job_id}`

Download results as ZIP file.

**Request:**
```bash
curl "http://localhost:8000/download/a3b8d1b6-0b3b-4b1a-9c1a-1a2b3c4d5e6f" \
  -o results.zip
```

**Response:**
ZIP archive containing all output PDB files (one per offset value):
- `output_offset_0.pdb`
- `output_offset_-12.pdb`
- `output_offset_-15.pdb`

### GET `/health`

Server health check.

**Request:**
```bash
curl "http://localhost:8000/health"
```

**Response:**
```json
{
  "status": "healthy",
  "jobs_directory": "jobs",
  "jobs_directory_exists": true
}
```

---

## 📄 File Format Guide

### PDB Format

Standard Protein Data Bank format with ATOM records:

```pdb
ATOM      1  N   MET A   1      20.154  30.229  25.123  1.00 10.50           N
ATOM      2  CA  MET A   1      21.289  29.345  25.456  1.00 11.20           C
```

**Note:** The B-factor column (columns 61-66) will be overwritten with density scores.

### bedGraph Format

Tab-separated values with chromosome, start, end, and density:

```bedgraph
mock_chrom	1	2	10.5
mock_chrom	2	3	12.3
mock_chrom	3	4	15.7
```

**Requirements:**
- Chromosome name must match FASTA header
- Coordinates are 0-based, half-open intervals
- Density values are floating-point numbers

### FASTA Format

Standard genomic sequence format:

```fasta
>mock_chrom
ATGAAAACCATAATAGCTCTGTCTTACATATTCTGCCTGGTGTTC
```

### GTF Format

Gene annotation format with CDS features:

```gtf
mock_chrom	test	CDS	1	45	.	+	0	gene_id "MOCK1"; transcript_id "MOCK1.1";
```

**Requirements:**
- Feature type must be `CDS`
- Must include `gene_id` and `transcript_id` attributes

---

## 🧪 Testing

```bash
# Run core module tests
pytest test_parser.py
pytest test_alignment.py
pytest test_processor.py
pytest test_injector.py

# Test server endpoints
python test_server.py

# Test cleanup functionality
python test_cleanup.py

# Test complete workflow
python test_complete_workflow.py
```

---

## 🏗️ Architecture

```
RiboStructMapper/
├── Core Modules
│   ├── parser.py         # FASTA/GTF/PDB/bedGraph parsing
│   ├── alignment.py      # Sequence translation and alignment
│   ├── processor.py      # Density aggregation and offset application
│   └── injector.py       # B-factor injection into PDB
│
├── Web Application
│   ├── server.py         # FastAPI server with REST API
│   └── static/
│       └── index.html    # Web interface
│
├── Integration
│   └── main_cli.py       # Command-line pipeline orchestration
│
└── Tests
    ├── test_*.py         # Unit and integration tests
    └── mock_data/        # Test fixtures
```

---

## 🔧 Configuration

### Offset Values

Offset values adjust the ribosome position on the mRNA:

- **0** - No adjustment (raw density)
- **-12** - Shift 12 nucleotides upstream (common for P-site)
- **-15** - Shift 15 nucleotides upstream (alternative P-site)

Multiple offsets generate separate output files for comparison.

### Aggregation Methods

The pipeline supports multiple density aggregation strategies:

- **mean** (default) - Average density across codon
- **max** - Maximum density value
- **sum** - Total density
- **median** - Median density value

Used in CLI mode: `--aggregation-method mean`

---

## 🐳 Docker Deployment

See [DOCKER.md](DOCKER.md) for detailed deployment instructions including:

- Production deployment
- Reverse proxy configuration
- Cloud platform deployment (AWS, GCP, Azure)
- Volume management
- Troubleshooting

---

## 📊 Example Workflow

```bash
# 1. Start the server
uvicorn server:app --reload

# 2. Submit a job via API
JOB_ID=$(curl -X POST "http://localhost:8000/submit_job" \
  -F "pdb_file=@data/protein.pdb" \
  -F "fasta_file=@data/genome.fasta" \
  -F "gtf_file=@data/annotation.gtf" \
  -F "density_file=@data/density.bedgraph" \
  -F "offsets=0,-12" | jq -r '.job_id')

# 3. Monitor status
watch -n 2 curl -s "http://localhost:8000/status/$JOB_ID" | jq

# 4. Download results when complete
curl "http://localhost:8000/download/$JOB_ID" -o results.zip

# 5. Extract and visualize in PyMOL/ChimeraX
unzip results.zip
pymol output_offset_-12.pdb
```

---

## 🛠️ Development

```bash
# Clone the repository
git clone <repository_url>
cd RiboStructMapper

# Create virtual environment
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Run tests
pytest

# Start development server
uvicorn server:app --reload --log-level debug
```

---

## 📝 Citation

If you use RiboStructMapper in your research, please cite:

```
[Your citation information here]
```

---

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

---

## 📜 License

[Your license information here]

---

## 🆘 Support

For issues, questions, or feature requests:

- Open an issue on GitHub
- See [SERVER_USAGE.md](SERVER_USAGE.md) for detailed usage instructions
- Check [DOCKER.md](DOCKER.md) for deployment help

---

## 👥 Authors

[Your author information here]

---

**Built with**: Python • FastAPI • Biopython • Docker

**Powered by**: Uvicorn • NumPy • Pandas • SciPy
