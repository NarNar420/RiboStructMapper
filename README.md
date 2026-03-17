# 🧬 RiboPrint

**Map ribosome density data onto protein structures for visualization**

RiboPrint is a bioinformatics pipeline that integrates ribosome profiling data with protein structures. It processes genomic sequences, ribosome density profiles, and PDB structures to generate visualization-ready output files with density scores encoded in the B-factor column.

[![Python](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/)
[![FastAPI](https://img.shields.io/badge/FastAPI-0.128-green.svg)](https://fastapi.tiangolo.com/)
[![Docker](https://img.shields.io/badge/docker-ready-blue.svg)](https://www.docker.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

---

## 📚 Documentation

For comprehensive documentation, see the [`docs/`](docs/) directory:

- **[Full Documentation](docs/README.md)** - Complete usage guide and API reference
- **[Docker Deployment](docs/DOCKER.md)** - Docker setup and deployment instructions
- **[Server Usage](docs/SERVER_USAGE.md)** - Detailed server usage guide
- **[Specifications](docs/specifications.md)** - Technical specifications and architecture

---

## 🚀 Quick Start

### Using Docker (Recommended)

```bash
docker-compose up -d
```

**Access the application**: http://riboprint.marx-group.edu:8000

### Manual Installation

```bash
pip install -r requirements.txt
uvicorn ribostruct.web.server:app --reload
```

**Access the application**: http://riboprint.marx-group.edu:8000

---

## 🏗️ Project Structure

```
RiboPrint/
├── ribostruct/              # Main package
│   ├── core/                # Core processing modules
│   │   ├── parser.py        # Input parsing
│   │   ├── alignment.py     # Sequence alignment
│   │   ├── processor.py     # Density processing
│   │   └── injector.py      # B-factor injection
│   ├── cli/                 # Command-line interface
│   │   └── main.py          # CLI orchestration
│   └── web/                 # Web application
│       ├── server.py        # FastAPI server
│       └── static/          # Frontend assets
│
├── tests/                   # Test suite
│   ├── unit/                # Unit tests
│   ├── integration/         # Integration tests
│   └── api/                 # API tests
│
├── docs/                    # Documentation
├── data/                    # Data files
│   └── mock/                # Mock/test data
├── scripts/                 # Utility scripts
└── jobs/                    # Runtime job storage
```

---

## ✨ Features

- 🧬 **Genomic Data Parsing** - Extracts CDS from FASTA sequence
- 🔬 **Structure Processing** - Parses PDB files and extracts amino acid sequences
- 🧮 **Sequence Alignment** - Global alignment with gap handling
- 📊 **Density Aggregation** - Nucleotide-to-amino acid mapping
- ⚙️ **Offset Correction** - Ribosome position phase adjustment
- 💉 **B-Factor Injection** - Writes density scores to PDB for visualization
- 🌐 **Web Interface** - Clean academic White/Green UI layout
- 🔎 **Interactive 3D Viewer** - Embedded browser viewing of mapped structures with dynamic color gradients
- 📡 **REST API** - Programmatic access
- 🐳 **Docker Ready** - One-command deployment

---

## 🧪 Testing

```bash
# Run all tests
pytest tests/

# Run specific test suites
pytest tests/unit/
pytest tests/integration/
pytest tests/api/
```

---

## 📖 Usage

### Web Interface

1. Open http://riboprint.marx-group.edu:8000
2. Upload required files (PDB, FASTA, bedGraph)
3. Enter offset values (e.g., `0,12`)
4. Click "Process Files"
5. Download results as ZIP

### Command-Line

```bash
python -m ribostruct.cli.main
```

### Python API

```python
from ribostruct.cli.main import run_pipeline

output_files = run_pipeline(
    pdb_path="protein.pdb",
    fasta_path="genome.fasta",
    bedgraph_path="density.bedgraph",
    offsets=[0, 12, 15]
)
```

---

## 🛠️ Development

```bash
# Clone repository
git clone https://github.com/NarNar420/RiboStructMapper.git
cd RiboPrint

# Create virtual environment
python -m venv .venv
source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Run development server
uvicorn ribostruct.web.server:app --reload --log-level debug
```

---

## 📜 License

This project is licensed under the **MIT License**. For full details, see the [LICENSE](LICENSE) file. 

The MIT License is a permissive free software license that allows anyone to do almost anything they want with the project (including commercial use, modification, and distribution), as long as they include the original copyright and license notice.

---

**Authors**: Einar Alhassid Gutkin, Ailie Marx (Marx Structural Biology Group)

---

**Built with**: Python • FastAPI • Biopython • Docker
