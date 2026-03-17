# Contributing to RiboPrint

We welcome contributions to RiboPrint! This guide will help you get started.

## Development Setup

```bash
# Clone the repository
git clone https://github.com/NarNar420/RiboStructMapper.git
cd RiboPrint

# Create virtual environment
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Run tests
pytest tests/
```

## Project Structure

```
RiboPrint/
├── ribostruct/          # Main package
│   ├── core/            # Core processing modules
│   ├── cli/             # Command-line interface
│   └── web/             # Web application
├── tests/               # Test suite
│   ├── unit/            # Unit tests
│   ├── integration/     # Integration tests
│   └── api/             # API tests
├── docs/                # Documentation
├── data/mock/           # Test data
└── scripts/             # Utility scripts
```

## Running Tests

```bash
# All tests
pytest tests/

# Specific test category
pytest tests/unit/
pytest tests/integration/
pytest tests/api/

# With coverage
pytest tests/ --cov=ribostruct
```

## Code Style

- Follow PEP 8 guidelines
- Use type hints where appropriate
- Add docstrings to all functions and classes
- Keep functions focused and single-purpose

## Pull Request Process

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-feature`)
3. Commit your changes with a clear message (`git commit -m 'feat: add new feature'`)
4. Push to your branch (`git push origin feature/my-feature`)
5. Open a Pull Request against `main` describing what was changed and why

## Questions?

Open an issue on GitHub or contact the maintainers.
