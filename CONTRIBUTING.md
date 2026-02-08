# Contributing to RiboStructMapper

We welcome contributions to RiboStructMapper! This guide will help you get started.

## Development Setup

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
pytest tests/
```

## Project Structure

```
RiboStructMapper/
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
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Run tests to ensure everything works
5. Commit your changes (`git commit -m 'Add amazing feature'`)
6. Push to the branch (`git push origin feature/amazing-feature`)
7. Open a Pull Request

## Questions?

Open an issue on GitHub or contact the maintainers.
