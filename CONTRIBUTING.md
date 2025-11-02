# Contributing to VCC-project

Thank you for your interest in contributing! 🎉

## How to Contribute

### Reporting Bugs

1. Check existing [Issues](https://github.com/ACTN3Bioinformatics/VCC-project/issues)
2. Create new issue with:
   - Clear description
   - Steps to reproduce
   - Expected vs actual behavior
   - System info (OS, Python version, etc.)
   - Relevant logs

### Suggesting Features

1. Open an issue with tag `enhancement`
2. Describe the feature and use case
3. Discuss implementation approach

### Submitting Code

1. **Fork the repository**
2. **Create feature branch**:
```bash
   git checkout -b feature/amazing-feature
```

3. **Make changes**:
   - Follow existing code style
   - Add tests for new features
   - Update documentation
   - Run tests locally

4. **Commit changes**:
```bash
   git commit -m "Add amazing feature"
```

5. **Push to branch**:
```bash
   git push origin feature/amazing-feature
```

6. **Open Pull Request**

## Development Setup
```bash
# Clone your fork
git clone https://github.com/YOUR_USERNAME/VCC-project.git
cd VCC-project

# Add upstream remote
git remote add upstream https://github.com/ACTN3Bioinformatics/VCC-project.git

# Create environment
conda env create -f environment.yml
conda activate vcc2025

# Install development dependencies
pip install -e .
```

## Code Style

- **Python**: Follow PEP 8
- **Formatting**: Use `black` for formatting
- **Linting**: Use `flake8`
- **Imports**: Use `isort`
```bash
# Format code
black scripts/

# Check linting
flake8 scripts/

# Sort imports
isort scripts/
```

## Testing
```bash
# Run all tests
pytest tests/

# Run with coverage
pytest --cov=scripts tests/

# Run specific test
pytest tests/test_qc.py::test_filter_cells -v
```

## Documentation

- Update `docs/` for major changes
- Use docstrings for functions (Google style)
- Update README.md if needed
- Add examples for new features

## Pull Request Checklist

- [ ] Code follows project style
- [ ] Tests pass locally
- [ ] New tests added for new features
- [ ] Documentation updated
- [ ] Commit messages are clear
- [ ] No merge conflicts

## Code Review Process

1. Maintainer reviews PR
2. Feedback addressed
3. CI/CD checks pass
4. Approved and merged

## Questions?

Open an issue or contact maintainers directly.

Thank you for contributing! 🚀