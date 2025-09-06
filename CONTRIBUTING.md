# Contributing to Codoff

Thank you for your interest in contributing to the Codoff project! This document provides guidelines and instructions for contributing to the codebase.

## Table of Contents

- [Development Setup](#development-setup)
- [Running Tests](#running-tests)
- [Code Style Guidelines](#code-style-guidelines)
- [Testing Guidelines](#testing-guidelines)
- [Submitting Changes](#submitting-changes)
- [Project Structure](#project-structure)

## Development Setup

### Prerequisites

- Python 3.8 or higher
- pip (Python package installer)

### Installation

1. **Fork the repository** on GitHub:
   - Go to the [codoff repository](https://github.com/Kalan-Lab/codoff)
   - Click the "Fork" button in the top-right corner
   - This creates a copy of the repository in your GitHub account

2. **Clone your fork** to your local machine:
```bash
git clone https://github.com/YOUR_USERNAME/codoff.git
cd codoff
```

3. **Add the upstream repository** as a remote:
```bash
git remote add upstream https://github.com/Kalan-Lab/codoff.git
```

4. **Install the package in development mode**:
```bash
pip install -e .
```

5. **Install development dependencies**:
```bash
pip install pytest pytest-cov pytest-mock numpy scipy biopython matplotlib seaborn pyrodigal
```

## Running Tests

The project uses pytest for testing. All tests are located in the `tests/` directory.

### Running All Tests

```bash
# Run all tests with verbose output
python3 -m pytest tests/ -v

# Run all tests with coverage report
python3 -m pytest tests/ --cov=src/codoff --cov-report=html

# Run tests and show coverage in terminal
python3 -m pytest tests/ --cov=src/codoff --cov-report=term-missing
```

### Running Specific Test Categories

```bash
# Run only unit tests
python3 -m pytest tests/unit/ -v

# Run only integration tests
python3 -m pytest tests/integration/ -v

# Run tests for specific functionality
python3 -m pytest tests/unit/test_caching.py -v
python3 -m pytest tests/unit/test_codon_analysis.py -v
python3 -m pytest tests/unit/test_threading.py -v
```

### Running Tests with Different Options

```bash
# Run tests in parallel (if pytest-xdist is installed)
python3 -m pytest tests/ -n auto

# Run tests and stop on first failure
python3 -m pytest tests/ -x

# Run tests with detailed output
python3 -m pytest tests/ -s

# Run only tests matching a pattern
python3 -m pytest tests/ -k "test_caching"
```

### Test Configuration

The project uses `pytest.ini` for configuration:
- Test discovery: `tests/` directory
- Coverage reporting: `src/codoff` module
- Custom markers: `unit`, `integration`
- Coverage threshold: 80%

## Code Style Guidelines

### Python Style

The project follows PEP 8 guidelines with the following specific requirements:

- **Line length**: Maximum 100 characters
- **Indentation**: 4 spaces (no tabs)
- **Function naming**: `snake_case`
- **Class naming**: `PascalCase`
- **Constants**: `UPPER_CASE`

### Type Hints

All functions must include type hints:

```python
def example_function(param1: str, param2: int) -> List[str]:
    """Example function with type hints."""
    return [param1] * param2
```

### Docstrings

Use Google-style docstrings for all functions:

```python
def example_function(param1: str, param2: int) -> List[str]:
    """
    Brief description of what the function does.
    
    Longer description if needed, explaining the function's purpose,
    behavior, or important implementation details.
    
    Args:
        param1: Description of the first parameter
        param2: Description of the second parameter
        
    Returns:
        Description of what the function returns
        
    Raises:
        ValueError: Description of when this exception is raised
    """
    return [param1] * param2
```

### Code Formatting

Use the following tools to maintain consistent formatting:

```bash
# Check code style
python3 -m flake8 src/codoff/ --max-line-length=100 --ignore=E501,W503

# Auto-format code (if black is installed)
black src/codoff/

# Sort imports (if isort is installed)
isort src/codoff/
```

## Testing Guidelines

### Writing Tests

1. **Test Structure**: Follow the `tests/unit/` and `tests/integration/` structure
2. **Test Naming**: Use descriptive names like `test_function_name_scenario`
3. **Test Coverage**: Aim for high test coverage of core functionality
4. **Test Data**: Use fixtures in `tests/conftest.py` for reusable test data

### Test Categories

- **Unit Tests** (`tests/unit/`): Test individual functions and methods
- **Integration Tests** (`tests/integration/`): Test complete workflows and interactions

### Example Test Structure

```python
import pytest
from codoff.codoff import function_to_test

class TestFunctionName:
    """Test class for function_name."""
    
    def test_function_name_valid_input(self):
        """Test function with valid input."""
        result = function_to_test("valid_input")
        assert result == expected_output
    
    def test_function_name_invalid_input(self):
        """Test function with invalid input."""
        with pytest.raises(ValueError):
            function_to_test("invalid_input")
```

### Running Tests Before Committing

Always run tests before submitting changes:

```bash
# Run all tests
python3 -m pytest tests/ -v

# Check code style
python3 -m flake8 src/codoff/ --max-line-length=100 --ignore=E501,W503

> [!NOTE] 
> F824 warnings about global `_codon_extraction_cache` and `_codon_counts_cache` are false positives - these variables are actively used for caching functionality throughout the codebase. These warnings can be safely ignored.

## Submitting Changes

### Before Submitting

1. **Run all tests** to ensure nothing is broken
2. **Check code style** with flake8
3. **Update documentation** if needed
4. **Add tests** for new functionality
5. **Update CHANGELOG.md** if applicable

### Development Workflow

1. **Keep your fork up to date**:
```bash
git fetch upstream
git checkout main
git merge upstream/main
```

2. **Create a feature branch**:
```bash
git checkout -b feature/your-feature-name
```

3. **Make your changes** and commit them:
```bash
git add .
git commit -m "feat: add new caching functionality"
```

4. **Push to your fork**:
```bash
git push origin feature/your-feature-name
```

### Commit Message Format

Use clear, descriptive commit messages:

```
feat: add new caching functionality
fix: resolve threading issue in simulation
docs: update function docstrings
test: add unit tests for codon analysis
refactor: improve code organization
```

### Pull Request Process

1. **Create a Pull Request** on GitHub:
   - Go to your fork on GitHub
   - Click "New Pull Request"
   - Select your feature branch
   - Write a clear description of your changes
   - Link any related issues

2. **Ensure all checks pass**:
   - All tests must pass
   - Code style checks must pass
   - Documentation is updated if needed

3. **Respond to feedback**:
   - Address any review comments
   - Make additional commits if needed
   - Keep the PR up to date with main

## Project Structure

```
codoff/
├── src/codoff/           # Main source code
│   ├── __init__.py
│   ├── codoff.py         # Core analysis functions
│   └── util.py           # Utility functions
├── tests/                # Test suite
│   ├── unit/             # Unit tests
│   ├── integration/      # Integration tests
│   ├── fixtures/         # Test data and fixtures
│   └── conftest.py       # Pytest configuration
├── bin/                  # Executable scripts
├── docs/                 # Documentation
├── pyproject.toml        # Project configuration
├── pytest.ini           # Pytest configuration
└── README.md            # Project overview
```

## Key Dependencies

- **Core**: `numpy`, `scipy`, `biopython`
- **Visualization**: `matplotlib`, `seaborn`
- **Gene Calling**: `pyrodigal`
- **Testing**: `pytest`, `pytest-cov`, `pytest-mock`

## Getting Help

- Check existing issues and pull requests
- Review the codebase, wiki, and existing tests
- Ask questions in discussions or issues

## License

By contributing to this project, you agree that your contributions will be licensed under the same license as the project.

---

Thank you for contributing to codoff! 🧬
