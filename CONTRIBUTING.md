# Contributing to codoff

Thank you for your interest in contributing to codoff! This document provides guidelines for contributing to the project, with a focus on the testing structure and development workflow.

## Development Setup

### Prerequisites

- Python 3.10 or higher
- Git

### Installation for Development

1. **Fork the repository** on GitHub:
   - Go to https://github.com/Kalan-Lab/codoff
   - Click the "Fork" button in the top-right corner
   - This creates your own copy of the repository

2. **Clone your fork** and set up the upstream remote:
   ```bash
   git clone https://github.com/YOUR_USERNAME/codoff.git
   cd codoff
   git remote add upstream https://github.com/Kalan-Lab/codoff.git
   ```

3. **Create a new branch** for your development work:
   ```bash
   git checkout -b feature/your-feature-name
   # or
   git checkout -b bugfix/issue-description
   ```
   
   **Branch naming conventions**:
   - `feature/description` - for new features (e.g., `feature/codon-caching`)
   - `bugfix/description` - for bug fixes (e.g., `bugfix/warning-messages`)
   - `hotfix/description` - for urgent fixes (e.g., `hotfix/critical-error`)
   - `docs/description` - for documentation updates (e.g., `docs/contributing-guide`)

4. **Create and activate a conda environment**:
   ```bash
   conda env create -f codoff_env.yml -n codoff_env
   conda activate codoff_env
   ```

5. **Install the package in development mode**:
   ```bash
   pip install -e .
   ```

6. **Keep your fork up to date**:
   ```bash
   git fetch upstream
   git checkout main
   git merge upstream/main
   ```

## Testing Structure

The project uses a comprehensive test suite organized into specialized test modules. All tests are located in the `tests/` directory and follow a structured approach to ensure code quality and reliability.

### Test Organization

The test suite is organized into the following categories:

#### 1. Core Functionality Tests (`test_codoff.py`)
- **Purpose**: Tests the main codoff functionality and core algorithms
- **Focus**: Codon caching, simulation logic, and basic statistical calculations
- **Key Test Classes**:
  - `TestCodonCaching`: Tests codon count caching functionality
  - `TestProportionalSampling`: Tests proportional sampling logic
  - `TestWarningMessages`: Tests warning message generation
  - `TestIntegration`: Tests simulation consistency and integration
  - `TestEdgeCases`: Tests edge cases like empty gene lists and single codon types

#### 2. Integration Tests (`test_integration.py`)
- **Purpose**: Tests complete workflows with realistic data
- **Focus**: End-to-end functionality with real-world scenarios
- **Key Test Classes**:
  - `TestIntegrationWorkflow`: Tests complete workflow with realistic genome data and edge cases

#### 3. antiSMASH Integration Tests (`test_antismash_codoff_integration.py`)
- **Purpose**: Tests integration with antiSMASH and caching functionality
- **Focus**: Data structure integrity, background calculations, and parameter handling
- **Key Test Classes**:
  - `TestAntismashCodoffIntegration`: Tests caching data structure, background calculation methods, and function signatures

#### 4. Warning System Tests (`test_warnings.py`)
- **Purpose**: Tests warning messages and error handling
- **Focus**: User feedback and error reporting
- **Key Test Classes**:
  - `TestWarningMessages`: Tests warning messages for missing locus tags, coordinate warnings, and warning conditions
  - `TestWarningMessageContent`: Tests warning message format and output capture

#### 5. Seed Reproducibility Tests (`test_seed_reproducibility.py`)
- **Purpose**: Tests that random seed parameter produces reproducible results
- **Focus**: Reproducibility and deterministic behavior
- **Key Test Classes**:
  - `TestSeedReproducibility`: Tests that identical seeds produce identical results and different seeds can produce different results

#### 6. Caching Functionality Tests (`test_caching.py`)
- **Purpose**: Tests genome data caching for performance optimization
- **Focus**: Data structure integrity and caching consistency
- **Key Test Classes**:
  - `TestCachingFunctionality`: Tests extract_genome_codon_data() structure, cached vs non-cached consistency, and simulation count parameter

### Running Tests

#### Using the Test Runner (Recommended)

The project includes a dedicated test runner script:

```bash
python run_tests.py
```

This will:
- Discover all test files in the `tests/` directory
- Run all tests with verbose output
- Return appropriate exit codes (0 for success, 1 for failure)

#### Using pytest (Alternative)

You can also use pytest directly:

```bash
# Run all tests
python -m pytest tests/ -v

# Run with coverage
python -m pytest tests/ --cov=src/codoff --cov-report=html

# Run specific test file
python -m pytest tests/test_codoff.py -v
```

**Note**: The project is configured to work with both unittest and pytest. No additional configuration files are needed.

#### Running Individual Test Modules

You can also run specific test modules:

```bash
# Run core functionality tests
python -m unittest tests.test_codoff

# Run integration tests
python -m unittest tests.test_integration

# Run antiSMASH integration tests
python -m unittest tests.test_antismash_codoff_integration

# Run warning system tests
python -m unittest tests.test_warnings

# Run seed reproducibility tests
python -m unittest tests.test_seed_reproducibility

# Run caching functionality tests
python -m unittest tests.test_caching

# Run utility function tests
python -m unittest tests.test_utils

# Run default behavior tests
python -m unittest tests.test_new_default
```

#### Running Specific Test Classes or Methods

```bash
# Run a specific test class
python -m unittest tests.test_codoff.TestCodonCaching

# Run a specific test method
python -m unittest tests.test_codoff.TestCodonCaching.test_gene_codons_structure
```

### Test Development Guidelines

#### Writing New Tests

When adding new tests, follow these guidelines:

1. **Test Organization**: Place tests in the appropriate module based on functionality:
   - Core algorithm tests → `test_codoff.py`
   - End-to-end workflow tests → `test_integration.py`
   - Processing antiSMASH results tests → `test_antismash_codoff_integration.py`
   - Warning/error tests → `test_warnings.py`
   - Reproducibility tests → `test_seed_reproducibility.py`
   - Caching tests → `test_caching.py`
   - Utility function tests → `test_utils.py`
   - Default behavior tests → `test_new_default.py`

2. **Test Naming**: Use descriptive test method names that explain what is being tested:
   ```python
   def test_codon_frequency_calculation_with_empty_gene_list(self):
       """Test that codon frequency calculation handles empty gene lists correctly."""
   ```

3. **Test Documentation**: Include docstrings that explain the test purpose:
   ```python
   def test_realistic_genome_simulation(self):
       """Test simulation with realistic genome data to ensure proper handling of real-world scenarios."""
   ```

4. **Test Data**: Use realistic test data that reflects actual biological scenarios when possible.

5. **Assertions**: Use specific assertions that test the exact behavior expected:
   ```python
   self.assertEqual(total_sampled, expected_count)
   self.assertAlmostEqual(proportion, expected_proportion, places=5)
   ```

#### Test Structure Requirements

- All test classes must inherit from `unittest.TestCase`
- Test methods must start with `test_`
- Use descriptive class and method names
- Include comprehensive docstrings
- Follow the existing code style and PEP 8 guidelines

### Code Quality Standards

#### Python Style Guidelines

- Follow PEP 8 style guidelines
- Use type hints for function parameters and return values
- Keep functions brief (preferably under 14 lines)
- Use meaningful variable and function names
- Add comprehensive docstrings for all public functions using NumPy-style format:
  ```python
  def function_name(param1: str, param2: int = 10) -> Dict[str, Any]:
      """
      Brief one-line description.
      
      Parameters
      ----------
      param1 : str
          Description of parameter
      param2 : int, optional
          Description with default value, by default 10
      
      Returns
      -------
      Dict[str, Any]
          Description of return value
      
      Notes
      -----
      Additional notes if needed.
      """
  ```

#### Simulation Processing

When working with simulation code:

- **Sequential Contiguous-Window Sampling**: The tool uses sequential contiguous-window sampling by default (as of v1.2.3), which randomly selects genomic windows of the same size as the focal region for more biologically realistic null distributions
- **Random Seeding**: Simulations use fixed random seeds (default: 42) for reproducible results. Users can specify custom seeds via `--seed/-x` parameter
- **Coordinate Information**: Genomic coordinates (gene positions and scaffold lengths) are always extracted and cached for efficient sequential sampling
- **Discordance Percentile**: Results report a discordance percentile (not p-value) indicating how unusual the focal region's codon usage is compared to similarly sized genomic windows
- **Testing**: Always test simulation consistency and reproducibility when modifying simulation code
- **Performance**: Sequential processing with coordinate-based sampling provides optimal performance
- **Compatibility**: Ensure changes work with both `codoff` and `antismash_codoff` scripts

#### Code Organization

- Place new functionality in appropriate modules within `src/codoff/`
- Update tests when adding new features
- Ensure all new code is covered by tests
- Maintain backward compatibility when possible

#### CLI Parameter Guidelines

When adding new CLI parameters:

- **Consistency**: Ensure parameters work in both `codoff` and `antismash_codoff` scripts
- **Defaults**: Use sensible defaults (e.g., 10000 simulations, seed=42, sequential sampling enabled)
- **Documentation**: Update help text, README.md, and docstrings
- **Testing**: Add tests for new parameters in the appropriate test modules
- **Backward Compatibility**: New parameters should be optional with sensible defaults
- **Coordinate Requirements**: Sequential sampling requires genomic coordinates, which are automatically extracted from GenBank files or generated via pyrodigal for FASTA files

### Pull Request Process

1. **Follow the development setup** above to fork and create a feature branch
2. **Make your changes** on your feature branch
3. **Write tests** for any new functionality
4. **Ensure all tests pass** using `python run_tests.py`
5. **Update documentation** if needed
6. **Commit your changes** with clear commit messages:
   ```bash
   git add .
   git commit -m "Add feature: brief description of changes"
   ```
7. **Push your branch** to your fork:
   ```bash
   git push origin feature/your-feature-name
   ```
8. **Submit a pull request** on GitHub with:
   - A clear title describing the changes
   - A detailed description of what was changed and why
   - Reference to any related issues
   - Screenshots or examples if applicable

### Test Coverage

The test suite provides comprehensive coverage across multiple specialized modules:
- Core algorithmic functionality and codon caching
- Edge cases and error conditions
- Integration scenarios with realistic data
- User-facing features and warnings
- Statistical accuracy and precision
- antiSMASH integration and genome data caching
- Seed-based reproducibility
- Sequential contiguous-window sampling
- Coordinate extraction and usage

### Continuous Integration

The project uses GitHub Actions for continuous integration, which automatically runs the test suite on:
- Pull requests
- Pushes to main branch
- Release tags

All tests must pass before code can be merged.

## Important Terminology

### Discordance Percentile (not P-value)

As of version 1.2.3, codoff reports a **discordance percentile** rather than an empirical P-value. This is an important distinction:

- **Discordance Percentile**: The proportion of simulations that show codon usage as or more discordant than the observed focal region, expressed as a percentage (0-100)
- **Interpretation**: Lower percentiles indicate more unusual/discordant codon usage
- **Example**: A percentile of 5.0 means the focal region is within the top 5% most discordant regions

When contributing code or documentation:
- Use "discordance percentile" (not "p-value" or "empirical p-value")
- The internal variable is `empirical_freq` (a proportion from 0-1)
- Output displays it as "Discordance Percentile" (multiplied by 100)

## Getting Help

If you have questions about contributing or need help with the codebase:

1. Check the existing issues on GitHub
2. Review the README.md and wiki for examples of how functionality is used
3. Open a new issue with specific questions

## License

By contributing to codoff, you agree that your contributions will be licensed under the same BSD 3-Clause License that covers the project.
